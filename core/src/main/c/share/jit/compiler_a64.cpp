/*******************************************************************************
 *     ___                  _   ____  ____
 *    / _ \ _   _  ___  ___| |_|  _ \| __ )
 *   | | | | | | |/ _ \/ __| __| | | |  _ \
 *   | |_| | |_| |  __/\__ \ |_| |_| | |_) |
 *    \__\_\\__,_|\___||___/\__|____/|____/
 *
 *  Copyright (c) 2014-2019 Appsicle
 *  Copyright (c) 2019-2024 QuestDB
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 ******************************************************************************/
#include "common.h"
#include "compiler.h"
#include "neon.h"
#include <asmjit/asmjit.h>
#include <asmjit/a64.h>
#include <iostream>

using namespace asmjit;

struct JitErrorHandler : public ErrorHandler {
    JitErrorHandler()
            : error(ErrorCode::kErrorOk) {}

    void handleError(Error err, const char *msg, BaseEmitter * /*origin*/) override {
        error = err;
        message.assign(msg);
    }

    asmjit::Error error;
    asmjit::String message;
};

struct JitGlobalContext {
    //rt allocator is thread-safe
    JitRuntime rt;
};

static JitGlobalContext gGlobalContext;

using CompiledFn = int64_t (*)(int64_t *cols, int64_t cols_count,
                               int64_t *varsize_indexes,
                               int64_t *vars, int64_t vars_count,
                               int64_t *rows, int64_t rows_count,
                               int64_t rows_start_offset);

struct Function {
    explicit Function(a64::Compiler &cc)
            : c(cc), zone(4094 - Zone::kBlockOverhead), allocator(&zone) {
        values.init(&allocator);

    };

    void compile(const instruction_t *istream, size_t size, uint32_t options) {
        enum type_size : uint32_t {
            scalar = 0,
            single_size = 1,
            mixed_size = 3,
        };

        uint32_t type_size = (options >> 1) & 7; // 0 - 1B, 1 - 2B, 2 - 4B, 3 - 8B, 4 - 16B
        uint32_t exec_hint = (options >> 4) & 3; // 0 - scalar, 1 - single size type, 2 - mixed size types, ...
        bool null_check = (options >> 6) & 1; // 1 - with null check
        int unroll_factor = 1;
        scalar_loop(istream, size, null_check, unroll_factor);
    };

    void scalar_tail(const instruction_t *istream, size_t size, bool null_check, const a64::Gp &stop, int unroll_factor = 1) {
        Label l_loop = c.newNamedLabel("loop");
        Label l_exit = c.newNamedLabel("exit");

        c.cmp(input_index, stop);
        c.b_gt(l_exit);

        c.bind(l_loop);

        for (int i = 0; i < unroll_factor; ++i) {
            questdb::neon::emit_code(c, istream, size, values, null_check, data_ptr, varsize_aux_ptr, vars_ptr, input_index);

            a64::Gp adjusted_id = c.newGp(TypeId::kInt64);
            c.add(adjusted_id, input_index, rows_id_start_offset);
            c.str(adjusted_id, ptr(rows_ptr, output_index, arm::Shift(arm::ShiftOp::kLSL, 3)));

            auto mask = values.pop();
            c.and_(mask.gp(), mask.gp(), 1);
            c.add(output_index, output_index, mask.gp().r64());
        }
        c.add(input_index, input_index, unroll_factor);

        c.cmp(input_index, stop);
        c.b_lt(l_loop); // input_index < stop
        c.bind(l_exit);
    }

    void scalar_loop(const instruction_t *istream, size_t size, bool null_check, int unroll_factor = 1) {
        if(unroll_factor > 1) {
            a64::Gp stop = c.newInt64("stop");
            c.sub(stop, stop, unroll_factor - 1);
            scalar_tail(istream, size, null_check, stop, unroll_factor);
            scalar_tail(istream, size, null_check, rows_size, 1);
        } else {
            scalar_tail(istream, size, null_check, rows_size, 1);
        }
        c.ret(output_index);
    }

    void begin_fn() {
        c.addFunc(FuncSignature::build<int64_t, int64_t *, int64_t, int64_t *, int64_t, int64_t *, int64_t, int64_t *, int64_t, int64_t>(
            CallConvId::kCDecl));
        data_ptr = c.newIntPtr("data_ptr");
        data_size = c.newInt64("data_size");

        c.setArg(0, data_ptr);
        c.setArg(1, data_size);

        varsize_aux_ptr = c.newIntPtr("varsize_aux_ptr");

        c.setArg(2, varsize_aux_ptr);

        vars_ptr = c.newIntPtr("vars_ptr");
        vars_size = c.newInt64("vars_size");

        c.setArg(3, vars_ptr);
        c.setArg(4, vars_size);

        rows_ptr = c.newIntPtr("rows_ptr");
        rows_size = c.newInt64("rows_size");

        c.setArg(5, rows_ptr);
        c.setArg(6, rows_size);

        rows_id_start_offset = c.newInt64("rows_id_start_offset");
        c.setArg(7, rows_id_start_offset);

        input_index = c.newInt64("input_index");
        c.mov(input_index, 0);

        output_index = c.newInt64("output_index");
        c.mov(output_index, 0);
    }

    void end_fn() {
        c.endFunc();
    }

    a64::Compiler &c;

    Zone zone;
    ZoneAllocator allocator;
    ZoneStack<jit_value_t> values;

    a64::Gp data_ptr;
    a64::Gp data_size;
    a64::Gp varsize_aux_ptr;
    a64::Gp vars_ptr;
    a64::Gp vars_size;
    a64::Gp rows_ptr;
    a64::Gp rows_size;
    a64::Gp input_index;
    a64::Gp output_index;
    a64::Gp rows_id_start_offset;
};


void fillJitErrorObject(JNIEnv *e, jobject error, uint32_t code, const char *msg) {

    if (!msg) {
        return;
    }

    jclass errorClass = e->GetObjectClass(error);
    if (errorClass) {
        jfieldID fieldError = e->GetFieldID(errorClass, "errorCode", "I");
        if (fieldError) {
            e->SetIntField(error, fieldError, static_cast<jint>(code));
        }
        jmethodID methodPut = e->GetMethodID(errorClass, "put", "(B)V");
        if (methodPut) {
            for (const char *c = msg; *c; ++c) {
                e->CallVoidMethod(error, methodPut, *c);
            }
        }
    }
}

JNIEXPORT jlong JNICALL
Java_io_questdb_jit_FiltersCompiler_compileFunction(JNIEnv *e,
                                                    jclass cl,
                                                    jlong filterAddress,
                                                    jlong filterSize,
                                                    jint options,
                                                    jobject error) {
    auto size = static_cast<size_t>(filterSize) / sizeof(instruction_t);

    if (filterAddress <= 0 || size <= 0) {
        fillJitErrorObject(e, error, ErrorCode::kErrorInvalidArgument, "Invalid argument passed");
        return 0;
    }

     CodeHolder code;
    code.init(gGlobalContext.rt.environment());
    FileLogger logger(stdout);
    bool debug = options & 1;
    if (debug) {
        logger.addFlags(FormatFlags::kRegCasts |
                        FormatFlags::kExplainImms);
        // logger.addFlags(FormatOptions::kFlagRegCasts |
        //                 FormatOptions::kFlagExplainImms |
        //                 FormatOptions::kFlagAnnotations);
        code.setLogger(&logger);
    }

    JitErrorHandler errorHandler;
    code.setErrorHandler(&errorHandler);

    a64::Compiler c(&code);
    // if (debug) {
        c.addDiagnosticOptions(DiagnosticOptions::kRAAnnotate);
    // }

    Function function(c);

    CompiledFn fn;

    std::cout << "---------------------------------------------------------------------" << std::endl;
    for (int i = 0; i < size; i++) {
        const instruction_t * instr = reinterpret_cast<const instruction_t *>(filterAddress);
        std::cout << "[" << opcode_to_string(instr[i].opcode) << "][" << instr[i].options << "] "
                  << "lo: " << instr[i].ipayload.lo << ", hi: " << instr[i].ipayload.hi
                  << ", dpayload: " << instr[i].dpayload << std::endl; 
    }

    function.begin_fn();
    function.compile(reinterpret_cast<const instruction_t *>(filterAddress), size, options);
    function.end_fn();

    Error err = errorHandler.error;

    if(err == ErrorCode::kErrorOk) {
        err = c.finalize();
    }

    if(err == ErrorCode::kErrorOk) {
        err = gGlobalContext.rt.add(&fn, &code);
    }

    fflush(logger.file());

    StringTmp<512> sb;
    Formatter::formatNodeList(sb, {}, &c);

    std::cout << "---------------------------------------------------------------------" << std::endl;
    std::cerr << sb.data() << '\n';
    std::cout << "---------------------------------------------------------------------" << std::endl;


    if(err != ErrorCode::kErrorOk) {
        fillJitErrorObject(e, error, err, errorHandler.message.data());
        return 0;
    }

    return reinterpret_cast<jlong>(fn);
}

JNIEXPORT void JNICALL
Java_io_questdb_jit_FiltersCompiler_freeFunction(JNIEnv *e, jclass cl, jlong fnAddress) {
    auto fn = reinterpret_cast<void *>(fnAddress);
    gGlobalContext.rt.release(fn);
}

JNIEXPORT jlong JNICALL Java_io_questdb_jit_FiltersCompiler_callFunction(JNIEnv *e,
                                                                         jclass cl,
                                                                         jlong fnAddress,
                                                                         jlong colsAddress,
                                                                         jlong colsSize,
                                                                         jlong varSizeIndexesAddress,
                                                                         jlong varsAddress,
                                                                         jlong varsSize,
                                                                         jlong rowsAddress,
                                                                         jlong rowsSize,
                                                                         jlong rowsStartOffset) {
    auto fn = reinterpret_cast<CompiledFn>(fnAddress);
    return fn(reinterpret_cast<int64_t *>(colsAddress),
              colsSize,
              reinterpret_cast<int64_t *>(varSizeIndexesAddress),
              reinterpret_cast<int64_t *>(varsAddress),
              varsSize,
              reinterpret_cast<int64_t *>(rowsAddress),
              rowsSize,
              rowsStartOffset);
}
