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
        Label l_loop = c.newLabel();
        Label l_exit = c.newLabel();

        c.cmp(input_index, stop);
        c.b_gt(l_exit);

        c.bind(l_loop);

        for (int i = 0; i < unroll_factor; ++i) {
            // questdb::a64::emit_code(c, istream, size, values, null_check, data_ptr, varsize_aux_ptr, vars_ptr, input_index);

            auto mask = values.pop();

            // a64::Gp adjusted_id = c.newInt64("input_index_+_rows_id_start_offset");
            // c.add(adjusted_id, ptr(input_index, rows_id_start_offset)); // input_index + rows_id_start_offset
            // c.mov(qword_ptr(rows_ptr, output_index, 3), adjusted_id);

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
            c.mov(stop, rows_size);
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

    std::cerr << "test\n";
    auto size = static_cast<size_t>(filterSize) / sizeof(instruction_t);

    std::cerr << "test2\n";
    if (filterAddress <= 0 || size <= 0) {
        fillJitErrorObject(e, error, ErrorCode::kErrorInvalidArgument, "Invalid argument passed");
        return 0;
    }

    CodeHolder code;
    code.init(gGlobalContext.rt.environment());

    std::cerr << "test3\n";
    FILE *file = fopen("/root/.questdb/log/asmjit_last.log", "w");

    std::cerr << "test4\n";
    FileLogger logger(file);

    std::cerr << "test5\n";
    bool debug = options & 1;
    debug = true;
    if (debug) {

        std::cerr << "test6\n";
        logger.addFlags(FormatFlags::kRegCasts |
                        FormatFlags::kExplainImms);
        // logger.addFlags(FormatOptions::kFlagRegCasts |
        //                 FormatOptions::kFlagExplainImms |
        //                 FormatOptions::kFlagAnnotations);

        std::cerr << "test7\n";
        code.setLogger(&logger);

        std::cerr << "test8\n";
    }

    JitErrorHandler errorHandler;
    std::cerr << "test9\n";
    code.setErrorHandler(&errorHandler);
    std::cerr << "test10\n";

    a64::Compiler c(&code);
    std::cerr << "test11\n";
    if (debug) {
        std::cerr << "test12\n";
        c.addDiagnosticOptions(DiagnosticOptions::kRAAnnotate);
    }

    std::cerr << "test13\n";
    Function function(c);

    std::cerr << "test14\n";
    CompiledFn fn;

    std::cerr << "test15\n";
    function.begin_fn();
    std::cerr << "test16\n";
    function.compile(reinterpret_cast<const instruction_t *>(filterAddress), size, options);
    std::cerr << "test17\n";
    function.end_fn();

    std::cerr << "test18\n";
    Error err = errorHandler.error;
    std::cerr << "test19\n";
    if(err == ErrorCode::kErrorOk) {
        std::cerr << "test20\n";
        err = c.finalize();
    }
    
    std::cerr << "test21\n";
    if(err == ErrorCode::kErrorOk) {
        std::cerr << "test22\n";
        err = gGlobalContext.rt.add(&fn, &code);
    }

   std::cerr << "test23\n";
    StringTmp<512> sb;
       std::cerr << "test24\n";
    Formatter::formatNodeList(sb, {}, &c);
       std::cerr << "test25\n";
    logger.log(sb);

   std::cerr << "test26\n";
    std::cerr << sb.data() << '\n';

    std::cerr << "test27\n";

   std::cerr << "test30\n";
    if(err != ErrorCode::kErrorOk) {
       std::cerr << "test31\n";
        fillJitErrorObject(e, error, err, errorHandler.message.data());
           std::cerr << "test32\n";
        return 0;
    }

   std::cerr << "test33\n";
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
