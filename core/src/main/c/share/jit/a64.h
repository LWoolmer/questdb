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

#ifndef QUESTDB_JIT_A64_H
#define QUESTDB_JIT_A64_H

#include "common.h"
#include "impl/a64.h"
#include <stdexcept>

namespace questdb::a64 {
    using namespace asmjit;
    using namespace asmjit::a64;
    using namespace asmjit::arm;

    // jit_value_t
    // read_vars_mem(Compiler &c, data_type_t type, int32_t idx, const Gp &vars_ptr) {
    //     auto shift = type_shift(type);
    //     auto type_size = 1 << shift;
    //     return {Mem(vars_ptr, 8 * idx, type_size), type, data_kind_t::kMemory};
    // }

    // // Reads length of variable size column with header stored in data vector (string, binary).
    // jit_value_t read_mem_varsize(Compiler &c,
    //                              uint32_t header_size,
    //                              int32_t column_idx,
    //                              const Gp &data_ptr,
    //                              const Gp &varsize_aux_ptr,
    //                              const Gp &input_index) {
    //     // Column has variable-size data with header stored in data vector.
    //     // First, we load this and the next data vector offsets from the aux vector.
    //     // When the offset difference is zero, it can indicate an empty  value (length 0)
    //     // or a NULL (length -1). In the zero difference case, we have to load the header
    //     // from the data vector. In the positive difference case, the difference is equal
    //     // to the length, so there is no need to do an extra load.
    //     Label l_nonzero = c.newLabel();
    //     Gp offset = c.newInt64("offset");
    //     Gp length = c.newInt64("length");
    //     Gp varsize_aux_address = c.newInt64("varsize_aux_address");
    //     Gp next_input_index = c.newInt64("next_input_index");
    //     c.mov(next_input_index, input_index);
    //     c.inc(next_input_index);
    //     c.mov(varsize_aux_address, ptr(varsize_aux_ptr, 8 * column_idx, 8));
    //     auto offset_shift = type_shift(data_type_t::i64);
    //     auto offset_size = 1 << offset_shift;
    //     c.mov(offset, ptr(varsize_aux_address, input_index, offset_shift, 0, offset_size));
    //     c.mov(length, ptr(varsize_aux_address, next_input_index, offset_shift, 0, offset_size));
    //     c.sub(length, offset);
    //     c.sub(length, header_size);
    //     // length now contains the length of the value. It can be zero for two reasons:
    //     // empty value or NULL value.
    //     c.jnz(l_nonzero);
    //     // If it's zero, we have to load the actual header value, which can be 0 or -1.
    //     Gp column_address = c.newInt64("column_address");
    //     c.mov(column_address, ptr(data_ptr, 8 * column_idx, 8));
    //     c.mov(length, ptr(column_address, offset, 0, 0, header_size));
    //     c.bind(l_nonzero);
    //     if (header_size == 4) {
    //         return {length.r32(), data_type_t::i32, data_kind_t::kMemory};
    //     }
    //     return {length, data_type_t::i64, data_kind_t::kMemory};
    // }

    // // Reads length part of the varchar header for aux vector.
    // // This part is stored in the lowest bytes of the header
    // // (see VarcharTypeDriver to understand the format).
    // //
    // // Note: unlike read_mem_varsize this method doesn't return the length,
    // //       so it can only be used in NULL checks.
    // jit_value_t read_mem_varchar_header(Compiler &c,
    //                                     int32_t column_idx,
    //                                     const Gp &varsize_aux_ptr,
    //                                     const Gp &input_index) {
    //     Gp varsize_aux_address = c.newInt64("varsize_aux_address");
    //     c.mov(varsize_aux_address, ptr(varsize_aux_ptr, 8 * column_idx, 8));

    //     Gp header_offset = c.newInt64("header_offset");
    //     c.mov(header_offset, input_index);
    //     auto header_shift = type_shift(data_type_t::i128);
    //     c.sal(header_offset, header_shift);

    //     Gp header = c.newInt64("header");
    //     c.mov(header, ptr(varsize_aux_address, header_offset, 0));

    //     return {header, data_type_t::i64, data_kind_t::kMemory};
    // }

    jit_value_t read_mem(
            Compiler &c, data_type_t type, int32_t column_idx, const Gp &data_ptr,
            const Gp &varsize_aux_ptr, const Gp &input_index
    ) {
        if (type == data_type_t::varchar_header) {
            throw std::runtime_error(
                "Reading varchar_header from memory is not supported on A64"
            );
            __builtin_unreachable();
            // return read_mem_varchar_header(c, column_idx, varsize_aux_ptr, input_index);
        }

        uint32_t header_size;
        switch (type) {
            case data_type_t::string_header:
                header_size = 4;
                break;
            case data_type_t::binary_header:
                header_size = 8;
                break;
            default:
                header_size = 0;
        }
        if (header_size != 0) {
            throw std::runtime_error(
                "Reading varsize from memory is not supported on A64"
            );
            __builtin_unreachable();
            // return read_mem_varsize(c, header_size, column_idx, data_ptr, varsize_aux_ptr, input_index);
        }

        // Simple case: column has fixed-length data.

        Gp column_address = c.newGpx();
        c.ldr(column_address, ptr(data_ptr, 8 * column_idx));

        Mem row_data = ptr(column_address, input_index, Shift(ShiftOp::kLSL, type_shift(type)));
        return {row_data, type, data_kind_t::kMemory};
    }

    jit_value_t mem2reg(Compiler &c, const jit_value_t &v) {
        auto type = v.dtype();
        auto mem = v.op().as<Mem>();
        Gp gp;
        Vec vec;
        switch (type) {
            case data_type_t::i8:
                gp = c.newInt32();
                c.ldrsb(gp, mem);
                return {gp, type, data_kind_t::kMemory};

            case data_type_t::i16:
                gp = c.newInt32();
                c.ldrsh(gp, mem);
                return {gp, type, data_kind_t::kMemory};

            case data_type_t::i32:
                gp = c.newInt32();
                c.ldr(gp, mem);
                return {gp, type, data_kind_t::kMemory};

            case data_type_t::i64:
                gp = c.newInt64();
                c.ldr(gp, mem);
                return {gp, type, data_kind_t::kMemory};

            case data_type_t::i128:
                vec = c.newVecQ();
                c.ldr(vec, mem);
                return {vec, type, data_kind_t::kMemory};

            case data_type_t::f32:
                vec = c.newVecS();
                c.ldr(vec, mem);
                return {vec, type, data_kind_t::kMemory};

            case data_type_t::f64:
                vec = c.newVecD();
                c.ldr(vec, mem);
                return {vec, type, data_kind_t::kMemory};

            default:
                throw std::runtime_error(
                    "Unsupported data type for mem2reg: " + std::string(data_type_to_string(type))
                );
                __builtin_unreachable();
        }
    }

    jit_value_t read_imm(Compiler &c, const instruction_t &instr) {
        auto type = static_cast<data_type_t>(instr.options);
        switch (type) {
            case data_type_t::i8:
            case data_type_t::i16:
            case data_type_t::i32:
            case data_type_t::i64:
                return {imm(instr.ipayload.lo), type, data_kind_t::kConst};

            // REVISIT 128bit imm
            // case data_type_t::i128: {
            //     return {
            //         c.newConst(ConstPoolScope::kLocal, &instr.ipayload, 16),
            //         type,
            //         data_kind_t::kMemory
            //     };
            // }

            case data_type_t::f32:
            case data_type_t::f64:
                return {imm(instr.dpayload), type, data_kind_t::kConst};

            default:
                throw std::runtime_error(
                    "Unsupported data type for immediate: " + std::string(data_type_to_string(type))
                );
                __builtin_unreachable();
        }
    }

    jit_value_t imm2reg(Compiler &c, data_type_t dst_type, const jit_value_t &v) {
        Imm imm = v.imm();
        Vec vec;
        Gp gp;
        Mem mem;

        if (imm.isInt()) {
            auto value = imm.valueAs<int64_t>();
            switch (dst_type) {
                case data_type_t::f32:
                    vec = c.newVecS();
                    mem = c.newFloatConst(ConstPoolScope::kLocal, static_cast<float>(value));
                    c.ldr(vec, mem);
                    return {vec, dst_type, data_kind_t::kConst};

                case data_type_t::f64:
                    vec = c.newVecD();
                    mem = c.newDoubleConst(ConstPoolScope::kLocal, static_cast<double>(value));
                    c.ldr(vec, mem);
                    return {vec, dst_type, data_kind_t::kConst};

                case data_type_t::i128:
                    vec = c.newVecQ();
                    c.movi(vec, value);
                    c.movi(vec.at(1), 0);
                    return {vec, dst_type, data_kind_t::kConst};

                case data_type_t::i64:
                    gp = c.newInt64();
                    c.mov(gp, value);
                    return {gp, dst_type, data_kind_t::kConst};

                default:
                    gp = c.newInt32();
                    c.mov(gp, static_cast<int32_t>(value));
                    return {gp, dst_type, data_kind_t::kConst};
            }
        } else {
            auto value = imm.valueAs<double>();
            switch (dst_type) {
                case data_type_t::f32:
                    vec = c.newVecS();
                    mem = c.newFloatConst(ConstPoolScope::kLocal, static_cast<float>(value));
                    c.ldr(vec, mem);
                    return {vec, dst_type, data_kind_t::kConst};

                case data_type_t::f64:
                    vec = c.newVecD();
                    mem = c.newDoubleConst(ConstPoolScope::kLocal, static_cast<double>(value));
                    c.ldr(vec, mem);
                    return {vec, dst_type, data_kind_t::kConst};

                case data_type_t::i128:
                    vec = c.newVecQ();
                    c.movi(vec, static_cast<int64_t>(value));
                    c.movi(vec.at(1), 0);
                    return {vec, dst_type, data_kind_t::kConst};

                case data_type_t::i64:
                    gp = c.newInt64();
                    c.mov(gp, value);
                    return {gp, dst_type, data_kind_t::kConst};

                default:
                    gp = c.newInt32();
                    c.mov(gp, static_cast<int32_t>(value));
                    return {gp, dst_type, data_kind_t::kConst};
            }
        }

        __builtin_unreachable();
    }

    jit_value_t load_register(Compiler &c, data_type_t dst_type, const jit_value_t &v) {
        if (v.op().isImm()) {
            return imm2reg(c, dst_type, v);
        } else if (v.op().isMem()) {
            return mem2reg(c, v);
        } else {
            return v;
        }
    }

    jit_value_t load_register(Compiler &c, const jit_value_t &v) {
        return load_register(c, v.dtype(), v);
    }

    std::pair<jit_value_t, jit_value_t> load_registers(Compiler &c, const jit_value_t &lhs, const jit_value_t &rhs) {
        data_type_t lt;
        data_type_t rt;
        if (lhs.op().isImm() && !rhs.op().isImm()) {
            lt = rhs.dtype();
            rt = rhs.dtype();
        } else if (rhs.op().isImm() && !lhs.op().isImm()) {
            lt = lhs.dtype();
            rt = lhs.dtype();
        } else {
            lt = lhs.dtype();
            rt = rhs.dtype();
        }
        jit_value_t l = load_register(c, lt, lhs);
        jit_value_t r = load_register(c, rt, rhs);
        return {l, r};
    }

    jit_value_t cmp(Compiler &c, const jit_value_t &lhs, const jit_value_t &rhs, CondCode cond, bool null_check) {
        auto dt = lhs.dtype();
        auto dk = dst_kind(lhs, rhs);
        switch (dt) {
            case data_type_t::i8:
            case data_type_t::i16:
            case data_type_t::i32:
            case data_type_t::string_header:
            case data_type_t::i64:
            case data_type_t::binary_header:
            case data_type_t::varchar_header:
                return {cmp(c, lhs.gp(), rhs.gp(), cond, null_check), data_type_t::i32, dst_kind(lhs, rhs)};
            case data_type_t::i128:
            case data_type_t::f32:
            case data_type_t::f64:
                return {cmp(c, lhs.vec(), rhs.vec(), cond), data_type_t::i32, dst_kind(lhs, rhs)};
            default:
                throw std::runtime_error(
                    "Unsupported data type for lhs of cmp: " + std::string(data_type_to_string(dt))
                );
                __builtin_unreachable();
        }
    }

    // jit_value_t neg(Compiler &c, const jit_value_t &lhs, bool null_check) {
    //     auto dt = lhs.dtype();
    //     auto dk = lhs.dkind();
    //     switch (dt) {
    //         case data_type_t::i8:
    //         case data_type_t::i16:
    //         case data_type_t::i32:
    //             return {int32_neg(c, lhs.gp().r32(), null_check), dt, dk};
    //         case data_type_t::i64:
    //             return {int64_neg(c, lhs.gp(), null_check), dt, dk};
    //         case data_type_t::f32:
    //             return {float_neg(c, lhs.xmm()), dt, dk};
    //         case data_type_t::f64:
    //             return {double_neg(c, lhs.xmm()), dt, dk};
    //         default:
    //             __builtin_unreachable();
    //     }
    // }

    jit_value_t bin_not(Compiler &c, const jit_value_t &lhs) {
        auto dt = lhs.dtype();
        auto dk = lhs.dkind();
        Gp lhs_gp = lhs.gp().r32();
        c.eor(lhs_gp, lhs_gp, -1);
        return {lhs_gp, dt, dk};
    }

    jit_value_t bin_and(Compiler &c, const jit_value_t &lhs, const jit_value_t &rhs) {
        auto dt = lhs.dtype();
        auto dk = dst_kind(lhs, rhs);
        Gp lhs_gp = lhs.gp().r32();
        Gp rhs_gp = rhs.gp().r32();
        c.and_(lhs_gp, lhs_gp, rhs_gp);
        return {lhs_gp, dt, dk};
    }

    jit_value_t bin_or(Compiler &c, const jit_value_t &lhs, const jit_value_t &rhs) {
        auto dt = lhs.dtype();
        auto dk = dst_kind(lhs, rhs);
        Gp lhs_gp = lhs.gp().r32();
        Gp rhs_gp = rhs.gp().r32();
        c.orr(lhs_gp, lhs_gp, rhs_gp);
        return {lhs_gp, dt, dk};
    }

    jit_value_t add(Compiler &c, const jit_value_t &lhs, const jit_value_t &rhs, bool null_check) {
        auto dt = lhs.dtype();
        auto dk = dst_kind(lhs, rhs);
        switch (dt) {
            case data_type_t::i8:
            case data_type_t::i16:
            case data_type_t::i32:
            case data_type_t::i64:
                return {add(c, lhs.gp(), rhs.gp(), null_check), dt, dk};
            case data_type_t::f32:
            case data_type_t::f64:
                return {add(c, lhs.vec(), rhs.vec(), null_check), dt, dk};
            default:
                __builtin_unreachable();
        }
    }

    jit_value_t sub(Compiler &c, const jit_value_t &lhs, const jit_value_t &rhs, bool null_check) {
        auto dt = lhs.dtype();
        auto dk = dst_kind(lhs, rhs);
        switch (dt) {
            case data_type_t::i8:
            case data_type_t::i16:
            case data_type_t::i32:
            case data_type_t::i64:
                return {sub(c, lhs.gp(), rhs.gp(), null_check), dt, dk};
            case data_type_t::f32:
            case data_type_t::f64:
                return {sub(c, lhs.vec(), rhs.vec(), null_check), dt, dk};
            default:
                __builtin_unreachable();
        }
    }

    // jit_value_t mul(Compiler &c, const jit_value_t &lhs, const jit_value_t &rhs, bool null_check) {
    //     auto dt = lhs.dtype();
    //     auto dk = dst_kind(lhs, rhs);
    //     switch (dt) {
    //         case data_type_t::i8:
    //         case data_type_t::i16:
    //         case data_type_t::i32:
    //             return {int32_mul(c, lhs.gp().r32(), rhs.gp().r32(), null_check), dt, dk};
    //         case data_type_t::i64:
    //             return {int64_mul(c, lhs.gp(), rhs.gp(), null_check), dt, dk};
    //         case data_type_t::f32:
    //             return {float_mul(c, lhs.xmm(), rhs.xmm()), dt, dk};
    //         case data_type_t::f64:
    //             return {double_mul(c, lhs.xmm(), rhs.xmm()), dt, dk};
    //         default:
    //             __builtin_unreachable();
    //     }
    // }

    // jit_value_t div(Compiler &c, const jit_value_t &lhs, const jit_value_t &rhs, bool null_check) {
    //     auto dt = lhs.dtype();
    //     auto dk = dst_kind(lhs, rhs);
    //     switch (dt) {
    //         case data_type_t::i8:
    //         case data_type_t::i16:
    //         case data_type_t::i32:
    //             return {int32_div(c, lhs.gp().r32(), rhs.gp().r32(), null_check), dt, dk};
    //         case data_type_t::i64:
    //             return {int64_div(c, lhs.gp(), rhs.gp(), null_check), dt, dk};
    //         case data_type_t::f32:
    //             return {float_div(c, lhs.xmm(), rhs.xmm()), dt, dk};
    //         case data_type_t::f64:
    //             return {double_div(c, lhs.xmm(), rhs.xmm()), dt, dk};
    //         default:
    //             __builtin_unreachable();
    //     }
    // }

    inline bool cvt_null_check(data_type_t type) {
        return !(type == data_type_t::i8 || type == data_type_t::i16);
    }

    inline std::pair<jit_value_t, jit_value_t>
    convert(Compiler &c, const jit_value_t &lhs, const jit_value_t &rhs, bool null_check) {
        switch (lhs.dtype()) {
            case data_type_t::i8:
            case data_type_t::i16:
            case data_type_t::i32:
                switch (rhs.dtype()) {
                    case data_type_t::i8:
                    case data_type_t::i16:
                    case data_type_t::i32:
                        return std::make_pair(lhs, rhs);
                    case data_type_t::i64:
                        return std::make_pair(
                                jit_value_t(
                                        to_int64(c, lhs.gp().r32(), null_check && cvt_null_check(lhs.dtype())),
                                        data_type_t::i64,
                                        lhs.dkind()), rhs);
                    case data_type_t::f32:
                        return std::make_pair(
                                jit_value_t(
                                        to_float(c, lhs.gp().r32(), null_check && cvt_null_check(lhs.dtype())),
                                        data_type_t::f32,
                                        lhs.dkind()), rhs);
                    case data_type_t::f64:
                        return std::make_pair(
                                jit_value_t(
                                        to_double(c, lhs.gp().r32(), null_check && cvt_null_check(lhs.dtype())),
                                        data_type_t::f64,
                                        lhs.dkind()), rhs);
                    default:
                        throw std::runtime_error(
                            "Unsupported data types for convert: " + std::string(data_type_to_string(lhs.dtype())) + " " + std::string(data_type_to_string(rhs.dtype()))
                        );
                        __builtin_unreachable();
                }
                break;
            case data_type_t::i64:
                switch (rhs.dtype()) {
                    case data_type_t::i8:
                    case data_type_t::i16:
                    case data_type_t::i32:
                        return std::make_pair(lhs,
                                              jit_value_t(
                                                      to_int64(c, rhs.gp().r32(),
                                                                     null_check && cvt_null_check(rhs.dtype())),
                                                      data_type_t::i64, rhs.dkind()));
                    case data_type_t::i64:
                        return std::make_pair(lhs, rhs);
                    case data_type_t::f32:
                        return std::make_pair(
                                jit_value_t(to_double(c, lhs.gp().r64(), null_check), data_type_t::f64,
                                            lhs.dkind()),
                                jit_value_t(to_double(c, rhs.vec().v32()), data_type_t::f64, rhs.dkind()));
                    case data_type_t::f64:
                        return std::make_pair(
                                jit_value_t(to_double(c, lhs.gp().r64(), null_check), data_type_t::f64, lhs.dkind()),
                                rhs);
                    default:
                        throw std::runtime_error(
                            "Unsupported data types for convert: " + std::string(data_type_to_string(lhs.dtype())) + " " + std::string(data_type_to_string(rhs.dtype()))
                        );
                        __builtin_unreachable();
                }
                break;
            case data_type_t::f32:
                switch (rhs.dtype()) {
                    case data_type_t::i8:
                    case data_type_t::i16:
                    case data_type_t::i32:
                        return std::make_pair(lhs,
                                              jit_value_t(
                                                      to_float(c, rhs.gp().r32(),
                                                                     null_check && cvt_null_check(rhs.dtype())),
                                                      data_type_t::f32, rhs.dkind()));
                    case data_type_t::i64:
                        return std::make_pair(jit_value_t(to_double(c, lhs.vec().v32()), data_type_t::f64, lhs.dkind()),
                                              jit_value_t(to_double(c, rhs.gp(), null_check), data_type_t::f64,
                                                          rhs.dkind()));
                    case data_type_t::f32:
                        return std::make_pair(lhs, rhs);
                    case data_type_t::f64:
                        return std::make_pair(jit_value_t(to_double(c, lhs.vec().v32()), data_type_t::f64, lhs.dkind()),
                                              rhs);
                    default:
                        throw std::runtime_error(
                            "Unsupported data types for convert: " + std::string(data_type_to_string(lhs.dtype())) + " " + std::string(data_type_to_string(rhs.dtype()))
                        );
                        __builtin_unreachable();
                }
                break;
            case data_type_t::f64:
                switch (rhs.dtype()) {
                    case data_type_t::i8:
                    case data_type_t::i16:
                    case data_type_t::i32:
                        return std::make_pair(lhs,
                                              jit_value_t(
                                                      to_double(c, rhs.gp().r32(),
                                                                      null_check && cvt_null_check(rhs.dtype())),
                                                      data_type_t::f64, rhs.dkind()));
                    case data_type_t::i64:
                        return std::make_pair(lhs,
                                              jit_value_t(to_double(c, rhs.gp().r64(), null_check), data_type_t::f64,
                                                          rhs.dkind()));
                    case data_type_t::f32:
                        return std::make_pair(lhs,
                                              jit_value_t(to_double(c, rhs.vec().v32()),
                                                          data_type_t::f64,
                                                          rhs.dkind()));
                    case data_type_t::f64:
                        return std::make_pair(lhs, rhs);
                    default:
                        throw std::runtime_error(
                            "Unsupported data types for convert: " + std::string(data_type_to_string(lhs.dtype())) + " " + std::string(data_type_to_string(rhs.dtype()))
                        );
                        __builtin_unreachable();
                }
                break;
            case data_type_t::i128:
            case data_type_t::string_header:
            case data_type_t::binary_header:
            case data_type_t::varchar_header:
                return std::make_pair(lhs, rhs);
            default:
                throw std::runtime_error(
                    "Unsupported data types for convert: " + std::string(data_type_to_string(lhs.dtype())) + " " + std::string(data_type_to_string(rhs.dtype()))
                );
                __builtin_unreachable();
        }
    }

    inline jit_value_t get_argument(Compiler &c, ZoneStack<jit_value_t> &values) {
        auto arg = values.pop();
        return load_register(c, arg);
    }

    inline std::pair<jit_value_t, jit_value_t>
    get_arguments(Compiler &c, ZoneStack<jit_value_t> &values, bool null_check) {
        auto lhs = values.pop();
        auto rhs = values.pop();
        auto args = load_registers(c, lhs, rhs);
        return convert(c, args.first, args.second, null_check);
    }

    void emit_bin_op(Compiler &c, const instruction_t &instr, ZoneStack<jit_value_t> &values, bool null_check) {
        auto args = get_arguments(c, values, null_check);
        jit_value_t lhs = args.first;
        jit_value_t rhs = args.second;

        char buf_instr[512];
        char buf_lhs[512];
        char buf_rhs[512];
        char buf_ret[512];
        lhs.to_string(c, buf_lhs);
        rhs.to_string(c, buf_rhs);
        instr.to_string(buf_instr);
        comment(c, "                                    >>> ; %s %s %s", buf_lhs, buf_instr, buf_rhs);

        jit_value_t ret;

        switch (instr.opcode) {
            case opcodes::And:
                ret = bin_and(c, lhs, rhs);
                break;
            case opcodes::Or:
                ret = bin_or(c, lhs, rhs);
                break;
            case opcodes::Eq:
                ret = cmp(c, lhs, rhs, CondCode::kEQ, null_check);
                break;
            case opcodes::Ne:
                ret = cmp(c, lhs, rhs, CondCode::kNE, null_check);
                break;
            case opcodes::Gt:
                ret = cmp(c, lhs, rhs, CondCode::kGT, null_check);
                break;
            case opcodes::Ge:
                ret = cmp(c, lhs, rhs, CondCode::kGE, null_check);
                break;
            case opcodes::Lt:
                ret = cmp(c, lhs, rhs, CondCode::kLT, null_check);
                break;
            case opcodes::Le:
                ret = cmp(c, lhs, rhs, CondCode::kLE, null_check);
                break;
            case opcodes::Add:
                ret = add(c, lhs, rhs, null_check);
                break;
            case opcodes::Sub:
                ret = sub(c, lhs, rhs, null_check);
                break;
            // case opcodes::Mul:
            //     ret = mul(c, lhs, rhs, null_check);
            //     break;
            // case opcodes::Div:
            //     ret = div(c, lhs, rhs, null_check);
            //     break;
            default:
                throw std::runtime_error(
                    "Unsupported operation: " + std::string(opcode_to_string(instr.opcode))
                );
                __builtin_unreachable();
        }

        ret.to_string(c, buf_ret);
        comment(c, "                                    <<< ; %s", buf_ret);
        values.append(ret);
    }

    void
    emit_code(Compiler &c, const instruction_t *istream, size_t size, ZoneStack<jit_value_t> &values,
              bool null_check,
              const Gp &data_ptr,
              const Gp &varsize_aux_ptr,
              const Gp &vars_ptr,
              const Gp &input_index) {

        for (size_t i = 0; i < size; ++i) {
            auto &instr = istream[i];
            auto type = static_cast<data_type_t>(instr.options);
            auto column_idx  = static_cast<int32_t>(instr.ipayload.lo);

            char buf[512];
            auto err = std::runtime_error("Unsupported operation: " + std::string(opcode_to_string(instr.opcode)));
            jit_value_t ret;

            switch (instr.opcode) {
                case opcodes::Inv:
                    return; // todo: throw exception

                case opcodes::Ret:
                    return;

                case opcodes::Var:
                    throw err;
                    // auto type = static_cast<data_type_t>(instr.options);
                    // auto idx  = static_cast<int32_t>(instr.ipayload.lo);
                    // values.append(read_vars_mem(c, type, idx, vars_ptr));
                    // break;

                case opcodes::Mem:
                    instr.to_string(buf);
                    comment(c, "                                    >>> ; %s", buf);
                    ret = read_mem(c, type, column_idx, data_ptr, varsize_aux_ptr, input_index);
                    ret.to_string(c, buf);
                    comment(c, "                                    <<< ; %s", buf);
                    values.append(ret);
                    break;

                case opcodes::Imm:
                    instr.to_string(buf);
                    comment(c, "                                    >>> ; %s", buf);
                    ret = read_imm(c, instr);
                    ret.to_string(c, buf);
                    comment(c, "                                    <<< ; %s", buf);
                    values.append(ret);
                    break;

                case opcodes::Neg:
                    throw err;
                    // values.append(neg(c, get_argument(c, values), null_check));
                    // break;

                case opcodes::Not:
                    throw err;
                    // values.append(bin_not(c, get_argument(c, values)));
                    // break;

                default:
                    emit_bin_op(c, instr, values, null_check);
                    break;
            }
        }
    }
}

#endif //QUESTDB_JIT_A64_H
