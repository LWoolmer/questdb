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

    // Reg
    // read_vars_mem(Compiler &c, data_type_t type, int32_t idx, const Gp &vars_ptr) {
    //     auto shift = type_shift(type);
    //     auto type_size = 1 << shift;
    //     return {Mem(vars_ptr, 8 * idx, type_size), type, data_kind_t::kMemory};
    // }

    // // Reads length of variable size column with header stored in data vector (string, binary).
    // Reg read_mem_varsize(Compiler &c,
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
    // Reg read_mem_varchar_header(Compiler &c,
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

    Reg read_mem(
            Compiler &c, data_type_t type, int32_t column_idx, const Gp &data_ptr,
            const Gp &varsize_aux_ptr, const Gp &input_index
    ) {
        if (type == data_type_t::varchar_header) {
            UNREACHABLE("Reading varchar_header from memory is not supported on A64");
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
            UNREACHABLE("Reading varsize from memory is not supported on A64");
            // return read_mem_varsize(c, header_size, column_idx, data_ptr, varsize_aux_ptr, input_index);
        }

        // Simple case: column has fixed-length data.

        Gp column_address = c.newGpx();
        c.ldr(column_address, ptr(data_ptr, 8 * column_idx));

        Mem row_data = ptr(column_address, input_index, Shift(ShiftOp::kLSL, type_shift(type)));
        Gp gp;
        Vec vec;
        switch (type) {
            case data_type_t::i8:
                gp = c.newInt32();
                c.ldrsb(gp, row_data);
                return gp;

            case data_type_t::i16:
                gp = c.newInt32();
                c.ldrsh(gp, row_data);
                return gp;

            case data_type_t::i32:
                gp = c.newInt32();
                c.ldr(gp, row_data);
                return gp;

            case data_type_t::i64:
                gp = c.newInt64();
                c.ldr(gp, row_data);
                return gp;

            case data_type_t::i128:
                vec = c.newVecQ();
                c.ldr(vec, row_data);
                return vec;

            case data_type_t::f32:
                vec = c.newVecS();
                c.ldr(vec, row_data);
                return vec;

            case data_type_t::f64:
                vec = c.newVecD();
                c.ldr(vec, row_data);
                return vec;

            default:
                UNREACHABLE("Unsupported data type for immediate: %s", data_type_to_string(type));
        }
        __builtin_unreachable();
    }

    Reg read_imm(Compiler &c, const instruction_t &instr) {
        auto type = static_cast<data_type_t>(instr.options);
        Gp gp;
        Gp gphi;
        Vec vec;
        switch (type) {
            case data_type_t::i8:
            case data_type_t::i16:
            case data_type_t::i32:
                gp = c.newInt32();
                c.mov(gp, instr.ipayload.lo);
                return gp;

            case data_type_t::i64:
                gp = c.newInt64();
                c.mov(gp, instr.ipayload.lo);
                return gp;

            case data_type_t::i128:
                vec = c.newVecQ();
                int64_t buffer[2];
                buffer[0] = instr.ipayload.lo;
                buffer[1] = instr.ipayload.hi;
                c.ldr(vec, c.newConst(ConstPoolScope::kLocal, buffer, 16));
                return vec;

            case data_type_t::f32:
                vec = c.newVecS();
                c.ldr(vec, c.newFloatConst(ConstPoolScope::kLocal, instr.dpayload));
                return vec;

            case data_type_t::f64:
                vec = c.newVecD();
                c.ldr(vec, c.newDoubleConst(ConstPoolScope::kLocal, instr.dpayload));
                return vec;

            default:
                UNREACHABLE("Unsupported data type for immediate: %s", data_type_to_string(type));
        }
        __builtin_unreachable();
    }

    std::pair<Reg, Reg> coerce_operands(Compiler &c, const Reg& lhs, const Reg &rhs, bool null_check) {
        RegType  lhs_type  = lhs.type();
        RegType  rhs_type  = rhs.type();
        RegGroup lhs_group = lhs.group();
        RegGroup rhs_group = rhs.group();

        switch(lhs_type) {
            case RegType::kGp32:
                switch (rhs_type) {
                    case RegType::kGp32:
                        return std::make_pair(lhs, rhs);
                    case RegType::kGp64:
                        return std::make_pair(to_int64(c, lhs.as<Gp>(), null_check), rhs);
                    case RegType::kVec32:
                        return std::make_pair(to_float(c, lhs.as<Gp>(), null_check), rhs);
                    case RegType::kVec64:
                        return std::make_pair(to_double(c, lhs.as<Gp>(), null_check), rhs);
                    default:
                        UNREACHABLE("Unsupported type coercion: %s --> %s", reg_type_to_string(lhs_type), reg_type_to_string(rhs_type));
                }
            case RegType::kGp64:
                switch (rhs_type) {
                    case RegType::kGp32:
                        return std::make_pair(lhs, to_int64(c, rhs.as<Gp>(), null_check));
                    case RegType::kGp64:
                        return std::make_pair(lhs, rhs);
                    case RegType::kVec32:
                        return std::make_pair(to_float(c, lhs.as<Gp>(), null_check), rhs);
                    case RegType::kVec64:
                        return std::make_pair(to_double(c, lhs.as<Gp>(), null_check), rhs);
                    default:
                        UNREACHABLE("Unsupported type coercion: %s --> %s", reg_type_to_string(lhs_type), reg_type_to_string(rhs_type));
                }
            case RegType::kVec32:
                switch (rhs_type) {
                    case RegType::kGp32:
                        return std::make_pair(lhs, to_float(c, rhs.as<Gp>(), null_check));
                    case RegType::kGp64:
                        return std::make_pair(lhs, to_float(c, rhs.as<Gp>(), null_check));
                    case RegType::kVec32:
                        return std::make_pair(lhs, rhs);
                    case RegType::kVec64:
                        return std::make_pair(to_double(c, lhs.as<Vec>(), null_check), rhs);
                    default:
                        UNREACHABLE("Unsupported type coercion: %s --> %s", reg_type_to_string(lhs_type), reg_type_to_string(rhs_type));
                }
            case RegType::kVec64:
                switch (rhs_type) {
                    case RegType::kGp32:
                        return std::make_pair(lhs, to_double(c, rhs.as<Gp>(), null_check));
                    case RegType::kGp64:
                        return std::make_pair(lhs, to_double(c, rhs.as<Gp>(), null_check));
                    case RegType::kVec32:
                        return std::make_pair(lhs, to_double(c, rhs.as<Vec>(), null_check));
                    case RegType::kVec64:
                        return std::make_pair(lhs, rhs);
                    default:
                        UNREACHABLE("Unsupported type coercion: %s --> %s", reg_type_to_string(lhs_type), reg_type_to_string(rhs_type));
                }
            case RegType::kVec128:
                return std::make_pair(lhs, rhs);
            default:
                UNREACHABLE("Unsupported type coercion: %s --> %s", reg_type_to_string(lhs_type), reg_type_to_string(rhs_type));
        }
        __builtin_unreachable();
    }

    inline Reg dispatch_bin_not(Compiler &c, const Reg &lhs) {
        if (lhs.isGp()) {
            return bin_not(c, lhs.as<Gp>());
        } else {
            UNREACHABLE("bin_not called with an argument that is not a Gp");
        }
    }

    inline Reg dispatch_bin_and(Compiler &c, const Reg &lhs, const Reg &rhs) {
        if (lhs.group() != rhs.group()) {
            UNREACHABLE("bin_and called with different groups");
        }

        if (lhs.isGp()) {
            return bin_and(c, lhs.as<Gp>(), rhs.as<Gp>());
        } else {
            UNREACHABLE("bin_and called with an argument that is not a Gp");
        }
    }

    inline Reg dispatch_bin_or(Compiler &c, const Reg &lhs, const Reg &rhs) {
        if (lhs.group() != rhs.group()) {
            UNREACHABLE("bin_or called with different groups");
        }

        if (lhs.isGp()) {
            return bin_or(c, lhs.as<Gp>(), rhs.as<Gp>());
        } else {
            UNREACHABLE("bin_or called with an argument that is not a Gp");
        }
    }

    inline Reg dispatch_neg(Compiler &c, const Reg &lhs, bool check_null) {
        if (lhs.isGp()) {
            return neg(c, lhs.as<Gp>(), check_null);
        } else if (lhs.isVec()) {
            return neg(c, lhs.as<Vec>(), check_null);
        } else {
            UNREACHABLE("neg called with an argument that is not a Gp or Vec");
        }
    }

    inline Reg dispatch_add(Compiler &c, const Reg &lhs, const Reg &rhs, bool check_null) {
        if (lhs.group() != rhs.group()) {
            UNREACHABLE("add called with different groups");
        }

        if (lhs.isGp() && rhs.isGp()) {
            return add(c, lhs.as<Gp>(), rhs.as<Gp>(), check_null);
        } else if  (lhs.isVec() && rhs.isVec()) {
            return add(c, lhs.as<Vec>(), rhs.as<Vec>(), check_null);
        } else {
            UNREACHABLE("add called with an argument that is not a Gp or Vec");
        }
    }

    inline Reg dispatch_sub(Compiler &c, const Reg &lhs, const Reg &rhs, bool check_null) {
        if (lhs.group() != rhs.group()) {
            UNREACHABLE("sub called with different groups");
        }

        if (lhs.isGp() && rhs.isGp()) {
            return sub(c, lhs.as<Gp>(), rhs.as<Gp>(), check_null);
        } else if  (lhs.isVec() && rhs.isVec()) {
            return sub(c, lhs.as<Vec>(), rhs.as<Vec>(), check_null);
        } else {
            UNREACHABLE("sub called with an argument that is not a Gp or Vec");
        }
    }

    inline Reg dispatch_cmp(Compiler &c, const Reg &lhs, const Reg &rhs, CondCode cond, bool check_null) {
        if (lhs.group() != rhs.group()) {
            UNREACHABLE("cmp called with different groups");
        }

        fprintf(stderr, "cmp called with lhs: %s, rhs: %s\n",
                reg_type_to_string(lhs.type()), reg_type_to_string(rhs.type()));

        if (lhs.isGp() && rhs.isGp()) {
            return cmp(c, lhs.as<Gp>(), rhs.as<Gp>(), cond, check_null);
        } else if (lhs.isVec() && rhs.isVec()) {
            return cmp(c, lhs.as<Vec>(), rhs.as<Vec>(), cond);
        } else {
            UNREACHABLE("cmp called with an argument that is not a Gp or Vec");
        }
    }

    // Reg mul(Compiler &c, const Reg &lhs, const Reg &rhs, bool null_check) {
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

    // Reg div(Compiler &c, const Reg &lhs, const Reg &rhs, bool null_check) {
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

    void emit_bin_op(Compiler &c, const instruction_t &instr, ZoneStack<Reg> &values, bool null_check) {
        std::pair<Reg, Reg> coerced = coerce_operands(c, values.pop(), values.pop(), null_check);
        Reg lhs = coerced.first;
        Reg rhs = coerced.second;

        char buf_instr[512];
        char buf_lhs[512];
        char buf_rhs[512];
        char buf_ret[512];
        reg_to_string(c, lhs, buf_lhs);
        reg_to_string(c, rhs, buf_rhs);
        instr.to_string(buf_instr);
        comment(c, "                                    >>> ; %s %s %s", buf_lhs, buf_instr, buf_rhs);

        Reg ret;

        switch (instr.opcode) {
            case opcodes::And:
                ret = dispatch_bin_and(c, lhs, rhs);
                break;
            case opcodes::Or:
                ret = dispatch_bin_or(c, lhs, rhs);
                break;
            case opcodes::Eq:
                ret = dispatch_cmp(c, lhs, rhs, CondCode::kEQ, null_check);
                break;
            case opcodes::Ne:
                ret = dispatch_cmp(c, lhs, rhs, CondCode::kNE, null_check);
                break;
            case opcodes::Gt:
                ret = dispatch_cmp(c, lhs, rhs, CondCode::kGT, null_check);
                break;
            case opcodes::Ge:
                ret = dispatch_cmp(c, lhs, rhs, CondCode::kGE, null_check);
                break;
            case opcodes::Lt:
                ret = dispatch_cmp(c, lhs, rhs, CondCode::kLT, null_check);
                break;
            case opcodes::Le:
                ret = dispatch_cmp(c, lhs, rhs, CondCode::kLE, null_check);
                break;
            case opcodes::Add:
                ret = dispatch_add(c, lhs, rhs, null_check);
                break;
            case opcodes::Sub:
                ret = dispatch_sub(c, lhs, rhs, null_check);
                break;
            // case opcodes::Mul:
            //     ret = mul(c, lhs, rhs, null_check);
            //     break;
            // case opcodes::Div:
            //     ret = div(c, lhs, rhs, null_check);
            //     break;
            default:
                UNREACHABLE("Unsupported operation: %s", opcode_to_string(instr.opcode));
        }

        reg_to_string(c, ret, buf_ret);
        comment(c, "                                    <<< ; %s", buf_ret);
        values.append(ret);
    }

    void
    emit_code(Compiler &c, const instruction_t *istream, size_t size, ZoneStack<Reg> &values,
              bool null_check,
              const Gp &data_ptr,
              const Gp &varsize_aux_ptr,
              const Gp &vars_ptr,
              const Gp &input_index) {

        for (size_t i = 0; i < size; ++i) {
            auto &instr = istream[i];
            auto type = static_cast<data_type_t>(instr.options);
            auto column_idx  = static_cast<int32_t>(instr.ipayload.lo);

            char buf_ret[512];
            Reg ret;

            fprintf(stderr, "Processing instruction %zu: %s\n", i, opcode_to_string(instr.opcode));

            switch (instr.opcode) {
                case opcodes::Inv:
                    return; // todo: throw exception

                case opcodes::Ret:
                    return;

                case opcodes::Var:
                    UNREACHABLE("Unsupported operation: %s", opcode_to_string(instr.opcode));
                    // auto type = static_cast<data_type_t>(instr.options);
                    // auto idx  = static_cast<int32_t>(instr.ipayload.lo);
                    // values.append(read_vars_mem(c, type, idx, vars_ptr));
                    // break;

                case opcodes::Mem:
                    instr.to_string(buf_ret);
                    comment(c, "                                    >>> ; %s", buf_ret);
                    ret = read_mem(c, type, column_idx, data_ptr, varsize_aux_ptr, input_index);
                    reg_to_string(c, ret, buf_ret);
                    comment(c, "                                    <<< ; %s", buf_ret);
                    values.append(ret);
                    break;

                case opcodes::Imm:
                    instr.to_string(buf_ret);
                    comment(c, "                                    >>> ; %s", buf_ret);
                    ret = read_imm(c, instr);
                    reg_to_string(c, ret, buf_ret);
                    comment(c, "                                    <<< ; %s", buf_ret);
                    values.append(ret);
                    break;

                case opcodes::Neg:
                    values.append(dispatch_neg(c, values.pop(), null_check));
                    break;

                case opcodes::Not:
                    UNREACHABLE("Unsupported operation: %s", opcode_to_string(instr.opcode));
                    // values.append(bin_not(c, values.pop()));
                    // break;

                default:
                    emit_bin_op(c, instr, values, null_check);
                    break;
            }
        }
    }
}

#endif //QUESTDB_JIT_A64_H
