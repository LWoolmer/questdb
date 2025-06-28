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

#ifndef QUESTDB_JIT_COMMON_H
#define QUESTDB_JIT_COMMON_H

#include <asmjit/asmjit.h>
#include <asmjit/a64.h>

enum class data_type_t : uint8_t {
    i8,
    i16,
    i32,
    f32,
    i64,
    f64,
    i128,
    string_header,
    binary_header,
    varchar_header
};

#define DATA_TYPE_STRING(op) \
    case data_type_t::op: return #op;

inline const char *data_type_to_string(data_type_t types) {
    switch (types) {
        DATA_TYPE_STRING(i8)
        DATA_TYPE_STRING(i16)
        DATA_TYPE_STRING(i32)
        DATA_TYPE_STRING(f32)
        DATA_TYPE_STRING(i64)
        DATA_TYPE_STRING(f64)
        DATA_TYPE_STRING(i128)
        DATA_TYPE_STRING(string_header)
        DATA_TYPE_STRING(binary_header)
        DATA_TYPE_STRING(varchar_header)
        default:
            return "Unknown";
    }
}

enum class data_kind_t : uint8_t {
    kMemory,
    kConst,
};

enum class opcodes : int32_t {
    Inv = -1,
    Ret = 0,
    Imm = 1,
    Mem = 2,
    Var = 3,
    Neg = 4,
    Not = 5,
    And = 6,
    Or = 7,
    Eq = 8,
    Ne,
    Lt,
    Le,
    Gt,
    Ge,
    Add,
    Sub,
    Mul,
    Div,
    Rem,
};

#define OPCODE_STRING(op) \
    case opcodes::op: return #op;

inline const char *opcode_to_string(opcodes opcode) {
    switch (opcode) {
        OPCODE_STRING(Inv)
        OPCODE_STRING(Ret)
        OPCODE_STRING(Imm)
        OPCODE_STRING(Mem)
        OPCODE_STRING(Var)
        OPCODE_STRING(Neg)
        OPCODE_STRING(Not)
        OPCODE_STRING(And)
        OPCODE_STRING(Or)
        OPCODE_STRING(Eq)
        OPCODE_STRING(Ne)
        OPCODE_STRING(Lt)
        OPCODE_STRING(Le)
        OPCODE_STRING(Gt)
        OPCODE_STRING(Ge)
        OPCODE_STRING(Add)
        OPCODE_STRING(Sub)
        OPCODE_STRING(Mul)
        OPCODE_STRING(Div)
        OPCODE_STRING(Rem)
        default:
            return "Unknown";
    }
}

struct instruction_t {
    opcodes opcode;
    int32_t options;
    union {
        struct {
            int64_t lo;
            int64_t hi;
        } ipayload;
        double dpayload;
    };

    void to_string(char* buf) const {
        data_type_t type = static_cast<data_type_t>(options);
        if (opcode == opcodes::Imm) {
            if (type == data_type_t::f32 || type == data_type_t::f64) {
                sprintf(buf, "[%s] (%s) %f", opcode_to_string(opcode), data_type_to_string(type), dpayload);
            } else {
                sprintf(buf, "[%s] (%s) %ld %ld", opcode_to_string(opcode), data_type_to_string(type), ipayload.hi, ipayload.lo);
            }
        } else if ((opcode == opcodes::Mem) || (opcode == opcodes::Var)) {
            sprintf(buf, "[%s] (%s) %ld", opcode_to_string(opcode), data_type_to_string(type), ipayload.lo);
        } else {
            sprintf(buf, "[%s]", opcode_to_string(opcode));
        }
    }
};

struct jit_value_t {

    inline jit_value_t() noexcept
            : op_(), type_(), kind_() {}

    inline jit_value_t(asmjit::Operand op, data_type_t type, data_kind_t kind) noexcept
            : op_(op), type_(type), kind_(kind) {}

    inline jit_value_t(const jit_value_t &other) noexcept = default;

    inline jit_value_t &operator=(const jit_value_t &other) noexcept = default;

#ifndef __aarch64__
    inline const asmjit::x86::Ymm &ymm() const noexcept { return op_.as<asmjit::x86::Ymm>(); }

    inline const asmjit::x86::Xmm &xmm() const noexcept { return op_.as<asmjit::x86::Xmm>(); }

    inline const asmjit::x86::Gpq &gp() const noexcept { return op_.as<asmjit::x86::Gpq>(); }
#else
    inline bool                   isVec() const noexcept { return op_.isVec(); }
    inline const asmjit::a64::Vec &vec() const noexcept { return op_.as<asmjit::a64::Vec>(); }

    inline bool                   isGp() const noexcept { return op_.isGp(); }
    inline const asmjit::a64::Gp &gp() const noexcept { return op_.as<asmjit::a64::Gp>(); }

    inline bool                   isImm() const noexcept { return op_.isImm(); }
    inline const asmjit::Imm      &imm() const noexcept { return op_.as<asmjit::Imm>(); }

    inline bool                   isMem() const noexcept { return op_.isMem(); }
    inline const asmjit::a64::Mem &mem() const noexcept { return op_.as<asmjit::a64::Mem>(); }
#endif

    inline data_type_t dtype() const noexcept { return type_; }

    inline data_kind_t dkind() const noexcept { return kind_; }

    inline const asmjit::Operand &op() const noexcept { return op_; }

    inline const char *to_string(asmjit::BaseCompiler &c) const {
        asmjit::StringTmp<512> sb;
        asmjit::Formatter::formatOperand(sb, {}, &c, asmjit::Arch::kAArch64, op_);
        return sb.data();
    }

private:
    asmjit::Operand op_;
    data_type_t type_;
    data_kind_t kind_;
};

inline uint32_t type_shift(data_type_t type) {
    switch (type) {
        case data_type_t::i8:
            return 0;
        case data_type_t::i16:
            return 1;
        case data_type_t::i32:
        case data_type_t::f32:
            return 2;
        case data_type_t::i64:
        case data_type_t::f64:
            return 3;
        case data_type_t::i128:
            return 4;
        default:
            __builtin_unreachable();
    }
}

inline data_kind_t dst_kind(const jit_value_t &lhs, const jit_value_t &rhs) {
    auto dk = (lhs.dkind() == data_kind_t::kConst && rhs.dkind() == data_kind_t::kConst) ? data_kind_t::kConst
                                                                                         : data_kind_t::kMemory;
    return dk;
}

inline void comment(asmjit::BaseCompiler &c, const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    char buffer[256];
    vsnprintf(buffer, sizeof(buffer), fmt, args);
    c.comment(buffer);
    va_end(args);
}

#endif //QUESTDB_JIT_COMMON_H
