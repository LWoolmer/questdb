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

#ifndef QUESTDB_JIT_IMPL_A64_H
#define QUESTDB_JIT_IMPL_A64_H

#include "consts.h"
#include <cassert>

namespace questdb::a64 {
    using namespace asmjit;
    using namespace asmjit::a64;
    using namespace asmjit::arm;

    inline Gp to_int64(Compiler &c, const Gp &gp, bool check_null) {
        Gp gp64 = c.newInt64();
        c.sxtw(gp64, gp);
        return gp64;

        // REVISIT int->int conversion check_null
    }

    inline Gp to_int64(Compiler &c, const Vec &fp) {
        Gp gp64 = c.newInt64();
        c.fcvtzs(gp64, fp);
        return gp64;

        // REVISIT float->int conversion check_null
    }

    inline Vec to_float(Compiler &c, const Gp &gp, bool check_null) {
        Vec fp = c.newVecS();
        c.scvtf(fp, gp);
        return fp;

        // REVISIT int->float conversion check_null
    }

    inline Vec to_double(Compiler &c, const Gp &gp, bool check_null) {
        Vec fp = c.newVecD();
        c.scvtf(fp, gp);
        return fp;
        // REVISIT int->double conversion check_null
    }

    inline Vec to_double(Compiler &c, const Vec &fps) {
        Vec fpd = c.newVecD();
        c.fcvt(fpd, fps);
        return fpd;
        // REVISIT float->double conversion check_null
    }

    // inline void check_int32_null(Compiler &c, const Gp &dst, const Gpd &lhs, const Gpd &rhs) {
    //     c.cmp(lhs, INT_NULL);
    //     c.cmove(dst, lhs);
    //     c.cmp(rhs, INT_NULL);
    //     c.cmove(dst, rhs);
    // }

    // inline Gpd int32_neg(Compiler &c, const Gpd &rhs, bool check_null) {
    //     c.comment("int32_neg");

    //     Gp r = c.newInt32();
    //     c.mov(r, rhs);
    //     c.neg(r);
    //     if (check_null) {
    //         Gp t = c.newInt32();
    //         c.mov(t, INT_NULL);
    //         c.cmp(rhs, t);
    //         c.cmove(r, t);
    //     }
    //     return r.as<Gpd>();
    // }

    // inline Gpq int64_neg(Compiler &c, const Gpq &rhs, bool check_null) {
    //     c.comment("int64_neg");

    //     Gp r = c.newInt64();
    //     c.mov(r, rhs);
    //     c.neg(r);
    //     if (check_null) {
    //         Gp t = c.newInt64();
    //         c.movabs(t, LONG_NULL);
    //         c.cmp(rhs, t);
    //         c.cmove(r, rhs);
    //     }
    //     return r.as<Gpq>();
    // }

    // inline Gpd int32_add(Compiler &c, const Gpd &lhs, const Gpd &rhs, bool check_null) {
    //     c.comment("int32_add");

    //     Gp r = c.newInt64();
    //     c.lea(r, ptr(lhs, rhs));
    //     if (check_null) check_int32_null(c, r, lhs, rhs);
    //     return r.as<Gpd>();
    // }

    // inline Gpd int32_sub(Compiler &c, const Gpd &lhs, const Gpd &rhs, bool check_null) {
    //     c.comment("int32_sub");

    //     Gp r = c.newInt32();
    //     c.mov(r, lhs);
    //     c.sub(r, rhs);
    //     if (check_null) check_int32_null(c, r, lhs, rhs);
    //     return r.as<Gpd>();
    // }

    // inline Gpd int32_mul(Compiler &c, const Gpd &lhs, const Gpd &rhs, bool check_null) {
    //     c.comment("int32_mul");

    //     Gp r = c.newInt32();
    //     c.mov(r, lhs);
    //     c.imul(r, rhs);
    //     if (check_null) check_int32_null(c, r, lhs, rhs);
    //     return r.as<Gpd>();
    // }

    // inline Gpd int32_div(Compiler &c, const Gpd &lhs, const Gpd &rhs, bool check_null) {
    //     c.comment("int32_div");

    //     Label l_null = c.newLabel();
    //     Label l_exit = c.newLabel();

    //     Gp r = c.newInt32();
    //     Gp t = c.newInt32();

    //     if (!check_null) {
    //         c.mov(r, lhs);
    //         c.test(rhs, rhs);
    //         c.je(l_null);
    //         c.cdq(t, r);
    //         c.idiv(t, r, rhs);
    //         c.jmp(l_exit);
    //         c.bind(l_null);
    //         c.mov(r, INT_NULL);
    //         c.bind(l_exit);
    //         return r.as<Gpd>();
    //     }
    //     c.mov(r, INT_NULL);
    //     c.test(rhs, 2147483647); //INT_NULL - 1
    //     c.je(l_null);
    //     c.cmp(lhs, INT_NULL);
    //     c.je(l_null);
    //     c.mov(r, lhs);
    //     c.cdq(t, r);
    //     c.idiv(t, r, rhs);
    //     c.bind(l_null);
    //     return r.as<Gpd>();
    // }

    // inline void check_int64_null(Compiler &c, const Gp &dst, const Gp &lhs, const Gp &rhs) {
    //     c.comment("check_int64_null");
    //     Gp n = c.newGpq();
    //     c.movabs(n, LONG_NULL);
    //     c.cmp(lhs, n);
    //     c.cmove(dst, lhs);
    //     c.cmp(rhs, n);
    //     c.cmove(dst, rhs);
    // }

    // inline Gpq int64_add(Compiler &c, const Gpq &lhs, const Gpq &rhs, bool check_null) {
    //     c.comment("int64_add");

    //     Gp r = c.newInt64();
    //     c.lea(r, ptr(lhs, rhs));
    //     if (check_null) check_int64_null(c, r, lhs, rhs);
    //     return r.as<Gpq>();
    // }

    // inline Gpq int64_sub(Compiler &c, const Gpq &lhs, const Gpq &rhs, bool check_null) {
    //     c.comment("int64_sub");
    //     Gp r = c.newInt64();
    //     c.mov(r, lhs);
    //     c.sub(r, rhs);
    //     if (check_null) check_int64_null(c, r, lhs, rhs);
    //     return r.as<Gpq>();
    // }

    // inline Gpq int64_mul(Compiler &c, const Gpq &lhs, const Gpq &rhs, bool check_null) {
    //     c.comment("int64_mul");
    //     Gp r = c.newInt64();
    //     c.mov(r, lhs);
    //     c.imul(r, rhs);
    //     if (check_null) check_int64_null(c, r, lhs, rhs);
    //     return r.as<Gpq>();
    // }

    // inline Gpq int64_div(Compiler &c, const Gpq &lhs, const Gpq &rhs, bool check_null) {
    //     c.comment("int64_div");

    //     Label l_null = c.newLabel();
    //     Label l_exit = c.newLabel();

    //     Gp r = c.newInt64();
    //     Gp t = c.newInt64();
    //     if (!check_null) {
    //         c.mov(r, lhs);
    //         c.test(rhs, rhs);
    //         c.je(l_null);
    //         c.cqo(t, r);
    //         c.idiv(t, r, rhs);
    //         c.jmp(l_exit);
    //         c.bind(l_null);
    //         c.movabs(r, LONG_NULL);
    //         c.bind(l_exit);
    //         return r.as<Gpq>();
    //     }
    //     c.mov(t, rhs);
    //     c.mov(r, lhs);
    //     c.btr(t, 63);
    //     c.test(t, t);
    //     c.je(l_null);
    //     c.movabs(t, LONG_NULL);
    //     c.cmp(lhs, t);
    //     c.je(l_null);
    //     c.cqo(t, r);
    //     c.idiv(t, r, rhs);
    //     c.jmp(l_exit);

    //     c.bind(l_null);
    //     c.movabs(r, LONG_NULL);
    //     c.bind(l_exit);
    //     return r.as<Gpq>();
    // }

    // inline Xmm float_neg(Compiler &c, const Xmm &rhs) {
    //     int32_t array[4] = {INT_NULL, 0, 0, 0};
    //     Mem mem = c.newConst(ConstPoolScope::kLocal, &array, 32);
    //     c.xorps(rhs, mem);
    //     return rhs;
    // }

    // inline Xmm double_neg(Compiler &c, const Xmm &rhs) {
    //     int32_t array[4] = {0, INT_NULL, 0, 0};
    //     Mem mem = c.newConst(ConstPoolScope::kLocal, &array, 32);
    //     c.xorpd(rhs, mem);
    //     return rhs;
    // }

    // inline Xmm float_add(Compiler &c, const Xmm &lhs, const Xmm &rhs) {
    //     c.addss(lhs, rhs);
    //     return lhs;
    // }

    // inline Xmm float_sub(Compiler &c, const Xmm &lhs, const Xmm &rhs) {
    //     c.subss(lhs, rhs);
    //     return lhs;
    // }

    // inline Xmm float_mul(Compiler &c, const Xmm &lhs, const Xmm &rhs) {
    //     c.mulss(lhs, rhs);
    //     return lhs;
    // }

    // inline Xmm float_div(Compiler &c, const Xmm &lhs, const Xmm &rhs) {
    //     c.divss(lhs, rhs);
    //     return lhs;
    // }

    // inline Xmm double_add(Compiler &c, const Xmm &lhs, const Xmm &rhs) {
    //     c.addsd(lhs, rhs);
    //     return lhs;
    // }

    // inline Xmm double_sub(Compiler &c, const Xmm &lhs, const Xmm &rhs) {
    //     c.subsd(lhs, rhs);
    //     return lhs;
    // }

    // inline Xmm double_mul(Compiler &c, const Xmm &lhs, const Xmm &rhs) {
    //     c.mulsd(lhs, rhs);
    //     return lhs;
    // }

    // inline Xmm double_div(Compiler &c, const Xmm &lhs, const Xmm &rhs) {
    //     c.divsd(lhs, rhs);
    //     return lhs;
    // }

    inline Gp cmp(Compiler &c, const Gp &lhs, const Gp &rhs, CondCode cond, bool check_null) {
        c.cmp(lhs, rhs);

        if (check_null) {
            // REVISIT non-int null comparison
            c.ccmp(lhs, INT_NULL, 0, CondCode::kNE);
            c.ccmp(rhs, INT_NULL, 0, CondCode::kNE);
        }

        Gp cmp = c.newSimilarReg(lhs);
        c.cset(cmp, cond);

        return cmp;
    }

    inline Gp cmp(Compiler &c, const Vec &lhs, const Vec &rhs, CondCode cond, bool check_null) {
        assert(false);
        __builtin_unreachable();
        // REVISIT cmp floating point and 128b
    }

}
#endif //QUESTDB_JIT_IMPL_A64_H
