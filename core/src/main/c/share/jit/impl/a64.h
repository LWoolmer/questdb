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
#include <stdexcept>

namespace questdb::a64 {
    using namespace asmjit;
    using namespace asmjit::a64;
    using namespace asmjit::arm;

    inline void cmp_null(Compiler &c, const Gp &gp) {
        comment(c, "cmp_null(gp)");
        if (gp.isGp64()) {
            Gp null = c.newInt64();
            c.mov(null, LONG_NULL);
            c.cmp(gp, null);
        } else {
            Gp null = c.newInt32();
            c.mov(null, INT_NULL);
            c.cmp(gp, null);
        }
    }

    inline void cmp_null(Compiler &c, const Vec &fp) {
        comment(c, "cmp_null(vec)");
        if (fp.isVec32()) {
            Gp gp = c.newInt32();
            Gp nan = c.newInt32();
            c.mov(nan, FLOAT_NAN);
            c.fmov(gp, fp);
            c.and_(gp, gp, nan);
            c.cmp(gp, nan);
        } else if (fp.isVec64()) {
            Gp gp = c.newInt64();
            Gp nan = c.newInt64();
            c.mov(nan, DOUBLE_NAN);
            c.fmov(gp, fp);
            c.and_(gp, gp, nan);
            c.cmp(gp, nan);
        } else {
            UNREACHABLE("Unsupported type for cmp_null(vec)");
        }
    }

    inline void cmp_eq_epsilon(Compiler &c, const Vec &lhs, const Vec &rhs) {
        comment(c, "cmp(eq, +-eps)");

        Mem eps;
        if (lhs.isVec32()) {
            eps = c.newFloatConst(ConstPoolScope::kLocal, FLOAT_EPSILON);
        } else {
            eps = c.newDoubleConst(ConstPoolScope::kLocal, DOUBLE_EPSILON);
        }

        Vec eps_vec = c.newSimilarReg(lhs);
        c.ldr(eps_vec, eps);

        Vec diff = c.newSimilarReg(lhs);
        c.fsub(diff, lhs, rhs);
        c.fabs(diff, diff);
        c.fcmp(diff, eps_vec);
    }

    inline Gp to_int64(Compiler &c, const Gp &gp, bool check_null) {
        Gp gp64 = c.newInt64();
        c.sxtw(gp64, gp);
        if (check_null) {
            Gp null = c.newInt64();
            c.mov(null, LONG_NULL);
            cmp_null(c, gp);
            c.csel(gp64, gp64, null, CondCode::kNE);
        }
        return gp64;
    }

    inline Gp to_int64(Compiler &c, const Vec &fp, bool check_null) {
        Gp gp64 = c.newInt64();
        c.fcvtzs(gp64, fp);
        return gp64;
    }

    inline Vec to_float(Compiler &c, const Gp &gp, bool check_null) {
        Vec fp = c.newVecS();
        c.scvtf(fp, gp);
        if (check_null) {
            Gp gpnan = c.newInt32();
            Vec nan = c.newVecS();
            c.mov(gpnan, FLOAT_NAN);
            c.fmov(nan, gpnan);
            cmp_null(c, gp);
            c.fcsel(fp, fp, nan, CondCode::kNE);
        }
        return fp;
    }

    inline Vec to_double(Compiler &c, const Gp &gp, bool check_null) {
        Vec fp = c.newVecD();
        c.scvtf(fp, gp);
        if (check_null) {
            Gp gpnan = c.newInt64();
            Vec nan = c.newVecD();
            c.mov(gpnan, DOUBLE_NAN);
            c.fmov(nan, gpnan);
            cmp_null(c, gp);
            c.fcsel(fp, fp, nan, CondCode::kNE);
        }
        return fp;
    }

    inline Vec to_double(Compiler &c, const Vec &fps, bool check_null) {
        Vec fpd = c.newVecD();
        c.fcvt(fpd, fps);
        return fpd;
    }

    inline Gp bin_not(Compiler &c, const Gp &lhs) {
        Gp result = c.newSimilarReg(lhs);
        c.eor(result, lhs, -1);
        return result;
    }

    inline Gp bin_and(Compiler &c, const Gp &lhs, const Gp &rhs) {
        Gp result = c.newSimilarReg(lhs);
        c.and_(result, lhs, rhs);
        return result;
    }

    inline Gp bin_or(Compiler &c, const Gp &lhs, const Gp &rhs) {
        Gp result = c.newSimilarReg(lhs);
        c.orr(result, lhs, rhs);
        return result;
    }

    inline Gp neg(Compiler &c, const Gp &lhs, bool check_null) {
        comment(c, "add(gp)");
        Gp result = c.newSimilarReg(lhs);
        c.neg(result, lhs);

        if (check_null) {
            cmp_null(c, lhs);
            c.csel(result, result, lhs, CondCode::kNE);
        }
        return result;
    }

    inline Vec neg(Compiler &c, const Vec &lhs, bool check_null) {
        comment(c, "add(gp)");
        if (lhs.type() != lhs.type()) {
            UNREACHABLE("neg called with different types");
        }

        Vec result = c.newSimilarReg(lhs);
        c.fneg(result, lhs);

        if (check_null) {
            cmp_null(c, lhs);
            c.fcsel(result, result, lhs, CondCode::kNE);
        }
        return result;
    }

    inline Gp add(Compiler &c, const Gp &lhs, const Gp &rhs, bool check_null) {
        comment(c, "add(gp)");
        if (lhs.type() != rhs.type()) {
            UNREACHABLE("add(gp) called with different types");
        }

        Gp result = c.newSimilarReg(lhs);
        c.add(result, lhs, rhs);

        if (check_null) {
            cmp_null(c, lhs);
            c.csel(result, result, lhs, CondCode::kNE);
            cmp_null(c, rhs);
            c.csel(result, result, rhs, CondCode::kNE);
        }
        return result;
    }

    inline Vec add(Compiler &c, const Vec &lhs, const Vec &rhs, bool check_null) {
        comment(c, "add(vec)");
        if (lhs.type() != rhs.type()) {
            UNREACHABLE("add(vec) called with different types");
        }

        Vec result = c.newSimilarReg(lhs);
        c.fadd(result, lhs, rhs);

        if (check_null) {
            cmp_null(c, lhs);
            c.fcsel(result, result, lhs, CondCode::kNE);
            cmp_null(c, rhs);
            c.fcsel(result, result, rhs, CondCode::kNE);
        }
        return result;
    }

    inline Gp sub(Compiler &c, const Gp &lhs, const Gp &rhs, bool check_null) {
        comment(c, "sub(gp)");
        if (lhs.type() != rhs.type()) {
            UNREACHABLE("sub(gp) called with different types");
        }

        Gp result = c.newSimilarReg(lhs);
        c.sub(result, lhs, rhs);

        if (check_null) {
            cmp_null(c, lhs);
            c.csel(result, result, lhs, CondCode::kNE);
            cmp_null(c, rhs);
            c.csel(result, result, rhs, CondCode::kNE);
        }
        return result;
    }

    inline Vec sub(Compiler &c, const Vec &lhs, const Vec &rhs, bool check_null) {
        comment(c, "sub(vec)");
        if (lhs.type() != rhs.type()) {
            UNREACHABLE("sub(vec) called with different types");
        }

        Vec result = c.newSimilarReg(lhs);
        c.fsub(result, lhs, rhs);

        if (check_null) {
            cmp_null(c, lhs);
            c.fcsel(result, result, lhs, CondCode::kNE);
            cmp_null(c, rhs);
            c.fcsel(result, result, rhs, CondCode::kNE);
        }
        return result;
    }


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

    // inline Xmm float_mul(Compiler &c, const Xmm &lhs, const Xmm &rhs) {
    //     c.mulss(lhs, rhs);
    //     return lhs;
    // }

    // inline Xmm float_div(Compiler &c, const Xmm &lhs, const Xmm &rhs) {
    //     c.divss(lhs, rhs);
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

    // REVISIT not optimized
    inline Gp cmp(Compiler &c, const Gp &lhs, const Gp &rhs, CondCode cond, bool check_null) {
        comment(c, "cmp(gp)");
        if (lhs.type() != rhs.type()) {
            UNREACHABLE("cmp called with different types");
        }

        Gp result = c.newInt64();
        c.cmp(lhs, rhs);
        c.cset(result, cond);

        if (!check_null || (cond == CondCode::kEQ) || (cond == CondCode::kNE)) {
            return result;
        }

        // check_null && (LT || LE || GT || GE)

        Gp lhs_null = c.newInt64();
        Gp rhs_null = c.newInt64();

        cmp_null(c, lhs);
        c.cset(lhs_null, CondCode::kEQ);
        cmp_null(c, rhs);
        c.cset(rhs_null, CondCode::kEQ);

        c.eor(lhs_null, lhs_null, rhs_null);
        c.mvn(lhs_null, lhs_null);

        c.and_(result, result, lhs_null);
        return result;
    }

    // REVISIT not optimized
    inline Gp cmp(Compiler &c, const Vec &lhs, const Vec &rhs, CondCode cond) {
        comment(c, "cmp(vec)");
        if (lhs.type() != rhs.type()) {
            UNREACHABLE("cmp called with different types");
        }
        // REVISIT very very not optimized
        if (lhs.isVec128()) {
            Gp lhs_lo = c.newInt64();
            Gp rhs_lo = c.newInt64();
            c.mov(lhs_lo, lhs.v64().at(0));
            c.mov(rhs_lo, rhs.v64().at(0));
            c.subs(lhs_lo, lhs_lo, rhs_lo);

            Gp rhs_hi = c.newInt64();
            Gp lhs_hi = c.newInt64();
            c.mov(lhs_hi, lhs.v64().at(1));
            c.mov(rhs_hi, rhs.v64().at(1));
            c.sbcs(lhs_hi, lhs_hi, rhs_hi);

            Gp result = c.newInt64();

            if (cond == CondCode::kEQ || cond == CondCode::kNE) {
                c.cset(result, cond);
                return result;
            } else {
                UNREACHABLE("128bit non-eq cmp not supported");
            }
            return result;
        }
        Gp result = c.newInt64();
        Gp eq_epsilon = c.newInt64();

        cmp_eq_epsilon(c, lhs, rhs);
        c.cset(eq_epsilon, CondCode::kLE);

        // This is a bit tricky:
        //  - if EQ then use the epsilon comparison
        //  - if NE then use the negation of the epsilon comparison
        //  - if LE/GE, then use the float comparison OR epsilon comparison
        //  - if LT/GT, then do the float comparison AND NOT epsilon comparison
        if (cond == CondCode::kEQ) {
            result = eq_epsilon;
        } else if (cond == CondCode::kNE) {
            c.mvn(result, eq_epsilon);
        } else {
            c.fcmp(lhs, rhs);
            c.cset(result, cond);
            if ((cond == CondCode::kLE) || (cond == CondCode::kGE)) {
                c.orr(result, result, eq_epsilon);
            } else {
                c.mvn(eq_epsilon, eq_epsilon);
                c.and_(result, result, eq_epsilon);
            }
        }

        Gp lhs_null = c.newInt64();
        Gp rhs_null = c.newInt64();

        cmp_null(c, lhs);
        c.cset(lhs_null, CondCode::kEQ);
        cmp_null(c, rhs);
        c.cset(rhs_null, CondCode::kEQ);

        c.eor(lhs_null, lhs_null, rhs_null);

        if (cond == CondCode::kNE) {
            c.orr(result, result, lhs_null);
        } else {
            c.mvn(lhs_null, lhs_null);
            c.and_(result, result, lhs_null);
        }
        return result;
    }

}
#endif //QUESTDB_JIT_IMPL_A64_H
