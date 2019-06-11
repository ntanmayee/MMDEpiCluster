//
// Created by David Helekal on 23/10/2018.
//

#ifndef MMDIFF3_MMDIFF_3_CPP
#define MMDIFF3_MMDIFF_3_CPP

#import "mmdiff3.hpp"

#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

#include <stdlib.h>
#include <tuple>
#include <vector>

#include "mmd.hpp"
#include "rbf_joint_discrete_kernel.hpp"


using namespace mmdiff3;

SEXP compute_jmmd(SEXP a1, SEXP a2, SEXP b1, SEXP b2, SEXP maxval, SEXP LUT) {
        int *ra1 = INTEGER(a1);
        int *ia2 = INTEGER(a2);

        int *rb1 = INTEGER(b1);
        int *ib2 = INTEGER(b2);

        double *lut = REAL(LUT);

        auto imaxv = INTEGER(maxval)[0];


        size_t maxv = imaxv;

        auto ker = rbf_joint_discrete_kernel(lut, maxv);

        mmd<std::tuple<int, int> > run_mmd;
        std::vector<std::tuple<int, int>> vec1, vec2;

        for (int i = 0; i < length(a1); ++i) {
            vec1.emplace_back(std::make_tuple(ra1[i], ia2[i]));
        }

        for (int i = 0; i < length(b1); ++i) {
            vec2.emplace_back(std::make_tuple(rb1[i], ib2[i]));
        }

        SEXP result;
        result = PROTECT(allocVector(REALSXP, 1));
        double* rresult = REAL(result);
        double mmd_res = run_mmd.compute_mmd(vec1, vec2, ker);
        rresult[0] = mmd_res;
        UNPROTECT(1);
        return result;
}

SEXP compute_LUT(SEXP maxval, SEXP sigma) {
        auto imax = INTEGER(maxval)[0];
        const double rsigma = REAL(sigma)[0];
        size_t maxv = imax;

        SEXP result;
        result = PROTECT(allocVector(REALSXP, maxv+1));
        double* rresult = REAL(result);
        rbf_joint_discrete_kernel::compute_LUT(maxv, rsigma, rresult);
        UNPROTECT(1);
        return result;
}
#endif //MMDIFF3_MMDIFF_3_HPP
