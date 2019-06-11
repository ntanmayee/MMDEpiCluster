//
// Created by David Helekal on 2019-01-03.
//

#include "rbf_joint_discrete_kernel.hpp"
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif


namespace mmdiff3 {

    rbf_joint_discrete_kernel::~rbf_joint_discrete_kernel() {
    }

    rbf_joint_discrete_kernel::rbf_joint_discrete_kernel(double* LUT, size_t maxval) {
        this->max_dist = maxval;
        this->lookup = LUT;
    }

    inline double rbf_joint_discrete_kernel::compute_kernel(std::tuple<int, int> &a, std::tuple<int, int> &b) const {
        int pos_a;
        int cat_a;

        int pos_b;
        int cat_b;

        double result = 0.0;

        std::tie(pos_a, cat_a) = a;
        std::tie(pos_b, cat_b) = b;

        size_t loc = abs(pos_a - pos_b);
        assert(loc <= this->max_dist);

        if (cat_a == cat_b) result = this->lookup[loc];

        return result;
    }

    void rbf_joint_discrete_kernel::compute_LUT(size_t maxval, double sigma, double* lut) {

        size_t i = 0;
        double val;

#if defined(_OPENMP)
#pragma omp parallel for \
        default(shared) \
        private(i, val) \
        schedule(dynamic, 100)

        for(i=0; i <= maxval; ++i) {
            assert(!std::isnan(sigma));
            val = exp((-1/(2*pow(sigma,2))) * pow(i, 2));
            lut[i]=val;
        }
#else
        for(i=0; i <= maxval; ++i) {
            assert(!std::isnan(sigma));
            val = exp((-1/(2*pow(sigma,2))) * pow(i, 2));
            lut[i]=val;
        }
#endif
    }
}