//
// Created by David Helekal on 2019-01-03.
//

#include "rbf_joint_kernel.hpp"
#include <cmath>

namespace mmdiff3 {
    rbf_joint_kernel::~rbf_joint_kernel() {

    }

    rbf_joint_kernel::rbf_joint_kernel(double sigma) {
        this->sigma = sigma;
    }

    inline double rbf_joint_kernel::compute_kernel(std::tuple<double, int> &a, std::tuple<double, int> &b) const {
        double pos_a;
        int cat_a;

        double pos_b;
        int cat_b;

        double result = 0.0;

        std::tie(pos_a, cat_a) = a;
        std::tie(pos_b, cat_b) = b;

        if (cat_a == cat_b) result = exp(-sigma * pow(pos_a - pos_b, 2));

        return result;
    }
}