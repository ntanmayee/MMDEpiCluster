//
// Created by David Helekal on 2019-01-03.
//

#ifndef MMDIFF3_rbf_joint_kernel_HPP
#define MMDIFF3_rbf_joint_kernel_HPP

#include <tuple>
#include "kernel_function.hpp"

namespace mmdiff3 {
    class rbf_joint_kernel : public kernel_function<std::tuple<double, int>> {
    public:
        double compute_kernel(std::tuple<double, int> &a, std::tuple<double, int> &b) const override;

        explicit rbf_joint_kernel(double sigma);

        ~rbf_joint_kernel();

    private:
        double sigma;
    };
}


#endif //MMDIFF3_KERNEL_HPP
