//
// Created by David Helekal on 2019-01-03.
//

#ifndef MMDIFF3_rbf_joint_discrete_kernel_HPP
#define MMDIFF3_rbf_joint_discrete_kernel_HPP

#include <tuple>
#include <vector>
#include "kernel_function.hpp"

namespace mmdiff3 {
    class rbf_joint_discrete_kernel : public kernel_function<std::tuple<int, int>> {
    public:

        static void compute_LUT(size_t maxval, double sigma, double* lut);

        double compute_kernel(std::tuple<int, int> &a, std::tuple<int, int> &b) const override;

        explicit rbf_joint_discrete_kernel(double* LUT, size_t maxval);

        ~rbf_joint_discrete_kernel();



    private:
        size_t max_dist;
        double *lookup;
    };
}


#endif //MMDIFF3_KERNEL_HPP
