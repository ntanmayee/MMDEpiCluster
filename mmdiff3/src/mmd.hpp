//
// Created by David Helekal on 2018-12-24.
//
#ifndef MMDIFF3_MMD_HPP
#define MMDIFF3_MMD_HPP

#include "vector"
#include "kernel_function.hpp"


namespace mmdiff3 {
    template<class T>
    class mmd {
    public:
        mmd();

        ~mmd();

        double compute_mmd(std::vector<T> &x, std::vector<T> &y, kernel_function<T> &k);

    private:
        inline double kernel_sum(std::vector<T> &x,
                std::vector<T> &y,
                kernel_function<T> &k,
                bool no_diag = false);
    };
}
#endif //MMDIFF3_MMD_HPP
