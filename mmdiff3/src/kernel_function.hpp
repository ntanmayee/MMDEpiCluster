//
// Created by David Helekal on 2019-01-04.
//

#ifndef MMDIFF3_KERNEL_FUNCTION_HPP
#define MMDIFF3_KERNEL_FUNCTION_HPP


#include <tuple>
namespace mmdiff3 {
    template<class T>
    class kernel_function {

    public:
        virtual double compute_kernel(T &a, T &b) const = 0;
    };
}
#endif //MMDIFF3_KERNEL_FUNCTION_HPP
