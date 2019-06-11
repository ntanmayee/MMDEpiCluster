//
// Created by David Helekal on 2019-01-16.
//

#ifndef MMDIFF3_MMDIFF3_HPP
#define MMDIFF3_MMDIFF3_HPP
#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>
extern "C" {
SEXP compute_jmmd(SEXP a1, SEXP a2, SEXP b1, SEXP b2, SEXP maxval, SEXP LUT);
SEXP compute_LUT(SEXP maxval, SEXP sigma);

}
#endif //MMDIFF3_MMDIFF3_HPP
