// we only include RcppArrayFire.h which pulls Rcpp.h in for us
#include "RcppArrayFire.h"
#include <arrayfire.h>
#include <cstdio>
#include <cstdlib>


// via the depends attribute we tell Rcpp to create hooks for
// RcppFire so that the build process will know what to do
//
// [[Rcpp::depends(RcppArrayFire)]]

// RcppArrayFire needs C++11
// add the following comment when you export your
// C++ function to R via Rcpp::SourceCpp()
// [[Rcpp::plugins(cpp11)]]

// simple example of creating two matrices and
// returning the result of an operation on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//


// [[Rcpp::export]]
af::array rcpparrayfire_hello_world() {
	af::array m1 = af::identity(3, 3, f32);
	af::array m2 = af::randn(3, 3, f32);

    return m1 + 3 * (m1 + m2);
}

using namespace af;

void testBackend()
{
  af::info();
  af_print(af::randu(5, 4));
}

//' @export
// [[Rcpp::export]]
int rcpparrayfire_test_backends()
  {
    try {
      printf("Trying CPU Backend\n");
      af::setBackend(AF_BACKEND_CPU);
      testBackend();
    } catch (af::exception& e) {
      printf("Caught exception when trying CPU backend\n");
      fprintf(stderr, "%s\n", e.what());
    }
    try {
      printf("Trying CUDA Backend\n");
      af::setBackend(AF_BACKEND_CUDA);
      testBackend();
    } catch (af::exception& e) {
      printf("Caught exception when trying CUDA backend\n");
      fprintf(stderr, "%s\n", e.what());
    }
    try {
      printf("Trying OpenCL Backend\n");
      af::setBackend(AF_BACKEND_OPENCL);
      testBackend();
    } catch (af::exception& e) {
      printf("Caught exception when trying OpenCL backend\n");
      fprintf(stderr, "%s\n", e.what());
    }
    return 0;
  }

