#ifndef GPU_KERNELS_H
#define GPU_KERNELS_H

#include "parser.h" 
#include <cuda_runtime.h>
#include <cusparse.h>
#include <type_traits>
#include <cstdio>
#include <cstdlib>


template <typename T>
static inline void cuda_check(T status) {
    if constexpr (std::is_same_v<T, cudaError_t>) {
        if (status != cudaSuccess) {
            fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(status));
            std::exit(EXIT_FAILURE);
        }
    } else if constexpr (std::is_same_v<T, cusparseStatus_t>) {
        if (status != CUSPARSE_STATUS_SUCCESS) {
            fprintf(stderr, "cuSPARSE error: status=%d\n", static_cast<int>(status));
            std::exit(EXIT_FAILURE);
        }
    } else {
        if (static_cast<int>(status) != 0) {
            fprintf(stderr, "Error: status=%d\n", static_cast<int>(status));
            std::exit(EXIT_FAILURE);
        }
    }
}

template <typename MatrixT>
void SpMV_cuda_basic(int warmup, int niter, const MatrixT& matrix,
                     const float* d_random_vector,
                     float* d_results,
                     double* runs_time);

extern template void SpMV_cuda_basic<COO>(
    int, int, const COO&, const float*, float*, double*
);

extern template void SpMV_cuda_basic<CSR>(
    int, int, const CSR&, const float*, float*, double*
);

template <typename MatrixT>
void SpMV_cuda_custom_1(int warmup, int niter, const MatrixT& matrix,
                        const float* d_random_vector,
                        float* d_results,
                        double* runs_time);

extern template void SpMV_cuda_custom_1<COO>(
    int, int, const COO&, const float*, float*, double*
);

extern template void SpMV_cuda_custom_1<CSR>(
    int, int, const CSR&, const float*, float*, double*
);

void SpMV_cuda_custom_2(int warmup, int niter, const CSR& matrix,
                        const float* d_random_vector, float* d_results,
                        double* runs_time, double* setup_time);

template <typename MatrixT>
void SpMV_cuda_cuSparse(int warmup, int niter, const MatrixT& matrix,
                        const float* d_random_vector,
                        float* d_results,
                        double* runs_time);

extern template void SpMV_cuda_cuSparse<COO>(
    int, int, const COO&, const float*, float*, double*
);

extern template void SpMV_cuda_cuSparse<CSR>(
    int, int, const CSR&, const float*, float*, double*
);


COO coo_upload_to_device(const COO& host_matrix);
CSR csr_upload_to_device(const CSR& host_matrix);
#endif
