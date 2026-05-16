#ifndef UTILS_H
#define UTILS_H

#include <math.h>
#include <string>
#include <filesystem>
#include "parser.h"

struct results_struct {
    std::string name;
    double working_set_bytes;
    double parsing_time_ns;

    double coo_time_ns;
    double coo_bw_effective_gbs;
    double coo_parallel_time_ns;
    double coo_parallel_bw_effective_gbs;

    double coo_cuda_basic_time_ns;
    double coo_cuda_basic_bw_gbs;
    double coo_cuda_custom_1_time_ns;
    double coo_cuda_custom_1_bw_gbs;
    double coo_cuda_cuSparse_time_ns;
    double coo_cuda_cuSparse_bw_gbs;

    double coo_to_csr_time_ns;

    double csr_time_ns;
    double csr_bw_effective_gbs;
    double csr_parallel_time_ns;
    double csr_parallel_bw_effective_gbs;
    
    double csr_cuda_basic_time_ns;
    double csr_cuda_basic_bw_gbs;
    double csr_cuda_custom_1_time_ns;
    double csr_cuda_custom_1_bw_gbs;
    double csr_cuda_custom_2_time_ns;
    double csr_custom_2_setup_time_us;
    double csr_cuda_custom_2_bw_gbs;
    double csr_cuda_cuSparse_time_ns;
    double csr_cuda_cuSparse_bw_gbs;

    double total_time_ns;
};

struct results_vectors_struct {
    float* COO_results;
    float* COO_parallel_results;
    float* CSR_results;
    float* CSR_parallel_results;
    float* COO_results_cuda_basic;
    float* COO_results_cuda_custom_1;
    float* COO_results_cuda_cuSparse;
    float* CSR_results_cuda_basic;
    float* CSR_results_cuda_custom_1;
    float* CSR_results_cuda_custom_2;
    float* CSR_results_cuda_cuSparse;
};



double geometric_mean(const double *v, int len);
double arithmetic_mean(const double *v, int len);
double sigma_fn_sol(double *v, double mu, int len);
double working_set_bytes(const COO& m);
double ns_to_gbs(double bytes, double runtime_ns);
double coo_effective_bytes(const COO& m);
double csr_effective_bytes(const CSR& m);
CSR coo_to_csr(const COO& COO_matrix, double* convertion_time);
void results_comparing(int n_rows, const char* v1, const float* COO_results, 
                        const char* v2, const float* CSR_results);
void save_matrix_number(std::string input, int matrix_number, double workset_bytes,
                        const std::filesystem::path& output_path);
void append_spmv_summary_csv(const std::filesystem::path& output_path,
                            const results_struct& summary);
void data_printing(results_struct results, results_vectors_struct& res_vector, const COO& m);

#endif
