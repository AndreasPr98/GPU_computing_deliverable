#include "utils.h"
#include "time_management.h"
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <filesystem>
#include <sys/time.h>

double arithmetic_mean(const double *v, int len) {

    double mu = 0.0;
    for (int i=0; i<len; i++)
        mu += (double)v[i];
    mu /= (double)len;

    return(mu);
}

double geometric_mean(const double *v, int len) {
    
    double mu = 1.0;
    for (int i=0; i<len; i++) {
        mu *= (v[i] > 0) ? ((double)v[i]) : 1;
    }
    mu = pow(mu, 1.0 / len);
    
    return(mu);
}

double sigma_fn_sol(double *v, double mu, int len) {

    double sigma = 0.0;
    for (int i=0; i<len; i++) {
        sigma += ((double)v[i] - mu)*((double)v[i] - mu);
    }
    sigma /= (double)len;

    return(sigma);
}

double working_set_bytes(const COO& m){
    return (double)m.nnz * (2.0 * sizeof(int) + sizeof(float)) +
            (double)m.n_cols * sizeof(float) + (double)m.n_rows * sizeof(float);
}

double ns_to_gbs(double bytes, double runtime_ns) {
    const double runtime_s = runtime_ns * 1e-9;
    if (runtime_s <= 0.0) { return 0.0; }
    return (bytes / runtime_s) / 1e9;
}

double coo_effective_bytes(const COO& m) {
    return (double)m.nnz * ((2.0 * sizeof(int)) + (4.0 * sizeof(float)));
}

double csr_effective_bytes(const CSR& m) {
    const double row_bytes = (double)m.n_rows * (2.0 * sizeof(int) + 1.0 * sizeof(float));
    const double nnz_bytes = (double)m.nnz * (1.0 * sizeof(int) + 4.0 * sizeof(float));
    return row_bytes + nnz_bytes;
}

// COO matrix to CSR matrix convertion
CSR coo_to_csr(const COO& COO_matrix, double* convertion_time){

    TIMER_DEF(convertion_timer);
    TIMER_START(convertion_timer);

    CSR CSR_matrix;
    CSR_matrix.n_rows = COO_matrix.n_rows;
    CSR_matrix.n_cols = COO_matrix.n_cols;
    CSR_matrix.nnz = COO_matrix.nnz;

    CSR_matrix.row_indices = (int*)calloc((CSR_matrix.n_rows + 1), sizeof(int));
    CSR_matrix.col_indices = (int*)malloc(CSR_matrix.nnz * sizeof(int));
    CSR_matrix.values = (float*)malloc(CSR_matrix.nnz * sizeof(float));
    if(!CSR_matrix.row_indices || !CSR_matrix.col_indices || !CSR_matrix.values) {
        fprintf(stderr, "COO_malloc_matrix allocation failed");
        return {};
    }

    for (int i = 0; i < CSR_matrix.nnz; i++) {
        CSR_matrix.row_indices[COO_matrix.row_indices[i] + 1]++;
    }
    for (int i = 0; i < CSR_matrix.n_rows; i++) {
        CSR_matrix.row_indices[i + 1] += CSR_matrix.row_indices[i];
    }

    for (int i = 0; i < CSR_matrix.nnz; i++) {
        CSR_matrix.col_indices[i] = COO_matrix.col_indices[i];
        CSR_matrix.values[i] = COO_matrix.values[i];
    }

    TIMER_STOP(convertion_timer);
    *convertion_time = TIMER_ELAPSED(convertion_timer);
    return CSR_matrix;
}

void results_comparing(int n_rows, const char* v1, const float* COO_results, 
                        const char* v2, const float* CSR_results){
    int equal = 1;
    
    const double abs_tol = 1e-2;
    const double rel_tol = 1e-2;
    for(int i = 0; i < n_rows; i++){
        double a = (double)COO_results[i];
        double b = (double)CSR_results[i];
        double diff = fabs(a - b);
        double scale = fmax(fabs(a), fabs(b));

        // Mixed tolerance is robust for both near-zero and large-magnitude values.
        if(diff > (abs_tol + rel_tol * scale)) {
            equal = 0; 
            printf(" Mismatch at row %d: %f vs %f", i, a, b);
            break;
        }
    }
    if(equal == 1) {printf(" %s = %s", v1, v2);} 
    else {printf(" %s != %s", v1, v2);}
}

void save_matrix_number(std::string input, int matrix_number, double workset_bytes,
                        const std::filesystem::path& output_path) {
    std::filesystem::create_directories(output_path.parent_path());
    std::ofstream outFile(output_path, std::ios::app);
    if (outFile.is_open()) {
        outFile << matrix_number << " -> " << input << " " 
                << std::fixed << std::setprecision(0) << workset_bytes << "\n";
        outFile.close();
    }
}

void append_spmv_summary_csv(const std::filesystem::path& output_path,
                            const results_struct& summary) {
    std::filesystem::create_directories(output_path.parent_path());
    const bool file_exists = std::filesystem::exists(output_path);
    const bool write_header = !file_exists || std::filesystem::file_size(output_path) == 0;

    std::ofstream out_file(output_path, std::ios::app);
    if (!out_file.is_open()) {
        fprintf(stderr, "Cannot open CSV output file: %s\n", output_path.string().c_str());
        return;
    }

    if (write_header) {
        out_file << "name,working_bytes,parsing_time_ms,"
                 << "coo_time_us,coo_bw_effective_gbs,"
                 << "coo_parallel_time_us,coo_parallel_bw_effective_gbs,"
                 << "coo_cuda_basic_time_us,coo_cuda_basic_bw_gbs,"
                 << "coo_cuda_custom_1_time_us,coo_cuda_custom_1_bw_gbs,"
                 << "coo_cuda_cuSparse_time_us,coo_cuda_cuSparse_bw_gbs,"
                 << "coo_to_csr_time_us,"
                 << "csr_time_us,csr_bw_effective_gbs,"
                 << "csr_parallel_time_us,csr_parallel_bw_effective_gbs,"
                 << "csr_cuda_basic_time_us,csr_cuda_basic_bw_gbs,"
                 << "csr_cuda_custom_1_time_us,csr_cuda_custom_1_bw_gbs,"
                 << "csr_cuda_custom_2_time_us,csr_cuda_custom_2_bw_gbs,"
                 << "csr_cuda_cuSparse_time_us,csr_cuda_cuSparse_bw_gbs,"
                 << "total_time_ms\n";
    }

    out_file << summary.name << ","
             << std::fixed << std::setprecision(0) << summary.working_set_bytes << ","
             << std::setprecision(6) << (summary.parsing_time_ns / 1e6) << ","
             << (summary.coo_time_ns / 1e3) << ","
             << summary.coo_bw_effective_gbs << ","
             << (summary.coo_parallel_time_ns / 1e3) << ","
             << summary.coo_parallel_bw_effective_gbs << ","
             << (summary.coo_cuda_basic_time_ns / 1e3) << ","
             << summary.coo_cuda_basic_bw_gbs << ","
             << (summary.coo_cuda_custom_1_time_ns / 1e3) << ","
             << summary.coo_cuda_custom_1_bw_gbs << ","
             << (summary.coo_cuda_cuSparse_time_ns / 1e3) << ","
             << summary.coo_cuda_cuSparse_bw_gbs << ","
             << (summary.coo_to_csr_time_ns / 1e3) << ","
             << (summary.csr_time_ns / 1e3) << ","
             << summary.csr_bw_effective_gbs << ","
             << (summary.csr_parallel_time_ns / 1e3) << ","
             << summary.csr_parallel_bw_effective_gbs << ","
             << (summary.csr_cuda_basic_time_ns / 1e3) << ","
             << summary.csr_cuda_basic_bw_gbs << ","
             << (summary.csr_cuda_custom_1_time_ns / 1e3) << ","
             << summary.csr_cuda_custom_1_bw_gbs << ","
             << (summary.csr_cuda_custom_2_time_ns / 1e3) << ","
             << summary.csr_cuda_custom_2_bw_gbs << ","
             << (summary.csr_cuda_cuSparse_time_ns / 1e3) << ","
             << summary.csr_cuda_cuSparse_bw_gbs << ","
             << (summary.total_time_ns / 1e6) << "\n";
}

void data_printing(results_struct results, results_vectors_struct& res_vector, const COO& m){
    fprintf(stdout, "%s:\n", results.name.c_str());
    fprintf(stdout, "Parsing time: %.1f us (%.0f ms),  ", 
            results.parsing_time_ns/1e3, results.parsing_time_ns/1e6);
    fprintf(stdout, "stats: rows=%d cols=%d nnz=%d workset~%.2f MB\n",
            m.n_rows, m.n_cols, m.nnz, results.working_set_bytes / 1e6);
    fprintf(stdout, "COO SpMV: - serial: (%.1f us) (%.1f GB/s)\n", results.coo_time_ns/1000, results.coo_bw_effective_gbs);
    fprintf(stdout, "          - parallel: (%.1f us) (%.1f GB/s)\n", results.coo_parallel_time_ns/1000, results.coo_parallel_bw_effective_gbs);
    fprintf(stdout, "          - GPU basic: (%.1f us) (%.1f GB/s)\n", results.coo_cuda_basic_time_ns/1000, results.coo_cuda_basic_bw_gbs);
    fprintf(stdout, "          - GPU custom 1: (%.1f us) (%.1f GB/s)\n", results.coo_cuda_custom_1_time_ns/1000, results.coo_cuda_custom_1_bw_gbs);
    fprintf(stdout, "          - GPU cuSParse: (%.1f us) (%.1f GB/s)\n", results.coo_cuda_cuSparse_time_ns/1000, results.coo_cuda_cuSparse_bw_gbs);
    fprintf(stdout, "COO to CSR conversion time: %.1f us\n", results.coo_to_csr_time_ns/1000);
    fprintf(stdout, "CSR SpMV: - serial: (%.1f us) (%.1f GB/s)\n", results.csr_time_ns/1000, results.csr_bw_effective_gbs);
    fprintf(stdout, "          - parallel: (%.1f us) (%.1f GB/s)\n", results.csr_parallel_time_ns/1000, results.csr_parallel_bw_effective_gbs);
    fprintf(stdout, "          - GPU basic: (%.1f us) (%.1f GB/s)\n", results.csr_cuda_basic_time_ns/1000, results.csr_cuda_basic_bw_gbs);
    fprintf(stdout, "          - GPU custom 1: (%.1f us) (%.1f GB/s)\n", results.csr_cuda_custom_1_time_ns/1000, results.csr_cuda_custom_1_bw_gbs);
    fprintf(stdout, "          - GPU custom 2: (%.1f us) (%.1f GB/s) [%.1f us setup time]\n", results.csr_cuda_custom_2_time_ns/1000, results.csr_cuda_custom_2_bw_gbs, results.csr_custom_2_setup_time_us);
    fprintf(stdout, "          - GPU cuSParse: (%.1f us) (%.1f GB/s)\n", results.csr_cuda_cuSparse_time_ns/1000, results.csr_cuda_cuSparse_bw_gbs);
    fprintf(stdout, "Results check: ");
    results_comparing(m.n_rows, "COO", res_vector.COO_results, "COO_parallel", res_vector.COO_parallel_results);
    results_comparing(m.n_rows, "CSR", res_vector.CSR_results, "CSR_parallel", res_vector.CSR_parallel_results);
    results_comparing(m.n_rows, "COO", res_vector.COO_results, "CSR\n", res_vector.CSR_results);
    results_comparing(m.n_rows, "COO", res_vector.COO_results, "COO_cuda_basic", res_vector.COO_results_cuda_basic);
    results_comparing(m.n_rows, "COO", res_vector.COO_results, "COO_cuda_custom_1", res_vector.COO_results_cuda_custom_1);
    results_comparing(m.n_rows, "COO", res_vector.COO_results, "COO_cuda_cuSparse", res_vector.COO_results_cuda_cuSparse);
    results_comparing(m.n_rows, "CSR", res_vector.CSR_results, "CSR_cuda_basic", res_vector.CSR_results_cuda_basic);
    results_comparing(m.n_rows, "CSR", res_vector.CSR_results, "CSR_cuda_custom_1", res_vector.CSR_results_cuda_custom_1);
    results_comparing(m.n_rows, "CSR", res_vector.CSR_results, "CSR_cuda_custom_2", res_vector.CSR_results_cuda_custom_2);
    results_comparing(m.n_rows, "CSR", res_vector.CSR_results, "CSR_cuda_cuSparse", res_vector.CSR_results_cuda_cuSparse);
}






