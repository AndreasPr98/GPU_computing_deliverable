#include "time_management.h"
#include "utils.h"
#include "parser.h"

static double coo_dram_lower_bound_bytes(const COO& m) {
    // Minimal compulsory traffic approximation (best-case cache reuse):
    // COO arrays once + x once + y once.
    return (double)m.nnz * ((2.0 * sizeof(int)) + (1.0 * sizeof(float))) +
           (double)m.n_cols * sizeof(float) +
           (double)m.n_rows * sizeof(float);
}

static double csr_dram_lower_bound_bytes(const CSR& m) {
    // Minimal compulsory traffic approximation (best-case cache reuse):
    // row ptr once + col/val once + x once + y once.
    return (double)(m.n_rows + 1) * sizeof(int) +
           (double)m.nnz * (sizeof(int) + sizeof(float)) +
           (double)m.n_cols * sizeof(float) +
           (double)m.n_rows * sizeof(float);
}

void COO_data_printing(const double* runs_time, int niter, const std::filesystem::directory_entry& entry, const COO& COO_matrix){
    
    // Print the COO runtime
    double COO_run_time = arithmetic_mean(runs_time, niter);
    fprintf(stdout, "COO matrix (%s) time: %.1f us", 
        entry.path().filename().c_str(), COO_run_time/1000);
    
    const double coo_bw_dram_lb = ns_to_gbs(coo_dram_lower_bound_bytes(COO_matrix), COO_run_time);
    const double coo_bw_effective = ns_to_gbs(coo_effective_bytes(COO_matrix), COO_run_time);
    fprintf(stdout, "      COO BW_dram_lb: %.4f GB/s  BW_effective: %.4f GB/s",
        coo_bw_dram_lb, coo_bw_effective);

}

void CSR_data_printing(const double* runs_time, int niter, const std::filesystem::directory_entry& entry, const CSR& CSR_matrix){

    // Print the CSR runtime
    double CSR_run_time = arithmetic_mean(runs_time, niter);
    fprintf(stdout, "CSR matrix (%s) time: %.1f us", entry.path().filename().c_str(), CSR_run_time/1000);
    
    const double csr_bw_dram_lb = ns_to_gbs(csr_dram_lower_bound_bytes(CSR_matrix), CSR_run_time);
    const double csr_bw_effective = ns_to_gbs(csr_effective_bytes(CSR_matrix), CSR_run_time);
    fprintf(stdout, "      CSR BW_dram_lb: %.4f GB/s  BW_effective: %.4f GB/s",
        csr_bw_dram_lb, csr_bw_effective); 
}
