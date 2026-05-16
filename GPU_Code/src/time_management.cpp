#include "time_management.h"
#include "utils.h"
#include "parser.h"

void COO_data_printing(const double* runs_time, int niter, const std::filesystem::directory_entry& entry, const COO& COO_matrix){
    
    // Print the COO runtime
    double COO_run_time = arithmetic_mean(runs_time, niter);
    fprintf(stdout, "COO matrix (%s) time: %.1f us", 
        entry.path().filename().c_str(), COO_run_time/1000);
    
    const double coo_gflops = ns_to_gflops(spmv_flop_count(COO_matrix.nnz), COO_run_time);
    fprintf(stdout, "      COO throughput: %.4f GFLOP/s", coo_gflops);

}

void CSR_data_printing(const double* runs_time, int niter, const std::filesystem::directory_entry& entry, const CSR& CSR_matrix){

    // Print the CSR runtime
    double CSR_run_time = arithmetic_mean(runs_time, niter);
    fprintf(stdout, "CSR matrix (%s) time: %.1f us", entry.path().filename().c_str(), CSR_run_time/1000);
    
    const double csr_gflops = ns_to_gflops(spmv_flop_count(CSR_matrix.nnz), CSR_run_time);
    fprintf(stdout, "      CSR throughput: %.4f GFLOP/s", csr_gflops); 
}
