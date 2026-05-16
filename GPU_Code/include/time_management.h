#ifndef TIME_MANAGEMENT_H
#define TIME_MANAGEMENT_H

#include <sys/time.h>
#include <stdio.h>
#include <filesystem>
#include "utils.h"
#include "parser.h"

#define TIMER_DEF(n)	 struct timespec time_1_##n={0,0}, time_2_##n={0,0}
#define TIMER_START(n)	 clock_gettime(CLOCK_MONOTONIC, &time_1_##n)
#define TIMER_STOP(n)	 clock_gettime(CLOCK_MONOTONIC, &time_2_##n)
#define TIMER_ELAPSED(n) ((time_2_##n.tv_sec-time_1_##n.tv_sec)*1.e9+(time_2_##n.tv_nsec-time_1_##n.tv_nsec))
#define TIMER_PRINT(n) \
    do { \
        printf("Timer elapsed: %.9f\n", TIMER_ELAPSED(n)/1e9);\
        fflush(stdout);\
    } while (0);

void COO_data_printing(const double* runs_time, int niter, 
    const std::filesystem::directory_entry& entry, 
    const COO& COO_matrix);
void CSR_data_printing(const double* runs_time, int niter, 
    const std::filesystem::directory_entry& entry, 
    const CSR& CSR_matrix);

#endif
