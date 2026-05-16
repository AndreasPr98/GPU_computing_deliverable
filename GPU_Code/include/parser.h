#ifndef PARSER_H
#define PARSER_H

#include <vector>
#include <filesystem>
#include <string.h>

struct COO {
    int n_rows;
    int n_cols;
    int nnz;           
    
    int* row_indices;
    int* col_indices;
    float* values;
};

struct CSR {
    int n_rows;
    int n_cols;
    int nnz;           
    
    int* row_indices;
    int* col_indices;
    float* values;
};

struct MatrixEntry {
    int r, c;
    float v;
};


COO parser(const std::string& file_path, double* parsing_time);

#endif
