# **SpMV Cache Analysis Tool**

This project provides an automated pipeline for analyzing the performance of Sparse Matrix-Vector Multiplication (SpMV) using Coordinate (COO) and Compressed Sparse Row (CSR) formats. It is designed to run on Linux/HPC systems and uses Valgrind Cachegrind for cache profiling.

## Requirements

- A C++ compiler with C++17 support (e.g. `g++` >= 7).  
- `make` (GNU Make).  
- `valgrind` and `cg_annotate` available in `PATH`.  
- `python3` (for `valgrind_to_csv.py`).  
- `bash` (for `generate_report.sh` and cluster scripts).  
- An HPC scheduler if you plan to run batch jobs (PBS or Slurm).
- The time for all matrixes completion is around 1.5-2 hours so be sure to have enough server time

> Note: The project is intended for Linux/HPC environments. The provided `Makefile` and shell scripts assume a Unix-like shell and POSIX tools — native Windows is not supported without modifications.

## 1. Environment Setup (Path Configuration)

The C++ source code uses an absolute `base_path` for file management. Before building or running on a different machine, update the `base_path` variable in `cache_miss_SpMV.cpp` to point to your project/matrices directory.

**Action**: Open [cache_miss_SpMV.cpp](cache_miss_SpMV.cpp) and update the `base_path` definition near the top of `main()`:

```cpp
// Change base path with the path where you have the code and the matrices folder
const std::filesystem::path base_path = "/home/youruser/GPU_computing";
const std::filesystem::path target_dir = base_path / "Matrixes";
const std::filesystem::path matrix_numbers_path = base_path / "Cache_miss_Code" / "program_output" / "matrix_numbers.txt";
```

## 2. Compilation

From the project root run:

```bash
cd Cache_miss_Code
make
```

* **Binary**: This creates the executable `bin/cache_miss_SpMV` (the Makefile builds `cache_miss_SpMV`).
* **Folders**: The build process automatically ensures `obj/`, `bin/`, `valgrind_output/`, and `program_output/` exist.

If your compiler doesn't default to C++17 you may need to add `-std=c++17` to `CXXFLAGS` in the [makefile](makefile).

## 3. Running (local)

To run the program on a single matrix (matrix number based on alphabetical order inside `Matrixes`):

```bash
./bin/cache_miss_SpMV 1
```

Profiling a single run with Cachegrind manually:

```bash
valgrind --tool=cachegrind --cache-sim=yes --cachegrind-out-file=valgrind_output/cache_miss.m1 --branch-sim=yes ./bin/cache_miss_SpMV 1
```

## 4. Submitting Profiling Jobs (batch)

The repository includes a PBS job script [cachegrind.pbs](cachegrind.pbs). A Slurm equivalent is provided as [cachegrind.slurm](cachegrind.slurm).

If you run [SpMV.slurm](SpMV.slurm) instead of Cachegrind, the code appends measured performance data to `program_output/spmv_results.csv`, including actual runtime and effective bandwidth for each SpMV stage/function (COO/CSR serial and parallel, plus parsing/conversion timings).

**PBS:**

```bash
qsub cachegrind.pbs
```

**Slurm:**

```bash
sbatch cachegrind.slurm
```

Both scripts will produce `valgrind_output/cache_miss.mX` files for each matrix ID.

## 5. Post-Processing & Reporting

Once the jobs are finished and the raw `cache_miss.mX` files are present in `valgrind_output/`, execute the reporting pipeline:

```bash
make report
```

This command runs `cg_annotate` (via `generate_report.sh`) to produce `valgrind_output/annotate_miss.mX` and then runs `valgrind_to_csv.py` to generate the final `spmv_cache_analysis.csv`.

## 6. Evaluating Results

Analyze the generated `spmv_cache_analysis.csv` to compare the two formats.

- **L1 Misses (D1mr)**: Measures primary cache efficiency.  
- **LLC Misses (DLmr)**: Measures high-latency memory traffic.

**Notes:**

- If a matrix fits entirely in the server's L3 cache (small matrices), the LLC misses can be zero.  
- For larger datasets the `CSR` format often reduces cache misses vs `COO`.

## Summary of Makefile Targets

| Command | Description |
| :---- | :---- |
| `make` | Compiles the SpMV program. |
| `make report` | Generates text annotations and the CSV analysis. |
| `make clean` | Removes compiled binaries and objects. |
| `make clean_reports` | Deletes generated annotation files. |

