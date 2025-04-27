# CellCDecon

CellCDecon is a program that attempts to find tissue cell composition effects from high throughput experiments in complex tissues such as whole blood or adipose tissue. It was designed for Illumina methylation beta values ranging from 0 to 1 but potentially has other applications in gene expression or other high throughput studies as well.

*Last updated: November 6, 2014*  
*Author: James Wagner (jamesrwagner@gmail.com)*

## Project Overview

CellCDecon takes as input a matrix of values and a parameter k (number of cell types), and outputs three files that describe:
1. Cell type proportions for each sample
2. Cell type-specific means and variances for each probe
3. Residuals after accounting for cell composition

## Project Structure

The project has been reorganized to follow modern C++ conventions:

```
CellCDecon/
├── include/               # Header files
├── src/                   # Source files
├── test/                  # Test files
│   ├── unit/              # Unit tests
│   ├── integration/       # Integration tests (TODO)
│   └── data/              # Test data
├── data/                  # Data directories
│   ├── synthetic_data/    # Synthetic datasets
│   └── evaluation_results/# Evaluation output
└── tools/                 # Supporting tools
    ├── analysis/          # Analysis scripts
    ├── data_generation/   # Test data generation
    └── scripts/           # Utility scripts
```

## Input File Format

The input file is tab or space separated with each row corresponding to one probe and each column to one sample/individual in the study. A header file is required listing the sample ids, and at least one column with probe ids.

- Probe ids and sample ids can be any valid string
- Missing values should be indicated by "-1"
- All non-missing measurement values should be non-negative

**Example input file** for a study with 3 samples and 2 probes:

```
id Samp_1 Samp_2 Samp_3
cg00001 0 0.2 0.3
cg00002 0.7 0.8 0.9
```

## Building and Running the Project

### Building with CMake (Recommended)

```bash
# Create a build directory
mkdir build && cd build

# Configure and build
cmake ..
make

# Install (optional)
make install
```

### Running the Main Executable

```bash
# If built with CMake
./build/cellcdecon_app -k <k> -n <nsamp> -f <filename> [-c <columns>] [-m <maxProbes>] [-s <seed>] [-g <gamma>]

```

#### Required Arguments

- `-k <k>`: An integer ≥ 1 specifying the number of cell types to infer in the algorithm. In practice, it is recommended to try several values of k and observe results.
- `-n <nsamp>`: An integer ≥ 1 specifying the number of samples/individuals in your study (i.e., the number of columns in your input file that are not probe identifiers or other probe metadata).
- `-f <filename>`: An absolute or relative path to your input file.

#### Optional Arguments

- `-c <columns>`: An integer ≥ 1 specifying how many columns at the beginning of each row correspond to probe meta-data. Defaults to 1 if not specified.
- `-m <maxProbes>`: An integer ≥ 1 specifying how many probes (i.e., rows excluding the header row) are to be read in. Default: 500,000.
- `-s <seed>`: An integer seed used for random number generation. Default: current time.
- `-g <gamma>`: Regularization parameter. Default: 0.0.

## Running Evaluation Scripts

The project includes Python scripts for data generation and evaluation. These scripts are located in the `tools` directory:

### Generating Test Data

```bash
# Generate synthetic test data
python tools/data_generation/generate_test_data.py [options]
```

### Evaluating Results

```bash
# Evaluate deconvolution results
python tools/analysis/evaluate_results.py [options]
```

### Simulation Scripts

```bash
# Run simulations
./tools/scripts/run_simulation.sh
```

## Running Tests

### Running Unit Tests with CMake

```bash
# From the build directory
make test

# Or to run all tests
ctest
```

### Running Unit Tests Directly

```bash
cd test
g++ -isystem ${GTEST_DIR}/include -I${GTEST_DIR} -pthread -c ${GTEST_DIR}/src/gtest-all.cc
ar -rv libgtest.a gtest-all.o
g++ -g -isystem ${GTEST_DIR}/include -pthread unit/test_cellc.cpp ../src/CellCDecon.cpp libgtest.a -o test_cellc
./test_cellc
```

## Output Files

Three output files of the form: `<input_file>.k<k>.seed<seed>.<suffix>` will appear in the same directory as the input file:

1. **`.w`**: A weight file indicating the composition of each sample. Each row is one sample and each column is one inferred "cell type". Values correspond to the proportion of that individual's cells that belong to a given cell type.

2. **`.meanvar`**: A file with the mean and variance of each probe, for each cell type. Each row is one probe.
   - The first column is the probe id
   - If a value greater than 1 was given for the `-c` parameter, these additional columns will appear at the start of each row
   - The following columns alternate between the mean and variance inferred for each of the k cell types
   - The final two columns correspond to the mean and variance *observed* for that probe in the experiment

3. **`.resid`**: This file contains residuals obtained by subtracting the weighted sum of means for a given probe from the observation. The weighting is done from each sample based on the cell types.
   
   **Example**: If for k = 2, an individual i has an observed methylation value for probe j of 0.5, individual i has cell composition weights of 0.75 and 0.25, and probe j has inferred mean methylation values for the two cell types of 0.8 and 0.4, then the residual will be 0.5 - (0.75*0.8 + 0.25*0.4) = -0.2. This indicates that this individual may have had some sort of effect in one or more of its constituent cell types that led to a relatively hypo-methylated phenotype at this probe, after taking into account his or her cell composition.

## Algorithm Design and Workflow

- `src/main.cpp` contains the main function which parses the command line arguments, allocates memory for data structures, and calls other functions.
- `initialize_weight()` and `initialize_meanvar()` start with random assignments for the weights, mean and variance that are to be optimized.
- `update_weight()` and `update_meanvar()` are called iteratively (1000 times). Each individual's weight vector and each probe's mean and variance vectors are perturbed randomly, with changes accepted if they result in a higher log likelihood of the observed experimental values.

## Future Work

Cell composition deconvolution is a rich area of future work and research given the plethora of high throughput platforms being developed to measure values such as DNA methylation in multiple sites in complex tissue obtained from multiple individuals.

Possible work to be done include:
1. Integrating prior knowledge about cell compositions or methylation/expression levels in pure cell types to render the algorithm more semi-supervised
2. Automatic determination of k
3. Parallelization of inference
4. Replacing random perturbations with gradient-based approach
5. Adding integration tests for end-to-end workflow validation

## License

See the [LICENSE](License.txt) file for details.

## Citation

If you use this software in your research, please cite:

```
Wagner, J.R. (2014). CellCDecon: A tool for cell type deconvolution from epigenomic data.
``` 