README for CellCDecon program, written by James Wagner, jamesrwagner@gmail.com. Last updated November 6, 2014

CellCDecon is a program that attempts to find tissue cell composition effects from high throughput experiments in complex tissues such as whole blood or adipose tissue.
It was designed for Illumina methylation beta values ranging from 0 to 1 but potentially has other applications in gene expression or other high throughput studies as well.
It takes as input a matrix of values and a parameter k, the output are three files described below. 


1) Input file
The input file is tab or space separated with each row corresponding to one probe and each column to one sample/individual in the study. A header file is required
listing the sample ids, and at least one column with probe ids. Probe ids and sample ids can be any valid string. Missing values should be indicated by "-1", and all non-missing measurement values should be non-negative.
Here is a sample input file for a study with 3 samples and 2 probes:

id Samp_1 Samp_2 Samp_3
cg00001 0 0.2 0.3
cg00002 0.7 0.8 0.9

2) Building
Download Makefile, the two .cpp files and the .h file into the same directory. Building requires make and g++ and makes use of no other specialized libraries, type "make" on a Unix command line to build the program. There should
now be an executable called "CellCDecon" in the same directory.

3) Running
CellCDecon requires at this time, three arguments, with three additional optional arguments: The basic form is ./CellCDecon -k <k> -n <nsamp> -f <filename>
a) k is an integer >= 1 specifying the number of cell types to infer in the algorithm. In practice working with real world datasets it is recommended to try several values of k and observe results.
b) nsamp is an integer >= 1 specifying the number of samples/individuals in your study (i.e. the number of columns in your input file that are not probe identifiers or other probe metadata)
c) filename is an absolute or relative path to your input file
Between 0 and 3 additional optional arguments that can be added if needed or wanted are specified as -c <columns>  -m <maxProbes> -s <seed>
d) columns is an integer >= 1 specifying how many columns at the beginning of each row correspond to probe meta-data (for example if each row starts with a "<probe_id> <chromosome> <genomic_coordinate>", columns would be 3)
defaults to 1 if not specified
e) maxProbes is an integer >= 1 specifying how many probes (i.e. rows excluding the header row) are to be read in. The program will stop reading in input either when maxProbes rows have been read, 
or it reaches the end of the input file, whichever comes first. Default: 500,000
f) seed is an integer (non-negative)? that is used in srand as at various points the algorithm generates random numbers which may lead to  different outputs with different seeds. Default: current time 

Arguments can be specified in any order, however a space must separate the parameter name and the parameter value (for example "-k 3" and NOT "-k3").

4) Output
Three output files of the form: <input_file>.k<k>.seed<seed>.<suffix> will appear in the same directory as the input file. <input_file>, <k> and <seed> were user-specified parameters, the suffixes are as follows
1: w. This a weight file indicating the composition of each sample. Each row is one sample and each column is one inferred "cell type". A given value will therefore correspond to the 
proportion of that individual's cell's that correspond to a given cell type.
2: meanvar. This is a file with the mean and variance of each probe, for each cell type. Each row will be one probe. The first column will be the probe id. If a value greater than 1
was given for the "-c" parameter, these additional columns will appear at the start of each row. The following columns will alternate between the mean and variance inferred for each of the k cell types.
The final two columns correspond to the mean and variance *observed* for that probe in the experiment.
3: .resid This is a file obtained by subtracting the weighted sum of means for a given probe from the observation. The weighting is done from each sample based on the cell types.
For example if for k = 2, an individual i has an observed methylation value for probe j of 0.5, individual i has cell composition weights of 0.75 and 0.25, and probe j has inferred mean methylation values
for the two cell types of 0.8 and 0.4, then the residual will be 0.5 - (0.75*0.8 + 0.25*0.4) = -0.2. This indicates that this individual may have had some sort of effect in one or more of its constituent 
cell types that led to a relatively hypo-methylated phenotype at this probe, after taking into account his or her cell composition.
Each row in the resid file will correspond to one probe. The first column(s) will correspond to the probe metadata in the original input files (as specified by the -c argument if greater than 1)
followed by the residual values. There is also a header line corresponding to the header line of the original input file. Note that as of now output file are all space-separated. 

 

5) Algorithm design and work flow

CellCMain.cpp contains the main function which will parse the command line arguments, allocate the memory necessary for the data structures used in the program, and then call other functions.
initliaze_weight() and initialize_meanvar() start with random assignments for the weights, mean and variance that are to be optimized by this approach. Mean and variance are determined based on some
properties of the observed means for a given probe.
update_weight() and update_meanvar() are then called iteratively (1000 times at the time of writing this README). Each individuals weight vector, and each probes mean and variance vectors will be 
perturbed in a random fashion, and the changes accepted if they result in an overall higher log likelihood of the observed experimental values. 

6) Unit testing 
unit test files are located in the test directory and were designed using gtest (https://code.google.com/p/googletest/) system. Please see gtest instructions for more details on how to build test files, 
however the following usage taken from gtest's README was sufficient to build a test executable :

cd test

g++ -isystem ${GTEST_DIR}/include -I${GTEST_DIR} -pthread -c ${GTEST_DIR}/src/gtest-all.cc
ar -rv libgtest.a gtest-all.o
g++ -g -isystem ${GTEST_DIR}/include -pthread test_cellc.cc ../CellCDecon.cpp libgtest.a -o out

7) Future work

Cell composition deconvolution is a rich area of future work and research given the plethora of high throughput platforms being developed to measure
values such as DNA methylation in multiple sites in complex tissue obtained from multiple individuals.
Possible work to be done include:
a) integrating prior knowledge about cell compositions or methylation/expression levels in pure cell types to render the algorithm more semi-supervised
b) Automatic determination of k
c) parallelization of inference
d) replacing random perturbations with gradient-based approach
