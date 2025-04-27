#include "CellCDecon.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <vector>
#include <iostream>
#include <string>
#include <chrono>

int parseArgs(int argc, char* argv[], std::string &filename, int &k, int &seed, int &maxProbes, int &nsamp, int &colskip, float &gamma) {
  for(int i = 1; i < argc; i+=2) {
    if(i + 1 >= argc) {
      return -1; // No value for the last parameter
    }
    
    if(!std::strcmp(argv[i], "-f")) {
      filename = std::string(argv[i+1]);
    }
    else if(!std::strcmp(argv[i], "-k")) {
      k = std::atoi(argv[i+1]);
    }
    else if(!std::strcmp(argv[i], "-s")) {
      seed = std::atoi(argv[i+1]);
    }
    else if(!std::strcmp(argv[i], "-m")) {
      maxProbes = std::atoi(argv[i+1]);
    }
    else if(!std::strcmp(argv[i], "-n")) {
      nsamp = std::atoi(argv[i+1]);
    }
    else if(!std::strcmp(argv[i], "-c")) {
      colskip = std::atoi(argv[i+1]);
    }
    else if (!std::strcmp(argv[i], "-g")) {
      gamma = std::atof(argv[i+1]);
    }
  }
  
  if(k < 1 || nsamp < 1 || filename.empty()) {
    return -1;
  }
  
  return 1;
}

int main(int argc, char *argv[]) {
  // Parse command line arguments
  std::string filename = "";
  int k = 0;
  int seed = std::time(nullptr);
  int maxProbes = 5e5;
  int nsamp = 0;
  int colskip = 1;
  float gamma = 0.0f;
  
  if(parseArgs(argc, argv, filename, k, seed, maxProbes, nsamp, colskip, gamma) < 1) {
    std::fprintf(stderr, "Requires arguments -k <number of cell types> -n <number of samples> -f <input filename>\n");
    std::exit(1);
  }

  // Seed random number generator
  std::srand(seed);

  // Data containers
  std::vector<std::string> sample_ids(nsamp);
  std::vector<std::string> probe_ids;
  std::vector<std::vector<float>> observations;
  std::vector<float> obs_mean;
  std::vector<float> obs_var;
  std::vector<float> obs_min;
  std::vector<float> obs_max;
  
  // Read input file
  std::string sample_prefix;
  int nprobe = 0;
  
  // Read input file
  if (!CellCDeconIO::processFile(filename, nsamp, colskip, sample_prefix, 
                                sample_ids, probe_ids, observations, 
                                obs_mean, obs_var, obs_min, obs_max, 
                                nprobe, maxProbes)) {
    std::cerr << "Failed to process input file: " << filename << std::endl;
    std::exit(1);
  }
  
  std::cout << "Read " << observations.size() << " samples with " << nprobe << " probes." << std::endl;
  
  // Configure the deconvolution algorithm
  CellCDecon::Deconvolution::Config config;
  config.max_unchanged = 10;
  config.max_unconsidered = 100;
  config.max_unchanged_ind = 20;
  config.max_unconsidered_ind = 75;
  
  // Create and initialize the deconvolution object
  CellCDecon::Deconvolution deconvolution(observations, k, gamma, config);
  
  float original_gamma = gamma;
  
  // Timing info
  auto start_time = std::chrono::high_resolution_clock::now();
  
  // Run first 500 iterations with gamma, then 500 with gamma=0
  std::cout << "Running first 500 iterations with gamma = " << gamma << std::endl;
  deconvolution.runIterations(500);
  
  // Set gamma to 0 for the second half
  std::cout << "Running next 500 iterations with gamma = 0" << std::endl;
  deconvolution.setGamma(0.0f);
  float final_likelihood = deconvolution.runIterations(500);
  
  // Timing info
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
  
  std::cout << "Deconvolution completed in " << duration.count() / 1000.0 
            << " seconds with final log likelihood: " << final_likelihood << std::endl;
  
  // Get results
  const auto& weights = deconvolution.getWeights();
  const auto& means = deconvolution.getMeans();
  const auto& vars = deconvolution.getVars();
  
  // Write results to files
  CellCDeconIO::writeFiles(filename, k, seed, weights, nsamp, nprobe, 
                          probe_ids, observations, means, vars, 
                          obs_mean, obs_var, sample_prefix, sample_ids, original_gamma);
  
  std::cout << "Results written to files with prefix: " << filename 
            << ".k" << k << ".seed" << seed << ".gamma" << original_gamma << ".SUMSQ" << std::endl;
  
  return 0;
}


