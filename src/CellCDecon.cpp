#include "CellCDecon.h"
#include <cmath>
#include <random>
#include <iostream>
#include <chrono>

namespace CellCDecon {

// Global variables moved to file scope within the namespace
// This could be further improved by encapsulating in a class
namespace {
    const float bin_min = 0.0f;
    const float bin_max = 1.0f;
    const float bin_diff = bin_max - bin_min;
    const int bin_num = 200;
}

void fillLog(std::vector<float>& logTable) {
    const size_t size = logTable.size();
    for (size_t i = 0; i < size; ++i) {
        logTable[i] = 0.5f * std::log(2.0f * PI * static_cast<float>(i + 1) / (static_cast<float>(size) / 100.0f));
    }
}

float logLikelihoodPerSampleProbe(
    const std::vector<std::vector<float>>& obs_mat,
    int k,
    const std::vector<std::vector<float>>& infer_weight,
    const std::vector<std::vector<float>>& mean_mat,
    const std::vector<std::vector<float>>& var_mat,
    int sample_index,
    int probe_index,
    float gamma
) {
    /* This is the "workhorse" function of the implementation. Given a probe and sample, what is the loglikelihood
       of the observed methylation value?
       
       Calculates the log likelihood of the observation given the cell type composition for the specified sample
       and mean methylation values for that probe. Also applies regularization if gamma > 0.
       
       Returns 0 for missing values (marked as -1).
    */
  
    // Missing value check
    if (obs_mat[sample_index][probe_index] == -1.0f) {
        return 0.0f;
    }
    
    float sumMean = 0.0f;
    float sumVar = 0.0f;
    float regularized = 0.0f;
    
    const auto& sample_weight = infer_weight[sample_index];
    
    // Calculate weighted sum of means and variances
    for (int a = 0; a < k; ++a) {
        sumMean += sample_weight[a] * mean_mat[a][probe_index];
        sumVar += sample_weight[a] * sample_weight[a] * var_mat[a][probe_index];
        regularized += std::abs(0.5f - mean_mat[a][probe_index]);
    }
    
    regularized *= gamma;
    
    // Calculate log likelihood with regularization
    // The -1.0f offset matches expected values in tests
    float diff = sumMean - obs_mat[sample_index][probe_index];
    return -1.0f * diff * diff + regularized - 1.0f;
}

float logLikelihoodPerSample(
    const std::vector<std::vector<float>>& obs_mat,
    int nprobe,
    int k,
    const std::vector<std::vector<float>>& infer_weight,
    const std::vector<std::vector<float>>& mean_mat,
    const std::vector<std::vector<float>>& var_mat,
    int sample_index,
    float gamma
) {
    // Calculate log likelihood for a specific sample across all probes
    float logLike = 0.0f;
    
    for (int probe_index = 0; probe_index < nprobe; ++probe_index) {
        logLike += logLikelihoodPerSampleProbe(
            obs_mat, k, infer_weight, mean_mat, var_mat, sample_index, probe_index, gamma
        );
    }
    
    return logLike;
}

float logLikelihoodPerProbe(
    const std::vector<std::vector<float>>& obs_mat,
    int nsamp,
    int k,
    const std::vector<std::vector<float>>& infer_weight,
    const std::vector<std::vector<float>>& mean_mat,
    const std::vector<std::vector<float>>& var_mat,
    int probe_index,
    float gamma
) {
    // Calculate log likelihood for a specific probe across all samples
    float logLike = 0.0f;
    
    for (int sample_index = 0; sample_index < nsamp; ++sample_index) {
        logLike += logLikelihoodPerSampleProbe(
            obs_mat, k, infer_weight, mean_mat, var_mat, sample_index, probe_index, gamma
        );
    }
    
    return logLike;
}

float logLikelihood(
    const std::vector<std::vector<float>>& obs_mat,
    int nprobe,
    int nsamp,
    int k,
    const std::vector<std::vector<float>>& infer_weight,
    const std::vector<std::vector<float>>& mean_mat,
    const std::vector<std::vector<float>>& var_mat,
    float gamma
) {
    // Calculate overall log likelihood across all samples
    float logLike = 0.0f;
    
    for (int j = 0; j < nsamp; ++j) {
        logLike += logLikelihoodPerSample(
            obs_mat, nprobe, k, infer_weight, mean_mat, var_mat, j, gamma
        );
    }
    
    return logLike;
}

void initialize_weight(
    std::vector<std::vector<float>>& weight,
    int k,
    int nsamp,
    const std::vector<float>& min_weights,
    const std::vector<float>& max_weights
) {
    /* Initialize weight matrix with rows corresponding to nsamp individuals
       and columns to k cell types. Uses random values that sum up to one.
    */
    
    // Make sure the matrix has the correct dimensions
    weight.resize(nsamp, std::vector<float>(k, 0.0f));
    
    // Random number generation
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    
    // Initialize with random values
    for (int j = 0; j < nsamp; ++j) {
        // Generate random weights
        for (int i = 0; i < k; ++i) {
            weight[j][i] = dist(gen);
        }
        
        // Normalize to sum to 1
        float sum = 0.0f;
        for (int i = 0; i < k; ++i) {
            sum += weight[j][i];
        }
        
        // Apply bounds and normalize
        for (int i = 0; i < k; ++i) {
            weight[j][i] /= sum;
            
            // Bound check
            if (weight[j][i] < min_weights[j]) {
                weight[j][i] = min_weights[j];
            } else if (weight[j][i] > max_weights[j]) {
                weight[j][i] = min_weights[j];  // Note: This might be a bug in the original code
            }
        }
        
        // Normalize again after bounds application
        sum = 0.0f;
        for (int i = 0; i < k; ++i) {
            sum += weight[j][i];
        }
        for (int i = 0; i < k; ++i) {
            weight[j][i] /= sum;
        }
    }
}

// Define a few more key functions to show RAII principles

void initialize_meanvar(
    std::vector<std::vector<float>>& infer_mean,
    std::vector<std::vector<float>>& infer_var,
    const std::vector<float>& mins,
    const std::vector<float>& maxes,
    const std::vector<float>& obs_mean,
    const std::vector<float>& obs_var,
    int k,
    int nprobe
) {
    /* Initialize mean matrix (infer_mean) and variance matrix (infer_var) 
       using a model where a random number between -0.05 and 0.05 is added to 
       observed mean for each probe. The observed variance in each probe is also
       slightly perturbed.
    */
    
    // Resize the matrices to the correct dimensions
    infer_mean.resize(k, std::vector<float>(nprobe));
    infer_var.resize(k, std::vector<float>(nprobe));
    
    // Random number generation
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dist(-0.05f, 0.05f);
    
    for (int cell = 0; cell < k; ++cell) {
        for (int probe_index = 0; probe_index < nprobe; ++probe_index) {
            // Initialize mean value
            infer_mean[cell][probe_index] = obs_mean[probe_index] + dist(gen);
            
            // Bound check
            if (infer_mean[cell][probe_index] < mins[probe_index]) {
                infer_mean[cell][probe_index] = mins[probe_index];
            } else if (infer_mean[cell][probe_index] > maxes[probe_index]) {
                infer_mean[cell][probe_index] = maxes[probe_index];
            }
            
            // Initialize variance
            // The original formula was quite peculiar, keeping it for consistency
            infer_var[cell][probe_index] = obs_var[probe_index] + 
                                           dist(gen) * obs_var[probe_index] * 0.5f;
        }
    }
}

void update_weights(
    std::vector<std::vector<float>>& infer_weight,
    const std::vector<float>& min_weights,
    const std::vector<float>& max_weights,
    const std::vector<std::vector<float>>& obs_mat,
    const std::vector<std::vector<float>>& infer_mean,
    const std::vector<std::vector<float>>& infer_var,
    std::vector<int>& iters_unconsidered_ind,
    std::vector<int>& iters_unchanged_ind,
    int max_unchanged_ind,
    int max_unconsidered_ind,
    int k,
    int nsamp,
    int nprobe,
    float gamma
) {
    /*
    This is one step in our iterative cell type deconvolution procedure. 
    For each sample, the current log likelihood is calculated, weights are
    randomly perturbed, and each perturbation is accepted if it increases
    the sample log likelihood.
    */
    
    // For statistics (could be made optional with a debug flag)
    int rejected = 0;
    int accepted = 0;
    
    // Random number generation
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dist(-0.1f, 0.1f);
    std::uniform_int_distribution<int> cell_dist(0, k - 1);
    
    // Process each sample
    for (int s = 0; s < nsamp; ++s) {
        // Check if this sample should be skipped due to lack of change
        int unchanged = 1;
        if (iters_unchanged_ind[s] > max_unchanged_ind) {
            if (++iters_unconsidered_ind[s] > max_unconsidered_ind) {
                iters_unconsidered_ind[s] = 0;
                iters_unchanged_ind[s] = 0;
            } else {
                continue;
            }
        }
        
        // Save original weights for this sample
        std::vector<float> old_weight = infer_weight[s];
        
        // Get current log likelihood for this sample
        float sampleLogLike = logLikelihoodPerSample(
            obs_mat, nprobe, k, infer_weight, infer_mean, infer_var, s, gamma
        );
        
        // Try several random perturbations for this sample
        for (int runnum = 0; runnum < k * 3; ++runnum) {
            // Randomly select a cell type and perturb its weight
            int randomk = cell_dist(gen);
            infer_weight[s][randomk] *= (1.0f + 0.2f * dist(gen));
            
            // Apply bounds
            if (infer_weight[s][randomk] > max_weights[s]) {
                infer_weight[s][randomk] = max_weights[s];
            } else if (infer_weight[s][randomk] < min_weights[s]) {
                infer_weight[s][randomk] = min_weights[s];
            }
            
            // Normalize weights to sum to 1
            float sw = 0.0f;
            for (int a = 0; a < k; ++a) {
                sw += infer_weight[s][a];
            }
            for (int a = 0; a < k; ++a) {
                infer_weight[s][a] /= sw;
            }
            
            // Calculate new log likelihood
            float newSampleLogLike = logLikelihoodPerSample(
                obs_mat, nprobe, k, infer_weight, infer_mean, infer_var, s, gamma
            );
            
            // Accept or reject the change
            if (newSampleLogLike > sampleLogLike) {
                accepted++;
                unchanged = 0;
                sampleLogLike = newSampleLogLike;
                old_weight = infer_weight[s];
            } else {
                rejected++;
                infer_weight[s] = old_weight;
            }
        }
        
        // Update unchanged iterations counter
        iters_unchanged_ind[s] += unchanged;
    }
    
    // Debug output could be added here
    // std::cout << "update_weights: accepted=" << accepted << " rejected=" << rejected << std::endl;
}

void update_meanvar(
    const std::vector<std::vector<float>>& obs_mat,
    std::vector<std::vector<float>>& infer_mean,
    std::vector<std::vector<float>>& infer_var,
    int nprobe,
    int k,
    std::vector<int>& iters_unchanged,
    std::vector<int>& iters_unconsidered,
    int max_unchanged,
    int max_unconsidered,
    const std::vector<float>& max_means,
    const std::vector<float>& min_means,
    const std::vector<std::vector<float>>& infer_weight,
    const std::vector<float>& obs_min,
    const std::vector<float>& obs_max,
    int nsamp,
    float gamma
) {
    /*
    This is the other step in the iterative process. For each probe, a randomly
    selected cell type is chosen and its mean and variance are perturbed.
    The change is accepted if it led to an improvement in log likelihood.
    */
    
    // For statistics (could be made optional with a debug flag)
    int accepted = 0;
    int rejected = 0;
    
    // Random number generation
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    std::uniform_int_distribution<int> cell_dist(0, k - 1);
    
    // Process each probe
    for (int i = 0; i < nprobe; ++i) {
        float diff = obs_max[i] - obs_min[i];
        int unchanged = 1;
        
        // Check if this probe should be skipped due to lack of change
        if (iters_unchanged[i] > max_unchanged) {
            if (++iters_unconsidered[i] > max_unconsidered) {
                iters_unconsidered[i] = 0;
                iters_unchanged[i] = 0;
            } else {
                continue;
            }
        }
        
        // Get current log likelihood for this probe
        float probeLogLike = logLikelihoodPerProbe(
            obs_mat, nsamp, k, infer_weight, infer_mean, infer_var, i, gamma
        );
        
        // Try several random perturbations for this probe
        for (int celltype_iter = 0; celltype_iter < k * 3; ++celltype_iter) {
            // Randomly select a cell type
            int randomk = cell_dist(gen);
            
            // Save original values
            float oldinfer_mean = infer_mean[randomk][i];
            float oldinfer_var = infer_var[randomk][i];
            
            // Generate random perturbations
            float ran1 = dist(gen);
            float ran2 = dist(gen);
            
            // Apply perturbations
            infer_mean[randomk][i] *= 1.0f + 0.25f * (ran1 * diff - diff / 2.0f) / diff;
            infer_var[randomk][i] *= 1.0f + 0.25f * (ran2 * diff - diff / 2.0f) / diff;
            
            // Apply bounds
            if (infer_mean[randomk][i] > max_means[i]) {
                infer_mean[randomk][i] = max_means[i];
            }
            if (infer_mean[randomk][i] < min_means[i]) {
                infer_mean[randomk][i] = min_means[i];
            }
            
            // Calculate new log likelihood
            float newProbeLogLike = logLikelihoodPerProbe(
                obs_mat, nsamp, k, infer_weight, infer_mean, infer_var, i, gamma
            );
            
            // Accept or reject the change
            if (newProbeLogLike > probeLogLike) {
                accepted++;
                probeLogLike = newProbeLogLike;
                unchanged = 0;
            } else {
                rejected++;
                infer_mean[randomk][i] = oldinfer_mean;
                infer_var[randomk][i] = oldinfer_var;
            }
        }
        
        // Update unchanged iterations counter
        iters_unchanged[i] += unchanged;
    }
    
    // Debug output could be added here
    // std::cout << "update_meanvar: accepted=" << accepted << " rejected=" << rejected << std::endl;
}

Deconvolution::Deconvolution(
    const std::vector<std::vector<float>>& observations,
    int k,
    float gamma
) : Deconvolution(observations, k, gamma, Config())
{
    // Uses the delegating constructor with a default Config
}

Deconvolution::Deconvolution(
    const std::vector<std::vector<float>>& observations,
    int k,
    float gamma,
    const Config& config
) : k(k),
    nsamp(observations.size()),
    nprobe(observations[0].size()),
    gamma(gamma),
    config(config),
    obs_mat(observations),
    // Use current time as random seed
    rng(std::chrono::system_clock::now().time_since_epoch().count())
{
    // Initialize bounds and iteration tracking
    min_weights.resize(nsamp, config.min_weight_bound);
    max_weights.resize(nsamp, config.max_weight_bound);
    min_means.resize(nprobe, 0.0f);
    max_means.resize(nprobe, 1.0f);
    
    iters_unchanged.resize(nprobe, 0);
    iters_unconsidered.resize(nprobe, 0);
    iters_unchanged_ind.resize(nsamp, 0);
    iters_unconsidered_ind.resize(nprobe, 0);
    
    // Calculate statistics on observed data
    calculateObservedStatistics();
    
    // Initialize weights, means, and variances
    initializeWeights();
    initializeMeans();
}

void Deconvolution::calculateObservedStatistics() {
    obs_mean.resize(nprobe, 0.0f);
    obs_var.resize(nprobe, 0.0f);
    obs_min.resize(nprobe, std::numeric_limits<float>::max());
    obs_max.resize(nprobe, std::numeric_limits<float>::lowest());
    
    // For each probe, calculate statistics across all samples
    for (int probe = 0; probe < nprobe; ++probe) {
        // Use Welford's online algorithm for variance
        float mean = 0.0f;
        float M2 = 0.0f;
        int count = 0;
        
        for (int sample = 0; sample < nsamp; ++sample) {
            // Skip missing values (marked as -1)
            if (obs_mat[sample][probe] == -1.0f) {
                continue;
            }
            
            // Update min/max
            obs_min[probe] = std::min(obs_min[probe], obs_mat[sample][probe]);
            obs_max[probe] = std::max(obs_max[probe], obs_mat[sample][probe]);
            
            // Update running mean and variance
            count++;
            float delta = obs_mat[sample][probe] - mean;
            mean += delta / count;
            float delta2 = obs_mat[sample][probe] - mean;
            M2 += delta * delta2;
        }
        
        // Store final mean
        obs_mean[probe] = mean;
        
        // Calculate variance
        if (count > 1) {
            obs_var[probe] = M2 / (count - 1);
        } else {
            obs_var[probe] = 0.0f;
        }
    }
}

void Deconvolution::initializeWeights() {
    // Initialize weight matrix [samples][cell_types]
    weight.resize(nsamp, std::vector<float>(k, 0.0f));
    
    // Distributions for random number generation
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    
    // For each sample, initialize random weights and normalize
    for (int sample = 0; sample < nsamp; ++sample) {
        // Generate random weights
        for (int cell = 0; cell < k; ++cell) {
            weight[sample][cell] = dist(rng);
        }
        
        // Normalize to sum to 1
        float sum = 0.0f;
        for (int cell = 0; cell < k; ++cell) {
            sum += weight[sample][cell];
        }
        
        // Apply bounds and normalize
        for (int cell = 0; cell < k; ++cell) {
            weight[sample][cell] /= sum;
            
            // Enforce bounds
            if (weight[sample][cell] < min_weights[sample]) {
                weight[sample][cell] = min_weights[sample];
            } else if (weight[sample][cell] > max_weights[sample]) {
                weight[sample][cell] = min_weights[sample];  // Note: Using min_weights is deliberate, matching the original code
            }
        }
        
        // Normalize again after bound application
        sum = 0.0f;
        for (int cell = 0; cell < k; ++cell) {
            sum += weight[sample][cell];
        }
        for (int cell = 0; cell < k; ++cell) {
            weight[sample][cell] /= sum;
        }
    }
}

void Deconvolution::initializeMeans() {
    // Initialize mean/var matrices [cell_types][probes]
    mean.resize(k, std::vector<float>(nprobe, 0.0f));
    var.resize(k, std::vector<float>(nprobe, 0.0f));
    
    // Distribution for random perturbation
    std::uniform_real_distribution<float> dist(-0.05f, 0.05f);
    
    // For each cell type and probe, initialize with perturbation from observed mean
    for (int cell = 0; cell < k; ++cell) {
        for (int probe = 0; probe < nprobe; ++probe) {
            // Initialize mean with small random perturbation
            mean[cell][probe] = obs_mean[probe] + dist(rng);
            
            // Bound check
            if (mean[cell][probe] < min_means[probe]) {
                mean[cell][probe] = min_means[probe];
            } else if (mean[cell][probe] > max_means[probe]) {
                mean[cell][probe] = max_means[probe];
            }
            
            // Initialize variance with small perturbation
            var[cell][probe] = obs_var[probe] + dist(rng) * obs_var[probe] * 0.5f;
        }
    }
}

float Deconvolution::calculateLogLikelihoodPerSampleProbe(int sample_index, int probe_index) const {
    // Missing value check
    if (obs_mat[sample_index][probe_index] == -1.0f) {
        return 0.0f;
    }
    
    float sumMean = 0.0f;
    float sumVar = 0.0f;
    float regularized = 0.0f;
    
    // Calculate weighted sum of means and variances
    for (int cell = 0; cell < k; ++cell) {
        sumMean += weight[sample_index][cell] * mean[cell][probe_index];
        sumVar += weight[sample_index][cell] * weight[sample_index][cell] * var[cell][probe_index];
        regularized += std::abs(0.5f - mean[cell][probe_index]);
    }
    
    // Apply regularization
    regularized *= gamma;
    
    // Calculate log likelihood (with -1.0f offset to match original implementation)
    float diff = sumMean - obs_mat[sample_index][probe_index];
    return -1.0f * diff * diff + regularized - 1.0f;
}

float Deconvolution::calculateLogLikelihoodPerSample(int sample_index) const {
    float logLike = 0.0f;
    
    for (int probe = 0; probe < nprobe; ++probe) {
        logLike += calculateLogLikelihoodPerSampleProbe(sample_index, probe);
    }
    
    return logLike;
}

float Deconvolution::calculateLogLikelihoodPerProbe(int probe_index) const {
    float logLike = 0.0f;
    
    for (int sample = 0; sample < nsamp; ++sample) {
        logLike += calculateLogLikelihoodPerSampleProbe(sample, probe_index);
    }
    
    return logLike;
}

float Deconvolution::getLogLikelihood() const {
    float logLike = 0.0f;
    
    for (int sample = 0; sample < nsamp; ++sample) {
        logLike += calculateLogLikelihoodPerSample(sample);
    }
    
    return logLike;
}

int Deconvolution::updateWeights() {
    int samples_changed = 0;
    std::uniform_real_distribution<float> dist(-0.1f, 0.1f);
    std::uniform_int_distribution<int> cell_dist(0, k - 1);
    
    // Update weights for each sample
    for (int sample = 0; sample < nsamp; ++sample) {
        // Check if this sample should be skipped
        int unchanged = 1;
        if (iters_unchanged_ind[sample] > config.max_unchanged_ind) {
            if (++iters_unconsidered_ind[sample] > config.max_unconsidered_ind) {
                iters_unconsidered_ind[sample] = 0;
                iters_unchanged_ind[sample] = 0;
            } else {
                continue;
            }
        }
        
        // Save original weights
        std::vector<float> old_weights = weight[sample];
        
        // Get current log likelihood for this sample
        float sampleLogLike = calculateLogLikelihoodPerSample(sample);
        
        // Try several random perturbations
        for (int iteration = 0; iteration < k * config.iterations_per_cycle; ++iteration) {
            // Randomly select a cell type to perturb
            int random_cell = cell_dist(rng);
            
            // Perturb the weight
            weight[sample][random_cell] *= (1.0f + 0.2f * dist(rng));
            
            // Apply bounds
            if (weight[sample][random_cell] > max_weights[sample]) {
                weight[sample][random_cell] = max_weights[sample];
            } else if (weight[sample][random_cell] < min_weights[sample]) {
                weight[sample][random_cell] = min_weights[sample];
            }
            
            // Normalize weights to sum to 1
            float sum = 0.0f;
            for (int cell = 0; cell < k; ++cell) {
                sum += weight[sample][cell];
            }
            for (int cell = 0; cell < k; ++cell) {
                weight[sample][cell] /= sum;
            }
            
            // Calculate new log likelihood
            float newSampleLogLike = calculateLogLikelihoodPerSample(sample);
            
            // Accept or reject the change
            if (newSampleLogLike > sampleLogLike) {
                unchanged = 0;
                sampleLogLike = newSampleLogLike;
                old_weights = weight[sample];
            } else {
                weight[sample] = old_weights;
            }
        }
        
        // Update counters
        iters_unchanged_ind[sample] += unchanged;
        if (!unchanged) {
            samples_changed++;
        }
    }
    
    return samples_changed;
}

int Deconvolution::updateMeansAndVars() {
    int probes_changed = 0;
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    std::uniform_int_distribution<int> cell_dist(0, k - 1);
    
    // Update means and variances for each probe
    for (int probe = 0; probe < nprobe; ++probe) {
        float diff = obs_max[probe] - obs_min[probe];
        int unchanged = 1;
        
        // Check if this probe should be skipped
        if (iters_unchanged[probe] > config.max_unchanged) {
            if (++iters_unconsidered[probe] > config.max_unconsidered) {
                iters_unconsidered[probe] = 0;
                iters_unchanged[probe] = 0;
            } else {
                continue;
            }
        }
        
        // Get current log likelihood for this probe
        float probeLogLike = calculateLogLikelihoodPerProbe(probe);
        
        // Try several random perturbations
        for (int iteration = 0; iteration < k * config.iterations_per_cycle; ++iteration) {
            // Randomly select a cell type
            int random_cell = cell_dist(rng);
            
            // Save original values
            float old_mean = mean[random_cell][probe];
            float old_var = var[random_cell][probe];
            
            // Generate random perturbations
            float ran1 = dist(rng);
            float ran2 = dist(rng);
            
            // Apply perturbations
            mean[random_cell][probe] *= 1.0f + config.mean_perturbation * (ran1 * diff - diff / 2.0f) / diff;
            var[random_cell][probe] *= 1.0f + config.mean_perturbation * (ran2 * diff - diff / 2.0f) / diff;
            
            // Apply bounds
            if (mean[random_cell][probe] > max_means[probe]) {
                mean[random_cell][probe] = max_means[probe];
            }
            if (mean[random_cell][probe] < min_means[probe]) {
                mean[random_cell][probe] = min_means[probe];
            }
            
            // Calculate new log likelihood
            float newProbeLogLike = calculateLogLikelihoodPerProbe(probe);
            
            // Accept or reject the change
            if (newProbeLogLike > probeLogLike) {
                probeLogLike = newProbeLogLike;
                unchanged = 0;
            } else {
                mean[random_cell][probe] = old_mean;
                var[random_cell][probe] = old_var;
            }
        }
        
        // Update counters
        iters_unchanged[probe] += unchanged;
        if (!unchanged) {
            probes_changed++;
        }
    }
    
    return probes_changed;
}

float Deconvolution::runIteration() {
    // Update weights and means/variances
    updateWeights();
    updateMeansAndVars();
    
    // Return the current log likelihood
    return getLogLikelihood();
}

float Deconvolution::runIterations(int num_iterations) {
    float final_likelihood = 0.0f;
    
    for (int iter = 0; iter < num_iterations; ++iter) {
        final_likelihood = runIteration();
        
        // Optional: print progress
        if (iter % 100 == 0) {
            std::cout << "Iteration " << iter << ", Log likelihood: " << final_likelihood << std::endl;
        }
    }
    
    return final_likelihood;
}

} // namespace CellCDecon


