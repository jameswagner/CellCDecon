#pragma once

#include <cmath>
#include <vector>
#include <string>
#include <memory>
#include <algorithm>
#include <limits>
#include <random>
#include "CellCDeconIO.h"

namespace CellCDecon {

constexpr double PI = 3.14159265358979323846;

/**
 * @brief Cell type deconvolution algorithm implementation
 * 
 * This class encapsulates the entire cell type deconvolution algorithm,
 * handling all state and operations internally to eliminate parameter-heavy
 * function calls and improve code organization.
 */
class Deconvolution {
public:
    /**
     * @brief Configuration for the algorithm
     */
    struct Config {
        int max_unchanged = 10;          ///< Max iterations where probe doesn't change before taking a break
        int max_unconsidered = 100;      ///< Number of iterations a probe "takes a break"
        int max_unchanged_ind = 20;      ///< Max iterations where sample doesn't change before taking a break
        int max_unconsidered_ind = 75;   ///< Number of iterations a sample "takes a break"
        float min_weight_bound = 0.05f;  ///< Minimum allowed weight for any cell type
        float max_weight_bound = 0.95f;  ///< Maximum allowed weight for any cell type
        float mean_perturbation = 0.25f; ///< Amount to perturb means by in optimization
        int iterations_per_cycle = 3;    ///< Number of iterations to try per cell type
    };

    /**
     * @brief Construct a new Deconvolution object
     * 
     * @param observations Matrix of observed values [samples][probes]
     * @param k Number of cell types to deconvolve
     * @param gamma Regularization parameter
     */
    Deconvolution(
        const std::vector<std::vector<float>>& observations,
        int k,
        float gamma = 0.0f
    );

    /**
     * @brief Construct a new Deconvolution object with custom configuration
     * 
     * @param observations Matrix of observed values [samples][probes]
     * @param k Number of cell types to deconvolve
     * @param gamma Regularization parameter
     * @param config Configuration parameters
     */
    Deconvolution(
        const std::vector<std::vector<float>>& observations,
        int k,
        float gamma,
        const Config& config
    );

    /**
     * @brief Get the number of samples
     * @return int Number of samples
     */
    int getNumSamples() const { return nsamp; }

    /**
     * @brief Get the number of probes
     * @return int Number of probes
     */
    int getNumProbes() const { return nprobe; }

    /**
     * @brief Get the number of cell types
     * @return int Number of cell types
     */
    int getNumCellTypes() const { return k; }

    /**
     * @brief Set the regularization parameter
     * @param gamma Regularization parameter
     */
    void setGamma(float gamma) { this->gamma = gamma; }

    /**
     * @brief Get the regularization parameter
     * @return float Regularization parameter
     */
    float getGamma() const { return gamma; }

    /**
     * @brief Run a single iteration of the algorithm
     * @return float The current log likelihood
     */
    float runIteration();

    /**
     * @brief Run multiple iterations of the algorithm
     * @param num_iterations Number of iterations to run
     * @return float The final log likelihood
     */
    float runIterations(int num_iterations);

    /**
     * @brief Get the inferred cell type weights
     * @return const std::vector<std::vector<float>>& Weights matrix [samples][cell_types]
     */
    const std::vector<std::vector<float>>& getWeights() const { return weight; }

    /**
     * @brief Get the inferred means
     * @return const std::vector<std::vector<float>>& Means matrix [cell_types][probes]
     */
    const std::vector<std::vector<float>>& getMeans() const { return mean; }

    /**
     * @brief Get the inferred variances
     * @return const std::vector<std::vector<float>>& Variances matrix [cell_types][probes]
     */
    const std::vector<std::vector<float>>& getVars() const { return var; }

    /**
     * @brief Get the observed data
     * @return const std::vector<std::vector<float>>& Observations matrix [samples][probes]
     */
    const std::vector<std::vector<float>>& getObservations() const { return obs_mat; }

    /**
     * @brief Calculate the overall log likelihood with current parameters
     * @return float Log likelihood
     */
    float getLogLikelihood() const;

private:
    // Core algorithm parameters
    int k;                                         ///< Number of cell types
    int nsamp;                                     ///< Number of samples
    int nprobe;                                    ///< Number of probes
    float gamma;                                   ///< Regularization parameter
    Config config;                                 ///< Algorithm configuration

    // Data matrices
    std::vector<std::vector<float>> obs_mat;       ///< Observed values [samples][probes]
    std::vector<std::vector<float>> weight;        ///< Inferred weights [samples][cell_types]
    std::vector<std::vector<float>> mean;          ///< Inferred means [cell_types][probes]
    std::vector<std::vector<float>> var;           ///< Inferred variances [cell_types][probes]

    // Statistics about observed data
    std::vector<float> obs_mean;                   ///< Mean of each probe
    std::vector<float> obs_var;                    ///< Variance of each probe
    std::vector<float> obs_min;                    ///< Minimum value of each probe
    std::vector<float> obs_max;                    ///< Maximum value of each probe

    // Bounds for optimization
    std::vector<float> min_weights;                ///< Minimum weight for each sample
    std::vector<float> max_weights;                ///< Maximum weight for each sample
    std::vector<float> min_means;                  ///< Minimum mean for each probe
    std::vector<float> max_means;                  ///< Maximum mean for each probe

    // Iteration tracking
    std::vector<int> iters_unchanged;              ///< Count of iterations each probe unchanged
    std::vector<int> iters_unconsidered;           ///< Count of iterations each probe unconsidered
    std::vector<int> iters_unchanged_ind;          ///< Count of iterations each sample unchanged
    std::vector<int> iters_unconsidered_ind;       ///< Count of iterations each sample unconsidered

    // Random number generation
    std::mt19937 rng;                              ///< Random number generator

    /**
     * @brief Calculate statistics for observed data
     */
    void calculateObservedStatistics();

    /**
     * @brief Initialize weights randomly
     */
    void initializeWeights();

    /**
     * @brief Initialize means and variances
     */
    void initializeMeans();

    /**
     * @brief Update weights for all samples in one iteration
     * @return int Number of samples that changed
     */
    int updateWeights();

    /**
     * @brief Update means and variances for all probes in one iteration
     * @return int Number of probes that changed
     */
    int updateMeansAndVars();

    /**
     * @brief Calculate log likelihood for a specific sample and probe
     * 
     * @param sample_index Index of the sample
     * @param probe_index Index of the probe
     * @return float Log likelihood
     */
    float calculateLogLikelihoodPerSampleProbe(int sample_index, int probe_index) const;

    /**
     * @brief Calculate log likelihood for a specific sample across all probes
     * 
     * @param sample_index Index of the sample
     * @return float Log likelihood
     */
    float calculateLogLikelihoodPerSample(int sample_index) const;

    /**
     * @brief Calculate log likelihood for a specific probe across all samples
     * 
     * @param probe_index Index of the probe
     * @return float Log likelihood
     */
    float calculateLogLikelihoodPerProbe(int probe_index) const;
};

} // namespace CellCDecon
