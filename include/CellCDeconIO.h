#pragma once

#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <limits>

/**
 * @brief Functions for file I/O operations in CellCDecon
 */
namespace CellCDeconIO {

/**
 * @brief Process an input data file
 * 
 * @param filename Path to the input file
 * @param nsamp Number of samples in the data
 * @param colskip Number of columns to skip (metadata columns)
 * @param samplePrefix Reference to store sample prefix information
 * @param sample_ids Vector to store sample IDs
 * @param probe_ids Vector to store probe IDs
 * @param obs_mat Matrix to store observation data
 * @param obs_mean Vector to store mean values for each probe
 * @param obs_var Vector to store variance values for each probe
 * @param obs_min Vector to store minimum values for each probe
 * @param obs_max Vector to store maximum values for each probe
 * @param nprobe Reference to store the actual number of probes read
 * @param maxProbe Maximum number of probes to read
 * @return bool True if file was processed successfully
 */
bool processFile(
    const std::string& filename,
    int nsamp,
    int colskip,
    std::string& samplePrefix,
    std::vector<std::string>& sample_ids,
    std::vector<std::string>& probe_ids,
    std::vector<std::vector<float>>& obs_mat,
    std::vector<float>& obs_mean,
    std::vector<float>& obs_var,
    std::vector<float>& obs_min,
    std::vector<float>& obs_max,
    int& nprobe,
    int maxProbe
);

/**
 * @brief Write results to output files
 * 
 * @param filename Base filename for output
 * @param k Number of cell types
 * @param seed Random seed used
 * @param infer_weight Matrix of inferred weights
 * @param nsamp Number of samples
 * @param nprobe Number of probes
 * @param probe_ids Vector of probe IDs
 * @param obs_mat Matrix of observations
 * @param infer_mean Matrix of inferred means
 * @param infer_var Matrix of inferred variances
 * @param obs_mean Vector of observed means
 * @param obs_var Vector of observed variances
 * @param samplePrefix Sample prefix string
 * @param sample_ids Vector of sample IDs
 * @param gamma Regularization parameter
 * @return bool True if files were written successfully
 */
bool writeFiles(
    const std::string& filename,
    int k,
    int seed,
    const std::vector<std::vector<float>>& infer_weight,
    int nsamp,
    int nprobe,
    const std::vector<std::string>& probe_ids,
    const std::vector<std::vector<float>>& obs_mat,
    const std::vector<std::vector<float>>& infer_mean,
    const std::vector<std::vector<float>>& infer_var,
    const std::vector<float>& obs_mean,
    const std::vector<float>& obs_var,
    const std::string& samplePrefix,
    const std::vector<std::string>& sample_ids,
    float gamma
);

} // namespace CellCDeconIO 