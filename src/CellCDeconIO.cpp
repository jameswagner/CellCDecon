#include "CellCDeconIO.h"
#include <iostream>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <iomanip>

namespace CellCDeconIO {

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
) {
    // Clear and resize output vectors
    sample_ids.resize(nsamp);
    probe_ids.clear();
    obs_mean.clear();
    obs_var.clear();
    obs_min.clear();
    obs_max.clear();
    
    // Make sure the observation matrix has the right dimensions
    obs_mat.resize(nsamp);
    for (auto& row : obs_mat) {
        row.clear();
        row.reserve(maxProbe);
    }
    
    // Open input file
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return false;
    }
    
    // Read header line
    std::string line;
    if (!std::getline(inputFile, line)) {
        std::cerr << "Error: Empty file or could not read header" << std::endl;
        return false;
    }
    
    std::istringstream headerStream(line);
    std::string token;
    
    // Process column headers that don't correspond to sample IDs
    samplePrefix.clear();
    for (int i = 0; i < colskip; ++i) {
        if (!(headerStream >> token)) {
            std::cerr << "Error: Not enough columns in header" << std::endl;
            return false;
        }
        samplePrefix += token + " ";
    }
    
    // Process sample IDs
    for (int i = 0; i < nsamp; ++i) {
        if (!(headerStream >> token)) {
            std::cerr << "Error: Not enough sample columns in header" << std::endl;
            return false;
        }
        sample_ids[i] = token;
    }
    
    // Read data rows
    nprobe = 0;
    while (std::getline(inputFile, line) && nprobe < maxProbe) {
        std::istringstream lineStream(line);
        
        // Read probe ID and metadata
        std::string probeId;
        for (int i = 0; i < colskip; ++i) {
            if (!(lineStream >> token)) {
                std::cerr << "Error: Not enough metadata columns in row " << nprobe + 1 << std::endl;
                return false;
            }
            
            if (i == 0) {
                probeId = token;
            } else {
                probeId += " " + token;
            }
        }
        probe_ids.push_back(probeId);
        
        // Variables for calculating mean and variance
        float line_mean = 0.0f;
        float M2 = 0.0f;
        int non_Na_count = 0;
        float line_min = std::numeric_limits<float>::max();
        float line_max = std::numeric_limits<float>::lowest();
        
        // Read observations for this probe
        for (int i = 0; i < nsamp; ++i) {
            if (!(lineStream >> token)) {
                std::cerr << "Error: Not enough data columns in row " << nprobe + 1 << std::endl;
                return false;
            }
            
            // Handle NA values
            if (token == "NA") {
                obs_mat[i].push_back(-1.0f);
            } else {
                try {
                    float obs = std::stof(token);
                    
                    // Update min/max values
                    line_min = std::min(line_min, obs);
                    line_max = std::max(line_max, obs);
                    
                    // Store observation
                    obs_mat[i].push_back(obs);
                    
                    // Update running statistics (Welford's algorithm)
                    non_Na_count++;
                    float delta = obs - line_mean;
                    line_mean += delta / non_Na_count;
                    float delta2 = obs - line_mean;
                    M2 += delta * delta2;
                } catch (const std::invalid_argument&) {
                    std::cerr << "Error: Invalid number format '" << token << "' in row " << nprobe + 1 << std::endl;
                    return false;
                }
            }
        }
        
        // Store statistics for this probe
        if (non_Na_count > 1) {
            obs_var.push_back(M2 / (non_Na_count - 1));
        } else {
            // Handle the case where there's only one or zero valid observations
            obs_var.push_back(0.0f);
        }
        
        obs_mean.push_back(line_mean);
        obs_min.push_back(line_min);
        obs_max.push_back(line_max);
        
        // Increment probe counter
        nprobe++;
    }
    
    // Check if we loaded any probes
    if (nprobe == 0) {
        std::cerr << "Warning: No probe data was loaded from file" << std::endl;
    }
    
    return true;
}

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
) {
    // Generate base filename
    std::ostringstream baseFilename;
    baseFilename << filename << ".k" << k << ".seed" << seed << ".gamma" << gamma << ".SUMSQ.";
    
    // Create the three output files
    std::string weightFile = baseFilename.str() + "w";
    std::string meanVarFile = baseFilename.str() + "meanvar";
    std::string residFile = baseFilename.str() + "resid";
    
    // Open weight file
    std::ofstream weightStream(weightFile);
    if (!weightStream.is_open()) {
        std::cerr << "Error: Could not open weight file for writing: " << weightFile << std::endl;
        return false;
    }
    
    // Write weights
    weightStream << std::fixed << std::setprecision(6);
    for (int s = 0; s < nsamp; ++s) {
        for (int a = 0; a < k; ++a) {
            weightStream << infer_weight[s][a] << "\t";
        }
        weightStream << "\n";
    }
    weightStream.close();
    
    // Open mean/var file
    std::ofstream meanVarStream(meanVarFile);
    if (!meanVarStream.is_open()) {
        std::cerr << "Error: Could not open mean/var file for writing: " << meanVarFile << std::endl;
        return false;
    }
    
    // Write inferred probe means and variances
    meanVarStream << std::fixed << std::setprecision(6);
    for (int i = 0; i < nprobe; ++i) {
        meanVarStream << probe_ids[i] << " ";
        for (int a = 0; a < k; ++a) {
            meanVarStream << infer_mean[a][i] << " " << infer_var[a][i] << " ";
        }
        meanVarStream << obs_mean[i] << " " << obs_var[i] << "\n";
    }
    meanVarStream.close();
    
    // Open residuals file
    std::ofstream residStream(residFile);
    if (!residStream.is_open()) {
        std::cerr << "Error: Could not open residuals file for writing: " << residFile << std::endl;
        return false;
    }
    
    // Write header for residuals file
    residStream << std::fixed << std::setprecision(6);
    residStream << samplePrefix << " ";
    for (int i = 0; i < nsamp; ++i) {
        residStream << sample_ids[i] << " ";
    }
    residStream << "\n";
    
    // Write residuals (observed - predicted)
    for (int i = 0; i < nprobe; ++i) {
        residStream << probe_ids[i];
        for (int ii = 0; ii < nsamp; ++ii) {
            // Calculate predicted value
            float pred = 0.0f;
            for (int a = 0; a < k; ++a) {
                pred += infer_weight[ii][a] * infer_mean[a][i];
            }
            
            // Write residual if observation is not missing
            if (obs_mat[ii][i] != -1.0f) {
                residStream << " " << obs_mat[ii][i] - pred;
            } else {
                residStream << " NA";
            }
        }
        residStream << "\n";
    }
    residStream.close();
    
    return true;
}

} // namespace CellCDeconIO 