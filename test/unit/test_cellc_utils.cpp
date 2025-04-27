#include "test_cellc.h"
#include "CellCDeconIO.h"

// Implementation for TestDataLoader
std::vector<float> TestDataLoader::LoadSampleLogLikeData(const std::string& filename) {
    std::vector<float> result;
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    
    float value;
    while (infile >> value) {
        result.push_back(value);
    }
    
    return result;
}

std::vector<float> TestDataLoader::LoadProbeLogLikeData(const std::string& filename) {
    return LoadSampleLogLikeData(filename); // Same implementation
}

std::vector<std::vector<float>> TestDataLoader::LoadSampleProbeLogLikeData(
    const std::string& filename, int numSamples, int numProbes) {
    
    std::vector<std::vector<float>> result(numSamples, std::vector<float>(numProbes));
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    
    for (int i = 0; i < numSamples; i++) {
        for (int j = 0; j < numProbes; j++) {
            if (!(infile >> result[i][j])) {
                throw std::runtime_error("Failed to read expected data from " + filename);
            }
        }
    }
    
    return result;
}

float TestDataLoader::LoadLogLikeData(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    
    float result;
    if (!(infile >> result)) {
        throw std::runtime_error("Failed to read log like data from " + filename);
    }
    
    return result;
}

std::vector<std::vector<float>> TestDataLoader::LoadTrueWeights(
    const std::string& filename, int nsamp, int k) {
    
    std::vector<std::vector<float>> result(nsamp, std::vector<float>(k));
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    
    for (int i = 0; i < nsamp; i++) {
        for (int j = 0; j < k; j++) {
            if (!(infile >> result[i][j])) {
                throw std::runtime_error("Failed to read true weights from " + filename);
            }
        }
    }
    
    return result;
}

std::vector<std::vector<float>> TestDataLoader::LoadTrueMeans(
    const std::string& filename, int k, int nprobe) {
    
    std::vector<std::vector<float>> result(k, std::vector<float>(nprobe));
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    
    // Skip header line
    std::string line;
    std::getline(infile, line);
    
    // Read data
    for (int i = 0; i < nprobe; i++) {
        // Skip probe ID
        std::string junk;
        infile >> junk;
        
        for (int j = 0; j < k; j++) {
            if (!(infile >> result[j][i])) {
                throw std::runtime_error("Failed to read true means from " + filename);
            }
        }
    }
    
    return result;
}

void TestDataLoader::LoadObservationData(
    const std::string& filename, 
    int nsamp, 
    int colskip,
    std::string& samplePrefix,
    std::vector<std::string>& sample_ids,
    std::vector<std::string>& probe_ids,
    float **obs_mat, 
    std::vector<float>& obs_mean,
    std::vector<float>& obs_var,
    std::vector<float>& obs_min,
    std::vector<float>& obs_max,
    int& nprobe, 
    int maxProbe) {
    
    // Create temporary 2D vector for observation matrix
    std::vector<std::vector<float>> obs_mat_vec(nsamp, std::vector<float>(maxProbe));
    
    // Initialize sample_ids and probe_ids to proper sizes
    sample_ids.resize(nsamp);
    probe_ids.resize(maxProbe);
    
    // Call the CellCDeconIO::processFile with vectors
    bool success = CellCDeconIO::processFile(
        filename, nsamp, colskip, samplePrefix,
        sample_ids, probe_ids, obs_mat_vec,
        obs_mean, obs_var, obs_min, obs_max,
        nprobe, maxProbe
    );
    
    if (!success) {
        throw std::runtime_error("Failed to process file: " + filename);
    }
    
    // Copy data from 2D vector to 2D array
    for (int i = 0; i < nsamp; i++) {
        for (int j = 0; j < nprobe; j++) {
            obs_mat[i][j] = obs_mat_vec[i][j];
        }
    }
}

// Implementation for TestMatrixUtils
float** TestMatrixUtils::CreateMatrix(int rows, int cols) {
    float** matrix = new float*[rows];
    for (int i = 0; i < rows; i++) {
        matrix[i] = new float[cols]();
    }
    return matrix;
}

void TestMatrixUtils::FreeMatrix(float** matrix, int rows) {
    for (int i = 0; i < rows; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;
}

void TestMatrixUtils::PrintMatrix(float** matrix, int rows, int cols, const std::string& label) {
    printf("Matrix %s (%d x %d):\n", label.c_str(), rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%8.4f ", matrix[i][j]);
        }
        printf("\n");
    }
} 