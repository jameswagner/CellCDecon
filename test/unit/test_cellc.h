#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <gtest/gtest.h>
#include "CellCDecon.h"

// Forward declarations of functions we're going to test are no longer needed
// since they're already included from CellCDecon.h

// Helper class for loading test data
class TestDataLoader {
public:
    // Load expected values for different test types
    static std::vector<float> LoadSampleLogLikeData(const std::string& filename);
    static std::vector<float> LoadProbeLogLikeData(const std::string& filename);
    static std::vector<std::vector<float>> LoadSampleProbeLogLikeData(const std::string& filename, int numSamples, int numProbes);
    static float LoadLogLikeData(const std::string& filename);
    
    // Load true weights and means for verification
    static std::vector<std::vector<float>> LoadTrueWeights(const std::string& filename, int nsamp, int k);
    static std::vector<std::vector<float>> LoadTrueMeans(const std::string& filename, int k, int nprobe);
    
    // Load observation data
    static void LoadObservationData(const std::string& filename, 
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
                                  int maxProbe);
};

// Utility functions for test matrices
class TestMatrixUtils {
public:
    // Create a 2D matrix with the specified dimensions
    static float** CreateMatrix(int rows, int cols);
    
    // Free a 2D matrix
    static void FreeMatrix(float** matrix, int rows);
    
    // Print a matrix for debugging
    static void PrintMatrix(float** matrix, int rows, int cols, const std::string& label);
};

// Constants for test file paths
namespace TestConstants {
    const std::string TEST_DATA_DIR = "../test/data/samples/";
    const std::string TEST_BETA_FILE = TEST_DATA_DIR + "testbeta.txt";
    const std::string TRUE_MEAN_FILE = TEST_DATA_DIR + "testbeta.txt.means";
    const std::string TRUE_WEIGHT_FILE = TEST_DATA_DIR + "testbeta.txt.cellcomp";
    const std::string SAMPLEPROBELOG_FILE = TEST_DATA_DIR + "R_output.sampleprobelog.txt";
    const std::string SAMPLELOG_FILE = TEST_DATA_DIR + "R_output.samplelog.txt";
    const std::string PROBELOG_FILE = TEST_DATA_DIR + "R_output.probelog.txt";
    const std::string LOG_FILE = TEST_DATA_DIR + "R_output.log.txt";
}

// A base class for CellCDecon tests
class CellCDeconTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Common setup for all tests
    }
    
    void TearDown() override {
        // Common cleanup for all tests
    }
    
    // Helper method to compare matrices
    void ExpectMatricesNear(float** expected, float** actual, int rows, int cols, float tolerance) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                EXPECT_NEAR(expected[i][j], actual[i][j], tolerance) 
                    << "Matrices differ at [" << i << "][" << j << "]";
            }
        }
    }
    
    // Helper method to compare vectors
    void ExpectVectorsNear(const std::vector<float>& expected, const std::vector<float>& actual, float tolerance) {
        ASSERT_EQ(expected.size(), actual.size()) << "Vector sizes differ";
        for (size_t i = 0; i < expected.size(); i++) {
            EXPECT_NEAR(expected[i], actual[i], tolerance) 
                << "Vectors differ at index " << i;
        }
    }
};

// Function declarations for test utility functions
float calculateCorrelation(const std::vector<float>& x, const std::vector<float>& y);

// Base fixture for deconvolution tests
class DeconBaseTest : public ::testing::Test {
protected:
    // Files for test data
    const std::string test_beta_file = "../test/data/samples/testbeta.txt";
    const std::string true_mean_file = "../test/data/samples/testbeta.txt.means";
    const std::string true_weight_file = "../test/data/samples/testbeta.txt.cellcomp";
    
    // Regularization parameter
    float gamma = 0.0f;
    
    // Helper method to read expected log likelihood values from file
    std::vector<float> ReadLogLikelihoodFile(const std::string& filename);
    
    // Helper method to read overall log likelihood from file
    float ReadSingleLogLikelihood(const std::string& filename);
    
    // Read 2D matrix of expected log likelihood values
    std::vector<std::vector<float>> ReadSampleProbeLogLikelihood(
        const std::string& filename, int nsamp, int nprobe);
    
    // Read true means from file
    std::vector<std::vector<float>> ReadTrueMeans(
        const std::string& filename, int k, int nprobe);
    
    // Read true weights from file
    std::vector<std::vector<float>> ReadTrueWeights(
        const std::string& filename, int nsamp, int k);
};
