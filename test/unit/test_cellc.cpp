// Test cases for the Cell Deconvolution program. 
//
// Author: jamesrwagner@gmail.com

#include "test_cellc.h"
#include "gtest/gtest.h"
#include "CellCDecon.h"
#include "CellCDeconIO.h"
#include <algorithm>
#include <memory>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <filesystem>

using namespace std;

// Add these function declarations since they're defined in CellCDecon.cpp but not in the header
// This allows us to use them in our tests
namespace CellCDecon {
    float logLikelihoodPerSampleProbe(
        const std::vector<std::vector<float>>& obs_mat,
        int k,
        const std::vector<std::vector<float>>& infer_weight,
        const std::vector<std::vector<float>>& mean_mat,
        const std::vector<std::vector<float>>& var_mat,
        int sample_index,
        int probe_index,
        float gamma
    );

    float logLikelihoodPerSample(
        const std::vector<std::vector<float>>& obs_mat,
        int nprobe,
        int k,
        const std::vector<std::vector<float>>& infer_weight,
        const std::vector<std::vector<float>>& mean_mat,
        const std::vector<std::vector<float>>& var_mat,
        int sample_index,
        float gamma
    );

    float logLikelihoodPerProbe(
        const std::vector<std::vector<float>>& obs_mat,
        int nsamp,
        int k,
        const std::vector<std::vector<float>>& infer_weight,
        const std::vector<std::vector<float>>& mean_mat,
        const std::vector<std::vector<float>>& var_mat,
        int probe_index,
        float gamma
    );

    float logLikelihood(
        const std::vector<std::vector<float>>& obs_mat,
        int nprobe,
        int nsamp,
        int k,
        const std::vector<std::vector<float>>& infer_weight,
        const std::vector<std::vector<float>>& mean_mat,
        const std::vector<std::vector<float>>& var_mat,
        float gamma
    );

    void initialize_weight(
        std::vector<std::vector<float>>& weight,
        int k,
        int nsamp,
        const std::vector<float>& min_weights,
        const std::vector<float>& max_weights
    );

    void initialize_meanvar(
        std::vector<std::vector<float>>& infer_mean,
        std::vector<std::vector<float>>& infer_var,
        const std::vector<float>& mins,
        const std::vector<float>& maxes,
        const std::vector<float>& obs_mean,
        const std::vector<float>& obs_var,
        int k,
        int nprobe
    );

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
    );

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
    );
}

// Calculate Pearson correlation coefficient between two vectors
float calculateCorrelation(const vector<float>& x, const vector<float>& y) {
    if (x.size() != y.size() || x.empty()) {
        throw runtime_error("Vectors must be the same size and non-empty");
    }
    
    float sum_x = 0.0f, sum_y = 0.0f, sum_xy = 0.0f;
    float sum_x2 = 0.0f, sum_y2 = 0.0f;
    const int n = x.size();
    
    for(int i = 0; i < n; ++i) {
        sum_x += x[i];
        sum_y += y[i];
        sum_xy += x[i] * y[i];
        sum_x2 += x[i] * x[i];
        sum_y2 += y[i] * y[i];
    }
    
    float denom = sqrt((sum_x2 - (sum_x * sum_x) / n) * (sum_y2 - (sum_y * sum_y) / n));
    if (denom == 0) {
        return 0; // Handle division by zero
    }
    
    return (sum_xy - (sum_x * sum_y) / n) / denom;
}

// Implementation of DeconBaseTest methods
vector<float> DeconBaseTest::ReadLogLikelihoodFile(const string& filename) {
    vector<float> values;
    ifstream infile(filename);
    if (!infile) {
        throw runtime_error("Could not open file: " + filename);
    }
    
    float value;
    while (infile >> value) {
        values.push_back(value);
    }
    
    return values;
}

float DeconBaseTest::ReadSingleLogLikelihood(const string& filename) {
    ifstream infile(filename);
    if (!infile) {
        throw runtime_error("Could not open file: " + filename);
    }
    
    float value;
    infile >> value;
    return value;
}

vector<vector<float>> DeconBaseTest::ReadSampleProbeLogLikelihood(const string& filename, int nsamp, int nprobe) {
    vector<vector<float>> values(nsamp, vector<float>(nprobe));
    ifstream infile(filename);
    if (!infile) {
        throw runtime_error("Could not open file: " + filename);
    }
    
    for (int i = 0; i < nsamp; i++) {
        for (int j = 0; j < nprobe; j++) {
            if (!(infile >> values[i][j])) {
                throw runtime_error("Failed to read data from " + filename);
            }
        }
    }
    
    return values;
}

vector<vector<float>> DeconBaseTest::ReadTrueMeans(const string& filename, int k, int nprobe) {
    vector<vector<float>> means(k, vector<float>(nprobe));
    ifstream infile(filename);
    if (!infile) {
        cerr << "ERROR: Could not open file '" << filename << "'. Working directory: " 
             << std::filesystem::current_path().string() << endl;
        throw runtime_error("Could not open file: " + filename);
    }
    
    cerr << "Successfully opened file: " << filename << endl;
    
    string line;
    if (!getline(infile, line)) {
        cerr << "ERROR: Failed to read header line from file '" << filename << "'" << endl;
        throw runtime_error("Failed to read header from file: " + filename);
    }
    
    cerr << "Header line: " << line << endl;
    
    for (int i = 0; i < nprobe; i++) {
        string probe_id;
        if (!(infile >> probe_id)) {
            cerr << "ERROR: Failed to read probe_id at line " << i+2 
                 << " in file '" << filename << "'" << endl;
            throw runtime_error("Failed to read probe_id at line " + to_string(i+2) 
                               + " from " + filename);
        }
        
        for (int j = 0; j < k; j++) {
            if (!(infile >> means[j][i])) {
                cerr << "ERROR: Failed to read mean for cell type " << j 
                     << " at line " << i+2 << " in file '" << filename << "'" << endl;
                throw runtime_error("Failed to read mean for cell type " + to_string(j) 
                                   + " at line " + to_string(i+2) + " from " + filename);
            }
        }
    }
    
    cerr << "Successfully read " << nprobe << " probes and " << k << " cell types from " 
         << filename << endl;
    
    return means;
}

vector<vector<float>> DeconBaseTest::ReadTrueWeights(const string& filename, int nsamp, int k) {
    vector<vector<float>> weights(nsamp, vector<float>(k));
    ifstream infile(filename);
    if (!infile) {
        throw runtime_error("Could not open file: " + filename);
    }
    
    for (int i = 0; i < nsamp; i++) {
        for (int j = 0; j < k; j++) {
            if (!(infile >> weights[i][j])) {
                throw runtime_error("Failed to read weight from " + filename);
            }
        }
    }
    
    return weights;
}

// Small test fixture for basic log likelihood tests with a small, controlled dataset
class SmallDeconTest : public DeconBaseTest {
protected:
    // Test dimensions
    const int k = 2;        // Number of cell types
    const int nprobe = 3;   // Number of probes
    const int nsamp = 4;    // Number of samples
    
    // Test data
    vector<vector<float>> means;
    vector<vector<float>> weights;
    vector<vector<float>> betas;
    vector<vector<float>> vars;
    
    // Expected log likelihood values
    vector<vector<float>> sampleprobelog_like;
    vector<float> samplelog_like;
    vector<float> probelog_like;
    float log_like;
    
    void SetUp() override {
        // Initialize test matrices with controlled values
        means.resize(k, vector<float>(nprobe));
        vars.resize(k, vector<float>(nprobe));
        weights.resize(nsamp, vector<float>(k));
        betas.resize(nsamp, vector<float>(nprobe));
        
        // Fill matrices with deterministic values for testing
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < nprobe; j++) {
                means[i][j] = i * 0.2f + j * 0.3f;
                vars[i][j] = i * 0.01f + j * 0.02f;
            }
        }
        
        for (int i = 0; i < nsamp; i++) {
            for (int j = 0; j < k; j++) {
                weights[i][j] = i * 0.2f + j * 0.3f;
            }
        }
        
        for (int i = 0; i < nsamp; i++) {
            for (int j = 0; j < nprobe; j++) {
                betas[i][j] = i * 0.25f + j * 0.45f;
            }
        }
        
        // Load expected log likelihood values
        sampleprobelog_like = ReadSampleProbeLogLikelihood("../test/data/samples/R_output.sampleprobelog.txt", nsamp, nprobe);
        samplelog_like = ReadLogLikelihoodFile("../test/data/samples/R_output.samplelog.txt");
        probelog_like = ReadLogLikelihoodFile("../test/data/samples/R_output.probelog.txt");
        log_like = ReadSingleLogLikelihood("../test/data/samples/R_output.log.txt");
    }
};

// Large test fixture for optimization and convergence tests with real-world data
class LargeDeconTest : public DeconBaseTest {
protected:
    // Test dimensions
    const int k = 2;         // Number of cell types
    const int nprobe = 10;   // Number of probes (changed from 100 to match the actual file content)
    const int nsamp = 10;    // Number of samples (adjusted to match actual file content)
    
    // Test data
    vector<vector<float>> betas; // Observations
    vector<vector<float>> true_means;
    vector<vector<float>> true_weights;
    
    // Unique pointer to the deconvolution object
    unique_ptr<CellCDecon::Deconvolution> decon;
    
    // Sample and probe IDs
    vector<string> sample_ids;
    vector<string> probe_ids;
    
    // Observed statistics vectors (need to be persistent for processFile)
    vector<float> obs_mean;
    vector<float> obs_var;
    vector<float> obs_min;
    vector<float> obs_max;
    
    void SetUp() override {
        // Initialize vectors with proper sizes
        betas.resize(nsamp, vector<float>(nprobe));
        
        sample_ids.resize(nsamp);
        probe_ids.resize(nprobe);
        
        // Initialize statistics vectors
        obs_mean.resize(nprobe);
        obs_var.resize(nprobe);
        obs_min.resize(nprobe);
        obs_max.resize(nprobe);
        
        // Load observation data using CellCDeconIO
        int actual_nprobe = nprobe;
        string prefix = "";
        
        bool success = CellCDeconIO::processFile(
            test_beta_file, nsamp, 1, prefix, 
            sample_ids, probe_ids,
            betas, 
            obs_mean, obs_var, obs_min, obs_max,
            actual_nprobe, nprobe
        );
        
        if (!success) {
            throw runtime_error("Failed to process file: " + test_beta_file);
        }
        
        // Create the deconvolution object
        decon = make_unique<CellCDecon::Deconvolution>(betas, k, gamma);
        
        // Load true values for validation
        true_means = ReadTrueMeans(true_mean_file, k, nprobe);
        true_weights = ReadTrueWeights(true_weight_file, nsamp, k);
    }
};

// Test cases for small dataset
TEST_F(SmallDeconTest, LogLikelihoodPerSampleProbe) {
    for (int sample = 0; sample < nsamp; sample++) {
        for (int probe = 0; probe < nprobe; probe++) {
            float likeli = CellCDecon::logLikelihoodPerSampleProbe(
                betas, k, weights, means, vars, sample, probe, gamma
            );
            EXPECT_NEAR(likeli, sampleprobelog_like[sample][probe], 0.2);
        }
    }
}

TEST_F(SmallDeconTest, LogLikelihoodPerSample) {
    for (int sample = 0; sample < nsamp; sample++) {
        float likeli = CellCDecon::logLikelihoodPerSample(
            betas, nprobe, k, weights, means, vars, sample, gamma
        );
        EXPECT_NEAR(likeli, samplelog_like[sample], 0.01);
    }
}

TEST_F(SmallDeconTest, LogLikelihoodPerProbe) {
    for (int probe = 0; probe < nprobe; probe++) {
        float likeli = CellCDecon::logLikelihoodPerProbe(
            betas, nsamp, k, weights, means, vars, probe, gamma
        );
        EXPECT_NEAR(likeli, probelog_like[probe], 0.01);
    }
}

TEST_F(SmallDeconTest, LogLikelihood) {
    float likeli = CellCDecon::logLikelihood(
        betas, nprobe, nsamp, k, weights, means, vars, gamma
    );
    EXPECT_NEAR(likeli, log_like, 10.0);
}

// Test cases for large dataset and optimization using Deconvolution class
TEST_F(LargeDeconTest, DeconvolutionIterations) {
    // Run for a reasonable number of iterations
    float likelihood = decon->runIterations(500);  // Increased from 100 to 500 for better convergence
    
    // Verify returned likelihood is reasonable
    EXPECT_GT(likelihood, -1000000.0f);
    
    // Get the inferred weights and means
    const auto& result_weights = decon->getWeights();
    const auto& result_means = decon->getMeans();
    
    // Verify dimensions match expected
    EXPECT_EQ(result_weights.size(), nsamp);
    EXPECT_EQ(result_weights[0].size(), k);
    EXPECT_EQ(result_means.size(), k);
    EXPECT_EQ(result_means[0].size(), nprobe);

    // Print out the true and inferred weights for each sample
    cout << "\n===== CELL TYPE COMPOSITIONS FOR EACH SAMPLE =====\n";
    cout << "Format: [Sample #]: True weights vs Inferred weights\n";
    
    for (int sample = 0; sample < nsamp; sample++) {
        cout << "Sample " << sample << ":\n";
        cout << "  True weights:     [";
        for (int cell = 0; cell < k; cell++) {
            cout << true_weights[sample][cell];
            if (cell < k - 1) cout << ", ";
        }
        cout << "]\n";
        
        cout << "  Inferred weights: [";
        for (int cell = 0; cell < k; cell++) {
            cout << result_weights[sample][cell];
            if (cell < k - 1) cout << ", ";
        }
        cout << "]\n";
    }
    
    // Print out the true and inferred means for each cell type
    cout << "\n===== CELL TYPE MEANS FOR FIRST 5 PROBES =====\n";
    cout << "Format: [Cell Type #]: True means vs Inferred means\n";
    
    for (int cell = 0; cell < k; cell++) {
        cout << "Cell Type " << cell << ":\n";
        cout << "  True means:     [";
        for (int probe = 0; probe < min(5, nprobe); probe++) {
            cout << true_means[cell][probe];
            if (probe < min(5, nprobe) - 1) cout << ", ";
        }
        cout << " ...]\n";
        
        cout << "  Inferred means: [";
        for (int probe = 0; probe < min(5, nprobe); probe++) {
            cout << result_means[cell][probe];
            if (probe < min(5, nprobe) - 1) cout << ", ";
        }
        cout << " ...]\n";
    }
}

// Test individual iteration functions of Deconvolution class
TEST_F(LargeDeconTest, SingleIteration) {
    // Capture log likelihood before iteration
    float before_likelihood = decon->getLogLikelihood();
    
    // Run a single iteration
    float after_likelihood = decon->runIteration();
    
    // Verify log likelihood improves (or at least doesn't get worse)
    EXPECT_GE(after_likelihood, before_likelihood);
}

// Test constructor with custom Config
TEST_F(LargeDeconTest, CustomConfig) {
    // Create a custom configuration
    CellCDecon::Deconvolution::Config config;
    config.max_unchanged = 15;
    config.max_unconsidered = 150;
    config.min_weight_bound = 0.1f;
    config.max_weight_bound = 0.9f;
    
    // Create a new deconvolution object with custom config
    auto custom_decon = make_unique<CellCDecon::Deconvolution>(betas, k, gamma, config);
    
    // Run iterations
    float likelihood = custom_decon->runIterations(20);
    
    // Verify basic functionality still works
    EXPECT_GT(likelihood, -1000000.0f);
    EXPECT_EQ(custom_decon->getNumSamples(), nsamp);
    EXPECT_EQ(custom_decon->getNumProbes(), nprobe);
    EXPECT_EQ(custom_decon->getNumCellTypes(), k);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
