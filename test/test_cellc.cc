// Test cases for the Cell Deconvolution program. 
//
// Author: jamesrwagner@gmail.com


#include "test_cellc.cc"
#include "gtest/gtest.h"
#include <fstream>
#include <time.h>
#include <algorithm>

using ::testing::EmptyTestEventListener;
using ::testing::InitGoogleTest;
using ::testing::Test;
using ::testing::TestCase;
using ::testing::TestEventListeners;
using ::testing::TestInfo;
using ::testing::TestPartResult;
using ::testing::UnitTest;


class DeconTest : public testing::Test {

#include "sample15.h"
 protected:  
  virtual void SetUp() {
    printf("Start of SetUp\n");
    test_beta_file = "samples/testbeta.txt";
    small_k = 2;
    small_nprob = 3;
    small_nsamp = 4;
    small_means = new float*[small_k];
    small_vars = new float*[small_k];
    small_weights = new float*[small_nsamp];
    small_betas = new float*[small_nsamp];
    logTable = new float[10000000]();
    fillLog(logTable, 10000000);
    

    large_k = 2;
    large_nprob = 100;
    large_nsamp = 33;
    large_means = new float*[large_k];
    true_means = new float*[large_k];
    large_vars = new float*[large_k];
    large_weights = new float*[large_nsamp];
    true_weights = new float*[large_nsamp];
    large_betas = new float*[large_nsamp];
    for(int i = 0; i < large_k; i++) {
      large_means[i] = new float[large_nprob];
      true_means[i] = new float[large_nprob];
      large_vars[i] = new float[large_nprob];
    }
    for(int i = 0; i < large_nsamp; i++) {
      large_weights[i] = new float[large_k];
      true_weights[i] = new float[large_k];
    }
    for(int i = 0; i < large_nsamp; i++) {
      large_betas[i] = new float[large_nprob];
    }
    sample_ids = new string[large_nsamp];
    probe_ids = new string[large_nprob];
    obs_mean = new float[large_nprob];
    obs_var = new float[large_nprob];
    obs_min = new float[large_nprob];
    obs_max = new float[large_nprob];
    maxProbe = large_nprob;
    prefix = "";

    max_means = new float[large_nprob]();
    min_means = new float[large_nprob]();

    max_weights = new float[large_nsamp]();
    min_weights = new float[large_nsamp]();


    
    for(int i = 0; i < large_nprob; i++) {
      max_means[i] = 0.9;
      min_means[i] = 0.1;
    }
    
    for(int i = 0; i < large_nsamp; i++) {
      max_weights[i] = 0.9;
      min_weights[i] = 0.1;
    }
    processFile(test_beta_file, large_nsamp, 1, prefix, sample_ids, probe_ids, large_betas, obs_mean, obs_var,  obs_min, obs_max, large_nprob, maxProbe) ;
    initialize_weight(large_weights, large_k, large_nsamp, min_weights, max_weights);
    initialize_meanvar(large_means, large_vars, obs_mean, obs_var, min_means, max_means, large_k, large_nprob);
    for(int i = 0; i < small_k; i++) {
      small_means[i] = new float[small_nprob];
      small_vars[i] = new float[small_nprob];
    }
    for(int i = 0; i < small_k; i++) {
      for(int j = 0; j < small_nprob; j++) {
	small_means[i][j] = i * 0.2 + j * 0.3;
	small_vars[i][j] = i * 0.01 + j * 0.02;
      }
    }
    

    for(int i = 0; i < small_nsamp; i++) {
      small_weights[i] = new float[small_k];
    }

    for(int i =0; i < small_nsamp; i++) {
      for(int j= 0; j < small_k;j++) {
	small_weights[i][j] = i * 0.2 + j * 0.3;
      }
    }
    

    for(int i = 0; i < small_nsamp; i++) {
      small_betas[i] = new float[small_nprob];
    }

    for(int i =0; i < small_nsamp; i++) {
      for(int j= 0; j < small_nprob;j++) {
	small_betas[i][j] = i*0.25 + j * 0.45;
      }
    }


    
    sampleprobelog_like = new float*[small_nsamp]();
    samplelog_like = new float[small_nsamp]();
    probelog_like = new float[small_nprob]();
    log_like = 0.0;

    for(int i = 0; i < small_nsamp; i++) {
      sampleprobelog_like[i] = new float[small_nprob]();
    }
    ifstream infile;
    string R_sampleprobelog_like = "samples/R_output.sampleprobelog.txt";
    infile.open(R_sampleprobelog_like.c_str());
    for(int i = 0; i < small_nsamp; i++) {
      for(int j = 0; j < small_nprob; j++) {
	infile >> sampleprobelog_like[i][j];
      }
    }
    
    infile.close();

    infile.open("samples/R_output.samplelog.txt");
    for(int i = 0; i < small_nsamp; i++) {
      infile >> samplelog_like[i];
    }
    infile.close();



    infile.open("samples/R_output.probelog.txt");
    for(int i = 0; i < small_nprob; i++) {
      infile >> probelog_like[i];
    }
    infile.close();

    infile.open("samples/R_output.log.txt");
    infile >> log_like;
    infile.close();


    true_mean_file = "samples/testbeta.txt.means";
    true_weight_file = "samples/testbeta.txt.cellcomp";
    infile.open(true_mean_file.c_str());
    char junk[2000] ;
    infile.getline(junk,2000); //ignore header line 

    //Get column headers in first line of file that do not correspond to sample ids
    for(int i=0; i < large_nprob; i++) {
      infile >> junk; // discard first column (probe id)
      for(int j = 0; j < large_k; j++) {
	infile >> true_means[j][i];
      }
    }
    infile.close();
    infile.open(true_weight_file.c_str());
    for(int i =0; i < large_nsamp; i++) {
      for(int j = 0; j < large_k; j++) {
	
	infile >> true_weights[i][j];

      }
    }


    iters_unconsidered_ind  = new int[large_nsamp]();
    iters_unchanged_ind  = new int[large_nsamp]();
    iters_unconsidered  = new int[large_nprob]();
    iters_unchanged  = new int[large_nprob]();


  }
  


    int small_k;
  float **small_means;
  float **small_weights; 
  float **small_betas;
  float **small_vars;
  int small_nprob;
  int small_nsamp;


  int large_k;
  float **large_means;
  float **large_weights; 
  float **large_betas;
  float **large_vars;
  int large_nprob;
  int large_nsamp;


  string *sample_ids ;
  string *probe_ids ;
  float *obs_mean ;
  float *obs_var ;
  float *obs_min ;
  float *obs_max ;
  int maxProbe ;
  string prefix;

  float **true_means;
  float **true_weights;
  string test_beta_file;
  string true_mean_file;
  string true_weight_file;
  
  float **sampleprobelog_like;
  float *samplelog_like;
  float *probelog_like;
  float log_like;
  float *logTable;
  
  int *iters_unconsidered_ind;
  int *iters_unchanged_ind;
  int *iters_unconsidered;
  int *iters_unchanged;

  float *max_weights;
  float *min_weights;
  
  float *max_means;
  float *min_means;
};

TEST_F(DeconTest, loglikelihoodPerSampleProbe) {

  for (int sample_index = 0; sample_index < small_nsamp; sample_index++) {
    for(int probe_index = 0; probe_index < small_nprob; probe_index++) {
      float likeli = logLikelihoodPerSampleProbe(small_betas,  small_k, small_weights, small_means, small_vars, sample_index, probe_index, logTable) ;
      EXPECT_NEAR(likeli,sampleprobelog_like[sample_index][probe_index], 0.01);
     }
  }


}


TEST_F(DeconTest, loglikelihoodPerSample) {
  for (int sample_index = 0; sample_index < small_nsamp; sample_index++) {
    float likeli = logLikelihoodPerSample(small_betas, small_nprob, small_k, small_weights, small_means, small_vars, sample_index,  logTable) ;
    EXPECT_NEAR(likeli,samplelog_like[sample_index], 0.01);
  }
}


TEST_F(DeconTest, loglikelihoodPerProbe) {
  for (int probe_index = 0; probe_index < small_nprob; probe_index++) {
    float likeli = logLikelihoodPerProbe(small_betas, small_nsamp, small_k, small_weights, small_means, small_vars, probe_index,  logTable) ;
    EXPECT_NEAR(likeli,probelog_like[probe_index], 0.01);
  }
}
 



TEST_F(DeconTest, loglikelihood) {
  float likeli = logLikelihood(small_betas, small_nprob,  small_nsamp, small_k, small_weights,    small_means, small_vars,  logTable) ;
  EXPECT_NEAR(likeli, log_like, 0.05);
}

TEST_F(DeconTest, update_weight) {
  for(int i =0 ; i < 100 ; i++) {
    float oldlike = logLikelihood(large_betas, large_nprob,  large_nsamp, large_k, large_weights,    large_means, large_vars,  logTable) ;
    update_weights(large_weights, min_weights,  max_weights, large_betas, large_means, large_vars, iters_unconsidered_ind,  iters_unchanged_ind,  1000, 1000,  large_k,  large_nsamp,  large_nprob, logTable ) ;
    float newlike = logLikelihood(large_betas, large_nprob,  large_nsamp, large_k, large_weights,    large_means, large_vars,  logTable) ;
    EXPECT_GE(newlike, oldlike);
    update_meanvar(large_betas, large_means, large_vars,  large_nprob,  large_k, iters_unchanged,  iters_unconsidered, 1000,  1000,  max_means, min_means, large_weights, obs_min, obs_max, large_nsamp, logTable) ;
  }
 
}


TEST_F(DeconTest, update_meanvar) {
  for(int i =0 ; i < 100 ; i++) {
    float oldlike = logLikelihood(large_betas, large_nprob,  large_nsamp, large_k, large_weights,    large_means, large_vars,  logTable) ;
    update_meanvar(large_betas, large_means, large_vars,  large_nprob,  large_k, iters_unchanged,  iters_unconsidered, 1000,  1000,  max_means, min_means, large_weights, obs_min, obs_max, large_nsamp, logTable) ;    

    float newlike = logLikelihood(large_betas, large_nprob,  large_nsamp, large_k, large_weights,    large_means, large_vars,  logTable) ;
    EXPECT_GE(newlike, oldlike);
    update_weights(large_weights, min_weights,  max_weights, large_betas, large_means, large_vars, iters_unconsidered_ind,  iters_unchanged_ind,  1000, 1000,  large_k,  large_nsamp,  large_nprob, logTable ) ;
  }
}


TEST_F(DeconTest, run_iteratively) {
  for(int i = 0; i < 500 ; i++) {
    update_meanvar(large_betas, large_means, large_vars,  large_nprob,  large_k, iters_unchanged,  iters_unconsidered, 1000,  1000,  max_means, min_means, large_weights, obs_min, obs_max, large_nsamp, logTable) ;    
    update_weights(large_weights, min_weights,  max_weights, large_betas, large_means, large_vars, iters_unconsidered_ind,  iters_unchanged_ind,  1000, 1000,  large_k,  large_nsamp,  large_nprob, logTable ) ;
  }
  for(int j =0 ; j < large_nprob; j++) {
    float *large_mean_vec = new float[large_k];
    float *true_mean_vec = new float[large_k];
    for(int i = 0; i < large_k; i++) {
      large_mean_vec[i] = large_means[i][j];
      true_mean_vec[i] = true_means[i][j];
    }
    std::sort(large_mean_vec, large_mean_vec + large_k);
    std::sort(true_mean_vec, true_mean_vec + large_k);
    for(int i =0; i < large_k; i++) {
    

      printf( "%d %d: %f %f\n", i,j,large_mean_vec[i], true_mean_vec[i]);
      EXPECT_NEAR(large_mean_vec[i], true_mean_vec[i], 0.1);
    }

    }

    for(int i = 0; i < large_nsamp; i++) {


      std::sort(large_weights[i], large_weights[i] + large_k);
      std::sort(true_weights[i], true_weights[i] + large_k);

      for(int j =0; j < large_k; j++) {
	printf("%d %d: %f %f\n", i,j, large_weights[i][j], true_weights[i][j]);
	EXPECT_NEAR(large_weights[i][j], true_weights[i][j], 0.05);
	
      }
    }
}








int main(int argc, char **argv) {
  //srand (11);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
