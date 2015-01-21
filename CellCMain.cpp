#include "CellCDecon.h"
#include <stdio.h>
#include <time.h>


int parseArgs(int argc, char* argv[], string &filename, int &k, int &seed, int &maxProbes, int &nsamp, int &colskip, float &gamma) {
  for(int i = 1; i < argc; i+=2) {
    if(!strcmp(argv[i], "-f")) {
      filename = string(argv[i+1]);
    }
    else if(!strcmp(argv[i], "-k")) {
      k = atoi(argv[i+1]);
    }
    else if(!strcmp(argv[i], "-s")) {
      seed = atoi(argv[i+1]);
    }
    else if(!strcmp(argv[i], "-m")) {
      maxProbes = atoi(argv[i+1]);
    }
    else if(!strcmp(argv[i], "-n")) {
      nsamp = atoi(argv[i+1]);
    }
    else if(!strcmp(argv[i], "-c")) {
      colskip = atoi(argv[i+1]);
    }
    else if (!strcmp(argv[i], "-g")) {
      gamma = atof(argv[i+1]);
    }

  }
  if(k < 1 || nsamp < 1 || !strcmp(filename.c_str(), "")) {
    return -1;
  }
  return 1;

}


int main(int argc, char *argv[]) {



  string filename = "";
  int k = 0;
  int seed = time(NULL);
  int maxProbes = 5e5;
  int nsamp = 0;
  int colskip = 1;
  float gamma = 0;
  if(parseArgs(argc, argv, filename, k, seed, maxProbes, nsamp, colskip, gamma) < 1) {
    fprintf(stderr, "Requires arguments -k <number of cell types> -n <number of samples> -f <input filename>\n");
    exit(1);
  }

  srand(seed);

  int log_tableSize = 10000000;
  float *logTable = new float[10000000]();
  fillLog(logTable, log_tableSize);
 

  string *sample_ids = new string[nsamp];
  string *probe_ids = new string[maxProbes];

  float *obs_mean =  new float[maxProbes];
  float *obs_var= new float[maxProbes];

  float *obs_min =  new float[maxProbes];
  float *obs_max= new float[maxProbes];
  float **obs_mat = new float*[nsamp];
  float *diff = new float[maxProbes];
  float *mins =  new float[maxProbes];

  float *maxes =  new float[maxProbes];
  string sample_prefix ;


  int *iters_unchanged = new int[maxProbes]();
  int *iters_unconsidered = new int[maxProbes]();
  int max_unchanged = 10;
  int max_unconsidered = 100;


  int *iters_unchanged_ind = new int[nsamp]();
  int *iters_unconsidered_ind = new int[maxProbes]();
  int max_unchanged_ind = 20;
  int max_unconsidered_ind = 75;




  float **infer_mean = new float*[k];
  float **infer_var = new float*[k];
  float **newmean = new float*[k];
  float **newvar = new float*[k];
  float **infer_weight = new float*[nsamp];



  float *min_means = new float[maxProbes]();
  float *max_means = new float[maxProbes]();
  float *min_weights = new float[nsamp]();
  float *max_weights = new float[nsamp]();

  


  for (int i=0;i<nsamp;i++) {
    obs_mat[i]= new float[maxProbes];
    infer_weight[i] = new float[k]();
    min_weights[i] = 0.05;
    max_weights[i] = 0.95;


  }
  int nprobe = maxProbes;


  processFile(filename, nsamp, colskip, sample_prefix, sample_ids,  probe_ids, obs_mat, obs_mean, obs_var, obs_min, obs_max, nprobe, maxProbes);
  

  for (int i = 0; i < nprobe; i++) {
    min_means[i] = 0;
    max_means[i] = 1.0;
  }

  initialize_weight(infer_weight, k, nsamp, min_weights, max_weights);
  initialize_meanvar(infer_mean, infer_var, min_means, max_means,  obs_mean, obs_var, k, nprobe);
  float original_gamma = gamma;

    
    for (int iter=0;iter<1000;iter++) {
           printf("iter %d\n", iter);
	   if(iter==500) {

             gamma=0;
           }
	   update_weights(infer_weight, min_weights,  max_weights, obs_mat, infer_mean, infer_var, iters_unconsidered_ind,  iters_unchanged_ind,  max_unchanged_ind,  max_unconsidered_ind,  k,  nsamp,  nprobe, logTable, gamma ) ;
	   update_meanvar(obs_mat, infer_mean, infer_var,  nprobe,  k, iters_unchanged,  iters_unconsidered, max_unchanged,  max_unconsidered,  max_means, min_means, infer_weight, obs_min, obs_max, nsamp, logTable, gamma) ;
	   if(iter % 100 == 99) {
	     write_files( filename,  k, seed,  infer_weight,  nsamp, nprobe, probe_ids,  obs_mat, infer_mean, infer_var, obs_mean, obs_var,  sample_prefix, sample_ids, original_gamma) ;
	   }
    }
    write_files( filename,  k, seed,  infer_weight,  nsamp, nprobe, probe_ids,  obs_mat, infer_mean, infer_var, obs_mean, obs_var,  sample_prefix, sample_ids, original_gamma) ;
    
    
  
}


