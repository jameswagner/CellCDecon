#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sstream>
#include <iostream>
#include <time.h>
#include <float.h>
#define PI           3.14159265358979323846

using namespace std;

void fillLog(float *logTable, int size) ;

inline float logLikelihoodPerSampleProbe(float **obs_mat,  int k, float **infer_weight, float **mean_mat, float **var_mat, int sample_index, int probe_index, float *logTable) ;
float logLikelihoodPerSample(float **obs_meth_mat, int nprobe, int k, float **infer_weight,  float **mean_mat, float **var_mat, int sample_index, float *logTable) ;
float logLikelihoodPerProbe(float **obs_meth_mat,   int nsamp, int k, float **infer_weight, float **mean_mat, float **var_mat, int probe_index, float *logTable) ;
float logLikelihood(float **obs_meth_mat, int nprobe, int nsamp, int k, float **infer_weight,    float **mean_mat, float **var_mat, float *logTable) ;

void write_files(string filename, int k, int seed, float** infer_weight, int nsamp, int nprobe, string *probe_ids, float** obs_mat, float **infer_mean, float **infer_var, float *obs_mean, float *obs_var, string samplePrefix, string *sample_ids) ;

void processFile(string filename, int nsamp, int colskip, string &samplePrefix, string *sample_ids, string *probe_ids, float **obs_mat, float *obs_mean, float *obs_var, float* obs_min, float* obs_max, int &nprobe, int maxProbe) ;
void initialize_weight(float **weight, int k, int nsamp, float* min_weights, float* max_weights) ;
void initialize_meanvar(float **infer_mean, float** infer_var, float* mins, float* maxes, float* obs_mean, float* obs_var, int k, int nprobe) ;
void update_weights(float **infer_weight, float* max_weights, float* min_weights, float **obs_mat, float **infer_mean, float **infer_var, int *iters_unconsidered_ind, int* iters_unchanged_ind, int max_unchanged_ind, int max_unconsidered_ind, int k, int nsamp, int nprobe, float *logTable ) ;

void update_meanvar(float **obs_mat, float **infer_mean, float **infer_var,  int nprobe, int k, int *iters_unchanged, int *iters_unconsidered, int max_unchanged, int max_unconsidered, float *max_means, float *min_means, float **infer_weight, float *obs_min, float *obs_max, int nsamp, float *logTable) ;
