#include "CellCDecon.h"

float *logTable ;
const float bin_min = 0.0;
const float bin_max = 1.0;
const float bin_diff = bin_max - bin_min;
const int bin_num = 200;
void fillLog(float *logTable, int size) {

  for (int i=0;i<size;i++) {
    logTable[i]=0.5*log(2*PI*(i+1)/(size/100.0));
  }
}

inline float logLikelihoodPerSampleProbe(float **obs_mat,  int k, float **infer_weight, float **mean_mat, float **var_mat, int sample_index, int probe_index, float *logTable) {
  /* This is the "workhorse" function of the current implementation. Given a probe and sample, what is the loglikehood of a  the observed methylation value ?
     input arguments:
     obs_mat: a matrix of observed values with the samples corresponding to rows and the probes corresponding to columns
     k: The number of cell types
     infer_weight: a matrix of weights for individuals, with each row being one individual (sample), and each column corresponding to one of k cell types
     mean_mat: a matrix of the mean methylation values for each probe, with each row being one of k cell types and each column being one of the probes
     var_mat: a matrix of mean methylation values, same dimension as mean
     sample_index: the index of the sample for the log likelihood
     probe_index: the index of the probe for the loglikelihood
     
     return value: a float corresponding to the log likelihood of the observation given the cell type composition for the specified sample and mean methylation values 
     for that probe
     assumption: warning that this function does not check that sample_index is within the bounds of the meth, mean or var matrices, nor that probe_index
     is within the bounds of the mean or infer_weightrices
  */
  
  float sumMean = 0.;
  float sumVar = 0.;
  //  printf("called with sample_index %d\n", sample_index); 
  float *sample_weight = infer_weight[sample_index];

  if(obs_mat[sample_index][probe_index] == -1 ) { //missing value
    return 0 ;
  }
  for (int a=0;a<k;a++) {
    sumMean+=sample_weight[a]*mean_mat[a][probe_index];
    sumVar+=sample_weight[a]*sample_weight[a]*var_mat[a][probe_index];

  }
  float t=obs_mat[sample_index][probe_index]-sumMean;
  //printf("came up with summean %f sumvar %f observed is %f t is %f\n", sumMean, sumVar, obs_mat[sample_index][probe_index], t);
  if(sumVar*100000 > 10000000) {
    return -t*t/sumVar/2- 0.5*log(6.283185*(sumVar*100000+1)/100000.0); // exceeds bounds of log table
  }

  return   (-t*t/sumVar/2-logTable[(int)(sumVar*100000)]); // lookup in logtable
  
}




float logLikelihoodPerSample(float **obs_meth_mat, int nprobe, int k, float **infer_weight,  float **mean_mat, float **var_mat, int sample_index, float *logTable) {
  /*
    Takes log likelihood for given sample, given the observations and the current weights and means 
 input arguments:
     obs_meth_mat: a matrix of observed values with the samples corresponding to rows and the probes corresponding to columns
     nprobe: total number of probes (columns) in obs_meth_mat and mean_mat
     k: The number of cell types
     infer_weight: a matrix of weights for individuals, with each row being one individual (sample), and each column corresponding to one of k cell types
     mean_mat: a matrix of the mean methylation values for each probe, with each row being one of k cell types and each column being one of the probes
     var_mat: a matrix of mean methylation values, same dimension as mean
     sample_index: the index of the sample for the log likelihood
    

     return value: a float corresponding to the log likelihood of the observations given the mean/variance for the full set of probes and cell composition for the individual specified in sample_index
     assumption: warning that this function does not check that probe_index is within the bounds of the obs_meth_mat, mean_mat or var_mat matrices, nor that sample_index
     is within the bounds of the obs_meth_mat or infer_weight matrices

   */
  float logLike=0;
  for (int probe_index = 0; probe_index<nprobe; probe_index++) {
    logLike+=logLikelihoodPerSampleProbe(obs_meth_mat,  k, infer_weight, mean_mat, var_mat, sample_index, probe_index, logTable);
  }
  return logLike;
}

float logLikelihoodPerProbe(float **obs_meth_mat,   int nsamp, int k, float **infer_weight, float **mean_mat, float **var_mat, int probe_index, float *logTable) {
  /*
    Takes log likelihood for given sample, given the observations and the current weights and means
 input arguments:
     obs_meth_mat: a matrix of observed values with the samples corresponding to rows and the probes corresponding to columns
     k: The number of cell types
     nsamp: The number of samples (rows) in obs_meth_mat and infer_weight
     infer_weight: a matrix of weights for individuals, with each row being one individual (sample), and each column corresponding to one of k cell types
     mean_mat: a matrix of the mean methylation values for each probe, with each row being one of k cell types and each column being one of the probes
     var_mat: a matrix of mean methylation values, same dimension as mean
     sample_index: the index of the sample for the log likelihood


     return value: a float corresponding to the log likelihood of the observations given the mean methylation values for the specified probe and mean methylation values for the
     set of individuals
     assumption: warning that this function does not check that  that sample_index
     is within the bounds of the obs_meth_mat or infer_weight matrices.


  */

 float logLike=0;
  for (int sample_index=0; sample_index<nsamp; sample_index++) {
    logLike+=logLikelihoodPerSampleProbe(obs_meth_mat,  k, infer_weight,mean_mat, var_mat,sample_index, probe_index, logTable);
  }

  return logLike;
}


float logLikelihood(float **obs_meth_mat, int nprobe, int nsamp, int k, float **infer_weight,    float **mean_mat, float **var_mat, float *logTable) {
/*
    Takes log likelihood for given sample, given the observations and the current weights and means
 input arguments:
     obs_meth_mat: a matrix of observed values with the samples corresponding to rows and the probes corresponding to columns
     k: The number of cell types
     nprobe: total number of probes (columns) in obs_meth_mat and mean_mat
     nsamp: The number of samples (rows) in obs_meth_mat and infer_weight
     infer_weight: a matrix of weights for individuals, with each row being one individual (sample), and each column corresponding to one of k cell types
     mean_mat: a matrix of the mean methylation values for each probe, with each row being one of k cell types and each column being one of the probes
     var_mat: a matrix of mean methylation values, same dimension as mean
     sample_index: the index of the sample for the log likelihood


     return value: a float corresponding to the log likelihood of the observations given the cell type composition for the full set of individuals probe and mean methylation values for the
     full set of probes
     
  */
  float logLike=0;
  for (int j=0;j<nsamp;j++) {
    logLike+=logLikelihoodPerSample(obs_meth_mat, nprobe, k, infer_weight, mean_mat, var_mat, j, logTable);
  }
  return logLike;
}



void write_files(string filename, int k, int seed, float** infer_weight, int nsamp, int nprobe, string *probe_ids, float** obs_mat, float **infer_mean, float **infer_var, float *obs_mean, float *obs_var, string samplePrefix, string *sample_ids) {

  //File names and opening of three files
  stringstream sstm;
  sstm << filename << ".k" << k  << ".seed" << seed << ".";
  string file_base = sstm.str();
  string weight_file = file_base + "w";
  string mean_file = file_base + "meanvar";
  string resid_file = file_base + "resid"; 
  
  FILE *fw=fopen(weight_file.c_str(),"w");
  FILE *fp=fopen(mean_file.c_str(),"w");
  FILE *fr=fopen(resid_file.c_str(),"w");
  

  // print weights
  for (int s=0;s<nsamp;s++) {
    for (int a=0;a<k;a++) {
      fprintf(fw,"%f\t", infer_weight[s][a]);
    }
    fprintf(fw, "\n");
  }

  //print inferred probe mean for component 1, probe var for component 1, probe mean for component 2... etc and then observed mean and var for that probe
  for(int i=0;i<nprobe;i++) {
    fprintf(fp,"%s ", probe_ids[i].c_str());
    for (int a=0;a<k;a++) {
      fprintf(fp,"%f %f ", infer_mean[a][i], infer_var[a][i]);
    }
    fprintf(fp,"%f %f\n", obs_mean[i], obs_var[i]);
  } 
  

  //print residuals (observed mean - weighted sum of inferred means), each row is one probe each column one sample, includes column header row as originally in input file.
 
  
  fprintf(fr, "%s ", samplePrefix.c_str());
  for(int i=0; i<nsamp;i++) {
    fprintf(fr, "%s ", sample_ids[i].c_str());
  }
  fprintf(fr, "\n");
  for (int i=0;i<nprobe;i++) {
    fprintf(fr, "%s ", probe_ids[i].c_str());
    for (int ii=0;ii<nsamp;ii++) {
      float pred =0.;
      for (int a=0;a<k;a++) {
	
	pred+=infer_weight[ii][a]*infer_mean[a][i];

      }
      
      fprintf(fr, "%f ", obs_mat[ii][i]-pred);
    }
    fprintf(fr, "\n");
  }



  fclose(fr);
  fclose(fw);
  fclose(fp);
}




void processFile(string filename, int nsamp, int colskip, string &samplePrefix, string *sample_ids, string *probe_ids, float **obs_mat, float *obs_mean, float *obs_var, float* obs_min, float* obs_max, int &nprobe, int maxProbe) {
  /* Read input file, stores information about input measurement in various matrices to me used throughout the program
     Input:
     
     filename - path of file location, assume space or tab separated with a header line for sample ids
     nsamp - number of samples measured. 
     colskip - the number of columns at the beginning of each row that contain probe meta data (probe ID, chromosome, location)
     sample_prefix - address for string that will store the prefix on the first line of files (i.e. column headers for the probe metadata such as id, chromosome, location)
     sample_ids - pointer to strings that will each correspond to one sample id
     probe_ids - pointer to strings that will each correspond to one probe id (if colskip is > 1, this means the probe id will correspond to other probe data stored in column 2 or more)
     obs_mat - float pointer for float values to store the observed (methylation, expression, etc ) values to deconvolve
     obs_mean - mean values of observed experimental values
     obs_var - variance value of observed experimental values
     obs_min - min values of observed experimental values
     obs_vax - maximum value of observed experimental values
   
     nprobe - the actual number of probes (rows) in file, note that this is distinct from maxProbe which is user specified
     maxProbe - user specified maximum number of rows to read from filename
     
     note that all pointers must be initialized and sufficient memory allocated prior to calling the function. Input file must be well formatted with exactly colskip metadata values for each probe
     followed by exactly nsamp measurements. Missing values are specified by -1 and otherwise observed values are assumed to be non-negative. Metadata stored in the first colskip columns may
     be of any format.
  */
  //printf("Called %s\n", filename.c_str());
  
  FILE *f=fopen(filename.c_str(),"r");
  char j[200] ;
  /*Get column headers in first line of file that do not correspond to sample ids*/
  for(int i=0; i < colskip; i++) {
    fscanf(f,"%s", j);
    samplePrefix = samplePrefix + j;
    //printf("prefix %s\n", samplePrefix.c_str());
  }


  /* Get sample IDs from first line of file */
  for (int i=0;i<nsamp;i++) 
    {
      j[0]=0;
      fscanf(f,"%s",j);
      if (j[0]) {
    	  string sample_id(j);
    	  sample_ids[i] = sample_id;
      }
      else 
	break;
    }
  

  // read in data, one probe per row, with colskip columns of probe metadata before data starts. Space or tab separated
  // obtain min, max, mean and var during this pass as well (single pass variance adapted from wikipedia)
  nprobe=0;
  while (fscanf(f,"%s", j) != EOF) {  
    string probe_id(j);     
   
    for (int i = 1 ; i <  colskip; i++) {

      fscanf(f, "%s", j);
      probe_id += j;
    }

    probe_ids[nprobe] = probe_id;
    
    float line_mean = 0;
    float M2 = 0;
    int non_Na_count = 0;
    float line_min = FLT_MAX;
    float line_max = FLT_MIN;
      

    
    for (int i=0;i<nsamp;i++) {
      fscanf(f,"%s",j);
      
      if(strcmp("NA", j)==0) {
	obs_mat[i][nprobe] = -1.0;
      }
      else {
	
	float obs = atof(j);
	if(obs < line_min) {
	  line_min = obs;
	}
	else if(obs > line_max) {
	  line_max = obs;
	}
	obs_mat[i][nprobe] = obs;
	float delta = obs - line_mean;
	non_Na_count++;
	line_mean += delta/non_Na_count;
	M2 += delta *(obs-line_mean);
      }
    }
    obs_var[nprobe] = M2 / (non_Na_count - 1);
    
    obs_mean[nprobe] = line_mean;

    obs_min[nprobe] = line_min;
    obs_max[nprobe] = line_max;
    //    printf("observed mean var %f %f %f %f\n", obs_mean[nprobe], obs_var[nprobe], obs_min[nprobe], obs_max[nprobe]);
    nprobe++;
    if (nprobe>=maxProbe) break;
  } // end of fscanf


}
 
  



void initialize_weight(float **weight, int k, int nsamp, float* min_weights, float* max_weights) {
  /* Initialize weight matrix
     with rows corresponding to nsamp individuals and columns to k cell types. For now use random values that sum up to one 
     in future integrate ways of using prior knowledge
     note that weight matrix must be allocated prior to calling function. 

   */
  for (int i=0;i<k;i++) {
    for (int j=0;j<nsamp;j++) weight[j][i]=(float)rand()/RAND_MAX;
  }
  for (int j=0;j<nsamp;j++) {
    float s=0;
    for (int i=0;i<k;i++) {
      s+=weight[j][i];
    }
    for (int i=0;i<k;i++) {
      weight[j][i]/=s;    //printf("original weight: %f\n",weight[j][i]);
      if(weight[j][i] < min_weights[j]) {
	weight[j][i] = min_weights[j];
      }
      else if (weight[j][i] > max_weights[j]) {
	weight[j][i] = min_weights[j];
      }
    }

  }
}
void initialize_meanvar(float **infer_mean, float** infer_var, float* mins, float* maxes, float* obs_mean, float* obs_var, int k, int nprobe) {
  /* Initialize mean matrix (infer_mean) and variance matrix (infer_var) 
currently using a model where a random number between -0.05 and 0.05 is added to observed mean for each probe
observed variance in each probe is also slightly perturbed
Input:
infer_mean matrix of values to store the mean, must be allocated prior to calling with k rows and nprobe columns
infer_var same as infer_mean but for variances to be inferred
mins the minimum value of each mean, any random initial value that goes below this will be automatically set to this
maxes same as mins but for the max!
obs_mean observed mean of each probe
obs_var observed variance of each probe
k - number of cell types
nprobe - number of probes in study.

Note: infer_mean and infer_var matrices must have memory allocated before calling initialize_meanvar
  */
  for (int i=0;i<k;i++) {


    infer_mean[i]=new float[nprobe];
    infer_var[i]=new float[nprobe];

  }
  for(int cell = 0; cell < k; cell++) {

	  for(int probe_index = 0; probe_index < nprobe; probe_index++) {
    
    	infer_mean[cell][probe_index] = obs_mean[probe_index] + ((float)rand()/RAND_MAX-0.5)*0.1;
    	if(infer_mean[cell][probe_index] < mins[probe_index]) {
	  infer_mean[cell][probe_index] = mins[probe_index];
	}
	else if (infer_mean[cell][probe_index] > maxes[probe_index]) {
	  infer_mean[cell][probe_index] = maxes[probe_index] ;
	}
	infer_var[cell][probe_index]=obs_var[probe_index]+((float)rand()/RAND_MAX-obs_var[probe_index]*0.5)*0.1;
	//printf(" inferred var %d %d %f vs inferred mean %f obs %f %f %f\n", cell, probe_index, infer_var[cell][probe_index], infer_mean[cell][probe_index], obs_mean[probe_index], mins[probe_index], maxes[probe_index]);
	  }
  }
}
 
  

   
void update_weights(float **infer_weight, float* min_weights, float* max_weights, float **obs_mat, float **infer_mean, float **infer_var, int *iters_unconsidered_ind, int* iters_unchanged_ind, int max_unchanged_ind, int max_unconsidered_ind, int k, int nsamp, int nprobe, float *logTable ) {
  /*This is one step in our iterative cell type deconvolution procedure. For each sample, the current log likelihood will be calculated, weights randomly perturbed and for
each perturbation the change accepted if it increases the sample log likelihood

input: 
infer_weight - matrix of weights that will be perturbed and optimized. Each row is one sample each column one cell type
max_weights - maximum proportion allowable for each sample's weights
min_weights - minimum proporation allowable for each sample's weights
obs_mat - observed experimental values
infer_mean - the currently inferred mean values (no values are changed in this function but is needed to determine log likelihood)
infer_var - same as infer_mean but for var
iters_unconsidered_ind - if a sample's weights have not changed for a certain number of iterations, it will not be considered for a certain number of iterations. this variable keeps track of how many iterations have passed since a given sample's weights have been considered for perturbation
iters_unchanged_ind - this keeps track of how many iterations an individual has been considered but no proposed change accepted
max_unchanged_ind - the number of iterations a sample has to be considered with no change accepted before it "takes a break" and is not considered
max_unconsidered_ind - the number of iterations a sample will "take a break" before being considered again
k - number of cell types
nsamp - number of samples
nprobe - number of probes
   */
  float *old_weight = new float [k]();
  for (int s=0;s<nsamp;s++) {
      // pick a component at random and adjust its weight                                                 
      //      printf("Adjusting weight for sample %d\n",s);
      int unchanged = 1;
      if(iters_unchanged_ind[s] > max_unchanged_ind) {
        if(++iters_unconsidered_ind[s] > max_unconsidered_ind) {
          iters_unconsidered_ind[s] = 0;
          iters_unchanged_ind[s] = 0;
        }
        else {
          continue;
        }
      }

      for (int celltype_iter = 0; celltype_iter<k; celltype_iter++) { 
    	  old_weight[celltype_iter]=infer_weight[s][celltype_iter];
      }

      
      float sampleLogLike = logLikelihoodPerSample(obs_mat, nprobe, k, infer_weight, infer_mean, infer_var, s, logTable);
      for (int runnum=0;runnum<k*3;runnum++) {
	

	
	
	

	
	int randomk=rand()%k;
	float sumrest  = 1.0 - infer_weight[s][randomk];
	// revise weights                                                                                   
	infer_weight[s][randomk] *=  1+0.2*((float)rand()/RAND_MAX-0.5);
	//	printf(" Revised %f %f from %f %f\n", infer_weight[s][0], infer_weight[s][1], old_weight[0], old_weight[1]);

	if(infer_weight[s][randomk] / (infer_weight[s][randomk] + sumrest) < min_weights[s]) {
	  infer_weight[s][randomk] = sumrest / (1/min_weights[s] - 1);
	}
	else if (infer_weight[s][randomk] / (infer_weight[s][randomk] + sumrest) > max_weights[s]) {
	  infer_weight[s][randomk] = sumrest / (1/max_weights[s] - 1);
	}
        
	float sw=0;
	for (int a=0;a<k;a++) {
	  sw+=infer_weight[s][a];
	}
	for (int a=0;a<k;a++) {
	  infer_weight[s][a] /=sw;
	}
	
	//	printf(" Revised %f %f from %f %f\n", infer_weight[s][0], infer_weight[s][1], old_weight[0], old_weight[1]);

	
      

      
	//float logLikelihoodPerSample(float **obs_meth_mat, int nprobe, int k, float **infer_weight,  float **mean_mat, float **var_mat, int sample_index) {

	float newSampleLogLike=logLikelihoodPerSample(obs_mat, nprobe, k, infer_weight, infer_mean, infer_var, s, logTable);
      
	if (newSampleLogLike>sampleLogLike) {
	  //printf("changed !\n");
	  unchanged = 0;
	  sampleLogLike=newSampleLogLike;
	  for (int celltype_iter = 0; celltype_iter<k; celltype_iter++) { 
	    old_weight[celltype_iter]=infer_weight[s][celltype_iter];
	  }
	}
	else {
	  for(int celltype_iter = 0; celltype_iter < k; celltype_iter++) {
	    infer_weight[s][celltype_iter] = old_weight[celltype_iter];
	  }
	}
	
      }
      iters_unchanged_ind[s] += unchanged;
  }
}

void update_meanvar(float **obs_mat, float **infer_mean, float **infer_var,  int nprobe, int k, int *iters_unchanged, int *iters_unconsidered, int max_unchanged, int max_unconsidered, float *max_means, float *min_means, float **infer_weight, float *obs_min, float *obs_max, int nsamp, float *logTable) {
  /*This is the other step in the iterative process. In this case for each probe a randomly selected one of the k cell types is selected and its mean and variance perturbed
    the change is accepted if this led to an improvement */



  for (int i=0;i<nprobe;i++) {
    float diff = obs_max[i] - obs_min[i];
    int unchanged = 1;
    if(iters_unchanged[i] > max_unchanged) {
      if(++iters_unconsidered[i] > max_unconsidered) {
	iters_unconsidered[i] = 0;
	iters_unchanged[i] = 0;
      }
      else {
	continue;
      }
    }
    float probeLogLike = logLikelihoodPerProbe(obs_mat, nsamp, k, infer_weight, infer_mean, infer_var,i, logTable);
    for (int celltype_iter=0;celltype_iter<k*3;celltype_iter++) {
      // pick a class randomly
      int randomk=rand()%k;
      float oldinfer_mean=infer_mean[randomk][i];
      float oldinfer_var=infer_var[randomk][i];
      //infer_mean[randomk][i]*=1+0.1*((float)rand()/RAND_MAX-0.5);
      //	infer_var[randomk][i]*=1+0.1*((float)rand()/RAND_MAX-0.5);
	
      float ran1 = (float)rand()/RAND_MAX;
      float ran2 = (float)rand()/RAND_MAX;
	
      infer_mean[randomk][i] *= 1 + 0.25*(ran1*diff-diff/2)/diff;
      infer_var[randomk][i]   *= 1 + 0.25*(ran2*diff-diff/2)/diff;
      if (infer_mean[randomk][i]>max_means[i]) {infer_mean[randomk][i]=max_means[i]; }
      if (infer_mean[randomk][i]<min_means[i]) {infer_mean[randomk][i]=min_means[i];        }
      
      float newProbeLogLike=logLikelihoodPerProbe(obs_mat, nsamp, k, infer_weight, infer_mean, infer_var,i, logTable);
      if (newProbeLogLike>probeLogLike) {

	probeLogLike=newProbeLogLike;
	  unchanged = 0;
      }
      else {

	infer_mean[randomk][i]=oldinfer_mean;
	infer_var[randomk][i]=oldinfer_var;
      }
      
    }
    iters_unchanged[i] += unchanged;
  }



  
}


