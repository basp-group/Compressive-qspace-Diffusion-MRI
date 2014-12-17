


#include "mex.h"  /* Always include this */
#include <math.h>
#include <matrix.h>

/**
 * L1 solver.
 *
 * Usage: [ W ] = fast_compute_weights( X, NEIGH_M, NEIGH_N, MASK) 

*/


 void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{

  int i,j, iaux, iin, iout, id, iv, v, d;

  int Mx, Nx, M_NeighM, N_NeighM, M_NeighN, N_NeighN, M_m, N_m;
  double *X=NULL, *W=NULL, *dir=NULL, *vox=NULL;
  // double tau;
  double w;

  mxArray *pdir=NULL, *pvox=NULL;
  mxLogical *pm=NULL;
  size_t ndir, nvox;


// Check number of arguments. 
  if (nrhs != 4) {
    mexErrMsgIdAndTxt("compute_weights_mex:InvalidInput:nrhs",
          "Require 5 inputs.");
  }
  if (nlhs != 1) {
    mexErrMsgIdAndTxt("compute_weights_mex:InvalidOutput:nlhs",
          "Require one output.");
  }

  // Parse X.
  iin = 0;
  mexPrintf("Parsing the %i input element\n",iin);
  Mx = mxGetM(prhs[iin]); /* Get the dimensions of X*/
  Nx = mxGetN(prhs[iin]);  
  if (!mxIsDouble(prhs[iin]))
    mexErrMsgIdAndTxt("compute_weights_mex:InvalidInput:y",
          "X must be a real 2D double array.");
  X = mxGetPr(prhs[iin]);  /* Get the pointer to the data of X */
  mexPrintf("Matrix X read: %i rows and %i columns\n", Mx, Nx);

  
  // Parse NEIGH_M
  iin++;
  mexPrintf("Parsing the %i input element\n",iin);
  M_NeighM = mxGetM(prhs[iin]); /* Get the dimensions of Neigh_M */
  N_NeighM = mxGetN(prhs[iin]); 
  if (!mxIsCell(prhs[iin]))
    mexErrMsgIdAndTxt("compute_weights_mex:InvalidInput:y",
          "NEIGH_M must be a cell array.");
  
  mexPrintf("Matrix NEIGH_M is of %s class: %i rows and %i columns\n",mxGetClassName(prhs[iin]), M_NeighM, N_NeighM);

  pdir = mxGetCell(prhs[iin],0);  /* Get contents of the zero cell */
  mexPrintf("get cell NEIGH_N done \n");
  ndir= mxGetNumberOfElements(pdir);
  mexPrintf("Direction 0 is of class %s has %i neighbors:\n",mxGetClassName(pdir), ndir);
  if (!mxIsDouble(pdir))
    mexErrMsgIdAndTxt("compute_weights_mex:InvalidInput:y",
          "Contents of cell NEIGH_M must be a double array.");
  dir=mxGetPr(pdir);
  for(i=0;i<ndir;i++)
  	mexPrintf("%f, ",dir[i]);
  mexPrintf("\n");

// Parse NEIGH_N
  iin++;
  mexPrintf("Parsing the %i input element\n",iin);
  M_NeighN = mxGetM(prhs[iin]); /* Get the dimensions of NEIGH_N */
  N_NeighN = mxGetN(prhs[iin]);  
  if (!mxIsCell(prhs[iin]))
    mexErrMsgIdAndTxt("compute_weights_mex:InvalidInput:y",
          "NEIGH_N must be a cell array.");
  
  mexPrintf("Matrix NEIGH_N is of %s class: %i rows and %i columns\n",mxGetClassName(prhs[iin]), M_NeighN, N_NeighN);

  pvox = mxGetCell(prhs[iin],0);  /* TEST: Get contents of the zero cell */
  mexPrintf("get cell NEIGH_N done \n");
  nvox= mxGetNumberOfElements(pvox);
  mexPrintf("Voxel 0 is of class %s has %i neighbors:\n",mxGetClassName(pvox), nvox);
  if (!mxIsDouble(pvox))
    mexErrMsgIdAndTxt("compute_weights_mex:InvalidInput:y",
          "Contents of cell NEIGH_M must be a double array.");
  vox=mxGetPr(pvox);
  for(i=0;i<nvox;i++)
  	mexPrintf("%f, ",vox[i]);
  mexPrintf("\n");

// Parse MASK
  iin++;
  mexPrintf("Parsing the %i input element\n",iin);
  M_m= mxGetM(prhs[iin]); /* Get the dimensions of NEIGH_N */
  N_m = mxGetN(prhs[iin]);  
  if (!mxIsLogical(prhs[iin]))
    mexErrMsgIdAndTxt("compute_weights_mex:InvalidInput:y",
          "MASK must be a logical array.");
  
  mexPrintf("Matrix MASK is of %s class: %i rows and %i columns\n",mxGetClassName(prhs[iin]), M_m, N_m);

  pm = mxGetLogicals(prhs[iin]);  /* TEST: Get data */

  i=1025;
  mexPrintf("get logical array MASK done. TEST: Element %i of the mask is %s.\n", i, (pm[i])?"open":"closed");
  

  iout=0;

  plhs[iout]=mxCreateDoubleMatrix(Mx,Nx,mxREAL); /* Create the output matrix */

  W = mxGetPr(plhs[iout]); /* Create the pointer to the output data */

  // Compute all the weights (loop for every voxel)

  i=0; j=0;

  for(i=0;i<Mx;i++) { // for all directions

    for(j=0;j<Nx;j++) {// for all voxels

  

      pvox = mxGetCell(prhs[2],j);  /* Get contents of the j cell */
      nvox= mxGetNumberOfElements(pvox);
      vox=mxGetPr(pvox);
      
      pdir = mxGetCell(prhs[1],i);  /* Get contents of the i cell */
      ndir= mxGetNumberOfElements(pdir);
      dir=mxGetPr(pdir);
      
      w=0;
      for(id=0;id<ndir;id++){
        for(iv=0;iv<nvox;iv++){

          d=dir[id]-1;
          v=vox[iv]-1;
          w+=X[d+Mx*v]; 
          }
      }

      
      w=w/(double)(nvox);
      W[i + Mx*j]=w;

  }
  } 
 return;


}