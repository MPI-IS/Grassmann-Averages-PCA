// Mex wrapper file for robust PCA


#include "mex.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  // arguments checking
  if (nrhs != 1) 
  {
	  mexErrMsgTxt("One input argument required.");
  }
  else if (nlhs > 1) 
  {
	  mexErrMsgTxt("Too many output arguments.");
  }

  const mxArray* const X = prhs[0];
  assert(X);

  // checking the format of the data
  if (!mxIsDouble(X) && !mxIsSingle(X))
  {
	  mexErrMsgTxt("Unsupported format");
  }

  if(mxIsComplex(X))
  {
    mexErrMsgTxt("Unsupported format (should be scalar)");
  }



  

  // Check the dimensions of inputs
  const size_t rows = mxGetM(X);
  const size_t columns = mxGetN(X);
  const size_t &dimension = columns; // an alias

  // for the library to work, we need to allocate some temporary storage
  // we also allocate the output if given



  plhs[0] = mxCreateDoubleMatrix(dimension , dimension , mxREAL); 

}