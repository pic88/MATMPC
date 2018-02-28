#include "mex.h"
#include "string.h"

#include "blas.h"

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{

    double *Q = mxGetPr( mxGetField(prhs[0], 0, "Q") );
    double *S = mxGetPr( mxGetField(prhs[0], 0, "S") );
    double *R = mxGetPr( mxGetField(prhs[0], 0, "R") );
    double *A = mxGetPr( mxGetField(prhs[0], 0, "A") );
    double *B = mxGetPr( mxGetField(prhs[0], 0, "B") );
    double *Cx = mxGetPr( mxGetField(prhs[0], 0, "Cx") );
    double *Cu = mxGetPr( mxGetField(prhs[0], 0, "Cu") );
    double *CN = mxGetPr( mxGetField(prhs[0], 0, "CN") );
    double *gx = mxGetPr( mxGetField(prhs[0], 0, "gx") );
    double *a = mxGetPr( mxGetField(prhs[0], 0, "a") );
    double *ds0 = mxGetPr( mxGetField(prhs[0], 0, "ds0") );
    
    double *mu_vec = mxGetPr(prhs[2]);
    
    mwSize nx = mxGetScalar( mxGetField(prhs[1], 0, "nx") );
    mwSize nu = mxGetScalar( mxGetField(prhs[1], 0, "nu") );
    mwSize nc = mxGetScalar( mxGetField(prhs[1], 0, "nc") );
    mwSize ncN = mxGetScalar( mxGetField(prhs[1], 0, "ncN") );
    mwSize N = mxGetScalar( mxGetField(prhs[1], 0, "N") );
       
    double *dx = mxGetPr( mxGetField(prhs[0], 0, "dx") );
    double *du = mxGetPr( mxGetField(prhs[0], 0, "du") );
    double *lambda = mxGetPr( mxGetField(prhs[0], 0, "lambda_new") );
    double *mu = mxGetPr( mxGetField(prhs[0], 0, "mu_new") );
    double *muN = mxGetPr( mxGetField(prhs[0], 0, "muN_new") );
   
    memcpy(mu, mu_vec, nc*N*sizeof(double));
    memcpy(muN, mu_vec+N*nc, ncN*sizeof(double));
 
    int i;
        
    char *nTrans = "N", *Trans="T";
    double one_d = 1.0, zero = 0.0;
    mwSignedIndex one_i = 1;
       
    memcpy(dx, ds0, nx*sizeof(double)); 
    
    for (i=0;i<N;i++){        
        memcpy(dx+(i+1)*nx, a+i*nx, nx*sizeof(double));          
        dgemv(nTrans,&nx,&nx,&one_d,A+i*nx*nx,&nx,dx+i*nx,&one_i,&one_d,dx+(i+1)*nx,&one_i);
        dgemv(nTrans,&nx,&nu,&one_d,B+i*nx*nu,&nx,du+i*nu,&one_i,&one_d,dx+(i+1)*nx,&one_i);
    }
    
    memcpy(lambda+N*nx, gx+N*nx, nx*sizeof(double));
    dgemv(nTrans,&nx,&nx,&one_d,Q+N*nx*nx,&nx,dx+N*nx,&one_i,&one_d,lambda+N*nx,&one_i);
    
    if (ncN>0){
        dgemv(Trans,&ncN,&nx,&one_d,CN,&ncN,mu_vec+N*nc,&one_i,&one_d,lambda+N*nx,&one_i);
    }
    for (i=N-1;i>-1;i--){
        memcpy(lambda+i*nx,gx+i*nx, nx*sizeof(double));
        dgemv(nTrans,&nx,&nx,&one_d,Q+i*nx*nx,&nx,dx+i*nx,&one_i,&one_d,lambda+i*nx,&one_i);
        dgemv(nTrans,&nx,&nu,&one_d,S+i*nx*nu,&nx,du+i*nu,&one_i,&one_d,lambda+i*nx,&one_i);
        dgemv(Trans,&nx,&nx,&one_d,A+i*nx*nx,&nx,lambda+(i+1)*nx,&one_i,&one_d,lambda+i*nx,&one_i);
        
        if (nc>0)
            dgemv(Trans,&nc,&nx,&one_d,Cx+i*nc*nx,&nc,mu_vec+i*nc,&one_i,&one_d,lambda+i*nx,&one_i);       
    }
    
}