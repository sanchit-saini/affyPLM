/*********************************************************************
 **
 ** file: rlm_se.c
 **
 ** Aim: implement computation of standard errors for robust linear models.
 **
 ** Copyright (C) 2003 Ben Bolstad
 **
 ** created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
 ** 
 ** created on: Jan 27, 2003
 **
 ** Last modified: Feb 11, 2003
 **
 ** We provide implemementations of all three methods of computing standard errors
 ** given in Huber (1981) Robust Statistic. Wiley. In particular equations
 ** (6.5), (6.6) and (6.7). Finally we will implement the standard errors in terms
 ** of the standard errors from a weighted linear regression. Note that Huber strongly advises
 ** against the use of this last method as it is "not robust in general". We provide it to agree 
 ** with previous implementations in R code.
 **  
 ** In particular we implement functions
 ** RLM_SE_Method_1 (6.5)
 ** RLM_SE_Method_2 (6.6)
 ** RLM_SE_Method_3 (6.7)
 ** RLM_SE_Method_4 (weighted regression method)
 **
 **
 ** History
 **
 ** Jan 27, 2003 - Initial version using weighted least squares type se
 ** Jan 28, 2003 - More work, try to get huber recommendations working.
 **                Realize that the W matrix is sometimes not of full 
 **                rank. to get around this problem recommend using
 **                generalized inverse (compute using svd). <-- TO DO
 **                Add in better checking to the inversion routines (Choleski).
 ** Jan 31, 2003 - make Choleski inverse take parameter to return only
 **                upper triangle of inverse matrix
 **                Actually make se routines check that nothing bad happens 
 **                when taking the inverse.
 ** Feb 02, 2003 - Test if choleski failed, if did use an svd to compute 
 **                a generalize invese and use this instead.
 **                code for testing failure to se.type 2,3 added
 **                SVD approach to inverses added and tested.
 ** Feb 08, 2003 - Move a code block which will improve efficiency for
 **                se.type = 4
 **                TO BE DONE: replace linpack routines with Lapack routines to further improve speed
 ** Feb 10, 2003 - Change sumpsi, sumpsi2 to user option 2 (the actual psi function) in psi_huber.
 **                this will fix a bug in se.type=1,2,3. we would agree completely
 **                with summary.rlm in R except that we estimate the scale parameter using
 **                residuals from the final fit and R uses a scale from one step back.
 **                A mechanism has been added to switch between LAPACK and LINPACK with the default 
 **                being LAPACK.
 ** Feb 11, 2003 - Make LINPACK DEFAULT, comment out LAPACK in chol, svd routines, solves linking 
 **                problems on windows, will get fixed later
 ** Feb 24, 2003 - comment out a few declared but unused variables, some of these might be used later
 ** Mar 28, 2003 - uncomment LAPACK code. LAPACK will be default using R 1.7.0 and later, which
 **                will be a requirement for AffyExtensions 0.5-1 and later
 ** Jun 11, 2003 - Modify Standard error routines to handle different psi_fns.
 ** Jul 23, 2003 - remove the warning about moduleCdynload by including appropriate header file
 ** Sep 06, 2003 - now we return the whole variance covariance matrix and residual SE with appropriate
 **                DF
 ** Sep 07, 2003 - variance matrix from method == 4 now returned
 ** Sep 08, 2003 - variance matrix from method == 1, 2, 3 returned
 ** Sep 13, 2003 - copy only upper triangle of variance matrix into output.
 **                Also if the variance matrix is the NULL pointer don't store anything
 **                at all.
 **                
 ********************************************************************/

#include "rlm.h"
#include "rlm_se.h"
#include "psi_fns.h"

#include <R_ext/Rdynload.h>
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/********************************************************************
 **
 ** two static global variables to define if lapack library loaded
 ** and which library to use. 0 is LINPACK, 1 is LAPACK
 **
 ** Kind of bad, but will do it anyway
 **
 ********************************************************************/

static int Lapack_initialized = 0;
static int use_lapack = 1;



/********************************************************************
 **
 ** external declarations for Choleski routines (LINPACK)
 **
 **
 *******************************************************************/

extern int dpofa_(double *x, int *lda, int *n, int *j);
extern int dpodi_(double *x, int *lda, int *n, double *d, int *j);


/********************************************************************
 **
 ** external declarations for Choleski routines (LAPACK)
 **
 **
 *******************************************************************/

extern int dpotrf_(const char *uplo, const int *n, double* a, const int *lda, int *info);
extern int dpotri_(const char *uplo, const int *n, double* a, const int *lda, int *info);


/*****************************************************************
 * svd routine - LINPACK
 *****************************************************************/

extern int dsvdc_(double *x, int *ldx, int *n, int *p, double *s, double *e, double *u, int *ldu,
		 double *v, int *ldv, double *work, int *job, int *info);

/*****************************************************************
 * svd routine - LAPACK
 *****************************************************************/

extern int dgesdd_(const char *jobz,
                      const int *m, const int *n,
                      double *a, const int *lda, double *s,
                      double *u, const int *ldu,
                      double *vt, const int *ldvt,
                      double *work, const int *lwork, int *iwork, int *info);




/*********************************************************************
 **
 ** static void Lapack_Init(void)
 **
 ** this function loads the Lapack library if it has not already been
 ** loaded and sets the use_lapack variable to 1 so that LAPACK
 ** is used (for Choleski and SVD routines)
 **
 **
 **
 ********************************************************************/

static void Lapack_Init(void)
{
    int res = moduleCdynload("lapack", 0, 1);
    Lapack_initialized = -1;
    if(!res) return;

    /* Initializing LAPACK */
    use_lapack = 1;
    Lapack_initialized = 1;
    return;
}




/*********************************************************************
 **
 ** int Choleski_decompose(double *X, double *L, int n)
 **
 ** double *X - a square matrix 
 ** double *L - space already allocated to store Cholesky decomposition
 ** int n - matrix dimension
 ** int lapack - if 0 use LINPACK otherwise LAPACK
 **
 ** RETURNS integer code indicating success or failure 0 means no error, 
 **      non zero indicates failure ie not positive definite
 **
 ** this function choleski decomposes a positive definite symmetric matrix,
 ** on output L will contain the Choleski decomposition in the upper triangle
 **
 **
 *********************************************************************/


int Choleski_decompose(double *X, double *L, int n, int lapack){
  int i,j,error_code;
  char upper = 'U';

  for (i=0; i < n; i++){
    for (j=0; j < n; j++){
      if (i > j)
	L[j*n+i] = 0.0;
      else {
	L[j*n+i] = X[j*n + i];
      }
    }
  }
  if (!lapack){
    dpofa_(L,&n,&n,&error_code);
  } else {
    dpotrf_(&upper,&n,L,&n,&error_code);
  }
    


  return error_code;
}

/***********************************************************************
 **
 ** int Choleski_2_inverse(double *X, double *Xinv, int n)
 **
 ** double *X - matrix containing choleski decomposition in upper triangle
 ** double *Xinv - on output will contain the inverse
 ** int n - dimension of matrix
 ** int upperonly - if non zero return only the upper triangle of the inverse.
 ** int lapack - use LINPACK if 0 otherwise LAPACK
 **
 ** RETURNS integer code, indicating success 0  or error (non zero) 
 **
 ** this function uses the choleski decomposition of a 
 ** matrix to compute the inverse of a matrix.
 ** typically it would be used in conjunction with the choleski_decompose
 ** function above.
 **
 **
 **********************************************************************/

int Choleski_2_inverse(double *X, double *Xinv, int n,int upperonly, int lapack){
  
  int i,j ,error_code=0,inverseonly;
  double d =0.0;
  char upper = 'U';
  
  for (i=0; i < n; i++){ 
    /* check for a zero or close to zero diagonal element */ 
    if(fabs(X[i*n+ i]) < 1e-06){
      error_code = 1;
      return error_code;
    }

    for (j=i; j < n; j++){
      Xinv[j*n + i] = X[j*n + i];
    }
  }

  inverseonly = 1;
  if (!lapack){
    dpodi_(Xinv,&n,&n,&d,&inverseonly);
  } else {
    dpotri_(&upper,&n,Xinv,&n,&error_code);
  }
  if (!upperonly){
    for (i=0; i < n; i++){
      for (j=0; j <= i-1; j++){
	Xinv[j*n+i] = Xinv[i*n+j];
      }
    }
  }
  return error_code;

}


/***********************************************************************
 **
 ** int Choleski_inverse(double *X, double *Xinv, double *work, int n)
 **
 ** double *X - matrix containing choleski decomposition in upper triangle
 ** double *Xinv - on output will contain the inverse
 ** double *work - working space n*n dimension
 ** int n - dimension of matrix
 ** int upperonly - if non zero return only upper triangle of inverse.
 **
 ** RETURNS integer code, indicating success 0  or error (non zero) 
 **
 ** This function will compute the inverse of a positive definite symmetric
 ** matrix using choleski decomposition.
 **
 **********************************************************************/

int Choleski_inverse(double *X, double *Xinv, double *work, int n, int upperonly){

  int error_code;
  
  error_code = Choleski_decompose(X, work, n,use_lapack);
  if (!error_code){
    error_code = Choleski_2_inverse(work, Xinv, n,upperonly,use_lapack);
  }
  return error_code;

}



/***************************************************************
 **
 ** int SVD_compute()
 **
 **
 ** Computes the singular value decomposition of a matrix. Current
 ** implemtnation uses a linpack routine, but this will later be transitioned
 ** to a lapack routine (which is faster)
 **
 ***************************************************************/

int SVD_compute(double *X, int n, double *s, double *u, double *v,int lapack){
  
  int i,j, error_code;
  int lwork = 7*n*n + 4*n;
  int job = 21;
  char jobz = 'A';
  double *Xcopy= calloc(n*n,sizeof(double));              //Calloc(n*n,double);
  double *e =    calloc(n,sizeof(double));              // Calloc(n,double);
  double *work =  calloc(n,sizeof(double));             //Calloc(n,double);
  double *work2 =  calloc(lwork,sizeof(double));
  int *iwork = calloc(8*n,sizeof(int));


  for (i=0; i < n; i++){
    for (j=0; j < n; j++){
      Xcopy[j*n + i] = X[j*n+i];
    }
  }
  if (!lapack){
    dsvdc_(Xcopy,&n,&n,&n,s,e,u,&n,v,&n,work,&job,&error_code);
  } else {
    dgesdd_(&jobz,&n,&n,Xcopy,&n,s,u,&n,v,&n,work2,&lwork,iwork,&error_code);
  }
    

  free(iwork);
  free(work2);
  free(work);
  free(e);
  free(Xcopy);
  
  return error_code;

}

/***************************************************************
 **
 ** int SVD_2_inverse(double *Xinv, int n, double *s, double *u, double *v,int lapack)
 **
 ** double *Xinv - on exit contains the inverse
 ** int n - Xinv is n by n
 ** double *s - SVD components length n
 ** double *u - SVD components n by n
 ** double *v - SVD components n by n
 ** int lapack - non zero if the decomposition was done by a LAPACK routine (implies v is the transpose)
 **
 ** given a Singular value decomposition of a matrix compute
 ** the generalized inverse.
 **
 ***************************************************************/

int SVD_2_inverse(double *Xinv, int n, double *s, double *u, double *v,int lapack){
  double tolerance = 1e-7; //1.490116e-08;
  int i,j,k;
  int nonzero =n;
  
  
  for (i = 0; i < n; i++){
    if (s[i] < tolerance*s[0]){
      nonzero = i;
      //printf("nonzero %d",i);
      break;
    }
  }


  /* for all columns where $d is not to small do */
  /*  svdu$v %*% (t(svdu$u)* 1/svdu$d); */
    
  for (i = 0; i < n; i++){
    for (j = 0; j < nonzero; j++){
      u[j*n + i] = u[j*n+i] * 1.0/s[j];
    }
  }
  if (!lapack){
    for (i = 0; i < n; i++){
      for (j = 0; j < n; j++){
	Xinv[j*n+i] = 0.0;
	for (k=0; k < nonzero; k++){
	  Xinv[j*n+i]+= v[k*n+i] * u[k*n+j];
	}
      }
    }
  } else {
    /* lapack so v is transposed */
    for (i = 0; i < n; i++){
      for (j = 0; j < n; j++){
	Xinv[j*n+i] = 0.0;
	for (k=0; k < nonzero; k++){
	  Xinv[j*n+i]+= v[i*n+k] * u[k*n+j];
	}
      }
    }
  }



  return 0;
}


/***************************************************************
 **
 ** int SVD_inverse(double *X, double *Xinv, int n)
 **
 ** double *X -  the matrix to be inverted
 ** double *Xinv - on exit contains inverse
 ** int n - X and Xinv are n by n
 **
 **
 ** given an n by n matrix compute its inverse by singular value decomposition
 ** this is particularly useful when dealling with singular or near singular
 ** matrices since if the regular inverse can not be computed a generalized
 ** inverse will be returned instead, rather than erroneous results.
 **
 ** Note that  we will use linpack routine for SVD at some point this will
 ** be transitioned to a lapack routine. (This has now been done).
 **
 **
 **
 **************************************************************/

int SVD_inverse(double *X, double *Xinv, int n){

  int error_code=0;
  double *s = Calloc(n+1,double);
  double *v = Calloc(n*n,double);
  double *u = Calloc(n*n,double);

  error_code = SVD_compute(X, n, s, u, v,use_lapack);
  SVD_2_inverse(Xinv,n, s, u, v,use_lapack);

  return error_code;

  Free(s);
  Free(v);
  Free(u);
}


/* void R_SVD_compute(double *X, int *n, double *s, double *u, double *v){
  SVD_compute(X, *n,s,u, v);
  } */


void R_SVD_inverse(double *X, double *Xinv, int *n){
  SVD_inverse(X, Xinv,*n);
}








/*************************************************************************
 **
 ** void RLM_SE_Method_1(double residvar, double *XTX, int p, double *se_estimates)
 **
 ** double residvar - residual variance estimate
 ** double *XTX - t(Design matrix)%*% Design Matrix
 ** double p - number of parameters
 ** double *se_estimates - on output contains standard error estimates for each of
 **                        the parametes
 **
 ** this function computes the parameter standard errors using the first
 ** method described in Huber (1981)
 ** 
 ** ie k^2 (sum psi^2/(n-p))/(sum psi'/n)^2 *(XtX)^(-1)
 **
 **
 ************************************************************************/


void RLM_SE_Method_1(double residvar, double *XTX, int p, double *se_estimates,double *varcov){
  int i,j;
  double *XTXinv = Calloc(p*p,double);
  double *work = Calloc(p*p,double);

  if (!Choleski_inverse(XTX,XTXinv,work,p,1)){
    for (i =0; i < p; i++){
      se_estimates[i] = sqrt(residvar*XTXinv[i*p + i]);
    }
  } else {
    printf("Singular matrix in SE inverse calculation");    
  }    


  if (varcov != NULL)
    for (i =0; i < p; i++){
      for (j = i; j < p; j++){
	varcov[j*p +i]= residvar*XTXinv[j*p +i];
      }
    }
  
  Free(work);
  Free(XTXinv);
}


/*************************************************************************
 **
 ** void RLM_SE_Method_2(double residvar, double *W, int p, double *se_estimates)
 **
 ** double residvar - residual variance estimate
 ** double *XTX - t(Design matrix)%*% Design Matrix
 ** double p - number of parameters
 ** double *se_estimates - on output contains standard error estimates for each of
 **                        the parametes
 **
 ** this function computes the parameter standard errors using the second
 ** method described in Huber (1981)
 ** 
 ** ie K*(sum psi^2/(n-p))/(sum psi'/n) *(W)^(-1)
 **
 **
 ************************************************************************/

void RLM_SE_Method_2(double residvar, double *W, int p, double *se_estimates,double *varcov){
  int i,j; /* l,k; */
  double *Winv = Calloc(p*p,double);
  double *work = Calloc(p*p,double);

  if (!Choleski_inverse(W,Winv,work,p,1)){
    for (i =0; i < p; i++){
      se_estimates[i] = sqrt(residvar*Winv[i*p + i]);
      /* printf("%f ", se_estimates[i]); */
    }
  } else {
    //printf("Using a G-inverse\n");
    SVD_inverse(W, Winv,p);
    for (i =0; i < p; i++){
      se_estimates[i] = sqrt(residvar*Winv[i*p + i]);
      /* printf("%f ", se_estimates[i]); */
    }
  }

  if (varcov != NULL)
    for (i =0; i < p; i++){
      for (j = i; j < p; j++){
	varcov[j*p +i]= residvar*Winv[j*p +i];
      }
    }
  
  
  Free(work);
  Free(Winv);

}

/*************************************************************************
 **
 ** void RLM_SE_Method_3(double residvar, double *XTX, double *W, int p, double *se_estimates)
 **
 ** double residvar - residual variance estimate
 ** double *XTX - t(Design matrix)%*% Design Matrix
 ** double p - number of parameters
 ** double *se_estimates - on output contains standard error estimates for each of
 **                        the parametes
 **
 ** this function computes the parameter standard errors using the third
 ** method described in Huber (1981)
 ** 
 ** ie 1/(K)*(sum psi^2/(n-p))*(W)^(-1)(XtX)W^(-1)
 **
 **
 ************************************************************************/

int RLM_SE_Method_3(double residvar, double *XTX, double *W, int p, double *se_estimates,double *varcov){
  int i,j,k;   /* l; */
  int rv;

  double *Winv = Calloc(p*p,double);
  double *work = Calloc(p*p,double);
  

  /***************** 

  double *Wcopy = Calloc(p*p,double);

  
  for (i=0; i <p; i++){
    for (j=0; j < p; j++){
      Wcopy[j*p + i] = W[j*p+i];
    }
    } **********************/
  
  if(Choleski_inverse(W,Winv,work,p,1)){
    SVD_inverse(W, Winv,p);
  }
  
  /*** want W^(-1)*(XtX)*W^(-1) ***/

  /*** first Winv*(XtX) ***/

  for (i=0; i <p; i++){
    for (j=0; j < p; j++){
      work[j*p + i] = 0.0;
      for (k = 0; k < p; k++){
	work[j*p+i]+= Winv[k*p +i] * XTX[j*p + k];
      }
    }
  }
 
  /* now again by W^(-1) */
  
   for (i=0; i <p; i++){
    for (j=0; j < p; j++){
      W[j*p + i] =0.0;
      for (k = 0; k < p; k++){
	W[j*p+i]+= work[k*p +i] * Winv[j*p + k];
      }
    }
   }
  
   for (i =0; i < p; i++){
     se_estimates[i] = sqrt(residvar*W[i*p + i]);
     // printf("%f ", se_estimates[i]);
   }
   
   rv = 0;
   if (varcov != NULL)
     for (i =0; i < p; i++){
       for (j = i; j < p; j++){
	 varcov[j*p +i]= residvar*W[j*p +i];
       }
   }
  
  Free(work);
  Free(Winv);

  return rv;

}



/*********************************************************************
 **
 ** void rlm_compute_se(double *X,double *Y, int n, int p, double *beta, double *resids,double *weights,double *se_estimates, int method)
 **
 ** given a robust linear model fit to data, compute the standard errors of the parameter estimates
 **
 ** double *X
 ** double *Y
 ** int n
 ** int p
 ** double *beta
 ** double *resids
 ** double *weights
 ** double *se_estimates
 ** int method
 **
 **
 **
 **
 **
 ** Note that we compute Kappa using a simplification for Huber Psi
 **
 ** ie Kappa = 1 + p/n* var(psi')/(E psi')^2
 **
 ** simplifies to
 **
 ** Kappa = 1 + p*(1-m)/(n*m)
 **
 ** where m is the proportion of psi' = 1 (ie the unwinsorized observations)
 **
 ** note that W_jk = sum_i (psi'(r_i)*x_{ij}*x_{ik}) ie a weighted form of XTX
 ** 
 **  
 ** 
 **
 **
 *********************************************************************/


void rlm_compute_se(double *X,double *Y, int n, int p, double *beta, double *resids,double *weights,double *se_estimates, double *varcov, double *residSE, int method,double (* PsiFn)(double, double, int), double psi_k){
  
  int i,j,k; /* counter/indexing variables */
  double k1 = psi_k;   /*  was 1.345; */
  double sumpsi2=0.0;  /* sum of psi(r_i)^2 */
  /*  double sumpsi=0.0; */
  double sumderivpsi=0.0; /* sum of psi'(r_i) */
  double Kappa=0.0;      /* A correction factor */
  double scale=0.0;

  double *XTX = Calloc(p*p,double);
  double *W = Calloc(p*p,double);
  double *work = Calloc(p*p,double);
  double RMSEw = 0.0;
  double vs=0.0,m;  /* varpsiprime=0.0; */

  /* Initialize Lapack library */
  if(!Lapack_initialized) Lapack_Init();

  if (method == 4){
    for (i=0; i < n; i++){
      RMSEw+= weights[i]*resids[i]*resids[i];
    }
    
    RMSEw = sqrt(RMSEw/(double)(n-p));

    residSE[0] =  RMSEw;

    for (j =0; j < p;j++){
      for (k=0; k < p; k++){
	W[k*p + j] = 0.0;
	for (i=0; i < n; i++){
	  W[k*p + j]+=  weights[i]*X[j*n +i]*X[k*n + i];
	}
      }
    }
    if (!Choleski_inverse(W, XTX, work, p,1)){
      
      for (i =0; i < p; i++){
	se_estimates[i] = RMSEw*sqrt(XTX[i*p + i]);
      }
    } else {
      printf("Singular matrix in SE inverse: Method 4\n");
      
    }


    if (varcov != NULL)
      for (i = 0; i < p; i++)
	for (j = i; j < p; j++)
	  varcov[j*p + i] =  RMSEw*RMSEw*XTX[j*p + i];
  } else {

    scale = med_abs(resids,n)/0.6745;
    
    residSE[0] =  scale;
    
    /* compute most of what we will need to do each of the different standard error methods */
    for (i =0; i < n; i++){
      sumpsi2+= PsiFn(resids[i]/scale,k1,2)*PsiFn(resids[i]/scale,k1,2); 
      /* sumpsi += psi_huber(resids[i]/scale,k1,2); */
      sumderivpsi+= PsiFn(resids[i]/scale,k1,1);
    }
    
    m = (sumderivpsi/(double) n);
    Kappa = 1.0 + (double)p/(double)n * (1.0-m)/(m);
    
    /* prepare XtX and W matrices */
  
    for (j =0; j < p;j++){
      for (k=0; k < p; k++){
	for (i = 0; i < n; i++){
	  XTX[k*p+j]+=  X[j*n +i]*X[k*n + i];
	  W[k*p + j]+=  PsiFn(resids[i]/scale,k1,1)*X[j*n +i]*X[k*n + i];	
	}
      }
    }

    if (method==1) {
      Kappa = Kappa*Kappa;
      vs = scale*scale*sumpsi2/(double)(n-p);
      Kappa = Kappa*vs/(m*m);
      RLM_SE_Method_1(Kappa, XTX, p, se_estimates,varcov);
    } else if (method==2){
      vs = scale*scale*sumpsi2/(double)(n-p);
      Kappa = Kappa*vs/m;
      RLM_SE_Method_2(Kappa, W, p, se_estimates,varcov);
      
    } else if (method==3){
      
      vs = scale*scale*sumpsi2/(double)(n-p);
      Kappa = 1.0/Kappa*vs;
      i = RLM_SE_Method_3(Kappa, XTX, W, p, se_estimates,varcov);
      if (i){
	for (i=0; i <n; i++){
	  printf("%2.1f ", PsiFn(resids[i]/scale,k1,1));
	} 
	printf("\n");
      }
    } 
  }
  Free(work);
  Free(XTX);
  Free(W);
}
