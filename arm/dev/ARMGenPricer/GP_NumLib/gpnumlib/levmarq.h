#ifndef _INGPNUMLIB_LEVMARQ_H
#define _INGPNUMLIB_LEVMARQ_H

///GP Base
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/ostringstream.h"
#include "gpbase/warningkeeper.h"
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/port.h"
#include "gpbase/gpvector.h"
#include "gpbase/typedef.h"
#include "gpbase/removenagwarning.h"

#include <glob/expt.h>
#include <cmath>
#include <iomanip> /// for setprecision()

#include <glob/armglob.h>

using namespace std;

CC_BEGIN_NAMESPACE( ARM )

/////////////////////////////
//
//		Encapsulation des fonctions pour le Levenberg - Marquardt
//
/////////////////////////////

class ARM_LEVMARQFunc
{
public:
	virtual void operator()(double p[], double hx[], int m, int n, void * adata = 0) const = 0;
};


#ifdef __cplusplus
extern "C" {
#endif

#ifndef MIN
#define MIN(x,y) (x < y ? x : y)
#endif

#ifndef MAX
#define MAX(x,y) (x	>= y ? x : y)
#endif

#ifndef FABS
#define FABS(x) (((x)>=0.0)? (x) : -(x))
#endif

/* work arrays size for LM with & without jacobian, should be multiplied by sizeof(double)
 * or sizeof(float) to be converted to bytes
 */
#define LM_DER_WORKSZ(npar, nmeas) (2*(nmeas) + 4*(npar) + (nmeas)*(npar) + (npar)*(npar))
#define LM_DIF_WORKSZ(npar, nmeas) (3*(nmeas) + 4*(npar) + (nmeas)*(npar) + (npar)*(npar))

#define LM_OPTS_SZ    	 5 /* max(4, 5) */
#define LM_INFO_SZ    	 9
#define LM_INIT_MU    	 1E-05
#define LM_STOP_THRESH	 1E-17
#define LM_DIFF_DELTA    1E-10
#define LM_VERSION       "2.1 (Apr. 2005)"

/* double precision LM, with & without jacobian */
/* unconstrained minimization */

extern int dlevmar_der(
		const ARM_LEVMARQFunc& func,
		const ARM_LEVMARQFunc& fjac,
		double * p, double * x, int m, int n, int itmax, double *opts,
		double *info, double *work, double *covar, void *adata);

extern int dlevmar_dif(
		const ARM_LEVMARQFunc& func,
		double * p, double * x, int m, int n, int itmax, double *opts,
		double *info, double *work, double *covar, void *adata);

/* box-constrained minimization */
extern int dlevmar_bc_der(
		const ARM_LEVMARQFunc& func,
		const ARM_LEVMARQFunc& fjac,
		double * p, double * x, int m, int n, double * lb, double * ub,
		int itmax, double *opts, double *info, double *work, double *covar, void *adata);

extern int dlevmar_bc_dif(
		const ARM_LEVMARQFunc& func,
		double * p, double * x, int m, int n, double * lb, double * ub,
		int itmax, double *opts, double *info, double *work, double *covar, void *adata);

#ifdef HAVE_LAPACK
/* linear equation constrained minimization */
//extern int dlevmar_lec_der(
//		const CLMFunc& func,
//		const CLMFunc& fjac,
//		double *p, double *x, int m, int n, double *A, double *b, int k,
//		int itmax, double *opts, double *info, double *work, double *covar, void *adata);
//
//extern int dlevmar_lec_dif(
//		const CLMFunc& func,
//		double *p, double *x, int m, int n, double *A, double *b, int k,
//		int itmax, double *opts, double *info, double *work, double *covar, void *adata);
#endif /* HAVE_LAPACK */

/* linear system solvers */
#ifdef HAVE_LAPACK
extern int dAx_eq_b_QR(double *A, double *B, double *x, int m);
extern int dAx_eq_b_QRLS(double *A, double *B, double *x, int m, int n);
extern int dAx_eq_b_Chol(double *A, double *B, double *x, int m);
extern int dAx_eq_b_LU(double *A, double *B, double *x, int m);
extern int dAx_eq_b_SVD(double *A, double *B, double *x, int m);

#else // no LAPACK
extern int dAx_eq_b_LU_noLapack(double *A, double *B, double *x, int n);

#endif /* HAVE_LAPACK */

/* jacobian verification, double & single precision */
extern void dlevmar_chkjac(
		const ARM_LEVMARQFunc& func,
		const ARM_LEVMARQFunc& fjac,
		ARM_GP_Vector & p, int m, int n, void *adata, double *err);

#ifdef __cplusplus
}
#endif

#define LCAT_(a, b)    #a b
#define LCAT(a, b)    LCAT_(a, b) // force substitution
#define RCAT_(a, b)    a #b
#define RCAT(a, b)    RCAT_(a, b) // force substitution

#define __BLOCKSZ__       32 /* block size for cache-friendly matrix-matrix multiply. It should be
                              * such that __BLOCKSZ__^2*sizeof(LM_REAL) is smaller than the CPU (L1)
                              * data cache size. Notice that a value of 32 when LM_REAL=double assumes
                              * an 8Kb L1 data cache (32*32*8=8K). This is a concervative choice since
                              * newer Pentium 4s have a L1 data cache of size 16K, capable of holding
                              * up to 45x45 double blocks.
                              */
#define __BLOCKSZ__SQ    (__BLOCKSZ__)*(__BLOCKSZ__)

#ifdef _MSC_VER
#define inline __inline //MSVC
#elif !defined(__GNUC__)
#define inline //other than MSVC, GCC: define empty
#endif

/* add a prefix in front of a token */
#define LM_CAT__(a, b) a ## b
#define LM_CAT_(a, b) LM_CAT__(a, b) // force substitution
#define LM_ADD_PREFIX(s) LM_CAT_(LM_PREFIX, s)

#ifdef __cplusplus
extern "C" {
#endif

extern void dtrans_mat_mat_mult(double *a, double *b, int n, int m, int bsize);

extern void dfdif_forw_jac_approx(const ARM_LEVMARQFunc& func,
					double * p, double * hx, double * hxx, double delta,
					double * jac, int m, int n, void *adata);

extern void dfdif_cent_jac_approx(const ARM_LEVMARQFunc& func,
          double * p, double * hxm, double * hxp, double delta,
          double * jac, int m, int n, void *adata);

/* covariance of LS fit */
extern int dlevmar_covar(double *JtJ, double *C, double sumsq, int m, int n);

#ifdef __cplusplus
}
#endif

struct LMBC_DIF_DATA{
  const ARM_LEVMARQFunc * func;
  double * hx, * hxx;
  void *adata;
  double delta;
};

class ARM_LEVMARQDIFFunc :
	public ARM_LEVMARQFunc
{
public:
	void operator()(double p[], double hx[], int m, int n, void * adata) const
	{
		struct LMBC_DIF_DATA *dta = (struct LMBC_DIF_DATA *)adata;

		/* call user-supplied function passing it the user-supplied data */
		(*dta->func)(p, hx, m, n, dta->adata);
	}
};


class ARM_LEVMARQDIFJACFunc :
	public ARM_LEVMARQFunc
{
public:
	void operator()(double p[], double jac[], int m, int n, void * adata) const
	{
		struct LMBC_DIF_DATA *dta = (struct LMBC_DIF_DATA *)adata;

		(*dta->func)(p, dta->hx, m, n, dta->adata);

		dfdif_forw_jac_approx(*dta->func, p, dta->hx, dta->hxx, dta->delta, jac, m, n, dta->adata);
	}
};

inline int LEVMARQMin(ARM_LEVMARQFunc& objectiveFunction,
					 ARM_LEVMARQFunc * gradientFunction,
					 ARM_GP_Vector& firstGuess,
					 ARM_GP_Vector& functionAtfirstGuess,
					 ARM_GP_Vector * lowerBound,
					 ARM_GP_Vector * upperBound,
					 double * info = NULL,
					 int MaxIter = 250,
					 double * opts = NULL,
					 double * work = NULL,
					 double * covar = NULL,
					 void * adata = NULL)
{
	int xsize = firstGuess.size();
	int hxsize = functionAtfirstGuess.size();

	double * x = new double [xsize];
	double * hx = new double [hxsize];
    int k;
	for(k = 0; k < xsize; k++) x[k] = firstGuess[k];
	for(k = 0; k < hxsize; k++) hx[k] = functionAtfirstGuess[k];

	int status;
	double * popts = opts;
	if(popts == NULL)
	{
		popts = new double [LM_OPTS_SZ];
		popts[0] = LM_INIT_MU;
		popts[1] = 1E-10;
		popts[2] = 1E-10;
		popts[3] = 1E-12;
		popts[4] = LM_DIFF_DELTA;
	}

	if(gradientFunction == NULL)
	{
		if(lowerBound == NULL && upperBound == NULL)
		{
			try
			{
				status = dlevmar_dif(objectiveFunction, x, hx, xsize, hxsize, MaxIter, popts, info, work, covar, adata);
			}
			catch(Exception& e)
			{
				if(opts == NULL) delete [] popts;
				delete [] x;
				delete [] hx;

				throw(e);
			}
		}
		else
		{
			double * lbound = new double [xsize];
			double * ubound = new double [xsize];

			int lowerBoundSize = lowerBound == NULL ? 0 : lowerBound->size();

			for(k = 0; k < MIN(lowerBoundSize, xsize); k++) lbound[k] = (*lowerBound)[k];
			for(k = lowerBoundSize; k < xsize; k++) lbound[k] = - 1.e10;

			int upperBoundSize = upperBound == NULL ? 0 : upperBound->size();

			for(k = 0; k < MIN(upperBoundSize, xsize); k++) ubound[k] = (*upperBound)[k];
			for(k = upperBoundSize; k < xsize; k++) ubound[k] = 1.e10;
			
			try
			{
				status = dlevmar_bc_dif(objectiveFunction, x, hx, xsize, hxsize, lbound, ubound, MaxIter, popts, info, work, covar, adata);
			}
			catch(Exception& e)
			{
				if(opts == NULL) delete [] popts;
				delete [] lbound;
				delete [] ubound;
				delete [] x;
				delete [] hx;

				throw(e);
			}

			delete [] lbound;
			delete [] ubound;
		}
	}
	else
	{
		if(lowerBound == NULL && upperBound == NULL)
		{
			try
			{
				status = dlevmar_der(objectiveFunction, *gradientFunction, x, hx, xsize, hxsize, MaxIter, popts, info, work, covar, adata);
			}
			catch(Exception& e)
			{
				if(opts == NULL) delete [] popts;
				delete [] x;
				delete [] hx;

				throw(e);
			}
		}
		else
		{
			double * lbound = new double [xsize];
			double * ubound = new double [xsize];

			int lowerBoundSize = lowerBound == NULL ? 0 : lowerBound->size();

			for(k = 0; k < MIN(lowerBoundSize, xsize); k++) lbound[k] = (*lowerBound)[k];
			for(k = lowerBoundSize; k < xsize; k++) lbound[k] = - 1.e10;

			int upperBoundSize = upperBound == NULL ? 0 : upperBound->size();

			for(k = 0; k < MIN(upperBoundSize, xsize); k++) ubound[k] = (*upperBound)[k];
			for(k = upperBoundSize; k < xsize; k++) ubound[k] = 1.e10;
			
			try
			{
				status = dlevmar_bc_der(objectiveFunction, *gradientFunction, x, hx, xsize, hxsize, lbound, ubound, MaxIter, popts, info, work, covar, adata);
			}
			catch(Exception& e)
			{
				if(opts == NULL) delete [] popts;
				delete [] x;
				delete [] hx;
				delete [] lbound;
				delete [] ubound;

				throw(e);
			}

			delete [] lbound;
			delete [] ubound;
		}
	}

	for(k = 0; k < xsize; k++) firstGuess[k] = x[k];
	for(k = 0; k < hxsize; k++) functionAtfirstGuess[k] = hx[k];

	delete [] x;
	delete [] hx;
	if(opts == NULL) delete [] popts;

	return status;
}

inline int LEVMARQMinization_WithDerivatives(ARM_LEVMARQFunc& objectiveFunction,
											 ARM_LEVMARQFunc& gradientFunction,
											 ARM_GP_Vector& firstGuess,
											 ARM_GP_Vector& functionAtfirstGuess,
											 double * info = NULL,
											 int MaxIter = 250,
											 double * opts = NULL,
											 double * work = NULL,
											 double * covar = NULL,
											 void * adata = NULL
											 )
{
	return LEVMARQMin(objectiveFunction, &gradientFunction, firstGuess, functionAtfirstGuess, 
				NULL, NULL, info, MaxIter, opts, work, covar, adata);
}

inline int LEVMARQMinization_WithNumDerivatives(ARM_LEVMARQFunc& objectiveFunction,
											 ARM_GP_Vector& firstGuess,
											 ARM_GP_Vector& functionAtfirstGuess,
											 double * info = NULL,
											 int MaxIter = 250,
											 double * opts = NULL,
											 double * work = NULL,
											 double * covar = NULL,
											 void * adata = NULL
											 )
{
	return LEVMARQMin(objectiveFunction, NULL, firstGuess, functionAtfirstGuess, 
				NULL, NULL, info, MaxIter, opts, work, covar, adata);
}

inline int LEVMARQConstrainedMinization_WithDerivatives(ARM_LEVMARQFunc& objectiveFunction,
											 ARM_LEVMARQFunc& gradientFunction,
											 ARM_GP_Vector& firstGuess,
											 ARM_GP_Vector& functionAtfirstGuess,
											 ARM_GP_Vector& lowerBound,
											 ARM_GP_Vector& upperBound,
											 double * info = NULL,
											 int MaxIter = 250,
											 double * opts = NULL,
											 double * work = NULL,
											 double * covar = NULL,
											 void * adata = NULL
											 )
{
	return LEVMARQMin(objectiveFunction, &gradientFunction, firstGuess, functionAtfirstGuess, 
				&lowerBound, &upperBound, info, MaxIter, opts, work, covar, adata);
}

inline int LEVMARQConstrainedMinization_WithNumDerivatives(ARM_LEVMARQFunc& objectiveFunction,
											 ARM_GP_Vector& firstGuess,
											 ARM_GP_Vector& functionAtfirstGuess,
											 ARM_GP_Vector& lowerBound,
											 ARM_GP_Vector& upperBound,
											 double * info = NULL,
											 int MaxIter = 250,
											 double * opts = NULL,
											 double * work = NULL,
											 double * covar = NULL,
											 void * adata = NULL
											 )
{
	return LEVMARQMin(objectiveFunction, NULL, firstGuess, functionAtfirstGuess, 
				&lowerBound, &upperBound, info, MaxIter, opts, work, covar, adata);
}

CC_END_NAMESPACE()

#endif
