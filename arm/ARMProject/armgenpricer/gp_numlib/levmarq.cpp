
#include "gpnumlib/levmarq.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "armglob.h"

CC_BEGIN_NAMESPACE( ARM )

/////////////////////////////////////////////////////////////////////////////////
// 
//  Levenberg - Marquardt non-linear minimization algorithm
//  Copyright (C) 2004  Manolis Lourakis (lourakis@ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
/////////////////////////////////////////////////////////////////////////////////

#define EPSILON       1E-12
#define ONE_THIRD     0.3333333334 /* 1.0/3.0 */
#define SUBCNST(x) x##F
#define CNST(x) SUBCNST(x) // force substitution

/* 
 * This function seeks the parameter vector p that best describes the measurements vector x.
 * More precisely, given a vector function  func : R^m --> R^n with n>=m,
 * it finds p s.t. func(p) ~= x, i.e. the squared second order (i.e. L2) norm of
 * e=x-func(p) is minimized.
 *
 * This function requires an analytic jacobian. In case the latter is unavailable,
 * use LEVMAR_DIF() bellow
 *
 * Returns the number of iterations (>=0) if successfull, -1 if failed
 *
 * For more details, see H.B. Nielsen's (http://www.imm.dtu.dk/~hbn) IMM/DTU
 * tutorial at http://www.imm.dtu.dk/courses/02611/nllsq.pdf
 */

int dlevmar_der(
	const ARM_LEVMARQFunc& func,		// fonction du type (void)(double * p, double * hx, int m, int n, void * adata) de Rm dans Rn
	const ARM_LEVMARQFunc& fjac,		// jacobien 
	double * p,					// first guess, est modifié avec la solution
	double * x,
	int m,						// dimension de p
	int n,						// dimension de x
	int itmax,
	double opts[4],				// I: minim. options [\mu, \epsilon1, \epsilon2, \epsilon3]. 
								// Respectively the scale factor for initial \mu,
								// stopping thresholds for ||J^T e||_inf, ||Dp||_2 and ||e||_2. 
								// Set to NULL for defaults to be used
	double info[LM_INFO_SZ],	//	O: information regarding the minimization. Set to NULL if don't care
								//	info[0]= ||e||_2 at initial p.
								//	info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||Dp||_2, mu/max[J^T J]_ii ], all computed at estimated p.
								//	info[5]= # iterations,
								//	info[6]=reason for terminating: 1 - stopped by small gradient J^T e
								//									2 - stopped by small Dp
								//									3 - stopped by itmax
								//									4 - singular matrix. Restart from current p with increased mu 
								//									5 - no further error reduction is possible. Restart with increased mu
								//									6 - stopped by small ||e||_2
								//	info[7]= # function evaluations
								//	info[8]= # jacobian evaluations
	double * work,
	double * covar,
	void * adata)
{
int i, j, k, l, kk;
int worksz, freework=0, issolved;
/* temp work arrays */
double *e,          /* nx1 */
       *hx,         /* \hat{x}_i, nx1 */
       *jacTe,      /* J^T e_i mx1 */
       *jac,        /* nxm */
       *jacTjac,    /* mxm */
       *Dp,         /* mx1 */
   *diag_jacTjac,   /* diagonal of J^T J, mx1 */
       *pDp;        /* p + Dp, mx1 */

double mu,  /* damping constant */
                tmp; /* mainly used in matrix & vector multiplications */
double p_eL2, jacTe_inf, pDp_eL2; /* ||e(p)||_2, ||J^T e||_inf, ||e(p+Dp)||_2 */
double p_L2, Dp_L2=DBL_MAX, dF, dL;
double tau, eps1, eps2, eps2_sq, eps3;
double init_p_eL2;
int nu=2, nu2, stop, nfev, njev=0;
const int nm=n*m;

  mu=jacTe_inf=0.0; /* -Wall */

  if(n<m){
 // fprintf(stderr, LCAT(LEVMAR_DER, "(): cannot solve a problem with fewer measurements [%d] than unknowns [%d]\n"), n, m);
 //	exit(1);
	  ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"LEVMAR_DER() : cannot solve a problem with fewer measurements than unknown");
  }

/*  if(!fjac){
    fprintf(stderr, RCAT("No function specified for computing the jacobian in ", LEVMAR_DER)
        RCAT("().\nIf no such function is available, use ", LEVMAR_DIF) RCAT("() rather than ", LEVMAR_DER) "()\n");
    exit(1);
  }
*/

  if(opts){
	  tau=opts[0];
	  eps1=opts[1];
	  eps2=opts[2];
	  eps2_sq=opts[2]*opts[2];
    eps3=opts[3];
  }
  else{ // use default values
	  tau=CNST(LM_INIT_MU);
	  eps1=CNST(LM_STOP_THRESH);
	  eps2=CNST(LM_STOP_THRESH);
	  eps2_sq=CNST(LM_STOP_THRESH)*CNST(LM_STOP_THRESH);
    eps3=CNST(LM_STOP_THRESH);
  }

  if(!work){
    worksz=LM_DER_WORKSZ(m, n); //2*n+4*m + n*m + m*m;
    work=(double *)malloc(worksz*sizeof(double)); /* allocate a big chunk in one step */
	
    if(!work){
//      fprintf(stderr, LCAT(LEVMAR_DER, "(): memory allocation request failed\n"));
//      exit(1);
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"LEVMAR_DER() : memory allocation request failed");
    }
    freework=1;
  }

  /* set up work arrays */
  for (kk=0;kk<worksz;kk++) {
      work[kk] = 0.0;
  }
  e=work;
  hx=e + n;
  jacTe=hx + n;
  jac=jacTe + m;
  jacTjac=jac + nm;
  Dp=jacTjac + m*m;
  diag_jacTjac=Dp + m;
  pDp=diag_jacTjac + m;

  /* compute e=x - f(p) and its L2 norm */
  func(p, hx, m, n, adata); nfev=1;
  for(i=0, p_eL2=0.0; i<n; ++i){
    e[i]=tmp=x[i]-hx[i];
    p_eL2+=tmp*tmp;
  }
  init_p_eL2=p_eL2;

  for(k=stop=0; k<itmax && !stop; ++k){
    /* Note that p and e have been updated at a previous iteration */

    if(p_eL2<=eps3){ /* error is small */
      stop=6;
      break;
    }

    /* Compute the jacobian J at p,  J^T J,  J^T e,  ||J^T e||_inf and ||p||^2.
     * Since J^T J is symmetric, its computation can be speeded up by computing
     * only its upper triangular part and copying it to the lower part
     */

    fjac(p, jac, m, n, adata); ++njev;

    /* J^T J, J^T e */
    if(nm<__BLOCKSZ__SQ){ // this is a small problem
      /* This is the straightforward way to compute J^T J, J^T e. However, due to
       * its noncontinuous memory access pattern, it incures many cache misses when
       * applied to large minimization problems (i.e. problems involving a large
       * number of free variables and measurements), in which J is too large to
       * fit in the L1 cache. For such problems, a cache-efficient blocking scheme
       * is preferable.
       *
       * Thanks to John Nitao of Lawrence Livermore Lab for pointing out this
       * performance problem.
       *
       * On the other hand, the straightforward algorithm is faster on small
       * problems since in this case it avoids the overheads of blocking. 
       */

      for(i=0; i<m; ++i){
        for(j=i; j<m; ++j){
          int lm;

          for(l=0, tmp=0.0; l<n; ++l){
            lm=l*m;
            tmp+=jac[lm+i]*jac[lm+j];
          }

		      /* store tmp in the corresponding upper and lower part elements */
          jacTjac[i*m+j]=jacTjac[j*m+i]=tmp;
        }

        /* J^T e */
        for(l=0, tmp=0.0; l<n; ++l)
          tmp+=jac[l*m+i]*e[l];
        jacTe[i]=tmp;
      }
    }
    else{ // this is a large problem
      /* Cache efficient computation of J^T J based on blocking
       */
      dtrans_mat_mat_mult(jac, jacTjac, n, m, __BLOCKSZ__);

      /* cache efficient computation of J^T e */
      for(i=0; i<m; ++i)
        jacTe[i]=0.0;

      for(i=0; i<n; ++i){
        register double *jacrow;

        for(l=0, jacrow=jac+i*m, tmp=e[i]; l<m; ++l)
          jacTe[l]+=jacrow[l]*tmp;
      }
    }

	  /* Compute ||J^T e||_inf and ||p||^2 */
    for(i=0, p_L2=jacTe_inf=0.0; i<m; ++i){
      if(jacTe_inf < (tmp=FABS(jacTe[i]))) jacTe_inf=tmp;

      diag_jacTjac[i]=jacTjac[i*m+i]; /* save diagonal entries so that augmentation can be later canceled */
      p_L2+=p[i]*p[i];
    }
    //p_L2=sqrt(p_L2);

#if 0
if(!(k%100)){
  printf("Current estimate: ");
  for(i=0; i<m; ++i)
    printf("%.9g ", p[i]);
  printf("-- errors %.9g %0.9g\n", jacTe_inf, p_eL2);
}
#endif

    /* check for convergence */
    if((jacTe_inf <= eps1)){
      Dp_L2=0.0; /* no increment for p in this case */
      stop=1;
      break;
    }

   /* compute initial damping factor */
    if(k==0){
      for(i=0, tmp=DBL_MIN; i<m; ++i)
        if(diag_jacTjac[i]>tmp) tmp=diag_jacTjac[i]; /* find max diagonal element */
      mu=tau*tmp;
    }

    /* determine increment using adaptive damping */
    while(1){
      /* augment normal equations */
      for(i=0; i<m; ++i)
        jacTjac[i*m+i]+=mu;

      /* solve augmented equations */
#ifdef HAVE_LAPACK
      /* 5 alternatives are available: LU, Cholesky, 2 variants of QR decomposition and SVD.
       * Cholesky is the fastest but might be inaccurate; QR is slower but more accurate;
       * SVD is the slowest but most accurate; LU offers a tradeoff between accuracy and speed
       */

      issolved=dAx_eq_b_LU(jacTjac, jacTe, Dp, m);
      //issolved=AX_EQ_B_CHOL(jacTjac, jacTe, Dp, m);
      //issolved=AX_EQ_B_QR(jacTjac, jacTe, Dp, m);
      //issolved=AX_EQ_B_QRLS(jacTjac, jacTe, Dp, m, m);
      //issolved=AX_EQ_B_SVD(jacTjac, jacTe, Dp, m);

#else
      /* use the LU included with levmar */
      issolved=dAx_eq_b_LU_noLapack(jacTjac, jacTe, Dp, m);
#endif /* HAVE_LAPACK */

      if(issolved){
        /* compute p's new estimate and ||Dp||^2 */
        for(i=0, Dp_L2=0.0; i<m; ++i){
          pDp[i]=p[i] + (tmp=Dp[i]);
          Dp_L2+=tmp*tmp;
        }
        //Dp_L2=sqrt(Dp_L2);

        if(Dp_L2<=eps2_sq*p_L2){ /* relative change in p is small, stop */
        //if(Dp_L2<=eps2*(p_L2 + eps2)){ /* relative change in p is small, stop */
          stop=2;
          break;
        }

       if(Dp_L2>=(p_L2+eps2)/(CNST(EPSILON)*CNST(EPSILON))){ /* almost singular */
       //if(Dp_L2>=(p_L2+eps2)/CNST(EPSILON)){ /* almost singular */
         stop=4;
         break;
       }

        func(pDp, hx, m, n, adata); ++nfev; /* evaluate function at p + Dp */
        for(i=0, pDp_eL2=0.0; i<n; ++i){ /* compute ||e(pDp)||_2 */
          hx[i]=tmp=x[i]-hx[i];
          pDp_eL2+=tmp*tmp;
        }

        for(i=0, dL=0.0; i<m; ++i)
          dL+=Dp[i]*(mu*Dp[i]+jacTe[i]);

        dF=p_eL2-pDp_eL2;

        if(dL>0.0 && dF>0.0){ /* reduction in error, increment is accepted */
          tmp=(CNST(2.0)*dF/dL-CNST(1.0));
          tmp=CNST(1.0)-tmp*tmp*tmp;
          mu=mu*( (tmp>=CNST(ONE_THIRD))? tmp : CNST(ONE_THIRD) );
          nu=2;

          for(i=0 ; i<m; ++i) /* update p's estimate */
            p[i]=pDp[i];

          for(i=0; i<n; ++i) /* update e and ||e||_2 */
            e[i]=hx[i];
          p_eL2=pDp_eL2;
          break;
        }
      }

      /* if this point is reached, either the linear system could not be solved or
       * the error did not reduce; in any case, the increment must be rejected
       */

      mu*=nu;
      nu2=nu<<1; // 2*nu;
      if(nu2<=nu){ /* nu has wrapped around (overflown). Thanks to Frank Jordan for spotting this case */
        stop=5;
        break;
      }
      nu=nu2;

      for(i=0; i<m; ++i) /* restore diagonal J^T J entries */
        jacTjac[i*m+i]=diag_jacTjac[i];
    } /* inner loop */
  }

  if(k>=itmax) stop=3;

  for(i=0; i<m; ++i) /* restore diagonal J^T J entries */
    jacTjac[i*m+i]=diag_jacTjac[i];

  if(info){
    info[0]=init_p_eL2;
    info[1]=p_eL2;
    info[2]=jacTe_inf;
    info[3]=Dp_L2;
    for(i=0, tmp=DBL_MIN; i<m; ++i)
      if(tmp<jacTjac[i*m+i]) tmp=jacTjac[i*m+i];
    info[4]=mu/tmp;
    info[5]=(double)k;
    info[6]=(double)stop;
    info[7]=(double)nfev;
    info[8]=(double)njev;
  }

  /* covariance matrix */
  if(covar){
    dlevmar_covar(jacTjac, covar, p_eL2, m, n);
  }

  if(freework) free(work);

  return (stop!=4)?  k : -1;
}

int dlevmar_dif(
	const ARM_LEVMARQFunc& func,		// fonction du type (void)(double * p, double * hx, int m, int n, void * adata) de Rm dans Rn
	double * p,					// first guess, est modifié avec la solution
	double * x,
	int m,						// dimension de p
	int n,						// dimension de x
	int itmax,
	double opts[5],    /* I: opts[0-4] = minim. options [\mu, \epsilon1, \epsilon2, \epsilon3, \delta]. Respectively the
                       * scale factor for initial \mu, stopping thresholds for ||J^T e||_inf, ||Dp||_2 and ||e||_2 and
                       * the step used in difference approximation to the jacobian. Set to NULL for defaults to be used.
                       * If \delta<0, the jacobian is approximated with central differences which are more accurate
                       * (but slower!) compared to the forward differences employed by default. 
                       */
 	double info[LM_INFO_SZ],	//	O: information regarding the minimization. Set to NULL if don't care
								//	info[0]= ||e||_2 at initial p.
								//	info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||Dp||_2, mu/max[J^T J]_ii ], all computed at estimated p.
								//	info[5]= # iterations,
								//	info[6]=reason for terminating: 1 - stopped by small gradient J^T e
								//									2 - stopped by small Dp
								//									3 - stopped by itmax
								//									4 - singular matrix. Restart from current p with increased mu 
								//									5 - no further error reduction is possible. Restart with increased mu
								//									6 - stopped by small ||e||_2
								//	info[7]= # function evaluations
								//	info[8]= # jacobian evaluations
	double * work,
	double * covar,
	void * adata)
{
	int i, j, k, l;
	int worksz, freework=0, issolved;

	// temp work arrays
	double * e, * hx, * jacTe, * jac, * jacTjac, * Dp, * diag_jacTjac, * pDp, * wrk;

	int using_ffdif = 1;
	double * wrk2 = NULL; // nx1, used for differentiating with central differences only

	double mu, tmp; 
	double p_eL2, jacTe_inf, pDp_eL2; // ||e(p)||_2, ||J^T e||_inf, ||e(p+Dp)||_2 
	double p_L2, Dp_L2 = DBL_MAX, dF, dL;
	double tau, eps1, eps2, eps2_sq, eps3, delta;
	double init_p_eL2;
	int nu, nu2, stop, nfev, njap = 0, K = (m >= 10)? m : 10, updjac, updp = 1, newjac;
	const int nm = n * m;

	mu = jacTe_inf = p_L2 = 0.0; /* -Wall */
	stop = updjac = newjac = 0; /* -Wall */

	if(n < m)
	{
		return -1; // plus d'inconnues que d'équations
	}

	if(opts)
	{
		tau		= opts[0];
		eps1	= opts[1];
		eps2	= opts[2];
		eps2_sq	= opts[2] * opts[2];
		eps3	= opts[3];
		delta	= opts[4];

		if(delta<0.0)
		{
			delta = - delta; /* make positive */
			using_ffdif = 0; /* use central differencing */
			wrk2 = (double *)malloc(n*sizeof(double));
			// assert(wrk2 != 0);
		}
	}
	else
	{ // use default values
		tau		= (LM_INIT_MU);
		eps1	= (LM_STOP_THRESH);
		eps2	= (LM_STOP_THRESH);
		eps2_sq	= (LM_STOP_THRESH) * (LM_STOP_THRESH);
		eps3	= (LM_STOP_THRESH);
		delta	= (LM_DIFF_DELTA);
	}

	if(!work)
	{
		worksz	= LM_DIF_WORKSZ(m, n); //3*n+4*m + n*m + m*m;
		work	= (double*)malloc(worksz*sizeof(double)); /* allocate a big chunk in one step */
		//assert(work != 0);
		freework= 1;
	}

	/* set up work arrays */
	e				= work;
	hx				= e + n;
	jacTe			= hx + n;
	jac				= jacTe + m;
	jacTjac			= jac + nm;
	Dp				= jacTjac + m * m;
	diag_jacTjac	= Dp + m;
	pDp				= diag_jacTjac + m;
	wrk				= pDp + m;

	/* compute e=x - f(p) and its L2 norm */
	func(p, hx, m, n, adata);
	nfev = 1;

	for(i=0, p_eL2 = 0.0; i < n; ++i)
	{
		e[i] = tmp = x[i] - hx[i];
		p_eL2 += tmp * tmp;
	}

	init_p_eL2 = p_eL2;

	nu = 20; /* force computation of J */

	for(k = 0; k < itmax; ++ k)
	{
		/* Note that p and e have been updated at a previous iteration */

		if(p_eL2 <= eps3)
		{ /* error is small */
			stop=6;
			break;
		}

		/* Compute the jacobian J at p,  J^T J,  J^T e,  ||J^T e||_inf and ||p||^2.
		* The symmetry of J^T J is again exploited for speed
		*/

		if((updp && nu > 16) || updjac == K)
		{ /* compute difference approximation to J */
			if(using_ffdif)
			{ /* use forward differences */
				dfdif_forw_jac_approx(func, p, hx, wrk, delta, jac, m, n, adata);
				++ njap; 
				nfev += m;
			}
			else
			{ /* use central differences */
				dfdif_cent_jac_approx(func, p, wrk, wrk2, delta, jac, m, n, adata);
				++ njap; 
				nfev += 2 * m;
			}
		
			nu		= 2; 
			updjac	= 0; 
			updp	= 0; 
			newjac	= 1;
		}

		if(newjac)
		{ /* jacobian has changed, recompute J^T J, J^t e, etc */

			newjac = 0;

			/* J^T J, J^T e */
			if(nm <=__BLOCKSZ__SQ)
			{ // this is a small problem
			/* This is the straightforward way to compute J^T J, J^T e. However, due to
				* its noncontinuous memory access pattern, it incures many cache misses when
				* applied to large minimization problems (i.e. problems involving a large
				* number of free variables and measurements), in which J is too large to
				* fit in the L1 cache. For such problems, a cache-efficient blocking scheme
				* is preferable.
				*
				* Thanks to John Nitao of Lawrence Livermore Lab for pointing out this
				* performance problem.
				*
				* On the other hand, the straightforward algorithm is faster on small
				* problems since in this case it avoids the overheads of blocking. 
				*/

				for(i = 0; i < m; ++i)
				{
					for(j = i; j < m; ++j)
					{
						int lm;

						for(l = 0, tmp = 0.0; l < n; ++l)
						{
							lm	= l * m;
							tmp	+= jac[lm+i] * jac[lm+j];
						}

						jacTjac[i*m+j] = jacTjac[j*m+i] = tmp;
					}

					/* J^T e */
					for(l=0, tmp = 0.0; l < n; ++l)
					{
						tmp += jac[l*m+i] * e[l];
					}

					jacTe[i] = tmp;
				}
			}
			else
			{
				// this is a large problem
				// Cache efficient computation of J^T J based on blocking
				
				dtrans_mat_mat_mult(jac, jacTjac, n, m, __BLOCKSZ__);

				// cache efficient computation of J^T e 
				for(i = 0; i < m; ++i) jacTe[i] = 0.0;

				for(i = 0; i < n; ++i)
				{
					double * jacrow;

					for(l = 0, jacrow = jac+i*m, tmp = e[i]; l < m; ++l) jacTe[l] += jacrow[l] * tmp;
				}
			}

			// Compute ||J^T e||_inf and ||p||^2 
			for(i = 0, p_L2 = jacTe_inf = 0.0; i < m; ++i)
			{
				if(jacTe_inf < (tmp = FABS(jacTe[i]))) jacTe_inf = tmp;

				diag_jacTjac[i] = jacTjac[i*m+i]; // save diagonal entries so that augmentation can be later canceled
				p_L2 += p[i] * p[i];
			}
			//p_L2=sqrt(p_L2);
		}

		/* check for convergence */
		if((jacTe_inf <= eps1))
		{
			Dp_L2 = 0.0; /* no increment for p in this case */
			stop = 1;
			break;
		}

		/* compute initial damping factor */
		if(k == 0)
		{
			for(i = 0, tmp = DBL_MIN; i < m; ++i)
			{
				if(diag_jacTjac[i] > tmp) tmp = diag_jacTjac[i]; /* find max diagonal element */
			}

			mu = tau * tmp;
		}

		/* determine increment using adaptive damping */

		/* augment normal equations */
		for(i = 0; i < m; ++i) jacTjac[i*m+i] += mu;

		/* solve augmented equations */
		#ifdef HAVE_LAPACK
		/* 5 alternatives are available: LU, Cholesky, 2 variants of QR decomposition and SVD.
			* Cholesky is the fastest but might be inaccurate; QR is slower but more accurate;
			* SVD is the slowest but most accurate; LU offers a tradeoff between accuracy and speed
			*/

		issolved = dAx_eq_b_LU(jacTjac, jacTe, Dp, m);
		//issolved=AX_EQ_B_CHOL(jacTjac, jacTe, Dp, m);
		//issolved=AX_EQ_B_QR(jacTjac, jacTe, Dp, m);
		//issolved=AX_EQ_B_QRLS(jacTjac, jacTe, Dp, m, m);
		//issolved=AX_EQ_B_SVD(jacTjac, jacTe, Dp, m);
		#else
		/* use the LU included with levmar */
		issolved = dAx_eq_b_LU_noLapack(jacTjac, jacTe, Dp, m);
		#endif /* HAVE_LAPACK */

		if(issolved)
		{
		/* compute p's new estimate and ||Dp||^2 */
			for(i = 0, Dp_L2 = 0.0; i < m; ++i)
			{
				pDp[i] = p[i] + (tmp = Dp[i]);
				Dp_L2 += tmp * tmp;
			}
			//Dp_L2=sqrt(Dp_L2);

			if(Dp_L2 <= eps2_sq * p_L2)
			{ /* relative change in p is small, stop */
			//if(Dp_L2<=eps2*(p_L2 + eps2)){ /* relative change in p is small, stop */
				stop=2;
				break;
			}

			if(Dp_L2 >= (p_L2 + eps2) / ((EPSILON) * (EPSILON)))
			{ /* almost singular */
			//if(Dp_L2>=(p_L2+eps2)/CNST(EPSILON)){ /* almost singular */
				stop=4;
				break;
			}

			func(pDp, wrk, m, n, adata); /* evaluate function at p + Dp */
			++ nfev; 

			for(i = 0, pDp_eL2 = 0.0; i < n; ++i)
			{ /* compute ||e(pDp)||_2 */
				tmp = x[i] - wrk[i];
				pDp_eL2 += tmp * tmp;
			}

			dF = p_eL2 - pDp_eL2;
			if(updp || dF > 0)
			{ /* update jac */
				for(i = 0; i < n; ++i)
				{
					for(l = 0, tmp = 0.0; l < m; ++l) tmp += jac[i*m+l] * Dp[l]; /* (J * Dp)[i] */
					
					tmp = (wrk[i] - hx[i] - tmp) / Dp_L2; /* (f(p+dp)[i] - f(p)[i] - (J * Dp)[i])/(dp^T*dp) */
					
					for(j = 0; j < m; ++j) jac[i*m+j] += tmp * Dp[j];
				}

				++ updjac;
				newjac = 1;
			}

			for(i = 0, dL = 0.0; i < m; ++i) dL += Dp[i] * (mu * Dp[i] + jacTe[i]);

			if(dL > 0.0 && dF > 0.0)
			{ /* reduction in error, increment is accepted */
				dF	= ((2.0)*dF/dL - (1.0));
				tmp	= dF * dF * dF;
				tmp	= (1.0) - tmp * tmp * dF;
				mu	= mu *( (tmp >= (ONE_THIRD))? tmp : (ONE_THIRD) );
				nu	= 2;

				/* update p's estimate */
				for(i = 0 ; i < m; ++i) p[i] = pDp[i];

				/* update e, hx and ||e||_2 */
				for(i=0; i<n; ++i)
				{ 
					e[i] = x[i] - wrk[i];
					hx[i]= wrk[i];
				}
				p_eL2 = pDp_eL2;
				updp = 1;
				continue;
			}
		}

		/* if this point is reached, either the linear system could not be solved or
			* the error did not reduce; in any case, the increment must be rejected
			*/

		mu	*= nu;
		nu2	= nu<<1; // 2*nu;

		if(nu2 <= nu)
		{ /* nu has wrapped around (overflown). Thanks to Frank Jordan for spotting this case */
			stop = 5;
			break;
		}

		nu = nu2;

		/* restore diagonal J^T J entries */
		for(i = 0; i < m; ++i) jacTjac[i*m+i] = diag_jacTjac[i];
	}  

	if(k >= itmax) stop=3;

	/* restore diagonal J^T J entries */
	for(i = 0; i< m; ++i) jacTjac[i*m+i] = diag_jacTjac[i];

	if(info)
	{
		info[0]	= init_p_eL2;
		info[1]	= p_eL2;
		info[2]	= jacTe_inf;
		info[3]	= Dp_L2;
		for(i = 0, tmp = DBL_MIN; i < m; ++i)
		{
			if(tmp < jacTjac[i*m+i]) tmp = jacTjac[i*m+i];
		}
		info[4]	= mu / tmp;
		info[5]	=(double)k;
		info[6]	=(double)stop;
		info[7]	=(double)nfev;
		info[8]	=(double)njap;
	}

	/* covariance matrix */
	if(covar)
	{
		dlevmar_covar(jacTjac, covar, p_eL2, m, n);
	}

	                                                            
	if(freework) free(work);

	if(wrk2) free(wrk2);

	return (stop!=4)?  k : -1;
}

/////////////////////////////////////////////////////////////////////////////////
// 
//  Solution of linear systems involved in the Levenberg - Marquardt
//  minimization algorithm
//  Copyright (C) 2004  Manolis Lourakis (lourakis@ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
/////////////////////////////////////////////////////////////////////////////////

/******************************************************************************** 
 * LAPACK-based implementations for various linear system solvers. The same core
 * code is used with appropriate #defines to derive single and double precision
 * solver versions, see also Axb_core.c
 ********************************************************************************/

/////////////////////////////////////////////////////////////////////////////////
// 
//  Solution of linear systems involved in the Levenberg - Marquardt
//  minimization algorithm
//  Copyright (C) 2004  Manolis Lourakis (lourakis@ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
/////////////////////////////////////////////////////////////////////////////////


/* Solvers for the linear systems Ax=b. Solvers should NOT modify their A & B arguments! */

#ifdef LINSOLVERS_RETAIN_MEMORY
#define __STATIC__ static
#else
#define __STATIC__ // empty
#endif /* LINSOLVERS_RETAIN_MEMORY */

#ifdef HAVE_LAPACK

/* prototypes of LAPACK routines */

#define GEQRF dgeqrf_
#define ORGQR dorgqr_
#define TRTRS dtrtrs_
#define POTF2 dpotf2_
#define POTRF dpotrf_
#define GETRF dgetrf_
#define GETRS dgetrs_
#define GESVD dgesvd_
#define GESDD dgesdd_

/* QR decomposition */
extern int GEQRF(int *m, int *n, double *a, int *lda, double *tau, double *work, int *lwork, int *info);
extern int ORGQR(int *m, int *n, int *k, double *a, int *lda, double *tau, double *work, int *lwork, int *info);

/* solution of triangular systems */
extern int TRTRS(char *uplo, char *trans, char *diag, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, int *info);

/* cholesky decomposition */
extern int POTF2(char *uplo, int *n, double *a, int *lda, int *info);
extern int POTRF(char *uplo, int *n, double *a, int *lda, int *info); /* block version of dpotf2 */

/* LU decomposition and systems solution */
extern int GETRF(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
extern int GETRS(char *trans, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

/* Singular Value Decomposition (SVD) */
extern int GESVD(char *jobu, char *jobvt, int *m, int *n, double *a, int *lda, double *s, double *u, int *ldu,
                   double *vt, int *ldvt, double *work, int *lwork, int *info);

/* lapack 3.0 new SVD routine, faster than xgesvd().
 * In case that your version of LAPACK does not include them, use the above two older routines
 */
extern int GESDD(char *jobz, int *m, int *n, double *a, int *lda, double *s, double *u, int *ldu, double *vt, int *ldvt,
                   double *work, int *lwork, int *iwork, int *info);

/* precision-specific definitions */
#define AX_EQ_B_QR dAx_eq_b_QR
#define AX_EQ_B_QRLS dAx_eq_b_QRLS
#define AX_EQ_B_CHOL dAx_eq_b_Chol
#define AX_EQ_B_LU dAx_eq_b_LU
#define AX_EQ_B_SVD dAx_eq_b_SVD

/*
 * This function returns the solution of Ax = b
 *
 * The function is based on QR decomposition with explicit computation of Q:
 * If A=Q R with Q orthogonal and R upper triangular, the linear system becomes
 * Q R x = b or R x = Q^T b.
 * The last equation can be solved directly.
 *
 * A is mxm, b is mx1
 *
 * The function returns 0 in case of error, 1 if successfull
 *
 * This function is often called repetitively to solve problems of identical
 * dimensions. To avoid repetitive malloc's and free's, allocated memory is
 * retained between calls and free'd-malloc'ed when not of the appropriate size.
 * A call with NULL as the first argument forces this memory to be released.
 */
int AX_EQ_B_QR(double *A, double *B, double *x, int m)
{
__STATIC__ double *buf=NULL;
__STATIC__ int buf_sz=0;

double *a, *qtb, *tau, *r, *work;
int a_sz, qtb_sz, tau_sz, r_sz, tot_sz;
register int i, j;
int info, worksz, nrhs=1;
register double sum;

#ifdef LINSOLVERS_RETAIN_MEMORY
    if(!A){
      if(buf) free(buf);
      buf_sz=0;
      return 1;
    }
#endif /* LINSOLVERS_RETAIN_MEMORY */
   
    /* calculate required memory size */
    a_sz=m*m;
    qtb_sz=m;
    tau_sz=m;
    r_sz=m*m; /* only the upper triangular part really needed */
    worksz=3*m; /* this is probably too much */
    tot_sz=a_sz + qtb_sz + tau_sz + r_sz + worksz;

#ifdef LINSOLVERS_RETAIN_MEMORY
    if(tot_sz>buf_sz){ /* insufficient memory, allocate a "big" memory chunk at once */
      if(buf) free(buf); /* free previously allocated memory */

      buf_sz=tot_sz;
      buf=(double *)malloc(buf_sz*sizeof(double));
      if(!buf){
//        fprintf(stderr, RCAT("memory allocation in ", AX_EQ_B_QR) "() failed!\n");
//        exit(1);
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"memory allocation in AX_EQ_B_QR failed");
      }
    }
#else
      buf_sz=tot_sz;
      buf=(double *)malloc(buf_sz*sizeof(double));
      if(!buf){
//        fprintf(stderr, RCAT("memory allocation in ", AX_EQ_B_QR) "() failed!\n");
//        exit(1);
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"memory allocation in AX_EQ_B_QR failed");
      }
#endif /* LINSOLVERS_RETAIN_MEMORY */

    a=buf;
    qtb=a+a_sz;
    tau=qtb+qtb_sz;
    r=tau+tau_sz;
    work=r+r_sz;

  /* store A (column major!) into a */
	for(i=0; i<m; i++)
		for(j=0; j<m; j++)
			a[i+j*m]=A[i*m+j];

  /* QR decomposition of A */
  GEQRF((int *)&m, (int *)&m, a, (int *)&m, tau, work, (int *)&worksz, (int *)&info);
  /* error treatment */
  if(info!=0){
    if(info<0){
//      fprintf(stderr, RCAT(RCAT("LAPACK error: illegal value for argument %d of ", GEQRF) " in ", AX_EQ_B_QR) "()\n", -info);
//		exit(1);
	  ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"LAPACK error: illegal value for argument");
    }
    else{
      fprintf(stderr, RCAT(RCAT("Unknown LAPACK error %d for ", GEQRF) " in ", AX_EQ_B_QR) "()\n", info);
#ifndef LINSOLVERS_RETAIN_MEMORY
      free(buf);
#endif

      return 0;
    }
  }

  /* R is stored in the upper triangular part of a; copy it in r so that ORGQR() below won't destroy it */ 
  for(i=0; i<r_sz; i++)
    r[i]=a[i];

  /* compute Q using the elementary reflectors computed by the above decomposition */
  ORGQR((int *)&m, (int *)&m, (int *)&m, a, (int *)&m, tau, work, (int *)&worksz, (int *)&info);
  if(info!=0){
    if(info<0){
      //fprintf(stderr, RCAT(RCAT("LAPACK error: illegal value for argument %d of ", ORGQR) " in ", AX_EQ_B_QR) "()\n", -info);
      //exit(1);
	  ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"LAPACK error: illegal value for argument");
    }
    else{
      fprintf(stderr, RCAT("Unknown LAPACK error (%d) in ", AX_EQ_B_QR) "()\n", info);
#ifndef LINSOLVERS_RETAIN_MEMORY
      free(buf);
#endif

      return 0;
    }
  }

  /* Q is now in a; compute Q^T b in qtb */
  for(i=0; i<m; i++){
    for(j=0, sum=0.0; j<m; j++)
      sum+=a[i*m+j]*B[j];
    qtb[i]=sum;
  }

  /* solve the linear system R x = Q^t b */
  TRTRS("U", "N", "N", (int *)&m, (int *)&nrhs, r, (int *)&m, qtb, (int *)&m, &info);
  /* error treatment */
  if(info!=0){
    if(info<0){
//      fprintf(stderr, RCAT(RCAT("LAPACK error: illegal value for argument %d of ", TRTRS) " in ", AX_EQ_B_QR) "()\n", -info);
//      exit(1);
	  ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"LAPACK error : illegal value for argument");
    }
    else{
      fprintf(stderr, RCAT("LAPACK error: the %d-th diagonal element of A is zero (singular matrix) in ", AX_EQ_B_QR) "()\n", info);
#ifndef LINSOLVERS_RETAIN_MEMORY
      free(buf);
#endif

      return 0;
    }
  }

	/* copy the result in x */
	for(i=0; i<m; i++)
    x[i]=qtb[i];

#ifndef LINSOLVERS_RETAIN_MEMORY
  free(buf);
#endif

	return 1;
}

/*
 * This function returns the solution of min_x ||Ax - b||
 *
 * || . || is the second order (i.e. L2) norm. This is a least squares technique that
 * is based on QR decomposition:
 * If A=Q R with Q orthogonal and R upper triangular, the normal equations become
 * (A^T A) x = A^T b  or (R^T Q^T Q R) x = A^T b or (R^T R) x = A^T b.
 * This amounts to solving R^T y = A^T b for y and then R x = y for x
 * Note that Q does not need to be explicitly computed
 *
 * A is mxn, b is mx1
 *
 * The function returns 0 in case of error, 1 if successfull
 *
 * This function is often called repetitively to solve problems of identical
 * dimensions. To avoid repetitive malloc's and free's, allocated memory is
 * retained between calls and free'd-malloc'ed when not of the appropriate size.
 * A call with NULL as the first argument forces this memory to be released.
 */
int AX_EQ_B_QRLS(double *A, double *B, double *x, int m, int n)
{
__STATIC__ double *buf=NULL;
__STATIC__ int buf_sz=0;

double *a, *atb, *tau, *r, *work;
int a_sz, atb_sz, tau_sz, r_sz, tot_sz;
register int i, j;
int info, worksz, nrhs=1;
register double sum;
   
#ifdef LINSOLVERS_RETAIN_MEMORY
    if(!A){
      if(buf) free(buf);
      buf_sz=0;
      return 1;
    }
#endif /* LINSOLVERS_RETAIN_MEMORY */
   
    if(m<n){
//		  fprintf(stderr, RCAT("Normal equations require that the number of rows is greater than number of columns in ", AX_EQ_B_QRLS) "() [%d x %d]! -- try transposing\n", m, n);
//		  exit(1);
		  ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"Normal equations require that the number of rows is greater than number of columns");
	  }
      
    /* calculate required memory size */
    a_sz=m*n;
    atb_sz=n;
    tau_sz=n;
    r_sz=n*n;
    worksz=3*n; /* this is probably too much */
    tot_sz=a_sz + atb_sz + tau_sz + r_sz + worksz;

#ifdef LINSOLVERS_RETAIN_MEMORY
    if(tot_sz>buf_sz){ /* insufficient memory, allocate a "big" memory chunk at once */
      if(buf) free(buf); /* free previously allocated memory */

      buf_sz=tot_sz;
      buf=(double *)malloc(buf_sz*sizeof(double));
      if(!buf){
//        fprintf(stderr, RCAT("memory allocation in ", AX_EQ_B_QRLS) "() failed!\n");
//        exit(1);
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"memory allocation in AX_EQ_B_QRLS failed");
      }
    }
#else
      buf_sz=tot_sz;
      buf=(double *)malloc(buf_sz*sizeof(double));
      if(!buf){
//        fprintf(stderr, RCAT("memory allocation in ", AX_EQ_B_QRLS) "() failed!\n");
//        exit(1);
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"memory allocation in AX_EQ_B_QRLS failed");
      }
#endif /* LINSOLVERS_RETAIN_MEMORY */

    a=buf;
    atb=a+a_sz;
    tau=atb+atb_sz;
    r=tau+tau_sz;
    work=r+r_sz;

  /* store A (column major!) into a */
	for(i=0; i<m; i++)
		for(j=0; j<n; j++)
			a[i+j*m]=A[i*n+j];

  /* compute A^T b in atb */
  for(i=0; i<n; i++){
    for(j=0, sum=0.0; j<m; j++)
      sum+=A[j*n+i]*B[j];
    atb[i]=sum;
  }

  /* QR decomposition of A */
  GEQRF((int *)&m, (int *)&n, a, (int *)&m, tau, work, (int *)&worksz, (int *)&info);
  /* error treatment */
  if(info!=0){
    if(info<0){
      //fprintf(stderr, RCAT(RCAT("LAPACK error: illegal value for argument %d of ", GEQRF) " in ", AX_EQ_B_QRLS) "()\n", -info);
      //exit(1);
	  ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"LAPACK error: illegal value for argument");
    }
    else{
      fprintf(stderr, RCAT(RCAT("Unknown LAPACK error %d for ", GEQRF) " in ", AX_EQ_B_QRLS) "()\n", info);
#ifndef LINSOLVERS_RETAIN_MEMORY
      free(buf);
#endif

      return 0;
    }
  }

  /* R is stored in the upper triangular part of a. Note that a is mxn while r nxn */
  for(j=0; j<n; j++){
    for(i=0; i<=j; i++)
      r[i+j*n]=a[i+j*m];

    /* lower part is zero */
    for(i=j+1; i<n; i++)
      r[i+j*n]=0.0;
  }

  /* solve the linear system R^T y = A^t b */
  TRTRS("U", "T", "N", (int *)&n, (int *)&nrhs, r, (int *)&n, atb, (int *)&n, &info);
  /* error treatment */
  if(info!=0){
    if(info<0){
//      fprintf(stderr, RCAT(RCAT("LAPACK error: illegal value for argument %d of ", TRTRS) " in ", AX_EQ_B_QRLS) "()\n", -info);
//      exit(1);
	  ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"LAPACK error: illegal value for argument");
    }
    else{
      fprintf(stderr, RCAT("LAPACK error: the %d-th diagonal element of A is zero (singular matrix) in ", AX_EQ_B_QRLS) "()\n", info);
#ifndef LINSOLVERS_RETAIN_MEMORY
      free(buf);
#endif

      return 0;
    }
  }

  /* solve the linear system R x = y */
  TRTRS("U", "N", "N", (int *)&n, (int *)&nrhs, r, (int *)&n, atb, (int *)&n, &info);
  /* error treatment */
  if(info!=0){
    if(info<0){
//      fprintf(stderr, RCAT(RCAT("LAPACK error: illegal value for argument %d of ", TRTRS) " in ", AX_EQ_B_QRLS) "()\n", -info);
//      exit(1);
	  ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"LAPACK error: illegal value for argument");
    }
    else{
      fprintf(stderr, RCAT("LAPACK error: the %d-th diagonal element of A is zero (singular matrix) in ", AX_EQ_B_QRLS) "()\n", info);
#ifndef LINSOLVERS_RETAIN_MEMORY
      free(buf);
#endif

      return 0;
    }
  }

	/* copy the result in x */
	for(i=0; i<n; i++)
    x[i]=atb[i];

#ifndef LINSOLVERS_RETAIN_MEMORY
  free(buf);
#endif

	return 1;
}

/*
 * This function returns the solution of Ax=b
 *
 * The function assumes that A is symmetric & postive definite and employs
 * the Cholesky decomposition:
 * If A=U^T U with U upper triangular, the system to be solved becomes
 * (U^T U) x = b
 * This amount to solving U^T y = b for y and then U x = y for x
 *
 * A is mxm, b is mx1
 *
 * The function returns 0 in case of error, 1 if successfull
 *
 * This function is often called repetitively to solve problems of identical
 * dimensions. To avoid repetitive malloc's and free's, allocated memory is
 * retained between calls and free'd-malloc'ed when not of the appropriate size.
 * A call with NULL as the first argument forces this memory to be released.
 */
int AX_EQ_B_CHOL(double *A, double *B, double *x, int m)
{
__STATIC__ double *buf=NULL;
__STATIC__ int buf_sz=0;

double *a, *b;
int a_sz, b_sz, tot_sz;
register int i, j;
int info, nrhs=1;
   
#ifdef LINSOLVERS_RETAIN_MEMORY
    if(!A){
      if(buf) free(buf);
      buf_sz=0;
      return 1;
    }
#endif /* LINSOLVERS_RETAIN_MEMORY */
   
    /* calculate required memory size */
    a_sz=m*m;
    b_sz=m;
    tot_sz=a_sz + b_sz;

#ifdef LINSOLVERS_RETAIN_MEMORY
    if(tot_sz>buf_sz){ /* insufficient memory, allocate a "big" memory chunk at once */
      if(buf) free(buf); /* free previously allocated memory */

      buf_sz=tot_sz;
      buf=(double *)malloc(buf_sz*sizeof(double));
      if(!buf){
//        fprintf(stderr, RCAT("memory allocation in ", AX_EQ_B_CHOL) "() failed!\n");
//        exit(1);
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"memory allocation in AX_EQ_B_CHOL");
      }
    }
#else
      buf_sz=tot_sz;
      buf=(double *)malloc(buf_sz*sizeof(double));
      if(!buf){
//        fprintf(stderr, RCAT("memory allocation in ", AX_EQ_B_CHOL) "() failed!\n");
//        exit(1);
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"memory allocation in AX_EQ_B_CHOL");
      }
#endif /* LINSOLVERS_RETAIN_MEMORY */

    a=buf;
    b=a+a_sz;

  /* store A (column major!) into a anb B into b */
	for(i=0; i<m; i++){
		for(j=0; j<m; j++)
			a[i+j*m]=A[i*m+j];

    b[i]=B[i];
  }

  /* Cholesky decomposition of A */
  POTF2("U", (int *)&m, a, (int *)&m, (int *)&info);
  /* error treatment */
  if(info!=0){
    if(info<0){
//      fprintf(stderr, RCAT(RCAT("LAPACK error: illegal value for argument %d of ", POTF2) " in ", AX_EQ_B_CHOL) "()\n", -info);
//      exit(1);
	  ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"LAPACK error: illegal value for argument");
    }
    else{
      fprintf(stderr, RCAT(RCAT("LAPACK error: the leading minor of order %d is not positive definite,\nthe factorization could not be completed for ", POTF2) " in ", AX_EQ_B_CHOL) "()\n", info);
#ifndef LINSOLVERS_RETAIN_MEMORY
      free(buf);
#endif

      return 0;
    }
  }

  /* solve the linear system U^T y = b */
  TRTRS("U", "T", "N", (int *)&m, (int *)&nrhs, a, (int *)&m, b, (int *)&m, &info);
  /* error treatment */
  if(info!=0){
    if(info<0){
//      fprintf(stderr, RCAT(RCAT("LAPACK error: illegal value for argument %d of ", TRTRS) " in ", AX_EQ_B_CHOL) "()\n", -info);
//      exit(1);
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"LAPACK error: illegal value for argument");
    }
    else{
      fprintf(stderr, RCAT("LAPACK error: the %d-th diagonal element of A is zero (singular matrix) in ", AX_EQ_B_CHOL) "()\n", info);
#ifndef LINSOLVERS_RETAIN_MEMORY
      free(buf);
#endif

      return 0;
    }
  }

  /* solve the linear system U x = y */
  TRTRS("U", "N", "N", (int *)&m, (int *)&nrhs, a, (int *)&m, b, (int *)&m, &info);
  /* error treatment */
  if(info!=0){
    if(info<0){
//      fprintf(stderr, RCAT(RCAT("LAPACK error: illegal value for argument %d of ", TRTRS) "in ", AX_EQ_B_CHOL) "()\n", -info);
//      exit(1);
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"LAPACK error: illegal value for argument");
    }
    else{
      fprintf(stderr, RCAT("LAPACK error: the %d-th diagonal element of A is zero (singular matrix) in ", AX_EQ_B_CHOL) "()\n", info);
#ifndef LINSOLVERS_RETAIN_MEMORY
      free(buf);
#endif

      return 0;
    }
  }

	/* copy the result in x */
	for(i=0; i<m; i++)
    x[i]=b[i];

#ifndef LINSOLVERS_RETAIN_MEMORY
  free(buf);
#endif

	return 1;
}

/*
 * This function returns the solution of Ax = b
 *
 * The function employs LU decomposition:
 * If A=L U with L lower and U upper triangular, then the original system
 * amounts to solving
 * L y = b, U x = y
 *
 * A is mxm, b is mx1
 *
 * The function returns 0 in case of error,
 * 1 if successfull
 *
 * This function is often called repetitively to solve problems of identical
 * dimensions. To avoid repetitive malloc's and free's, allocated memory is
 * retained between calls and free'd-malloc'ed when not of the appropriate size.
 * A call with NULL as the first argument forces this memory to be released.
 */
int AX_EQ_B_LU(double *A, double *B, double *x, int m)
{
__STATIC__ double *buf=NULL;
__STATIC__ int buf_sz=0;

int a_sz, ipiv_sz, b_sz, work_sz, tot_sz;
register int i, j;
int info, *ipiv, nrhs=1;
double *a, *b, *work;
   
#ifdef LINSOLVERS_RETAIN_MEMORY
    if(!A){
      if(buf) free(buf);
      buf_sz=0;
      return 1;
    }
#endif /* LINSOLVERS_RETAIN_MEMORY */
   
    /* calculate required memory size */
    ipiv_sz=m;
    a_sz=m*m;
    b_sz=m;
    work_sz=100*m; /* this is probably too much */
    tot_sz=ipiv_sz + a_sz + b_sz + work_sz; // ipiv_sz counted as double here, no harm is done though

#ifdef LINSOLVERS_RETAIN_MEMORY
    if(tot_sz>buf_sz){ /* insufficient memory, allocate a "big" memory chunk at once */
      if(buf) free(buf); /* free previously allocated memory */

      buf_sz=tot_sz;
      buf=(double *)malloc(buf_sz*sizeof(double));
      if(!buf){
//        fprintf(stderr, RCAT("memory allocation in ", AX_EQ_B_LU) "() failed!\n");
//        exit(1);
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"memory allocation in AX_EQ_B_LU");
      }
    }
#else
      buf_sz=tot_sz;
      buf=(double *)malloc(buf_sz*sizeof(double));
      if(!buf){
//        fprintf(stderr, RCAT("memory allocation in ", AX_EQ_B_LU) "() failed!\n");
//        exit(1);
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"memory allocation in AX_EQ_B_LU");
      }
#endif /* LINSOLVERS_RETAIN_MEMORY */

    ipiv=(int *)buf;
    a=(double *)(ipiv + ipiv_sz);
    b=a+a_sz;
    work=b+b_sz;

    /* store A (column major!) into a and B into b */
	  for(i=0; i<m; i++){
		  for(j=0; j<m; j++)
        a[i+j*m]=A[i*m+j];

      b[i]=B[i];
    }

  /* LU decomposition for A */
	GETRF((int *)&m, (int *)&m, a, (int *)&m, ipiv, (int *)&info);  
	if(info!=0){
		if(info<0){
//			fprintf(stderr, RCAT(RCAT("argument %d of ", GETRF) " illegal in ", AX_EQ_B_LU) "()\n", -info);
//			exit(1);
			ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"argument illegal in AX_EQ_B_LU");
		}
		else{
      fprintf(stderr, RCAT(RCAT("singular matrix A for ", GETRF) " in ", AX_EQ_B_LU) "()\n");
#ifndef LINSOLVERS_RETAIN_MEMORY
      free(buf);
#endif

			return 0;
		}
	}

  /* solve the system with the computed LU */
  GETRS("N", (int *)&m, (int *)&nrhs, a, (int *)&m, ipiv, b, (int *)&m, (int *)&info);
	if(info!=0){
		if(info<0){
//			fprintf(stderr, RCAT(RCAT("argument %d of ", GETRS) " illegal in ", AX_EQ_B_LU) "()\n", -info);
//			exit(1);
			ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"argument illegal in AX_EQ_B_LU");
		}
		else{
			fprintf(stderr, RCAT(RCAT("unknown error for ", GETRS) " in ", AX_EQ_B_LU) "()\n");
#ifndef LINSOLVERS_RETAIN_MEMORY
      free(buf);
#endif

			return 0;
		}
	}

	/* copy the result in x */
	for(i=0; i<m; i++){
		x[i]=b[i];
	}

#ifndef LINSOLVERS_RETAIN_MEMORY
  free(buf);
#endif

	return 1;
}

/*
 * This function returns the solution of Ax = b
 *
 * The function is based on SVD decomposition:
 * If A=U D V^T with U, V orthogonal and D diagonal, the linear system becomes
 * (U D V^T) x = b or x=V D^{-1} U^T b
 * Note that V D^{-1} U^T is the pseudoinverse A^+
 *
 * A is mxm, b is mx1.
 *
 * The function returns 0 in case of error, 1 if successfull
 *
 * This function is often called repetitively to solve problems of identical
 * dimensions. To avoid repetitive malloc's and free's, allocated memory is
 * retained between calls and free'd-malloc'ed when not of the appropriate size.
 * A call with NULL as the first argument forces this memory to be released.
 */
int AX_EQ_B_SVD(double *A, double *B, double *x, int m)
{
__STATIC__ double *buf=NULL;
__STATIC__ int buf_sz=0;
static double eps=(-1.0);

register int i, j;
double *a, *u, *s, *vt, *work;
int a_sz, u_sz, s_sz, vt_sz, tot_sz;
double thresh, one_over_denom;
register double sum;
int info, rank, worksz, *iwork, iworksz;
   
#ifdef LINSOLVERS_RETAIN_MEMORY
    if(!A){
      if(buf) free(buf);
      buf_sz=0;
      return 1;
    }
#endif /* LINSOLVERS_RETAIN_MEMORY */
   
  /* calculate required memory size */
  worksz=16*m; /* more than needed */
  iworksz=8*m;
  a_sz=m*m;
  u_sz=m*m; s_sz=m; vt_sz=m*m;

  tot_sz=iworksz*sizeof(int) + (a_sz + u_sz + s_sz + vt_sz + worksz)*sizeof(double);

#ifdef LINSOLVERS_RETAIN_MEMORY
  if(tot_sz>buf_sz){ /* insufficient memory, allocate a "big" memory chunk at once */
    if(buf) free(buf); /* free previously allocated memory */

    buf_sz=tot_sz;
    buf=(double *)malloc(buf_sz);
    if(!buf){
//      fprintf(stderr, RCAT("memory allocation in ", AX_EQ_B_SVD) "() failed!\n");
//      exit(1);
	  ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"memory allocation in AX_EQ_B_SVD");
    }
  }
#else
    buf_sz=tot_sz;
    buf=(double *)malloc(buf_sz);
    if(!buf){
//      fprintf(stderr, RCAT("memory allocation in ", AX_EQ_B_SVD) "() failed!\n");
//      exit(1);
	  ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"memory allocation in AX_EQ_B_SVD");
      }
#endif /* LINSOLVERS_RETAIN_MEMORY */

  iwork=(int *)buf;
  a=(double *)(iwork+iworksz);
  /* store A (column major!) into a */
  for(i=0; i<m; i++)
    for(j=0; j<m; j++)
      a[i+j*m]=A[i*m+j];

  u=a + a_sz;
  s=u+u_sz;
  vt=s+s_sz;
  work=vt+vt_sz;

  /* SVD decomposition of A */
  GESVD("A", "A", (int *)&m, (int *)&m, a, (int *)&m, s, u, (int *)&m, vt, (int *)&m, work, (int *)&worksz, &info);
  //GESDD("A", (int *)&m, (int *)&m, a, (int *)&m, s, u, (int *)&m, vt, (int *)&m, work, (int *)&worksz, iwork, &info);

  /* error treatment */
  if(info!=0){
    if(info<0){
//      fprintf(stderr, RCAT(RCAT(RCAT("LAPACK error: illegal value for argument %d of ", GESVD), "/" GESDD) " in ", AX_EQ_B_SVD) "()\n", -info);
//      exit(1);
	  ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"LAPACK error: illegal value for argument");
      }
    else{
      fprintf(stderr, RCAT("LAPACK error: dgesdd (dbdsdc)/dgesvd (dbdsqr) failed to converge in ", AX_EQ_B_SVD) "() [info=%d]\n", info);
#ifndef LINSOLVERS_RETAIN_MEMORY
      free(buf);
#endif

      return 0;
    }
  }

  if(eps<0.0){
    double aux;

    /* compute machine epsilon */
    for(eps=(1.0); aux=eps+(1.0), aux-(1.0)>0.0; eps*=(0.5))
                                          ;
    eps*=(2.0);
  }

  /* compute the pseudoinverse in a */
	for(i=0; i<a_sz; i++) a[i]=0.0; /* initialize to zero */
  for(rank=0, thresh=eps*s[0]; rank<m && s[rank]>thresh; rank++){
    one_over_denom=(1.0)/s[rank];

    for(j=0; j<m; j++)
      for(i=0; i<m; i++)
        a[i*m+j]+=vt[rank+i*m]*u[j+rank*m]*one_over_denom;
  }

	/* compute A^+ b in x */
	for(i=0; i<m; i++){
	  for(j=0, sum=0.0; j<m; j++)
      sum+=a[i*m+j]*B[j];
    x[i]=sum;
  }

#ifndef LINSOLVERS_RETAIN_MEMORY
  free(buf);
#endif

	return 1;
}

/* undefine all. IT MUST REMAIN IN THIS POSITION IN FILE */
#undef AX_EQ_B_QR
#undef AX_EQ_B_QRLS
#undef AX_EQ_B_CHOL
#undef AX_EQ_B_LU
#undef AX_EQ_B_SVD

#undef GEQRF
#undef ORGQR
#undef TRTRS
#undef POTF2
#undef POTRF
#undef GETRF
#undef GETRS
#undef GESVD
#undef GESDD

#else // no LAPACK

/* precision-specific definitions */
#define AX_EQ_B_LU dAx_eq_b_LU_noLapack

/*
 * This function returns the solution of Ax = b
 *
 * The function employs LU decomposition followed by forward/back substitution (see 
 * also the LAPACK-based LU solver above)
 *
 * A is mxm, b is mx1
 *
 * The function returns 0 in case of error,
 * 1 if successfull
 *
 * This function is often called repetitively to solve problems of identical
 * dimensions. To avoid repetitive malloc's and free's, allocated memory is
 * retained between calls and free'd-malloc'ed when not of the appropriate size.
 * A call with NULL as the first argument forces this memory to be released.
 */
int AX_EQ_B_LU(double *A, double *B, double *x, int m)
{
__STATIC__ void *buf=NULL;
__STATIC__ int buf_sz=0;

int i, j, k;
int *idx, maxi=-1, idx_sz, a_sz, work_sz, tot_sz;
double *a, *work, max, sum, tmp;

#ifdef LINSOLVERS_RETAIN_MEMORY
    if(!A){
      if(buf) free(buf);
      buf_sz=0;
      return 1;
    }
#endif /* LINSOLVERS_RETAIN_MEMORY */
   
  /* calculate required memory size */
  idx_sz=m;
  a_sz=m*m;
  work_sz=m;
  tot_sz=idx_sz*sizeof(int) + (a_sz+work_sz)*sizeof(double);

#ifdef LINSOLVERS_RETAIN_MEMORY
  if(tot_sz>buf_sz){ /* insufficient memory, allocate a "big" memory chunk at once */
    if(buf) free(buf); /* free previously allocated memory */

    buf_sz=tot_sz;
    buf=(void *)malloc(tot_sz);
    if(!buf){
//      fprintf(stderr, RCAT("memory allocation in ", AX_EQ_B_LU) "() failed!\n");
//      exit(1);
	  ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"memory allocation in AX_EQ_B_LU");
    }
  }
#else
    buf_sz=tot_sz;
    buf=(void *)malloc(tot_sz);
    if(!buf){
//      fprintf(stderr, RCAT("memory allocation in ", AX_EQ_B_LU) "() failed!\n");
//      exit(1);
	  ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"memory allocation in AX_EQ_B_LU");
    }
#endif /* LINSOLVERS_RETAIN_MEMORY */

  idx=(int *)buf;
  a=(double *)(idx + idx_sz);
  work=a + a_sz;

  /* avoid destroying A, B by copying them to a, x resp. */
  for(i=0; i<m; ++i){ // B & 1st row of A
    a[i]=A[i];
    x[i]=B[i];
  }
  for(  ; i<a_sz; ++i) a[i]=A[i]; // copy A's remaining rows
  /****
  for(i=0; i<m; ++i){
    for(j=0; j<m; ++j)
      a[i*m+j]=A[i*m+j];
    x[i]=B[i];
  }
  ****/

  /* compute the LU decomposition of a row permutation of matrix a; the permutation itself is saved in idx[] */
	for(i=0; i<m; ++i){
		max=0.0;
		for(j=0; j<m; ++j)
			if((tmp=FABS(a[i*m+j]))>max)
        max=tmp;
		  if(max==0.0){
        fprintf(stderr, RCAT("Singular matrix A in ", AX_EQ_B_LU) "()!\n");
#ifndef LINSOLVERS_RETAIN_MEMORY
        free(buf);
#endif

        return 0;
      }
		  work[i]=(1.0)/max;
	}

	for(j=0; j<m; ++j){
		for(i=0; i<j; ++i){
			sum=a[i*m+j];
			for(k=0; k<i; ++k)
        sum-=a[i*m+k]*a[k*m+j];
			a[i*m+j]=sum;
		}
		max=0.0;
		for(i=j; i<m; ++i){
			sum=a[i*m+j];
			for(k=0; k<j; ++k)
        sum-=a[i*m+k]*a[k*m+j];
			a[i*m+j]=sum;
			if((tmp=work[i]*FABS(sum))>=max){
				max=tmp;
				maxi=i;
			}
		}
		if(j!=maxi){
			for(k=0; k<m; ++k){
				tmp=a[maxi*m+k];
				a[maxi*m+k]=a[j*m+k];
				a[j*m+k]=tmp;
			}
			work[maxi]=work[j];
		}
		idx[j]=maxi;
		if(a[j*m+j]==0.0)
      a[j*m+j]=DBL_EPSILON;
		if(j!=m-1){
			tmp=(1.0)/(a[j*m+j]);
			for(i=j+1; i<m; ++i)
        a[i*m+j]*=tmp;
		}
	}

  /* The decomposition has now replaced a. Solve the linear system using
   * forward and back substitution
   */
	for(i=k=0; i<m; ++i){
		j=idx[i];
		sum=x[j];
		x[j]=x[i];
		if(k!=0)
			for(j=k-1; j<i; ++j)
        sum-=a[i*m+j]*x[j];
		else
      if(sum!=0.0)
			  k=i+1;
		x[i]=sum;
	}

	for(i=m-1; i>=0; --i){
		sum=x[i];
		for(j=i+1; j<m; ++j)
      sum-=a[i*m+j]*x[j];
		x[i]=sum/a[i*m+i];
	}

#ifndef LINSOLVERS_RETAIN_MEMORY
  free(buf);
#endif

  return 1;
}

/* undefine all. IT MUST REMAIN IN THIS POSITION IN FILE */
#undef AX_EQ_B_LU

#endif /* HAVE_LAPACK */


/////////////////////////////////////////////////////////////////////////////////
// 
//  Levenberg - Marquardt non-linear minimization algorithm
//  Copyright (C) 2004-05  Manolis Lourakis (lourakis@ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
/////////////////////////////////////////////////////////////////////////////////

/******************************************************************************** 
 * Box-constrained Levenberg-Marquardt nonlinear minimization. The same core code
 * is used with appropriate #defines to derive single and double precision versions,
 * see also lmbc_core.c
 ********************************************************************************/


#define EPSILON       1E-12
#define ONE_THIRD     0.3333333334 /* 1.0/3.0 */

#define FUNC_STATE LM_ADD_PREFIX(func_state)
#define LNSRCH LM_ADD_PREFIX(lnsrch)
#define BOXPROJECT LM_ADD_PREFIX(boxProject)
#define BOXCHECK LM_ADD_PREFIX(boxCHECK)

#define LEVMAR_BC_DER dlevmar_bc_der
#define LEVMAR_BC_DIF dlevmar_bc_dif //CHECKME
#define FDIF_FORW_JAC_APPROX dfdif_forw_jac_approx
#define FDIF_CENT_JAC_APPROX dfdif_cent_jac_approx
#define TRANS_MAT_MAT_MULT dtrans_mat_mat_mult
#define LEVMAR_COVAR dlevmar_covar

#ifdef HAVE_LAPACK
#define AX_EQ_B_LU dAx_eq_b_LU
#define AX_EQ_B_CHOL dAx_eq_b_Chol
#define AX_EQ_B_QR dAx_eq_b_QR
#define AX_EQ_B_QRLS dAx_eq_b_QRLS
#define AX_EQ_B_SVD dAx_eq_b_SVD
#else
#define AX_EQ_B_LU dAx_eq_b_LU_noLapack
#endif /* HAVE_LAPACK */

/* find the median of 3 numbers */
#define __MEDIAN3(a, b, c) ( ((a) >= (b))?\
        ( ((c) >= (a))? (a) : ( ((c) <= (b))? (b) : (c) ) ) : \
        ( ((c) >= (b))? (b) : ( ((c) <= (a))? (a) : (c) ) ) )

#define _POW_ (2.1)

struct FUNC_STATE{
  int n, *nfev;
  double *hx, *x;
  void *adata;
};

static void
LNSRCH(int m, double *x, double f, double *g, double *p, double alpha, double *xpls,
       double *ffpls, const ARM_LEVMARQFunc& func, struct FUNC_STATE state,
       int *mxtake, int *iretcd, double stepmx, double steptl, double *sx)
{
/* Find a next newton iterate by backtracking line search.
 * Specifically, finds a \lambda such that for a fixed alpha<0.5 (usually 1e-4),
 * f(x + \lambda*p) <= f(x) + alpha * \lambda * g^T*p
 *
 * Translated (with minor changes) from Schnabel, Koontz & Weiss uncmin.f,  v1.3

 * PARAMETERS :

 *	m       --> dimension of problem (i.e. number of variables)
 *	x(m)    --> old iterate:	x[k-1]
 *	f       --> function value at old iterate, f(x)
 *	g(m)    --> gradient at old iterate, g(x), or approximate
 *	p(m)    --> non-zero newton step
 *	alpha   --> fixed constant < 0.5 for line search (see above)
 *	xpls(m) <--	 new iterate x[k]
 *	ffpls   <--	 function value at new iterate, f(xpls)
 *	func    --> name of subroutine to evaluate function
 *	state   <--> information other than x and m that func requires.
 *			    state is not modified in xlnsrch (but can be modified by func).
 *	iretcd  <--	 return code
 *	mxtake  <--	 boolean flag indicating step of maximum length used
 *	stepmx  --> maximum allowable step size
 *	steptl  --> relative step size at which successive iterates
 *			    considered close enough to terminate algorithm
 *	sx(m)	  --> diagonal scaling matrix for x, can be NULL

 *	internal variables

 *	sln		 newton length
 *	rln		 relative length of newton step
*/

    int i;
    int firstback = 1;
    double disc;
    double a3, b;
    double t1, t2, t3, lambda, tlmbda, rmnlmb;
    double scl, rln, sln, slp;
    double tmp1, tmp2;
    double fpls, pfpls = 0., plmbda = 0.; /* -Wall */

    f*=(0.5);
    *mxtake = 0;
    *iretcd = 2;
    tmp1 = 0.;
    if(!sx) /* no scaling */
      for (i = 0; i < m; ++i)
        tmp1 += p[i] * p[i];
    else
      for (i = 0; i < m; ++i)
        tmp1 += sx[i] * sx[i] * p[i] * p[i];
    sln = (double)sqrt(tmp1);
    if (sln > stepmx) {
	  /*	newton step longer than maximum allowed */
	    scl = stepmx / sln;
      for(i=0; i<m; ++i) /* p * scl */
        p[i]*=scl;
	    sln = stepmx;
    }
    for(i=0, slp=0.; i<m; ++i) /* g^T * p */
      slp+=g[i]*p[i];
    rln = 0.;
    if(!sx) /* no scaling */
      for (i = 0; i < m; ++i) {
	      tmp1 = (FABS(x[i])>=(1.))? FABS(x[i]) : (1.);
	      tmp2 = FABS(p[i])/tmp1;
	      if(rln < tmp2) rln = tmp2;
      }
    else
      for (i = 0; i < m; ++i) {
	      tmp1 = (FABS(x[i])>=(1.)/sx[i])? FABS(x[i]) : (1.)/sx[i];
	      tmp2 = FABS(p[i])/tmp1;
	      if(rln < tmp2) rln = tmp2;
      }
    rmnlmb = steptl / rln;
    lambda = (1.0);

    /*	check if new iterate satisfactory.  generate new lambda if necessary. */

    while(*iretcd > 1) {
	    for (i = 0; i < m; ++i)
	      xpls[i] = x[i] + lambda * p[i];

      /* evaluate function at new point */
      func(xpls, state.hx, m, state.n, state.adata);
      for(i=0, tmp1=0.0; i<state.n; ++i){
        state.hx[i]=tmp2=state.x[i]-state.hx[i];
        tmp1+=tmp2*tmp2;
      }
      fpls=(0.5)*tmp1; *ffpls=tmp1; ++(*(state.nfev));

	    if (fpls <= f + slp * alpha * lambda) { /* solution found */
	      *iretcd = 0;
	      if (lambda == (1.) && sln > stepmx * (.99)) *mxtake = 1;
	      return;
	    }

	    /* else : solution not (yet) found */

      /* First find a point with a finite value */

	    if (lambda < rmnlmb) {
	      /* no satisfactory xpls found sufficiently distinct from x */

	      *iretcd = 1;
	      return;
	    }
	    else { /*	calculate new lambda */

	      /* modifications to cover non-finite values */
	      if (fpls >= DBL_MAX) {
		      lambda *= (0.1);
		      firstback = 1;
	      }
	      else {
		      if (firstback) { /*	first backtrack: quadratic fit */
		        tlmbda = -lambda * slp / ((fpls - f - slp) * (2.));
		        firstback = 0;
		      }
		      else { /*	all subsequent backtracks: cubic fit */
		        t1 = fpls - f - lambda * slp;
		        t2 = pfpls - f - plmbda * slp;
		        t3 = (1.) / (lambda - plmbda);
		        a3 = (3.) * t3 * (t1 / (lambda * lambda)
				      - t2 / (plmbda * plmbda));
		        b = t3 * (t2 * lambda / (plmbda * plmbda)
			          - t1 * plmbda / (lambda * lambda));
		        disc = b * b - a3 * slp;
		        if (disc > b * b)
			    /* only one positive critical point, must be minimum */
			        tlmbda = (-b + ((a3 < 0)? -(double)sqrt(disc): (double)sqrt(disc))) /a3;
		        else
			    /* both critical points positive, first is minimum */
			        tlmbda = (-b + ((a3 < 0)? (double)sqrt(disc): -(double)sqrt(disc))) /a3;

		        if (tlmbda > lambda * (.5))
			        tlmbda = lambda * (.5);
		    }
		    plmbda = lambda;
		    pfpls = fpls;
		    if (tlmbda < lambda * (.1))
		      lambda *= (.1);
		    else
		      lambda = tlmbda;
      }
	  }
  }
} /* LNSRCH */

static void BOXPROJECT(double *p, double *lb, double *ub, int m)
{
register int i;

  if(!lb){ /* no lower bounds */
    if(!ub) /* no upper bounds */
      return;
    else{ /* upper bounds only */
      for(i=0; i<m; ++i)
        if(p[i]>ub[i]) p[i]=ub[i];
    }
  }
  else
    if(!ub){ /* lower bounds only */
      for(i=0; i<m; ++i)
        if(p[i]<lb[i]) p[i]=lb[i];
    }
    else /* box bounds */
      for(i=0; i<m; ++i)
        p[i]=__MEDIAN3(lb[i], p[i], ub[i]);
}


static int BOXCHECK(double *lb, double *ub, int m)
{
register int i;

  if(!lb || !ub) return 1;

  for(i=0; i<m; ++i)
    if(lb[i]>ub[i]) return 0;

  return 1;
}


/* 
 * This function seeks the parameter vector p that best describes the measurements
 * vector x under box constraints.
 * More precisely, given a vector function  func : R^m --> R^n with n>=m,
 * it finds p s.t. func(p) ~= x, i.e. the squared second order (i.e. L2) norm of
 * e=x-func(p) is minimized under the constraints lb[i]<=p[i]<=ub[i].
 * If no lower bound constraint applies for p[i], use -DBL_MAX/-FLT_MAX for lb[i];
 * If no upper bound constraint applies for p[i], use DBL_MAX/FLT_MAX for ub[i].
 *
 * This function requires an analytic jacobian. In case the latter is unavailable,
 * use LEVMAR_BC_DIF() bellow
 *
 * Returns the number of iterations (>=0) if successfull, -1 if failed
 *
 * For details, see C. Kanzow, N. Yamashita and M. Fukushima: "Levenberg-Marquardt
 * methods for constrained nonlinear equations with strong local convergence properties",
 * Journal of Computational and Applied Mathematics 172, 2004, pp. 375-397.
 * Also, see H.B. Nielsen's (http://www.imm.dtu.dk/~hbn) IMM/DTU tutorial on
 * unconrstrained Levenberg-Marquardt at http://www.imm.dtu.dk/courses/02611/nllsq.pdf
 */

int LEVMAR_BC_DER(
  const ARM_LEVMARQFunc& func, /* functional relation describing measurements. A p \in R^m yields a \hat{x} \in  R^n */
  const ARM_LEVMARQFunc& jacf,  /* function to evaluate the jacobian \part x / \part p */ 
  double *p,         /* I/O: initial parameter estimates. On output has the estimated solution */
  double *x,         /* I: measurement vector */
  int m,              /* I: parameter vector dimension (i.e. #unknowns) */
  int n,              /* I: measurement vector dimension */
  double *lb,        /* I: vector of lower bounds. If NULL, no lower bounds apply */
  double *ub,        /* I: vector of upper bounds. If NULL, no upper bounds apply */
  int itmax,          /* I: maximum number of iterations */
  double opts[4],    /* I: minim. options [\mu, \epsilon1, \epsilon2, \epsilon3]. Respectively the scale factor for initial \mu,
                       * stopping thresholds for ||J^T e||_inf, ||Dp||_2 and ||e||_2. Set to NULL for defaults to be used.
                       * Note that ||J^T e||_inf is computed on free (not equal to lb[i] or ub[i]) variables only.
                       */
  double info[LM_INFO_SZ],
					           /* O: information regarding the minimization. Set to NULL if don't care
                      * info[0]= ||e||_2 at initial p.
                      * info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||Dp||_2, mu/max[J^T J]_ii ], all computed at estimated p.
                      * info[5]= # iterations,
                      * info[6]=reason for terminating: 1 - stopped by small gradient J^T e
                      *                                 2 - stopped by small Dp
                      *                                 3 - stopped by itmax
                      *                                 4 - singular matrix. Restart from current p with increased mu 
                      *                                 5 - no further error reduction is possible. Restart with increased mu
                      *                                 6 - stopped by small ||e||_2
                      * info[7]= # function evaluations
                      * info[8]= # jacobian evaluations
                      */
  double *work,     /* working memory, allocate if NULL */
  double *covar,    /* O: Covariance matrix corresponding to LS solution; mxm. Set to NULL if not needed. */
  void *adata)       /* pointer to possibly additional data, passed uninterpreted to func & jacf.
                      * Set to NULL if not needed
                      */
{
register int i, j, k, l;
int worksz, freework=0, issolved;
/* temp work arrays */
double *e,          /* nx1 */
       *hx,         /* \hat{x}_i, nx1 */
       *jacTe,      /* J^T e_i mx1 */
       *jac,        /* nxm */
       *jacTjac,    /* mxm */
       *Dp,         /* mx1 */
   *diag_jacTjac,   /* diagonal of J^T J, mx1 */
       *pDp;        /* p + Dp, mx1 */

register double mu,  /* damping constant */
                tmp; /* mainly used in matrix & vector multiplications */
double p_eL2, jacTe_inf, pDp_eL2; /* ||e(p)||_2, ||J^T e||_inf, ||e(p+Dp)||_2 */
double p_L2, Dp_L2=DBL_MAX, dF, dL;
double tau, eps1, eps2, eps2_sq, eps3;
double init_p_eL2;
int nu=2, nu2, stop, nfev, njev=0;
const int nm=n*m;

/* variables for constrained LM */
struct FUNC_STATE fstate;
double alpha=(1e-4), beta=(0.9), gamma=(0.99995), gamma_sq=gamma*gamma, rho=(1e-8);
double t, t0;
double steptl=(1e3)*(double)sqrt(DBL_EPSILON), jacTeDp;
double tmin=(1e-12), tming=(1e-18); /* minimum step length for LS and PG steps */
const double tini=(1.0); /* initial step length for LS and PG steps */
int nLMsteps=0, nLSsteps=0, nPGsteps=0, gprevtaken=0;
int numactive;

  mu=jacTe_inf=t=0.0;  tmin=tmin; /* -Wall */

  if(n<m){
//    fprintf(stderr, LCAT(LEVMAR_BC_DER, "(): cannot solve a problem with fewer measurements [%d] than unknowns [%d]\n"), n, m);
//    exit(1);
	  ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"LEVMAR_BC_DER : cannot solve a problem with fewer measurements than unknowns");
  }

  //if(!jacf){
  //  fprintf(stderr, RCAT("No function specified for computing the jacobian in ", LEVMAR_BC_DER)
  //      RCAT("().\nIf no such function is available, use ", LEVMAR_BC_DIF) RCAT("() rather than ", LEVMAR_BC_DER) "()\n");
  //  exit(1);
  //}

  if(!BOXCHECK(lb, ub, m)){
//    fprintf(stderr, LCAT(LEVMAR_BC_DER, "(): at least one lower bound exceeds the upper one\n"));
//    exit(1);
	  ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"LEVMAR_BC_DER : at least one lower bound exceeds the upper one");
  }

  if(opts){
	  tau=opts[0];
	  eps1=opts[1];
	  eps2=opts[2];
	  eps2_sq=opts[2]*opts[2];
	  eps3=opts[3];
  }
  else{ // use default values
	  tau=(LM_INIT_MU);
	  eps1=(LM_STOP_THRESH);
	  eps2=(LM_STOP_THRESH);
	  eps2_sq=(LM_STOP_THRESH)*(LM_STOP_THRESH);
	  eps3=(LM_STOP_THRESH);
  }

  if(!work){
    worksz=LM_DER_WORKSZ(m, n); //2*n+4*m + n*m + m*m;
    work=(double *)malloc(worksz*sizeof(double)); /* allocate a big chunk in one step */
    if(!work){
//      fprintf(stderr, LCAT(LEVMAR_BC_DER, "(): memory allocation request failed\n"));
//      exit(1);
	  ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"LEVMAR_BC_DER : memory allocation request failed");
    }
    freework=1;
  }

  /* set up work arrays */
  e=work;
  hx=e + n;
  jacTe=hx + n;
  jac=jacTe + m;
  jacTjac=jac + nm;
  Dp=jacTjac + m*m;
  diag_jacTjac=Dp + m;
  pDp=diag_jacTjac + m;

  fstate.n=n;
  fstate.hx=hx;
  fstate.x=x;
  fstate.adata=adata;
  fstate.nfev=&nfev;
  
  /* see if starting point is within the feasile set */
  for(i=0; i<m; ++i)
    pDp[i]=p[i];
  BOXPROJECT(p, lb, ub, m); /* project to feasible set */
  for(i=0; i<m; ++i)
    if(pDp[i]!=p[i])
      fprintf(stderr, RCAT("Warning: component %d of starting point not feasible in ", LEVMAR_BC_DER) "()! [%g projected to %g]\n",
                      i, p[i], pDp[i]);

  /* compute e=x - f(p) and its L2 norm */
  func(p, hx, m, n, adata); nfev=1;
  for(i=0, p_eL2=0.0; i<n; ++i){
    e[i]=tmp=x[i]-hx[i];
    p_eL2+=tmp*tmp;
  }
  init_p_eL2=p_eL2;

  for(k=stop=0; k<itmax && !stop; ++k){
 //printf("%d  %.15g\n", k, 0.5*p_eL2);
    /* Note that p and e have been updated at a previous iteration */

    if(p_eL2<=eps3){ /* error is small */
      stop=6;
      break;
    }

    /* Compute the jacobian J at p,  J^T J,  J^T e,  ||J^T e||_inf and ||p||^2.
     * Since J^T J is symmetric, its computation can be speeded up by computing
     * only its upper triangular part and copying it to the lower part
     */

    jacf(p, jac, m, n, adata); ++njev;

    /* J^T J, J^T e */
    if(nm<__BLOCKSZ__SQ){ // this is a small problem
      /* This is the straightforward way to compute J^T J, J^T e. However, due to
       * its noncontinuous memory access pattern, it incures many cache misses when
       * applied to large minimization problems (i.e. problems involving a large
       * number of free variables and measurements), in which J is too large to
       * fit in the L1 cache. For such problems, a cache-efficient blocking scheme
       * is preferable.
       *
       * Thanks to John Nitao of Lawrence Livermore Lab for pointing out this
       * performance problem.
       *
       * On the other hand, the straightforward algorithm is faster on small
       * problems since in this case it avoids the overheads of blocking. 
       */

      for(i=0; i<m; ++i){
        for(j=i; j<m; ++j){
          int lm;

          for(l=0, tmp=0.0; l<n; ++l){
            lm=l*m;
            tmp+=jac[lm+i]*jac[lm+j];
          }

		      /* store tmp in the corresponding upper and lower part elements */
          jacTjac[i*m+j]=jacTjac[j*m+i]=tmp;
        }

        /* J^T e */
        for(l=0, tmp=0.0; l<n; ++l)
          tmp+=jac[l*m+i]*e[l];
        jacTe[i]=tmp;
      }
    }
    else{ // this is a large problem
      /* Cache efficient computation of J^T J based on blocking
       */
      dtrans_mat_mat_mult(jac, jacTjac, n, m, __BLOCKSZ__);

      /* cache efficient computation of J^T e */
      for(i=0; i<m; ++i)
        jacTe[i]=0.0;

      for(i=0; i<n; ++i){
        register double *jacrow;

        for(l=0, jacrow=jac+i*m, tmp=e[i]; l<m; ++l)
          jacTe[l]+=jacrow[l]*tmp;
      }
    }

	  /* Compute ||J^T e||_inf and ||p||^2. Note that ||J^T e||_inf
     * is computed for free (i.e. inactive) variables only. 
     * At a local minimum, if p[i]==ub[i] then g[i]>0;
     * if p[i]==lb[i] g[i]<0; otherwise g[i]=0 
     */
    for(i=j=numactive=0, p_L2=jacTe_inf=0.0; i<m; ++i){
      if(ub && p[i]==ub[i]){ ++numactive; if(jacTe[i]>0.0) ++j; }
      else if(lb && p[i]==lb[i]){ ++numactive; if(jacTe[i]<0.0) ++j; }
      else if(jacTe_inf < (tmp=FABS(jacTe[i]))) jacTe_inf=tmp;

      diag_jacTjac[i]=jacTjac[i*m+i]; /* save diagonal entries so that augmentation can be later canceled */
      p_L2+=p[i]*p[i];
    }
    //p_L2=sqrt(p_L2);

#if 0
if(!(k%100)){
  printf("Current estimate: ");
  for(i=0; i<m; ++i)
    printf("%.9g ", p[i]);
  printf("-- errors %.9g %0.9g, #active %d [%d]\n", jacTe_inf, p_eL2, numactive, j);
}
#endif

    /* check for convergence */
    if(j==numactive && (jacTe_inf <= eps1)){
      Dp_L2=0.0; /* no increment for p in this case */
      stop=1;
      break;
    }

   /* compute initial damping factor */
    if(k==0){
      if(!lb && !ub){ /* no bounds */
        for(i=0, tmp=DBL_MIN; i<m; ++i)
          if(diag_jacTjac[i]>tmp) tmp=diag_jacTjac[i]; /* find max diagonal element */
        mu=tau*tmp;
      }
      else 
        mu=(0.5)*tau*p_eL2; /* use Kanzow's starting mu */
    }

    /* determine increment using a combination of adaptive damping, line search and projected gradient search */
    while(1){
      /* augment normal equations */
      for(i=0; i<m; ++i)
        jacTjac[i*m+i]+=mu;

      /* solve augmented equations */
#ifdef HAVE_LAPACK
      /* 5 alternatives are available: LU, Cholesky, 2 variants of QR decomposition and SVD.
       * Cholesky is the fastest but might be inaccurate; QR is slower but more accurate;
       * SVD is the slowest but most accurate; LU offers a tradeoff between accuracy and speed
       */

      issolved=AX_EQ_B_LU(jacTjac, jacTe, Dp, m);
      //issolved=AX_EQ_B_CHOL(jacTjac, jacTe, Dp, m);
      //issolved=AX_EQ_B_QR(jacTjac, jacTe, Dp, m);
      //issolved=AX_EQ_B_QRLS(jacTjac, jacTe, Dp, m, m);
      //issolved=AX_EQ_B_SVD(jacTjac, jacTe, Dp, m);

#else
      /* use the LU included with levmar */
      issolved=AX_EQ_B_LU(jacTjac, jacTe, Dp, m);
#endif /* HAVE_LAPACK */

      if(issolved){
        for(i=0; i<m; ++i)
          pDp[i]=p[i] + Dp[i];

        /* compute p's new estimate and ||Dp||^2 */
        BOXPROJECT(pDp, lb, ub, m); /* project to feasible set */
        for(i=0, Dp_L2=0.0; i<m; ++i){
          Dp[i]=tmp=pDp[i]-p[i];
          Dp_L2+=tmp*tmp;
        }
        //Dp_L2=sqrt(Dp_L2);

        if(Dp_L2<=eps2_sq*p_L2){ /* relative change in p is small, stop */
          stop=2;
          break;
        }

       if(Dp_L2>=(p_L2+eps2)/((EPSILON)*(EPSILON))){ /* almost singular */
         stop=4;
         break;
       }

        func(pDp, hx, m, n, adata); ++nfev; /* evaluate function at p + Dp */
        for(i=0, pDp_eL2=0.0; i<n; ++i){ /* compute ||e(pDp)||_2 */
          hx[i]=tmp=x[i]-hx[i];
          pDp_eL2+=tmp*tmp;
        }

        if(pDp_eL2<=gamma_sq*p_eL2){
          for(i=0, dL=0.0; i<m; ++i)
            dL+=Dp[i]*(mu*Dp[i]+jacTe[i]);

#if 1
          if(dL>0.0){
            dF=p_eL2-pDp_eL2;
            tmp=((2.0)*dF/dL-(1.0));
            tmp=(1.0)-tmp*tmp*tmp;
            mu=mu*( (tmp>=(ONE_THIRD))? tmp : (ONE_THIRD) );
          }
          else
            mu=(mu>=pDp_eL2)? pDp_eL2 : mu; /* pDp_eL2 is the new pDp_eL2 */
#else

          mu=(mu>=pDp_eL2)? pDp_eL2 : mu; /* pDp_eL2 is the new pDp_eL2 */
#endif

          nu=2;

          for(i=0 ; i<m; ++i) /* update p's estimate */
            p[i]=pDp[i];

          for(i=0; i<n; ++i) /* update e and ||e||_2 */
            e[i]=hx[i];
          p_eL2=pDp_eL2;
          ++nLMsteps;
          gprevtaken=0;
          break;
        }
      }
      else{

      /* the augmented linear system could not be solved, increase mu */

        mu*=nu;
        nu2=nu<<1; // 2*nu;
        if(nu2<=nu){ /* nu has wrapped around (overflown). Thanks to Frank Jordan for spotting this case */
          stop=5;
          break;
        }
        nu=nu2;

        for(i=0; i<m; ++i) /* restore diagonal J^T J entries */
          jacTjac[i*m+i]=diag_jacTjac[i];

        continue; /* solve again with increased nu */
      }

      /* if this point is reached, the LM step did not reduce the error;
       * see if it is a descent direction
       */

      /* negate jacTe (i.e. g) & compute g^T * Dp */
      for(i=0, jacTeDp=0.0; i<m; ++i){
        jacTe[i]=-jacTe[i];
        jacTeDp+=jacTe[i]*Dp[i];
      }

      if(jacTeDp<=-rho*pow(Dp_L2, _POW_/(2.0))){
        /* Dp is a descent direction; do a line search along it */
        int mxtake, iretcd;
        double stepmx;

        tmp=(double)sqrt(p_L2); stepmx=(1e3)*( (tmp>=(1.0))? tmp : (1.0) );

#if 1
        /* use Schnabel's backtracking line search; it requires fewer "func" evaluations */
        LNSRCH(m, p, p_eL2, jacTe, Dp, alpha, pDp, &pDp_eL2, func, fstate,
               &mxtake, &iretcd, stepmx, steptl, NULL); /* NOTE: LNSRCH() updates hx */
        if(iretcd!=0) goto gradproj; /* rather inelegant but effective way to handle LNSRCH() failures... */
#else
        /* use the simpler (but slower!) line search described by Kanzow */
        for(t=tini; t>tmin; t*=beta){
          for(i=0; i<m; ++i){
            pDp[i]=p[i] + t*Dp[i];
            //pDp[i]=__MEDIAN3(lb[i], pDp[i], ub[i]); /* project to feasible set */
          }

          func(pDp, hx, m, n, adata); ++nfev; /* evaluate function at p + t*Dp */
          for(i=0, pDp_eL2=0.0; i<n; ++i){ /* compute ||e(pDp)||_2 */
            hx[i]=tmp=x[i]-hx[i];
            pDp_eL2+=tmp*tmp;
          }
          //if((0.5)*pDp_eL2<=(0.5)*p_eL2 + t*alpha*jacTeDp) break;
          if(pDp_eL2<=p_eL2 + (2.0)*t*alpha*jacTeDp) break;
        }
#endif
        ++nLSsteps;
        gprevtaken=0;

        /* NOTE: new estimate for p is in pDp, associated error in hx and its norm in pDp_eL2.
         * These values are used below to update their corresponding variables 
         */
      }
      else{
gradproj: /* Note that this point can also be reached via a goto when LNSRCH() fails */

        /* jacTe is a descent direction; make a projected gradient step */

        /* if the previous step was along the gradient descent, try to use the t employed in that step */
        /* compute ||g|| */
        for(i=0, tmp=0.0; i<m; ++i)
          tmp=jacTe[i]*jacTe[i];
        tmp=(double)sqrt(tmp);
        tmp=(100.0)/((1.0)+tmp);
        t0=(tmp<=tini)? tmp : tini; /* guard against poor scaling & large steps; see (3.50) in C.T. Kelley's book */

        for(t=(gprevtaken)? t : t0; t>tming; t*=beta){
          for(i=0; i<m; ++i)
            pDp[i]=p[i] - t*jacTe[i];
          BOXPROJECT(pDp, lb, ub, m); /* project to feasible set */
          for(i=0; i<m; ++i)
            Dp[i]=pDp[i]-p[i];

          func(pDp, hx, m, n, adata); ++nfev; /* evaluate function at p - t*g */
          for(i=0, pDp_eL2=0.0; i<n; ++i){ /* compute ||e(pDp)||_2 */
            hx[i]=tmp=x[i]-hx[i];
            pDp_eL2+=tmp*tmp;
          }
          for(i=0, tmp=0.0; i<m; ++i) /* compute ||g^T * Dp|| */
            tmp+=jacTe[i]*Dp[i];

          if(gprevtaken && pDp_eL2<=p_eL2 + (2.0)*(0.99999)*tmp){ /* starting t too small */
            t=t0;
            gprevtaken=0;
            continue;
          }
          //if((0.5)*pDp_eL2<=(0.5)*p_eL2 + alpha*tmp) break;
          if(pDp_eL2<=p_eL2 + (2.0)*alpha*tmp) break;
        }

        ++nPGsteps;
        gprevtaken=1;
        /* NOTE: new estimate for p is in pDp, associated error in hx and its norm in pDp_eL2 */
      }

      /* update using computed values */

      for(i=0, Dp_L2=0.0; i<m; ++i){
        tmp=pDp[i]-p[i];
        Dp_L2+=tmp*tmp;
      }
      //Dp_L2=sqrt(Dp_L2);

      if(Dp_L2<=eps2_sq*p_L2){ /* relative change in p is small, stop */
        stop=2;
        break;
      }

      for(i=0 ; i<m; ++i) /* update p's estimate */
        p[i]=pDp[i];

      for(i=0; i<n; ++i) /* update e and ||e||_2 */
        e[i]=hx[i];
      p_eL2=pDp_eL2;
      break;
    } /* inner loop */
  }

  if(k>=itmax) stop=3;

  for(i=0; i<m; ++i) /* restore diagonal J^T J entries */
    jacTjac[i*m+i]=diag_jacTjac[i];

  if(info){
    info[0]=init_p_eL2;
    info[1]=p_eL2;
    info[2]=jacTe_inf;
    info[3]=Dp_L2;
    for(i=0, tmp=DBL_MIN; i<m; ++i)
      if(tmp<jacTjac[i*m+i]) tmp=jacTjac[i*m+i];
    info[4]=mu/tmp;
    info[5]=(double)k;
    info[6]=(double)stop;
    info[7]=(double)nfev;
    info[8]=(double)njev;
  }

  /* covariance matrix */
  if(covar){
    dlevmar_covar(jacTjac, covar, p_eL2, m, n);
  }
                                                               
  if(freework) free(work);

#if 0
printf("%d LM steps, %d line search, %d projected gradient\n", nLMsteps, nLSsteps, nPGsteps);
#endif

  return (stop!=4)?  k : -1;
}


/* No jacobian version of the LEVMAR_BC_DER() function above: the jacobian is approximated with 
 * the aid of finite differences (forward or central, see the comment for the opts argument)
 * Ideally, this function should be implemented with a secant approach. Currently, it just calls
 * LEVMAR_BC_DER()
 */
int LEVMAR_BC_DIF(
  const ARM_LEVMARQFunc& func, /* functional relation describing measurements. A p \in R^m yields a \hat{x} \in  R^n */
  double *p,         /* I/O: initial parameter estimates. On output has the estimated solution */
  double *x,         /* I: measurement vector */
  int m,              /* I: parameter vector dimension (i.e. #unknowns) */
  int n,              /* I: measurement vector dimension */
  double *lb,        /* I: vector of lower bounds. If NULL, no lower bounds apply */
  double *ub,        /* I: vector of upper bounds. If NULL, no upper bounds apply */
  int itmax,          /* I: maximum number of iterations */
  double opts[5],    /* I: opts[0-4] = minim. options [\mu, \epsilon1, \epsilon2, \epsilon3, \delta]. Respectively the
                       * scale factor for initial \mu, stopping thresholds for ||J^T e||_inf, ||Dp||_2 and ||e||_2 and
                       * the step used in difference approximation to the jacobian. Set to NULL for defaults to be used.
                       * If \delta<0, the jacobian is approximated with central differences which are more accurate
                       * (but slower!) compared to the forward differences employed by default. 
                       */
  double info[LM_INFO_SZ],
					           /* O: information regarding the minimization. Set to NULL if don't care
                      * info[0]= ||e||_2 at initial p.
                      * info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||Dp||_2, mu/max[J^T J]_ii ], all computed at estimated p.
                      * info[5]= # iterations,
                      * info[6]=reason for terminating: 1 - stopped by small gradient J^T e
                      *                                 2 - stopped by small Dp
                      *                                 3 - stopped by itmax
                      *                                 4 - singular matrix. Restart from current p with increased mu 
                      *                                 5 - no further error reduction is possible. Restart with increased mu
                      *                                 6 - stopped by small ||e||_2
                      * info[7]= # function evaluations
                      * info[8]= # jacobian evaluations
                      */
  double *work,     /* working memory, allocate if NULL */
  double *covar,    /* O: Covariance matrix corresponding to LS solution; mxm. Set to NULL if not needed. */
  void *adata)       /* pointer to possibly additional data, passed uninterpreted to func.
                      * Set to NULL if not needed
                      */
{
struct LMBC_DIF_DATA data;
int ret;

    //fprintf(stderr, RCAT("\nWarning: current implementation of ", LEVMAR_BC_DIF) "() does not use a secant approach!\n\n");

    data.func=&func;
    data.hx=(double *)malloc(2*n*sizeof(double)); /* allocate a big chunk in one step */
    if(!data.hx){
//      fprintf(stderr, LCAT(LEVMAR_BC_DIF, "(): memory allocation request failed\n"));
//      exit(1);
	  ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"LEVMAR_BC_DIF : memory allocation request failed");
    }
    data.hxx=data.hx+n;
    data.adata=adata;
    data.delta=(opts)? FABS(opts[4]) : (double)LM_DIFF_DELTA; // no central differences here...

	ARM_LEVMARQDIFFunc nfunc;
	ARM_LEVMARQDIFJACFunc njacf;

    ret=LEVMAR_BC_DER(nfunc, njacf, p, x, m, n, lb, ub, itmax, opts, info, work, covar, (void *)&data);

    if(info) /* correct the number of function calls */
      info[7]+=info[8]*(m+1); /* each jacobian evaluation costs m+1 function calls */

    free(data.hx);

    return ret;
}

#undef AX_EQ_B_LU
#undef AX_EQ_B_CHOL
#undef AX_EQ_B_QR
#undef AX_EQ_B_QRLS
#undef AX_EQ_B_SVD

/////////////////////////////////////////////////////////////////////////////////
// 
//  Levenberg - Marquardt non-linear minimization algorithm
//  Copyright (C) 2004-05  Manolis Lourakis (lourakis@ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
/////////////////////////////////////////////////////////////////////////////////

#ifdef HAVE_LAPACK
#define LEVMAR_PSEUDOINVERSE dlevmar_pseudoinverse
static int LEVMAR_PSEUDOINVERSE(double *A, double *B, int m);

/* LAPACK SVD routines */
#define GESVD LM_ADD_PREFIX(gesvd_)
#define GESDD LM_ADD_PREFIX(gesdd_)
extern int GESVD(char *jobu, char *jobvt, int *m, int *n, double *a, int *lda, double *s, double *u, int *ldu,
                 double *vt, int *ldvt, double *work, int *lwork, int *info);

/* lapack 3.0 new SVD routine, faster than xgesvd() */
extern int GESDD(char *jobz, int *m, int *n, double *a, int *lda, double *s, double *u, int *ldu, double *vt, int *ldvt,
                 double *work, int *lwork, int *iwork, int *info);

#else
#define LEVMAR_LUINVERSE dlevmar_LUinverse_noLapack

static int LEVMAR_LUINVERSE(double *A, double *B, int m);
#endif /* HAVE_LAPACK */

/* blocked multiplication of the transpose of the nxm matrix a with itself (i.e. a^T a)
 * using a block size of bsize. The product is returned in b.
 * Since a^T a is symmetric, its computation can be speeded up by computing only its
 * upper triangular part and copying it to the lower part.
 *
 * More details on blocking can be found at 
 * http://www-2.cs.cmu.edu/afs/cs/academic/class/15213-f02/www/R07/section_a/Recitation07-SectionA.pdf
 */

void dtrans_mat_mat_mult(double *a, double *b, int n, int m, int bsize)
{
	int i, j, k, jj, kk;
	double sum, *bim, *akm;

	/* compute upper triangular part using blocking */
	for(jj = 0; jj < m; jj += bsize)
	{
		for(i = 0; i < m; ++i)
		{
			bim = b+i*m;
			for(j = MAX(jj, i); j < MIN(jj+bsize, m); ++j) bim[j]=0.0; //b[i*m+j]=0.0;
		}

		for(kk = 0; kk < n; kk += bsize)
		{
			for(i = 0; i < m; ++i)
			{
				bim = b + i*m;
				for(j=MAX(jj, i); j < MIN(jj+bsize, m); ++j)
				{
					sum=0.0;
					for(k = kk; k < MIN(kk+bsize, n); ++k)
					{
						akm =a + k * m;
						sum += akm[i] * akm[j]; //a[k*m+i]*a[k*m+j];
					}
					
					bim[j] += sum; //b[i*m+j]+=sum;
				}
			}
		}
	}

	/* copy upper triangular part to the lower one */
	for(i = 0; i < m; ++i)
	{
		for(j = 0; j < i; ++j) b[i*m+j]=b[j*m+i];
	}
}

void dfdif_forw_jac_approx(
    const ARM_LEVMARQFunc& func,
													   /* function to differentiate */
    double *p,              /* I: current parameter estimate, mx1 */
    double *hx,             /* I: func evaluated at p, i.e. hx=func(p), nx1 */
    double *hxx,            /* W/O: work array for evaluating func(p+delta), nx1 */
    double delta,           /* increment for computing the jacobian */
    double *jac,            /* O: array for storing approximated jacobian, nxm */
    int m,
    int n,
    void *adata)
{
int i, j;
double tmp;
double d;

  for(j=0; j<m; ++j){
    /* determine d=max(1E-04*|p[j]|, delta), see HZ */
    d=(1E-04)*p[j]; // force evaluation
    d=FABS(d);
    if(d<delta)
      d=delta;

    tmp=p[j];
    p[j]+=d;

    func(p, hxx, m, n, adata);

    p[j]=tmp; /* restore */

    d=(1.0)/d; /* invert so that divisions can be carried out faster as multiplications */
    for(i=0; i<n; ++i){
      jac[i*m+j]=(hxx[i]-hx[i])*d;
    }
  }
}

void dfdif_cent_jac_approx(
    const ARM_LEVMARQFunc& func,
													   /* function to differentiate */
    double *p,              /* I: current parameter estimate, mx1 */
    double *hxm,            /* W/O: work array for evaluating func(p-delta), nx1 */
    double *hxp,            /* W/O: work array for evaluating func(p+delta), nx1 */
    double delta,           /* increment for computing the jacobian */
    double *jac,            /* O: array for storing approximated jacobian, nxm */
    int m,
    int n,
    void *adata)
{
	int i, j;
	double tmp;
	double d;

	for(j = 0; j < m; ++j)
	{
		/* determine d=max(1E-04*|p[j]|, delta), see HZ */
		d = (1E-04) * p[j]; // force evaluation
		d = FABS(d);
		if(d < delta) d = delta;

		tmp = p[j];
		p[j] -= d;
		func(p, hxm, m, n, adata);

		p[j] = tmp + d;
		func(p, hxp, m, n, adata);
		p[j] = tmp; /* restore */

		d = (0.5)/d; /* invert so that divisions can be carried out faster as multiplications */
		for(i = 0; i < n; ++i)
		{
			jac[i*m+j] = (hxp[i] - hxm[i]) * d;
		}
	}
}

/* 
 * Check the jacobian of a n-valued nonlinear function in m variables
 * evaluated at a point p, for consistency with the function itself.
 *
 * Based on fortran77 subroutine CHKDER by
 * Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
 * Argonne National Laboratory. MINPACK project. March 1980.
 *
 *
 * func points to a function from R^m --> R^n: Given a p in R^m it yields hx in R^n
 * jacf points to a function implementing the jacobian of func, whose correctness
 *     is to be tested. Given a p in R^m, jacf computes into the nxm matrix j the
 *     jacobian of func at p. Note that row i of j corresponds to the gradient of
 *     the i-th component of func, evaluated at p.
 * p is an input array of length m containing the point of evaluation.
 * m is the number of variables
 * n is the number of functions
 * adata points to possible additional data and is passed uninterpreted
 *     to func, jacf.
 * err is an array of length n. On output, err contains measures
 *     of correctness of the respective gradients. if there is
 *     no severe loss of significance, then if err[i] is 1.0 the
 *     i-th gradient is correct, while if err[i] is 0.0 the i-th
 *     gradient is incorrect. For values of err between 0.0 and 1.0,
 *     the categorization is less certain. In general, a value of
 *     err[i] greater than 0.5 indicates that the i-th gradient is
 *     probably correct, while a value of err[i] less than 0.5
 *     indicates that the i-th gradient is probably incorrect.
 *
 *
 * The function does not perform reliably if cancellation or
 * rounding errors cause a severe loss of significance in the
 * evaluation of a function. therefore, none of the components
 * of p should be unusually small (in particular, zero) or any
 * other value which may cause loss of significance.
 */

void dlevmar_chkjac(const ARM_LEVMARQFunc& func, const ARM_LEVMARQFunc& funcjac, double *p, int m, int n, void *adata, double *err)
{
	double factor = (100.0);
	double one = (1.0);
	double zero = (0.0);
	double *fvec, *fjac, *pp, *fvecp, *buf;

	int i, j;
	double eps, epsf, temp, epsmch;
	double epslog;
	int fvec_sz = n, fjac_sz = n * m, pp_sz = m, fvecp_sz = n;

	epsmch = DBL_EPSILON;
	eps = (double)sqrt(epsmch);

	buf = (double *)malloc((fvec_sz + fjac_sz + pp_sz + fvecp_sz)*sizeof(double));
//	assert(buf != 0);
	fvec = buf;
	fjac = fvec + fvec_sz;
	pp	 = fjac + fjac_sz;
	fvecp= pp + pp_sz;

	/* compute fvec=func(p) */
	func(p, fvec, m, n, adata);

	/* compute the jacobian at p */
	funcjac(p, fjac, m, n, adata);

	/* compute pp */
	for(j = 0; j < m; ++j)
	{
		temp = eps*FABS(p[j]);
		if(temp==0.) temp = eps;
		pp[j] = p[j] + temp;
	}

	/* compute fvecp=func(pp) */
	func(pp, fvecp, m, n, adata);

	epsf = factor * epsmch;
	epslog = (double)log10(eps);

	for(i = 0; i < n; ++i) err[i] = zero;

	for(j=0; j<m; ++j)
	{
		temp = FABS(p[j]);
		if(temp == zero) temp = one;

		for(i=0; i<n; ++i) err[i] += temp * fjac[i*m+j];
	}

	for(i=0; i<n; ++i)
	{
		temp = one;
		if(fvec[i] != zero && fvecp[i] != zero && FABS(fvecp[i] - fvec[i]) >= epsf * FABS(fvec[i]))
			temp = eps * FABS((fvecp[i] - fvec[i]) / eps - err[i]) / (FABS(fvec[i]) + FABS(fvecp[i]));

		err[i] = one;
		
		if(temp>epsmch && temp<eps)
			err[i] = ((double)log10(temp) - epslog)/epslog;

		if(temp >= eps) err[i]=zero;
	}

	free(buf);

	return;
}

#ifdef HAVE_LAPACK
/*
 * This function computes the pseudoinverse of a square matrix A
 * into B using SVD. A and B can coincide
 * 
 * The function returns 0 in case of error (e.g. A is singular),
 * the rank of A if successfull
 *
 * A, B are mxm
 *
 */
static int LEVMAR_PSEUDOINVERSE(double *A, double *B, int m)
{
	double *buf = NULL;
	int buf_sz=0;
	static double eps = (-1.0);

	int i, j;
	double *a, *u, *s, *vt, *work;
	int a_sz, u_sz, s_sz, vt_sz, tot_sz;
	double thresh, one_over_denom;
	int info, rank, worksz, *iwork, iworksz;
   
	/* calculate required memory size */
	worksz = 16 * m; /* more than needed */
	iworksz = 8 * m;
	a_sz = m * m;
	u_sz = m * m; s_sz = m; vt_sz = m * m;

	tot_sz = iworksz * sizeof(int) + (a_sz + u_sz + s_sz + vt_sz + worksz)*sizeof(double);

    buf_sz = tot_sz;
    buf=(double *)malloc(buf_sz);
	assert(buf != NULL);

	iwork = (int *)buf;
	a=(double *)(iwork+iworksz);
	/* store A (column major!) into a */
	for(i=0; i<m; i++)
	for(j=0; j<m; j++)
		a[i+j*m]=A[i*m+j];

	u=a + a_sz;
	s=u+u_sz;
	vt=s+s_sz;
	work=vt+vt_sz;

	/* SVD decomposition of A */
	GESVD("A", "A", (int *)&m, (int *)&m, a, (int *)&m, s, u, (int *)&m, vt, (int *)&m, work, (int *)&worksz, &info);
  //GESDD("A", (int *)&m, (int *)&m, a, (int *)&m, s, u, (int *)&m, vt, (int *)&m, work, (int *)&worksz, iwork, &info);

  /* error treatment */
	if(info!=0)
	{
		if(info<0)
		{
			//fprintf(stderr, RCAT(RCAT(RCAT("LAPACK error: illegal value for argument %d of ", GESVD), "/" GESDD) " in ", LEVMAR_PSEUDOINVERSE) "()\n", -info);
			//exit(1);
		}
		else
		{
			//fprintf(stderr, RCAT("LAPACK error: dgesdd (dbdsdc)/dgesvd (dbdsqr) failed to converge in ", LEVMAR_PSEUDOINVERSE) "() [info=%d]\n", info);
			free(buf);

			return 0;
		}
	}

	if(eps<0.0)
	{
		double aux;

		/* compute machine epsilon */
		for(eps=(1.0); aux=eps+(1.0), aux-(1.0)>0.0; eps*=(0.5));
		eps*=(2.0);
	}

	/* compute the pseudoinverse in B */
	for(i=0; i<a_sz; i++) B[i]=0.0; /* initialize to zero */
	
	for(rank=0, thresh=eps*s[0]; rank<m && s[rank]>thresh; rank++)
	{
		one_over_denom=(1.0)/s[rank];

		for(j=0; j<m; j++)
			for(i=0; i<m; i++)
				B[i*m+j]+=vt[rank+i*m]*u[j+rank*m]*one_over_denom;
	}

	free(buf);

	return rank;
}
#else // no LAPACK

/*
 * This function computes the inverse of A in B. A and B can coincide
 *
 * The function employs LAPACK-free LU decomposition of A to solve m linear
 * systems A*B_i=I_i, where B_i and I_i are the i-th columns of B and I.
 *
 * A and B are mxm
 *
 * The function returns 0 in case of error,
 * 1 if successfull
 *
 */
static int LEVMAR_LUINVERSE(double *A, double *B, int m)
{
	void *buf=NULL;
	int buf_sz=0;

	int i, j, k, l;
	int *idx, maxi=-1, idx_sz, a_sz, x_sz, work_sz, tot_sz;
	double *a, *x, *work, max, sum, tmp;

	/* calculate required memory size */
	idx_sz=m;
	a_sz=m*m;
	x_sz=m;
	work_sz=m;
	tot_sz=idx_sz*sizeof(int) + (a_sz+x_sz+work_sz)*sizeof(double);

	buf_sz=tot_sz;
	buf=(void *)malloc(tot_sz);
//	assert(buf != NULL);

	idx=(int *)buf;
	a=(double *)(idx + idx_sz);
	x=a + a_sz;
	work=x + x_sz;

	/* avoid destroying A by copying it to a */
	for(i=0; i<a_sz; ++i) a[i]=A[i];

	/* compute the LU decomposition of a row permutation of matrix a; the permutation itself is saved in idx[] */
	for(i=0; i<m; ++i)
	{
		max=0.0;
		for(j=0; j<m; ++j)
			if((tmp=FABS(a[i*m+j]))>max) max=tmp;

		if(max==0.0)
		{
			// matrice singulière
			free(buf);

			return 0;
		}

		work[i]=(1.0)/max;
	}

	for(j=0; j<m; ++j)
	{
		for(i=0; i<j; ++i)
		{
			sum=a[i*m+j];
			for(k=0; k<i; ++k)
        sum-=a[i*m+k]*a[k*m+j];
			a[i*m+j]=sum;
		}
		max=0.0;
		for(i=j; i<m; ++i){
			sum=a[i*m+j];
			for(k=0; k<j; ++k)
        sum-=a[i*m+k]*a[k*m+j];
			a[i*m+j]=sum;
			if((tmp=work[i]*FABS(sum))>=max){
				max=tmp;
				maxi=i;
			}
		}
		if(j!=maxi){
			for(k=0; k<m; ++k){
				tmp=a[maxi*m+k];
				a[maxi*m+k]=a[j*m+k];
				a[j*m+k]=tmp;
			}
			work[maxi]=work[j];
		}
		idx[j]=maxi;
		if(a[j*m+j]==0.0)
      a[j*m+j]=DBL_EPSILON;
		if(j!=m-1){
			tmp=(1.0)/(a[j*m+j]);
			for(i=j+1; i<m; ++i)
        a[i*m+j]*=tmp;
		}
	}

  /* The decomposition has now replaced a. Solve the m linear systems using
   * forward and back substitution
   */
  for(l=0; l<m; ++l){
    for(i=0; i<m; ++i) x[i]=0.0;
    x[l]=(1.0);

	  for(i=k=0; i<m; ++i){
		  j=idx[i];
		  sum=x[j];
		  x[j]=x[i];
		  if(k!=0)
			  for(j=k-1; j<i; ++j)
          sum-=a[i*m+j]*x[j];
		  else
        if(sum!=0.0)
			    k=i+1;
		  x[i]=sum;
	  }

	  for(i=m-1; i>=0; --i){
		  sum=x[i];
		  for(j=i+1; j<m; ++j)
        sum-=a[i*m+j]*x[j];
		  x[i]=sum/a[i*m+i];
	  }

    for(i=0; i<m; ++i)
      B[i*m+l]=x[i];
  }

  free(buf);

  return 1;
}
#endif /* HAVE_LAPACK */

/*
 This function computes in C the covariance matrix corresponding to a least
 squares fit. JtJ is the approximate Hessian at the solution (i.e. J^T*J, where
 J is the jacobian at the solution), sumsq is the sum of squared residuals
 (i.e. goodnes of fit) at the solution, m is the number of parameters (variables)
 and n the number of observations. JtJ can coincide with C.
 
 if JtJ is of full rank, C is computed as sumsq/(n-m)*(JtJ)^-1
 otherwise and if LAPACK is available, C=sumsq/(n-r)*(JtJ)^+
 where r is JtJ's rank and ^+ denotes the pseudoinverse
 The diagonal of C is made up from the estimates of the variances
 of the estimated regression coefficients.
 See the documentation of routine E04YCF from the NAG fortran lib

 The function returns the rank of JtJ if successful, 0 on error

 A and C are mxm

 */
int dlevmar_covar(double *JtJ, double *C, double sumsq, int m, int n)
{
int i;
int rnk;
double fact;

#ifdef HAVE_LAPACK
   rnk=LEVMAR_PSEUDOINVERSE(JtJ, C, m);
   if(!rnk) return 0;
#else
#ifdef _MSC_VER
#pragma message("LAPACK not available, LU will be used for matrix inversion when computing the covariance; this might be unstable at times")
#else
#warning LAPACK not available, LU will be used for matrix inversion when computing the covariance; this might be unstable at times
#endif // _MSC_VER

   rnk=LEVMAR_LUINVERSE(JtJ, C, m);
   if(!rnk) return 0;

   rnk=m; /* assume full rank */
#endif /* HAVE_LAPACK */

   fact=sumsq/(double)(n-rnk);
   for(i=0; i<m*m; ++i)
     C[i]*=fact;

   return rnk;
}

/* undefine everything. THIS MUST REMAIN AT THE END OF THE FILE */
#undef LEVMAR_PSEUDOINVERSE
#undef LEVMAR_LUINVERSE
#undef LEVMAR_COVAR
#undef LEVMAR_CHKJAC
#undef FDIF_FORW_JAC_APPROX
#undef FDIF_CENT_JAC_APPROX
#undef TRANS_MAT_MAT_MULT


CC_END_NAMESPACE()
