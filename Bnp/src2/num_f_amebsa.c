/* ==========================================================
   FILENAME :        num_f_amebsa.c

   PURPOSE:          for Simplex
   ========================================================== */

#define NRANSI
#include "utallhdr.h"
#include "math.h"

#define GET_PSUM \
					for (n=1;n<=ndim;n++) {\
					for (sum=0.0,m=1;m<=mpts;m++) sum += p[m][n];\
					psum[n]=sum;\
}


extern long sa_seed;    /* Defined and initialized in main */


/* ------------------------------------------------------------------------ 
	Extrapolates by a factor fac through the face of simplex across from 
	the high point, tries it, and replaces the high point if new point is 
	better.
   ------------------------------------------------------------------------ */

static double amotsa(
		double **p,                  /* Simplex points */
		double y[],                  /* Value of funk at the simplex points */
		double psum[],               
		int ndim, 
		double pb[],
		double *yb, 
		double (*funk)(double []),
		int ihi,
		double *yhi, 
		double fac,
		double temptr)
{
	double uniform_fast(long *seed);
	int j;
	double fac1,fac2,yflu,ytry,*ptry;

	ptry = dvector(1,ndim);
	fac1 = (1.0-fac)/ndim;
	fac2 = fac1-fac;
	for (j=1; j<=ndim; j++)
	{
		ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
	}
	ytry=(*funk)(ptry);
	if (ytry <= *yb)
	{
		for (j=1;j<=ndim;j++) 
			pb[j]=ptry[j];
		*yb=ytry;
	}
/* We added a thermal fluctuation to all the current points,
	but we substract it here: we want the algorithm to accept a priori 
	any suggested change */	
	yflu = ytry  - temptr * (-log(uniform_fast(&sa_seed)));
	if (yflu < *yhi)
	{
		y[ihi]=ytry;
		*yhi=yflu;
		for (j=1;j<=ndim;j++)
		{
			psum[j] += ptry[j]-p[ihi][j];
			p[ihi][j]=ptry[j];
		}
	}
	free_dvector(ptry,1,ndim);
	return yflu;
}

/* ------------------------------------------------------------------------
                            SIMULATED ANNEALING 
   ------------------------------------------------------------------------ */

void amebsa(
		double **p,                   /* Starting points of the simplex */
		double y[],                   /* Values at the starting points */
		int ndim,                     /* Number of dimension */
		double pb[],                  /* Best point ever encountered */
		double *yb,                   /* Best function value encountered */
		double ftol,                  /* Fractional convergence tolerance */
		double (*funk)(double []), 
		int *iter,                    /* Number of function evaluations */
		double temptr)                /* Temperature */
{
	double amotsa(double **p, double y[], double psum[], int ndim, double pb[],
		double *yb, double (*funk)(double []), int ihi, double *yhi, 
		double fac, double temptr);
	double uniform_fast(long *seed);
	int i,ihi,ilo,j,m,n,mpts=ndim+1;
	double rtol,sum,swap,yhi,ylo,ynhi,ysave,yt,ytry,*psum;

	psum = dvector(1,ndim);
	GET_PSUM
	for (;;)
	{
	/* Determine which point is the worst, second worst and best */
		ilo=1;
		ihi=2;
	/* Each point of the simplex gets a positive random thermal fluctuation */
		ynhi=ylo=y[1] + temptr * (-log(uniform_fast(&sa_seed)));
		yhi=y[2] + temptr * (-log(uniform_fast(&sa_seed)));
		if (ylo > yhi) 
		{
			ihi=1;
			ilo=2;
			ynhi=yhi;
			yhi=ylo;
			ylo=ynhi;
		}
	/* Loop over the points in the simplex */		
		for (i=3;i<=mpts;i++)
		{
		/* More positive random thermal fluctuations */
			yt=y[i] + temptr * (-log(uniform_fast(&sa_seed)));
			if (yt <= ylo) 
			{
				ilo=i;
				ylo=yt;
			}
			if (yt > yhi)
			{
				ynhi=yhi;
				ihi=i;
				yhi=yt;
			} 
			else 
			if (yt > ynhi) 
			{
				ynhi=yt;
			}
		}
	/* Compute the fractional range from worst to best */
		rtol=2.0*fabs(yhi-ylo)/(fabs(yhi)+fabs(ylo));
		if (rtol < ftol || *iter < 0) 
		{
		/* If criteria satisfactory put best point and value in slot 1 */
			swap=y[1];
			y[1]=y[ilo];
			y[ilo]=swap;
			for (n=1;n<=ndim;n++)
			{
				swap=p[1][n];
				p[1][n]=p[ilo][n];
				p[ilo][n]=swap;
			}
			break;
		}
	/* Begins a new iteration: reflects the simplex from the high point */
		*iter -= 2;
		ytry = amotsa(p,y,psum,ndim,pb,yb,funk,ihi,&yhi,-1.0, temptr);
		if (ytry <= ylo)
		{
		/* Gives a result better than the best point: extrapolation by 2*/
			ytry = amotsa(p,y,psum,ndim,pb,yb,funk,ihi,&yhi,2.0, temptr);
		}
		else 
		if (ytry >= ynhi)
		{
		/* Reflected point worse than second worse: one dimensional contract. */
			ysave=yhi;
			ytry = amotsa(p,y,psum,ndim,pb,yb,funk,ihi,&yhi,0.5, temptr);
			if (ytry >= ysave)
			{
			/* Cannot get rid of tha high point: contract around best point */
				for (i=1;i<=mpts;i++)
				{
					if (i != ilo) 
					{
						for (j=1;j<=ndim;j++) 
						{
							psum[j]=0.5*(p[i][j]+p[ilo][j]);
							p[i][j]=psum[j];
						}
						y[i]=(*funk)(psum);
					}
				}
				*iter -= ndim;
			/* Recompute psum */
				GET_PSUM
			}
		}
		else 
	/* Correct the evaluation count */
		++(*iter);
	}
	free_dvector(psum,1,ndim);
}
#undef GET_PSUM
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software +135[)6=. */
