/* ==========================================================
   FILENAME :        num_f_amoeba.c
   ========================================================== */
#define NRANSI
#include "utallhdr.h"
#include "math.h"
/*#define NMAX 5000*/
#define GET_PSUM \
				for (j=1;j<=ndim;j++) {\
				for (sum=0.0,i=1;i<=mpts;i++) sum += p[i][j];\
				psum[j]=sum;}

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

static double amotry(double **p, double y[], double psum[], int ndim,
	double (*funk)(double []), int ihi, double fac)
{
	int j;
	double fac1,fac2,ytry,*ptry;

	ptry=vector(1,ndim);
	fac1=(1.0-fac)/ndim;
	fac2=fac1-fac;
	for (j=1;j<=ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
	ytry=(*funk)(ptry);
	if (ytry < y[ihi]) {
		y[ihi]=ytry;
		for (j=1;j<=ndim;j++) {
			psum[j] += ptry[j]-p[ihi][j];
			p[ihi][j]=ptry[j];
		}
	}
	free_vector(ptry,1,ndim);
	return ytry;
}


Err amoeba(
		double **p, 
		double y[], 
		int ndim, 
		double ftol,
		int maxiter, 
		double (*funk)(double []), 
		int *nfunk)
{
	int i,ihi,ilo,inhi,j,mpts=ndim+1;
	double rtol,sum,swap,ysave,ytry,*psum;

	psum = dvector(1,ndim);
	*nfunk=0;
	GET_PSUM
	for (;;)
	{
		ilo=1;
	/* Determine wich points are the best, worst, second worst */		
		ihi = y[1]>y[2] ? (inhi=2,1) : (inhi=1,2);
		for (i=1;i<=mpts;i++)
		{
			if (y[i] <= y[ilo]) ilo=i;
			if (y[i] > y[ihi])
			{
				inhi=ihi;
				ihi=i;
			} 
			else 
			if (y[i] > y[inhi] && i != ihi) 
				inhi=i;
		}
		rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]));
	/* Compute fractionnal range from highest to lowest and return if good */
		if ((rtol < ftol) || (*nfunk >= maxiter))
		{
			SWAP(y[1],y[ilo])
			for (i=1;i<=ndim;i++) SWAP(p[1][i],p[ilo][i])
			break;
		}
		*nfunk += 2;
	/* Begin a new iteration */
		ytry=amotry(p,y,psum,ndim,funk,ihi,-1.0);
		if (ytry <= y[ilo])
		{
		/* Gives a result better than the best point */
			ytry = amotry(p,y,psum,ndim,funk,ihi,2.0);
		}
		else
		if (ytry >= y[inhi])
		{
		/* Worst than second highest: do a one-dimensionnal contraction */
			ysave=y[ihi];
			ytry=amotry(p,y,psum,ndim,funk,ihi,0.5);
			if (ytry >= ysave)
			{
			/* Cannot get rid of the high point: see around best point */
				for (i=1;i<=mpts;i++)
				{
					if (i != ilo)
					{
						for (j=1;j<=ndim;j++)
							p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
						y[i]=(*funk)(psum);
					}
				}
				*nfunk += ndim;
				GET_PSUM
			}
		} 
		else --(*nfunk);
	}
	free_dvector(psum,1,ndim);

/* Return a success message */
	return NULL;
}

#undef SWAP
#undef GET_PSUM
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software +135[)6=. */
