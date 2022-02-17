/* ============================================================

  FILENAME:	    num_f_fdjac.c

  
  PURPOSE:	    diagonalise_symmetric_matrix
  ============================================================= */
#define NRANSI
#include "utallhdr.h"
#include "math.h"
#define EPS_FDJAC 1.0e-4

void fdjac(int n, double x[], double fvec[], double **df,
	void (*vecfunc)(int, double [], double []))
{
	int i,j;
	double h,temp,*f;

	f=dvector(1,n);
	for (j=1;j<=n;j++) {
		temp=x[j];
		h=EPS_FDJAC*fabs(temp);
		if (h == 0.0) h=EPS_FDJAC;
		x[j]=temp+h;
		h=x[j]-temp;
		(*vecfunc)(n,x,f);
		x[j]=temp;
		for (i=1;i<=n;i++) df[i][j]=(f[i]-fvec[i])/h;
	}
	free_dvector(f,1,n);
}
#undef EPS_FDJAC
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software +135[)6=. */
