/* ==========================================================================
   FILE_NAME:	num_f_qrupdt.c

   PURPOSE:         
   ========================================================================== */
#define NRANSI
#include "utallhdr.h"
#include "math.h"

void qrupdt(double **r, double **qt, int n, double u[], double v[])
{
	void rotate(double **r, double **qt, int n, int i, double a, double b);
	int i,j,k;

	for (k=n;k>=1;k--) {
		if (u[k]) break;
	}
	if (k < 1) k=1;
	for (i=k-1;i>=1;i--) {
		rotate(r,qt,n,i,u[i],-u[i+1]);
		if (u[i] == 0.0) u[i]=fabs(u[i+1]);
		else if (fabs(u[i]) > fabs(u[i+1]))
			u[i]=fabs(u[i])*sqrt(1.0+SQR(u[i+1]/u[i]));
		else u[i]=fabs(u[i+1])*sqrt(1.0+SQR(u[i]/u[i+1]));
	}
	for (j=1;j<=n;j++) r[1][j] += u[1]*v[j];
	for (i=1;i<k;i++)
		rotate(r,qt,n,i,r[i][i],-r[i+1][i]);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software +135[)6=. */
