/****************************************************************************/
/*      Utility routines for Newton-Raphson                                 */
/****************************************************************************/
/*      UTIL_NR.c                                                           */
/****************************************************************************/


/*
$Header: /nasdev/export2/home/drdev/cvsadmin/cvs/libs/fix3/src/util.c,v 1.22 2005/02/04 20:01:49 skuzniar Exp $
*/


#include <stdio.h>
#include <stdlib.h>  
#include <math.h>
#include <string.h>
#include "fix123head.h"

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

/*******************************************************************/
int  gaussjORI(double **a, int n, double **b, int m)
{
    int  status = FAILURE;

	int *indxc = NULL,
        *indxr = NULL,
        *ipiv  = NULL;

	int    i,icol=0,irow=0,j,k,l,ll;
	double big,dum,pivinv,temp;

	indxc = DR_Array(INT, 1, n);
	indxr = DR_Array(INT, 1, n);
	ipiv  = DR_Array(INT, 1, n);
    if ((indxc == NULL) || (indxr == NULL) || (ipiv == NULL)) goto RETURN;

	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} 
                    else if (ipiv[k] > 1)
                    {
                        DR_Error("gaussj: Singular Matrix-1");
                        goto RETURN;
                    }
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
			for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0)
        {
            DR_Error("gaussj: Singular Matrix-2");
            goto RETURN;
        }
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *= pivinv;
		for (l=1;l<=m;l++) b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n;l>=1;l--) {
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}

    status = SUCCESS;

RETURN:

    Free_DR_Array(ipiv,  INT, 1, n);
    Free_DR_Array(indxr, INT, 1, n);
    Free_DR_Array(indxc, INT, 1, n);

    return(status);
}


/*****  MatrixInverse *******************************************************/
/*
 *   Finds the inverse of A and return it in C.
 *   On entry, memory must be allocated for C.
 *   This routine allows A and C to be the same variable (i.e. share the 
 *   same memory space)
 *
 */
int   MatrixInverse(long      n,     /* (I)   dimensions                 */
                    double  **A,     /* (I)   A[0..n-1][0..n-1]          */ 
                    double  **C)     /* (O)   inv (mem alloc'd on input) */ 
{
    int       status = FAILURE;

    double  **ALocal = NULL;
    double  **b      = NULL;

    int       i, j;

    if (n < 1L) goto RETURN;

    if (n == 1L)
    {
        C[0][0] = 1.0/A[0][0];
        return (SUCCESS);
    }

    ALocal = (double **) DR_Matrix(DOUBLE, 1, n, 1, n);
    b      = (double **) DR_Matrix(DOUBLE, 1, n, 1, 1);
    if ((ALocal == NULL) || (b == NULL)) goto RETURN;

    for (i=1; i<=n; i++)
    {
        b[i][1] = 0.0;
        for (j=1; j<=n; j++)
        {
            ALocal[i][j] = A[i-1][j-1];
        }
    }

    if (gaussjORI(ALocal, n, b, 1) == FAILURE) goto RETURN;

    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
        {
            C[i][j] = ALocal[i+1][j+1];
        }
    }

    status = SUCCESS;

RETURN:

    Free_DR_Matrix(ALocal, DOUBLE, 1, n, 1, n);
    Free_DR_Matrix(b,      DOUBLE, 1, n, 1, n);

    return(status);

} /* MatrixInverse */

