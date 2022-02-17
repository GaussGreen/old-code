/************************************************************************
 * Module:	DCU - LINQ
 * Function:	Eigen stuff
 * Author:	C. Daher
 * REMARK:	<<< call to NAG Fortran library >>
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */
#include "drlerr.h"		/* DrlErrMsg */

#include <math.h>

#include "drlmem.h"

#include "drllineq.h"		/* prototype consistency */

#if defined(USEIMSL)
# include "imsl.h"
#endif

#if defined(USENAG)
extern	void	c_f04axf(int *, ...) ;
#endif


static	int	jacobi(double **a, int n, double d[], double **v, int *nrot);



/*f---------------------------------------------------------------------
 * Matrix eigenvalues and vector.
 * 
 * <br><br>
 * This routines computes the eigen values and vectors
 * of the square matrix <i> a</i> of size n.
 * The <i>n</i> eigen values are placed in the vector <i> vap</i> and
 * the eigen vectors in <i> vep</i>  (<i> vep[i]</i> is the vector 
 * for eigen value <i> vap[i]</i>).
 * Both <i> vap</i> and <i> vep</i> should be allocated <i> before</i> the 
 * routine call.
 * The routine returns a failure if one of the eigenvalues
 * has a non-zero imaginary part.
 * Returns 0 on success.<br>
 * <b> Remark: Must be linked with IMSL</b>.
 */

int
DrlMatrixRealEigenVect(
	int n,		/* (I) dimension */
	double **a,	/* (I) matrix */
	double *vap,	/* (O) array of eigen values */
	double **vep)	/* (O) array of eigen vectors */
{
static	char	routine[] = "DrlMatrixRealEigenVect";

#if defined(USENAG)
	int	i, j, ifail;
	double	*vfar,	/* contains Re a[i][j] */
		*vfai,	/* Im a[i][j] */
		*vfvr,	/* Re v_i[j] (i-th eigen vec) */
		*vfvi ;	/* Im v_i[j] (i-th eigen vec) */

	ifail = 0 ;
	vfar = NULL;
	vfai = NULL;
	vfvr = NULL;
	vfvi = NULL;

	/*
	 * Allocate memory
	 */
	if ( ((vfar = DoubleVectAlloc(0, n*n-1)) == NULL)
	  || ((vfai = DoubleVectAlloc(0, n*n-1)) == NULL)
	  || ((vfvr = DoubleVectAlloc(0, n*n-1)) == NULL)
	  || ((vfvi = DoubleVectAlloc(0, n*n-1)) == NULL)) {
		ifail = -4 ;
		goto done;
	}


	/*
	 * transfer matrix
	 */
	for (i=0; i<=n-1; i++)
	for (j=0; j<=n-1; j++) {
		vfar[i*n+j] = a[i][j] ; /* Re a[i][j] */
		vfai[i*n+j] = 0.0 ;	/* Im a[i][j] */
	}


	/*
	 * Appel routine fortran
	 */
	c_f04axf(&n, vfar, vfai, vap, vfvr, vfvi, &ifail) ;


	/*
	 * Transfer vectors
	 */
	for (i=0; i<=n-1; i++)
	for (j=0; j<=n-1; j++) {
		vep[i][j] = vfvr[i*n+j] ;
	}


	/*
	 * FREE memory
	 */

done:
	if (vfar != NULL) DoubleVectFree(0, vfar-1);
	if (vfai != NULL) DoubleVectFree(0, vfai-1);
	if (vfvr != NULL) DoubleVectFree(0, vfvr-1);
	if (vfvi != NULL) DoubleVectFree(0, vfvi-1);

	if (ifail != 0) {
		DrlErrMsg("%s: Failed\n", routine);
	}

	return(ifail);
#elif defined(USEIMSL)
	int		status = FAILURE;
	int		i, j;
	double		*avec = NULL;
	d_complex	*evalu,
			*evecu,
			*eval;

	if ((avec = DrlDoubleVectAlloc(0, n*n-1)) == NULL)
		goto done;
	for (i=0; i<=n-1; i++)
	for (j=0; j<=n-1; j++)
		avec[n*i+j] = a[i][j];
		
	if ((evalu = NEW_ARRAY(d_complex, n)) == NULL)
		goto done;
	if ((evecu = NEW_ARRAY(d_complex, n*n)) == NULL)
		goto done;


	eval = imsl_d_eig_gen(n, avec,
		IMSL_VECTORS_USER, evecu,
		IMSL_RETURN_USER, evalu,
		0);
	if (eval == NULL) {
		DrlErrMsg("%s: IMSL failure.\n", routine);
		goto done;
	}


	for (i=0; i<=n-1; i++) {
	    if (!IS_ALMOST_ZERO(evalu[i].im)) {
		DrlErrMsg("%s: eigen val #%d (%lf,%lf) has im part.\n",
			routine, i,
	    		evalu[i].re,
	    		evalu[i].im);
		goto done;
	    }
	    vap[i] = evalu[i].re;
	}


	for (i=0; i<=n-1; i++)
	for (j=0; j<=n-1; j++) {
		vep[i][j] = evecu[n*i+j].re;
	}


	status = SUCCESS;
done:
	if (evalu != NULL) FREE(evalu);
	if (evecu != NULL) FREE(evecu);
	DrlDoubleVectFree(avec, 0, n*n-1);
	if (status != SUCCESS) {
		DrlErrMsg("%s: failed.\n", routine);
	}
	return(status);
#else

	DrlErrMsg("%s: not implemented\n", routine);
	return(-1);
#endif
}



/*f---------------------------------------------------------------------
 * Matrix eigenvalues and vector for symmetric matrix using Jacobi..
 * 
 * <br><br>
 * This routines computes the eigen values and vectors
 * of the square matrix <i> a</i> of size n.
 * The <i>n</i> eigen values are placed in the vector <i> vap</i> and
 * the eigen vectors in <i> vep</i>  (<i> vep[i]</i> is the vector 
 * for eigen value <i> vap[i]</i>).
 * Both <i> vap</i> and <i> vep</i> should be allocated <i> before</i> the 
 * routine call.
 * The routine returns a failure if the matrix os not symmetric
 * or convergence is not achieved.
 * Returns 0 on success.<br>
 */

int
DrlMatrixRealSymEigen(
	int n,		/* (I) dimension */
	double **a,	/* (I) matrix */
	double *vap,	/* (O) array of eigen values */
	double **vep)	/* (O) array of eigen vectors */
{
static	char	routine[] = "DrlMatrixRealSymEigen";
	int	status = FAILURE;
	int	i, j;
	double	**a1 = NULL;
	int	nrot;


	/* Check symmetric */
	for (i=0; i<n; i++)
	for (j=i+1; j<n; j++) {
	    if (!IS_ALMOST_ZERO(a[i][j] - a[j][i])) {
		DrlErrMsg("%s: a[%d][%d] (%lf) != a[%d][%d] (%lf).\n",
			routine, i, j, a[i][j], j, i, a[j][i]);
		DrlErrMsg("%s: matrix not symmetric.\n", routine);
		goto done;
	    }
	}


	/*
	 * Use static numerical recipe routine
	 */

	if (DrlMatrixCopy(&a1, a, n, n) != SUCCESS)
		goto done;
	DrlDoubleMatrCToNr(&a1, 1, n, 1, n);


	DrlDoubleMatrCToNr(&vep, 1, n, 1, n);
	vap -= 1;

	IF_FAILED_DONE ( jacobi(
		a1,
		n,
		vap,
		vep,
		&nrot));


	DrlDoubleMatrNrToC(&a1, 1, n, 1, n);
	vap += 1;
	DrlDoubleMatrNrToC(&vep, 1, n, 1, n);

	/*
	 * Reorder the eigenvalues
	 */
	{
		double	a,
			*b = NULL;
		int	k;

		ASSERT_OR_DONE((b = DrlDoubleVectAlloc(0, n-1)) != NULL);

		for (j=1;j<=n-1;j++) {
			a = vap[j];
			for (k=0; k<n; k++) 
				b[k] = vep[k][j];
			i=j-1;
			while (i >= 0 && vap[i] < a) {
				vap[i+1] = vap[i];
				for (k=0; k<n; k++) 
					vep[k][i+1] = vep[k][i];
				i--;
			}
			vap[i+1] = a;
			for (k=0; k<n; k++) 
				vep[k][i+1] = b[k];
		}
		DrlDoubleVectFree(b, 0, n-1);
	}


	status = SUCCESS;
	/* made it through OK */
done:
	DrlMatrixFree(a1, n, n);
	if (status != SUCCESS) {
		DrlErrMsg("%s: failed\n", routine);
	}
	return(status);
}


/*----------------------------------------------------------------------
 * Support numerical recipes routines. See use guide.
 *
 */


#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);

static	int	jacobi(
	double **a,
	int n,
	double d[],
	double **v,
	int *nrot)
{
static	char	routine[] = "jacobi";
    int j,iq,ip,i;
    double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

    b=DrlDoubleVectAlloc(1,n);
    ASSERT_OR_DONE(b != NULL);
    z=DrlDoubleVectAlloc(1,n);
    ASSERT_OR_DONE(z != NULL);

    for (ip=1;ip<=n;ip++) {
        for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
        v[ip][ip]=1.0;
    }
    for (ip=1;ip<=n;ip++) {
        b[ip]=d[ip]=a[ip][ip];
        z[ip]=0.0;
    }
    *nrot=0;
    for (i=1;i<=50;i++) {
        sm=0.0;
        for (ip=1;ip<=n-1;ip++) {
            for (iq=ip+1;iq<=n;iq++)
                sm += fabs(a[ip][iq]);
        }
        if (sm == 0.0) {
            DrlDoubleVectFree(z,1,n);
            DrlDoubleVectFree(b,1,n);
            return(0);
        }
        if (i < 4)
            tresh=0.2*sm/(n*n);
        else
            tresh=0.0;
        for (ip=1;ip<=n-1;ip++) {
            for (iq=ip+1;iq<=n;iq++) {
                g=100.0*fabs(a[ip][iq]);
                if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
                    && (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
                    a[ip][iq]=0.0;
                else if (fabs(a[ip][iq]) > tresh) {
                    h=d[iq]-d[ip];
                    if ((double)(fabs(h)+g) == (double)fabs(h))
                        t=(a[ip][iq])/h;
                    else {
                        theta=0.5*h/(a[ip][iq]);
                        t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                        if (theta < 0.0) t = -t;
                    }
                    c=1.0/sqrt(1+t*t);
                    s=t*c;
                    tau=s/(1.0+c);
                    h=t*a[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    a[ip][iq]=0.0;
                    for (j=1;j<=ip-1;j++) {
                        ROTATE(a,j,ip,j,iq)
                    }
                    for (j=ip+1;j<=iq-1;j++) {
                        ROTATE(a,ip,j,j,iq)
                    }
                    for (j=iq+1;j<=n;j++) {
                        ROTATE(a,ip,j,iq,j)
                    }
                    for (j=1;j<=n;j++) {
                        ROTATE(v,j,ip,j,iq)
                    }
                    ++(*nrot);
                }
            }
        }
        for (ip=1;ip<=n;ip++) {
            b[ip] += z[ip];
            d[ip]=b[ip];
            z[ip]=0.0;
        }
    }

    DrlErrMsg("jacobi: too many iterations.\n");
done:
    DrlErrMsg("jacobi: failed.\n");
    return(-1);
}
#undef ROTATE











