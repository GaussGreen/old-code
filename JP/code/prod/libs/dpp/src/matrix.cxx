/****************************************************************
 * Module:	PenGuin
 * Submodule:	
 * File:	modpar.c
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 ****************************************************************/
#include "kmatrix.h"
#include "kutilios.h"

extern	"C" {
#include "drlio.h"
#include "drlstr.h"
#include "drllineq.h"
#include "drlmem.h"
};

#if !defined (__SUNPRO_CC) || (__SUNPRO_CC != 0x0500) 
//===============================================================
//
//
//===============================================================


static	double**
mNew(const KDMatrix& mat)
{
	int	i, j;
	double	**a = NULL;
	ASSERT_OR_THROW((a = DrlMatrixNew(mat.rows(), mat.cols())) != NULL);
	for (i=0; i<mat.rows(); i++)
	for (j=0; j<mat.cols(); j++) {
		a[i][j] = mat[i][j];
	}
	return (a);
}

static	void
mFree(double **a, int nRows, int nCols, KDMatrix *mat)
{
	int	i, j;
	if (mat) {
		mat->resize(nRows, nCols);
		for (i=0; i<mat->rows(); i++)
		for (j=0; j<mat->cols(); j++) {
			(*mat)[i][j] = a[i][j];
		}
	}

	DrlMatrixFree(a, nRows, nCols);
}



//---------------------------------------------------------------
//

KDMatrix::KDMatrix(int nRows, KDMatrixType type)
	: KMatrix<double>(nRows, nRows, NULL)
{
	int	i, j;
	for (i=0; i<rows(); i++)
	for (j=0; j<cols(); j++)
		(*this)[i][j] = (i == j ? 1e0 : 0e0);

}

//---------------------------------------------------------------
//


istream&
operator>>(istream& is, KDMatrix& mat)
{
static	char	routine[] = "operator>>(istream&, KDMatrix&)";
	int	nRows, nCols, i, j;
    try {
	nRows = getDouble(is, "nRows");
	nCols = getDouble(is, "nCols");

	mat.resize(nRows, nCols);

	for (i=0; i<nRows; i++)
	for (j=0; j<nCols; j++) {
		mat[i][j] = getDouble(is, "matrix element");
	}
	return(is);
    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}


//---------------------------------------------------------------

KDVector operator*(const KDMatrix &aMat, const KDVector &bVect)
{
static	char	routine[] = "operator*(const KDMatrix&, const KDVector &)";
	int	i, j;
	KDVector retVal;


	if (aMat.cols() != bVect.size()) {
		throw KFailure("%s: can't multiply: a-cols %d, b-rows %d.\n",
			routine, aMat.cols(), bVect.size());
	}
	retVal.resize(aMat.rows());

	for(i=0; i<aMat.rows(); i++) {
		retVal[i] = 0e0;
		for(j=0; j<bVect.size(); j++) {
			retVal[i] += aMat[i][j]*bVect[j];
		}
	}
	return(retVal);
}



//---------------------------------------------------------------
//

KDMatrix
KDMatrix::operator*(const KDMatrix &bMat) const
{
static	char	routine[] = "KDMatrix::operator*";
	int	i, j, k;
	double	x;
	const KDMatrix &aMat = *this;
	KDMatrix	prodMat;


	if (aMat.cols() != bMat.rows()) {
		throw KFailure("%s: can't multiply: a-cols %d, b-rows %d.\n",
			routine, aMat.cols(), bMat.rows());
	}


	prodMat.resize(aMat.rows(), bMat.cols());

	for(i=0; i<aMat.rows(); i++)
	for(j=0; j<bMat.cols(); j++) {
		x = 0e0;
		for(k=0; k<aMat.cols(); k++) {
			x += aMat[i][k]*bMat[k][j];
		}
		prodMat[i][j] = x;
	}

	return(prodMat);
}








//---------------------------------------------------------------
//

KDMatrix
KDMatrix::inverse() const
{
static	char	routine[] = "KDMatrix::inverse";

	KDMatrix invMat;
	double	**amat = NULL,
		**bmat = NULL;
	int	n=rows();

    try {

	if (rows() != cols()) {
		throw KFailure("%s: rows %d != cols %d.\n", routine,
			rows(), cols());
	}
	amat = mNew(*this);

	IF_FAILED_THROW( DrlMatrixInvert(
		&bmat, 
		amat,
		n));

	mFree(amat, n, n, NULL);
	mFree(bmat, n, n, &invMat);

	return(invMat);

    }
    catch (KFailure) {
	mFree(amat, n, n, NULL);
	mFree(bmat, n, n, &invMat);
	throw KFailure("%s: failed.\n", routine);
    }

}



//---------------------------------------------------------------
// Computes eigen values and vectors.


void
KDMatrix::eigen(KDVector &eigVal, KDMatrix &eigVec) const
{
static	char	routine[] = "KDMatrix::EigenValues";


	double	*vap = NULL,
		**vep = NULL,
		**am= NULL;
	int	n, i, j;

    try {

	if (cols() != rows()) {
		throw KFailure("%s: non diagonal matrix (%d x %d).\n",
			routine, rows(), cols());
	}
	n = rows();


	ASSERT_OR_THROW((vap = DrlDoubleVectAlloc(0, n-1)) != NULL);
	ASSERT_OR_THROW((vep = DrlDoubleMatrAlloc(0, n-1, 0, n-1)) != NULL);
	ASSERT_OR_THROW((am  = DrlDoubleMatrAlloc(0, n-1, 0, n-1)) != NULL);

	for(i=0; i<=n-1; i++) 
	for(j=0; j<=n-1; j++) 
		am[i][j] = (*this)[i][j];

	IF_FAILED_THROW( DrlMatrixRealEigenVect(
		n,
		am,
		vap,
		vep));

	//
	eigVal.resize(n);
	for(i=0; i<=n-1; i++) 
		eigVal[i] = vap[i];
	eigVec.resize(n, n);
	for(i=0; i<=n-1; i++)
	for(j=0; j<=n-1; j++)
		eigVec[i][j] = vep[i][j];


	DrlDoubleVectFree(vap, 0, n-1);
	DrlDoubleMatrFree(vep, 0, n-1, 0, n-1);
	DrlDoubleMatrFree(am, 0, n-1, 0, n-1);

    }
    catch (KFailure) {
	DrlDoubleVectFree(vap, 0, n-1);
	DrlDoubleMatrFree(vep, 0, n-1, 0, n-1);
	DrlDoubleMatrFree(am, 0, n-1, 0, n-1);
	throw KFailure("%s: failed.\n", routine);
    }
}


//---------------------------------------------------------------
//


KDMatrix
directsum(const KDMatrix& aMat, const KDMatrix &bMat)
{
static	char	routine[] = "KDMatrix::directsum";
	int	i, j, nar, nac;
	KDMatrix	prodMat;

	prodMat.resize(
		aMat.rows() + bMat.rows(),
		aMat.cols() + bMat.cols());

	prodMat = 0e0;

	for(i=0; i<aMat.rows(); i++)
	for(j=0; j<aMat.cols(); j++) {
		prodMat[i][j] = aMat[i][j];
	}

	nar = aMat.rows();
	nac = aMat.cols();
	for(i=0; i<bMat.rows(); i++)
	for(j=0; j<bMat.cols(); j++) {
		prodMat[i+nar][j+nac] = aMat[i][j];
	}


	return(prodMat);
}

//---------------------------------------------------------------
// Computes the pseudo inverse of a matrix
// obtained from an SVD. The arguments "regType", "alpha"
// and "param" determine the type and amount
// of regularization applied. Assume
// (w_0, ..., w_N) are the eigen values of the input matrix, then
// If "regType" is 0, no regulariztion is performed, i.e.
// lambda_i = 1/w_i if w_i < 0, 0 otherwise.
// If "regType" is 1, performs Tikhonov regulariztion, i.e.
// lambda_i = (w_i / alpha w^2 + w_i^2}$$
// where w=max|w_i|.
// If "regType" is 2, performs Landweber regulariztion, i.e.
// lambda_i = {1 - (1- (1 - param w_i^2)^(1/alpha) / w_i.
// The constant is dtermined by
// a = param / sup w_i^2.
// If "regType" is 3, performs cutoff regulariztion, i.e. 
// lambda_i = 1/w_i if w_i >= alpha w <> 0, otherwise 0.

KDMatrix
KDMatrix::pseudoinverse(
	int regType,		// (I) see below
	double alpha,		// (I) see below
	double param)		// (I) see below
 		const
{
static	char	routine[] = "KDMatrix::pseudoinverse";
	double	**u = NULL,	/* [0..m-1][0..n-1] */
		*w = NULL,	/* [0..n-1] */
		**v = NULL;	/* [0..n-1][0..n-1] */
	int	m = rows(),
		n = cols();
	double	**a = (double**) *this;
	KDMatrix invMat;

	int	i, j, k;


    try {

	IF_FAILED_THROW( DrlMatrixSvdDecomp(m, n, a, &u, &w, &v));

	// compute pseudo inverse of w 
	IF_FAILED_THROW( DrlMatrixSvdRegularizeEV(
		n, w, w, regType, alpha, param));



	// compute X = V W^-1 tU b 
	//
	invMat.resize(n, m);
	for (k=0; k<=n-1; k++)
	for (j=0; j<=m-1; j++) {
		invMat[k][j] = 0e0;
		for (i=0; i<=n-1; i++)  {
			invMat[k][j] += v[k][i] * w[i] * u[j][i];
		}
	}

	DrlMatrixFree(u, m, n);
	DrlDoubleVectFree(w, 0, n-1);
	DrlMatrixFree(v, n, n);

	return(invMat);

    }
    catch (KFailure) {
	DrlMatrixFree(u, m, n);
	DrlDoubleVectFree(w, 0, n-1);
	DrlMatrixFree(v, n, n);

	throw KFailure("%s: failed.\n", routine);
    }
}

#endif
