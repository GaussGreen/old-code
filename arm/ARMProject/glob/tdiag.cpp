
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "error.h"
#include "LinAlg.h"



/* *********************************** /SOH/ ***********************************
	SYNOPSIS	ARM_TDiag::ARM_TDiag(int size)
	Constructor
* *********************************** /EOH/ ***********************************/
ARM_TDiag::ARM_TDiag(int size)
{
	if (!Resize(size)) SetError(kErrMemoryAlloc, "ARM_TDiag::ARM_TDiag", NULL);
}


/* *********************************** /SOH/ ***********************************
	SYNOPSIS	ARM_TDiag::ARM_TDiag(ARM_TDiag &t)
	Constructor (copy)
* *********************************** /EOH/ ***********************************/
ARM_TDiag::ARM_TDiag(ARM_TDiag &t) : ARM_GenMatrix(t)
{
	int	i;
	
	if (!Resize(t.GetNumLines())) {	
		SetError(kErrMemoryAlloc, "ARM_TDiag::ARM_TDiag(ARM_TDiag &)", NULL);
		return;
	}
	
	for (i=1; i<t.GetNumLines(); i++) Elt(i-1, i) = t.Elt(i, i-1);
	for (i=0; i<t.GetNumLines(); i++) Elt(i, i) = t.Elt(i, i);
	for (i=0; i<t.GetNumLines()-1; i++) Elt(i, i+1) = t.Elt(i, i+1);
}


/* *********************************** /SOH/ ***********************************
	SYNOPSIS	ARM_TDiag & ARM_TDiag::operator = (ARM_TDiag &t)

* *********************************** /EOH/ ***********************************/
ARM_TDiag & ARM_TDiag::operator = (ARM_TDiag &t)
{	
	(*this).ARM_Obj::operator = (t);
	
	itsSize = t.itsSize;
	itsLower = t.itsLower;
	itsDiag = t.itsDiag;
	itsUpper = t.itsUpper;
	
	if (!itsLower.GetErrorStatus() || !itsDiag.GetErrorStatus() ||!itsUpper.GetErrorStatus())
		SetError(kErrAssign, "ARM_TDiag::operator = (ARM_TDiag &)", NULL);
	return(*this);
}



/* *********************************** /SOH/ ***********************************
	SYNOPSIS	int ARM_TDiag::Resize(int size)

	Resize tridiagonal matrix. Elements are set to 0. Return 0 if failed, 1 if not.
* *********************************** /EOH/ ***********************************/
int ARM_TDiag::Resize(int size)
{
	int	i, j;
	
	ASSERT_ERROR_STATUS("ARM_TDiag::Resize(int)", 0);

	if (!itsLower.Resize(size) || !itsDiag.Resize(size) || !itsUpper.Resize(size)) {
		SetError(kErrMemoryAlloc,"ARM_TDiag::Resize(int)", NULL);
		return(0);
	}
	itsSize = size;	
	return(1);
}


/* *********************************** /SOH/ ***********************************
	SYNOPSIS	double &ARM_TDiag::Elt(int i, int j)

	Returns reference to element (i,j)
* *********************************** /EOH/ ***********************************/
double &ARM_TDiag::Elt(int i, int j)
{
	static	double	e0, reterr;
	reterr = kHugeDouble;

	ASSERT_ERROR_STATUS("ARM_TDiag::Elt(int, int)", reterr);
	if (i<0 || j<0 || i>itsSize-1 || j>itsSize-1) {	
		SetError(kErrInvalidData, "ARM_TDiag::Elt(int,int)", NULL);
		return(reterr);
	}

	switch (i-j) {
		case (1) : return(itsLower[i]);
		case (0) : return(itsDiag[i]);
		case (-1) : return(itsUpper[i]);
		default: 
			e0 = 0.0;
			return(e0);
	}
}


/* *********************************** /SOH/ ***********************************
	SYNOPSIS	double &ARM_TDiag::LowerElt(int i)

	Returns reference to lower element i.
* *********************************** /EOH/ ***********************************/
double &ARM_TDiag::LowerElt(int i)
{
	static	double	reterr;
	reterr = 0.0;
	ASSERT_ERROR_STATUS("ARM_TDiag::LowerElt(int)", reterr);
	if (i<0 || i+1>itsSize) {
		SetError(kErrInvalidData, "ARM_TDiag::LowerElt(int)", NULL);
		return(reterr);
	}
	return(itsLower[i]);
}

/* *********************************** /SOH/ ***********************************
	SYNOPSIS	double &ARM_TDiag::DiagElt(int i)

	Returns reference to diagonal element i.
* *********************************** /EOH/ ***********************************/
double &ARM_TDiag::DiagElt(int i)
{
	static	double	reterr;
	reterr = 0.0;
	ASSERT_ERROR_STATUS("ARM_TDiag::DiagElt(int)", reterr);
	if (i<0 || i+1>itsSize) {
		SetError(kErrInvalidData, "ARM_TDiag::DiagElt(int)", NULL);
		return(reterr);
	}
	return(itsDiag[i]);
}

/* *********************************** /SOH/ ***********************************
	SYNOPSIS	double &ARM_TDiag::UpperElt(int i)

	Returns reference to upper element i.
* *********************************** /EOH/ ***********************************/
double &ARM_TDiag::UpperElt(int i)
{
	static	double	reterr;
	reterr = 0.0;
	ASSERT_ERROR_STATUS("ARM_TDiag::UpperElt(int)", reterr);
	if (i<0 || i+1>itsSize) {
		SetError(kErrInvalidData, "ARM_TDiag::UpperElt(int)", NULL);
		return(reterr);
	}
	return(itsUpper[i]);
}


/* *********************************** /SOH/ ***********************************
	SYNOPSIS	ARM_TDiag & ARM_TDiag::operator -()

* *********************************** /EOH/ ***********************************/
ARM_TDiag & ARM_TDiag::operator -()
{
	ASSERT_ERROR_STATUS("ARM_TDiag::operator -()", *this);
	itsLower = -itsLower;
	itsDiag = -itsDiag;
	itsUpper = -itsUpper;
	return(*this);
}



/* *********************************** /SOH/ ***********************************
	SYNOPSIS	int ARM_TDiag::Transpose(void)

	Transpose tri diagonal matrix
* *********************************** /EOH/ ***********************************/
int ARM_TDiag::Transpose(void)
{
	double	tmp;
	int	i;
	ASSERT_ERROR_STATUS("ARM_TDiag::Transpose(void)", 0);
	for (i=0; i<itsSize-1; i++) {
		tmp = itsUpper[i];
		itsUpper[i] = itsLower[i+1];
		itsLower[i+1] = tmp;
	}
	return(1);
}


/* *********************************** /SOH/ ***********************************
	SYNOPSIS	double ARM_TDiag::Det(void)

	Returns determinant of matrix, kHugeDouble if failed.
* *********************************** /EOH/ ***********************************/
double ARM_TDiag::Det(void)
{
	int	i;
	double	det, r;
	static	ARM_TDiag	tmp;
	
	ASSERT_ERROR_STATUS("ARM_TDiag::Det(void)", 0);

	tmp.ClearError();
	if (!tmp.Resize(itsSize)) {
		setError(kErrMemoryAlloc, "ARM_TDiag::Det", NULL);
		return(kHugeDouble);
	}
	
	tmp = *this;
	
	for (i=itsSize-2; i>=0; i--) {
		r = tmp.UpperElt(i) / tmp.DiagElt(i+1);
		tmp.DiagElt(i) -= r * tmp.LowerElt(i+1);
	}

	det = 1.0;
	for (i=0; i<itsSize; i++) det *= tmp.DiagElt(i);
	
	return(det);

}



/* *********************************** /SOH/ ***********************************
	SYNOPSIS	int ARM_TDiag::Invert(ARM_Matrix *invert, double &det)

	Returns invert of matrix into invert and determinant in det. Returns
	0 if failed, 1 if not.
* *********************************** /EOH/ ***********************************/
int ARM_TDiag::Invert(ARM_Matrix *invert, double &det)
{
	int	i, j;
	double	d, x;
	static	ARM_TDiag	tmp;
	
	ASSERT_ERROR_STATUS("ARM_TDiag::Invert(ARM_Matrix *, double &)", 0);
	if (invert->GetNumLines() != itsSize || invert->GetNumCols() != itsSize) {
		SetError(kErrInvalidData, "ARM_TDiag::Invert(ARM_Matrix *, double &)", NULL);
		return(0);
	}
	
	tmp.ClearError();
	tmp = *this;
	if (!tmp.GetErrorStatus()) {
		SetError(kErrMemoryAlloc, "ARM_TDiag::Invert(ARM_Matrix *, double &)", NULL);
		return(0);
	}
	
	invert->SetToId();
	

	for (i=itsSize-2; i>=0; i--) {
		d = tmp.UpperElt(i) / tmp.DiagElt(i+1);
	
		tmp.DiagElt(i) -= d * tmp.LowerElt(i+1);
		
		for (j=0; j<itsSize; j++) {
			x = invert->Elt(i, j) - d * invert->Elt(i+1, j);
			invert->Elt(i, j) = x;
		}
	}

	for (j=0; j<itsSize; j++) {
		x = invert->Elt(0,j) / tmp.DiagElt(0);
		invert->Elt(0, j) = x;
	}

	for (i=1; i<itsSize; i++) {
		d = tmp.itsLower[i] / tmp.DiagElt(i);
		for (j=0; j<itsSize; j++) {
			x = invert->Elt(i,j) / tmp.DiagElt(i) - d * invert->Elt(i-1, j);
			invert->Elt(i, j) = x;
		}
	}

	det = 1.0;
	for (i=0; i<itsSize; i++) det *= tmp.DiagElt(i);

	return(1);
}  



/* *********************************** /SOH/ ***********************************
	SYNOPSIS	int ARM_TDiag::LinSolve(ARM_Vector *x, double &det)

	Solves linear system with x as right side. Returns result in x
	and determinant in det. Returns 0 if failed, 1 if not.
* *********************************** /EOH/ ***********************************/
int ARM_TDiag::LinSolve(ARM_Vector *x, double &det)
{
	int	i;
	double	r, c;
	static	ARM_TDiag	tmp;
	
	ASSERT_ERROR_STATUS("ARM_TDiag::LinSolve(ARM_Vector *, double &)", 0);
	if (x->GetSize() != itsSize) {
		SetError(kErrInvalidData, "ARM_TDiag::LinSolve(ARM_Vector *, double &)", NULL);
		return(0);
	}
	
	tmp.ClearError();
	tmp = *this;
	if (!tmp.GetErrorStatus()) {
		SetError(kErrAssign, "ARM_TDiag::LinSolve(ARM_Vector *, double &)", NULL);
		return(0);
	}
	
	

	for (i=itsSize-2; i>=0; i--) {
		r = tmp.UpperElt(i) / tmp.DiagElt(i+1);
		tmp.DiagElt(i) -= r * tmp.LowerElt(i+1);
		(*x)[i] -=  r * (*x)[i+1];
	}

	(*x)[0] /= tmp.DiagElt(0);
	
	for (i=1; i<itsSize; i++) {
		c = (*x)[i] - tmp.itsLower[i] * (*x)[i-1];
		c /= tmp.DiagElt(i);
		(*x)[i] = c;
	}	
		
	det = 1.0;
	for (i=0; i<itsSize; i++) det *= tmp.DiagElt(i);
	
	return(1);

}









