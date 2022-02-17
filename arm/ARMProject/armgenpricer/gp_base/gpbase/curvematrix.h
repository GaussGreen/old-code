/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file curvematrix.h
 *
 *  \brief files to define a curve matrix
 *
 *	\author  Richard Guillemot
 *	\version 1.0
 *	\date December 2006
 */


#ifndef _INGPBASE_CURVEMATRIX_H
#define _INGPBASE_CURVEMATRIX_H

/// this header comes firts as it includes some preprocessor constants!
#include "removeidentifiedwarning.h"

/// use our macro for namespace
#include "gpbase/assignop.h"

#include "port.h"
#include "curvetypedef.h"
#include "curve.h"


CC_BEGIN_NAMESPACE( ARM )

/// class for conversion of string to the corresponding enum!
class ARM_CurveMatrix : public ARM_RootObject
{
protected:
	ARM_GP_T_Vector<ARM_Curve> itsCurves;
	size_t itsNbRows;
	size_t itsNbCols;
public:
	ARM_CurveMatrix(){}
	ARM_CurveMatrix(const ARM_GP_Matrix& matrix);
	ARM_CurveMatrix(const vector<ARM_Curve*>& curves, size_t nbRows, size_t nbCols);

	ASSIGN_OPERATOR(ARM_CurveMatrix)

	ARM_CurveMatrix(const ARM_CurveMatrix& rhs){ Copy(rhs); }

	virtual void Copy( const ARM_CurveMatrix& rhs){
	itsCurves = rhs.itsCurves;
	itsNbRows = rhs.itsNbRows;
	itsNbCols = rhs.itsNbCols;
	}

	virtual ~ARM_CurveMatrix() {};

	virtual ARM_Object* Clone() const { return new ARM_CurveMatrix(*this); }

	ARM_GP_Matrix Interpolate(double lag) const;

	std::vector<double> GetTimes() const;

	bool CheckCorrelMatrix() const;
	size_t rows() const { return itsNbRows; }
	size_t cols() const { return itsNbCols; }
	bool empty() const { return ((itsNbRows==0)&&(itsNbCols)); }

	virtual string toString(const string& indent="", const string& nextIndent="") const;
};

struct ARM_CurveWithLabel{

public:
	ARM_CurveWithLabel(){}
	ARM_CurveWithLabel(const ARM_Curve& cur, const string& lab){ itsCurve= cur; itsLabel = lab;	}
	ARM_CurveWithLabel( const ARM_CurveWithLabel & rhs):itsCurve(rhs.itsCurve),itsLabel(rhs.itsLabel){}
	virtual ~ARM_CurveWithLabel() {};
	ASSIGN_OPERATOR(ARM_CurveWithLabel)

	ARM_Curve itsCurve;
	string itsLabel;
};

class ARM_CurveCorrelMatrix : public ARM_CurveMatrix
{
protected:
	vector<ARM_CurveWithLabel> itsCorrelCurves;

public:
	ARM_CurveCorrelMatrix( ){ }
	ARM_CurveCorrelMatrix(const ARM_GP_Matrix& matrix, const vector<string>& lab):ARM_CurveMatrix(matrix){	SetLabel(lab);	}
	ARM_CurveCorrelMatrix(const vector<ARM_Curve*>& curves, size_t nbRows, size_t nbCols, const vector<string>& lab):ARM_CurveMatrix(curves, nbRows, nbCols){	SetLabel(lab);	}

	ASSIGN_OPERATOR(ARM_CurveCorrelMatrix)

	virtual void Copy( const ARM_CurveCorrelMatrix& rhs){
		ARM_CurveMatrix::Copy(rhs);
		itsCorrelCurves = rhs.itsCorrelCurves;
	}

	ARM_CurveCorrelMatrix(const ARM_CurveCorrelMatrix& rhs){ Copy(rhs); }
	virtual ~ARM_CurveCorrelMatrix() {};
	virtual ARM_Object* Clone() const { return new ARM_CurveCorrelMatrix(*this); }

	void							SetLabel( const vector<string>& lab);
	vector<string>					GetLabel( ) const;
	void							SetCorrelCurves( const vector<ARM_CurveWithLabel> & vec ){ itsCorrelCurves = vec; }
	vector<ARM_CurveWithLabel>		GetCorrelCurves( ) const { return itsCorrelCurves; }



/*	ARM_Curve		GetLabel( const string & ) const;
	void			SetLabel( const string & , const ARM_Curve & ) const;*/
};




CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
