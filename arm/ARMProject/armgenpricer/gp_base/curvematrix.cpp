/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file curvematrix.cpp
 *
 *  \brief
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date December 2006
 */


#include "gpbase/curvematrix.h"
#include "gpbase/curve.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/vectormanip.h"
#include "expt.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_CurveMatrix
///	Routine: Constructor
///	Returns: built the object
///	Action : Create a constant curve matrix
////////////////////////////////////////////////////

ARM_CurveMatrix::ARM_CurveMatrix(const ARM_GP_Matrix& matrix)
{
	itsNbRows = matrix.GetRowsNb();
	itsNbCols = matrix.GetColsNb();

	itsCurves.resize(itsNbRows*itsNbCols);

	for (size_t i = 0; i < itsNbRows; ++i)
		for (size_t j = 0; j < itsNbCols; ++j)
		{
			std::vector<double> abs(1,0.0);
			std::vector<double> ord(1,matrix(i,j));
			ARM_Curve crv;//(abs, ord,new ARM_LinInterpCstExtrapolDble);

			itsCurves(i*itsNbRows+j) = crv;
		}
}

////////////////////////////////////////////////////
///	Class  : ARM_CurveMatrix
///	Routine: Constructor
///	Returns: built the object
///	Action : Create a curve matrix with a vector 
/// of matrix
////////////////////////////////////////////////////

ARM_CurveMatrix::ARM_CurveMatrix(const vector<ARM_Curve*>& curves, size_t nbRows, size_t nbCols)
{
	if (nbRows*nbCols != curves.size())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "The curves vector don't have the good size." );

	itsNbRows = nbRows;
	itsNbCols = nbCols;

	itsCurves.resize(itsNbRows*itsNbCols);

	for (size_t i = 0; i < itsNbRows; ++i)
		for (size_t j = 0; j < itsNbCols; ++j)
		{
			itsCurves[i*itsNbRows+j] = *(curves[i*itsNbRows+j]);
		}
}

////////////////////////////////////////////////////
///	Class  : ARM_CurveMatrix
///	Routine: Interpolate
///	Returns: ARM_GP_Matrix
///	Action : Interpolate on the matrix curve
////////////////////////////////////////////////////

ARM_GP_Matrix ARM_CurveMatrix::Interpolate(double lag) const
{
	ARM_GP_Matrix mat(itsNbRows,itsNbCols);

	for (size_t i = 0; i < itsNbRows; ++i)
		for (size_t j = 0; j < itsNbCols; ++j)
			mat(i,j) = itsCurves[i*itsNbRows+j].Interpolate(lag);

	return mat;
}

////////////////////////////////////////////////////
///	Class  : ARM_CurveMatrix
///	Routine: GetTimes
///	Returns: std::vector<double>
///	Action : Interpolate on the matrix curve
////////////////////////////////////////////////////

std::vector<double> ARM_CurveMatrix::GetTimes() const
{
	std::vector<double> prevVec(0);
	/*std::vector<double>* newVec;

	for (size_t i = 0; i < itsNbRows*itsNbCols;++i)
	{
		std::vector<double> tmpVec = itsCurves[i].GetAbscisses();

		newVec = MergeSortedVectorNoDuplicates(&prevVec,&tmpVec);
		prevVec = *newVec;
		delete newVec;
	}*/

	return prevVec;
}

////////////////////////////////////////////////////
///	Class  : ARM_CurveMatrix
///	Routine: toString
///	Returns: string
///	Action : Interpolate on the matrix curve
////////////////////////////////////////////////////

string ARM_CurveMatrix::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;

	std::vector<double> times = GetTimes();
	ARM_GP_Matrix matrix;

	for(size_t i = 0; i < times.size(); ++i)
	{
		os << "Time " << times[i] << "\n";
		matrix = Interpolate(times[i]);
			
		os << matrix.toString();
		os << "\n";
	}

	return os.str();
}

////////////////////////////////////////////////////
///	Class  : ARM_CurveMatrix
///	Routine: CheckCorrelMatrix
///	Returns: bool
///	Action : Check the correlation matrix
////////////////////////////////////////////////////

bool ARM_CurveMatrix::CheckCorrelMatrix() const
{
	if (itsNbRows != itsNbCols)
		return false;

	return true;
}

////////////////////////////////////////////////////
///	Class  : ARM_CurveCorrelMatrix
///	Routine: SetLabel
///	Returns: void
///	Action : assign a label to each curve of the class
////////////////////////////////////////////////////
void ARM_CurveCorrelMatrix::SetLabel( const vector<string>& lab){

	std::vector<double> tmp; 

	if( itsNbRows!= itsNbCols)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "NbCols should be equal to NbRows to provide a correlation matrix." );

	if( itsNbRows!= lab.size() )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "No correspondance can, be made from this correlation matrix and label vector" );

	int Index=0;
	itsCorrelCurves.resize( itsNbRows*itsNbRows);

	for( int i = 0; i < itsNbRows; i++){
		for ( int j = 0; j < itsNbRows; j++){
			itsCorrelCurves[i*itsNbRows+j] = ARM_CurveWithLabel (itsCurves(i*itsNbRows+j), lab[i]+"|"+lab[j] );
			Index++;
		}
		//tmp = itsCorrelCurves[i*itsNbRows+i].itsCurve.GetOrdinates();
	
		for(int j= 0; j<tmp.size(); j++){
			if( tmp[j] != 1)
				ARMTHROW(ERR_INVALID_ARGUMENT, "the digonal element of a correlation matrix should be equal to 1" );
		}
	}
}

vector<string>	ARM_CurveCorrelMatrix::GetLabel( ) const{

	vector<string>	str;
	string			tmp;
	int				pos;

	for( int i = 0; i < itsNbRows; i++){
		tmp = itsCorrelCurves[i*itsNbRows+i].itsLabel;
		pos = tmp.find("|");
	
		if( pos < 0 || pos > tmp.size() )
			ARMTHROW(ERR_INVALID_ARGUMENT, "label has been correctly inserted in correlation matrix" );

		str.push_back(tmp.substr(0,pos));
	}
	return str;
}
/*

ARM_Curve		ARM_CurveCorrelMatrix::GetLabel( const string & str ) const{

	ARM_Curve tmp;
	for( int i =0; i< itsNbRows*itsNbRows; i++){
		if( itsCorrelCurves[i].itsLabel == str){
			tmp =itsCorrelCurves[i];
			break;
		}
	}
	if (tmp.size()==0)
		ARMTHROW(ERR_INVALID_ARGUMENT," this is not an element of CurveCorrelMatrix");

	return tmp;
}

void	ARM_CurveCorrelMatrix::SetLabel( const string & str, const ARM_Curve & curve ) const{
	
	for( int i =0; i< itsNbRows; i++){
		for( int j =i+1; j< itsNbRows; j++){
			if( itsCorrelCurves[i*itsNbRows+j].itsLabel == str || itsCorrelCurves[j*itsNbRows+i].itsLabel == str){
				itsCorrelCurves[i*itsNbRows+j]=curve;
				itsCorrelCurves[j*itsNbRows+i]=curve;
				break;
			}
		}
	}
}

*/

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

