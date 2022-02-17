/*!
 *
 * Copyright (c) CDC IXIS CM July 2004 Paris
 *
 *	\file surfacelistmodelparam.cpp
 *
 *  \brief 
 *
 *	\author  A. Rajaona
 *	\version 1.0
 *	\date December 2004
 */


/// this headers has to come first
/// as firsttoinc defines pragma for warning to avoid

#include "gpinfra/surfacelistmodelparam.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/surface.h"
#include "gpbase/surfacetypedef.h"
#include "gpbase/numericconstant.h"


/// ARM Kernel
#include <glob/expt.h>			/// necessary for exception throwing

/// STL
#include <iomanip> /// for setprecision()
#include <cmath>

CC_BEGIN_NAMESPACE(ARM)

#define K_STD_DOUBLE_TOL 1e-14

const int K_MAX_ROW_SURFACE = 120;
const int K_MAX_COLUMN_SURFACE = 120;
const int INTERMEDIATESTEPS = 2;

////////////////////////////////////////////////////
///	Class  : ARM_SurfaceListModelParam
///	Routine: Constructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_SurfaceListModelParam::ARM_SurfaceListModelParam( 
	ParamType type, 
	const std::vector<double>& index,
	const ARM_SurfacePtrVector& surfacelist, 
	const string& paramName,
	double lowerBound,
	double upperBound,
	bool adviseBreakPointTimes )
:
	ARM_ModelParam(type,adviseBreakPointTimes),
	itsIndex(index),
	itsSurfaceList(surfacelist),
	itsLowerBound(lowerBound),
    itsUpperBound(upperBound)
{
	if( itsIndex.size() != itsSurfaceList.size() )
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME + ": itsIndex.size() != itsSurfaceList.size()!" );
	if( !itsIndex.size() )
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME + ": itsIndex.size() ==0!" );

	CheckVectorStrictlyIncreasing(itsIndex,"Index","ARM_SurfaceListModelParam::ARM_SurfaceListModelParam",__LINE__,__FILE__);
	
	if(		lowerBound == -ARM_NumericConstants::ARM_BIGGEST_POSITIVE_NUMBER
		&&  upperBound ==  ARM_NumericConstants::ARM_BIGGEST_POSITIVE_NUMBER )
	{
		itsLowerBound = ARM_CalibParamCst::ModelParamStdLowerBound( GetType() );
		itsUpperBound = ARM_CalibParamCst::ModelParamStdUpperBound( GetType() );
	}
};
	

////////////////////////////////////////////////////
///	Class  : ARM_SurfaceListModelParam
///	Routine: Copy constructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////

ARM_SurfaceListModelParam::ARM_SurfaceListModelParam( const ARM_SurfaceListModelParam& rhs )
:	ARM_ModelParam( rhs ), 
	itsIndex( rhs.itsIndex ), 
	itsLowerBound( rhs.itsLowerBound ), 
	itsUpperBound( rhs.itsUpperBound ) 
{
	// surface list deep copy
	itsSurfaceList.resize ( rhs.itsSurfaceList.size() );

	for (size_t i(0); i<itsSurfaceList.size(); i++)
		itsSurfaceList[i] = ARM_SurfacePtr ((ARM_Surface*)rhs.itsSurfaceList[i]->Clone());

}

////////////////////////////////////////////////////
///	Class  : ARM_SurfaceListModelParam
///	Routine: assignment operator
///	Returns: 
///	Action : 
////////////////////////////////////////////////////

ARM_SurfaceListModelParam& ARM_SurfaceListModelParam::operator =( const ARM_SurfaceListModelParam& rhs )
{
	if (&rhs != this)
	{ 
		this->~ARM_SurfaceListModelParam();
		new (this) ARM_SurfaceListModelParam (rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_SurfaceListModelParam
///	Routine: ~ARM_SurfaceListModelParam
///	Returns: 
///	Action : destructor
////////////////////////////////////////////////////

ARM_SurfaceListModelParam::~ARM_SurfaceListModelParam()
{}


///////////////////////////////////////
///	Class  : ARM_SurfaceListModelParam
///	Routine: Clone
///	Returns: ARM_Object*
///	Action : Clone the object
///////////////////////////////////////

ARM_Object* ARM_SurfaceListModelParam::Clone() const
{
	return new ARM_SurfaceListModelParam(*this); 
}


////////////////////////////////////////////////////
///	Class   : ARM_SurfaceListModelParam
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_SurfaceListModelParam::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
	os << "\n"<< indent << "Model Param Surface List Type: " << GetTypeString() << "\n";
	os << " lower bound = " << itsLowerBound << "\n";
	os << " upper bound = " << itsUpperBound << "\n";

	for (size_t i = 0; i < itsSurfaceList.size(); ++i)
	{
		os << "Index: " << itsIndex[i] << std::endl;		
		int rowSize = itsSurfaceList[i]->GetX1().size(); 
		int columnSize = itsSurfaceList[i]->GetX2().size();
		int affichRowSize = rowSize < K_MAX_ROW_SURFACE ? rowSize :  K_MAX_ROW_SURFACE;
		int affichColumnSize = columnSize < K_MAX_COLUMN_SURFACE ? columnSize :  K_MAX_COLUMN_SURFACE;
		bool isDumped=false;
		if( ((ARM_SurfaceWithInterpol*) (&(itsSurfaceList[i]))) && (rowSize!=0) && (columnSize!=0)  )
		{
// FIXMEFRED: mig.vc8 (31/05/2007 10:41:45): floor(int/int) doesnt exist
			int stepRow = /*static_cast<int>(floor(*/rowSize/affichRowSize/*))*/;
			int stepCol = /*static_cast<int>(floor(*/columnSize/affichColumnSize/*))*/;
			if((affichRowSize!=rowSize) || (affichColumnSize!=columnSize))
			{
				os << "here is a sample of the surface : "<< itsIndex[i] << std::endl;
				os << " Model Param Surface with rows       = " << rowSize << " cols = " << columnSize << "\n";
				os << " whereas the displayed one with rows = " << affichRowSize << " cols = " << affichColumnSize << "\n";

				
				ARM_GP_Vector x1(affichRowSize);
				ARM_GP_Vector x2(affichColumnSize);
				ARM_GP_Matrix x3(affichRowSize,affichColumnSize);
							
				ARM_InterpolType type = ARM_InterpolationType::linear_row_extrapoleCst;
				double defaultvalue = 0.0;//((itsSurfaceList[i]))->GetDefaultValue();

				for( size_t j=0; j<affichRowSize; ++j )
				{
					x1[j] = itsSurfaceList[i]->GetX1()[j*(stepRow)];
					for( size_t k=0; k<affichColumnSize; ++k )
					{
						x2[k] = itsSurfaceList[i]->GetX2()[k*(stepCol)];
						x3(j,k) = itsSurfaceList[i]->GetX3()(j*(stepRow),k*(stepCol));
					}
				}
				
				ARM_SurfaceWithInterpol surface (x1,
										 x2, 
										 x3,
										 type,defaultvalue);
				os << surface.toString();
				isDumped=true;
			}
		}
		if(!isDumped)
		{
			/// Check if necessary to dump all values
			if(rowSize > 0 && columnSize > 0)
			{
				double theVal = itsSurfaceList[i]->InterpolateAtPoint(0,0);
				if(*(itsSurfaceList[i]) == theVal)
					os << " constant values = " << theVal << "\n";
				else
					os << itsSurfaceList[i]->toString();
			}
			else
				os << " empty surface" << "\n";
		}
	}
    return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_SurfaceListModelParam
///	Routine: rows
///	Returns: size_t
///	Action : returns the size of the model param
////////////////////////////////////////////////////
size_t ARM_SurfaceListModelParam::rows() const
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "rows with only one argument unimplemented for a surface!" );
}

////////////////////////////////////////////////////
///	Class  : ARM_SurfaceListModelParam
///	Routine: cols
///	Returns: size_t
///	Action : returns the size of the model param
////////////////////////////////////////////////////
size_t ARM_SurfaceListModelParam::cols() const
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "cols with only one argument unimplemented for a surface!" );
}

////////////////////////////////////////////////////
///	Class  : ARM_SurfaceListModelParam
///	Routine: size
///	Returns: size_t
///	Action : returns the size of the model param
////////////////////////////////////////////////////
size_t ARM_SurfaceListModelParam::size() const
{
	return itsSurfaceList.size();  
}

////////////////////////////////////////////////////
///	Class  : ARM_SurfaceListModelParam
///	Routine: rows
///	Returns: size_t
///	Action : returns the size of the model param
////////////////////////////////////////////////////
size_t ARM_SurfaceListModelParam::rows(size_t surfIdx) const
{
	return GetSurface(surfIdx)->GetX1().size();  
}

////////////////////////////////////////////////////
///	Class  : ARM_SurfaceListModelParam
///	Routine: cols
///	Returns: size_t
///	Action : returns the size of the model param
////////////////////////////////////////////////////
size_t ARM_SurfaceListModelParam::cols(size_t surfIdx) const
{
	return GetSurface(surfIdx)->GetX2().size();  
}


////////////////////////////////////////////////////
///	Class  : ARM_SurfaceListModelParam
///	Routine: GetValue
///	Returns: double
///	Action : Various getvalue functions
////////////////////////////////////////////////////

double ARM_SurfaceListModelParam::GetValue( double x1 ) const
{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "GetValue with only one argument unimplemented for a surface!" ); }

double ARM_SurfaceListModelParam::GetValue( double x1, double x2 ) const
{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "GetValue with only two arguments unimplemented for a list of surface!" ); }

double ARM_SurfaceListModelParam::GetValue( double x1, double x2, double x3 ) const
{	return Interpolate(x1,x2,x3); }

double ARM_SurfaceListModelParam::GetValue( size_t i, double x1, double x2 ) const
{	return itsSurfaceList[i]->Interpolate(x1,x2); }

double ARM_SurfaceListModelParam::GetValueAtPoint( size_t i ) const
{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "GetValue with only one argument unimplemented for a surface!" ); }

double ARM_SurfaceListModelParam::GetTimeAtPoint( size_t i ) const
{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "GetValue with only one argument unimplemented for a surface!" ); }
double ARM_SurfaceListModelParam::GetTenorAtPoint( size_t i) const
{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "GetTenorAtPoint(i) unimplemented for a curve!" ); }

double ARM_SurfaceListModelParam::GetValueAtPoint( size_t i, size_t j ) const
{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "GetValue with only two arguments unimplemented for a list of surface!" ); }

double ARM_SurfaceListModelParam::GetTimeAtPoint( size_t i, size_t j ) const
{	return itsSurfaceList[i]->GetX1()[j]; }

double ARM_SurfaceListModelParam::GetTenorAtPoint( size_t i, size_t j ) const
{	return itsSurfaceList[i]->GetX2()[j]; }

double ARM_SurfaceListModelParam::GetValueAtPoint( size_t i, size_t j, size_t k) const
{	return itsSurfaceList[i]->InterpolateAtPoint(j,k); }


////////////////////////////////////////////////////
///	Class  : ARM_SurfaceListModelParam
///	Routine: SetValue
///	Returns: void
///	Action : Various SetValue functions
////////////////////////////////////////////////////

void ARM_SurfaceListModelParam::SetValue( double x1, double value )
{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "SetValue with only one argument unimplemented for a surface list!" ); }

void ARM_SurfaceListModelParam::SetValue( double x1, double x2, double value  )
{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "SetValue with only two arguments unimplemented for a surface list!" ); }

void ARM_SurfaceListModelParam::SetValue( double x1, double x2, double x3, double value  )
{	
	// itsSurfaceList->insert( x1, x2, x3, value ); 
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "SetValue unimplemented for a surface list!" );
}

void ARM_SurfaceListModelParam::SetValue( size_t i, double x1, double x2, double value  )
{	
	itsSurfaceList[i]->insert (x1, x2, value);
}

void ARM_SurfaceListModelParam::SetValueAtPoint( size_t i, double value )
{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "SetValue with only one argument unimplemented for a surface list!" ); }

void ARM_SurfaceListModelParam::SetValueAtPoint( size_t i, size_t j, double value )
{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "SetValue with only two arguments unimplemented for a surface list!" ); }

void ARM_SurfaceListModelParam::SetValueAtPoint( size_t i, size_t j, size_t k, double value ) 
{	itsSurfaceList[i]->insertAtPoint(j,k,value); }


///////////////////////////////////////////////////
///	Class  : ARM_SurfaceListModelParam
///	Routine: Interpolate
///	Returns: 
///	Action : Direct implementation of interpolation
///////////////////////////////////////////////////
	
double ARM_SurfaceListModelParam::Interpolate( double time, double tenor, double K_F ) const
{

	int i= CC_NS(std,lower_bound)(itsIndex.begin(), itsIndex.end(), tenor ) - itsIndex.begin();
	
	/// in the range
	if( i>0 && i<itsIndex.size() )
	{
		if( fabs( itsIndex[i]-tenor) < K_STD_DOUBLE_TOL )
			return itsSurfaceList[i]->Interpolate(time, K_F);
		else
		{
			double v_p	= itsSurfaceList[i-1]->Interpolate(time, K_F);
			double v_n	= itsSurfaceList[i]->Interpolate(time, K_F);
			double x_p	= itsIndex[i-1];
			double x_n	= itsIndex[i];
			double v_i	= v_p + (tenor - x_p) * (v_n - v_p)/(x_n - x_p);
			return v_i;
		}
	}
	else if( i==0)
		return itsSurfaceList[0]->Interpolate(time, K_F);
	else /// if( i==itsIndex.size())
		return itsSurfaceList[itsIndex.size()-1]->Interpolate(time, K_F);

/*
	size_t k=0;
	size_t maxIdxSize, i_p=0, i_n=0, i=0;
	double v_p=0., v_n=0., v_i=0., x_n=0., x_p=0.;

	maxIdxSize = itsIndex.size()-1;

	// Find i, (ip, x1p) and (in, x1n) satisfying x1p <= x1 <= x1n
	// with special case : i = 0 if  x1<= itsIndex[0] 
	// and i = itsIndex->Size()-1 if x1 >= itsIndex[maxSize] with maxSize = itsIndex->GetSize() 

	x_p = itsIndex[0];
	x_n = x_p;
	i_n = i_p;

	while (i < maxIdxSize)
	{
		if (itsIndex[i] < x1)
		{
			i_p = i;
			i_n = i_p + 1;
			x_n = itsIndex[i_n];
			x_p = itsIndex[i_p];
			i++;
		}
		else
			break;
	}

	if (x_n <= x1)
	{
		x_p = x_n;
		i_p = i_n;
	}

	if (i_p != i_n)
	{
		v_p = itsSurfaceList[i_p]->Interpolate(x2, x3);
		v_n = itsSurfaceList[i_n]->Interpolate(x2, x3);
		v_i = v_p + (x1 - x_p) * (v_n - v_p)/(x_n - x_p);
	}
	else
		v_i = itsSurfaceList[i_p]->Interpolate(x2, x3);

	return(v_i);
*/
}


////////////////////////////////////////////////////////////
///	Class  : ARM_SurfaceListModelParam
///	Routine: InterpolateAtPoint
///	Returns: 
///	Action : Direct implementation of interpolation at point
////////////////////////////////////////////////////////////

double ARM_SurfaceListModelParam::InterpolateAtPoint( size_t i, size_t j, size_t k ) const
{
	// i is index of itsIndex
	// j, k are index of itsSurfaceList[i] 
	return itsSurfaceList[i]->InterpolateAtPoint(j, k);
}



///////////////////////////////////////
///	Class  : ARM_SurfaceListModelParam
///	Routine: SetAndUpDate
///	Returns: 
///	Action : Replace old values using interpolate if necessairy
///////////////////////////////////////
	
void ARM_SurfaceListModelParam::SetAndUpDate(std::vector<double>*values)
{
	for (size_t i=0, k=0; i<itsSurfaceList.size(); i++)
	{
		for(size_t j=0; j<itsSurfaceList[i]->GetX3().rows(); ++j )
			for(size_t p=0; p<itsSurfaceList[i]->GetX3().cols(); ++p, ++k)
				itsSurfaceList[i]->insert(itsSurfaceList[i]->GetX1()[j],itsSurfaceList[i]->GetX2()[p],(*values)[k]);
	}
}

///////////////////////////////////////
///	Class  : ARM_SurfaceListModelParam
///	Routine: MergeModelParam
///	Returns: 
///	Action : merge a model param with another model param of same type
///         Warning! this assumes that the newValue has the same type
///////////////////////////////////////

void ARM_SurfaceListModelParam::MergeModelParam(ARM_ModelParam* NewValue )
{
	ARM_SurfaceListModelParam* modelParam = dynamic_cast<ARM_SurfaceListModelParam*>(NewValue);

	if( !modelParam)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
		    "trying to merge model param that is not a surface list with a model param surface list!" );

/// theoretically this should not be tested 
/// but in strict validation mode
/// we test and retest!
#if defined( __GP_STRICT_VALIDATION )
    if(modelParam->GetType() != GetType() )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
		    "trying to merge model param of a different type!" );
#endif
    for(size_t i=0; i<modelParam->size(); ++i)
	{
		ARM_FlatSurface* surface = dynamic_cast<ARM_FlatSurface*>(&*(itsSurfaceList[i]));
		if(!surface)
			itsSurfaceList[i]->insert( &*(modelParam->GetSurface(i)) );
		else
			itsSurfaceList[i]= ARM_SurfacePtr((ARM_Surface*)(modelParam->GetSurface(i))->Clone());
	}
}

///////////////////////////////////////
///	Class  : ARM_SurfaceListModelParam
///	Routine: SetSurfaceList
///	Returns: 
///	Action : Set the surface of a surface model param!
///////////////////////////////////////

void ARM_SurfaceListModelParam::SetSurfaceList( const ARM_SurfacePtrVector& surfacelist )
{ 
	itsSurfaceList = surfacelist; 
}


void ARM_SurfaceListModelParam::SetSurfaceListIndex( const std::vector<double>& index )
{ 
	if( itsIndex.size() != index.size() )
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME + ": itsIndex.size() != index.size()!" );
	itsIndex=index; 
}


std::vector<double>* ARM_SurfaceListModelParam::GetData( ARM_DataType type, double index, long& rows, long& cols ) const
{
	return new std::vector<double>();
	/*size_t i = 0; 
	while (i < itsIndex.size())
	{
		if (index == itsIndex[i])
			break;
		else
			i++;
	}

	if (i  >= itsIndex.size())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
            "Unknown Index... (Check your Data)" );
	else
	{

	switch(type)
		{
		case ARM_ModelParamType::BreakPointTimes:
			{
				const ARM_GP_Vector& result = itsSurfaceList[i]->GetX1();
				rows=result.size();
				cols=1;
				return result.GetValues();
			}
        case ARM_ModelParamType::Tenors:
            {
				const std::vector<double>& result = itsSurfaceList[i]->GetX2();
				rows=result.size();
				cols=1;
				return new std::vector<double>(result);
			}
        case ARM_ModelParamType::Values:
			{
				const ARM_GP_Matrix& values = itsSurfaceList[i]->GetX3();
				rows = values.rows();
				cols = values.cols();
				std::vector<double>& result = new std::vector<double>( values.size() );
				CC_NS(std,copy)(values.begin(), values.end(), result->begin());
				return result;
			}

        default:
            throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
            "Unknown type... permitted are BreakPointTimes,Tenors or Values" );
		}
	}*/
}

#undef K_STD_DOUBLE_TOL

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

