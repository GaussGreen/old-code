
#include "gpinfra/gramfunctorarg.h"
#include "gpbase/ostringstream.h"
#include "gpbase/env.h"
#include "gpbase/gpvector.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/datestrip.h"

#include <glob/expt.h>


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArg
///	Routine: CheckType
///	Returns: void
///	Action : throw an exception if type different from 
///			itsType
////////////////////////////////////////////////////

void ARM_GramFctorArg::CheckType(Type t) const
{
#ifdef __GP_STRICT_VALIDATION
	if( itsType != t ) 
		{ 
			CC_Ostringstream os; 
			os << "looking for type " << TypeToString(t) << " but found " << TypeToString(itsType)
				<< " Please advise\n"; 
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() ); 
		}
#else
	/// nothing
#endif
}


////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArg
///	Routine: Constructor with a double
///	Returns: void
///	Action : builds the object.. 	by default an
///			object is a double with value 0
////////////////////////////////////////////////////

ARM_GramFctorArg::ARM_GramFctorArg( double d )
:	itsType(GFAT_DOUBLE_TYPE), itsDouble(d)
{}


///////////////////////////////////////////////////
///	Class  : ARM_GramFctorArg
///	Routine: Constructor with a string
///	Returns: 
///	Action : build the object
////////////////////////////////////////////////////

ARM_GramFctorArg::ARM_GramFctorArg( const string& s)
:	itsType(GFAT_STRING_TYPE), itsString(s)
{}


////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArg
///	Routine: Constructor with a date
///	Returns: 
///	Action : build the object with ARM_Date
////////////////////////////////////////////////////

ARM_GramFctorArg::ARM_GramFctorArg( const ARM_Date& d )
:	itsType(GFAT_DATE_TYPE), itsDate(d)
{}


////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArg
///	Routine: Constructor with a vector
///	Returns:  
///	Action : build the object with std::vector<double>
////////////////////////////////////////////////////
ARM_GramFctorArg::ARM_GramFctorArg( const ARM_VectorPtr& v )
:	itsType(GFAT_VECTOR_TYPE), itsVector(ARM_GP_VectorPtr(new ARM_GP_Vector(*v)))
{}

ARM_GramFctorArg::ARM_GramFctorArg( const ARM_GP_VectorPtr& v )
:	itsType(GFAT_VECTOR_TYPE), itsVector(v)
{}


////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArg
///	Routine: Constructor with a matrix
///	Returns:  
///	Action : build the object with ARM_GP_Matrix
////////////////////////////////////////////////////

ARM_GramFctorArg::ARM_GramFctorArg( const ARM_GP_MatrixPtr& m )
:	itsType(GFAT_MATRIX_TYPE), itsMatrix(m)
{}

////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArg
///	Routine: Constructor with a curve
///	Returns: void 
///	Action : build the object with ARM_GP_Curve
////////////////////////////////////////////////////

ARM_GramFctorArg::ARM_GramFctorArg( const ARM_GP_CurvePtr& c)
:	itsType(GFAT_CURVE_TYPE), itsCurve(c)
{}


////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArg
///	Routine: Constructor with a string vector
///	Returns: void 
///	Action : build the object with ARM_StringVectorPtr
////////////////////////////////////////////////////

ARM_GramFctorArg::ARM_GramFctorArg( const ARM_StringVectorPtr& sv)
:	itsType(GFAT_STRINGVEC_TYPE), itsStringVector(sv)
{}

////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArg
///	Routine: Constructor with a date strip
///	Returns: void 
///	Action : build the object with ARM_DateStripPtr
////////////////////////////////////////////////////

ARM_GramFctorArg::ARM_GramFctorArg( const ARM_DateStripPtr& sched)
:	itsType(GFAT_DATESTRIP_TYPE), itsDateStrip(sched)
{}


////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArg
///	Routine: CopyNoCleanUp
///	Returns: void
///	Action : Copy the various part of the class
////////////////////////////////////////////////////

void ARM_GramFctorArg::CopyNoCleanUp( const ARM_GramFctorArg& rhs )
{
	itsType = rhs.itsType;
	
	switch( itsType )
	{
	case GFAT_DOUBLE_TYPE:
		{
			itsDouble = rhs.itsDouble;
			break;
		}
	case GFAT_STRING_TYPE:
		{
			itsString = rhs.itsString;
			break;
		}
	case GFAT_DATE_TYPE:
		{
			itsDate = rhs.itsDate;
			break;
		}
	case GFAT_VECTOR_TYPE:
		{
			itsVector = rhs.itsVector;
			break;
		}
	case GFAT_MATRIX_TYPE:
		{
			itsMatrix = rhs.itsMatrix;
			break;
		}
    case GFAT_CURVE_TYPE:
		{
			itsCurve = rhs.itsCurve;
			break;
		}
	case GFAT_STRINGVEC_TYPE:
		{
			itsStringVector = rhs.itsStringVector;
			break;
		}
	case GFAT_STRINGVECTRANS_TYPE:
		{
			itsStringVector = rhs.itsStringVector;
			break;
		}

	case GFAT_DATESTRIP_TYPE:
		{
			itsDateStrip = rhs.itsDateStrip;
			break;
		}
	default:
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				"Unknown type... please advise\n" ); 
		}
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArg
///	Routine: GetVector
///	Returns: vector pointor value
///	Action : Check and get vector pointor value
////////////////////////////////////////////////////

ARM_GramFctorArg::ARM_GramFctorArg( const ARM_GramFctorArg& rhs )
: itsType( rhs.itsType ) 
{
	CopyNoCleanUp( rhs );
}


////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArg
///	Routine: GetVector
///	Returns: vector pointor value
///	Action : Check and get vector pointor value
////////////////////////////////////////////////////

ARM_GramFctorArg& ARM_GramFctorArg::operator=( const ARM_GramFctorArg& rhs )
{
	if( this != &rhs )
	{
		CopyNoCleanUp( rhs );
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArg
///	Routine: GetVector
///	Returns: vector pointor value
///	Action : Check and get vector pointor value
////////////////////////////////////////////////////
ARM_GramFctorArg::~ARM_GramFctorArg()
{
	/// nothing because the pointor are given from the outside world! 
}




////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArg
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////

string ARM_GramFctorArg::toString() const
{
	CC_Ostringstream os;
	size_t i;
	os << TypeToString(itsType) << " (=";
	
	switch( itsType )
	{
	case GFAT_DOUBLE_TYPE:
		{
			os << itsDouble;
			break;
		}
	case GFAT_STRING_TYPE:
		{
			os << itsString;
			break;
		}
	case GFAT_DATE_TYPE:
		{
			os << itsDate.toString();
			break;
		}
	case GFAT_STRINGVEC_TYPE:
		{
			if( itsStringVector != ARM_StringVectorPtr(NULL) )
				if ( size_t vecSize = itsStringVector->size() )
				{
					const size_t truncatedSize = 50;
					size_t displaySize;
					bool truncated	= false;

					os << "size=" << vecSize;
				
					if( vecSize > truncatedSize )
					{
						displaySize = truncatedSize;
						truncated	= true;
					}
					else
						displaySize = vecSize;
	
					os << " [" ;
					if( displaySize )
					{
						for(i=0; i<displaySize-1; ++i )
							os << (*itsStringVector)[i] << ", ";
						os << (*itsStringVector)[i];
					}	

					if( truncated )
						os << " ... " << (*itsStringVector)[itsStringVector->size()-1];
					os << "]";
				}
				else
					os << "empty string vector";
		}
	case GFAT_STRINGVECTRANS_TYPE:
		{
			if( itsStringVector != ARM_StringVectorPtr(NULL) )
				if ( size_t vecSize = itsStringVector->size() )
				{
					const size_t truncatedSize = 50;
					size_t displaySize;
					bool truncated	= false;

					os << "size=" << vecSize;
				
					if( vecSize > truncatedSize )
					{
						displaySize = truncatedSize;
						truncated	= true;
					}
					else
						displaySize = vecSize;
	
					os << " [" ;
					if( displaySize )
					{
						for(i=0; i<displaySize-1; ++i )
							os << (*itsStringVector)[i] << ", ";
						os << (*itsStringVector)[i];
					}	

					if( truncated )
						os << " ... " << (*itsStringVector)[itsStringVector->size()-1];
					os << "]";
				}
				else
					os << "empty string vector";
		}
	case GFAT_VECTOR_TYPE:
		{
			if( itsVector != ARM_GP_VectorPtr(NULL) )
				if ( size_t vecSize = itsVector->size() )
				{
					const size_t truncatedSize = 50;
					size_t displaySize;
					bool truncated	= false;

					os << "size=" << vecSize;
				
					if( vecSize > truncatedSize )
					{
						displaySize = truncatedSize;
						truncated	= true;
					}
					else
						displaySize = vecSize;
	
					os << " [" ;
					if( displaySize )
					{
						for(i=0; i<displaySize-1; ++i )
							os << (*itsVector)[i] << ", ";
						os << (*itsVector)[i];
					}	

					if( truncated )
						os << " ... " << (*itsVector)[itsVector->size()-1];
					os << "]";
				}
				else
					os << "empty vector";
			break;
		}
	case GFAT_MATRIX_TYPE:
		{
			if ( itsMatrix != ARM_GP_MatrixPtr(NULL) )
			{
				size_t matNbRows = itsMatrix->GetRowsNb();
				size_t matNbCols = itsMatrix->GetColsNb();
				if (matNbRows && matNbCols)
				{
					const size_t truncatedSize = 50;
					size_t displaySize;
					bool truncated	= false;

					os << "nbRows=" << matNbRows;
					os << " nbCols=" << matNbCols;
				
					if( matNbCols > truncatedSize )
					{
						displaySize = truncatedSize;
						truncated	= true;
					}
					else
						displaySize = matNbCols;

					os << " [" ;
					if( displaySize )
					{
						for(i=0; i<displaySize-1; ++i )
							os << (*itsMatrix)(0,i) << ", ";
						os << (*itsMatrix)(0,i);
					}

					if( truncated )
						os << " ... " << (*itsMatrix)(0,matNbCols-1);
					os << "]";


					if (matNbRows > 1)
					{
						os << " ... ";

						os << " [" ;
						for(size_t i=0; i<displaySize; ++i )
						{
							os << (*itsMatrix)(matNbRows-1,i);
							if (i < matNbRows-1)
								os << ", ";
						}

						if( truncated )
							os << " ... " << (*itsMatrix)(matNbRows-1,matNbCols-1);
						os << "]";
					}
				}
			}
			break; 
		}

    case GFAT_CURVE_TYPE:
		{
			if( itsCurve != ARM_GP_CurvePtr(NULL) ) 
				if ( size_t curveSize = itsCurve->size() )
				{
					const size_t truncatedSize = 50;
					size_t displaySize;
					bool truncated	= false;
	
					os << "size=" << curveSize;
					
					if( curveSize > truncatedSize )
					{
							displaySize = truncatedSize;
						truncated	= true;
					}
					else
						displaySize = curveSize;
	
						os << " [" ;
					for(size_t i=0; i<displaySize; ++i )
					{
						os << " ";
						os << itsCurve->GetOrdinates()[i];
					}

					if( truncated )
						os << " ... " << itsCurve->GetOrdinates()[itsCurve->size()-1];
					os << "]";
				}
				else
					os << "empty curve";
			break;
		}

    case GFAT_DATESTRIP_TYPE:
		{
			if( itsDateStrip != ARM_DateStripPtr(NULL) ) 
			{
                os << itsDateStrip->toString();
            }
			else
				os << "empty date strip";

            break;
        }

	default:
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				"Unknown type... please advise\n" ); 
		}
	}

	os << ")\n";

	return os.str();
}



////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArg
///	Routine: Duplicate
///	Returns: an object similar to this one
///	Action : Clone 
////////////////////////////////////////////////////
ARM_GramFctorArg ARM_GramFctorArg::Duplicate() const
{
	switch( itsType )
	{
	case GFAT_DOUBLE_TYPE:
		{
			return ARM_GramFctorArg( itsDouble );
			break;
		}
	case GFAT_STRING_TYPE:
		{
			return ARM_GramFctorArg( itsString );
			break;
		}
	case GFAT_DATE_TYPE:
		{
			return ARM_GramFctorArg( itsDate );
			break;
		}
	/*case GFAT_VECTOR_TYPE:
		{
			return ARM_GramFctorArg( *itsVector );
			break;
		}*/
	case GFAT_STRINGVEC_TYPE:
		{
			return ARM_GramFctorArg( ARM_StringVectorPtr( new ARM_StringVector(*itsStringVector) ) );
			break;
		}
	case GFAT_STRINGVECTRANS_TYPE:
		{
			return ARM_GramFctorArg( ARM_StringVectorPtr( new ARM_StringVector(*itsStringVector) ) );
			break;
		}
	case GFAT_MATRIX_TYPE:
		{
			return ARM_GramFctorArg( ARM_GP_MatrixPtr( (ARM_GP_Matrix*) itsMatrix->Clone() ) );
			break;
		}
	default:
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				"Unknown type... please advise\n" ); 
		}
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArg
///	Routine: GetResultType
///	Returns: Type
///	Action : get the validation type 
////////////////////////////////////////////////////

GramFuncArgType ARM_GramFctorArg::GetResultType( GramFuncArgType t )
{
	switch( t )
	{

	case GFAT_DOUBLE_TYPE:			
		return GFAT_DOUBLE_TYPE;

	case GFAT_STRING_TYPE: case GFAT_MODEL_TYPE: case GFAT_MATURITY_TYPE: case GFAT_MULTITOKENSTRING_TYPE:
		return GFAT_STRING_TYPE;
	
	case GFAT_DATEORMATU_TYPE:
		return GFAT_DATEORMATU_TYPE;

	case GFAT_DATE_TYPE:			
		return GFAT_DATE_TYPE;

	case GFAT_DATE_OR_DOUBLE_TYPE:			
		return GFAT_DATE_OR_DOUBLE_TYPE;

	case GFAT_DATE_OR_VECTOR_TYPE:			
		return GFAT_DATE_OR_VECTOR_TYPE;

    case GFAT_VECTOR_TYPE: case GFAT_FUTUREREF_TYPE: case GFAT_VECTOR_OR_FUTUREREF_TYPE:	
		return GFAT_VECTOR_TYPE;

	case GFAT_MATRIX_TYPE:	
		return GFAT_MATRIX_TYPE;

	case GFAT_CURVE_TYPE:	
		return GFAT_CURVE_TYPE;
	
	case GFAT_VECTOR_OR_CURVE_TYPE:	
		return GFAT_VECTOR_OR_CURVE_TYPE;

	case GFAT_STRING_OR_CURVE_TYPE:	
		return GFAT_STRING_OR_CURVE_TYPE;

	case GFAT_DATESTRIP_TYPE:	
		return GFAT_DATESTRIP_TYPE;

	case GFAT_TERMINATOR: case GFAT_UNKNOWN_TYPE:
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"Should never EVER return the type " + TypeToString(t) );
	default:
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Unknonw type");

	}
}

////////////////////////////////////////////////////
///	Class  : ARM_GramFctorArg
///	Routine: Reset
///	Returns: void
///	Action : free all the memory of the gram functor!
////////////////////////////////////////////////////

void ARM_GramFctorArg::Reset()
{
	/// free the smart pointors!
	itsVector		= ARM_GP_VectorPtr(NULL);
	itsMatrix		= ARM_GP_MatrixPtr(NULL);
    itsCurve		= ARM_GP_CurvePtr(NULL);
	itsStringVector	= ARM_StringVectorPtr(NULL);
	itsString		= string("");
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

