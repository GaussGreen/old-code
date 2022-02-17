#include <ARM\libarm_local\firstToBeIncluded.h>
#include <functional>
#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_gp_pricer.h>
#include "ARM_xl_gp_gensecurity_local.h"
#include "ARM_xl_wrapper_local.h"
#include <GP_Infra\gpinfra\gramfunctorarg.h>
#include "ARM_xl_trycatch_local.h"
#include <gpbase\gpmatrix.h>

using ARM::ARM_GramFctorArg;
using ARM::GFAT_DOUBLE_TYPE;
using ARM::GFAT_VECTOR_TYPE;
using ARM::GFAT_MATRIX_TYPE;
using ARM::GFAT_STRING_TYPE;
using ARM::GFAT_STRINGVEC_TYPE;
using ARM::GFAT_STRINGVECTRANS_TYPE;
using ARM::GFAT_DATE_TYPE;

LPXLOPER ARM_GramFunctorToXLOPER(const ARM_GramFctorArg& gramFunctorArg, XLOPER& XL_result, ARM_result C_result, bool fromExcel,int index)
{
	switch( gramFunctorArg.GetType() )
	{
	case GFAT_DOUBLE_TYPE:
		{
			XL_result.xltype  = xltypeNum;
			XL_result.val.num = gramFunctorArg.GetDouble();
		}
		break;
	case GFAT_VECTOR_TYPE:
		{
			/// convert to a vector
			int sizeVector = gramFunctorArg.GetVector()->size();
			if (index>=sizeVector)
			{
				ARM_ERR();
			}
			sizeVector-= index;
			VECTOR<double> vectorResult(sizeVector);
			
			for( size_t i=0; i<sizeVector; ++i )
				vectorResult[i] = gramFunctorArg.GetVector()->Elt(i + index);
			
			/// add these additional lines 
			/// to display blank lines
			const int additionalLinesNb = 100;
			bool fillWithBlank = fromExcel;
			XL_writeNumVectorWithOptions( XL_result, vectorResult, " ARM_ERR: Could not get gramFunctorArg array data", C_result, additionalLinesNb, fillWithBlank );
		}
		break;
	case GFAT_STRINGVEC_TYPE:
		{
			/// convert to a vector
			int sizeVector = gramFunctorArg.GetStringVector()->size();
			if (index>=sizeVector)
			{
				ARM_ERR();
			}
			sizeVector-= index;
			VECTOR<CCString> vectorResult(sizeVector);
			char * tmp;

			for( size_t i=0; i<sizeVector; ++i ){
				tmp = const_cast<char*> ( (*gramFunctorArg.GetStringVector())[i + index].c_str()  );
				vectorResult[i] =  CCString( tmp  );
			
			}
			/// add these additional lines 
			/// to display blank lines
			const int additionalLinesNb = 100;
			bool fillWithBlank = fromExcel;
			XL_writeStrVectorWithOptions( XL_result, vectorResult, " ARM_ERR: Could not get gramFunctorArg array data", C_result, additionalLinesNb, fillWithBlank );
		}
		break;
	case GFAT_STRINGVECTRANS_TYPE:
		{
			/// convert to a vector
			int sizeVector = gramFunctorArg.GetStringVectorTrans()->size();
			if (index>=sizeVector)
			{
				ARM_ERR();
			}
			sizeVector-= index;
			VECTOR<CCString> vectorResult(sizeVector);
			char * tmp;

			for( size_t i=0; i<sizeVector; ++i ){
				tmp = const_cast<char*> ( (*gramFunctorArg.GetStringVectorTrans())[i + index].c_str()  );
				vectorResult[i] =  CCString( tmp  );
				}
			/// add these additional lines 
			/// to display blank lines
			const int additionalColsNb = 100;
			bool fillWithBlank = fromExcel;
			XL_writeStrLineVectorWithOptions( XL_result, vectorResult, " ARM_ERR: Could not get gramFunctorArg array data", C_result, additionalColsNb, fillWithBlank );
		}
		break;

	case GFAT_MATRIX_TYPE:
		{
			int nbRows = gramFunctorArg.GetMatrix()->GetRowsNb(), nbCols =  gramFunctorArg.GetMatrix()->GetColsNb();
			if (index>=nbCols)
			{
				ARM_ERR();
			}
			nbCols -= index;

			VECTOR<double> vectorResult(nbRows*nbCols);
			for ( size_t i=0; i < nbRows; ++i)
				for ( size_t j=0; j < nbCols; ++j)
					vectorResult[i*nbCols+j] = (*gramFunctorArg.GetMatrix())(i,index+j);

            /// add these additional lines 
			/// to display blank lines
			const int additionalLinesNb = 100;
			bool fillWithBlank = fromExcel;
			XL_writeNumMatrixSizeWithOptions( XL_result, vectorResult, nbRows, nbCols, " ARM_ERR: Could not set the num matrix", C_result,additionalLinesNb,fillWithBlank );
		}
		break;
	case GFAT_STRING_TYPE:
		{
			XL_result.xltype  = xltypeStr;
			XL_result.xltype |= xlbitDLLFree;
			XL_result.val.str = XL_StrC2StrPascal( gramFunctorArg.GetString().c_str() );
		}
		break;
	case GFAT_DATE_TYPE:
		{
			const double JULIANDATEADD	= 2415019.0;
			XL_result.xltype  = xltypeNum;
			XL_result.val.num = gramFunctorArg.GetDate().GetJulian() - JULIANDATEADD;
		}
		break;
	default:
		ARM_ERR();
	}

	return &XL_result;
}
