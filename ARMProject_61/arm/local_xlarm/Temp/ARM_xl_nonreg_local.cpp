#pragma warning(disable : 4786)
#pragma warning(disable : 4005)

#include <libCCxll\CCxll.h>


#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"
#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_util.h>
#include "ARM_xl_trycatch_local.h"

#include "util\tech_macro.h"

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_NumericalRegression(LPXLOPER XL_refValue,
																	LPXLOPER XL_newValue,
																	LPXLOPER XL_tolerance,
																	LPXLOPER XL_epsilon)
{
	ADD_LOG("ARM_NumericalRegression");
//	ARM_BEGIN();
	
	// return
	static XLOPER XL_result;
	LPXLOPER pxArray;

	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_refValue;
	double C_newValue;
	double C_tolerance;
	double C_tolerance_default = 10E-8;
	double C_epsilon;
	double C_epsilon_default = 10E-14;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(	XL_refValue,	C_refValue," ARM_ERR: refValue: numeric expected", C_result);
	XL_readNumCell(	XL_newValue,	C_newValue," ARM_ERR: newValue: numeric expected", C_result);
	XL_readNumCellWD(XL_tolerance,	C_tolerance,C_tolerance_default," ARM_ERR: tolerance: numeric expected",C_result);
	XL_readNumCellWD(XL_epsilon,	C_epsilon,	C_epsilon_default,	" ARM_ERR: epsilon: numeric expected",C_result);

	// gap calculation
	double gap = 0.0;
	if (C_refValue > C_epsilon)
	{
		gap = fabs( (C_refValue - C_newValue) / C_refValue);
	}
	else
	{
		gap = fabs( (C_refValue - C_newValue) / C_epsilon);
	}
	
	// regression status : 1 (OK no regression), 0 (KO regression)
	int res = 0;
	if (gap > C_tolerance)
		res = 0;
	else
		res = 1;

	int nbrows = 1;
	int nbcolumns = 1;
	
	// excel result construction
	FreeCurCellErr ();
	XL_result.xltype = xltypeMulti;
	XL_result.val.array.columns = nbcolumns;
	XL_result.val.array.rows = nbrows; 
	XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

	pxArray[XL_Coordonnate2Rank (0, 0, 1)].xltype = xltypeNum;
	pxArray[XL_Coordonnate2Rank (0, 0, 1)].val.num = res; 
	
//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetCorrelInst" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
