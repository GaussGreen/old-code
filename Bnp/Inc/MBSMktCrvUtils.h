#ifndef MBSMKTCRVUTILS_H
#define MBSMKTCRVUTILS_H

#include "Sort\inc\srt_h_types.h"
///////
	// C function
	#if defined __cplusplus
		#define _MBS_EXTERN_C extern "C"
	#else
		#define _MBS_EXTERN_C extern
	#endif

///
/////////////////////
//Vol Functions

_MBS_EXTERN_C
double MBS_GetVol(void * VC, //VolatilityCurve
				  long start,//(today + 2b) + term
				  long end,//start + tenor
				  double Strike//.05 for 5 %
				  );

_MBS_EXTERN_C
double MBS_GetVolTermTenor(void * VC, //VolatilityCurve
						   long today,
						   int TermNumUnits,SrtUnit TermUnitType,
						   int TenorNumUnits,SrtUnit TenorUnitType,
						   double Strike);

	/////////////////////////////////////////////////////////
	// swap rate 'functors' for C
	// NB:
	//	* functor represention by unique opaque handles
	//	* usage:
	//		+ SwapRate_Create: create handle for swap rate corresponding to given {yield curve, ref rate, basis, frequency}
	//		+ SwapRate_Evaluate: evaluate swaprate as bivariate function;
	//		+ SwapRate_Delete: delete handle when done

	_MBS_EXTERN_C 
	unsigned MBS_SwapRate_Create(	// NB: returns 0 on error
		SrtBasisCode eBasis_Fixed,
		SrtCompounding eFreq_Fixed,
		const char * szYC,
		const char * szRefRate
	);

	_MBS_EXTERN_C
	unsigned MBS_SwapRate_CreateFromYC(
				void * YC//IYieldCurve
				);

	_MBS_EXTERN_C 
	void MBS_SwapRate_Delete(
		unsigned hSwapRate	// handle to swaprate functor
	);

	_MBS_EXTERN_C 
	double MBS_SwapRate_Evaluate(	// NB: returns -1 on error
		unsigned hSwapRate,	// handle to swaprate functor
		long nDateStart,	// actual start date
		long nDateEnd	// theoretical end date
	);


	/////////////////////////////////////////////////////////
	// forward rate 'functors' for C
	// NB:
	//	* functor represention by unique opaque handles
	//	* usage:
	//		+ ForwardRate_Create: create handle for swap rate corresponding to given {yield curve, ref rate, basis, frequency}
	//		+ ForwardRate_Evaluate: evaluate swaprate as bivariate function;
	//		+ ForwardRate_Delete: delete handle when done

	_MBS_EXTERN_C 
	unsigned MBS_ForwardRate_Create(	// NB: returns 0 on error
		const char * szYC,
		const char * szRefRate
	);

	_MBS_EXTERN_C 
	void MBS_ForwardRate_Delete(
		unsigned hForwardRate	// handle to swaprate functor
	);

	_MBS_EXTERN_C 
	double MBS_ForwardRate_Evaluate(	// NB: returns -1 on error
		unsigned hForwardRate,	// handle to swaprate functor
		long nDateStart,	// actual start date
		long nDateEnd	// actual end date
	);

	_MBS_EXTERN_C
	double MBS_SwapRateEvaluate_TermTenor(
				void * YC_in,//IYieldCurve
				int TermNumUnits,SORT::SrtUnit TermUnitType,
				int TenorNumUnits,SORT::SrtUnit TenorUnitType
	);

#endif//MBSMKTCRVUTILS_H
