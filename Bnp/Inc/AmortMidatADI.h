// prevent multiple inclusions
#pragma once

//////////////////////
//	warnings
#pragma warning(disable : 4786)	//"identifier was truncated to '255' characters in the debug information"

// NB: force warnings for unused arguments and local parameters
#pragma warning (1 : 4100 4101)

#include "AmortMidatADIUtils.h"

/// signature of the coefficient functors needed by the 3D PDE 
typedef double (__stdcall * PFN_COEFF_3DPDE)(const void* pX1, const void *pX2,const void *pT,void*pComm);

void __stdcall _MULTD_PFN_BI(
	const double* pIn1,
	const double* pIn2, 
	const void *pComm,
	double *pProduct
	);


typedef struct _ADI_MultipleSlices_Storage
{
	// size at least (6*nSize_AuxX+5*nSize_X)*sizeof(double*)
	void **m_pm; // matrix
	// size at least 11*nSize_AuxX*nSize_X*sizeof(double))
	void *m_pdata;
} _ADI_Storage;

//-------------------------------------------------------------------------------------------------------
//	Description:  initializes the storage needed by the call to _ADI_MultipleSlices
//
//	Returns : storage
//-------------------------------------------------------------------------------------------------------
const char* Init_ADIStorage(_ADI_Storage *pstore,int nSize_X, int nSize_AuxX);

//-------------------------------------------------------------------------------------------------------
//	Description:  free storage
//
//	Returns :
//-------------------------------------------------------------------------------------------------------
void Free_ADIStorage(const _ADI_Storage *pStorage);

//-------------------------------------------------------------------------------------------------------
//	Description:   core PDE solving engine 
//						 solver traverses from pdT_Begin to pdT_End
//						  this version without mem storage in the interface
//	
//	Returns : Grid ppdU
//-------------------------------------------------------------------------------------------------------
const char *ADI_MultipleSlices(
		/// slice from pdT_Begin to pdT_End 
		/// IMPORTANT: by default backwards equation is solve
		/// pdT_Begin is a DESCENDING sequence of time knots
		const double *pdT_Begin,
		const double *pdT_End,
		// orthogonal states 
		const double *pdX_Begin, // dominant  state
		const double *pdX_End, 
		const double *pdAuxX_Begin, // secondary state
		const double *pdAuxX_End,		
		/// common parameters for updating vol and drift
		void *pFuncArg,
		// functor to update vol of the dominant variable
		// functor must accept the argument seqence (X,AuxX,T)
		PFN_COEFF_3DPDE _SIGMAX_Sqrd_FUNC_,
		// functor to update drift of the dominant variable
		// functor must accept the argument seqence (X,AuxX,T)
		PFN_COEFF_3DPDE _MUX_FUNC_,
		// functor to update vol of the secondary variable
		// functor must accept the argument seqence (AuxX,X,T)
		PFN_COEFF_3DPDE _SIGMAAuxX_Sqrd_FUNC_,
		// functor to update drift of the secondary variable
		// functor must accept the argument seqence (AuxX,X,T)
		PFN_COEFF_3DPDE _MUAuxX_FUNC_,
	    // On input :  Grid at time pdT_Begin
		// On output: Grid at time PdT_End-1
		// Row Indices = [nX_Begin,...] , Col Indices = [nAuxX_Begin, ...]
		double **ppdU,
		int nX_Begin,
		int nAuxX_Begin
		);


//-------------------------------------------------------------------------------------------------------
//	Description:   core PDE solving engine 
//						 solver traverses from pdT_Begin to pdT_End
//						  this version with mem storage in the interface
//	
//	Returns : Grid ppdU
//-------------------------------------------------------------------------------------------------------
void _ADI_MultipleSlices(
		/// slice from pdT_Begin to pdT_End 
		/// IMPORTANT: by default backwards equation is solve
		/// pdT_Begin is a DESCENDING sequence of time knots
		const double *pdT_Begin,
		const double *pdT_End,
		// orthogonal states 
		const double *pdX_Begin, // dominant  state
		const double *pdX_End, 
		const double *pdAuxX_Begin, // secondary state
		const double *pdAuxX_End,		
		/// common parameters for updating vol and drift
		void *pFuncArg,
		// functor to update vol of the dominant variable
		// functor must accept the argument seqence (X,AuxX,T)
		PFN_COEFF_3DPDE _SIGMAX_Sqrd_FUNC_,
		// functor to update drift of the dominant variable
		// functor must accept the argument seqence (X,AuxX,T)
		PFN_COEFF_3DPDE _MUX_FUNC_,
		// functor to update vol of the secondary variable
		// functor must accept the argument seqence (AuxX,X,T)
		PFN_COEFF_3DPDE _SIGMAAuxX_Sqrd_FUNC_,
		// functor to update drift of the secondary variable
		// functor must accept the argument seqence (AuxX,X,T)
		PFN_COEFF_3DPDE _MUAuxX_FUNC_,
	    // On input :  Grid at time pdT_Begin
		// On output: Grid at time PdT_End-1
		// Row Indices = [nX_Begin,...] , Col Indices = [nAuxX_Begin, ...]
		double **ppdU,
		int nX_Begin,
		int nAuxX_Begin,
		const _ADI_Storage *pStore
		);


#ifdef _DEBUG
	void _adi_tester_();
#endif 