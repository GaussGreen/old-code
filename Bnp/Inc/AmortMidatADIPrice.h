// prevent multiple inclusions
#pragma once

//////////////////////
//	warnings
#pragma warning(disable : 4786)	//"identifier was truncated to '255' characters in the debug information"
// NB: force warnings for unused arguments and local parameters
#pragma warning (1 : 4100 4101)
#include "AmortMidatADI.h"
#include "AmortMidatADIModel.h"


typedef struct _Precomputed_Quantities_TimeSequence_
{
	const double *m_pdT_Begin; /// time sequence in year fraction
	const double *m_pdT_End; 

	const double *m_pdDF_Begin; /// discount factors observed at today = 0 for the time sequence
	const double *m_pdLAMZero1_Begin; // LAMBDAZERO1 for the time sequence
	const double *m_pdLAMZero2_Begin;// LAMBDAZERO1 for the time sequence
} _PCQ_Seq;

typedef struct _Precomputed_Quantities_Exercise_
{
	const _PCQ_Seq	*m_pSeq; /// info related to exercise date sequence
	const _PCQ_Seq	*m_pExFeeSeq; /// info related to exercise date sequence
	
	const double *m_pdEIL1_Begin; /// EIL a.k.a. exp(integration of lambda)
	const double *m_pdEIL2_Begin;
	const double *m_pdPHI1_Begin;
	const double *m_pdPHI2_Begin;
	const double *m_pdPHI12_Begin;
	//const double *m_pdCoupon_Begin;
	const double *m_pdExFee_Begin;
} _Ex;

typedef struct _Precomputed_Quantities_FixedLegPay_
{
	const _PCQ_Seq	*m_pSeq;
	const double *const m_pdCoverage_Begin;
	const double *const m_pdNotlProdCvg_Begin;
	const _PCQ_Seq	*m_pFeeSeq;
	const double *const m_pFixFee_Begin; // fixed fee 
	const double *m_pdFixCoupon_Begin;

} _Fix;

typedef struct _Precomputed_Quantities_FloatingLegPay_
{
	const _PCQ_Seq	*m_pSeq;
	const double *const m_pdCoverage_Begin;
	const double *const m_pdPaymentAdj_Begin; // notional * coverage * (spread + margin)
	const double *const m_pdCouponProdCvg_Begin; // coverage * (spread + margin)
	const double *const m_pdNotional_Begin;
} _Flt;

typedef struct _Generic_MidAt_Product_
{
		const _Ex *m_pEx;
		const _PCQ_Seq *m_pExStart;
		
		const _Fix *m_pFixPay;
		const double *m_pdFixStart_Begin;
		const double *m_pdFixEnd_Begin;
		
		const _Flt *m_pFltPay;
		const double *m_pdFltStart_Begin;
		const double *m_pdFltEnd_Begin;

		const double m_dPayRec;

} _GenMidAt;



///  const long * for term structure
const char *Price_GenMidAt(
	///////Market	
	long lToday,
	const char *szYC,
	/// MidAt Specs
	const char *szPayRec, //"REC" or "PAY"
	// MidAt Exercise 
	const long *plEx_Begin,// exercise dates
	const long *plEx_End, 
	const long *plExStart_Begin, // exercise premium dats
	const double *pdExFee_Begin,// exercise fee
	// MidAt Fixed Leg 
	const long* plFixPay_Begin,
	const long* plFixPay_End,
	const long *plFixStart_Begin,
	const long* plFixEnd_Begin,
	const double *pdFixCvg_Begin,
	const double *pdFixNotional_Begin,
	const double *pdFixCoupon_Begin,// coupon
	const double *pdFixFee_Begin,
	// MidAt Funding/Floating leg 
	const long* plFltPay_Begin,
	const long* plFltPay_End,
	const long* plFltStart_Begin,
	const long* plFltEnd_Begin,
	const double *pdFltCvg_Begin,
	const double *pdMargin_Begin,
	const double *pdSpread_Begin,
	const double *pdFltNotional_Begin,
	// term structure
	const long *plSigDate_Begin,
	const long *plSigDate_End,
	const double *pdSig_Begin,
	const long *plLamDate_Begin,
	const long *plLamDate_End,
	const double *pdLam_Begin,
	const double *pdAlpha,
	const double *pdGamma,
	const double *pdRho,
	// grid
	int nNumPointT,
	int nNumPointX,
	/// results,
	double *pdOptionPV,
	double *pdFixLegPV,
	double *pdFltLegPV,
	double *pdExProb,
	double *pdExBoundary,
	int nUse_Backward
	);

//// const double * for term structure
const char* _Price_GenMidAt(
				/// grid
				int nNumPointT,
				int nNumPointX,
				// market
				const char *szYC, // yield curve id
				long lToday,
				// model
				const double *pdSigT_Begin,
				const double *pdSigT_End,
				const double *pdSig_Begin,
				const double *pdLamT_Begin,
				const double *pdLamT_End,
				const double *pdLam_Begin,
				const double *pdAlpha,
				const double *pdGamma,
				const double *pdRho,
				/// Call
				const long *plEx_Begin,
				const long *plEx_End,
				const long *plExFee_Begin,
				const double *pdExFee_Begin,
				//// Ex Start
				const long *plExStart_Begin,
				//// fixed leg
				const long *plFixPay_Begin,
				const long *plFixPay_End,
				const long *plFixStart_Begin,
				const long *plFixEnd_Begin,
				const long *plFixFee_Begin,
				const double *pdFixCoverage_Begin,
				const double *pdFixNotional_Begin,
				const double *pdFixCoupon_Begin,
				const double *pdFixFee_Begin,
				//// floating leg
				const long *plFltPay_Begin,
				const long *plFltPay_End,
				const long *plFltStart_Begin,
				const long *plFltEnd_Begin,
				const double *pdFltCoverage_Begin,
				const double *pdFltNotional_Begin,
				const double *pdMargin_Begin,
				const double *pdSpread_Begin,
				/// pay or receive
				double dPayRec,
				/// result
				double *pdFixPV,
				double *pdFltPV,
				double *pOptionPV,
				double *pdExProb,
				double *pdExBoundary,
				int nUse_Backward
				);

const char *_GenMidAt_OptionPV_Core(
				/// grid
				int nNumPointT,
				int nNumPointX,
				// option
				const _GenMidAt *pOption,
				// model
				const _Model *pModel,
				/// numeraire
				const _PCQ_Seq *pNumeraire, // 0 to use jumping numeraire
				/// result
				double *pOptionPV
				);

const char* _alloc_Option(
						int nExSize,
						int nFixSize,
						int nFltSize,
						/// output
						_GenMidAt *pOption
						);

void _preinit_Option(
			 // market
			const char *szYC,
			long lToday,
			/// Ex
			const long *plEx_Begin,
			const long *plEx_End,
			const long *plExFee_Begin,
			const double *pdExFee_Begin,
			/// Ex Start
			const long *plExStart_Begin,
			//// fixed leg
			const long *plFixPay_Begin,
			const long *plFixPay_End,
			const long *plFixStart_Begin,
			const long *plFixEnd_Begin,
			const long *plFixFee_Begin,
			const double *pdFixCoverage_Begin,
			const double *pdFixNotional_Begin,
			const double *pdFixCoupon_Begin,
			const double *pdFixFee_Begin,
			//// floating leg
			const long *plFltPay_Begin,
			const long *plFltPay_End,
			const long *plFltStart_Begin,
			const long *plFltEnd_Begin,
			const double *pdFltCoverage_Begin,
			const double *pdFltNotional_Begin,
			const double *pdMargin_Begin,
			const double *pdSpread_Begin,
			// pay rec
			double dPayRec,
			/// output
			_GenMidAt *pOption
			);

const char *_alloc_Seq(
					int nDateSize,
					/// output
					_PCQ_Seq *pSeq,
					const double **ppdEIL1_Begin,
					const double **ppdEIL2_Begin
					);
void _preinit_Seq(
					 const char *szYC,
					long lToday,
					const long *plDate_Begin,
					const long *plDate_End,
					/// output
					const _PCQ_Seq *pSeq
					);

void _free_Seq(
		  const _PCQ_Seq* pSeq
		  );

const char * _GenMidAt_OptionPV(
		const _GenMidAt *pOption,
		const _Model *pModel,
		const _PCQ_Seq *pNumeraire,
		int nNumPointX,
		int nNumPointT,
		double *pOptionPV,
		double *pdExProb,
		double *pdExBoundary,
		int nUse_Backward
		);


void _free_Option(_GenMidAt *pOption);
