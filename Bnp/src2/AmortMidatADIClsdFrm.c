/*-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

	Description	: implementation of the Closed Form Solution of the European Amortizing/Accreting Swaption.

				(C)	2004 BNP Paribas.  All rights reserved.

	Author		:	 Albert Wang

	Created		:	02.12.2004

	Changes		:

-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-*/

#include "AmortMidatADIClsdFrm.h"
#include "AmortMidatADIUtils.h"
#include "CPDCalib.h"

#ifdef _DEBUG
#include "crtdbg.h"
#endif 


static void _GenMidAt_CashFlow_(
				// market
				const char *szYC, 
				long lToday,
				//call
				double dCoupon,
				//// fixed leg
				const long *plFixPay_Begin,
				const long *plFixPay_End,
				const double *pdFixCoverage_Begin,
				const double *pdFixNotional_Begin,
				//// floating leg
				long lFltStart,
				const long *plFltPay_Begin,
				const long *plFltPay_End,
				const double *pdFltCoverage_Begin,
				const double *pdFltNotional_Begin,
				const double *pdMargin_Begin,
				const double *pdSpread_Begin,
				/// pay or receive
				double dPayRec, // 1 - receive fixed; 
				// optional Param_Caplet
				const long *plExFeePay,
				const double *pdExFee,
				const long *plFixFeePay,
				const double *pdFixFee,
				// results
				double *pdCashFlow_Begin,
				double *pdT_CashFlow_Begin
				)
{
	const int nSize_FixPay = plFixPay_End-plFixPay_Begin;
	const int nSize_FltPay = plFltPay_End-plFltPay_Begin;
	const int nSize_FltStart = 1;
	const int nSize_FixFee = plFixFeePay?1:0;
	const int nSize_ExFee = plExFeePay?1:0;
	const double *pdCashFlow_End 
		= pdCashFlow_Begin+nSize_FixPay+nSize_FltPay+nSize_FltStart+nSize_FixFee+nSize_ExFee;
		
	// exercise fee
	if(plExFeePay&&pdExFee)
	{
		const double dDF = swp_f_df(lToday,*plExFeePay,szYC);
		*pdCashFlow_Begin = -dDF*(*pdExFee);// always subtract the exercise fee
		_fill_tenor(lToday,plExFeePay,1+plExFeePay,pdT_CashFlow_Begin);
		++pdCashFlow_Begin;
		++pdT_CashFlow_Begin;
	}

	// fixed leg
	{
		for(;plFixPay_Begin < plFixPay_End; 
		++pdCashFlow_Begin,++pdT_CashFlow_Begin,++plFixPay_Begin,++pdFixCoverage_Begin,++pdFixNotional_Begin
			)
		{
			const double dDF = swp_f_df(lToday,*plFixPay_Begin,szYC);
			*pdCashFlow_Begin = (*pdFixCoverage_Begin)*(*pdFixNotional_Begin)*dCoupon*dDF ;
			*pdCashFlow_Begin *= dPayRec;
			_fill_tenor(lToday,plFixPay_Begin,1+plFixPay_Begin,pdT_CashFlow_Begin);
		}
	}
	
	// fixed fee
	if(plFixFeePay&&pdFixFee)
	{
		const double dDF = swp_f_df(lToday,*plFixFeePay,szYC);
		*pdCashFlow_Begin -= dDF*(*pdFixFee)*dPayRec;
		_fill_tenor(lToday,plFixFeePay,1+plFixFeePay,pdT_CashFlow_Begin);
		++pdCashFlow_Begin;
		++pdT_CashFlow_Begin;
	}

	// floating leg
	{
		// floating start
		const double dDF = swp_f_df(lToday,lFltStart,szYC);
		*pdCashFlow_Begin = -dPayRec*dDF*(*pdFltNotional_Begin);
		_fill_tenor(lToday,&lFltStart,1+&lFltStart,pdT_CashFlow_Begin);
		++pdCashFlow_Begin;
		++pdT_CashFlow_Begin;
		
		// cash flows in between float start and the last float cash flow
		for(;plFltPay_Begin<plFltPay_End-1; // NB: not inclusing the last cash flow 
		++pdCashFlow_Begin,++pdT_CashFlow_Begin,++plFltPay_Begin,++pdFltCoverage_Begin,++pdFltNotional_Begin,++pdMargin_Begin,++pdSpread_Begin
		)
		{
			const double dDF = swp_f_df(lToday,*plFltPay_Begin,szYC);
			*pdCashFlow_Begin = *(pdFltNotional_Begin+1)- *pdFltNotional_Begin;
			*pdCashFlow_Begin+= (*pdFltNotional_Begin)*(*pdFltCoverage_Begin)*(*pdMargin_Begin+*pdSpread_Begin);
			*pdCashFlow_Begin*= -dDF*dPayRec;
			_fill_tenor(lToday,plFltPay_Begin,1+plFltPay_Begin,pdT_CashFlow_Begin);
		}

		// last floating cash flow
		*pdCashFlow_Begin = swp_f_df(lToday,*plFltPay_Begin,szYC);
		*pdCashFlow_Begin *= (*pdFltNotional_Begin);
		*pdCashFlow_Begin *= dPayRec;
		_fill_tenor(lToday,plFltPay_Begin,1+plFltPay_Begin,pdT_CashFlow_Begin);
	}
}

const char* GenMidAt_ClsdFrm_1F(
				// market
				const char *szYC, 
				long lToday,
				// model
				const double *pdSigT_Begin,
				const double *pdSigT_End,
				const double *pdSig_Begin,
				const double *pdLamT_Begin,
				const double *pdLamT_End,
				const double *pdLam_Begin,
				/// Call
				long lEx,
				double dCoupon,
				//// fixed leg
				const long *plFixPay_Begin,
				const long *plFixPay_End,
				const double *pdFixCoverage_Begin,
				const double *pdFixNotional_Begin,
				//// floating leg
				long lFltStart,
				const long *plFltPay_Begin,
				const long *plFltPay_End,
				const double *pdFltCoverage_Begin,
				const double *pdFltNotional_Begin,
				const double *pdMargin_Begin,
				const double *pdSpread_Begin,
				/// pay or receive
				double dPayRec,
				/// optional Param_Caplet
				const long *plExFeePay,
				const double *pdExFee,
				const long *plFixFeePay,
				const double *pdFixFee,
				/// result
				double *pdOptionPV
				)
{
	const int nSize_FixPay = plFixPay_End-plFixPay_Begin;
	const int nSize_FltPay = plFltPay_End-plFltPay_Begin;
	const int nSize_FltStart = 1;
	const int nSize_FixFee = plFixFeePay&&pdFixFee?1:0;
	const int nSize_ExFee = plExFeePay&&pdExFee?1:0;
	const int nSize_CashFlow = 
		nSize_FixPay+nSize_FltPay+nSize_FltStart+nSize_FixFee+nSize_ExFee;
	
	// cash flows
	double *pdCashFlow = _alloca(nSize_CashFlow*sizeof(double));
	double *pdT_CashFlow = _alloca(nSize_CashFlow*sizeof(double));

	// Zeta
	double dEx = 0.;
	double dZeta=0.;

	// G1 and G2 at cash flow dates
	double *pdG = _alloca(nSize_CashFlow*sizeof(double));
	double dG_Ex=0.;

	const char *szErr = 0;

	/// create cash flows for the GenMidAt
	_GenMidAt_CashFlow_(
		szYC, 
		lToday,
		dCoupon,
		plFixPay_Begin,
		plFixPay_End,
		pdFixCoverage_Begin,
		pdFixNotional_Begin,
		lFltStart,
		plFltPay_Begin,
		plFltPay_End,
		pdFltCoverage_Begin,
		pdFltNotional_Begin,
		pdMargin_Begin,
		pdSpread_Begin,
		dPayRec,
		plExFeePay,
		pdExFee,
		plFixFeePay,
		pdFixFee,
		pdCashFlow,
		pdT_CashFlow
		);

	/// Covert Sigma to Zeta
	_fill_tenor(lToday,&lEx,1+&lEx,&dEx);
	szErr = Convert_SigToZeta_1F(
			pdSigT_Begin,
			pdSigT_End,
			pdSig_Begin,
			pdLamT_Begin,
			pdLamT_End,
			pdLam_Begin,
			&dEx,
			1+&dEx,
			&dZeta
			);
	if(szErr) return szErr;

	/// Covert Lambda to G
	Convert_LamToG_1F(
		pdLamT_Begin,
		pdLamT_End,
		pdLam_Begin,
		pdT_CashFlow,
		pdT_CashFlow+nSize_CashFlow,
		pdG // G at cash flow dates
		);

	Convert_LamToG_1F(
		pdLamT_Begin,
		pdLamT_End,
		pdLam_Begin,
		&dEx,
		1+&dEx,
		&dG_Ex // G at exercise date
		);

	/// price 
	*pdOptionPV = 
		lgmopval1F(nSize_CashFlow,pdCashFlow,pdG,dZeta,dG_Ex);

	return 0;
}


const char* GenMidAt_ClsdFrm_2F(
				// market
				const char *szYC, 
				long lToday,
				// model
				const double *pdSigT_Begin,
				const double *pdSigT_End,
				const double *pdSig_Begin,
				const double *pdLamT_Begin,
				const double *pdLamT_End,
				const double *pdLam_Begin,
				double dAlpha,
				double dGamma,
				double dRho,
				/// Call
				long lEx,
				double dCoupon,
				//// fixed leg
				const long *plFixPay_Begin,
				const long *plFixPay_End,
				const double *pdFixCoverage_Begin,
				const double *pdFixNotional_Begin,
				//// floating leg
				long lFltStart,
				const long *plFltPay_Begin,
				const long *plFltPay_End,
				const double *pdFltCoverage_Begin,
				const double *pdFltNotional_Begin,
				const double *pdMargin_Begin,
				const double *pdSpread_Begin,
				/// pay or receive
				double dPayRec,
				const long *plExFeePay,
				const double *pdExFee,
				const long *plFixFeePay,
				const double *pdFixFee,
				/// result
				double *pdOptionPV
				)
{
	const int nSize_FixPay = plFixPay_End-plFixPay_Begin;
	const int nSize_FltPay = plFltPay_End-plFltPay_Begin;
	const int nSize_FltStart = 1;
	const int nSize_FixFee = plFixFeePay&&pdFixFee?1:0;
	const int nSize_ExFee = plExFeePay&&pdExFee?1:0;
	const int nSize_CashFlow 
		= nSize_FixPay+nSize_FltPay+nSize_FltStart+nSize_FixFee+nSize_ExFee;
	
	// cash flows
	double *pdCashFlow = _alloca(nSize_CashFlow*sizeof(double));
	double *pdT_CashFlow = _alloca(nSize_CashFlow*sizeof(double));

	// Zeta
	double dEx = 0.;
	double dZeta1=0., dZeta2=0.,dZeta12=0.;

	// G1 and G2 at cash flow dates
	double *pdG1 = _alloca(nSize_CashFlow*sizeof(double));
	double *pdG2 = _alloca(nSize_CashFlow*sizeof(double));
	double dG1_Ex=0.,dG2_Ex=0.;

	const char *szErr = 0;

	/// create cash flows for the GenMidAt
	_GenMidAt_CashFlow_(
		szYC, 
		lToday,
		dCoupon,
		plFixPay_Begin,
		plFixPay_End,
		pdFixCoverage_Begin,
		pdFixNotional_Begin,
		lFltStart,
		plFltPay_Begin,
		plFltPay_End,
		pdFltCoverage_Begin,
		pdFltNotional_Begin,
		pdMargin_Begin,
		pdSpread_Begin,
		dPayRec,
		plExFeePay,
		pdExFee,
		plFixFeePay,
		pdFixFee,
		pdCashFlow,
		pdT_CashFlow
		);


	/// Covert Sigma to Zeta
	_fill_tenor(lToday,&lEx,1+&lEx,&dEx);
	szErr = Convert_SigToZeta_2F(
			pdSigT_Begin,
			pdSigT_End,
			pdSig_Begin,
			pdLamT_Begin,
			pdLamT_End,
			pdLam_Begin,
			&dEx,
			1+&dEx,
			dAlpha,
			dGamma,
			dRho,
			&dZeta1,
			&dZeta2,
			&dZeta12
			);
	if(szErr) return szErr;

	/// Covert Lambda to G
	Convert_LamToG_2F(
		pdLamT_Begin,
		pdLamT_End,
		pdLam_Begin,
		pdT_CashFlow,
		pdT_CashFlow+nSize_CashFlow,
		dGamma,
		pdG1, // G1 at cash flow dates
		pdG2  // G2 at cash flow dates
		);

	Convert_LamToG_2F(
		pdLamT_Begin,
		pdLamT_End,
		pdLam_Begin,
		&dEx,
		1+&dEx,
		dGamma,
		&dG1_Ex, // G1 at exercise date
		&dG2_Ex  // G2 at exercise date
		);

	/// price 
	*pdOptionPV = 
		lgmopval2F(nSize_CashFlow,pdCashFlow,pdG1,pdG2,dZeta1,dZeta2,dZeta12,dG1_Ex,dG2_Ex);

	return 0;
}

void Convert_LamToG_1F(
			const double *pdLamT_Begin,
			const double *pdLamT_End,
			const double *pdLam_Begin,
			const double *pdGT_Begin,
			const double *pdGT_End,
			// results
			double *pdG_Begin
			)
{
	_ASSERTE(pdLamT_End-pdLamT_Begin>=1);
	export_lgmsetupG_ts(
		pdLamT_End-pdLamT_Begin,
		(double*)pdLamT_Begin,
		(double*)pdLam_Begin,
		pdGT_End-pdGT_Begin,
		(double*)pdGT_Begin,
		pdG_Begin,
		0,
		0,
		0);

	
#ifdef _DEBUG
	// Check again export_lgmsetupG if there is a single lambda
	if(	pdLamT_End-pdLamT_Begin ==1)
	{
		double *pdG = _alloca((pdGT_End-pdGT_Begin)*sizeof(double));
		export_lgmsetupG(
			*pdLam_Begin,
			pdGT_End-pdGT_Begin,
			(double*)pdGT_Begin,
			pdG,
			0,
			0,
			0);
		
		for(;pdGT_Begin<pdGT_End;++pdGT_Begin,++pdG,++pdG_Begin)
			_ASSERTE(fabs(*pdG-*pdG_Begin)<1.e-15);
	}
#endif //_DEBUG
}

void Convert_LamToG_2F(
			const double *pdLamT_Begin,
			const double *pdLamT_End,
			const double *pdLam_Begin,
			const double *pdGT_Begin,
			const double *pdGT_End,
			double dGamma,
			// results
			double *pdG1_Begin,
			double *pdG2_Begin
			)
{
	_ASSERTE(pdLamT_End-pdLamT_Begin>=1);
	
	// fill in pdG1_Begin
	Convert_LamToG_1F(pdLamT_Begin,	pdLamT_End,	pdLam_Begin,pdGT_Begin,	pdGT_End,pdG1_Begin);
	export_lgmsetupG2_ts(
		pdLamT_End-pdLamT_Begin,
		(double*)pdLamT_Begin,
		(double*)pdLam_Begin,
		dGamma,
		pdGT_End-pdGT_Begin,
		(double*)pdGT_Begin,
		pdG2_Begin,
		0,
		0,
		0);

#ifdef _DEBUG
	// Check again export_lgmsetupG2 if there is a single lambda
	if(	pdLamT_End-pdLamT_Begin ==1)
	{
		double *pdG2 = _alloca((pdGT_End-pdGT_Begin)*sizeof(double));
		export_lgmsetupG2(
			*pdLam_Begin,
			dGamma,
			pdGT_End-pdGT_Begin,
			(double*)pdGT_Begin,
			pdG2,
			0,
			0,
			0);

		for(;pdGT_Begin<pdGT_End;++pdGT_Begin,++pdG2,++pdG2_Begin)
			_ASSERTE(fabs(*pdG2-*pdG2_Begin)<1.e-15);
	}
#endif //_DEBUG

}

static const char *_Convert_SigToZeta_(
			const double *pdSigT_Begin,
			const double *pdSigT_End,
			const double *pdSig_Begin,
			const double *pdLamT_Begin,
			const double *pdLamT_End,
			const double *pdLam_Begin,
			const double *pdZetaT_Begin,
			const double *pdZetaT_End,
			double dGamma,
			// results
			double *pdZeta_Begin
			)
{
	char *computeZeta_ts(double,double*,double *,int,double,double*,double*,int,double*);
	char *computeZeta(double, double, double *, double *, int, double*);

	const double *pdZetaT = pdZetaT_Begin;
	double *pdZeta = pdZeta_Begin;
	for(;pdZetaT<pdZetaT_End;++pdZetaT,++pdZeta)
	{
		const char *szErr = computeZeta_ts(
			*pdZetaT,
			(double*)pdLamT_Begin,
			(double*)pdLam_Begin,
			pdLamT_End-pdLamT_Begin,
			dGamma,
			(double*)pdSigT_Begin,
			(double*)pdSig_Begin,
			pdSigT_End-pdSigT_Begin,
			pdZeta
			);

		if(szErr) return szErr;
	}

#ifdef _DEBUG
	// Check again computeZeta if there is a single lambda
	if(pdLamT_End-pdLamT_Begin==1)
	{
		double *pdZeta = _alloca((pdSigT_End-pdSigT_Begin)*sizeof(double));
		const double *pdZetaT = pdZetaT_Begin;		
		for(;pdZetaT<pdZetaT_End;++pdZeta,++pdZetaT,++pdZeta_Begin)
		{
			const char *szErr = computeZeta(
				*pdZetaT,
				*pdLam_Begin,
				(double*)pdSigT_Begin,
				(double*)pdSig_Begin,
				pdSigT_End-pdSigT_Begin,
				pdZeta
				);

			if(szErr) return szErr;
			_ASSERTE(fabs(*pdZeta-*pdZeta_Begin)<1.e-15);
		}	
	}
#endif //_DEBUG

	return 0;
}

const char *Convert_SigToZeta_1F(
			const double *pdSigT_Begin,
			const double *pdSigT_End,
			const double *pdSig_Begin,
			const double *pdLamT_Begin,
			const double *pdLamT_End,
			const double *pdLam_Begin,
			const double *pdZetaT_Begin,
			const double *pdZetaT_End,
			// results
			double *pdZeta_Begin
			)
{
	const double dGamma = 0.;
	return _Convert_SigToZeta_(
			pdSigT_Begin,
			pdSigT_End,
			pdSig_Begin,
			pdLamT_Begin,
			pdLamT_End,
			pdLam_Begin,
			pdZetaT_Begin,
			pdZetaT_End,
			dGamma,
			pdZeta_Begin
			);
}

const char *Convert_SigToZeta_2F(
			const double *pdSigT_Begin,
			const double *pdSigT_End,
			const double *pdSig_Begin,
			const double *pdLamT_Begin,
			const double *pdLamT_End,
			const double *pdLam_Begin,
			const double *pdZetaT_Begin,
			const double *pdZetaT_End,
			double dAlpha,
			double dGamma,
			double dRho,
			// results
			double *pdZeta1_Begin,
			double *pdZeta2_Begin,
			double *pdZeta12_Begin
			)
{
	const char *szErr = 0;
	
	/// fill Zeta1
	szErr = _Convert_SigToZeta_(
			pdSigT_Begin,
			pdSigT_End,
			pdSig_Begin,
			pdLamT_Begin,
			pdLamT_End,
			pdLam_Begin,
			pdZetaT_Begin,
			pdZetaT_End,
			0.,
			pdZeta1_Begin
			);
	if(szErr) return szErr;
	
	/// fill Zeta2
	szErr = _Convert_SigToZeta_(
			pdSigT_Begin,
			pdSigT_End,
			pdSig_Begin,
			pdLamT_Begin,
			pdLamT_End,
			pdLam_Begin,
			pdZetaT_Begin,
			pdZetaT_End,
			dGamma,
			pdZeta2_Begin
			);
	if(szErr) return szErr;
	
	_multiply_vector2( // multiply by dAlpha^2
		pdZeta2_Begin,
		pdZeta2_Begin+(pdZetaT_End-pdZetaT_Begin),
		dAlpha*dAlpha,
		pdZeta2_Begin
		);

	/// fill Zeta12
	szErr = _Convert_SigToZeta_(
			pdSigT_Begin,
			pdSigT_End,
			pdSig_Begin,
			pdLamT_Begin,
			pdLamT_End,
			pdLam_Begin,
			pdZetaT_Begin,
			pdZetaT_End,
			0.5*dGamma,
			pdZeta12_Begin
			);
	if(szErr) return szErr;
	
	_multiply_vector2( // multiply by dAlpha*dRho
		pdZeta12_Begin,
		pdZeta12_Begin+(pdZetaT_End-pdZetaT_Begin),
		dRho*dAlpha,
		pdZeta12_Begin
		);

	return 0;
}




