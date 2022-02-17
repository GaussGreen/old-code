/*-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

	Description	: implementation of the calibration of Callable Amortizing/Accreting Swaption (CAS).

				(C)	2004 BNP Paribas.  All rights reserved.

	Author		:	 Albert Wang

	Created		:	02.12.2004

	Changes		:

-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-*/
#include "AmortMidatADICalib.h"
#include "AmortMidatADIUtils.h"
#include "AmortEuroCalib.h"
#include "AmortEuroPrice.h"
#include "DiagCalibDLM.h"
#include "AmortMidatCalib.h"
#include "AmortMidatADIAutoCal.h"

#if defined SIGN
#undef SIGN	// NB: SORT definition collides with Nag's
#endif //defined SIGN

#include "nag.h"
#include "nage04.h"

#ifdef _DEBUG
#include "crtdbg.h"
#endif 

typedef struct _Calibrate_LamTS_Common_
{
	const int m_nNumPointT;
	const int m_nNumPointX;
	const _GenMidAt *const m_pOption_Begin;
	const _GenMidAt *const m_pOption_End;
	const double *const m_pOptionPrice_Begin;// market price
	const _PCQ_Seq *m_pNumeraire_Begin;
	double *m_pdf;
	const  _Calib_LamTS_DiagEuro_Comm* m_pCalib_LamTS_DiagEuro_Begin;
	const double *m_pdAlpha;
	const double *m_pdGamma;
	const double *m_pdRho;
	double *m_pdLamT_Begin;
	double *m_pdLamT_End;
	const char *m_szErr;
	double *m_pdx;
} _Calib_LamTS_Comm;


typedef struct _GenEuropean_ // NB: const pointer just to be sure
{
	// exercise time
	const double m_dEx;
	// cash flow time
	const double *const m_pdT_CashFlow_Begin;
	const double *const m_pdT_CashFlow_End;

	// cash flows size
	const double *const m_pdCashFlow_Begin; 

	// term structures
	const double* const m_pdZeta_Ex;
	const double* const m_pdG_Ex;
	const double* const m_pdG_CashFlow_Begin;

	// market/target price
	const double m_dMktPrice;

} _GenEuropean;


typedef struct _Calib_GenEuropean_Arg_
{
	// vector of instruments to be calibrated
	const _GenEuropean** m_ppInst_Begin;
	
	// global term structure 
	double *const m_pdTS_Begin;
	const double *const m_pdTS_End;

	// beginning positions of Zeta and G on the pdTS
	const int m_nZeta_Begin;
	const int m_nG_Begin;

} _GenEuropeanArg;


typedef struct _Calib_GenMidAt_Arg_
{
	double a;
} _GenMidAt_Common;

void _init_Calib_LamTS_DiagEuro_Comm(
		long lToday,
		const char *szYC,
		const char * szVC,
		Err (*get_cash_vol)(char*,double,double,double,int,char*,double*,double*),
		long lTheoEnd,
		const char * szRefRate,
		const char * szFreq,
		const char * szBasis,
		const long *plEx_Begin,
		const long *plEx_End,
		const int nOnefequi,
		const double *pdCoupon_Begin,
		const long *plFixStart_Begin,
		const long *plFixStart_End,
		const long *plFixEnd_Begin,
		const long* plFixPay_Begin,
		const long *plFltStart_Begin,
		const long *plFltStart_End,
		const long *plFltEnd_Begin,
		const long* plFltPay_Begin,
		const double *pdMargin_Begin,
		const double *pdSpread_Begin,
		const cpd_diag_calib_param *pParam_Swaption,
		const cpd_diag_calib_param *pParam_Caplet,
		const diag_calib_lm_params *pCalibParam,
		/// output
		_Calib_LamTS_DiagEuro_Comm *pCalib_LamTS_DiagEuro
		)
{
	const int nSize_Ex = plEx_End-plEx_Begin;
	
	pCalib_LamTS_DiagEuro->m_pszSwapTenor = calloc(nSize_Ex,sizeof(char *));
	pCalib_LamTS_DiagEuro->m_pdSwapStrike = calloc(nSize_Ex,sizeof(double));
	pCalib_LamTS_DiagEuro->m_pszCapTenor=calloc(nSize_Ex,sizeof(char *));
	pCalib_LamTS_DiagEuro->m_pdCapStrike=calloc(nSize_Ex,sizeof(double));
	
	//pCalib_LamTS_DiagEuro->m_pParam_Swaption=calloc(1,sizeof(cpd_diag_calib_param));
	//pCalib_LamTS_DiagEuro->m_pParam_Caplet=calloc(1,sizeof(cpd_diag_calib_param));
	//pCalib_LamTS_DiagEuro->m_pCalibParam=calloc(1,sizeof(diag_calib_lm_params));
	
	pCalib_LamTS_DiagEuro->m_szYC = szYC ;
	pCalib_LamTS_DiagEuro->m_szVC= szVC;
	pCalib_LamTS_DiagEuro->m_get_cash_vol=get_cash_vol;
	*(long*)(&pCalib_LamTS_DiagEuro->m_lTheoEnd)=lTheoEnd;
	pCalib_LamTS_DiagEuro->m_szRefRate=szRefRate;
	pCalib_LamTS_DiagEuro->m_szFreq=szFreq;
	pCalib_LamTS_DiagEuro->m_szBasis=szBasis;
	*(const long**)(&pCalib_LamTS_DiagEuro->m_plEx_Begin)=plEx_Begin;
	*(const long**)(&pCalib_LamTS_DiagEuro->m_plEx_End)=plEx_End;
	*(int*)(&pCalib_LamTS_DiagEuro->m_nOnefequi)=nOnefequi;
	pCalib_LamTS_DiagEuro->m_pParam_Swaption=pParam_Swaption;
	pCalib_LamTS_DiagEuro->m_pParam_Caplet=pParam_Caplet;
	pCalib_LamTS_DiagEuro->m_pCalibParam=pCalibParam;
	
	////////
	Calibrate_MidAt_Preprocess(
	lToday,
	szYC,
	szRefRate, // Refrate
	szFreq, // frequency of the fixed leg -"M","Q","S","A" 
	pdCoupon_Begin,// coupon
	plEx_Begin,
	plEx_End,// exercise dates
	plFixStart_Begin,
	plFixStart_End,
	plFixEnd_Begin,
	plFixPay_Begin,
	plFltStart_Begin,
	plFltStart_End,
	plFltEnd_Begin,
	plFltPay_Begin,
	pdMargin_Begin,
	pdSpread_Begin,
	(char**)pCalib_LamTS_DiagEuro->m_pszSwapTenor,
	(double*)pCalib_LamTS_DiagEuro->m_pdSwapStrike,
	(char**)pCalib_LamTS_DiagEuro->m_pszCapTenor,
	(double*)pCalib_LamTS_DiagEuro->m_pdCapStrike
	);
}


void _free_Calib_LamTS_DiagEuro_Comm(
		_Calib_LamTS_DiagEuro_Comm *pCalib_LamTS_DiagEuro
		)
{
	const int nSize_Ex = pCalib_LamTS_DiagEuro->m_plEx_End-pCalib_LamTS_DiagEuro->m_plEx_Begin;
	const char **pszSwapTenor=pCalib_LamTS_DiagEuro->m_pszSwapTenor;
	const char **pszCapTenor=pCalib_LamTS_DiagEuro->m_pszCapTenor;
	int nI;

	free((double*)pCalib_LamTS_DiagEuro->m_pdSwapStrike);// = calloc(nSize_Ex,sizeof(double));
	free((double*)pCalib_LamTS_DiagEuro->m_pdCapStrike);//=calloc(nSize_Ex,sizeof(double));
	//free((cpd_diag_calib_param*)pCalib_LamTS_DiagEuro->m_pParam_Swaption);//=calloc(1,sizeof(cpd_diag_calib_param));
	//free((cpd_diag_calib_param*)pCalib_LamTS_DiagEuro->m_pParam_Caplet);//=calloc(1,sizeof(cpd_diag_calib_param));
	//free((diag_calib_lm_params*)pCalib_LamTS_DiagEuro->m_pCalibParam);//=calloc(1,sizeof(diag_calib_lm_params));

	for(nI=0;nI<nSize_Ex;++nI)
	{
		free((char*)pszSwapTenor[nI]);
		free((char*)pszCapTenor[nI]);
	}
	
	free((char**)pszSwapTenor);// = calloc(nSize_Ex,sizeof(char *));
	free((char**)pszCapTenor);//=calloc(nSize_Ex,sizeof(char *));
			
}



static int _initialize_ptr_TS_(
							 const double *pdIn_Begin, 
							 const double *pdIn_End,
							 double dKey
							 )
{
	const double *pd = find_min_greater_than(&dKey,pdIn_Begin,pdIn_End-pdIn_Begin,sizeof(dKey),dless);
	_ASSERTE(pd!=pdIn_End);// ASSERTE found!
	return (pd-pdIn_End);
}


// Map term structure in each _GenEuropean onto
// positions in global term structure
static void _initialize_ptr_TS(
								const double *pdZetaT_Begin,
								const double *pdZetaT_End,
								const double *pdGT_Begin,
								const double *pdGT_End,							   
								const double *pdTS, // 0 Based
								int nZeta_Begin,// beginning position of Zeta no pdTS 
								int nG_Begin,// beginning position of G no pdTS
								// results
								_GenEuropean* pInst
							   )
{
	const double *pdZeta = pdTS+nZeta_Begin;
	const double *pdG = pdTS+nG_Begin;
	const double *pdKey = &pInst->m_dEx;
	const double *pdG_CashFlow = pInst->m_pdG_CashFlow_Begin;

	// zeta at exercise
	*(const double**)(&pInst->m_pdZeta_Ex) = pdZeta+_initialize_ptr_TS_(pdZetaT_Begin,pdZetaT_End,*pdKey);
	// G at exercise
	*(const double**)(&pInst->m_pdG_Ex) = pdG + _initialize_ptr_TS_(pdGT_Begin,pdGT_End,*pdKey);

	// G at cash flow dates
	for(
		pdKey=pInst->m_pdT_CashFlow_Begin;
		pdKey<pInst->m_pdT_CashFlow_End;
		++pdKey,++pdG_CashFlow
		)
	{
		pdG_CashFlow =  pdG + _initialize_ptr_TS_(pdGT_Begin,pdGT_End,*pdKey);
	}

}

static const char* _Calibrate_(
				// calibration instruments
				const _GenEuropean** ppGenEuroBlock_Begin, 
				const _GenEuropean** ppGenEuroBlock_End,
				// Zeta T and GT
				const double *pdZetaT_Begin,
				const double *pdZetaT_End,
				const double *pdGT_Begin,
				const double *pdGT_End,
				// Results
				double *pdZeta_Begin,
				double *pdG_Begin
				)
{
	// local functor declaration
	void _SetNagFail_(NagError *pFail);
	void NAG_CALL _lsqfunc_(long,long,double*,double*,Nag_Comm*);

	// alias
	const _GenEuropean** ppInst = ppGenEuroBlock_Begin;
	const _GenEuropean** ppInst_end = ppGenEuroBlock_End;
	const _GenEuropean *pInst = 0;

	Nag_Comm comm;
	NagError fail;
	const int nNumZeta = pdZetaT_End-pdZetaT_Begin;
	const int nNumG = pdGT_End-pdGT_Begin;
	const int n= nNumZeta+nNumG; // number of unknowns
	const int m = ppInst_end-ppInst;// number of equations

	// global term structure
	double *pdTS = memset(_alloca(n*sizeof(double)),0,n*sizeof(double)), dfsumsq=0.;
	// begin positions of Zeta and G in the global term structrure
	const int nZeta_Begin = 0;
	const int nG_Begin = nNumZeta;
	
	// final residule
	double *pdf = _alloca(m*sizeof(double)); 
	// jacobian
	double *pdjac = _alloca(m*sizeof(double)); 
	
	// check that the problem is not over-determininistic
	_ASSERTE(m>=n);
	
	/// initialize term struture pointers in the calibration block
	for(;ppInst<ppInst_end;++ppInst)
	{
		for(pInst=*ppInst;pInst<*(ppInst+1);++pInst)
		{
			_initialize_ptr_TS(
				pdZetaT_Begin,
				pdZetaT_End,
				pdGT_Begin,
				pdGT_End,
				pdTS,
				nZeta_Begin,
				nG_Begin,
				(_GenEuropean *)pInst
				);
		}
	}

	// attach data to common Param_Caplet
	comm.p = (void*)ppGenEuroBlock_Begin;		

	// failure handling
	_SetNagFail_(&fail);

	// delegate to nag ... 
	nag_opt_lsq_no_deriv(
		m,
		n,
		_lsqfunc_,
		pdTS,
		&dfsumsq,
		pdf,
		pdjac,
		n,
		E04_DEFAULT,
		&comm,
		&fail);

	// error handling
	if(fail.code == NE_NOERROR)
		return fail.message;

	// assign to results
	memcpy(pdZeta_Begin,pdTS+nZeta_Begin,(nNumZeta-nZeta_Begin)*sizeof(double));

	memcpy(pdG_Begin,pdTS+nG_Begin,(nNumG-nG_Begin)*sizeof(double));
	
	return 0;	
}

static double _Price_GenEuroBlock(
								  const _GenEuropean* pGenEuroBlock_Begin,
								  const _GenEuropean* pGenEuroBlock_End
								  )
{
	const _GenEuropean *pInst = pGenEuroBlock_Begin;

	double dPrice = 0.;
	for(;pInst<pGenEuroBlock_End;++pInst)
	{
		dPrice += lgmopval1F(
			pInst->m_pdT_CashFlow_End-pInst->m_pdT_CashFlow_Begin,
			(double*)pInst->m_pdCashFlow_Begin,
			(double*)pInst->m_pdG_CashFlow_Begin,
			*pInst->m_pdZeta_Ex,
			*pInst->m_pdG_Ex
			);
	}

	return dPrice;
}


static void NAG_CALL _lsqfunc_(
									long m, // number of equations
									long n, // number of unknowns
									double xVec[], // x 
									double fVec[], 
									Nag_Comm* pComm
									)
{
	// declare local functor
	double _Price_GenEuroBlock(const _GenEuropean*,const _GenEuropean*);

	// alias
	const _GenEuropeanArg* pArg = (const _GenEuropeanArg*)(pComm->p);
	const _GenEuropean** ppInst = pArg->m_ppInst_Begin;
	const _GenEuropean** ppInst_end = ppInst+m;

	// overwrite the term structure
	_ASSERTE((pArg->m_pdTS_End-pArg->m_pdTS_Begin)==n);
	memcpy(pArg->m_pdTS_Begin,xVec,n*sizeof(*xVec));
	
	// price
	for(;ppInst<ppInst_end;++ppInst,++fVec)
		*fVec=_Price_GenEuroBlock(*ppInst,*(ppInst+1));
}

static void _SetNagFail_(NagError *pFail)
{
	SET_FAIL(*pFail);
	pFail->print = FALSE; // dont print!
}


static void _Equi_EuroSwapStrike(
		long lToday,
		const char *szYC,
		// exercise
		const long *plEx_Begin,
		const long *plEx_End,
		// fix
		const long *plFixStart_Begin,
		const long *plFixStart_End,
		const long *plFixEnd_Begin,
		const long *plFixPay_Begin,
		const double *pdCoupon_Begin,
		SrtBasisCode eFixBasis,
		// floating
		const long *plFltStart_Begin,
		const long *plFltStart_End,
		const long *plFltEnd_Begin,
		const long *plFltPay_Begin,
		const double *pdMargin_Begin,
		const double *pdSpread_Begin,
		SrtBasisCode eFltBasis,
		/// results
		double *pdStrike_Begin
		)
{
	const int nSize_Fix = plFixStart_End-plFixStart_Begin;
	const int nSize_Flt = plFltStart_End-plFltStart_Begin;
	double *pdFixNotional_Begin = _alloca(nSize_Fix*sizeof(double));
	//double *pdCoupon_Begin = _alloca(nSize_Fix*sizeof(double));
	double *pdFltNotional_Begin=_alloca(nSize_Flt*sizeof(double));
	const double dEqui_Notional = 1.;
		
	_memset(&dEqui_Notional,pdFixNotional_Begin,nSize_Fix,sizeof(dEqui_Notional));
	_memset(&dEqui_Notional,pdFltNotional_Begin,nSize_Flt,sizeof(dEqui_Notional));
	//_memset(&dCoupon,pdCoupon_Begin,nSize_Fix,sizeof(*pdCoupon_Begin));
	
	for(;plEx_Begin<plEx_End;++plEx_Begin,++pdStrike_Begin)
	{
		/// locate indeces 
		const long *plFixStart_Begin_ = find_min_greater_than(plEx_Begin,plFixStart_Begin,nSize_Fix,sizeof(*plEx_Begin),lless);
		int nFixOffSet = plFixStart_Begin_ - plFixStart_Begin;
		
		const long *plFltStart_Begin_ = find_min_greater_than(plEx_Begin,plFltStart_Begin,nSize_Flt,sizeof(*plEx_Begin),lless);
		int nFltOffSet = plFltStart_Begin_ - plFltStart_Begin;
		
		//_ASSERTE(nFixOffSet!=nSize_Fix);
		//_ASSERTE(nFltOffSet!=nSize_Flt);
		
		// NB: in case indeces are not located, e.g. the last call in call freq = q and fix freq = s
		*pdStrike_Begin = Equi_Strike(
			lToday,
			szYC,
			nFixOffSet!=nSize_Fix?plFixStart_Begin_:plEx_Begin,
			nFixOffSet!=nSize_Fix?plFixStart_End:plEx_Begin+1,
			nFixOffSet!=nSize_Fix?plFixEnd_Begin+nFixOffSet:plFixEnd_Begin+nSize_Fix-1,
			nFixOffSet!=nSize_Fix?plFixPay_Begin+nFixOffSet:plFixPay_Begin+nSize_Fix-1,
			nFixOffSet!=nSize_Fix?pdCoupon_Begin+nFixOffSet:pdCoupon_Begin+nSize_Fix-1,
			nFixOffSet!=nSize_Fix?pdFixNotional_Begin+nFixOffSet:pdFixNotional_Begin+nSize_Fix-1,
			eFixBasis,
			nFltOffSet!=nSize_Flt?plFltStart_Begin_:plEx_Begin,
			nFltOffSet!=nSize_Flt?plFltStart_End:plEx_Begin+1,
			nFltOffSet!=nSize_Flt?plFltEnd_Begin+nFltOffSet:plFltEnd_Begin+nSize_Flt-1 ,
			nFltOffSet!=nSize_Flt?plFltPay_Begin+nFltOffSet:plFltPay_Begin+nSize_Flt-1 ,
			nFltOffSet!=nSize_Flt?pdMargin_Begin+nFltOffSet:pdMargin_Begin+nSize_Flt-1 ,
			nFltOffSet!=nSize_Flt?pdSpread_Begin+nFltOffSet:pdSpread_Begin+nSize_Flt-1,
			nFltOffSet!=nSize_Flt?pdFltNotional_Begin+nFltOffSet:pdFltNotional_Begin+nSize_Flt-1,
			eFltBasis,
			dEqui_Notional
			);
	}
}

static void _ATM_3MStrike(
		long lToday,
		const char *szYC,
		const long *plEx_Begin,
		const long *plEx_End,
		const long *plFltStart_Begin,
		const long *plFltStart_End,
		const double *pdMargin_Begin,
		const double *pdSpread_Begin,
		SrtBasisCode eFltBasis,
		/// results
		double *pdStrike_Begin
		)
{
	const int nSize_Flt = plFltStart_End-plFltStart_Begin;

	for(;plEx_Begin<plEx_End;++plEx_Begin,++pdStrike_Begin)
	{
		const long lMat_3M = add_unit(*plEx_Begin,3,SRT_MONTH,NO_BUSDAY_CONVENTION);
		const long *plFltStart_Begin_ = find_min_greater_than(plEx_Begin,plFltStart_Begin,nSize_Flt,sizeof(*plEx_Begin),lless);
		const int nFltOffSet = plFltStart_Begin_ - plFltStart_Begin;

		*pdStrike_Begin = Swap_Rate(
			lToday,
			szYC,
			plEx_Begin,
			1+plEx_Begin,// size 1
			&lMat_3M,
			&lMat_3M,
			eFltBasis,
			plEx_Begin,
			1+plEx_Begin,
			&lMat_3M,
			&lMat_3M,
			pdMargin_Begin+nFltOffSet,
			pdSpread_Begin+nFltOffSet,
			eFltBasis
			);
	}
}

void Calibrate_MidAt_Preprocess(
	///////Market	
	long lToday,
	const char *szYC,
	/// MidAt General Specs
	const char *szRefRate, // Refrate
	const char *szFreq, // frequency of the fixed leg -"M","Q","S","A" 
	const double *pdCoupon_Begin,// coupon
	// MidAt Exercise
	const long *plEx_Begin,
	const long *plEx_End,// exercise dates
	// MidAt Fixed Leg 
	const long *plFixStart_Begin,
	const long *plFixStart_End,
	const long *plFixEnd_Begin,
	const long* plFixPay_Begin,
	// MidAt Funding/Floating leg 
	const long *plFltStart_Begin,
	const long *plFltStart_End,
	const long *plFltEnd_Begin,
	const long* plFltPay_Begin,
	const double *pdMargin_Begin,
	const double *pdSpread_Begin,
	// results
	char **pszSwapTenor,
	double *pdSwapStrike,
	char **pszCapTenor,
	double *pdCapStrike
	)
{
	const int nSize_Ex = plEx_End-plEx_Begin;

	// Swaption tenor
	const char *szSwapTenor = "DIAG";
	// cap tenor
	char *szCapTenor = "3M";
	
	int nI=0;

	SrtBasisCode eFixBasis,eFltBasis;
	SrtCompounding eFltFreq;
	interp_basis(szFreq,&eFixBasis);
	swp_f_get_ref_rate_details((char*)szRefRate,&eFltBasis,&eFltFreq);	
	
	// fill in swap and cap tenors
	for(nI=0;nI<nSize_Ex;++nI)
	{
		pszSwapTenor[nI] = calloc(1,sizeof("DIAG"));
		pszCapTenor[nI] = calloc(1,sizeof("3M"));
		_memset("DIAG",pszSwapTenor[nI],1,sizeof("DIAG"));
		_memset("3M",pszCapTenor[nI],1,sizeof("3M"));
	}
	
	// fill in equivalent swap strike
	_Equi_EuroSwapStrike(
		lToday,
		szYC,
		plEx_Begin,
		plEx_End,
		plFixStart_Begin,
		plFixStart_End,
		plFixEnd_Begin,
		plFixPay_Begin,
		pdCoupon_Begin,
		eFixBasis,
		plFltStart_Begin,
		plFltStart_End,
		plFltEnd_Begin,
		plFltPay_Begin,
		pdMargin_Begin,
		pdSpread_Begin,
		eFltBasis,
		pdSwapStrike 
		);

#ifdef _DEBUG
	//_print_vector(pdSwapStrike,plEx_End-plEx_Begin);
#endif 
	
	// fill in ATM Cap strike
	_ATM_3MStrike(
		lToday,
		szYC,
		plEx_Begin,
		plEx_End,
		plFltStart_Begin,
		plFltStart_End,
		pdMargin_Begin,
		pdSpread_Begin,
		eFltBasis,
		pdCapStrike
		);

#ifdef _DEBUG
	//_print_vector(pdCapStrike,plEx_End-plEx_Begin);
#endif 

}

const char *_Calibrate_MidAt(
	///////Market	
	const char *szYC,
	const char *szVC,
	Err (*get_cash_vol)(char*,double,double,double,int,char	*,double*,double*),	//	functor to get cash vol from the market
	/// MidAt General Specs
	long lTheoEnd,
	const char *szRefRate, // Refrate
	const char *szFreq, // frequency of the fixed leg -"M","Q","S","A" 
	const char *szBasis, // basis of the fixed leg
	// MidAt Exercise
	const long *plEx_Begin,
	const long *plEx_End,// exercise dates
	// calib params	
	const char **pszSwapTenor_Begin,
	const double *pdSwapStrike_Begin,
	const char **pszCapTenor_Begin, 
	const double *pdCapStrike_Begin,
	int nCalibTauToCap, // 1 - calibrate Tau to cap, 0 otherwise
	const double *pdAlpha,
	const double *pdGamma,
	const double *pdRho,
	const cpd_diag_calib_param *pParam_Swaption,
	const cpd_diag_calib_param *pParam_Caplet,
	const diag_calib_lm_params *pCalibParam,
	int nOnefequi,
	const double *pdLamT_Begin, // NB: might be changed!
	const double *pdLamT_End,// NB: might be changed!
	double *pdLam_Begin,// NB: might be changed!
	// results
	double **ppdSigT_Begin,
	double **ppdSigT_End,
	double **ppdSig
	)
{
	const int nSize_Ex = plEx_End-plEx_Begin;
	const int nSize_Lam = pdLamT_End-pdLamT_Begin;
	int nSize_Sig=0;
	const int nOne = 1;
	const double dZero = 0.;
	const int nNumFactor = pdAlpha&&pdGamma&&pdRho? 2:1;
	const char *szErr =0;

	szErr =  cpd_calib_diagonal_dlm(
		 (char*)szYC,
		 (char*)szVC,
		 get_cash_vol,
		 (char*)szRefRate,	
		 ///////////////////////
		 // Diag Swaptions
		 (char*)szFreq,
		 (char*)szBasis,
		 (char*)szRefRate,
		 nSize_Ex,
		 (long*)plEx_Begin,
		 // all 1's to use exercise dates as calibration dates
		 (int*)_memset(&nOne,_alloca(nSize_Ex*sizeof(int)), nSize_Ex,sizeof(int)),
		 (char**)pszSwapTenor_Begin,// "DIAG"
		lTheoEnd,
		(double*)pdSwapStrike_Begin, ///  	
		(cpd_diag_calib_param*)pParam_Swaption,
		///////////////////////
		// Caplets
		"Q", // Always 3M caplet ???
		"MM", // Always "MM"?
		(char*)szRefRate,	
		nSize_Ex,	
		(long*)plEx_Begin,
		// all 1's to use exercise dates as calibration dates
		(int*)_memset(&nOne,_alloca(nSize_Ex*sizeof(int)), nSize_Ex,sizeof(int)),
		(char**)pszCapTenor_Begin,
		lTheoEnd,
		(double*)pdCapStrike_Begin,
		0,//pdOne,//weights
		(cpd_diag_calib_param*)pParam_Caplet,
		nCalibTauToCap?0:1,
		nOnefequi,
		nSize_Lam,
		(double*)pdLamT_Begin,
		//nSize_Ex==1?(double*)_memset(&dLam_Euro,_alloca(nSize_Lam*sizeof(double)),nSize_Lam,sizeof(double)):pdLam_Begin,
		pdLam_Begin,
		(double*)_memset(&dZero,_alloca(nSize_Lam*sizeof(double)),nSize_Lam,sizeof(double)), // lambda shift
		nNumFactor,
		nNumFactor==1?1.e-5:*pdAlpha,
		nNumFactor==1?1.e-5:*pdGamma,
		nNumFactor==1?1.e-5:*pdRho,
		// Shift parameters
		0,0,0,0,0,0,
		//	Output 
		&nSize_Sig,
		ppdSigT_Begin,
		ppdSig,
		(diag_calib_lm_params*)pCalibParam,
		0);

	*ppdSigT_End = *ppdSigT_Begin+nSize_Sig;
		
	return szErr;
}

/// calibrate the diagnoal europeans and possibly caplets 
/// underlying a single regular midat
const char * Calibrate_MidAt(
	///////Market	
	long lToday,
	const char *szYC,
	const char *szVC,
	Err (*get_cash_vol)(char*,double,double,double,int,char	*,double*,double*),	//	functor to get cash vol from the market
	/// MidAt General Specs
	long lTheoEnd,
	const char *szRefRate, // Refrate
	const char *szFreq, // frequency of the fixed leg -"M","Q","S","A" 
	const char *szBasis, // basis of the fixed leg
	const double *pdCoupon_Begin,
	// MidAt Exercise
	const long *plEx_Begin,
	const long *plEx_End,// exercise dates
	// MidAt Fixed Leg 
	const long *plFixStart_Begin,
	const long *plFixStart_End,
	const long *plFixEnd_Begin,
	const long* plFixPay_Begin,
	// MidAt Funding/Floating leg 
	const long *plFltStart_Begin,
	const long *plFltStart_End,
	const long *plFltEnd_Begin,
	const long* plFltPay_Begin,
	const double *pdMargin_Begin,
	const double *pdSpread_Begin,
	// term structure
	int nCalibTauToCap, // 2 - calibrate lambda ATM, then prim ATS:  1 - calibrate Tau to cap, 0 otherwise
	const double *pdAlpha,
	const double *pdGamma,
	const double *pdRho,
	// calibration Param_Caplet
	const cpd_diag_calib_param *pParam_Swaption,
	const cpd_diag_calib_param *pParam_Caplet,
	const diag_calib_lm_params *pCalibParam,
	int nOnefequi,
	double *pdLamT_Begin, // NB: might be changed!
	double *pdLamT_End,// NB: might be changed!
	double *pdLam_Begin,// NB: might be changed!
	// results
	double **ppdSigT_Begin,
	double **ppdSigT_End,
	double **ppdSig_Begin
	)
{
// Variable declaration
	cpd_diag_calib_param copy_SwapParam, copy_CapParam;
	const char* err = 0;
	const int nSize_Ex = plEx_End-plEx_Begin;

	// Swaption tenor
	char **pszSwapTenor = _alloca(nSize_Ex*sizeof(char *));
	double *pdSwapStrike = _alloca(nSize_Ex*sizeof(double));
	
	// cap tenor
	char **pszCapTenor = _alloca(nSize_Ex*sizeof(char *));
	double *pdCapStrike = _alloca(nSize_Ex*sizeof(double));
	
// copy the parameters
	copy_SwapParam = *pParam_Swaption;
	copy_CapParam = *pParam_Caplet;

// check to see if we do a two-step calibration
	if ( nCalibTauToCap == 2 )
	{
// set the values to ATM
		copy_CapParam.strike_type = 5;
	}

	// delegate to prepare for the call to diagcalibdlm
	Calibrate_MidAt_Preprocess(
		lToday,
		szYC,
		szRefRate, // Refrate
		szFreq, // frequency of the fixed leg -"M","Q","S","A" 
		pdCoupon_Begin,// coupon
		plEx_Begin,
		plEx_End,
		plFixStart_Begin,
		plFixStart_End,
		plFixEnd_Begin,
		plFixPay_Begin,
		plFltStart_Begin,
		plFltStart_End,
		plFltEnd_Begin,
		plFltPay_Begin,
		pdMargin_Begin,
		pdSpread_Begin,
		pszSwapTenor,
		pdSwapStrike,
		pszCapTenor,
		pdCapStrike
		);
	
	// delegate
	return _Calibrate_MidAt(
		szYC,
		szVC,
		get_cash_vol,
		lTheoEnd,
		szRefRate, 
		szFreq, 
		szBasis,
		plEx_Begin,
		plEx_End,
		pszSwapTenor,
		pdSwapStrike,
		pszCapTenor, 
		pdCapStrike,
		nCalibTauToCap, 
		pdAlpha,
		pdGamma,
		pdRho,
		&copy_SwapParam,
		&copy_CapParam,
		pCalibParam,
		nOnefequi,
		pdLamT_Begin,
		pdLamT_End,
		pdLam_Begin,
		ppdSigT_Begin,
		ppdSigT_End,
		ppdSig_Begin
		);
}


/// calibrate the diagnoal europeans and possibly caplets 
/// underlying a single generic midat
const char * Calibrate_GenMidAt(
		/// Market
		long lToday,
		const char	*szYC,
		const char	*szVC,
		Err	(*get_cash_vol)(char*,double,double,double,int,char*,double*,double*),	
		const char	*szCorrel,
		Err	(*get_correl)(char*, double,double,double,double*),
		// GenMidAt general info
		int	eod_flag,//	EOD Flag 0: I, 1: E 
		long lTheoEnd,			
		const char	*szRefRate,	//ref rate
		const char *szFreq,//swap_freq, //swap freq 
		const char *szBasis, //	swap basis 
		// Exercise info		
		const long *plEx_Begin,
		const long *plEx_End,
		// Fixed leg
		const long *plFixStart_Begin,		
		const long *plFixStart_End,
		const long *plFixEnd_Begin,
		const long *plFixPay_Begin,
		const double *pdFixRate_Begin,
		const double *pdFixNotional_Begin,
		const double *pdFixFee_Begin,
		// Floating leg
		const long *plFltStart_Begin,		
		const long *plFltStart_End,
		const long *plFltEnd_Begin,
		const long *plFltPay_Begin,
		const double	*pdMargin_Begin,
		const double	*pdSpread_Begin,
		const double *pdFltNotional_Begin,
		// Given term structure 
		const double *pdLamT_Begin,
		const double *pdLamT_End,
		const double *pdLam_Begin,				
		const double *pdAlpha,		
		const double *pdGamma,	
		const double *pdRho,
		//	Calib params 
		const char	*default_ref,	//	ref rate
		const char	*default_swap_freq,//swap freq 
		const char	*default_swap_basis,//swap basis
		double		dMinTime,
		double		dMinInterval,
		long lNotice,
		double dMaxStdShort,
		int			nFixLam,//0: calib lambda to cap, 1: fix lambda calib	to diagonal
		int			n1FEqui,//1: 2F lambda will calibrate to the cap priced within calibrated 1F	with the given lambda 
		int			nSkipLast,//If 1, the last option is disregardedand the forward volatility is flat from option	n-1 
		int			nUseJump,
		double		dMaxVarJump,
		int			nStrikeType,
		int			nEuroModel,
		// Calibration results
		double **ppdSigT_Begin,
		double **ppdSigT_End,
		double **ppdSig_Begin,
		// market prices
		int *pnExBool,
		double *pdDiagPrice
		)
{
	const int nSize_Ex = plEx_End-plEx_Begin;
	const int nSize_Fix = plFixStart_End-plFixStart_Begin;
	const int nSize_Flt = plFltStart_End-plFltStart_Begin;
	const int	nUseVol = nStrikeType == 5?1:0;

	int nSize_Sig;
	const char *szFreq_Flt=_alloca(256*sizeof(char));
	const char *szBasis_Flt=_alloca(256*sizeof(char));

	const char *szErr=0;
	
	if(nSize_Ex <= 0 || (nSize_Ex == 1 && plEx_Begin[0] <= lToday + eod_flag)) 
		return 0;

	if(nEuroModel!=1)
		return "Only SMM model currently available for pricing amortizing/accreting european swaptions!";

	/// Delegate to compute diag generic european option prices ...
	szErr = Price_DiagGenEuro(
			lToday,
			szYC,
			szVC,
			get_correl,
			szCorrel,
			szRefRate,
			szFreq,
			szBasis,
			lTheoEnd,			
			*plEx_Begin,
			plFixStart_Begin,
			plFixStart_End,
			plFixEnd_Begin,
			pdFixNotional_Begin,
			pdFixRate_Begin,
			pdFixFee_Begin,
			plFltStart_Begin,
			plFltStart_End,
			plFltEnd_Begin,
			pdFltNotional_Begin,
			pdMargin_Begin,
			pdSpread_Begin,
			dMinTime, 
			dMinInterval, 
			nUseVol,
			pnExBool, 
			pdDiagPrice
			);

	if(szErr) return szErr;


	szErr =Decode_RefRate(szRefRate,(char**)(&szFreq_Flt),(char**)(&szBasis_Flt));
	if(szErr) return szErr;

	_ASSERTE((pdLamT_End- pdLamT_Begin) >=1);
	
	if( (pdLamT_End- pdLamT_Begin) ==1)
	{
		szErr = amortMidat_cpd_calib_diagonal_new(
					lNotice,
					(char*)szYC,
					(char*)szVC,
					(char*)default_ref,
					(char*)default_swap_basis,
					(char*)default_swap_freq,
					get_cash_vol,
					0,
					1, /// shift type
					pnExBool,
					pdDiagPrice,
					memcpy(_alloca(sizeof(double)*nSize_Fix),pdFixFee_Begin,sizeof(double)*nSize_Fix),
					(char*)szFreq,
					(char*)szBasis,
					nSize_Fix,
					memcpy(_alloca(sizeof(long)*nSize_Fix),plFixStart_Begin,sizeof(long)*nSize_Fix),
					memcpy(_alloca(sizeof(long)*nSize_Fix),plFixEnd_Begin,sizeof(long)*nSize_Fix),
					memcpy(_alloca(sizeof(long)*nSize_Fix),plFixPay_Begin,sizeof(long)*nSize_Fix),//(long*)plFixPay_Begin,
					memcpy(_alloca(sizeof(double)*nSize_Fix),pdFixRate_Begin,sizeof(double)*nSize_Fix),//(double*)pdFixRate_Begin,
					memcpy(_alloca(sizeof(double)*nSize_Fix),pdFixNotional_Begin,sizeof(double)*nSize_Fix),//(double*)pdFixNotional_Begin,
					(char*)szFreq_Flt,
					(char*)szBasis_Flt,
					nSize_Flt,
					memcpy(_alloca(sizeof(long)*nSize_Flt),plFltStart_Begin,sizeof(long)*nSize_Flt),
					memcpy(_alloca(sizeof(long)*nSize_Flt),plFltEnd_Begin,sizeof(long)*nSize_Flt),
					memcpy(_alloca(sizeof(long)*nSize_Flt),plFltPay_Begin,sizeof(long)*nSize_Flt),//(long*)plFltPay_Begin,
					memcpy(_alloca(sizeof(double)*nSize_Flt),pdMargin_Begin,sizeof(double)*nSize_Flt),//(double*)pdMargin_Begin,
					memcpy(_alloca(sizeof(double)*nSize_Flt),pdSpread_Begin,sizeof(double)*nSize_Flt),//(double*)pdSpread_Begin,
					memcpy(_alloca(sizeof(double)*nSize_Flt),pdFltNotional_Begin,sizeof(double)*nSize_Flt),//(double*)pdFltNotional_Begin,
					memcpy(_alloca(sizeof(double)*nSize_Fix),pdFixRate_Begin,sizeof(double)*nSize_Fix),//(double*)pdFixRate_Begin,
					nStrikeType,
					dMaxStdShort,
					nFixLam,
					n1FEqui,
					nSkipLast,
					nUseJump,
					dMaxVarJump,
					memcpy(_alloca((pdLamT_End-pdLamT_Begin)*sizeof(double)),pdLam_Begin,(pdLamT_End-pdLamT_Begin)*sizeof(double)),
					pdAlpha&&pdGamma&&pdRho?2:1,
					pdAlpha?*pdAlpha:1.e-5,
					pdGamma?*pdGamma:1.e-5,
					pdRho?*pdRho:1.e-5,
					&nSize_Sig,
					ppdSigT_Begin,
					ppdSig_Begin
					);

		*ppdSigT_End=*ppdSigT_Begin+nSize_Sig;
	}
	else
	{
		szErr = amortMidat_cpd_calib_diagonal_new_ts(
					lNotice,
					(char*)szYC,
					(char*)szVC,
					(char*)default_ref,
					(char*)default_swap_basis,
					(char*)default_swap_freq,
					get_cash_vol,
					0,
					1, /// shift type
					pnExBool,
					pdDiagPrice,
					(double*)pdFixFee_Begin,
					(char*)szFreq,
					(char*)szBasis,
					nSize_Fix,
					(long*)plFixStart_Begin,
					(long*)plFixEnd_Begin,
					(long*)plFixPay_Begin,
					(double*)pdFixRate_Begin,
					(double*)pdFixNotional_Begin,

					(char*)szFreq_Flt,
					(char*)szBasis_Flt,
					nSize_Flt,
					(long*)plFltStart_Begin,
					(long*)plFltEnd_Begin,
					(long*)plFltPay_Begin,
					(double*)pdMargin_Begin,
					(double*)pdSpread_Begin,
					(double*)pdFltNotional_Begin,

					(double*)pdFixRate_Begin,
					
					nStrikeType,
					
					dMaxStdShort,

					nFixLam,
					n1FEqui,
					
					nSkipLast,

					nUseJump,
					dMaxVarJump,
					

					pdLamT_End-pdLamT_Begin,
					(double*)pdLam_Begin,
					(double*)pdLamT_Begin,

					pdAlpha&&pdGamma&&pdRho?2:1,
					pdAlpha?*pdAlpha:1.e-5,
					pdGamma?*pdGamma:1.e-5,
					pdRho?*pdRho:1.e-5,
					
					&nSize_Sig,
					ppdSigT_Begin,
					ppdSig_Begin
					);
	}

	*ppdSigT_End = nSize_Sig+*ppdSigT_Begin;

	return szErr;
}

static const char* _nag_return(
	   const char *szErr,
	   Nag_Comm *pNagComm
	   )
{
	_ASSERTE(szErr);
	// set flag to some negative number to cause nag optimzer to terminate
	pNagComm->flag=-1;
	return szErr;
}

static void _stdcall _OBJFUN_LAMTS_CONVERTSIG(const double *pdx,int nEle,double *pdsin)
{
	const double *pdxend = pdx+nEle;
	
	for(;pdx<pdxend;++pdx,++pdsin)
	{
		*pdsin = *pdx;//)*(*pdx);//sin(*pdx);
		//*pdsin =(*pdx)*(*pdx);//)*sin(*pdx);//*sin(*pdx);//*sin(*pdx*1.e-1);
		//_ASSERTE(*pdsin>0.);
	}
}

static void _stdcall _OBJFUN_LAMTS_CONVERTLAM (const double *pdx,int nEle,double *pdsin)
{
	const double *pdxend = pdx+nEle;
	for(;pdx<pdxend;++pdx,++pdsin)
	{
		*pdsin = *pdx;
		//*pdsin=*pdx;//sin(*pdx)*cos(*pdx);
	}
}

static void __stdcall _OBJFUN_LAMTS_CONVERT(
	const double *pdx,
	_Model *pModel
	)
{
	const int nSize_Sig = pModel->m_pdSigT_End-pModel->m_pdSigT_Begin;
	const int nSize_Lam = pModel->m_pdLamT_End-pModel->m_pdLamT_Begin;
	/// convert pdx to sigma
	_OBJFUN_LAMTS_CONVERTSIG(pdx,nSize_Sig,(double*)pModel->m_pdSig_Begin);
	/// convert pdx to lambda
	_OBJFUN_LAMTS_CONVERTLAM (pdx+nSize_Sig,nSize_Lam,(double*)pModel->m_pdLam_Begin);
}

static const char *_OBJFUN_LAMTS_CALIBDIAG(
	const _Calib_LamTS_DiagEuro_Comm* pCalibToEuro,
	// model
	const double *pdLamT_Begin,
	const double *pdLamT_End,
	const double *pdLam_Begin, // unknowns
	const double *pdGamma,
	const double *pdAlpha,
	const double *pdRho,
	// output
	double**ppdSigT_Begin,
	double**ppdSigT_End,
	double**ppdSig_Begin
    )
{

	// first calibrate to its europeans
	return _Calibrate_MidAt(
		pCalibToEuro->m_szYC,
		pCalibToEuro->m_szVC,
		pCalibToEuro->m_get_cash_vol, /// get vol functor
		pCalibToEuro->m_lTheoEnd,
		pCalibToEuro->m_szRefRate, 
		pCalibToEuro->m_szFreq, 
		pCalibToEuro->m_szBasis,
		pCalibToEuro->m_plEx_Begin,
		pCalibToEuro->m_plEx_End,
		pCalibToEuro->m_pszSwapTenor,
		pCalibToEuro->m_pdSwapStrike,
		pCalibToEuro->m_pszCapTenor, 
		pCalibToEuro->m_pdCapStrike,
		0,//nCalibTauToCap, 
		pdAlpha,
		pdGamma,
		pdRho,
		pCalibToEuro->m_pParam_Swaption,
		pCalibToEuro->m_pParam_Caplet,
		pCalibToEuro->m_pCalibParam,
		pCalibToEuro->m_nOnefequi,
		pdLamT_Begin,
		pdLamT_End,
		(double*)pdLam_Begin,
		ppdSigT_Begin,
		ppdSigT_End,
		ppdSig_Begin
		);
}

static const char *_OBJFUN_LAMTS_BOUNDS_AUTOCAL(
	const _Calib_LamTS_DiagEuro_Comm* pCalibToEuro,
	const _GenMidAt *pOption,
	const _PCQ_Seq *pNumeraire,
	// model
	const double *pdLamT_Begin,
	const double *pdLamT_End,
	const double *pdLam_Begin, // unknowns
	const double *pdGamma,
	const double *pdAlpha,
	const double *pdRho,
	/// Grid,
	int nNumPointT,
	int nNumPointX,
	// result
	double *pPrice
    )
{
	double *pdSigT_Begin=0,*pdSigT_End=0,*pdSig_Begin=0;
	
	// calibrate sigmas
	const char *szErr = _OBJFUN_LAMTS_CALIBDIAG(
		pCalibToEuro,
		pdLamT_Begin,
		pdLamT_End,
		pdLam_Begin,
		pdGamma,
		pdAlpha,
		pdRho,
		&pdSigT_Begin,
		&pdSigT_End,
		&pdSig_Begin
		);

	// create model that contains the sigma and lambda
	const _Model MyModel = {
		pdSigT_Begin,
		pdSigT_End,
		pdSig_Begin,
		pdLamT_Begin,
		pdLamT_End,
		pdLam_Begin,
		pdAlpha&&pdGamma&&pdRho?2:1,
		pdAlpha?*pdAlpha:0,
		pdGamma?*pdGamma:0,
		pdRho?*pdRho:0,
		MyModel.m_dAlpha*MyModel.m_dRho,//const double m_dAlphaRho;
		MyModel.m_dAlpha*sqrt(1-MyModel.m_dRho*MyModel.m_dRho),//const double m_dAlpha_Sqrt1mRhoSqrd;
		MyModel.m_dGamma*MyModel.m_dRho/sqrt(1.-MyModel.m_dRho*MyModel.m_dRho),//const double m_dGammaRho_Over_Sqrt1mRhoSqrd;
		0
	};

	if(szErr)
	{
		free(pdSig_Begin);
		free(pdSigT_Begin);
		return szErr;
	}

	// fill in models
	if(0){_init_Model(
		pdSigT_Begin,
		pdSigT_End,
		pdSig_Begin,
		pdLamT_Begin,
		pdLamT_End,
		pdLam_Begin,
		pdAlpha,
		pdGamma,
		pdRho,
		0//&MyModel
		);
	}
	
	/// price
	*pPrice = 0.;
	if(szErr = _GenMidAt_OptionPV_Core(nNumPointT,nNumPointX,pOption,&MyModel,pNumeraire,pPrice))
	{
		//free(pdSig_Begin);
		//free(pdSigT_Begin);
		return szErr;
	}
		
	//free(pdSig_Begin);
	//free(pdSigT_Begin);
	return 0;

}

static void NAG_CALL  _OBJFUN_LAMTS_BOUNDS(
	long n,
	double *pdx,
	double *pdf,
	double *pdg,
	Nag_Comm *pNagComm
	)
{
	// unpack common parameters
	_Calib_LamTS_Comm* pComm = (_Calib_LamTS_Comm*)(pNagComm->p);
	
	// alias
	const int nNumPointT=pComm->m_nNumPointT;
	const int nNumPointX=pComm->m_nNumPointX;
	const _GenMidAt *pOption_Begin=pComm->m_pOption_Begin;
	const _GenMidAt *pOption_End=pComm->m_pOption_End;
	const int m=pOption_End-pOption_Begin;
	const double *pOptionPrice_Begin=pComm->m_pOptionPrice_Begin;
	const _PCQ_Seq *pNumeraire = pComm->m_pNumeraire_Begin;
	//const char *szErr = pComm->m_szErr;

	// more alias
	const _Calib_LamTS_DiagEuro_Comm* pCalib_LamTS_DiagEuro_Begin=pComm->m_pCalib_LamTS_DiagEuro_Begin;
	const double *pdAlpha=pComm->m_pdAlpha;
	const double *pdGamma=pComm->m_pdGamma;
	const double *pdRho=pComm->m_pdRho;

	double *pdLamT_Begin=pComm->m_pdLamT_Begin;
	double *pdLamT_End=pdLamT_Begin+n;
	double *pdf_=pComm->m_pdf;
	_ASSERTE(pComm->m_pdLamT_End==pdLamT_End);

	*pdf=0.;
	
	for(;pOption_Begin<pOption_End;
	++pOption_Begin,
	++pOptionPrice_Begin,
	++pNumeraire,
	++pCalib_LamTS_DiagEuro_Begin,
	++pdf_
	)
	{
		const char *szErr = _OBJFUN_LAMTS_BOUNDS_AUTOCAL(
			pCalib_LamTS_DiagEuro_Begin,
			pOption_Begin,
			pNumeraire,
			pdLamT_Begin,
			pdLamT_End,
			pdx, // unknowns,i.e. lambda
			pdGamma,
			pdAlpha,
			pdRho,
			nNumPointT,
			nNumPointX,
			// result
			pdf_
			);
		
		_ASSERTE(*pdf_!=0.);
		
		if(szErr) 	
		{
			pComm->m_szErr=memset(calloc(1,strlen(szErr)),0,strlen(szErr));
			strcpy((char*)pComm->m_szErr,szErr);
			pNagComm->flag=-1;//NE_USER_STOP;			
		}

		*pdf += (*pdf_-*pOptionPrice_Begin)*(*pdf_-*pOptionPrice_Begin);
	}

	*pdf/=m;// in BP
	*pdf=sqrt(*pdf);

#ifdef _DEBUG
	_print_vector(pdx,m);
	_print_vector(pdf,1);
#endif 

	if(*pdf<1.e-4)
	{
		pNagComm->flag=-1;//NE_USER_STOP;			
		_ASSERTE(*pdf!=0.);
	}
	
	memcpy(pComm->m_pdx,pdx,sizeof(double)*m);
	
	// not referenced
	pdg;
}

static void _Set_NagFail(
	NagError *pfail
	)
{
	SET_FAIL(*pfail);
	pfail->print = FALSE; // dont print!
}

static void NAG_CALL my_print_fun(
	const Nag_Search_State *st, 
	Nag_Comm *comm
	)
{
	/// do nothing
	st,comm;
}


static void _Set_NagOption(
	Nag_E04_Opt *poptions
	)
{
	nag_opt_init(poptions);
	poptions->list = FALSE;
	poptions->print_level = Nag_NoPrint;
	poptions->print_fun = my_print_fun;
	poptions->optim_tol = 1.e-2;
}

const char *_Calibrate_TS_(
	int nNumPointT,
	int nNumPointX,
	const _GenMidAt *pOption_Begin,
	const _GenMidAt *pOption_End,
	const _PCQ_Seq *pNumeraire_Begin,
	const double *pOptionPrice_Begin,
	const _Calib_LamTS_DiagEuro_Comm* pCalib_LamTS_DiagEuro_Begin,
	// model
	const double *pdLamT_Begin,
	const double *pdLamT_End,
	const double *pdAlpha,
	const double *pdGamma,
	const double *pdRho,
	/// output
	double *pdLam_InitGuess_Begin // lambda ts
	)
{
	const char *szErr = 0;
	const int n = pdLamT_End-pdLamT_Begin;// # of unknowns
	const int m = pOption_End-pOption_Begin;// # of equations

	// result of least squares optimization
	double df= 0.; 
		
	// initialize common params
	const _Calib_LamTS_Comm Comm = {
	nNumPointT,
	nNumPointX,
	pOption_Begin,
	pOption_End,
	pOptionPrice_Begin,// market price
	pNumeraire_Begin,
	memset(_alloca(sizeof(double)*m),0,sizeof(double)*m),
	pCalib_LamTS_DiagEuro_Begin,
	pdAlpha,
	pdGamma,
	pdRho,
	(double*)pdLamT_Begin,
	(double*)pdLamT_End,
	0,//memset(_alloca(225),0,225)
	memset(_alloca(sizeof(double)*n),0,sizeof(double)*n)
	};

	if(m<n) return "_Calibrate_LamTS(...) : over-deterministic system!";

	/// delegate to nag solver .... 
	{
		int nI=0;

		/// specify boundaries
		const double dLam_LB = -0.3;
		const double dLam_UB = 0.3;
		double *pdbl = _alloca(n*sizeof(double));
		double *pdbu = _alloca(n*sizeof(double));
		// nag failure handling
		NagError fail;
		// nag options
		Nag_E04_Opt options;
		// nag common params
		Nag_Comm commParams;
		// attach data to nag common parameters
		commParams.p = (_Calib_LamTS_Comm*) (&Comm);		
		_memset(&dLam_LB,pdbl,n,sizeof(double));
		_memset(&dLam_UB,pdbu,n,sizeof(double));
		// failure handling
		_Set_NagFail(&fail);
		//// disable nag print
		_Set_NagOption(&options);

		nag_opt_bounds_no_deriv(
			n,
			_OBJFUN_LAMTS_BOUNDS,
			Nag_Bounds,
			pdbl,
			pdbu,
			pdLam_InitGuess_Begin,
			&df,
			memset(_alloca(sizeof(double)*n),0,sizeof(double)*n),
			&options,
			&commParams,
			&fail
			);
		
		if(fail.code != NE_NOERROR && fail.code!=NE_USER_STOP) 
			return fail.message;
	}

	memcpy(pdLam_InitGuess_Begin,Comm.m_pdx,sizeof(double)*n);

#ifdef _DEBUG
	_print_vector(pdLam_InitGuess_Begin,n);
	_print_vector(&df,1);
	_print_vector(Comm.m_pOptionPrice_Begin,Comm.m_pOption_End-Comm.m_pOption_Begin);
	_print_vector(Comm.m_pdf,Comm.m_pOption_End-Comm.m_pOption_Begin);
#endif //_DEBUG

	return Comm.m_szErr;
}


static const char *_Calibrate_TS(
	int nNumPointT,
	int nNumPointX,
	const _GenMidAt *pOption_Begin,
	const _GenMidAt *pOption_End,
	const _PCQ_Seq *pNumeraire_Begin,
	const double *pdLam_Begin,
	const _Calib_LamTS_DiagEuro_Comm* pCalib_LamTS_DiagEuro_Begin,
	// model
	const double *pdLamT_Begin,
	const double *pdLamT_End,
	const double *pdAlpha,
	const double *pdGamma,
	const double *pdRho,
	/// output
	double *pOptionPrice_Begin, // market option prices
	double *pdLam_InitGuess_Begin // lambda ts
	)
{
	const char *szErr = 0;
	const int n = pdLamT_End-pdLamT_Begin;// # of unknowns
	const int m = pOption_End-pOption_Begin;// # of equations

	// result of least squares optimization
	double df= 0.; 
		
	// initialize common params
	const _Calib_LamTS_Comm Comm = {
	nNumPointT,
	nNumPointX,
	pOption_Begin,
	pOption_End,
	pOptionPrice_Begin,// market price
	pNumeraire_Begin,
	memset(_alloca(sizeof(double)*m),0,sizeof(double)*m),
	pCalib_LamTS_DiagEuro_Begin,
	pdAlpha,
	pdGamma,
	pdRho,
	(double*)pdLamT_Begin,
	(double*)pdLamT_End,
	0,//memset(_alloca(225),0,225)
	memset(_alloca(sizeof(double)*n),0,sizeof(double)*n)
	};
		
	if(m<n) return "_Calibrate_LamTS(...) : over-deterministic system!";

	/// delegate to nag solver .... 
	{
		int nI=0;

		/// specify boundaries
		const double dLam_LB = -0.3;
		const double dLam_UB = 0.3;
		double *pdbl = _alloca(n*sizeof(double));
		double *pdbu = _alloca(n*sizeof(double));
		// nag failure handling
		NagError fail;
		// nag options
		Nag_E04_Opt options;
		// nag common params
		Nag_Comm commParams;
		// attach data to nag common parameters
		commParams.p = (_Calib_LamTS_Comm*) (&Comm);		
		_memset(&dLam_LB,pdbl,n,sizeof(double));
		_memset(&dLam_UB,pdbu,n,sizeof(double));
		// failure handling
		_Set_NagFail(&fail);
		//// disable nag print
		_Set_NagOption(&options);

		// market option prices
		for(nI=0;nI<m;++nI)
		{
			const char *szErr = _OBJFUN_LAMTS_BOUNDS_AUTOCAL(
				pCalib_LamTS_DiagEuro_Begin+nI,
				pOption_Begin+nI,
				pNumeraire_Begin+nI,
				pdLamT_Begin+nI,
				pdLamT_Begin+nI+1,
				pdLam_Begin+nI, // unknowns,i.e. lambda
				pdGamma,
				pdAlpha,
				pdRho,
				nNumPointT,
				nNumPointX,
				// result
				pOptionPrice_Begin+nI
				);
		
			_ASSERTE(*(pOptionPrice_Begin+nI)!=0.);

		}

		nag_opt_bounds_no_deriv(
			n,
			_OBJFUN_LAMTS_BOUNDS,
			Nag_Bounds,
			pdbl,
			pdbu,
			pdLam_InitGuess_Begin,
			&df,
			memset(_alloca(sizeof(double)*n),0,sizeof(double)*n),
			&options,
			&commParams,
			&fail
			);
		
		if(fail.code != NE_NOERROR && fail.code!=NE_USER_STOP) 
			return fail.message;
	}

	memcpy(pdLam_InitGuess_Begin,Comm.m_pdx,sizeof(double)*n);

#ifdef _DEBUG
	_print_vector(pdLam_InitGuess_Begin,n);
	_print_vector(&df,1);
	_print_vector(Comm.m_pOptionPrice_Begin,Comm.m_pOption_End-Comm.m_pOption_Begin);
	_print_vector(Comm.m_pdf,Comm.m_pOption_End-Comm.m_pOption_Begin);
#endif //_DEBUG

	return Comm.m_szErr;
}

const char *Calibrate_TS_Diag_Euro_MidAt(
	// Grid
	int nNumPointT,
	int nNumPointX,
	 // market
	const char *szYC,
	const char *szVC,
	Err (*get_cash_vol)(char*,double,double,double,int,char	*,double*,double*),	//	functor to get cash vol from the market
	long lToday,
	const char *szFreq,
	const char *szBasis,
	const char *szRefRate,
	long lTheoEnd,
	int nOnefequi,
	// model
	const double *pdLambda_Begin,// lambda from the matrices
	const double *pdAlpha,
	const double *pdGamma,
	const double *pdRho,
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
	// PayRec
	double dPayRec,
	const cpd_diag_calib_param *pParam_Swaption,
	const cpd_diag_calib_param *pParam_Caplet,
	const diag_calib_lm_params *pCalibParam,
	// results
	double *pdOptionPV_Begin,
	double *pdCoupon_Begin,
	double *pdLam_InitGuess_Begin // term struture returned
	)
{
	const int m = plEx_End-plEx_Begin;
	
	/// options
	_GenMidAt *pOption_Begin = _alloca(m*sizeof(_GenMidAt)),*pOption_End=pOption_Begin+m, *pOption=pOption_Begin;
	
	/// numeraires
	_PCQ_Seq Numeraire={0},*pNumeraire_Begin=_alloca(m*sizeof(_PCQ_Seq));
	// set lambda T to be exercise dates
	const double *pdLamT_Begin =_fill_tenor(lToday,plEx_Begin,plEx_End,_alloca(m*sizeof(double)));
	// diag european swaption pricing
	const _Calib_LamTS_DiagEuro_Comm *pCalib_LamTS_DiagEuro_Begin=_alloca(sizeof(_Calib_LamTS_DiagEuro_Comm)*m); 
	const _Calib_LamTS_DiagEuro_Comm *pCalib_LamTS_DiagEuro_End = pCalib_LamTS_DiagEuro_Begin+m;
	_Calib_LamTS_DiagEuro_Comm *pCalib_LamTS_DiagEuro = (_Calib_LamTS_DiagEuro_Comm*)pCalib_LamTS_DiagEuro_Begin;

	const int nSize_Fix = plFixPay_End-plFixPay_Begin;
	const int nSize_Flt = plFltPay_End-plFltPay_Begin;

	// ATM rate 
	const char *szErr = _fill_swapcashrate_diag(szYC,plExStart_Begin,plExStart_Begin+m,lTheoEnd,szFreq,szBasis,szRefRate,pdCoupon_Begin);
	if(szErr)return szErr;
	
	
	// because we have a co-terminal system, all options share the same terminal measure numeraire
	if(szErr =_alloc_Seq(1,&Numeraire,0,0)) return szErr;
	{
		 const long lFltPayBack = *(plFltPay_End-1);
		 const long lFixPayBack = *(plFixPay_End-1);
		 const long lLastPay = lFixPayBack>lFltPayBack?lFixPayBack:lFltPayBack;		
		 _ASSERTE(lLastPay>lToday);
		 _ASSERTE(lLastPay>=*(plEx_End-1));
		 _preinit_Seq(szYC,lToday,&lLastPay,1+&lLastPay,&Numeraire);
		 _memset(&Numeraire,pNumeraire_Begin,m,sizeof(Numeraire));
	}
	

	// allocation option strucutures 
	for(;plEx_Begin<plEx_End;
	++plEx_Begin,
	++plExStart_Begin,
	++plExFee_Begin,
	++pdCoupon_Begin,
	++pdExFee_Begin,
	++pOption,
	++pCalib_LamTS_DiagEuro
	)
	{
		const char *szErr = "GenMidAt_AutoCal_Diag(...): internal dates error!";

		// locate indeces
		const long *plFixPay_nI_Begin = (const long*)find_min_greater_than(plExStart_Begin,plFixPay_Begin,nSize_Fix,sizeof(long),lless);
		const long *plFltPay_nI_Begin = (const long*)find_min_greater_than(plExStart_Begin,plFltPay_Begin,nSize_Flt,sizeof(long),lless);
		
		const int nOffSet_Fix = plFixPay_nI_Begin-plFixPay_Begin;
		const int nOffSet_Flt = plFltPay_nI_Begin -plFltPay_Begin;
		const int nSize_Fix_nI = plFixPay_End-plFixPay_nI_Begin;
		const int nSize_Flt_nI = plFltPay_End-plFltPay_nI_Begin;

		const double dOne=1.;

		// check that indeces are located correctly
		if(nOffSet_Fix== nSize_Fix) return szErr;
		if(nOffSet_Flt== nSize_Flt) return szErr;

		_ASSERTE(nOffSet_Fix<nSize_Fix&&nOffSet_Fix>=0);
		_ASSERTE(nOffSet_Flt<nSize_Flt&&nOffSet_Flt>=0);
		_ASSERTE(nOffSet_Fix+nSize_Fix_nI==nSize_Fix);
		_ASSERTE(nOffSet_Flt+nSize_Flt_nI==nSize_Flt);

		// allocate Diagonal MidAts 
		if(szErr = _alloc_Option(plEx_End-plEx_Begin,nSize_Fix_nI,nSize_Flt_nI,pOption)) return szErr;

		// pre-initialize Diagonal MidAts 
		_preinit_Option(
			szYC,
			lToday,
			plEx_Begin,
			plEx_End,
			plExFee_Begin,
			pdExFee_Begin,
			plExStart_Begin,
			plFixPay_Begin+nOffSet_Fix,
			plFixPay_End,
			plFixStart_Begin+nOffSet_Fix,
			plFixEnd_Begin+nOffSet_Fix,
			plFixFee_Begin+nOffSet_Fix,
			pdFixCoverage_Begin+nOffSet_Fix,
			pdFixNotional_Begin+nOffSet_Fix,
			(const double*)_memset(pdCoupon_Begin,_alloca(sizeof(double)*nSize_Fix_nI),nSize_Fix_nI,sizeof(double)),//pdFixNotional_Begin+nOffSet_Fix,
			pdFixFee_Begin+nOffSet_Fix,
			plFltPay_Begin+nOffSet_Flt,
			plFltPay_End,
			plFltStart_Begin+nOffSet_Flt,
			plFltEnd_Begin+nOffSet_Flt,
			pdFltCoverage_Begin+nOffSet_Flt,
			pdFltNotional_Begin+nOffSet_Flt,
			pdMargin_Begin+nOffSet_Flt,
			pdSpread_Begin+nOffSet_Flt,
			dPayRec,
			pOption
			);
		
			// european diag calibration
		_init_Calib_LamTS_DiagEuro_Comm(
			lToday,
			szYC,
			szVC,
			get_cash_vol,
			lTheoEnd,
			szRefRate,
			szFreq,
			szBasis,
			plEx_Begin,
			plEx_End,
			nOnefequi,
			
			pdCoupon_Begin,
			plFixStart_Begin+nOffSet_Fix,
			plFixStart_Begin+nOffSet_Fix+nSize_Fix_nI,
			plFixEnd_Begin+nOffSet_Fix,
			plFixPay_Begin+nOffSet_Fix,
			
			plFltStart_Begin+nOffSet_Flt,
			plFltStart_Begin+nOffSet_Flt+nSize_Flt_nI,
			plFltEnd_Begin+nOffSet_Flt,
			plFltPay_Begin+nOffSet_Flt,
			pdMargin_Begin+nOffSet_Flt,
			pdSpread_Begin+nOffSet_Flt,
			pParam_Swaption,
			pParam_Caplet,
			pCalibParam,
			pCalib_LamTS_DiagEuro
			);	
	}


	// delegate ...
	szErr = _Calibrate_TS(
		nNumPointT,
		nNumPointX,
		pOption_Begin,
		pOption_Begin+m,
		pNumeraire_Begin,
		pdLambda_Begin,
		pCalib_LamTS_DiagEuro_Begin,
		pdLamT_Begin,
		pdLamT_Begin+m,
		pdAlpha,
		pdGamma,
		pdRho,
		pdOptionPV_Begin,
		pdLam_InitGuess_Begin
		);

	/// free numeraire
	_free_Seq(&Numeraire);

	// free options
	for(pOption=pOption_Begin;pOption<pOption_End;++pOption)
		_free_Option(pOption);

	pCalib_LamTS_DiagEuro=(_Calib_LamTS_DiagEuro_Comm*)pCalib_LamTS_DiagEuro_Begin;
	for(;pCalib_LamTS_DiagEuro<pCalib_LamTS_DiagEuro_End;++pCalib_LamTS_DiagEuro)
		_free_Calib_LamTS_DiagEuro_Comm(pCalib_LamTS_DiagEuro);

	return szErr;
}
