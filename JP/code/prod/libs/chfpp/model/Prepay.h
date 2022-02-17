// Prepay.h : main header file for the PREPAY DLL
//

#if !defined(AFX_PREPAY_H__6A7521CB_E62B_11D2_9FE5_080009C16B4A__INCLUDED_)
#define AFX_PREPAY_H__6A7521CB_E62B_11D2_9FE5_080009C16B4A__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifndef __AFXWIN_H__
#error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"		// main symbols

#include "common.h"

//#define DEBUG_DUMP 1				// Uncomment to Create Dump.csv for Debugging.

//#define DUMMY_DLL  1				// Uncomment to Create a dummy dll that returns 13% CPR

// Defines Globals
#define	C_1M_LAG_MTG_RATE	0
#define	C_2M_LAG_MTG_RATE	1
#define	C_3M_LAG_MTG_RATE	2
#define	C_MIN_MTG_5Y	3
#define	C_1M_LAG_SW2Y	4
#define	C_2M_LAG_SW2Y	5
#define	C_3M_LAG_SW2Y	6
#define	C_1M_LAG_SW5Y	7
#define	C_2M_LAG_SW5Y	8
#define	C_3M_LAG_SW5Y	9
#define	C_1M_LAG_SW10Y	10
#define	C_2M_LAG_SW10Y	11
#define	C_3M_LAG_SW10Y	12
/*
#define	C_1M_LAG_VIX	13
#define	C_2M_LAG_VIX	14
#define	C_HOUSE_TURNOVER	15
*/
#define	C_MIN_MTG_LOW_CUM_COUNT 13
#define	C_MIN_MTG_COUNT 14
#define	C_SPOT_MTG_F30Y 15
#define	C_SPOT_MTG_F15Y 16
#define	C_SPOT_PMMS_F30Y 17
#define	C_SPOT_PMMS_F15Y 18
#define	C_BC_MTG_SPREAD 19
#define	C_HEL_MH_SPREAD 20
#define	C_DATA_MONTH    21

// Defines Model Constants
#define	FACTOR_30Y_SW2Y			0
#define	FACTOR_30Y_SW5Y			1
#define	FACTOR_30Y_SW10Y		2
#define	FACTOR_30Y_VOL			3
#define FACTOR_30Y_CONST		4
#define HT_TYPE					5
#define B_CURVE_OUT				6
#define B_USE_MORT_ARRAY		7
#define B_USE_CURVE_ARRAY		8
#define D_CPR_ADDON				9
#define SHORT_PP_MONTHS			10
#define REVERSION_START_TIME	11
#define REVERSION_END_TIME		12
#define	FACTOR_15Y_SW2Y			13
#define	FACTOR_15Y_SW5Y			14
#define	FACTOR_15Y_SW10Y		15
#define	FACTOR_15Y_VOL			16
#define FACTOR_15Y_CONST		17

// Defines Model Parameters
#define	NUM_MODEL_PARAMS 2
#define MPARAM_VIX		 0
#define MPARAM_HT		 1

// Defines Rate Override
#define OV_SWAP_2Y		 0
#define OV_SWAP_10Y		 1
#define OV_MTG_RATE		 2
#define OV_SWAP_10Y_NEW	 3

//Defines number of 15 Factor Prepay products
#define MAX_15FACTOR_PRODUCTS 3

// Defines FNMA 30Y
#define 	CF30Y_LAG2_INCENTIVE_MIN			0
#define 	CF30Y_LAG2_INCENTIVE_MAX			1
#define 	CF30Y_PICK_LOAN_AGE_MIN				2
#define 	CF30Y_PICK_LOAN_AGE_MAX				3
#define 	CF30Y_LOAN_AGE_MIN					4
#define 	CF30Y_LOAN_AGE_MAX					5
#define 	CF30Y_INCENTIVE_MIN					6
#define 	CF30Y_INCENTIVE_MAX					7
#define 	CF30Y_CONSTANT						8
#define 	CF30Y_LAG2_INCENTIVE				9
#define 	CF30Y_LAG1_INCENTIVE				10
#define 	CF30Y_LAG3_INCENTIVE				11
#define 	CF30Y_MAX_INCENTIVE					12
#define 	CF30Y_POOL_FACTOR					13
#define 	CF30Y_POOL_FACTOR_SQ				14
#define 	CF30Y_LOAN_AGE						15
#define 	CF30Y_LOAN_AGE_SQ					16
#define 	CF30Y_LOAN_AGE_LAG2_INCENTIVE		17
#define 	CF30Y_LOWEST_MORT_RATE_FLG			18
#define 	CF30Y_MORT_RATE_SPR_SWAP10			19
#define 	CF30Y_CURVE							20
#define 	CF30Y_VIX							21
#define 	CF30Y_HOUSING_TURN_OVER				22
#define 	CF30Y_MONTH1						23
#define 	CF30Y_MONTH2						24
#define 	CF30Y_MONTH3						25
#define 	CF30Y_MONTH4						26
#define 	CF30Y_MONTH5						27
#define 	CF30Y_MONTH6						28
#define 	CF30Y_MONTH7						29
#define 	CF30Y_MONTH8						30
#define 	CF30Y_MONTH9						31
#define 	CF30Y_MONTH10						32
#define 	CF30Y_MONTH11						33

#define 	CF30Y_MAX_PARAMS					34

// Defines for GNMA30Y
#define 	GF30Y_LAG2_INCENTIVE_MIN	0
#define 	GF30Y_LAG2_INCENTIVE_MAX	1
#define 	GF30Y_PICK_LOAN_AGE_MIN	2
#define 	GF30Y_PICK_LOAN_AGE_MAX	3
#define 	GF30Y_LOAN_AGE_MIN	4
#define 	GF30Y_LOAN_AGE_MAX	5
#define 	GF30Y_INCENTIVE_MIN	6
#define 	GF30Y_INCENTIVE_MAX	7
#define 	GF30Y_CONSTANT	8
#define 	GF30Y_LAG2_INCENTIVE	9
#define 	GF30Y_LAG1_INCENTIVE	10
#define 	GF30Y_LAG3_INCENTIVE	11
#define 	GF30Y_POOL_FACTOR	12
#define 	GF30Y_POOL_FACTOR_SQ	13
#define 	GF30Y_LOAN_AGE	14
#define 	GF30Y_LOAN_AGE_SQ	15
#define 	GF30Y_LOAN_AGE_FLG0_3	16
#define 	GF30Y_LOAN_AGE_FLG3_6	17
#define 	GF30Y_LOAN_AGE_FLG6_12	18
#define 	GF30Y_LOAN_AGE_LAG2_INCENTIVE	19
#define 	GF30Y_CURVE_LAG2_INCENTIVE	20
#define 	GF30Y_FACTOR_LAG2_INCENTIVE	21
#define 	GF30Y_LOWEST_MORT_RATE_FLG_0	22
#define 	GF30Y_LOWEST_MORT_RATE_FLG_25	23
#define 	GF30Y_LOWEST_MORT_RATE_FLG_50	24
#define 	GF30Y_MORT_RATE_SPR_SWAP10	25
#define 	GF30Y_CURVE	26
#define 	GF30Y_VIX	27
#define 	GF30Y_HOUSING_TURN_OVER	28
#define 	GF30Y_MONTH1	29
#define 	GF30Y_MONTH2	30
#define 	GF30Y_MONTH3	31
#define 	GF30Y_MONTH4	32
#define 	GF30Y_MONTH5	33
#define 	GF30Y_MONTH6	34
#define 	GF30Y_MONTH7	35
#define 	GF30Y_MONTH8	36
#define 	GF30Y_MONTH9	37
#define 	GF30Y_MONTH10	38
#define 	GF30Y_MONTH11	39
#define 	GF30Y_TREND	40
		
#define 	GF30Y_MAX_PARAMS	41

// Defines for FNMA15Y
#define 	CF15Y_LAG2_INCENTIVE_MIN	0
#define 	CF15Y_LAG2_INCENTIVE_MAX	1
#define 	CF15Y_PICK_LOAN_AGE_MIN	2
#define 	CF15Y_PICK_LOAN_AGE_MAX	3
#define 	CF15Y_LOAN_AGE_MIN	4
#define 	CF15Y_LOAN_AGE_MAX	5
#define 	CF15Y_INCENTIVE_MIN	6
#define 	CF15Y_INCENTIVE_MAX	7
#define 	CF15Y_CONSTANT	8
#define 	CF15Y_LAG2_INCENTIVE	9
#define 	CF15Y_LAG1_INCENTIVE	10
#define 	CF15Y_LAG3_INCENTIVE	11
#define 	CF15Y_POOL_FACTOR	12
#define 	CF15Y_POOL_FACTOR_SQ	13
#define 	CF15Y_LOAN_AGE	14
#define 	CF15Y_LOAN_AGE_SQ	15
#define 	CF15Y_LOAN_AGE_LAG2_INCENTIVE	16
#define 	CF15Y_LOWEST_MORT_RATE_FLG_0	17
#define 	CF15Y_LOWEST_MORT_RATE_FLG_50	18
#define 	CF15Y_MORT_RATE_SPR_SWAP10	19
#define 	CF15Y_CURVE	20
#define 	CF15Y_VIX	21
#define 	CF15Y_HOUSING_TURN_OVER	22
#define 	CF15Y_MONTH1	23
#define 	CF15Y_MONTH2	24
#define 	CF15Y_MONTH3	25
#define 	CF15Y_MONTH4	26
#define 	CF15Y_MONTH5	27
#define 	CF15Y_MONTH6	28
#define 	CF15Y_MONTH7	29
#define 	CF15Y_MONTH8	30
#define 	CF15Y_MONTH9	31
#define 	CF15Y_MONTH10	32
#define 	CF15Y_MONTH11	33
		
#define 	CF15Y_MAX_PARAMS	34

//MIAC Utility Functions
double  cvt_smm_to_cpr ( double mth_smm );
double  cvt_cpr_to_smm ( double ann_cpr );

typedef double PARRAY[CHF_MAX_FORWARD_RATES];

class chfError
{
public:
	virtual char* GetErrorText() {return "General Error"; }
};

class chfDLLVersionError: public chfError
{
public:
	virtual char* GetErrorText() { return "CHF Prepay: Wrong DLL Version"; }
};


typedef struct ROWDATA_TAG
{
	SQLINTEGER	lLoanNum;
	SQLINTEGER	ind_lLoanNum;
	SQLCHAR		szProduct[21];
	SQLINTEGER	ind_szProduct;
	SQL_TIMESTAMP_STRUCT	dtMaturityDate;
	SQLINTEGER	ind_dtMaturityDate;
	SQLDOUBLE	dBalance;
	SQLINTEGER	ind_dBalance;
} ROWDATA;

class CDBBase : public CObject
{
public:
	CDBBase() {}
	~CDBBase() {}
	virtual short ConnectDB();
	virtual void FreeHandles();
	virtual short LoadData() = 0;	
protected:
	CDatabase db;
	SQLHENV  henv;
	SQLHDBC  hdbc;
	SQLHSTMT hstmt;
	SQLRETURN  retcode;
	SQLCHAR  SqlState[6], Msg[SQL_MAX_MESSAGE_LENGTH];
	SQLINTEGER NativeError;
	SQLSMALLINT MsgLen;
	SQLRETURN  rc2;
};

typedef struct CONSTANTS_TAG
{
    SQLINTEGER  iID;
    SQLINTEGER  ind_iID;
    SQLCHAR     szName[51];
    SQLINTEGER  ind_iName;
    SQLDOUBLE   dValue;
    SQLINTEGER  ind_iValue;
} CONSTANTS;

class CConstants : public CDBBase
{
public:
    CConstants();
    ~CConstants();

    virtual short LoadData();
	virtual short UnLoad();
	
	double *GetConstantsPtr() { return m_pConstants; }
    
	CString		m_szSql;
	CONSTANTS   m_ConstantsRow;
	int			m_iNumConstants;
	double		*m_pConstants;
};

typedef struct MODELCONSTANTS_TAG
{
    SQLINTEGER  iID;
    SQLINTEGER  ind_iID;
    SQLCHAR     szName[51];
    SQLINTEGER  ind_iName;
    SQLDOUBLE   dValue;
    SQLINTEGER  ind_iValue;
} MODELCONSTANTS;

class CModelConstants : public CDBBase
{
public:
    CModelConstants();
    ~CModelConstants();

    virtual short LoadData();
	virtual short UnLoad();

	double *GetConstantsPtr() { return m_pModelConstants; }
    CString m_szSql;
	MODELCONSTANTS   m_ModelConstantsRow;
	int m_iNumModelConstants;
	double *m_pModelConstants;
};

typedef struct MODELPARAMETERS_TAG
{
    SQLINTEGER  iparameter_id;
    SQLINTEGER  ind_iparameter_id;
    SQLINTEGER  imonth;
    SQLINTEGER  ind_imonth;
    SQLCHAR     szparameter[11];
    SQLINTEGER  ind_iparameter;
    SQLDOUBLE   dvalue;
    SQLINTEGER  ind_ivalue;
} MODELPARAMETERS;

class CModelParameters : public CDBBase
{
public:
    CModelParameters();
    ~CModelParameters();
	
	PARRAY * GetPtr() { return m_dModelParameters; }
    virtual short LoadData();
	virtual short UnLoad();

    CString m_szSql;
	MODELPARAMETERS   ModelParametersRow;

	PARRAY *m_dModelParameters;
};

typedef struct LRATESMOOTHPARAM_TAG
{
    SQLINTEGER  iProduct_id;
	SQLDOUBLE   dDensity;
	SQLDOUBLE   dBeta1;
 	SQLDOUBLE   dAlpha1;
	SQLDOUBLE   dBeta2;
	SQLDOUBLE   dAlpha2;
  
} LRATESMOOTHPARAM;


class CLowRateSmoothingParam : public CDBBase
{
public:
    CLowRateSmoothingParam();
    ~CLowRateSmoothingParam();
	
	LRATESMOOTHPARAM * GetSmoothParam(int nProductID);
    virtual short LoadData();
	virtual short UnLoad();

    CString m_szSql;
	LRATESMOOTHPARAM   LRateSmoothParamRow;

	LRATESMOOTHPARAM m_LRateSmoothParam[MAX_15FACTOR_PRODUCTS];
};

typedef struct TRANCHEDATA_TAG
{
    SQLCHAR     szid[51];
    SQLINTEGER  ind_iid;
    SQLDOUBLE   dOriginalBalance;
    SQLINTEGER  ind_iOriginalBalance;
} TRANCHEDATA;

class CTranche : public CObject
{
public:
	CTranche(double dOrigBal)
		{ dOriginalBalance = dOrigBal; }

	double dOriginalBalance;
};

class CTrancheData : public CDBBase
{
public:
    CTrancheData();
    ~CTrancheData();
	
	virtual short UnLoad();
    virtual short LoadData();
    int m_iNumTranches;

	CTranche *GetTrancheData(CString szTrancheID);

    TRANCHEDATA   TrancheDataRow;
	
	CMapStringToOb m_TrancheMap;
};

typedef struct MODELCOEFFS_TAG
{
    SQLINTEGER  iid;
    SQLINTEGER  ind_iid;
    SQLINTEGER  iproduct_id;
    SQLINTEGER  ind_iproduct_id;
    SQLINTEGER  imodel_id;
    SQLINTEGER  ind_imodel_id;
    SQLINTEGER  iparameter_id;
    SQLINTEGER  ind_iparameter_id;
    SQLDOUBLE   dvalue;
    SQLINTEGER  ind_ivalue;
} MODELCOEFFS;

class CModelCoeffs : public CDBBase
{
public:
    CModelCoeffs();
    ~CModelCoeffs();
	
	virtual short UnLoad();
    virtual short LoadData();
	double * GetCoeffsPtr(int nProductType) 
	{ return pBaseCoeffs[iMapProducts[nProductType]]; }
	
    CString m_szSql;

	int			m_iNumCoeffs;
	int			m_iNumProducts;
	double		*m_pCoeffs;
	int			iMapProducts[PP_MAX_ENUM_LOANTYPE];
	double		*pBaseCoeffs[PP_MAX_ENUM_LOANTYPE];
    MODELCOEFFS	ModelCoeffRow;
};

typedef struct RATEOVEREIDE_TAG
{
    SQLINTEGER  iparameter_id;
    SQLINTEGER  ind_iparameter_id;
    SQLINTEGER  imonth;
    SQLINTEGER  ind_imonth;
    SQLCHAR     szparameter[11];
    SQLINTEGER  ind_iparameter;
    SQLDOUBLE   dvalue;
    SQLINTEGER  ind_ivalue;
} RATEOVEREIDE;

class CRateOverride : public CDBBase
{
public:
    CRateOverride();
    ~CRateOverride();
	
	PARRAY * GetPtr() { return m_dRateOverride; }
    virtual short LoadData();
	virtual short UnLoad();

    CString m_szSql;
	RATEOVEREIDE   RateOverrideRow;
	int m_iNumRates;
	PARRAY *m_dRateOverride;
};

typedef struct PPTABLES_TAG
{
    SQLSMALLINT  sID;
    SQLINTEGER  ind_iID;
    SQLINTEGER  iTerm;
    SQLINTEGER  ind_iTerm;
    SQLDOUBLE   dIncentive;
    SQLINTEGER  ind_iIncentive;
    SQLDOUBLE   dCPR;
    SQLINTEGER  ind_iCPR;
} PPTABLE_ROW;

class CPPTables : public CDBBase
{
public:
    CPPTables();
    ~CPPTables();

    virtual short LoadData();
	virtual short UnLoad() ;
	virtual double GetCPR(int ePPTable, float dIncentive, float dAge);
	int RowIndexCount() { return PP_MAXINCENTIVE; }
	int ColIndexCount() { return PP_MAXTERM; }
	int RowIndexLookup(int idxTable, int i) { return (int)PPTablesArray[idxTable][i+1][0]; }
	int ColIndexLookup(int idxTable, int i) { return (int)PPTablesArray[idxTable][0][i+1]; }
    double CPRLookup(int idxTable, int row, int col) { return PPTablesArray[idxTable][row+1][col+1]; }

	CString m_szSql;
    int m_iNumRows;
	int m_iNumTables;
    int m_iNumRowsFetched;
    PPTABLE_ROW    PPTableRow;
	int iMapPPTable[PP_MAX_ENUM_LOANTYPE];
	float (*PPTablesArray)[PP_MAXINCENTIVE+1][PP_MAXTERM+1];
};

class CForecast : public CObject
{
public:
	CForecast();
	~CForecast();
	
	CHF_RATE_STRUCTURE m_stFwdRates[CHF_MAX_FORWARD_RATES+1];
	double chf_dSwap_2y[CHF_MAX_FORWARD_RATES+1], chf_dSwap_10y[CHF_MAX_FORWARD_RATES+1];
	double chf_dMortgage_rate[CHF_MAX_FORWARD_RATES+1];
	double chf_dMortgage_rate_15Y[CHF_MAX_FORWARD_RATES+1];

	CHF_PREPAYMENTMODELSTRUCTURE m_stLocalParams;
	CHF_PREPAYMENTMODELSTRUCTURE *m_lpPrepayParams;
	
	CConstants *m_pConstantsObj;
	CModelConstants *m_pModelConstantsObj;
	CModelParameters *m_pModelParametersObj;
	CRateOverride *m_pRateOverrideObj;
	CModelCoeffs *m_pModelCoeffsObj;
	CPPTables *m_pPPTablesObj;
	CTrancheData *m_pTrancheDataObj;
	CTranche *pTranche;
	CLowRateSmoothingParam *m_pLowRateSmoothingParamObj;

	double *m_pCoeffs;
	double *m_pConstants;
	double *m_pModelConstants;
	PARRAY *m_pModelParameters;
	PARRAY *m_pRateOverride;
	LRATESMOOTHPARAM *m_pLRateSmoothParam;

	int m_bReloadRates;

	CTime m_dtValuationMonth;
	
	CTime AddMonths(CTime dt, long lMonths);

	void Load();
	double CalcPrepay();
	
	double GetVIX(int iSimMonth);
	double GetHousingTurnover(int iSimMonth);

	// For FNMA 30 Year
	double CalcPrepayFNMA30(int iSimMonth);

	// For GNMA 30 Year
	double CalcPrepayGNMA30(int iSimMonth);

	// For FNMA 15 Year
	double CalcPrepayFNMA15(int iSimMonth);

	void Unload();
	BYTE m_bLoaded;
private:
	double m_dNextWAC;
	double m_dNewLow;
	int m_nLowCumCount;
	int m_nLow2Now;
	double lastSw2Y, lastSw5Y, lastSw10Y, lastMtgRate;

	CString szPrevTranche_ID;

	double Max2(double x1, double x2);
	double Min2(double x1, double x2);
	double Max3(double x1, double x2, double x3);
	double Min3(double x1, double x2, double x3);
	double Limit(double x, double low, double high);
	
	void SetMortgageOverride();
	void SetSwapOverride();
	void DumpRates();
	void NewLoanAdjustment(double& mtgLag1,double& mtgLag2,double& mtgLag3,double loan_age, double dMortgage);
	void ReloadRates();
};

#define MAX_ROWS 5

/////////////////////////////////////////////////////////////////////////////
// CPrepayApp
// See Prepay.cpp for the implementation of this class
//

class CPrepayApp : public CWinApp
{
public:
	CPrepayApp();

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CPrepayApp)
	public:
	virtual BOOL InitInstance();
	//}}AFX_VIRTUAL

	//{{AFX_MSG(CPrepayApp)
		// NOTE - the ClassWizard will add and remove member functions here.
		//    DO NOT EDIT what you see in these blocks of generated code !
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};


/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_PREPAY_H__6A7521CB_E62B_11D2_9FE5_080009C16B4A__INCLUDED_)
