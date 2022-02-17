// Prepay.cpp : Defines the initialization routines for the DLL.
//

#include "stdafx.h"

#include <math.h>
#include "Prepay.h"
#include "Globals.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

CForecast ForecastObj;

//
//	Note!
//
//		If this DLL is dynamically linked against the MFC
//		DLLs, any functions exported from this DLL which
//		call into MFC must have the AFX_MANAGE_STATE macro
//		added at the very beginning of the function.
//
//		For example:
//
//		extern "C" BOOL PASCAL EXPORT ExportedFunction()
//		{
//			AFX_MANAGE_STATE(AfxGetStaticModuleState());
//			// normal function body here
//		}
//
//		It is very important that this macro appear in each
//		function, prior to any calls into MFC.  This means that
//		it must appear as the first statement within the 
//		function, even before any object variable declarations
//		as their constructors may generate calls into the MFC
//		DLL.
//
//		Please see MFC Technical Notes 33 and 58 for additional
//		details.
//

/////////////////////////////////////////////////////////////////////////////
// CPrepayApp

BEGIN_MESSAGE_MAP(CPrepayApp, CWinApp)
	//{{AFX_MSG_MAP(CPrepayApp)
		// NOTE - the ClassWizard will add and remove mapping macros here.
		//    DO NOT EDIT what you see in these blocks of generated code!
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CPrepayApp construction

CPrepayApp::CPrepayApp()
{
	// TODO: add construction code here,
	// Place all significant initialization in InitInstance
}

/////////////////////////////////////////////////////////////////////////////
// The one and only CPrepayApp object

CPrepayApp theApp;

/////////////////////////////////////////////////////////////////////////////
// CPrepayApp initialization

BOOL CPrepayApp::InitInstance()
{
	// Register all OLE server (factories) as running.  This enables the
	//  OLE libraries to create objects from other applications.
	COleObjectFactory::RegisterAll();

	return TRUE;
}

/////////////////////////////////////////////////////////////////////////////
// Special entry points required for inproc servers

STDAPI DllGetClassObject(REFCLSID rclsid, REFIID riid, LPVOID* ppv)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState());
	return AfxDllGetClassObject(rclsid, riid, ppv);
}

STDAPI DllCanUnloadNow(void)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState());
	return AfxDllCanUnloadNow();
}

// by exporting DllRegisterServer, you can use regsvr.exe
STDAPI DllRegisterServer(void)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState());
	COleObjectFactory::UpdateRegistryAll();
	return S_OK;
}


short CDBBase::ConnectDB()
{
	CString dbconn = "PrepayDB;"; // PrepayDB;

	if(szGblDBName.GetLength() > 0)
		dbconn += "DBQ=" + szGblDBName + ";";

  retcode = SQLAllocHandle(SQL_HANDLE_ENV, SQL_NULL_HANDLE, &henv);  
	if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
		goto connect_error;

	retcode = SQLSetEnvAttr(henv, SQL_ATTR_ODBC_VERSION, (void*)SQL_OV_ODBC3, 0);
	if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
		goto connect_error;

	if (!db.Open(dbconn, false, false, "ODBC;", false))
        return 1;

	hdbc = db.m_hdbc;

	if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
		goto connect_error;
	retcode = SQLAllocHandle(SQL_HANDLE_STMT, hdbc, &hstmt); 
	if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
		goto connect_error;
	
	return 0;
connect_error:
	FreeHandles();
	return 1;
}

void CDBBase::FreeHandles()
{
	if(hstmt)
		SQLFreeHandle(SQL_HANDLE_STMT, hstmt);
	hstmt = NULL;
	if(hdbc)
	{
		db.Close();
		hdbc = NULL;
	}
	if(henv)
		SQLFreeHandle(SQL_HANDLE_ENV, henv);
	henv = NULL;
}


CForecast::CForecast() 
{
	m_bLoaded = FALSE;
	m_bReloadRates = 0;
	memset(&m_stLocalParams, 0, sizeof(CHF_PREPAYMENTMODELSTRUCTURE));
	m_lpPrepayParams = &m_stLocalParams;
}

CForecast::~CForecast() 
{
	Unload();
}

void CForecast::Load()
{
	if(!m_bLoaded)
	{
		TRY
		{
			m_pConstantsObj = new CConstants();
			if(m_pConstantsObj->LoadData())
			{
				MessageBox(NULL, "Unable to Load Constants", "Prepay.DLL Error Message", MB_OK);
				delete m_pConstantsObj;
				return;
			}	
			m_pConstants = m_pConstantsObj->GetConstantsPtr();
			
			m_pModelConstantsObj = new CModelConstants();
			if(m_pModelConstantsObj->LoadData())
			{
				MessageBox(NULL, "Unable to Load Model Constants", "Prepay.DLL Error Message", MB_OK);
				return;
			}	
			m_pModelConstants = m_pModelConstantsObj->GetConstantsPtr();

			m_pModelParametersObj = new CModelParameters();
			if(m_pModelParametersObj->LoadData())
			{
				MessageBox(NULL, "Unable to Load Model Parameters", "Prepay.DLL Error Message", MB_OK);
				return;
			}	
			m_pModelParameters =  m_pModelParametersObj->GetPtr();

			m_pModelCoeffsObj = new CModelCoeffs();
			if(m_pModelCoeffsObj->LoadData())
			{
				MessageBox(NULL, "Unable to Load Coefficients", "Prepay.DLL Error Message", MB_OK);
				return;
			}

			m_pPPTablesObj = new CPPTables();
			if(m_pPPTablesObj->LoadData())
			{
				MessageBox(NULL, "Unable to Load Prepay Tables", "Prepay.DLL Error Message", MB_OK);
				return;
			}	
			m_pTrancheDataObj = new CTrancheData();
			if(m_pTrancheDataObj->LoadData())
			{
				MessageBox(NULL, "Unable to Load Tranche Level Data", "Prepay.DLL Error Message", MB_OK);
				return;
			}	

			m_pRateOverrideObj = new CRateOverride();
			if(m_pRateOverrideObj->LoadData())
			{
				MessageBox(NULL, "Unable to Load Override Rates", "Prepay.DLL Error Message", MB_OK);
				return;
			}	
			m_pRateOverride =  m_pRateOverrideObj->GetPtr();

			m_pLowRateSmoothingParamObj = new CLowRateSmoothingParam();
			if(m_pLowRateSmoothingParamObj->LoadData())
			{
				MessageBox(NULL, "Unable to Load Low Rate Smoothing Parameters", "Prepay.DLL Error Message", MB_OK);
				return;
			}
			
			if(m_pModelConstantsObj->m_iNumModelConstants <= REVERSION_END_TIME)
			{
				MessageBox(NULL, "Invalid Prepay Database Version", "Prepay.DLL Error Message", MB_OK);
				m_bLoaded = TRUE;
				Unload();
				return;
			}
		}
		CATCH(CException, pEx)
		{
			TCHAR   szCause[255];
			pEx->GetErrorMessage(szCause, 255);
			MessageBox(NULL, szCause, "Prepay.DLL Error Message", MB_OK);
		}
		END_CATCH

		m_bLoaded = TRUE;
	}
}

void CForecast::Unload()
{
	if(m_bLoaded)
	{
		m_pConstantsObj->UnLoad();
		m_pModelConstantsObj->UnLoad();
		m_pModelParametersObj->UnLoad();
		m_pModelCoeffsObj->UnLoad();
		m_pPPTablesObj->UnLoad();
		m_pTrancheDataObj->UnLoad();
		m_pRateOverrideObj->UnLoad();
		m_pLowRateSmoothingParamObj->UnLoad();
			
		delete m_pConstantsObj;
		delete m_pModelConstantsObj;
		delete m_pModelParametersObj;
		delete m_pModelCoeffsObj;
		delete m_pPPTablesObj;
		delete m_pTrancheDataObj;
		delete m_pRateOverrideObj;
		delete m_pLowRateSmoothingParamObj;

		m_bLoaded = FALSE;
	}
}

CTime CForecast::AddMonths(CTime dt, long lMonths)
{
	long mm,dd,yy;

	yy = dt.GetYear();
	mm=dt.GetMonth();
	dd=dt.GetDay();
	mm += lMonths-1;
	yy+=mm/12;
	mm=(long)fmod(mm,12)+1;

	return(CTime(yy,mm,dd,0,0,0));
}

double CForecast::Max2(double x1, double x2)
{
	return (x1>x2)?x1:x2;
}

double CForecast::Min2(double x1, double x2)
{
	return (x1<x2)?x1:x2;
}

double CForecast::Max3(double x1, double x2, double x3)
{
	double x = (x1>x2)?x1:x2;

	return (x>x3)?x:x3;
}

double CForecast::Min3(double x1, double x2, double x3)
{
	double x = (x1<x2)?x1:x2;

	return (x<x3)?x:x3;
}

double CForecast::Limit(double x, double low, double high)
{
	return (x>high)?high:(x<low)?low:x;
}

double CForecast::CalcPrepay()
{
	double dResult = 0;
	int iMonth_ind;
	int i;
	double dMortgageSpread=0;
	char* pszFileName = "c:\\Rates.csv";
	CFile cf;

	iMonth_ind = m_lpPrepayParams->simulation_year*12+m_lpPrepayParams->simulation_month-1;

	m_pLRateSmoothParam = m_pLowRateSmoothingParamObj->GetSmoothParam(m_lpPrepayParams->loan_type);

	if(ForecastObj.m_pRateOverride == 0)
	{
		ForecastObj.m_pModelConstants[B_USE_MORT_ARRAY] = 0;
		ForecastObj.m_pModelConstants[B_USE_CURVE_ARRAY] = 0;
	}

	i = iMonth_ind;

	// CHASE Mortgage Rate Model - Overwrites MIAC Mortgage Rate !!

	ReloadRates();

	switch (m_lpPrepayParams->loan_type)
	{
		// FNMA 30 Year
		case 0: 
			if(iMonth_ind < ForecastObj.m_pModelConstants[SHORT_PP_MONTHS])
			{
				m_pCoeffs = m_pModelCoeffsObj->GetCoeffsPtr((int) m_lpPrepayParams->loan_type);
				dResult = CalcPrepayFNMA30(iMonth_ind); 
			}
            else {
				dResult = m_pPPTablesObj->GetCPR(m_lpPrepayParams->loan_type, 
					(float) (m_lpPrepayParams->gross_wac - (chf_dMortgage_rate[iMonth_ind] + dMortgageSpread))*10000,
					(float) (m_lpPrepayParams->original_term - m_lpPrepayParams->remain_term))/100;
            }
			break;

		// GNMA 30 Year
		case 10: 
			if(iMonth_ind < ForecastObj.m_pModelConstants[SHORT_PP_MONTHS])
			{
				m_pCoeffs = m_pModelCoeffsObj->GetCoeffsPtr((int) m_lpPrepayParams->loan_type);
				dResult = CalcPrepayGNMA30(iMonth_ind); 
			}
			else
				dResult = m_pPPTablesObj->GetCPR(m_lpPrepayParams->loan_type, 
					(float) (m_lpPrepayParams->gross_wac - (chf_dMortgage_rate[iMonth_ind] + dMortgageSpread))*10000,
					(float) (m_lpPrepayParams->original_term - m_lpPrepayParams->remain_term))/100;
			break;

		// FNMA 15 Year
		case 2: 
			//dMortgageSpread = ForecastObj.m_pConstants[C_SPOT_MTG_F15Y]-ForecastObj.m_pConstants[C_SPOT_MTG_F30Y]; 
			if(iMonth_ind < ForecastObj.m_pModelConstants[SHORT_PP_MONTHS])
			{
				m_pCoeffs = m_pModelCoeffsObj->GetCoeffsPtr((int) m_lpPrepayParams->loan_type);
				dResult = CalcPrepayFNMA15(iMonth_ind); 
			}
			else
				dResult = m_pPPTablesObj->GetCPR(m_lpPrepayParams->loan_type, 
					(float) (m_lpPrepayParams->gross_wac - (chf_dMortgage_rate_15Y[iMonth_ind] + dMortgageSpread))*10000,
					(float) (m_lpPrepayParams->original_term - m_lpPrepayParams->remain_term))/100;

			break;
		
		// All These are 15 Year Products 
		case 	3	: // CBL7YR
		case 	4	: // CBL5YR
		case	7	: // FH15YR
		case 	13	: // GF15YR
		case 	19	: // JF15YR
		case 	27	: // EF15YR
		case 	37	: // P715
				//dMortgageSpread = ForecastObj.m_pConstants[C_SPOT_MTG_F15Y]-ForecastObj.m_pConstants[C_SPOT_MTG_F30Y]; 
				dResult = m_pPPTablesObj->GetCPR(m_lpPrepayParams->loan_type, 
					(float) (m_lpPrepayParams->gross_wac - (chf_dMortgage_rate_15Y[iMonth_ind] + dMortgageSpread))*10000,
					(float) (m_lpPrepayParams->original_term - m_lpPrepayParams->remain_term))/100;
				break;

		case 49:
		case 50:
				dMortgageSpread = ForecastObj.m_pConstants[C_BC_MTG_SPREAD]; 
				dResult = m_pPPTablesObj->GetCPR(m_lpPrepayParams->loan_type, 
					(float) (m_lpPrepayParams->gross_wac - (chf_dMortgage_rate[iMonth_ind] + dMortgageSpread))*10000,
					(float) (m_lpPrepayParams->original_term - m_lpPrepayParams->remain_term))/100;
				break;

		case 55:
		case 58:
				dMortgageSpread = ForecastObj.m_pConstants[C_HEL_MH_SPREAD]; 
				dResult = m_pPPTablesObj->GetCPR(m_lpPrepayParams->loan_type, 
					(float) (m_lpPrepayParams->gross_wac - (chf_dMortgage_rate[iMonth_ind] + dMortgageSpread))*10000,
					(float) (m_lpPrepayParams->original_term - m_lpPrepayParams->remain_term))/100;
				break;

		default:

			dResult = m_pPPTablesObj->GetCPR(m_lpPrepayParams->loan_type, 
				(float) (m_lpPrepayParams->gross_wac - (chf_dMortgage_rate[iMonth_ind] + dMortgageSpread))*10000,
				(float) (m_lpPrepayParams->original_term - m_lpPrepayParams->remain_term))/100;
	}

//	dResult *= m_lpPrepayParams->state_factor;

	dResult += ForecastObj.m_pModelConstants[D_CPR_ADDON];

	if (dResult >= 0.99) dResult = 0.99;
	if (dResult <= 0) dResult = 0;
	
	return dResult;
}

/*
Function Name: CalcPrepayFNMA30
Class: CForecast

Purpose: Calculates CPR for FNMA 30 Year Product
*/
double CForecast::CalcPrepayFNMA30(int iSimMonth)
{
	double dCpr[4];
	double wac, incent, incent1, incent2, incent3, incentMax, incentMin;
	double factor, loan_age;
	double mtgLag1, mtgLag2, mtgLag3, sw10YLag2, sw2YLag2;
	int iCurMonth, nModel, i;
	double dParamsArray[CF30Y_MAX_PARAMS];
	double calib_incent_min, calib_incent_max;
	

#ifdef DEBUG_DUMP
	char* pszFileName = "c:\\DUMP.csv";
	CFile cf;
	CFileException fileException;

	if ( !cf.Open( pszFileName, CFile::modeCreate |    
			  CFile::modeWrite | CFile::modeNoTruncate , &fileException ))
	{
		TRACE( "Can't open file %s, error = %u\n",
		   pszFileName, fileException.m_cause );
	}
	cf.SeekToEnd();
#endif

	memset(dParamsArray, 0, CF30Y_MAX_PARAMS*sizeof(double));
	
	if(szPrevTranche_ID != m_lpPrepayParams->szTranche_id)
	{
		pTranche = m_pTrancheDataObj->GetTrancheData(m_lpPrepayParams->szTranche_id);
		szPrevTranche_ID = m_lpPrepayParams->szTranche_id;
	}

	if(iSimMonth == 0)
		m_dNextWAC = m_lpPrepayParams->gross_wac;

	wac = m_dNextWAC;

	loan_age = m_lpPrepayParams->wala; //+ iSimMonth;
	
	// Gets the LAG Rates (NOTE: Order of the Constants 
	// is used here to determine the 1m, 2m, ... LAG rates.)
	// HENCE DO NOT CHANGE THE ORDER OF THE GLOBAL CONSTANTS.

	if(iSimMonth-1<0) 
		mtgLag1 = m_pConstants[C_1M_LAG_MTG_RATE];
	else
		mtgLag1 = chf_dMortgage_rate[iSimMonth-1] + m_pConstants[C_SPOT_PMMS_F30Y] - m_pConstants[C_SPOT_MTG_F30Y];

	if(iSimMonth-2<0)
	{
		mtgLag2 = m_pConstants[C_1M_LAG_MTG_RATE+1-iSimMonth];
		sw2YLag2 = m_pConstants[C_1M_LAG_SW2Y+1-iSimMonth];
		sw10YLag2 = m_pConstants[C_1M_LAG_SW10Y+1-iSimMonth];
	}
	else
	{
        mtgLag2 = chf_dMortgage_rate[iSimMonth-2]  + m_pConstants[C_SPOT_PMMS_F30Y] - m_pConstants[C_SPOT_MTG_F30Y];
		sw2YLag2 = chf_dSwap_2y[iSimMonth-2];
		sw10YLag2 = chf_dSwap_10y[iSimMonth-2];
	}

    
	if(iSimMonth-3<0) 
		mtgLag3 = m_pConstants[C_1M_LAG_MTG_RATE+2-iSimMonth];
	else
		mtgLag3 = chf_dMortgage_rate[iSimMonth-3] + m_pConstants[C_SPOT_PMMS_F30Y] - m_pConstants[C_SPOT_MTG_F30Y];
	
    // bmc 
//     fprintf (stderr, "# %d %lf %lf %lf %lf %lf\n",iSimMonth,mtgLag1,mtgLag2,mtgLag3,sw2YLag2,sw10YLag2);
    if (mtgLag2 <= 0. || mtgLag3 <= 0. || mtgLag1 <= 0. || sw2YLag2 <= 0. || sw10YLag2 <= 0.) {
        throw "bad rate"; 
    }

        double x1 = chf_dMortgage_rate[iSimMonth];
    double x2 = m_pConstants[C_SPOT_PMMS_F30Y];
    double x3 = m_pConstants[C_SPOT_MTG_F30Y];
    
	NewLoanAdjustment(mtgLag1, mtgLag2, mtgLag3, loan_age, chf_dMortgage_rate[iSimMonth] + m_pConstants[C_SPOT_PMMS_F30Y] - m_pConstants[C_SPOT_MTG_F30Y]);

	incent1 = sqrt(wac/mtgLag1)-1;
	incent2 = sqrt(wac/mtgLag2)-1;
	incent3 = sqrt(wac/mtgLag3)-1;
	incentMax = sqrt(wac/Min3(mtgLag1, mtgLag2, mtgLag3)) - 1;
	incentMin = sqrt(wac/Max3(mtgLag1, mtgLag2, mtgLag3)) - 1;

//Model Pick
	double dModel[4];

/*	nModel = 0;
	if(incent2 >= m_pCoeffs[CF30Y_LAG2_INCENTIVE_MIN] && loan_age <= m_pCoeffs[CF30Y_PICK_LOAN_AGE_MAX])
		nModel = 0;

	if(incent2 < m_pCoeffs[CF30Y_MAX_PARAMS+CF30Y_LAG2_INCENTIVE_MAX] && loan_age <= m_pCoeffs[CF30Y_MAX_PARAMS+CF30Y_PICK_LOAN_AGE_MAX])
		nModel = 1;

	if(loan_age > m_pCoeffs[CF30Y_MAX_PARAMS*2+CF30Y_PICK_LOAN_AGE_MIN] && loan_age < m_pCoeffs[CF30Y_MAX_PARAMS*2+CF30Y_PICK_LOAN_AGE_MAX])
		nModel = 2;

	if(loan_age >= m_pCoeffs[CF30Y_MAX_PARAMS*3+CF30Y_PICK_LOAN_AGE_MIN])
		nModel = 3;
*/		
	double p,q,r;
	if(fabs(incent2)<0.01) 
		p = (incent2+0.01)/(0.01 + 0.01);
	else if(incent2 >= 0.01)
			p = 1;
		 else
		 	p = 0;

/*	if(abs(incent2)<0.01) 
		dModel[1] = (0.01 - incent2)/(0.01 + 0.01);
	else if(incent2 <= -0.01)
			dModel[1] = 1;
		else if(incent2 < m_pCoeffs[CF30Y_MAX_PARAMS+CF30Y_LAG2_INCENTIVE_MAX])
				dModel[1] = 1;
			else
				dModel[1] = 0;
*/
	if(loan_age > 96 && loan_age < (96+36)) 
		q = (loan_age-96)/36.0;
	else if(loan_age > m_pCoeffs[CF30Y_MAX_PARAMS*2+CF30Y_PICK_LOAN_AGE_MIN])
		q = 1;
	else 
		q = 0;

	if(loan_age >= 210 && loan_age < (210+24)) 
		r = (loan_age-210)/24.0;
	else if(loan_age >= (m_pCoeffs[CF30Y_MAX_PARAMS*3+CF30Y_PICK_LOAN_AGE_MIN] + 24))
		r = 1;
	else 
		r = 0;

	dModel[0] = (1-q)*(1-r)*p;
	dModel[1] = (1-q)*(1-r)*(1-p);
	dModel[2] = q*(1-r);
	dModel[3] = r;

	for(nModel = 0; nModel < 4; nModel++)
	{
		if(dModel[nModel] == 0) 
		{
			dCpr[nModel] = 0;
			continue;
		}

	//Calibration (sets upper and lower limits for incentive)
		calib_incent_min = m_pCoeffs[CF30Y_MAX_PARAMS*nModel+CF30Y_INCENTIVE_MIN];
		calib_incent_max = m_pCoeffs[CF30Y_MAX_PARAMS*nModel+CF30Y_INCENTIVE_MAX];

		if(loan_age <= 6)
			dParamsArray[CF30Y_LOAN_AGE_LAG2_INCENTIVE] = 6*incent2;
		else if(loan_age <= 50)
			dParamsArray[CF30Y_LOAN_AGE_LAG2_INCENTIVE] = loan_age*incent2;
		else
			dParamsArray[CF30Y_LOAN_AGE_LAG2_INCENTIVE] = (50+log(loan_age-49))*incent2;

		incent1 = Limit(incent1, calib_incent_min, calib_incent_max);
		incent2 = Limit(incent2, calib_incent_min, calib_incent_max);
		incent3 = Limit(incent3, calib_incent_min, calib_incent_max);
		incentMax = Limit(incentMax, calib_incent_min, calib_incent_max);
		incentMin = Limit(incentMin, calib_incent_min, calib_incent_max);


	//  SET UP PARAMETERS
		
		dParamsArray[CF30Y_CONSTANT] = 1;
		incent = dParamsArray[CF30Y_LAG2_INCENTIVE] = incent2;
		dParamsArray[CF30Y_LAG1_INCENTIVE] = incent1;
		dParamsArray[CF30Y_LAG3_INCENTIVE] = incent3;
		dParamsArray[CF30Y_MAX_INCENTIVE] = incentMax;

		if(pTranche)
			factor = m_lpPrepayParams->remain_balance/pTranche->dOriginalBalance;
		else
			factor = 0.99; // Default if the Tranche Lookup Fails.
			
		dParamsArray[CF30Y_POOL_FACTOR] = log(factor/(1-Min2(factor, 0.99)));
		dParamsArray[CF30Y_POOL_FACTOR_SQ] = pow((dParamsArray[CF30Y_POOL_FACTOR]-0.65),2);
		dParamsArray[CF30Y_LOAN_AGE] = loan_age;
		dParamsArray[CF30Y_LOAN_AGE_SQ] = pow(dParamsArray[CF30Y_LOAN_AGE]-50, 2);
		
		// m_dNewLow is the lowest Mortgage Rate in last 5 Years
		if(iSimMonth == 0) 
		{
			m_dNewLow = m_pConstants[C_MIN_MTG_5Y];
			m_nLowCumCount = (int) m_pConstants[C_MIN_MTG_LOW_CUM_COUNT];
			m_nLow2Now = (int) m_pConstants[C_MIN_MTG_COUNT];
		}

	//		m_dNewLow = m_pConstants[C_MIN_MTG_5Y];
		
	/*	if( Min2(mtgLag1, mtgLag2) - m_dNewLow <= 0.005)
			dParamsArray[CF30Y_LOWEST_MORT_RATE_FLG] = 1;
		else
			dParamsArray[CF30Y_LOWEST_MORT_RATE_FLG] = 0;
	*/
		double dNewLowM = (Min2(mtgLag1, mtgLag2) - m_dNewLow) * 100;

		dParamsArray[CF30Y_LOWEST_MORT_RATE_FLG] = ((dNewLowM<=m_pLRateSmoothParam->dAlpha2) ? 1:0) * 
			 m_pLRateSmoothParam->dDensity * 
			(2/(1+exp(m_pLRateSmoothParam->dBeta2 * (dNewLowM-m_pLRateSmoothParam->dAlpha2)))-1) +
			1/(1+exp(m_pLRateSmoothParam->dBeta1 * (dNewLowM-m_pLRateSmoothParam->dAlpha1)));

		if ( Min2(mtgLag1, mtgLag2) - m_dNewLow <= 0.000 )
		{
			m_nLowCumCount=1;
			m_nLow2Now=1;
		}
		else 
		{
			m_nLowCumCount=m_nLowCumCount + ((Min2(mtgLag1, mtgLag2) - m_dNewLow <= 0.005) ? 1:0);
			m_nLow2Now = m_nLow2Now + 1; 
		}
		
		if(m_nLow2Now >4 && m_nLowCumCount >= m_nLow2Now) dParamsArray[CF30Y_LOWEST_MORT_RATE_FLG] = 0;

		//update m_dNewLow
		m_dNewLow = Min2(mtgLag1, m_dNewLow);  

		// Here the rate is multiplied by 100 since MIAC interest rates are decimals.

		dParamsArray[CF30Y_MORT_RATE_SPR_SWAP10] = (mtgLag2 - sw10YLag2)*100;
		dParamsArray[CF30Y_CURVE] = (sw10YLag2 - sw2YLag2) * 100;
		
		//Currently FIXED
	//	dParamsArray[CF30Y_VIX] = m_pConstants[C_1M_LAG_VIX];
		dParamsArray[CF30Y_VIX] = GetVIX(iSimMonth);
		
		//Currently FIXED
	//	dParamsArray[CF30Y_HOUSING_TURN_OVER] = m_pConstants[C_HOUSE_TURNOVER];
		dParamsArray[CF30Y_HOUSING_TURN_OVER] = GetHousingTurnover(iSimMonth);

		iCurMonth = (iSimMonth+m_dtValuationMonth.GetMonth()-1) % 12;
		
		if(iCurMonth < 11)
			dParamsArray[CF30Y_MONTH1+iCurMonth] = 1;

		
		// CPR Formula:
		dCpr[nModel] = 0;

		
#ifdef DEBUG_DUMP
		double a, b, c;
		char szFileBuf[256];
		sprintf(szFileBuf, "%d,%f,%f,%f,%d,%d,%f,", iSimMonth, 
				m_stFwdRates[iSimMonth].dSwap_5y,
				chf_dMortgage_rate[iSimMonth],
				m_dNewLow,
				m_nLowCumCount,
				m_nLow2Now,
				m_lpPrepayParams->remain_balance);

		cf.Write(szFileBuf, lstrlen(szFileBuf));
#endif

		for(i=CF30Y_CONSTANT; i<CF30Y_MAX_PARAMS; i++)
		{

#ifdef DEBUG_DUMP
			a = dParamsArray[i];
			b = m_pCoeffs[CF30Y_MAX_PARAMS*nModel+i];

			sprintf(szFileBuf, "%d %f, %f, ", i, a, b);
			
			cf.Write(szFileBuf, lstrlen(szFileBuf));
			c = a * b;
#endif

			dCpr[nModel] += dParamsArray[i] * m_pCoeffs[CF30Y_MAX_PARAMS*nModel+i];
		}
	}
	
	double dCprFinal=0;


	for(nModel=0; nModel<4; nModel++)
	{
		dCprFinal += dModel[nModel]*dCpr[nModel];
	}
	
	if(dCprFinal <= -100.0)
		dCprFinal = 0;
	else
		dCprFinal = 1 / (1 + exp(-dCprFinal));

	//m_dNextWAC = Min2(-0.00077812+1.00009*wac*100-0.001128*dCpr, wac*100)/100;

#ifdef DEBUG_DUMP	
    char szFileBuf[256];
	sprintf(szFileBuf, "%f\r\n", dCprFinal);
	cf.Write(szFileBuf, lstrlen(szFileBuf));
	cf.Close();
#endif

	return dCprFinal;
}

/*
Function Name: CalcPrepayGNMA30
Class: CForecast

Purpose: Calculates CPR for GNMA 30 Year Product
*/
double CForecast::CalcPrepayGNMA30(int iSimMonth)
{
	double dCpr[4];
	double wac, incent, incent1, incent2, incent3, incentMax, incentMin;
	double factor, loan_age;
	double mtgLag1, mtgLag2, mtgLag3, sw10YLag2, sw2YLag2;
	int iCurMonth, nModel, i;
	double dParamsArray[GF30Y_MAX_PARAMS];
	double calib_incent_min, calib_incent_max;

#ifdef DEBUG_DUMP
	char* pszFileName = "c:\\DUMP.csv";
	CFile cf;
	CFileException fileException;

	if ( !cf.Open( pszFileName, CFile::modeCreate |    
			  CFile::modeWrite | CFile::modeNoTruncate ), &fileException )
	{
		TRACE( "Can't open file %s, error = %u\n",
		   pszFileName, fileException.m_cause );
	}
	cf.SeekToEnd();
#endif

	memset(dParamsArray, 0, GF30Y_MAX_PARAMS*sizeof(double));
	
	if(szPrevTranche_ID != m_lpPrepayParams->szTranche_id)
	{
		pTranche = m_pTrancheDataObj->GetTrancheData(m_lpPrepayParams->szTranche_id);
		szPrevTranche_ID = m_lpPrepayParams->szTranche_id;
	}

	if(iSimMonth == 0)
		m_dNextWAC = m_lpPrepayParams->gross_wac;

	wac = m_dNextWAC;

	loan_age = m_lpPrepayParams->wala; //+ iSimMonth;
	
	// Gets the LAG Rates (NOTE: Order of the Constants 
	// is used here to determine the 1m, 2m, ... LAG rates.)
	// HENCE DO NOT CHANGE THE ORDER OF THE GLOBAL CONSTANTS.
	if(iSimMonth-1<0) 
		mtgLag1 = m_pConstants[C_1M_LAG_MTG_RATE];
	else
		mtgLag1 = chf_dMortgage_rate[iSimMonth-1] + m_pConstants[C_SPOT_PMMS_F30Y] - m_pConstants[C_SPOT_MTG_F30Y];

	if(iSimMonth-2<0)
	{
		mtgLag2 = m_pConstants[C_1M_LAG_MTG_RATE+1-iSimMonth];
		sw2YLag2 = m_pConstants[C_1M_LAG_SW2Y+1-iSimMonth];
		sw10YLag2 = m_pConstants[C_1M_LAG_SW10Y+1-iSimMonth];
	}
	else
	{
		mtgLag2 = chf_dMortgage_rate[iSimMonth-2] + m_pConstants[C_SPOT_PMMS_F30Y] - m_pConstants[C_SPOT_MTG_F30Y];
		sw2YLag2 = chf_dSwap_2y[iSimMonth-2];
		sw10YLag2 = chf_dSwap_10y[iSimMonth-2];
	}

	if(iSimMonth-3<0) 
		mtgLag3 = m_pConstants[C_1M_LAG_MTG_RATE+2-iSimMonth];
	else
		mtgLag3 = chf_dMortgage_rate[iSimMonth-3] + m_pConstants[C_SPOT_PMMS_F30Y] - m_pConstants[C_SPOT_MTG_F30Y];

	NewLoanAdjustment(mtgLag1, mtgLag2, mtgLag3, loan_age, chf_dMortgage_rate[iSimMonth] + m_pConstants[C_SPOT_PMMS_F30Y] - m_pConstants[C_SPOT_MTG_F30Y]);

	incent1 = sqrt(wac/mtgLag1)-1;
	incent2 = sqrt(wac/mtgLag2)-1;
	incent3 = sqrt(wac/mtgLag3)-1;
	incentMax = sqrt(wac/Min3(mtgLag1, mtgLag2, mtgLag3)) - 1;
	incentMin = sqrt(wac/Max3(mtgLag1, mtgLag2, mtgLag3)) - 1;

//Model Pick
	double dModel[4];
/*	
	nModel = 0;
	if(incent2 >= m_pCoeffs[GF30Y_LAG2_INCENTIVE_MIN] && loan_age <= m_pCoeffs[GF30Y_PICK_LOAN_AGE_MAX])
		nModel = 0;

	if(incent2 < m_pCoeffs[GF30Y_MAX_PARAMS+GF30Y_LAG2_INCENTIVE_MAX] && loan_age <= m_pCoeffs[GF30Y_MAX_PARAMS+GF30Y_PICK_LOAN_AGE_MAX])
		nModel = 1;

	if(loan_age > m_pCoeffs[GF30Y_MAX_PARAMS*2+GF30Y_PICK_LOAN_AGE_MIN] && loan_age < m_pCoeffs[GF30Y_MAX_PARAMS*2+GF30Y_PICK_LOAN_AGE_MAX])
		nModel = 2;

	if(loan_age >= m_pCoeffs[GF30Y_MAX_PARAMS*3+GF30Y_PICK_LOAN_AGE_MIN])
		nModel = 3;
*/
	double p,q,r;

	if(fabs(incent2)<0.01) 
		p = (incent2+0.01)/(0.01 + 0.01);
	else if(incent2 >= 0.01)
			p = 1;
		else
			p = 0;

	if(loan_age > 96 && loan_age < (96+36)) 
		q = (loan_age-96)/36.0;
	else if(loan_age > m_pCoeffs[GF30Y_MAX_PARAMS*2+GF30Y_PICK_LOAN_AGE_MIN])
		q = 1;
	else 
		q = 0;

	if(loan_age >= 210 && loan_age < (210+24)) 
		r = (loan_age-210)/24.0;
	else if(loan_age >= (m_pCoeffs[GF30Y_MAX_PARAMS*3+GF30Y_PICK_LOAN_AGE_MIN] + 24))
		r = 1;
	else 
		r = 0;

	dModel[0] = (1-q)*(1-r)*p;
	dModel[1] = (1-q)*(1-r)*(1-p);
	dModel[2] = q*(1-r);
	dModel[3] = r;

	for(nModel=0; nModel<4; nModel++)
	{
		if(dModel[nModel] == 0) 
		{
			dCpr[nModel] = 0;
			continue;
		}

	//Calibration (sets upper and lower limits for incentive)
		calib_incent_min = m_pCoeffs[GF30Y_MAX_PARAMS*nModel+GF30Y_INCENTIVE_MIN];
		calib_incent_max = m_pCoeffs[GF30Y_MAX_PARAMS*nModel+GF30Y_INCENTIVE_MAX];

		if(loan_age <= 6)
			dParamsArray[GF30Y_LOAN_AGE_LAG2_INCENTIVE] = 6*incent2;
		else if(loan_age <= 50)
			dParamsArray[GF30Y_LOAN_AGE_LAG2_INCENTIVE] = loan_age*incent2;
		else
			dParamsArray[GF30Y_LOAN_AGE_LAG2_INCENTIVE] = (50+log(loan_age-49))*incent2;

		incent1 = Limit(incent1, calib_incent_min, calib_incent_max);
		incent2 = Limit(incent2, calib_incent_min, calib_incent_max);
		incent3 = Limit(incent3, calib_incent_min, calib_incent_max);
		incentMax = Limit(incentMax, calib_incent_min, calib_incent_max);
		incentMin = Limit(incentMin, calib_incent_min, calib_incent_max);


	//  SET UP PARAMETERS
		
		dParamsArray[GF30Y_CONSTANT] = 1;
		incent = dParamsArray[GF30Y_LAG2_INCENTIVE] = incent2;
		dParamsArray[GF30Y_LAG1_INCENTIVE] = incent1;
		dParamsArray[GF30Y_LAG3_INCENTIVE] = incent3;
		dParamsArray[GF30Y_INCENTIVE_MAX] = incentMax;

		if(pTranche)
			factor = m_lpPrepayParams->remain_balance/pTranche->dOriginalBalance;
		else
			factor = 0.99; // Default if the Tranche Lookup Fails.
			
		dParamsArray[GF30Y_POOL_FACTOR] = log(factor/(1-Min2(factor, 0.99)));
		dParamsArray[GF30Y_POOL_FACTOR_SQ] = pow((dParamsArray[GF30Y_POOL_FACTOR]-0.65),2);
		dParamsArray[GF30Y_LOAN_AGE] = loan_age;
		dParamsArray[GF30Y_LOAN_AGE_SQ] = pow(dParamsArray[GF30Y_LOAN_AGE]-42, 2);
		
		// m_dNewLow is the lowest Mortgage Rate in last 5 Years
		if(iSimMonth==0)
			m_dNewLow =  m_pConstants[C_MIN_MTG_5Y];
		
	/*	if( Min2(mtgLag1, mtgLag2) - m_dNewLow <= 0.000)
			dParamsArray[GF30Y_LOWEST_MORT_RATE_FLG_0] = 1;
		else
			dParamsArray[GF30Y_LOWEST_MORT_RATE_FLG_0] = 0;
		if( Min2(mtgLag1, mtgLag2) - m_dNewLow <= 0.0025 && Min2(mtgLag1, mtgLag2) - m_dNewLow > 0.000)
			dParamsArray[GF30Y_LOWEST_MORT_RATE_FLG_25] = 1;
		else
			dParamsArray[GF30Y_LOWEST_MORT_RATE_FLG_25] = 0;
		if( Min2(mtgLag1, mtgLag2) - m_dNewLow <= 0.005 && Min2(mtgLag1, mtgLag2) - m_dNewLow > 0.0025)
			dParamsArray[GF30Y_LOWEST_MORT_RATE_FLG_50] = 1;
		else
			dParamsArray[GF30Y_LOWEST_MORT_RATE_FLG_50] = 0;

	*/
		double dNewLowM = (Min2(mtgLag1, mtgLag2) - m_dNewLow) * 100;

		dParamsArray[GF30Y_LOWEST_MORT_RATE_FLG_0] = 1/(1+exp(m_pLRateSmoothParam->dBeta2 * 
			(dNewLowM-m_pLRateSmoothParam->dAlpha2)));

		dParamsArray[GF30Y_LOWEST_MORT_RATE_FLG_25] = (1+m_pLRateSmoothParam->dDensity)/
			(1+exp(m_pLRateSmoothParam->dBeta1 * (dNewLowM-m_pLRateSmoothParam->dAlpha1)));

		if( Min2(mtgLag1, mtgLag2) - m_dNewLow <= 0.005 && Min2(mtgLag1, mtgLag2) - m_dNewLow > 0.0025)
			dParamsArray[GF30Y_LOWEST_MORT_RATE_FLG_50] = 1;
		else
			dParamsArray[GF30Y_LOWEST_MORT_RATE_FLG_50] = 0;


		// Age flag
		if( loan_age <= 3 )
			dParamsArray[GF30Y_LOAN_AGE_FLG0_3] = 1;
		else
			dParamsArray[GF30Y_LOAN_AGE_FLG0_3] = 0;
		if( loan_age <=6  && loan_age > 3)
			dParamsArray[GF30Y_LOAN_AGE_FLG3_6] = 1;
		else
			dParamsArray[GF30Y_LOAN_AGE_FLG3_6] = 0;
		if( loan_age <=12  && loan_age >6 )
			dParamsArray[GF30Y_LOAN_AGE_FLG6_12] = 1;
		else
			dParamsArray[GF30Y_LOAN_AGE_FLG6_12] = 0;

		//update m_dNewLow
		m_dNewLow = Min2(mtgLag1, m_dNewLow);  

		// Here the rate is multiplied by 100 since MIAC interest rates are decimals.

		dParamsArray[GF30Y_MORT_RATE_SPR_SWAP10] = (mtgLag2 - sw10YLag2)*100;
		dParamsArray[GF30Y_CURVE] = (sw10YLag2 - sw2YLag2) * 100;

		dParamsArray[GF30Y_CURVE_LAG2_INCENTIVE] = (sw10YLag2 - sw2YLag2) * 100 *incent;
		dParamsArray[GF30Y_FACTOR_LAG2_INCENTIVE] = 
		  incent * log(factor/(1-Min2(factor, 0.9999999)));        
		
		//Currently FIXED
		dParamsArray[GF30Y_VIX] = GetVIX(iSimMonth);
		
		//Currently FIXED
		dParamsArray[GF30Y_HOUSING_TURN_OVER] = GetHousingTurnover(iSimMonth);

		//Currently FIXED
		dParamsArray[GF30Y_TREND] = (( 2001-1985 )*12+12)/12;

		iCurMonth = (iSimMonth+m_dtValuationMonth.GetMonth()-1) % 12;
		
		if(iCurMonth < 11)
			dParamsArray[GF30Y_MONTH1+iCurMonth] = 1;

		
		// CPR Formula:
		dCpr[nModel] = 0;

	#ifdef DEBUG_DUMP
		double a, b, c;
		char szFileBuf[256];
		sprintf(szFileBuf, "%d,%f,%f,%f,%f,%f,%d,%d,%f,", iSimMonth, 
			chf_dSwap_2y[iSimMonth],
			m_stFwdRates[iSimMonth].dSwap_5y,
			chf_dSwap_10y[iSimMonth],
			chf_dMortgage_rate[iSimMonth],
			m_dNewLow,
			m_nLowCumCount,
			m_nLow2Now,
			m_lpPrepayParams->remain_balance);
		cf.Write(szFileBuf, lstrlen(szFileBuf));
	#endif

	  for(i=GF30Y_CONSTANT; i<GF30Y_MAX_PARAMS; i++)
		{

	#ifdef DEBUG_DUMP
			a = dParamsArray[i];
			b = m_pCoeffs[GF30Y_MAX_PARAMS*nModel+i];
			sprintf(szFileBuf, "%f, %f, ", a, b);
			cf.Write(szFileBuf, lstrlen(szFileBuf));
			c = a * b;
	#endif
			
			dCpr[nModel] += dParamsArray[i] * m_pCoeffs[GF30Y_MAX_PARAMS*nModel+i];
		}
	}
	double dCprFinal=0;

	for(nModel=0; nModel<4; nModel++)
	{
		dCprFinal += dModel[nModel]*dCpr[nModel];
	}

	if(dCprFinal <= -100.0)
		dCprFinal = 0;
	else
		dCprFinal = 1 / (1 + exp(-dCprFinal));

	//m_dNextWAC = Min2(-0.00077812+1.00009*wac*100-0.001128*dCpr, wac*100)/100;

#ifdef DEBUG_DUMP	
    char szFileBuf[256];
	sprintf(szFileBuf, "%f\r\n", dCprFinal);
	cf.Write(szFileBuf, lstrlen(szFileBuf));
	cf.Close();
#endif

  return dCprFinal;
}

/*
Function Name: CalcPrepayFNMA15
Class: CForecast

Purpose: Calculates CPR for FNMA 15 Year Product
*/
double CForecast::CalcPrepayFNMA15(int iSimMonth)
{
	double dCpr[4];
	double wac, incent, incent1, incent2, incent3;
	double factor, loan_age;
	double mtgLag1, mtgLag2, mtgLag3, sw10YLag2, sw2YLag2;
	int iCurMonth, nModel, i;
	double dParamsArray[CF15Y_MAX_PARAMS];
	double calib_incent_min, calib_incent_max, spread_15y;

#ifdef DEBUG_DUMP
	char* pszFileName = "c:\\DUMP.csv";
    char szFileBuf[256];
	CFile cf;
	char   szCause[255];
	CFileException fileException;

	if ( !cf.Open( pszFileName, CFile::modeCreate |    
			  CFile::modeWrite | CFile::modeNoTruncate , &fileException ))
	{
		fileException.GetErrorMessage(szCause, 255);
		
		TRACE( "Can't open file %s, error = %s\n",
		   pszFileName, szCause );
	}	
	cf.SeekToEnd();
#endif

	memset(dParamsArray, 0, CF15Y_MAX_PARAMS*sizeof(double));
	
	spread_15y = m_pConstants[C_SPOT_PMMS_F30Y]-m_pConstants[C_SPOT_PMMS_F15Y];
	
	if(szPrevTranche_ID != m_lpPrepayParams->szTranche_id)
	{
		pTranche = m_pTrancheDataObj->GetTrancheData(m_lpPrepayParams->szTranche_id);
		szPrevTranche_ID = m_lpPrepayParams->szTranche_id;
	}

	if(iSimMonth == 0)
		m_dNextWAC = m_lpPrepayParams->gross_wac;

	wac = m_dNextWAC;

	loan_age = m_lpPrepayParams->wala; //+ iSimMonth;
	
	// Gets the LAG Rates (NOTE: Order of the Constants 
	// is used here to determine the 1m, 2m, ... LAG rates.)
	// HENCE DO NOT CHANGE THE ORDER OF THE GLOBAL CONSTANTS.
	if(iSimMonth-1<0) 
		mtgLag1 = m_pConstants[C_1M_LAG_MTG_RATE];
	else
		mtgLag1 = chf_dMortgage_rate[iSimMonth-1] + m_pConstants[C_SPOT_PMMS_F15Y] - m_pConstants[C_SPOT_MTG_F15Y];

	if(iSimMonth-2<0)
	{
		mtgLag2 = m_pConstants[C_1M_LAG_MTG_RATE+1-iSimMonth];
		sw2YLag2 = m_pConstants[C_1M_LAG_SW2Y+1-iSimMonth];
		sw10YLag2 = m_pConstants[C_1M_LAG_SW10Y+1-iSimMonth];
	}
	else
	{
		mtgLag2 = chf_dMortgage_rate[iSimMonth-2] + m_pConstants[C_SPOT_PMMS_F15Y] - m_pConstants[C_SPOT_MTG_F15Y];
		sw2YLag2 = chf_dSwap_2y[iSimMonth-2];
		sw10YLag2 = chf_dSwap_10y[iSimMonth-2];
	}

	if(iSimMonth-3<0) 
		mtgLag3 = m_pConstants[C_1M_LAG_MTG_RATE+2-iSimMonth];
	else
		mtgLag3 = chf_dMortgage_rate[iSimMonth-3] + m_pConstants[C_SPOT_PMMS_F15Y] - m_pConstants[C_SPOT_MTG_F15Y];

	NewLoanAdjustment(mtgLag1, mtgLag2, mtgLag3, loan_age, chf_dMortgage_rate[iSimMonth] + m_pConstants[C_SPOT_PMMS_F15Y] - m_pConstants[C_SPOT_MTG_F15Y]);

	// Need to transform cf30y rate to cf15y rate?       
	incent1 = sqrt(wac/Max2(mtgLag1-spread_15y, 0.01))-1;
	incent2 = sqrt(wac/Max2(mtgLag2-spread_15y, 0.01))-1;
	incent3 = sqrt(wac/Max2(mtgLag3-spread_15y, 0.01))-1;
	

//Model Pick
/*	nModel = 0;
	if(incent2 >= m_pCoeffs[CF15Y_LAG2_INCENTIVE_MIN] 
	   && loan_age <= m_pCoeffs[CF15Y_PICK_LOAN_AGE_MAX])
		nModel = 0;

	if(incent2 < m_pCoeffs[CF15Y_MAX_PARAMS+CF15Y_LAG2_INCENTIVE_MAX] 
	   && loan_age <= m_pCoeffs[CF15Y_MAX_PARAMS+CF15Y_PICK_LOAN_AGE_MAX])
		nModel = 1;

	//Rag need to check this
	if(incent2 >= m_pCoeffs[CF15Y_MAX_PARAMS*2+CF15Y_LAG2_INCENTIVE_MIN] 
	   && loan_age > m_pCoeffs[CF15Y_MAX_PARAMS*2+CF15Y_PICK_LOAN_AGE_MIN])
		nModel = 2;

	if(incent2 < m_pCoeffs[CF15Y_MAX_PARAMS*2+CF15Y_LAG2_INCENTIVE_MAX] 
	   && loan_age >= m_pCoeffs[CF15Y_MAX_PARAMS*2+CF15Y_PICK_LOAN_AGE_MIN])
		nModel = 3;
*/
	double dModel[4];
	double p,q,r;
	
	if(fabs(incent2-m_pCoeffs[CF15Y_LAG2_INCENTIVE_MIN])<0.01) 
		p = (incent2-m_pCoeffs[CF15Y_LAG2_INCENTIVE_MIN]+0.01)/(0.01 + 0.01);
	else if(incent2 >= (0.01+m_pCoeffs[CF15Y_LAG2_INCENTIVE_MIN]))
			p = 1;
		else
			p = 0;

	if(fabs(incent2-m_pCoeffs[CF15Y_MAX_PARAMS*2+CF15Y_LAG2_INCENTIVE_MIN])<0.01) 
		q = (incent2-m_pCoeffs[CF15Y_MAX_PARAMS*2+CF15Y_LAG2_INCENTIVE_MIN]+0.01)/(0.01 + 0.01);
	else if(incent2 >= (0.01+m_pCoeffs[CF15Y_MAX_PARAMS*2+CF15Y_LAG2_INCENTIVE_MIN]))
			q = 1;
		else
			q = 0;

	if(loan_age > 84 && loan_age < (84+24)) 
		r = (loan_age-84)/24.0;
	else if(loan_age >= (m_pCoeffs[CF15Y_MAX_PARAMS*3+CF15Y_PICK_LOAN_AGE_MIN] + 24))
		r = 1;
	else 
		r = 0;

	dModel[0] = (1-r)*p;
	dModel[1] = (1-r)*(1-p);
	dModel[2] = q*r;
	dModel[3] = (1-q)*r;

	for(nModel=0; nModel<4; nModel++)
	{
		if(dModel[nModel] == 0) 
		{
			dCpr[nModel] = 0;
			continue;
		}


	//Calibration (sets upper and lower limits for incentive)
		calib_incent_min = m_pCoeffs[CF15Y_MAX_PARAMS*nModel+CF15Y_INCENTIVE_MIN];
		calib_incent_max = m_pCoeffs[CF15Y_MAX_PARAMS*nModel+CF15Y_INCENTIVE_MAX];

		if(loan_age <= 6)
			dParamsArray[CF15Y_LOAN_AGE_LAG2_INCENTIVE] = 6*incent2;
		else if(loan_age <= 60)
			dParamsArray[CF15Y_LOAN_AGE_LAG2_INCENTIVE] = loan_age*incent2;
		else
			dParamsArray[CF15Y_LOAN_AGE_LAG2_INCENTIVE] = (60+log(loan_age-59))*incent2;

		incent1 = Limit(incent1, calib_incent_min, calib_incent_max);
		incent2 = Limit(incent2, calib_incent_min, calib_incent_max);
		incent3 = Limit(incent3, calib_incent_min, calib_incent_max);

	//  SET UP PARAMETERS
		dParamsArray[CF15Y_CONSTANT] = 1;
		incent = dParamsArray[CF15Y_LAG2_INCENTIVE] = incent2;
		dParamsArray[CF15Y_LAG1_INCENTIVE] = incent1;
		dParamsArray[CF15Y_LAG3_INCENTIVE] = incent3;

		if(pTranche)
			factor = m_lpPrepayParams->remain_balance/pTranche->dOriginalBalance;
		else
			factor = 0.99; // Default if the Tranche Lookup Fails.
			
		dParamsArray[CF15Y_POOL_FACTOR] = log(factor/(1-Min2(factor, 0.99)));
		dParamsArray[CF15Y_POOL_FACTOR_SQ] = pow((dParamsArray[CF15Y_POOL_FACTOR]),2);
		dParamsArray[CF15Y_LOAN_AGE] = loan_age;
		dParamsArray[CF15Y_LOAN_AGE_SQ] = pow(dParamsArray[CF15Y_LOAN_AGE], 2);
		
		//should  change the above into 
		if(iSimMonth == 0) 
		{
			m_dNewLow = m_pConstants[C_MIN_MTG_5Y];
			m_nLowCumCount = (int) m_pConstants[C_MIN_MTG_LOW_CUM_COUNT];
			m_nLow2Now = (int) m_pConstants[C_MIN_MTG_COUNT];
		}

	/*
		if( Min2(mtgLag1, mtgLag2) - m_dNewLow <= 0.000)
			dParamsArray[CF15Y_LOWEST_MORT_RATE_FLG_0] = 1;
		else
			dParamsArray[CF15Y_LOWEST_MORT_RATE_FLG_0] = 0;
		if( Min2(mtgLag1, mtgLag2) - m_dNewLow <= 0.005 && Min2(mtgLag1, mtgLag2) - m_dNewLow > 0.000)
			dParamsArray[CF15Y_LOWEST_MORT_RATE_FLG_50] = 1;
		else
			dParamsArray[CF15Y_LOWEST_MORT_RATE_FLG_50] = 0;
		
		if ( Min2(mtgLag1, mtgLag2) - m_dNewLow <= 0.000 )
		{
			m_nLowCumCount=1;
			m_nLow2Now=1;
		}
		else 
		{
			m_nLowCumCount=m_nLowCumCount + ((Min2(mtgLag1, mtgLag2) - m_dNewLow <= 0.005) ? 1:0);
		  m_nLow2Now = m_nLow2Now + 1; 
		}
		
		if(m_nLow2Now >4 && m_nLowCumCount >= m_nLow2Now) dParamsArray[CF15Y_LOWEST_MORT_RATE_FLG_50] = 0;
	*/
		double dNewLowM = (Min2(mtgLag1, mtgLag2) - m_dNewLow) * 100;

		dParamsArray[CF15Y_LOWEST_MORT_RATE_FLG_0] = (1.0 + m_pLRateSmoothParam->dDensity) /
			(1+exp(m_pLRateSmoothParam->dBeta2 * (dNewLowM-m_pLRateSmoothParam->dAlpha2)));

		dParamsArray[CF15Y_LOWEST_MORT_RATE_FLG_50] = 1.0 /
			(1+exp(m_pLRateSmoothParam->dBeta1 * (dNewLowM-m_pLRateSmoothParam->dAlpha1)));


		//update m_dNewLow
		m_dNewLow = Min2(mtgLag1, m_dNewLow);  

		// Here the rate is multiplied by 100 since MIAC interest rates are decimals.

		dParamsArray[CF15Y_MORT_RATE_SPR_SWAP10] = (mtgLag2 - sw10YLag2 - spread_15y)*100;
		dParamsArray[CF15Y_CURVE] = (sw10YLag2 - sw2YLag2) * 100;

		
		//Currently FIXED
		dParamsArray[CF15Y_VIX] = GetVIX(iSimMonth);
		
		//Currently FIXED
		dParamsArray[CF15Y_HOUSING_TURN_OVER] = GetHousingTurnover(iSimMonth);


		iCurMonth = (iSimMonth+m_dtValuationMonth.GetMonth()-1) % 12;
		
		if(iCurMonth < 11)
			dParamsArray[CF15Y_MONTH1+iCurMonth] = 1;

		
		// CPR Formula:
		dCpr[nModel] = 0;

	#ifdef DEBUG_DUMP
		double a, b, c;
		sprintf(szFileBuf, "%d,%f,%f,%f,%f,%f,%d,%d,%f,", iSimMonth, 
			chf_dSwap_2y[iSimMonth],
			m_stFwdRates[iSimMonth].dSwap_5y,
			chf_dSwap_10y[iSimMonth],
			chf_dMortgage_rate[iSimMonth],
			m_dNewLow,
			m_nLowCumCount,
			m_nLow2Now,
			m_lpPrepayParams->remain_balance);
		cf.Write(szFileBuf, lstrlen(szFileBuf));
	#endif

	  for(i=CF15Y_CONSTANT; i<CF15Y_MAX_PARAMS; i++)
		{

	#ifdef DEBUG_DUMP
			a = dParamsArray[i];
			b = m_pCoeffs[CF15Y_MAX_PARAMS*nModel+i];
			sprintf(szFileBuf, "%f, %f, ", a, b);
			cf.Write(szFileBuf, lstrlen(szFileBuf));
			c = a * b;
	#endif
			
			dCpr[nModel] += dParamsArray[i] * m_pCoeffs[CF15Y_MAX_PARAMS*nModel+i];
		}
	}
	
	double dCprFinal=0;

	for(nModel=0; nModel<4; nModel++)
	{
		dCprFinal += dModel[nModel]*dCpr[nModel];
	}

	if(dCprFinal <= -100.0)
		dCprFinal = 0;
	else
		dCprFinal = 1 / (1 + exp(-dCprFinal));

	//m_dNextWAC = Min2(-0.00077812+1.00009*wac*100-0.001128*dCpr, wac*100)/100;

#ifdef DEBUG_DUMP	
	
	sprintf(szFileBuf, "%f\r\n",  dCpr);
	cf.Write(szFileBuf, lstrlen(szFileBuf));
	cf.Close();
#endif

  return dCprFinal;
}

double CForecast::GetHousingTurnover(int iSimMonth)
{
	double dValue;

	dValue = 0;

	dValue = m_pModelParameters[MPARAM_HT][iSimMonth];
	return dValue;
}

double CForecast::GetVIX(int iSimMonth)
{
	double dValue;

	dValue = 0;

	dValue = m_pModelParameters[MPARAM_VIX][iSimMonth];
	return dValue;
}

void CForecast::NewLoanAdjustment(double& mtgLag1,double& mtgLag2,double& mtgLag3,double loan_age, double dMortgage)
{
	if(loan_age <= 1)
		mtgLag1 = mtgLag2 = mtgLag3 = dMortgage;
	else if(loan_age <= 2)
		mtgLag2 = mtgLag3 = mtgLag1;
	else if(loan_age <= 3)
		mtgLag3 = mtgLag2;
}		

// Exported Functions
short WINAPI LoadConstants()
{
#ifndef DUMMY_DLL
	ForecastObj.Load();
	return ForecastObj.m_bLoaded;
#else
	MessageBox(NULL, "Called LoadConstants - PrepayDB", "CHF Prepayment Model", MB_OK);
	return TRUE;
#endif
}

void CForecast::ReloadRates()
{
	if(lastSw2Y != m_stFwdRates[12].dSwap_2y || 
		lastSw5Y !=  m_stFwdRates[12].dSwap_5y ||
		lastSw10Y != m_stFwdRates[12].dSwap_10y || 
		lastMtgRate != m_stFwdRates[12].dMortgage_rate ||
		!m_bReloadRates)
	{
		lastSw2Y = m_stFwdRates[12].dSwap_2y;
		lastSw5Y =  m_stFwdRates[12].dSwap_5y;
		lastSw10Y = m_stFwdRates[12].dSwap_10y;
		lastMtgRate = m_stFwdRates[12].dMortgage_rate;
	
		SetSwapOverride();
		SetMortgageOverride();
		//DumpRates();

		m_bReloadRates = 1;
	}
}

void CForecast::SetSwapOverride()
{
	int i;
	
	for(i=0; i<CHF_MAX_FORWARD_RATES; i++)
	{
		chf_dSwap_2y[i] = m_stFwdRates[i].dSwap_2y;
		chf_dSwap_10y[i] = m_stFwdRates[i].dSwap_10y;

		if(m_pModelConstants[B_USE_CURVE_ARRAY] == 1)
		{
			chf_dSwap_2y[i] =  m_pRateOverride[OV_SWAP_2Y][i] + (m_stFwdRates[i].dSwap_10y - m_pRateOverride[OV_SWAP_10Y_NEW][i]);
			chf_dSwap_10y[i] = m_pRateOverride[OV_SWAP_10Y][i] + (m_stFwdRates[i].dSwap_10y - m_pRateOverride[OV_SWAP_10Y_NEW][i]);
		}
	}
}

void CForecast::SetMortgageOverride()
{
	int i;
	int t2, t3;
	double dMtg15YSpread;
	
	double m2, m3, m22, m23;

	dMtg15YSpread = m_pConstants[C_SPOT_MTG_F15Y]-m_pConstants[C_SPOT_MTG_F30Y];

	t2 = (int) m_pModelConstants[REVERSION_START_TIME];
	t3 = (int) m_pModelConstants[REVERSION_END_TIME];
	
	t2 = (t2<0) ? 0:t2;
	t2 = (t2>CHF_MAX_FORWARD_RATES) ? CHF_MAX_FORWARD_RATES:t2;

	t3 = (t3<0) ? 0:t3;
	t3 = (t3>CHF_MAX_FORWARD_RATES) ? CHF_MAX_FORWARD_RATES:t3;

	m2 = m_stFwdRates[t2].dMortgage_rate;
	m3 = chf_dSwap_2y[t3] * m_pModelConstants[FACTOR_30Y_SW2Y] + 
		m_stFwdRates[t3].dSwap_5y * m_pModelConstants[FACTOR_30Y_SW5Y] +
		chf_dSwap_10y[t3] * m_pModelConstants[FACTOR_30Y_SW10Y] + 
		m_pModelConstants[FACTOR_30Y_CONST];
	
	m22 = m_stFwdRates[t2].dMortgage_rate + dMtg15YSpread;
	m23 = chf_dSwap_2y[t3] * m_pModelConstants[FACTOR_15Y_SW2Y] + 
		m_stFwdRates[t3].dSwap_5y * m_pModelConstants[FACTOR_15Y_SW5Y] +
		chf_dSwap_10y[t3] * m_pModelConstants[FACTOR_15Y_SW10Y] + 
		m_pModelConstants[FACTOR_15Y_CONST];

	for(i=0; i<CHF_MAX_FORWARD_RATES; i++)
	{
		chf_dMortgage_rate[i] = m_stFwdRates[i].dMortgage_rate;

		chf_dMortgage_rate_15Y[i] = m_stFwdRates[i].dMortgage_rate + dMtg15YSpread; 

		if(m_pModelConstants[B_USE_MORT_ARRAY] == 1)
		{
			chf_dMortgage_rate[i] = m_pRateOverride[OV_MTG_RATE][i] + (m_stFwdRates[i].dSwap_10y - m_pRateOverride[OV_SWAP_10Y_NEW][i]);
		}

        /*	

  COMMENTED OUT BY BMC; CHF TO MAKE THIS LOGIC OPTIONAL

		if(i>0)
		{
			if(m_pModelConstants[B_USE_MORT_ARRAY] != 1 && m_pModelConstants[FACTOR_30Y_CONST] >= 9999)
			{
				chf_dMortgage_rate[i] = 
					m_stFwdRates[i-1].dMortgage_rate + 
					(chf_dSwap_2y[i]-chf_dSwap_2y[i-1]) * m_pModelConstants[FACTOR_30Y_SW2Y] + 
					(m_stFwdRates[i].dSwap_5y-m_stFwdRates[i-1].dSwap_5y) * m_pModelConstants[FACTOR_30Y_SW5Y] + 
					(chf_dSwap_10y[i]-chf_dSwap_10y[i-1]) * m_pModelConstants[FACTOR_30Y_SW10Y] + 
					0 * m_pModelConstants[FACTOR_30Y_VOL];

				chf_dMortgage_rate_15Y[i] = 
					m_stFwdRates[i-1].dMortgage_rate + dMtg15YSpread + 
					(chf_dSwap_2y[i]-chf_dSwap_2y[i-1]) * m_pModelConstants[FACTOR_15Y_SW2Y] + 
					(m_stFwdRates[i].dSwap_5y-m_stFwdRates[i-1].dSwap_5y) * m_pModelConstants[FACTOR_15Y_SW5Y] + 
					(chf_dSwap_10y[i]-chf_dSwap_10y[i-1]) * m_pModelConstants[FACTOR_30Y_SW10Y] + 
					0 * m_pModelConstants[FACTOR_15Y_VOL];
			}
			
			if(m_pModelConstants[B_USE_MORT_ARRAY] != 1 && m_pModelConstants[FACTOR_30Y_CONST] < 9999)
			{
				
				if(i>t3)
				{
					chf_dMortgage_rate[i] = chf_dSwap_2y[i] * m_pModelConstants[FACTOR_30Y_SW2Y] + 
						m_stFwdRates[i].dSwap_5y * m_pModelConstants[FACTOR_30Y_SW5Y] +
						chf_dSwap_10y[i] * m_pModelConstants[FACTOR_30Y_SW10Y] + 
						m_pModelConstants[FACTOR_30Y_CONST];

					chf_dMortgage_rate_15Y[i] = chf_dSwap_2y[i] * m_pModelConstants[FACTOR_15Y_SW2Y] + 
						m_stFwdRates[i].dSwap_5y * m_pModelConstants[FACTOR_15Y_SW5Y] +
						chf_dSwap_10y[i] * m_pModelConstants[FACTOR_15Y_SW10Y] + 
						m_pModelConstants[FACTOR_15Y_CONST];
				}
				else if(i<t2)
				{
					chf_dMortgage_rate[i] = m_stFwdRates[i].dMortgage_rate;
					chf_dMortgage_rate_15Y[i] = m_stFwdRates[i].dMortgage_rate + dMtg15YSpread;
				}
				else if(t2!=t3)
				{
					chf_dMortgage_rate[i] = (m3-m2)/(t3-t2)*(i-t2) + m2;
					chf_dMortgage_rate_15Y[i] = (m23-m22)/(t3-t2)*(i-t2) + m22;
				}
			}
		}
        */
	}
}

void CForecast::DumpRates()
{
	char* pszFileName = "./Rates.csv";
	CFile cf;
	char szFileBuf[1024];
	int i;

	if(1 || m_pModelConstants[B_CURVE_OUT] == 1)
	{	
		// Dump Interest rates into a file
		CFileException fileException;

		if ( !cf.Open( pszFileName, CFile::modeCreate |    
				  CFile::modeWrite , &fileException ))  // | CFile::modeNoTruncate 
		{
			TRACE( "Can't open file %s, error = %u\n",
			   pszFileName, fileException.m_cause );
		}
		cf.SeekToEnd();
		sprintf(szFileBuf, "simMonth,1m,3m,6m,9m,12m,2y,3y,4y,5y,7y,10y,15y,20y,30y,t3m,t6m,t1y,t2y,t5y,t10y,t30y,mtgRate,mtg15Yr\n");
		cf.Write(szFileBuf, lstrlen(szFileBuf));

		for(i=0; i<CHF_MAX_FORWARD_RATES; i++)
		{

			sprintf(szFileBuf, "%d,%14.12f,%14.12f,%14.12f,%14.12f,%14.12f,%14.12f,%14.12f,%14.12f,%14.12f,%14.12f,%14.12f,%14.12f,%14.12f,%14.12f,%14.12f,%14.12f,%14.12f,%14.12f,%14.12f,%14.12f,%14.12f,%14.12f,%14.12f\n", 
				i, 
				m_stFwdRates[i].dSwap_1m,
				m_stFwdRates[i].dSwap_3m,
				m_stFwdRates[i].dSwap_6m,
				m_stFwdRates[i].dSwap_9m,
				m_stFwdRates[i].dSwap_12m,
				chf_dSwap_2y[i],
				m_stFwdRates[i].dSwap_3y,
				m_stFwdRates[i].dSwap_4y,
				m_stFwdRates[i].dSwap_5y,
				m_stFwdRates[i].dSwap_7y,
				chf_dSwap_10y[i],
				m_stFwdRates[i].dSwap_15y,
				m_stFwdRates[i].dSwap_20y,
				m_stFwdRates[i].dSwap_30y,
				m_stFwdRates[i].dTreasury_3m,
				m_stFwdRates[i].dTreasury_6m,
				m_stFwdRates[i].dTreasury_1y,
				m_stFwdRates[i].dTreasury_2y,
				m_stFwdRates[i].dTreasury_5y,
				m_stFwdRates[i].dTreasury_10y,
				m_stFwdRates[i].dTreasury_30y,
				chf_dMortgage_rate[i],
				chf_dMortgage_rate_15Y[i]);

			cf.Write(szFileBuf, lstrlen(szFileBuf));
		}
		cf.Close();
	}
}


/*
// OLD DLL Type
BOOL WINAPI CHFSetInitParams(CHF_RATE_STRUCTURE *pRateStruct, COleDateTime date)
{
	int yyyy,mm,dd;
	BOOL bLoaded = TRUE;
	
#ifndef DUMMY_DLL
	ForecastObj.Load();
	if(ForecastObj.m_bLoaded)
	{
		yyyy = (int)ForecastObj.m_pConstants[C_DATA_MONTH]/100;
		mm = (int)ForecastObj.m_pConstants[C_DATA_MONTH]-yyyy*100+1;
	
		if(mm==13)
		{
			yyyy+=1;
			mm=1;
		}

		dd=1;
		ForecastObj.m_dtValuationMonth = CTime(yyyy,mm,dd,0,0,0);
	}
	bLoaded = ForecastObj.m_bLoaded;
#endif

	return bLoaded;
}
*/
// New DLL Type
BOOL WINAPI CHFSetInitParams(CHF_RATE_STRUCTURE *pRateStruct, COleDateTime date, LPSTR path, BOOL reload)
{
	int yyyy,mm,dd;
	BOOL bLoaded = TRUE;
	
	if(reload)
	{	
		szGblDBName = path;
#ifndef DUMMY_DLL
		ForecastObj.Unload();
		ForecastObj.m_bReloadRates = 0; // Force reload
#endif
	}
#ifndef DUMMY_DLL
	ForecastObj.Load();
	if(ForecastObj.m_bLoaded)
	{
		yyyy = (int)ForecastObj.m_pConstants[C_DATA_MONTH]/100;
		mm = (int)ForecastObj.m_pConstants[C_DATA_MONTH]-yyyy*100+1;
	
		if(mm==13)
		{
			yyyy+=1;
			mm=1;
		}

		dd=1;
		ForecastObj.m_dtValuationMonth = CTime(yyyy,mm,dd,0,0,0);
	}
	bLoaded = ForecastObj.m_bLoaded;
#else
	MessageBox(NULL, "Called CHFInitParams - PrepayDB", "CHF Prepayment Model", MB_OK);
#endif

	return bLoaded;
}


CHF_RATE_STRUCTURE * WINAPI CHFGetForwardRatesPointer()
{
#ifdef DUMMY_DLL
	MessageBox(NULL, "Called CHFInitParams - PrepayDB", "CHF Prepayment Model", MB_OK);
#endif
	return &ForecastObj.m_stFwdRates[0];
}

double WINAPI CHFCalcPrepay(CHF_PREPAYMENTMODELSTRUCTURE *pPrepayStruct)
{
	double ret = 0.13;

	ForecastObj.m_lpPrepayParams = pPrepayStruct;
#ifndef DUMMY_DLL
	try 
	{
		ret = ForecastObj.CalcPrepay();
	}
	catch(chfDLLVersionError e)
	{
		MessageBox(NULL, e.GetErrorText(), "CHF Prepayment Model", MB_OK);
		pPrepayStruct->error_code = (enum e_CHF_error_code) 1;
	}
#endif

	return ret;
}

void WINAPI CHFSetDefaults()
{
	MessageBox(NULL, "Called Setup in - PrepayDB", "CHF Prepayment Model", MB_OK);
}

ProductTypeArray* WINAPI CHFGetProductList(int &iNumProducts)
{
	int i;
	productTypes.ClearAll();
	iNumProducts = 64;
	for (i=0; i<iNumProducts; i++)
	{
		productTypes.add(i, szGblProductTypes[i]);
	}
#ifdef DUMMY_DLL
	MessageBox(NULL, "Called CHFGetProductList - PrepayDB", "CHF Prepayment Model", MB_OK);
#endif
	return &productTypes;
}

//************ Monthly SMM to annual CPR 
double  cvt_smm_to_cpr ( double mth_smm )
{
	if(mth_smm >= 1.0)
		return 1.0;
	return (1.0 - exp( 12.0 * log(1.0 - mth_smm) )) * 100;
}

//************ Annual CPR to monthly SMM 
double  cvt_cpr_to_smm ( double ann_cpr )
{
	ann_cpr = ann_cpr/100;	 
	if(ann_cpr >= 1.0)
		 return 1.0;
	return 1.0 - exp( log(1.0 - ann_cpr)/12.0 );
}
