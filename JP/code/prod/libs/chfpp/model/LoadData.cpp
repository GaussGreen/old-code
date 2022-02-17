// LoadData.cpp : Contains Functions to load data from Prepay.mdb into memory
//

#include "stdafx.h"
#include <math.h>
#include "Prepay.h"
#include "extern.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

CConstants::CConstants()
{
    m_szSql = "select id, value from tblConstants order by id";
	henv = NULL;
    hdbc = NULL;
    hstmt = NULL;
}

CConstants::~CConstants()
{
	
}

short CConstants::UnLoad()
{
    if(m_pConstants)
		delete [] m_pConstants;
	m_pConstants = 0;
	return TRUE;
}

short CConstants::LoadData()
{
	int i;
	SQLINTEGER ind;
	CString m_szTempSql;
	
	if(ConnectDB()) 
	{
		MessageBox(NULL, "Unable to connect to Access ODBC Database - PrepayDB", "ODBC Error", MB_OK);
		return 1;
	}
	m_szTempSql = "select count(*) as NumConstants from tblConstants";
	retcode = SQLBindCol(hstmt, 1, SQL_C_SLONG, &m_iNumConstants, sizeof(m_iNumConstants), &ind);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;

	retcode = SQLExecDirect(hstmt, (SQLCHAR*)(LPCSTR)m_szTempSql, SQL_NTS);

    if (retcode == SQL_ERROR) {
        SQLCHAR state[5];
        SQLINTEGER errN=0;
        SQLCHAR errMsg[256];
        SQLSMALLINT bufSz=256;
        SQLSMALLINT errSz;

        SQLGetDiagRec(SQL_HANDLE_STMT,hstmt,1,state,&errN,errMsg,bufSz,&errSz);

    }
	retcode = SQLFetch(hstmt);

	SQLCloseCursor( hstmt);
    
	m_pConstants = new double[m_iNumConstants];

	retcode = SQLExecDirect(hstmt, (SQLCHAR*)(LPCSTR)m_szSql, SQL_NTS);
	if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;

    retcode = SQLBindCol(hstmt, 1, SQL_C_SLONG, &m_ConstantsRow.iID, sizeof(&m_ConstantsRow.iID), &m_ConstantsRow.ind_iID);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;
    retcode = SQLBindCol(hstmt, 2, SQL_C_DOUBLE, &m_ConstantsRow.dValue, sizeof(&m_ConstantsRow.dValue), &m_ConstantsRow.ind_iValue);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;

	for(i=0; i<m_iNumConstants; i++)
	{
		retcode = SQLFetch(hstmt);
		if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
			goto load_error;
		m_pConstants[i] = m_ConstantsRow.dValue;
	}

    FreeHandles();
    return 0;

load_error:
    FreeHandles();
    return 1;
}

CModelConstants::CModelConstants()
{
    m_szSql = "select id, value from tblModelConstants order by id";
	henv = NULL;
    hdbc = NULL;
    hstmt = NULL;
}

CModelConstants::~CModelConstants()
{
	
}

short CModelConstants::UnLoad()
{
    if(m_pModelConstants)
		delete [] m_pModelConstants;
	m_pModelConstants = 0;
	return TRUE;
}

short CModelConstants::LoadData()
{
	int i;
	SQLINTEGER ind;
	CString m_szTempSql;
	
	if(ConnectDB()) 
	{
		MessageBox(NULL, "Unable to connect to Access ODBC Database - PrepayDB", "ODBC Error", MB_OK);
		return 1;
	}
	m_szTempSql = "select count(*) as NumConstants from tblModelConstants";
	retcode = SQLBindCol(hstmt, 1, SQL_C_SLONG, &m_iNumModelConstants, sizeof(m_iNumModelConstants), &ind);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;

	retcode = SQLExecDirect(hstmt, (SQLCHAR*)(LPCSTR)m_szTempSql, SQL_NTS);
	retcode = SQLFetch(hstmt);

	SQLCloseCursor( hstmt);
    
	m_pModelConstants = new double[m_iNumModelConstants];

	retcode = SQLExecDirect(hstmt, (SQLCHAR*)(LPCSTR)m_szSql, SQL_NTS);
	if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;

    retcode = SQLBindCol(hstmt, 1, SQL_C_SLONG, &m_ModelConstantsRow.iID, sizeof(&m_ModelConstantsRow.iID), &m_ModelConstantsRow.ind_iID);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;
    retcode = SQLBindCol(hstmt, 2, SQL_C_DOUBLE, &m_ModelConstantsRow.dValue, sizeof(&m_ModelConstantsRow.dValue), &m_ModelConstantsRow.ind_iValue);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;

	for(i=0; i<m_iNumModelConstants; i++)
	{
		retcode = SQLFetch(hstmt);
		if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
			goto load_error;
		m_pModelConstants[i] = m_ModelConstantsRow.dValue;
	}

    FreeHandles();
    return 0;

load_error:
    FreeHandles();
    return 1;
}

CModelParameters::CModelParameters()
{
    m_szSql = "select a.* from tblModelParameters as a order by parameter_id, month";
    henv = NULL;
    hdbc = NULL;
    hstmt = NULL;
}

CModelParameters::~CModelParameters()
{

}

short CModelParameters::UnLoad()
{
	if(m_dModelParameters)
	delete[] m_dModelParameters;
	m_dModelParameters = 0;
	return TRUE;
}

short CModelParameters::LoadData()
{
	int i, j;

	if(ConnectDB()) 
	{
		MessageBox(NULL, "Unable to connect to Access ODBC Database - PrepayDB", "ODBC Error", MB_OK);
		return 1;
	}
	
	m_dModelParameters = new PARRAY[NUM_MODEL_PARAMS];

	retcode = SQLExecDirect(hstmt, (SQLCHAR*)(LPCSTR)m_szSql, SQL_NTS);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;
	
    retcode = SQLBindCol(hstmt, 1, SQL_C_SLONG, &ModelParametersRow.iparameter_id, sizeof(&ModelParametersRow.iparameter_id), &ModelParametersRow.ind_iparameter_id);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;
    retcode = SQLBindCol(hstmt, 2, SQL_C_SLONG, &ModelParametersRow.imonth, sizeof(&ModelParametersRow.imonth), &ModelParametersRow.ind_imonth);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;
    retcode = SQLBindCol(hstmt, 3, SQL_C_CHAR, &ModelParametersRow.szparameter, 11, &ModelParametersRow.ind_iparameter);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;
    retcode = SQLBindCol(hstmt, 4, SQL_C_DOUBLE, &ModelParametersRow.dvalue, sizeof(&ModelParametersRow.dvalue), &ModelParametersRow.ind_ivalue);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;

	for(i=0; i<NUM_MODEL_PARAMS; i++)
	{
		for(j=0; j<CHF_MAX_FORWARD_RATES; j++)
		{
			retcode = SQLFetch(hstmt);
			if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
				goto load_error;
			m_dModelParameters[i][j] = ModelParametersRow.dvalue;
		}
	}


    FreeHandles();
    return 0;

load_error:
    FreeHandles();
    return 1;
}

CModelCoeffs::CModelCoeffs()
{
//    m_szSql = "select a.id, a.product_id, a.model_id, a.parameter_id, a.value from tblModelCoefficients as a";
	m_szSql = "";
    henv = NULL;
    hdbc = NULL;
    hstmt = NULL;
}

CModelCoeffs::~CModelCoeffs()
{

}

short CModelCoeffs::UnLoad()
{
	if(m_pCoeffs)
		delete [] m_pCoeffs;
	m_pCoeffs = 0;
	return TRUE;
}

short CModelCoeffs::LoadData()
{
	SQLINTEGER ind1, ind2;
	CString m_szTempSql;
	int i, iProdCtr;

    if(ConnectDB()) return 1;

	m_szTempSql = "select a.NumProducts, a.NumCoeffs as num_c from qryNumProducts as a";
	retcode = SQLBindCol(hstmt, 1, SQL_C_SLONG, &m_iNumProducts, sizeof(m_iNumProducts), &ind1);
	retcode = SQLBindCol(hstmt, 2, SQL_C_SLONG, &m_iNumCoeffs, sizeof(m_iNumCoeffs), &ind2);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;

	retcode = SQLExecDirect(hstmt, (SQLCHAR*)(LPCSTR)m_szTempSql, SQL_NTS);
	retcode = SQLFetch(hstmt);
	
	SQLCloseCursor( hstmt);

	m_pCoeffs = new double[m_iNumCoeffs];

	m_szTempSql = "select a.id, a.product_id, a.model_id, a.parameter_id, a.value from tblModelCoefficients as a order by a.id";
	
	retcode = SQLExecDirect(hstmt, (SQLCHAR*)(LPCSTR)m_szTempSql, SQL_NTS);
	if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;

    retcode = SQLBindCol(hstmt, 1, SQL_C_SLONG, &ModelCoeffRow.iid, sizeof(&ModelCoeffRow.iid), &ModelCoeffRow.ind_iid);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;
    retcode = SQLBindCol(hstmt, 2, SQL_C_SLONG, &ModelCoeffRow.iproduct_id, sizeof(&ModelCoeffRow.iproduct_id), &ModelCoeffRow.ind_iproduct_id);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;
    retcode = SQLBindCol(hstmt, 3, SQL_C_SLONG, &ModelCoeffRow.imodel_id, sizeof(&ModelCoeffRow.imodel_id), &ModelCoeffRow.ind_imodel_id);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;
    retcode = SQLBindCol(hstmt, 4, SQL_C_SLONG, &ModelCoeffRow.iparameter_id, sizeof(&ModelCoeffRow.iparameter_id), &ModelCoeffRow.ind_iparameter_id);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;
    retcode = SQLBindCol(hstmt, 5, SQL_C_DOUBLE, &ModelCoeffRow.dvalue, sizeof(&ModelCoeffRow.dvalue), &ModelCoeffRow.ind_ivalue);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;
	
	iProdCtr = 0;

	for(i=0; i<m_iNumCoeffs; i++)
	{
		retcode = SQLFetch(hstmt);
		if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
			goto load_error;
		m_pCoeffs[i]=ModelCoeffRow.dvalue;
		if(ModelCoeffRow.iparameter_id == 0)
		{
			iMapProducts[ModelCoeffRow.iproduct_id] = iProdCtr;
			pBaseCoeffs[iProdCtr] = &(m_pCoeffs[i]);
			iProdCtr++;
		}
	}


    FreeHandles();
    return 0;

load_error:
    FreeHandles();
    return 1;
}

CTrancheData::CTrancheData()
{
//    m_szSql = "select a.id, a.OriginalBalance from tblTrancheData as a";
	henv = NULL;
    hdbc = NULL;
    hstmt = NULL;
}

CTrancheData::~CTrancheData()
{
}

short CTrancheData::UnLoad() 
{
	POSITION pos;
	CString key;
	CTranche* pTranche;

	// Iterate through the entire map to delete the contents
	for( pos = m_TrancheMap.GetStartPosition(); pos != NULL; )
	{
	   m_TrancheMap.GetNextAssoc( pos, key, (CObject*&)pTranche );
	   if(pTranche)
		delete pTranche;
	   pTranche = 0;
	}

	m_TrancheMap.RemoveAll();

	return TRUE;
}

short CTrancheData::LoadData()
{
	SQLINTEGER ind1;
	CString m_szTempSql;
	int i;

    if(ConnectDB()) return 1;

	m_szTempSql = "select count(*) as num_t from tblTrancheData  as a";
	retcode = SQLBindCol(hstmt, 1, SQL_C_SLONG, &m_iNumTranches, sizeof(m_iNumTranches), &ind1);

	retcode = SQLExecDirect(hstmt, (SQLCHAR*)(LPCSTR)m_szTempSql, SQL_NTS);
	retcode = SQLFetch(hstmt);
	
	SQLCloseCursor( hstmt);

	m_szTempSql = "select a.id, a.OriginalBalance from tblTrancheData as a order by a.id";

	retcode = SQLBindCol(hstmt, 1, SQL_C_CHAR, &TrancheDataRow.szid, 51, &TrancheDataRow.ind_iid);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;
    retcode = SQLBindCol(hstmt, 2, SQL_C_DOUBLE, &TrancheDataRow.dOriginalBalance, sizeof(&TrancheDataRow.dOriginalBalance), &TrancheDataRow.ind_iOriginalBalance);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;

	retcode = SQLExecDirect(hstmt, (SQLCHAR*)(LPCSTR)m_szTempSql, SQL_NTS);
	if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;

	for(i=0; i<m_iNumTranches; i++)
	{
		retcode = SQLFetch(hstmt);
		if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
			goto load_error;
		m_TrancheMap[(LPCSTR)TrancheDataRow.szid]= new CTranche(TrancheDataRow.dOriginalBalance);
	}

    FreeHandles();
    return 0;

load_error:
    FreeHandles();
    return 1;
}

CTranche *CTrancheData::GetTrancheData(CString szTrancheID)
{
	CTranche *pTranche;
	if(!m_TrancheMap.Lookup( szTrancheID, ( CObject*& ) pTranche ))
		pTranche = NULL;
	
	return pTranche;
}

CPPTables::CPPTables()
{
//    m_szSql = "select * from tblPrepayTables";
	m_szSql = "";
    m_iNumRows = 0;
    henv = NULL;
    hdbc = NULL;
    hstmt = NULL;
}

CPPTables::~CPPTables()
{
    
}

short CPPTables::UnLoad()
{
	if(PPTablesArray)
		delete [] PPTablesArray;
	PPTablesArray = 0;
	return TRUE;
}

short CPPTables::LoadData()
{
	SQLINTEGER ind;
	CString m_szTempSql;
	int i,j,k;

    if(ConnectDB()) return 1;
	
	m_szTempSql = "select num_tables from qryNumPPTables";
	retcode = SQLBindCol(hstmt, 1, SQL_C_SLONG, &m_iNumTables, sizeof(m_iNumTables), &ind);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;

	retcode = SQLExecDirect(hstmt, (SQLCHAR*)(LPCSTR)m_szTempSql, SQL_NTS);
	retcode = SQLFetch(hstmt);
	
	SQLCloseCursor( hstmt);
	PPTablesArray = new float[m_iNumTables][PP_MAXINCENTIVE+1][PP_MAXTERM+1];
	
	m_szTempSql = "SELECT p.ID, p.WAC_diff, p.Months, p.WAC_val FROM tblPrepayTables AS p order by p.ID, p.WAC_diff, p.Months";
	
	retcode = SQLExecDirect(hstmt, (SQLCHAR*)(LPCSTR)m_szTempSql, SQL_NTS);
	if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;

	retcode = SQLBindCol(hstmt, 1, SQL_C_SSHORT, &PPTableRow.sID, sizeof(&PPTableRow.sID), &PPTableRow.ind_iID);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;
    retcode = SQLBindCol(hstmt, 2, SQL_C_DOUBLE, &PPTableRow.dIncentive, sizeof(&PPTableRow.dIncentive), &PPTableRow.ind_iIncentive);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;
    retcode = SQLBindCol(hstmt, 3, SQL_C_SLONG, &PPTableRow.iTerm, sizeof(&PPTableRow.iTerm), &PPTableRow.ind_iTerm);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;
    retcode = SQLBindCol(hstmt, 4, SQL_C_DOUBLE, &PPTableRow.dCPR, sizeof(&PPTableRow.dCPR), &PPTableRow.ind_iCPR);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;
	
	for (i=0;i<m_iNumTables;i++)
	{
		for(j=1;j<=PP_MAXINCENTIVE;j++)
		{
			for(k=1;k<=PP_MAXTERM;k++)
			{
				retcode = SQLFetch(hstmt);
				if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
					goto load_error;

				iMapPPTable[PPTableRow.sID] = i;
				PPTablesArray[i][0][k]=(float)PPTableRow.iTerm;
				PPTablesArray[i][j][0]=(float)PPTableRow.dIncentive;
				PPTablesArray[i][j][k]=(float)PPTableRow.dCPR;
			}
		}
	}

    FreeHandles();
    return 0;

load_error:
    FreeHandles();
    return 1;
}

/*
double CPPTables::GetCPR(int ePPTable, float dIncentive, float dTerm)
{
	int idxTable, i, j;
	float dMin, dMax, dCpr;

	idxTable = iMapPPTable[ePPTable];

	for (i=1; i<PP_MAXINCENTIVE;i++)
	{
		if(PPTablesArray[idxTable][i][0]>=dIncentive) break;
	}

	if(i==PP_MAXINCENTIVE) dIncentive=PPTablesArray[idxTable][i][0];

	i--;

	if(i==0) 
	{
		i=1;
		dIncentive=PPTablesArray[idxTable][i][0];
	}


	for (j=1; j<PP_MAXTERM;j++)
	{
		if(PPTablesArray[idxTable][0][j]>=dTerm) break;
	}
	
	if(j==PP_MAXTERM) dTerm=PPTablesArray[idxTable][0][j];
	
	j--;
	if(j==0) 
	{
		j=1;
		dTerm=PPTablesArray[idxTable][0][j];
	}

	

	dMin = (PPTablesArray[idxTable][i][j+1]-PPTablesArray[idxTable][i][j]) / 
		(PPTablesArray[idxTable][0][j+1]-PPTablesArray[idxTable][0][j]) * 
		(dTerm-PPTablesArray[idxTable][0][j]) + PPTablesArray[idxTable][i][j];

	dMax = (PPTablesArray[idxTable][i+1][j+1]-PPTablesArray[idxTable][i+1][j]) / 
		(PPTablesArray[idxTable][0][j+1]-PPTablesArray[idxTable][0][j]) * 
		(dTerm-PPTablesArray[idxTable][0][j]) + PPTablesArray[idxTable][i+1][j];
	
	dCpr = (dMax-dMin) / 
		(PPTablesArray[idxTable][i+1][0]-PPTablesArray[idxTable][i][0]) * 
		(dIncentive-PPTablesArray[idxTable][i][0]) + dMin;
	
	if(dCpr >= (float)PP_MAX_CPR) dCpr = (float)PP_MAX_CPR;

	return dCpr;
}
*/

double CPPTables::GetCPR(int ePPTable, float dIncentive, float dAge)
{
	int row1, row2;
	int col1, col2;
	int idxTable;

	idxTable = iMapPPTable[ePPTable];

	if (dIncentive <= RowIndexLookup(idxTable, 0)) 
	{
		row1 = 0;
		row2 = 1;
	}
	else if (dIncentive > RowIndexLookup(idxTable, RowIndexCount()-1)) 
	{
	   row1 = RowIndexCount()-2;
	   row2 = RowIndexCount()-1;
	}
	else  
	{
	   for (int i = 0; dIncentive > RowIndexLookup(idxTable, i); ++i);
	   row1 = i-1;
	   row2 = i;
	}

	if (dAge <= ColIndexLookup(idxTable, 0))
	{
	   col1 = 0;
	   col2 = 1;
	}
	else if (dAge > ColIndexLookup(idxTable, ColIndexCount()-1))
	{
	   col1 = ColIndexCount()-2;
	   col2 = ColIndexCount()-1;
	}
	else  
	{
	   for (int i = 0; dAge > ColIndexLookup(idxTable, i); ++i);
	   col1 = i-1;
	   col2 = i;
	}

	double X0, X1, X2;
	double Y0, Y1, Y2;
	double Z1, Z2, Z3, Z4;
	X0 = dIncentive;
	Y0 = dAge;

	X1 = RowIndexLookup(idxTable, row1);
	Y1 = ColIndexLookup(idxTable, col1);

	X2 = RowIndexLookup(idxTable, row2);
	Y2 = ColIndexLookup(idxTable, col2);

	Z1 = cvt_cpr_to_smm(CPRLookup(idxTable, row1, col1));
	Z2 = cvt_cpr_to_smm(CPRLookup(idxTable, row1, col2));
	Z3 = cvt_cpr_to_smm(CPRLookup(idxTable, row2, col1));
	Z4 = cvt_cpr_to_smm(CPRLookup(idxTable, row2, col2));
	double Z0 = ((Z1 - Z2 - Z3 + Z4)*(Y0 - Y2)/(Y1 - Y2) + Z2 - Z4)*(X0 - X2)/(X1 - X2) + (Z3 - Z4)*(Y0 - Y2)/(Y1 - Y2) + Z4;
	return cvt_smm_to_cpr(Z0);
}


CRateOverride::CRateOverride()
{
    m_szSql = "select a.* from tblRateOverride as a order by parameter_id, month";
    henv = NULL;
    hdbc = NULL;
    hstmt = NULL;
	m_iNumRates = 0;
	m_dRateOverride = 0;
}

CRateOverride::~CRateOverride()
{

}

short CRateOverride::UnLoad()
{
	if(m_dRateOverride)
		delete[] m_dRateOverride;
	m_dRateOverride = 0;
	return TRUE;
}

short CRateOverride::LoadData()
{
	int i, j;
	CString m_szTempSql;
	SQLINTEGER ind;

	m_dRateOverride = 0;
	if(ConnectDB()) 
	{
		MessageBox(NULL, "Unable to connect to Access ODBC Database - PrepayDB", "ODBC Error", MB_OK);
		return 1;
	}
	
	m_szTempSql = "select ids from qryGetOverrideCount";
	retcode = SQLBindCol(hstmt, 1, SQL_C_SLONG, &m_iNumRates, sizeof(m_iNumRates), &ind);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;

	retcode = SQLExecDirect(hstmt, (SQLCHAR*)(LPCSTR)m_szTempSql, SQL_NTS);
	retcode = SQLFetch(hstmt);
	
	SQLCloseCursor( hstmt);

	if(m_iNumRates == 0)
		return 0;

	m_dRateOverride = new PARRAY[m_iNumRates];

	retcode = SQLExecDirect(hstmt, (SQLCHAR*)(LPCSTR)m_szSql, SQL_NTS);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;
	
    retcode = SQLBindCol(hstmt, 1, SQL_C_SLONG, &RateOverrideRow.iparameter_id, sizeof(&RateOverrideRow.iparameter_id), &RateOverrideRow.ind_iparameter_id);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;
    retcode = SQLBindCol(hstmt, 2, SQL_C_SLONG, &RateOverrideRow.imonth, sizeof(&RateOverrideRow.imonth), &RateOverrideRow.ind_imonth);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;
    retcode = SQLBindCol(hstmt, 3, SQL_C_CHAR, &RateOverrideRow.szparameter, 11, &RateOverrideRow.ind_iparameter);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;
    retcode = SQLBindCol(hstmt, 4, SQL_C_DOUBLE, &RateOverrideRow.dvalue, sizeof(&RateOverrideRow.dvalue), &RateOverrideRow.ind_ivalue);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;

	for(i=0; i<m_iNumRates; i++)
	{
		for(j=0; j<CHF_MAX_FORWARD_RATES; j++)
		{
			retcode = SQLFetch(hstmt);
			if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
				goto load_error;
			m_dRateOverride[i][j] = RateOverrideRow.dvalue;
		}
	}

load_error:
    FreeHandles();
    return 0;
}


CLowRateSmoothingParam::CLowRateSmoothingParam()
{
    m_szSql = "select a.PRODUCT_ID, a.DENSITY, a.BETA1, a.ALPHA1, a.BETA2, a.ALPHA2 from tblLowRateSmoothingParam AS a order by a.product_id";
	henv = NULL;
    hdbc = NULL;
    hstmt = NULL;
}

CLowRateSmoothingParam::~CLowRateSmoothingParam()
{
	
}

short CLowRateSmoothingParam::UnLoad()
{
	return TRUE;
}

short CLowRateSmoothingParam::LoadData()
{
	int i;
    SQLINTEGER  ind_iProduct_id;
    SQLINTEGER  ind_iDensity;
    SQLINTEGER  ind_iBeta1;
    SQLINTEGER  ind_iAlpha1;
    SQLINTEGER  ind_iBeta2;
    SQLINTEGER  ind_iAlpha2;

	if(ConnectDB())
	{
		MessageBox(NULL, "Unable to connect to Access ODBC Database - PrepayDB", "ODBC Error", MB_OK);
		return 1;
	}
   
	retcode = SQLExecDirect(hstmt, (SQLCHAR*)(LPCSTR)m_szSql, SQL_NTS);
	if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;

    retcode = SQLBindCol(hstmt, 1, SQL_C_SLONG, &LRateSmoothParamRow.iProduct_id, sizeof(&LRateSmoothParamRow.iProduct_id), &ind_iProduct_id);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;
    retcode = SQLBindCol(hstmt, 2, SQL_C_DOUBLE, &LRateSmoothParamRow.dDensity, sizeof(&LRateSmoothParamRow.dDensity), &ind_iDensity);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;
    retcode = SQLBindCol(hstmt, 3, SQL_C_DOUBLE, &LRateSmoothParamRow.dBeta1, sizeof(&LRateSmoothParamRow.dBeta1), &ind_iBeta1);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;
    retcode = SQLBindCol(hstmt, 4, SQL_C_DOUBLE, &LRateSmoothParamRow.dAlpha1, sizeof(&LRateSmoothParamRow.dAlpha1), &ind_iAlpha1);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;
    retcode = SQLBindCol(hstmt, 5, SQL_C_DOUBLE, &LRateSmoothParamRow.dBeta2, sizeof(&LRateSmoothParamRow.dBeta2), &ind_iBeta2);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;
    retcode = SQLBindCol(hstmt, 6, SQL_C_DOUBLE, &LRateSmoothParamRow.dAlpha2, sizeof(&LRateSmoothParamRow.dAlpha2), &ind_iAlpha2);
    if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
        goto load_error;

	for(i=0; i<MAX_15FACTOR_PRODUCTS; i++)
	{
		retcode = SQLFetch(hstmt);
		if (retcode != SQL_SUCCESS && retcode != SQL_SUCCESS_WITH_INFO) 
			goto load_error;
		m_LRateSmoothParam[i].iProduct_id = LRateSmoothParamRow.iProduct_id;
		m_LRateSmoothParam[i].dDensity = LRateSmoothParamRow.dDensity;
		m_LRateSmoothParam[i].dAlpha1 = LRateSmoothParamRow.dAlpha1;
		m_LRateSmoothParam[i].dBeta1 = LRateSmoothParamRow.dBeta1;
		m_LRateSmoothParam[i].dAlpha2 = LRateSmoothParamRow.dAlpha2;
		m_LRateSmoothParam[i].dBeta2 = LRateSmoothParamRow.dBeta2;
	}

    FreeHandles();
    return 0;

load_error:
    FreeHandles();
    return 1;
}

LRATESMOOTHPARAM * CLowRateSmoothingParam::GetSmoothParam(int nProductID)
{
	int idx=0;

	switch(nProductID)
	{
		
		case 0: // FNMA 30 Yr Product
			idx = 0;
			break;
		case 2:	// FNMA 15 Yr Product
			idx = 1;
			break;
		case 10:// GNMA 30 Yr Product
			idx = 2;
			break;
		default:
			break;
	}
	return &m_LRateSmoothParam[idx];
}
