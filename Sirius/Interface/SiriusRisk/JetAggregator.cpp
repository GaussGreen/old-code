//	JetAggregator.cpp : Handles aggregation of a results collection
//						to an Jet database.
//
//	Author :			David Cuin
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "SiriusRisk.h"
#include "JetAggregator.h"
#include "daotabledefex.h"


/*static*/ const DataTypeEnum CJetAggregator::s_avVartypeToDataTypeEnum[13] = {
	/*VT_EMPTY		= 0*/	(DataTypeEnum)0,
	/*VT_NULL		= 1*/	(DataTypeEnum)0,
	/*VT_I2			= 2*/	dbInteger,
	/*VT_I4			= 3*/	dbLong,
	/*VT_R4			= 4*/	dbSingle,
	/*VT_R8			= 5*/	dbDouble,
	/*VT_CY			= 6*/	dbCurrency,
	/*VT_DATE		= 7*/	dbDate,
	/*VT_BSTR		= 8*/	dbText,
	/*VT_DISPATCH	= 9*/	(DataTypeEnum)-1,
	/*VT_ERROR		= 10*/	(DataTypeEnum)-1,
	/*VT_BOOL		= 11*/	dbByte,
	/*VT_VARIANT	= 12*/	(DataTypeEnum)-1
};

STDMETHODIMP CJetAggregator::Aggregate(BSTR Name, IResults* pResults)
{
	begin_function
	Aggregate(estring(Name), pResults);
	end_function
}
void CJetAggregator::Aggregate(const std::string& szTableName, CComPtr<IResults> spResults) const
{
	try {
		CDaoDatabase					db;
		CDaoTableDefEx					td;
		long							nFields;
		table_fields					atf[MAX_TABLE_FIELDS];
		
		// Create + Open
		if (!m_szDatabaseName.size()){
			throw "An output database has not been specified";
		} else {			
			CFile			file;
			CFileStatus		fs;						
			if (!CFile::GetStatus(m_szDatabaseName.c_str(), fs)){		
				db.Create(m_szDatabaseName.c_str());
			}
		}		
		if (!db.IsOpen()){			
			db.Open(m_szDatabaseName.c_str());
		}	
								
		// Write the data and close
		GetTableFields(MAX_TABLE_FIELDS, spResults, reinterpret_cast<table_fields*>(atf), &nFields);
		OpenDestinationTable(db, &td, szTableName, atf, nFields, m_bClearTables);
		WriteResults(db, td, spResults);
		db.Close();
	} catch (CDaoException* pe){
		TCHAR szCause[255];
		pe->GetErrorMessage(szCause, 255);
		pe->Delete();
		throw szCause;
	}
}

HRESULT CJetAggregator::FinalConstruct(void)
{
	m_bRecreateTables = false;
	m_bClearTables = false;
	return S_OK;
}

STDMETHODIMP CJetAggregator::get_ClearTables(VARIANT_BOOL* pVal)
{
	*pVal = m_bClearTables ? VARIANT_TRUE : VARIANT_FALSE;
	return S_OK;
}

STDMETHODIMP CJetAggregator::get_DatabaseName(BSTR* pVal)
{
	return m_szDatabaseName.GetBSTR(pVal);
}

STDMETHODIMP CJetAggregator::get_RecreateTables(VARIANT_BOOL* pVal)
{
	*pVal = m_bRecreateTables ? VARIANT_TRUE : VARIANT_FALSE;
	return S_OK;
}

void CJetAggregator::GetTableFields(long nMaxFields, CComPtr<IResults> spResults, table_fields* atf, long* pnFields) const
//	nMaxFields - Maximum number of fields allowed in a table.
{
	CComVariant							v;
	CParameterMap						pm;
	
	if (spResults->GetFields(&v)) propagate_error;
	pm.SetValue(v);

	if (pm.GetRows() > nMaxFields) throw "You are attempting to create a database table with too many fields. You have " + estring(pm.GetRows()) + ", the limit is " + estring(nMaxFields);	
	*pnFields = pm.GetRows();

	for (long nField = 0; nField < *pnFields; nField++){
		std::string		szName;
		VARTYPE			vt;
		long			nSize;
		
		pm.GetValue(nField, 0, &szName);
		pm.GetValue(nField, 1, &vt);
		pm.GetValue(nField, 2, &nSize);
		
		atf[nField].szName = szName;
		atf[nField].nType = VariantTypeToDataTypeEnum(vt);
		atf[nField].lSize = nSize;
		atf[nField].bRequired = true;		// Optional field support goes here
		atf[nField].bPrimaryKey = false;	// Primary key support goes here
	}
}

//	Opens a table in the destination database, creating it if necessary.
//	We also optionally delete any existing records in the tables.
void CJetAggregator::OpenDestinationTable(CDaoDatabase& db, CDaoTableDefEx* ptd, const std::string& szTable, table_fields* ptf, int nFields, bool bClear) const
//	db - JET destination database
//	ptd - table definition to open
//	szTable - name of the table to open
//	ptf - definition of the table (supplied in case we have to create the table)
//	nFields - number of fields in ptd
//	bClear - true if we delete any existing records in the table
{		
	CDaoTableDefInfo					tdi;
	bool								bCreated;						// true if we needed to create the table		
	CDaoIndexFieldInfo					idx_fi[MAX_PRIMARY_KEY_FIELDS];	// array of fields making the primary key
		
	ptd->m_pDatabase = &db;
	bCreated = false;		
	if (m_bRecreateTables){																
		try {
			estring szSQL = "DROP TABLE [" + szTable + "]";
			db.Execute(szSQL.c_str(), dbSQLPassThrough);
		} catch (...){
			// Do nothing
		}
	}			
	
	try {
		db.GetTableDefInfo(szTable.c_str(), tdi, AFX_DAO_PRIMARY_INFO);
	} catch (...){
		// create the table		
		CDaoFieldInfo	fi;
		CStringList		slPrimary;		// list of fields making up the primary key
		
		bCreated = true;
		memset(&fi.m_nType, NULL, (size_t)&fi.m_strForeignName - (size_t)&fi.m_nType);	// initialise all the non-string members
		
		ptd->Create(szTable.c_str());
		
		for (int nElement = 0; nElement < nFields; nElement++){
			fi.m_strName = ptf[nElement].szName.c_str();
			fi.m_nType = ptf[nElement].nType;
			fi.m_lSize = ptf[nElement].lSize;
			fi.m_bRequired = ptf[nElement].bRequired;
			if (ptf[nElement].bPrimaryKey) slPrimary.AddTail(ptf[nElement].szName.c_str());
			ptd->CreateField(fi);
		}

		// set the primary key	
		if (slPrimary.GetCount() > MAX_PRIMARY_KEY_FIELDS){				
			throw "Fatal Error! Too many fields in the primary key. Maximum is " + estring(MAX_PRIMARY_KEY_FIELDS) + ", you have " + estring(slPrimary.GetCount());
		}									
		if (slPrimary.GetCount()){
			CDaoIndexInfo		idx;
			idx.m_strName = "PrimaryKey";
			idx.m_pFieldInfos = new CDaoIndexFieldInfo[slPrimary.GetCount()];
			
			POSITION pos = slPrimary.GetHeadPosition();
			int	nElement = 0;
			while (pos){
				idx_fi[nElement].m_bDescending = false;
				idx_fi[nElement].m_strName = slPrimary.GetNext(pos);
				nElement++;
			}
			idx.m_pFieldInfos = idx_fi;
			idx.m_nFields = slPrimary.GetCount();
			idx.m_bPrimary = TRUE;
			idx.m_bUnique = TRUE;
			idx.m_bClustered = FALSE;
			idx.m_bIgnoreNulls = FALSE;
			idx.m_bRequired = TRUE;
			idx.m_bForeign = FALSE;
			idx.m_lDistinctCount = 1;	// this is returned, not set			
			ptd->CreateIndex(idx);
		}		
		ptd->Append();
	}

	if (!bCreated){	
		// the table already existed - we need to open it			
		ptd->Open(szTable.c_str());
		if (bClear){
			// delete any records in the table							
			estring szSQL = "DELETE FROM [" + szTable + "]";				
			db.Execute(szSQL.c_str(), dbSQLPassThrough);
		}
	}			
}

STDMETHODIMP CJetAggregator::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IJetAggregator };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CJetAggregator::put_ClearTables(VARIANT_BOOL newVal)
{
	m_bClearTables = newVal ? true : false;
	return S_OK;
}

STDMETHODIMP CJetAggregator::put_DatabaseName(BSTR newVal)
{
	m_szDatabaseName = newVal;
	return S_OK;
}

STDMETHODIMP CJetAggregator::put_RecreateTables(VARIANT_BOOL newVal)
{
	m_bRecreateTables = newVal ? true : false;
	return S_OK;
}

DataTypeEnum CJetAggregator::VariantTypeToDataTypeEnum(VARTYPE vt) const
{
	if (vt < 0 || vt > 12) throw "Unknown or unsupported variant type in CResults::VariantTypeToDataTypeEnum";
	DataTypeEnum dt = s_avVartypeToDataTypeEnum[vt];
	if (dt == -1) throw "Unsupported variant type " + estring(vt);
	return dt;
}

void CJetAggregator::WriteResults(CDaoDatabase& db, CDaoTableDefEx& td, const CComPtr<IResults> spResults) const
//	db - JET destination database
//	td - destination table
//	spResults - results set
{	
	CDaoRecordset						rs;								// Destination recordset.
	long								nResults = 0L;
					
	rs.m_pDatabase = &db;
	rs.Open(&td, dbOpenDynaset);
	spResults->get_Count(&nResults);
	for (CComVariant vResult = 1L; vResult.lVal <= nResults; vResult.lVal++){	
		CComPtr<IResult>	spResult;
		
		spResults->get_Item(vResult, &spResult);
		ATLASSERT(spResult);
		
		// Create a new record.
		rs.AddNew();
		
		// Set the recordset fields.
		CComVariant v;
		spResult->get_Value(&v);
		CParameterMap pm;
		pm.SetValue(v);

		for (long nRow = 0; nRow < pm.GetRows(); nRow++){
			std::string szFieldName;
			std::string szFieldValue;			
			pm.GetValue(nRow, 0, &szFieldName);
			pm.GetValue(nRow, 1, &szFieldValue);						// We can use the string overload even for the numeric cases.
			rs.SetFieldValue(szFieldName.c_str(), szFieldValue.c_str());
		}

		// Commit the recordset.
		rs.Update();
	}

	// Close the recordset.
	rs.Close();	
}