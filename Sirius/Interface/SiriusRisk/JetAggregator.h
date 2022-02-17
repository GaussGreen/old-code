//	AggregatorJet.h : Handles aggregation of a results collection
//					  to an Jet (Joint Engine Technology) database.
//
//	Author :		  David Cuin
//////////////////////////////////////////////////////////////////////

#ifndef __JETAGGREGATOR_H_
#define __JETAGGREGATOR_H_

#include "resource.h"
class CDaoTableDefEx;

class ATL_NO_VTABLE CJetAggregator : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CJetAggregator, &CLSID_JetAggregator>,
	public ISupportErrorInfo,
	public IDispatchImpl<IJetAggregator, &IID_IJetAggregator, &LIBID_SiriusRisk>,
	public IDispatchImpl<IAggregatable, &IID_IAggregatable, &LIBID_SiriusRisk>
{
public:
	CJetAggregator(){}

	DECLARE_REGISTRY_RESOURCEID(IDR_JETAGGREGATOR)
	DECLARE_NOT_AGGREGATABLE(CJetAggregator)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CJetAggregator)
		COM_INTERFACE_ENTRY(IJetAggregator)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
		COM_INTERFACE_ENTRY2(IDispatch, IJetAggregator)
		COM_INTERFACE_ENTRY(IAggregatable)
	END_COM_MAP()

	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_DatabaseName)(/*[out, retval]*/ BSTR* pVal);
	STDMETHOD(put_DatabaseName)(/*[in]*/ BSTR newVal);	
	STDMETHOD(get_RecreateTables)(/*[out, retval]*/ VARIANT_BOOL* pVal);
	STDMETHOD(put_RecreateTables)(/*[in]*/ VARIANT_BOOL newVal);
	STDMETHOD(get_ClearTables)(/*[out, retval]*/ VARIANT_BOOL* pVal);
	STDMETHOD(put_ClearTables)(/*[in]*/ VARIANT_BOOL newVal);
	STDMETHOD(Aggregate)(/*[in]*/ BSTR Name, /*[in]*/ IResults* pResults);

	HRESULT FinalConstruct();

protected:
	#define MAX_PRIMARY_KEY_FIELDS		32								// Maximum number of fields in a destination table primary key.
	#define MAX_TABLE_FIELDS			256								// Maximum number of fields in a destination table. This number is arbitrary.
	
	typedef struct _table_fields {
		std::string						szName;
		short							nType;
		long							lSize;
		bool							bRequired;
		bool							bPrimaryKey;					// True if the field is part of the primary key.
	} table_fields;
		
	estring								m_szDatabaseName;
	bool								m_bRecreateTables;				// True if we delete a database table before writing records.
	bool								m_bClearTables;
	static const DataTypeEnum			s_avVartypeToDataTypeEnum[13];
	
	void								Aggregate(const std::string& szTableName, CComPtr<IResults> spResults) const;
	void								GetTableFields(long nMaxFields, CComPtr<IResults> spResults, table_fields* atf, long* pnFields) const;
	void								OpenDestinationTable(CDaoDatabase& db, CDaoTableDefEx* ptd, const std::string& szTable, table_fields* ptf, int nFields, bool bClear) const;
	DataTypeEnum						VariantTypeToDataTypeEnum(VARTYPE vt) const;
	void								WriteResults(CDaoDatabase& db, CDaoTableDefEx& td, const CComPtr<IResults> spResults) const;	
};

#endif
