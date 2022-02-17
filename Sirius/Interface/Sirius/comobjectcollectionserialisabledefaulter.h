//	ComObjectCollectionSerialisableDefaulter.h
//
//  Author :							David Cuin
//
//	This is used in CComObjectCollection::get_Item (via virtual functions)
//  when an object exports a name property and a load from a database is
//  attempted. This class is used to provide default values for key variables
//  (best understood by looking at the key structure in this class). The class
//  maintains a stack of keys. When an object of this class is constructed, 
//  the key is placed on a stack. The key is removed on destruction. The
//  function CSiriusComModule::_Load instantiates one of these objects which
//  therefore allows key data to propagate recursively.
//
/////////////////////////////////////////////////////////////////////////////

#ifndef _COMOBJECTCOLLECTIONSERIALISABLEDEFAULTER_H
#define _COMOBJECTCOLLECTIONSERIALISABLEDEFAULTER_H

#define pin_date(DataSource, Date)			CComObjectCollectionSerialisableDefaulter cd((DataSource), (Date));
class CComObjectCollectionSerialisableDefaulter
{	
protected:
	struct key
	{		
		DataSourceEnum					m_ds;
		std::string						m_szDataSource;					// if this is a valid DataSourceEnum then m_ds is set to this.
		long							m_nDate;
		bool							m_bDataSourceValid;				// true if m_ds represents m_szDataSource
	};
	static std::stack<key>				s_keys;	

public:
	explicit CComObjectCollectionSerialisableDefaulter(const std::string& szDataSource, long nDate);
	explicit CComObjectCollectionSerialisableDefaulter(DataSourceEnum ds, long nDate);
	virtual ~CComObjectCollectionSerialisableDefaulter(void);
	static DataSourceEnum				GetDataSource(void);
	static const std::string&			GetDataSourceStr(void);
	static long							GetDate(void);
};

#endif
