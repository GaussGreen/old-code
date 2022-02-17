/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: mktdatamanagerrep.h,v $
 * Revision 1.1  2004/03/01 07:52:17  ebenhamou
 * Initial revision
 *
 */


/*! \file mktdatamanagerrep.h
 *
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */


#ifndef _INGPCALCULATORS_MKTDATAMANAGERREP_H
#define _INGPCALCULATORS_MKTDATAMANAGERREP_H

#include "gpbase/port.h"

#include "gpbase/rootobject.h"

#include <string>
CC_USING_NS(std,string)

#include <set>
CC_USING_NS(std,set)

#include <glob/dates.h>
#include "typedef.h"

/// forward declaration in global namespace
class ARM_ZeroCurve;
class ARM_VolCurve;
class ARM_VolLInterpol;
class ARM_Model;

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
struct ARM_MarketData_ManagerRep;
class ARM_MarketData_ManagerImp;

typedef void (ARM_MarketData_ManagerRep::*MktDataManagerFunctionVoidToVoid)();
typedef void (ARM_MarketData_ManagerRep::*MktDataManagerFunctionDateToVoid)(const ARM_Date& );
typedef void (ARM_MarketData_ManagerRep::*MktDataManagerFunctionMktDataRepToVoid)( const ARM_MarketData_ManagerRep& rhs );
typedef ARM_MarketData_ManagerImp* (ARM_MarketData_ManagerRep::*MktDataManagerFunctionVoidToMktDataImpCst)() const;


/////////////////////////////////////////////////////////////////
/// \class ARM_MarketData_ManagerImp 
/// this class is responsible for hosting all the market data
///
/// the design is to have the object progressively growing when asked
/// market data... asking at demand market data source for market data
/// object if not already defined... the functions Register...blahblah
/// enables to change the market data, hence giving some flexibility
/// the market data manager keeps track of all the keys for the mkt data
/////////////////////////////////////////////////////////////////

struct ARM_MarketData_ManagerRepStrategy
{
	static MktDataManagerFunctionDateToVoid				CreateMktDataInstanceFunc;
	static MktDataManagerFunctionVoidToVoid				CleanMktDataInstanceFunc;
	static MktDataManagerFunctionMktDataRepToVoid		CopyMktDataInstanceFunc;
	static MktDataManagerFunctionVoidToMktDataImpCst	GetMktDataInstanceFunc;
};



struct ARM_MarketData_ManagerRep : public ARM_RootObject
{
private:
	set<string> itsKeys;
	void ThrowErrorIfKeyMissing( const string& key ) const;
	ARM_MarketData_ManagerImp* itsMktDataManager;

	/// redirector mechanism
	friend struct ARM_MarketData_ManagerRepStrategy;

	/// singleton feature
	void CreateMktDataInstance_Singleton( const ARM_Date& asOf);
	void CleanMktDataInstance_Singleton();
	void CopyMktDataInstance_Singleton( const ARM_MarketData_ManagerRep& rhs );
	ARM_MarketData_ManagerImp* GetMktDataInstance_Singleton() const;

	/// clone feature
	void CreateMktDataInstance_Clone( const ARM_Date& asOf );
	void CleanMktDataInstance_Clone();
	void CopyMktDataInstance_Clone( const ARM_MarketData_ManagerRep& rhs );
	ARM_MarketData_ManagerImp* GetMktDataInstance_Clone() const;

public:
	/// constructor,copy constructor,assignment operator and destructor
	ARM_MarketData_ManagerRep( const ARM_Date& asOf);
	ARM_MarketData_ManagerRep( const ARM_MarketData_ManagerRep& rhs );
	ARM_MarketData_ManagerRep& operator=( const ARM_MarketData_ManagerRep& rhs );
	virtual ~ARM_MarketData_ManagerRep();

	/// accessor method
	/// Zero curve part
	string RegisterZCCurve( const string& indexName, const string& ccy, const string& cvName, 
		const ARM_Date& asOf, const string& source, ARM_ZeroCurve* curve );
	ARM_ZeroCurve* GetZCCurve( const string& indexName, const string& ccy, const string& cvName, 
		const ARM_Date& asOf, const string& source ) const;
	ARM_ZeroCurve* GetZCCurveAndClone( const string& indexName, const string& ccy, const string& cvName, 
		const ARM_Date& asOf, const string& source ) const;

	/// vol market part
	string RegisterVolCurve( 	const string& indexName, const string& ccy, const string& cvName, const ARM_Date& asOf,
		const string& volMktType, const string& volType, const string& source, ARM_VolCurve* volCurve );
	ARM_VolCurve* GetVolCurve( const string& indexName, const string& ccy, const string& cvName, const ARM_Date& asOf,
		const string& volMktType, const string& volType, const string& source ) const;
	ARM_VolCurve* GetVolCurveAndClone( const string& indexName, const string& ccy, const string& cvName, const ARM_Date& asOf,
		const string& volMktType, const string& volType, const string& source ) const;

	/// vol mkt model part
	string RegisterVolMktModel( const string& indexName, const string& ccy, const string& cvName, 
		const ARM_Date& asOf, const string& volMktType,	const string& source, ARM_Model* volMktModel );
	ARM_Model* GetVolMktModel( const string& indexName, const string& ccy, const string& cvName, 
		const ARM_Date& asOf, const string& volMktType,	const string& source ) const;
	ARM_Model* GetVolMktModelAndClone( const string& indexName, const string& ccy, const string& cvName, 
		const ARM_Date& asOf, const string& volMktType,	const string& source ) const;

	/// general accessor
	inline ARM_MarketData_ManagerImp* GetMktDataManager() const { return itsMktDataManager;}

	/// function to get and set data
	string RegisterData( const string& objectKey, ARM_Object* object );
	ARM_Object* GetData(const string& objectKey) const;
	ARM_Object* GetDataAndClone( const string& objectKey ) const;
	void ResetAllData();	/// reset data of the central object
	void ResetMyData();		/// reset data of only the current object

	/// management of the asOfDate (accessor get and set)
    const ARM_Date& GetAsOfDate() const;
	void SetAsOfDate( const ARM_Date& newAsOf ) const;
	inline set<string> GetKeys() const { return itsKeys;}

	/// for detail mode
	bool SetDetailMode( bool value ) const;
	bool GetDetailMode() const;

	bool TestIfKeyMissing( const string& key ) const;

	typedef CC_NS(std,set)<string>::const_iterator const_iterator;
	typedef CC_NS(std,set)<string>::iterator iterator;

	/// iterator support
	iterator begin() { return itsKeys.begin(); }
	const_iterator begin() const { return itsKeys.begin(); }
	iterator end() { return itsKeys.end(); }
	const_iterator end() const { return itsKeys.end(); }

    /// standard ARM Object support
	/// because of the clone implementation in terms of the copy constructor
	/// there is no need to define BitwiseCopy and Copy!
	/// Fix fix there only tostring(0 must be redifined but kernel doesn't supported yet
	virtual void View(char* id, FILE* ficOut) const;
	virtual ARM_Object* Clone() const;
	virtual string toString(const string& indent="", const string& nextIndent="") const;
    virtual string ExportShortName() const { return "LMKTM"; }

};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

