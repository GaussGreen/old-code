
#ifndef _INGPCALCULATORS_MKTDATAMANAGER_H
#define _INGPCALCULATORS_MKTDATAMANAGER_H

#include "gpbase/port.h"


#include <string>
CC_USING_NS(std,string)

#include <map>
CC_USING_NS(std,map)
CC_USING_NS(std,less)

#include "typedef.h"
#include <glob/dates.h>
#include "gpbase/rootobject.h"

/// forward declaration in global namespace
class ARM_ZeroCurve;
class ARM_VolCurve;
class ARM_Model;

CC_BEGIN_NAMESPACE( ARM )


/// forward declaration
template <typename T> class ARM_SingletonHolder;
class ARM_MarketData_ManagerImp;
struct ARM_MarketData_ManagerRep;

extern ARM_SingletonHolder<ARM_MarketData_ManagerImp> ARM_TheMarketData_Manager;


/////////////////////////////////////////////////////////////////
/// \class ARM_MarketData_ManagerImp 
/// this class is responsible for hosting all the market data
///
/// the design is to have the object progressively growing when asked
/// market data... asking at demand market data source for market data
/// object if not already defined... the functions Register...blahblah
/// enables to change the market data, hence giving some flexibility
/////////////////////////////////////////////////////////////////

class ARM_MarketData_ManagerImp : public ARM_RootObject
{
private:
	typedef map<string,ARM_ObjectPtr,less<string> > stringARM_ObjectPtrMap;

	// as of represent the as of for all the market data
		ARM_Date itsAsOf;

		/// this content is a full dictionary of market data ... it is filled on demand
		/// if an object is not present when asked ... it is created with default..
		/// and stored as an ARM_ObjectPtr in itsContents
		/// the ARM_ObjectPtr allows us to be sure of deletion of objects!
		CC_IS_MUTABLE stringARM_ObjectPtrMap itsContents;

		/// for view with details
		bool itsDetailViewMode;

	/// function to state wheter a user is authorized to change mkt data
	static bool isAuthorizedToChangedMarketData;

	/// function to test if someone is authorized to change market data!
	static bool GetIsAuthorizedTochangeMarketData();

	/// for all object, throw an excption that an object could not be found
	void ThrowError_ObjectNoFound( const string& key, const string& type ) const;

	/// function to allow pull technology ... (with the help of the market data retriever!)
	string RegisterDataNoClone(	const string& objectKey, ARM_ObjectPtr object ) const;

	void ValidateAsOf( const ARM_Date& asOf ) const;

	static bool UseMktDataRetriever;

	/// copy and assignment operator to avoid duplication
	ARM_MarketData_ManagerImp( const ARM_MarketData_ManagerImp& rhs );
	ARM_MarketData_ManagerImp& operator=(const ARM_MarketData_ManagerImp& rhs);

	/// constructor and destructor private to allow only the singleton to access it!
	void CopyClonedData( const ARM_MarketData_ManagerImp& rhs );
	void CleanData();
	ARM_MarketData_ManagerImp( const ARM_Date& asOf = ARM_Date() );
	~ARM_MarketData_ManagerImp();


	friend class ARM_SingletonHolder<ARM_MarketData_ManagerImp>;
	friend struct ARM_MarketData_ManagerRep;

public:

// FIXMEFRED: mig.vc8 (25/05/2007 11:33:24):missing return type
	inline unsigned int size() const { return static_cast<unsigned int>(itsContents.size());}
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
		const ARM_Date& asOf, const string& volMktType,	const string& source, ARM_Model* volMktModel	);
	ARM_Model* GetVolMktModel( const string& indexName, const string& ccy, const string& cvName, 
		const ARM_Date& asOf, const string& volMktType,	const string& source ) const;
	ARM_Model* GetVolMktModelAndClone( const string& indexName, const string& ccy, const string& cvName, 
		const ARM_Date& asOf, const string& volMktType,	const string& source ) const;

	/// general accessor
	/// function to get and set data
	string RegisterData( const string& objectKey, ARM_Object* object );
	ARM_Object* GetData(const string& objectKey) const;
	ARM_Object* GetDataAndClone( const string& objectKey ) const;
	void ResetAllData();

	/// static method for the computation of the unique key!
	/// Zero curve part
	static string CoinZCObjKey( const string& indexName, const string& ccy, const string& cvName, 
		const ARM_Date& asOf, const string& source );
	/// vol market part
	static string CoinVolObjKey( const string& indexName, const string& ccy, 
		const string& cvName, const ARM_Date& asOf, const string& volMktType, const string& volType,
		const string& source );
	/// vol mkt model part
	static string CoinVolMktModelObjKey( const string& indexName, const string& ccy, const string& cvName, 
		const ARM_Date& asOf, const string& volMktType,	const string& source );

    /// accessor to the asOfDate
    const ARM_Date& GetAsOfDate() const {return itsAsOf;}
	void SetAsOfDate( const ARM_Date& newAsOf );

	/// for detail mode
	inline bool SetDetailMode( bool value ) { return itsDetailViewMode=value;}
	inline bool GetDetailMode() const {return itsDetailViewMode;}

	/// STL Like iterator support
	typedef stringARM_ObjectPtrMap::iterator iterator;
	iterator begin() { return itsContents.begin(); }
	iterator end() { return itsContents.end(); }

	typedef stringARM_ObjectPtrMap::const_iterator const_iterator;
	const_iterator begin() const { return itsContents.begin(); }
	const_iterator end() const { return itsContents.end(); }

	/// only the View to allow representant to use it!
	virtual ARM_Object* Clone() const;
	/// stringify for easy debugging!
	virtual string toString(const string& indent="", const string& nextIndent="") const;
	/// specitic view
	virtual void View(char* id, FILE* ficOut) const;
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

