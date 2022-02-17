#include "gpbase/removeidentifiedwarning.h"

#include "gpinfra/mktdatamanagerrep.h"
#include "gpinfra/mktdatamanager.h"
#include "gpbase/singleton.h"
#include "gpbase/ostringstream.h"

CC_BEGIN_NAMESPACE( ARM )

/////////////////////////////////////////////////////////////////
/// strategy
/////////////////////////////////////////////////////////////////
#define __GP_MKTDATAMGR_CLONE

#if defined(__GP_MKTDATAMGR_CLONE)
	MktDataManagerFunctionDateToVoid  ARM_MarketData_ManagerRepStrategy::CreateMktDataInstanceFunc			= &ARM_MarketData_ManagerRep::CreateMktDataInstance_Clone;
	MktDataManagerFunctionVoidToVoid  ARM_MarketData_ManagerRepStrategy::CleanMktDataInstanceFunc			= &ARM_MarketData_ManagerRep::CleanMktDataInstance_Clone;
	MktDataManagerFunctionMktDataRepToVoid  ARM_MarketData_ManagerRepStrategy::CopyMktDataInstanceFunc		= &ARM_MarketData_ManagerRep::CopyMktDataInstance_Clone;
	MktDataManagerFunctionVoidToMktDataImpCst  ARM_MarketData_ManagerRepStrategy::GetMktDataInstanceFunc	= &ARM_MarketData_ManagerRep::GetMktDataInstance_Clone;
#else
	MktDataManagerFunctionDateToVoid  ARM_MarketData_ManagerRepStrategy::CreateMktDataInstanceFunc			= ARM_MarketData_ManagerRep::CreateMktDataInstance_Singleton;
	MktDataManagerFunctionVoidToVoid  ARM_MarketData_ManagerRepStrategy::CleanMktDataInstanceFunc			= ARM_MarketData_ManagerRep::CleanMktDataInstance_Singleton;
	MktDataManagerFunctionMktDataRepToVoid  ARM_MarketData_ManagerRepStrategy::CopyMktDataInstanceFunc		= ARM_MarketData_ManagerRep::CopyMktDataInstance_Singleton;
	MktDataManagerFunctionVoidToMktDataImpCst  ARM_MarketData_ManagerRepStrategy::GetMktDataInstanceFunc	= ARM_MarketData_ManagerRep::GetMktDataInstance_Singleton;
#endif



/////////////////////////////////////////////////////////////////
// Singleton behavior
/////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerRep
///	Routine: CreateMktDataInstance_Singleton,CleanMktDataInstance_Singleton,
///				CopyMktDataInstance_Singleton,GetMktDataInstance_Singleton
///	Action : just redirect to the singleton to the singleton!
/////////////////////////////////////////////////////////////////

void ARM_MarketData_ManagerRep::CreateMktDataInstance_Singleton( const ARM_Date& asOf)
{
	((this->*ARM_MarketData_ManagerRepStrategy::GetMktDataInstanceFunc)())->SetAsOfDate( asOf );
}

void ARM_MarketData_ManagerRep::CleanMktDataInstance_Singleton()
{}

void ARM_MarketData_ManagerRep::CopyMktDataInstance_Singleton( const ARM_MarketData_ManagerRep& rhs )
{}

ARM_MarketData_ManagerImp* ARM_MarketData_ManagerRep::GetMktDataInstance_Singleton() const
{
	return ARM_TheMarketData_Manager.Instance(); 
}


/////////////////////////////////////////////////////////////////
// Clone behavior
/////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerRep
///	Routine: CreateMktDataInstance_Clone,CleanMktDataInstance_Clone,
///				CopyMktDataInstance_Clone,GetMktDataInstance_Clone
///	Action : create, copy and clone the mkt data implementation
/////////////////////////////////////////////////////////////////

void ARM_MarketData_ManagerRep::CreateMktDataInstance_Clone( const ARM_Date& asOf )
{ 
	itsMktDataManager = new ARM_MarketData_ManagerImp(asOf); 
}

void ARM_MarketData_ManagerRep::CleanMktDataInstance_Clone()
{
	delete itsMktDataManager; 
}

void ARM_MarketData_ManagerRep::CopyMktDataInstance_Clone( const ARM_MarketData_ManagerRep& rhs )
{	
	itsMktDataManager = rhs.itsMktDataManager? (ARM_MarketData_ManagerImp*) rhs.itsMktDataManager->Clone() : NULL; 
}

ARM_MarketData_ManagerImp* ARM_MarketData_ManagerRep::GetMktDataInstance_Clone() const
{
	return itsMktDataManager;
}




/////////////////////////////////////////////////////////////////
/// ARM_MarketData_ManagerRep
/////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////
///	Class  : 
///	Routine: constructor,copy constructor,assignment operator and destructor
///	Returns: 
///	Action : standard action for derived object
/////////////////////////////////////////////////////////////////

ARM_MarketData_ManagerRep::ARM_MarketData_ManagerRep( const ARM_Date& asOf)
:	ARM_RootObject(), itsKeys(), itsMktDataManager()
{ 	
	/// forced to use this as this refers to the current class object!
	(this->*ARM_MarketData_ManagerRepStrategy::CreateMktDataInstanceFunc)(asOf);
	CC_ARM_SETNAME(ARM_MKTDATAMANAGER); 
}

ARM_MarketData_ManagerRep::ARM_MarketData_ManagerRep( const ARM_MarketData_ManagerRep& rhs ) 
:	ARM_RootObject( rhs ), itsKeys( rhs.itsKeys )
{
	/// forced to use this as this refers to the current class object!
	(this->*ARM_MarketData_ManagerRepStrategy::CopyMktDataInstanceFunc)(rhs);
}

ARM_MarketData_ManagerRep& ARM_MarketData_ManagerRep::operator=( const ARM_MarketData_ManagerRep& rhs )
{
	if( this != &rhs )
	{
		ARM_RootObject::operator=(rhs);
		/// forced to use this as this refers to the current class object!
		(this->*ARM_MarketData_ManagerRepStrategy::CleanMktDataInstanceFunc)();
		(this->*ARM_MarketData_ManagerRepStrategy::CopyMktDataInstanceFunc)(rhs);
		itsKeys	= rhs.itsKeys;
	}
	return *this;
}


ARM_MarketData_ManagerRep::~ARM_MarketData_ManagerRep()
{
	(this->*ARM_MarketData_ManagerRepStrategy::CleanMktDataInstanceFunc)();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerRep
///	Routine: Clone
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_MarketData_ManagerRep::Clone() const
{
	return new ARM_MarketData_ManagerRep(*this);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerRep
///	Routine: ThrowErrorIfKeyMissing
///	Returns: void
///	Action : throw an exception if the key is missing!
/////////////////////////////////////////////////////////////////
void ARM_MarketData_ManagerRep::ThrowErrorIfKeyMissing( const string& key ) const
{
	if( itsKeys.find( key ) == itsKeys.end() )
       throw Exception(__LINE__, __FILE__,ARM_FRMMODEL, ARM_USERNAME + ": " + key + " missing! please advise!" );    
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerRep
///	Routine: TestIfKeyMissing
///	Returns: void
///	Action : test if the key is missing
/////////////////////////////////////////////////////////////////
bool ARM_MarketData_ManagerRep::TestIfKeyMissing( const string& key ) const
{
    if( itsKeys.find( key ) == itsKeys.end() )
        return true;
    else
        return false;
}


////////////////////////////////////////////////////
///	Class   : ARM_MarketData_ManagerRep
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_MarketData_ManagerRep::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;

	/// show keys
    os << "\n\n";
    os << indent << "ARM_MarketData_ManagerRep\n";
    os << indent << "-------------------------\n";
    os << indent << "Keys = ";
	CC_NS(std,set)<string>::const_iterator iter = itsKeys.begin();
	while( iter != itsKeys.end() )
	{
		os << (*iter) << "   ";
		++iter;
	}

	os << ((this->*ARM_MarketData_ManagerRepStrategy::GetMktDataInstanceFunc)())->toString(indent);
    return os.str();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_RootObject 
///	Routine: View
///	Returns: void
///	Action :
/////////////////////////////////////////////////////////////////
void ARM_MarketData_ManagerRep::View(char* id, FILE* ficOut) const
{
    FILE* fOut;
    char fOutName[200];
	
	/// first determine that the file is not already opened
    if ( NULL == ficOut )
    {
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		(void) unlink(fOutName);
		fOut = fopen(fOutName, "w");
    }
    else
		fOut = ficOut;

	/// use the method to string
	/// and just says the type and what is in it
	string allErrors = toString();
    fprintf(fOut, "%s", allErrors.c_str() );

	((this->*ARM_MarketData_ManagerRepStrategy::GetMktDataInstanceFunc)())->View(id,fOut);

	/// to allow to have nested view
    if ( NULL == ficOut )
       fclose(fOut);
}


////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerRep
///	Routine: 
///	Returns: setAsOfDate
///	Action : if we completely change the asOfDate, the contents are erased!
/////////////////////////////////////////////////////////////////
void ARM_MarketData_ManagerRep::SetAsOfDate( const ARM_Date& newAsOf ) const
{
	((this->*ARM_MarketData_ManagerRepStrategy::GetMktDataInstanceFunc)())->SetAsOfDate( newAsOf );
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerRep
///	Routine: RegisterData
///	Returns: register any ARM_Object data with a key
///	Action : 
/////////////////////////////////////////////////////////////////
string ARM_MarketData_ManagerRep::RegisterData(
	const string& objectKey, 
	ARM_Object* object )
{
	/// registers the key and delegates
	itsKeys.insert( objectKey );

	return ((this->*ARM_MarketData_ManagerRepStrategy::GetMktDataInstanceFunc)())->RegisterData( objectKey, object );
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerRep
///	Routine: GetData
///	Returns: get any ARM_Object data with a key
///	Action : 
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_MarketData_ManagerRep::GetData(
	const string& objectKey ) const
{
	/// validates the key
	ThrowErrorIfKeyMissing( objectKey );

	return ((this->*ARM_MarketData_ManagerRepStrategy::GetMktDataInstanceFunc)())->GetData( objectKey );
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerRep
///	Routine: GetData
///	Returns: get any ARM_Object data with a key and clone it
///	Action : 
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_MarketData_ManagerRep::GetDataAndClone( 
	const string& objectKey ) const
{
	/// validates the key
	ThrowErrorIfKeyMissing( objectKey );

	return ((this->*ARM_MarketData_ManagerRepStrategy::GetMktDataInstanceFunc)())->GetDataAndClone( objectKey );
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerRep
///	Routine: RegisterZCCurve
///	Returns: void
///	Action : Register a curve. a curve is registered by the following tag:
///		indexName_ccy_cvName__source_asOf
/////////////////////////////////////////////////////////////////
string ARM_MarketData_ManagerRep::RegisterZCCurve( 
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf, 
	const string& source, 
	ARM_ZeroCurve* curve )
{
	/// registers the key and delegates
	itsKeys.insert(ARM_MarketData_ManagerImp::CoinZCObjKey(indexName, ccy, cvName, asOf, source ) );
	
	return ((this->*ARM_MarketData_ManagerRepStrategy::GetMktDataInstanceFunc)())->RegisterZCCurve( 
		indexName, ccy, cvName, asOf, source, curve );
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerRep
///	Routine: GetZCCurve
///	Returns: get a curve
///	Action : non const to allow to create a default curve!
/////////////////////////////////////////////////////////////////
ARM_ZeroCurve* ARM_MarketData_ManagerRep::GetZCCurve( 
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf, 
	const string& source ) const
{
	/// validates the key
	ThrowErrorIfKeyMissing( ARM_MarketData_ManagerImp::CoinZCObjKey(indexName, ccy, cvName, asOf, source ) );
		
	return ((this->*ARM_MarketData_ManagerRepStrategy::GetMktDataInstanceFunc)())->GetZCCurve( 
		indexName, ccy, cvName, asOf, source );
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerRep
///	Routine: GetZCCurveAndClone
///	Returns: get a curve
///	Action : non const to allow to create a default curve!
/////////////////////////////////////////////////////////////////
ARM_ZeroCurve* ARM_MarketData_ManagerRep::GetZCCurveAndClone( 
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf, 
	const string& source ) const
{
	/// validates the key
	ThrowErrorIfKeyMissing( ARM_MarketData_ManagerImp::CoinZCObjKey(indexName, ccy, cvName, asOf, source ) );

	return ((this->*ARM_MarketData_ManagerRepStrategy::GetMktDataInstanceFunc)())->GetZCCurveAndClone( 
		indexName, ccy, cvName, asOf, source );
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerRep
///	Routine: RegisterVolCurve
///	Returns: register a volatility curve
///	Action : 
/////////////////////////////////////////////////////////////////
string ARM_MarketData_ManagerRep::RegisterVolCurve( 
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf,
	const string& volMktType,	/// IRG or SWOPT
	const string& volType,		/// ATM or Smile
	const string& source,
	ARM_VolCurve* volCurve	)
{
	/// registers the key and delegates
	itsKeys.insert( ARM_MarketData_ManagerImp::CoinVolObjKey( indexName, ccy, cvName, asOf, volMktType, volType, source ) );

	return ((this->*ARM_MarketData_ManagerRepStrategy::GetMktDataInstanceFunc)())->RegisterVolCurve( 
		indexName, ccy, cvName, asOf, volMktType, volType, source, volCurve );
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerRep
///	Routine: GetVolCurve
///	Returns: ARM_VolCurve*
///	Action : get a vol curve
/////////////////////////////////////////////////////////////////
ARM_VolCurve* ARM_MarketData_ManagerRep::GetVolCurve( 
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf,
	const string& volMktType,		/// IRG or SWOPT
	const string& volType,			/// ATM or Smile
	const string& source ) const
{
	/// validates the key
	ThrowErrorIfKeyMissing( ARM_MarketData_ManagerImp::CoinVolObjKey( indexName, ccy, cvName, asOf, volMktType, volType, source ) );

	return ((this->*ARM_MarketData_ManagerRepStrategy::GetMktDataInstanceFunc)())->GetVolCurve( 
		indexName, ccy, cvName, asOf, volMktType, volType, source );
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerRep
///	Routine: GetVolCurveAndClone
///	Returns: ARM_VolCurve*
///	Action : get a vol curve
/////////////////////////////////////////////////////////////////
ARM_VolCurve* ARM_MarketData_ManagerRep::GetVolCurveAndClone( 
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf,
	const string& volMktType,	/// IRG or SWOPT
	const string& volType,		/// ATM or Smile
	const string& source ) const
{
	/// validates the key
	ThrowErrorIfKeyMissing( ARM_MarketData_ManagerImp::CoinVolObjKey( indexName, ccy, cvName, asOf, volMktType, volType, source ) );

	return ((this->*ARM_MarketData_ManagerRepStrategy::GetMktDataInstanceFunc)())->GetVolCurveAndClone( 
		indexName, ccy, cvName, asOf, volMktType, volType, source );
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerRep
///	Routine: RegisterVolMktModel
///	Returns: register a volatility mkt model
///	Action : 
/////////////////////////////////////////////////////////////////
string ARM_MarketData_ManagerRep::RegisterVolMktModel( 
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf,
	const string& volMktType,	/// IRG or SWOPT
	const string& source,
	ARM_Model* volMktModel	)
{
	/// registers the key and delegates
	itsKeys.insert( ARM_MarketData_ManagerImp::CoinVolMktModelObjKey( indexName, ccy, cvName, asOf, volMktType, source ) );

	return ((this->*ARM_MarketData_ManagerRepStrategy::GetMktDataInstanceFunc)())->RegisterVolMktModel( 
		indexName, ccy, cvName, asOf, volMktType, source, volMktModel );
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerRep
///	Routine: GetVolMktModel
///	Returns: ARM_Model*
///	Action : return mkt vol model
/////////////////////////////////////////////////////////////////
ARM_Model* ARM_MarketData_ManagerRep::GetVolMktModel( 
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf,
	const string& volMktType,	/// IRG or SWOPT
	const string& source ) const
{
	/// validates the key
	ThrowErrorIfKeyMissing( ARM_MarketData_ManagerImp::CoinVolMktModelObjKey( indexName, ccy, cvName, asOf, volMktType, source ) );

	return ((this->*ARM_MarketData_ManagerRepStrategy::GetMktDataInstanceFunc)())->GetVolMktModel( 
		indexName, ccy, cvName, asOf, volMktType, source  );
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerRep
///	Routine: GetVolMktModelAndClone
///	Returns: ARM_Model*
///	Action : return mkt vol model
/////////////////////////////////////////////////////////////////
ARM_Model* ARM_MarketData_ManagerRep::GetVolMktModelAndClone( 
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf,
	const string& volMktType,	/// IRG or SWOPT
	const string& source ) const
{
	/// validates the key
	ThrowErrorIfKeyMissing( ARM_MarketData_ManagerImp::CoinVolMktModelObjKey( indexName, ccy, cvName, asOf, volMktType, source ) );

	return ((this->*ARM_MarketData_ManagerRepStrategy::GetMktDataInstanceFunc)())->GetVolMktModelAndClone( 
		indexName, ccy, cvName, asOf, volMktType, source  );
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerRep
///	Routine: SetDetailMode
///	Returns: bool the setted value
///	Action : change the view mode of the market data manager
/////////////////////////////////////////////////////////////////
bool ARM_MarketData_ManagerRep::SetDetailMode( bool value ) const
{
	return ((this->*ARM_MarketData_ManagerRepStrategy::GetMktDataInstanceFunc)())->SetDetailMode( value );
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerRep
///	Routine: GetDetailMode
///	Returns: bool the return value
///	Action : return the view mode of the market data manager
/////////////////////////////////////////////////////////////////
bool ARM_MarketData_ManagerRep::GetDetailMode() const
{
	return ((this->*ARM_MarketData_ManagerRepStrategy::GetMktDataInstanceFunc)())->GetDetailMode();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerRep
///	Routine: GetAsOfDate
///	Returns: const ARM_Date& 
///	Action : accessor to the asOfDate
/////////////////////////////////////////////////////////////////
const ARM_Date& ARM_MarketData_ManagerRep::GetAsOfDate() const
{
	return ((this->*ARM_MarketData_ManagerRepStrategy::GetMktDataInstanceFunc)())->GetAsOfDate();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerRep
///	Routine: ResetAllData
///	Returns: void
///	Action : reset all data
/////////////////////////////////////////////////////////////////
void ARM_MarketData_ManagerRep::ResetAllData()
{
	ResetMyData();
	((this->*ARM_MarketData_ManagerRepStrategy::GetMktDataInstanceFunc)())->ResetAllData();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerRep
///	Routine: ResetMyData
///	Returns: void
///	Action : reset only the current key
/////////////////////////////////////////////////////////////////
void ARM_MarketData_ManagerRep::ResetMyData()
{
	itsKeys = set<string>();
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

