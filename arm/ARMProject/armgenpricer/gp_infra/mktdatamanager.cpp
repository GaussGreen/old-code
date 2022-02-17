/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file mktdatamanager.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpinfra/mktdatamanager.h"
#include "gpinfra/advanceduserstable.h"
#include "gpinfra/argconvdefault.h"
#include "gpinfra/mktdataretriever.h"
#include "gpbase/eventviewer.h"

#include "gpbase/singleton.h"
#include "gpbase/ostringstream.h"

#include <ccy/currency.h>

#include <typeinfo>
#include <iomanip>

CC_USING_NS(std,pair)


CC_BEGIN_NAMESPACE( ARM )

/// static as this is for all users!
bool ARM_MarketData_ManagerImp::isAuthorizedToChangedMarketData = true; ////ARM_MarketData_ManagerImp::GetIsAuthorizedTochangeMarketData();


/// we do not allow to get data .. the user should be responsible for 
/// 
bool ARM_MarketData_ManagerImp::UseMktDataRetriever = false;

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: ARM_MarketData_ManagerImp
///	Returns: Constructor
///	Action : build object
/////////////////////////////////////////////////////////////////
ARM_MarketData_ManagerImp::ARM_MarketData_ManagerImp(const ARM_Date& asOf )
:	itsAsOf(asOf), itsContents(), itsDetailViewMode(false)
{}


////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: ~ARM_MarketData_ManagerImp
///	Returns: Destructor
///	Action : destroy the object
/////////////////////////////////////////////////////////////////
ARM_MarketData_ManagerImp::~ARM_MarketData_ManagerImp()
{
	/// nothing as everything is handled by the ARM_ObjectPtr
}

////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: 
///	Returns: setAsOfDate
///	Action : if we completely change the asOfDate, the contents are erased!
/////////////////////////////////////////////////////////////////
void ARM_MarketData_ManagerImp::SetAsOfDate( const ARM_Date& newAsOf )
{
	if( newAsOf != itsAsOf )
	{
		itsContents = stringARM_ObjectPtrMap();
		itsAsOf		= newAsOf;
	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: RegisterData
///	Returns: register any ARM_Object data with a key
///	Action : 
/////////////////////////////////////////////////////////////////
string ARM_MarketData_ManagerImp::RegisterData(
	const string& objectKey, 
	ARM_Object* object )
{
	stringARM_ObjectPtrMap::iterator found = itsContents.find(objectKey);

	/// does it already exist
	if(found!=itsContents.end())
	{
		if(ARM_MarketData_ManagerImp::isAuthorizedToChangedMarketData)
		{
			if(!object)
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": try to insert a null object!" );

			/// we clone to avoid crashing!
			(*found).second= ARM_ObjectPtr( object->Clone() );
		}
		else
		{
			CC_Ostringstream os;
			os << ARM_USERNAME << ": you 're not authorized to overwrite data";
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
		}
	}
	else
	{
		if(!object)
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": try to insert a null object!" );
		itsContents.insert( pair< const string, ARM_ObjectPtr>( objectKey, ARM_ObjectPtr( object->Clone() ) ) );
	}

	return objectKey;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: RegisterData
///	Returns: register any ARM_Object data with a key
///	Action : 
/////////////////////////////////////////////////////////////////
string ARM_MarketData_ManagerImp::RegisterDataNoClone(
	const string& objectKey, 
	ARM_ObjectPtr object ) const
{
	stringARM_ObjectPtrMap::iterator found = CC_MUTABLE( ARM_MarketData_ManagerImp, itsContents ).find(objectKey);

	/// does it already exist
	if(found!=itsContents.end())
	{
		if(ARM_MarketData_ManagerImp::isAuthorizedToChangedMarketData)
		{
			/// we clone to avoid crashing!
			(*found).second= object;
		}
		else
		{
			CC_Ostringstream os;
			os << ARM_USERNAME << ": you 're not authorized to overwrite data";
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
		}
	}
	else
		CC_MUTABLE( ARM_MarketData_ManagerImp, itsContents ).insert( pair< const string, ARM_ObjectPtr>( objectKey, object ) );

	return objectKey;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: GetData
///	Returns: get any ARM_Object data with a key
///	Action : 
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_MarketData_ManagerImp::GetData(
	const string& objectKey ) const
{
	stringARM_ObjectPtrMap::const_iterator found = itsContents.find(objectKey);
	if(found!=itsContents.end())
		return &*(*found).second;
	else
		return NULL;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: GetData
///	Returns: get any ARM_Object data with a key and clone it
///	Action : 
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_MarketData_ManagerImp::GetDataAndClone( 
	const string& objectKey ) const
{
	ARM_Object* myObject= GetData( objectKey );

	if( !myObject )
		ThrowError_ObjectNoFound( objectKey, "any ARM_Object*" );
	
	/// other cases
	return myObject->Clone();
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: CoinZCObjKey
///	Returns: string
///	Action : create the object key for a zc curve
/////////////////////////////////////////////////////////////////
string ARM_MarketData_ManagerImp::CoinZCObjKey( 
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf, 
	const string& source )
{
	return string("ZC_") + ccy + "_" + indexName + "_" + cvName + "." + asOf.toString() + "." + source;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: ValidateAsOf
///	Returns: void
///	Action : Check that the given asOf is valid, otherwise throw exception
/////////////////////////////////////////////////////////////////
void ARM_MarketData_ManagerImp::ValidateAsOf( const ARM_Date& asOf ) const
{
	if( asOf != itsAsOf )
	{
		CC_Ostringstream os;
		os	<< ARM_USERNAME << ": incompatible asOF! Mkt Data with " << itsAsOf.toString()
			<< " while data asked for " << asOf.toString();
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: RegisterZCCurve
///	Returns: void
///	Action : Register a curve. a curve is registered by the following tag:
///		indexName_ccy_cvName__source_asOf
/////////////////////////////////////////////////////////////////
string ARM_MarketData_ManagerImp::RegisterZCCurve( 
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf, 
	const string& source, 
	ARM_ZeroCurve* curve )
{
	ValidateAsOf( asOf );
	return RegisterData( CoinZCObjKey( indexName, ccy, cvName, asOf, source ),NULL);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: ThrowError_ObjectNoFound
///	Returns: throw an exception (hence void)
///	Action : tells that an object could not be found!
/////////////////////////////////////////////////////////////////
void ARM_MarketData_ManagerImp::ThrowError_ObjectNoFound( 
	const string& key, 
	const string& type  ) const
{
	CC_Ostringstream os;
	os << ARM_USERNAME << ": could not find data " << type << " with key: " << key;
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: GetZCCurve
///	Returns: get a curve
///	Action : non const to allow to create a default curve!
/////////////////////////////////////////////////////////////////
ARM_ZeroCurve* ARM_MarketData_ManagerImp::GetZCCurve( 
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf, 
	const string& source ) const
{
	ARM_ZeroCurve* zc = NULL;//static_cast<ARM_ZeroCurve*>( GetData( CoinZCObjKey( indexName, ccy, cvName, asOf, source ) ) );
	if(!zc)
	{
		if( ARM_MarketData_ManagerImp::UseMktDataRetriever && ARM_MarketData_ManagerImp::isAuthorizedToChangedMarketData)
		{
			string key	= CoinZCObjKey( indexName, ccy, cvName, asOf, "Summit" );

			if( ARM_EventViewerImp::VerboseIsOn )
			{
				/// notice that we will use default curve!
				CC_Ostringstream os;
				os << "\nthe Market Data Manager will use a default curve for the zc curve object with key: " << key;
				ARM_TheEventViewer.Instance()->AddToMessage( os.str() );
			}

			zc	= ARM_TheMarketDataRetriever.Instance()->GetZCCurveFromSummit( indexName, ccy, cvName, asOf );

			if( ARM_EventViewerImp::VerboseIsOn )
				ARM_TheEventViewer.Instance()->AddToMessage( "\nloading of the zc curve was successful!\n" );

			/// and obviously register it!
			//RegisterDataNoClone( key, ARM_ObjectPtr( zc ) );
		}
		else
			ThrowError_ObjectNoFound( CoinZCObjKey( indexName, ccy, cvName, asOf, source ), "ZC Curve" );
	}

	return zc;
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: GetZCCurveAndClone
///	Returns: get a curve
///	Action : non const to allow to create a default curve!
/////////////////////////////////////////////////////////////////
ARM_ZeroCurve* ARM_MarketData_ManagerImp::GetZCCurveAndClone( 
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf, 
	const string& source ) const
{
	ARM_ZeroCurve* myCurve = GetZCCurve( indexName, ccy, cvName, asOf, source );

	if( !myCurve )
		ThrowError_ObjectNoFound( CoinZCObjKey( indexName, ccy, cvName, asOf, source ), "ZC Curve" );
	
	/// other cases
	return NULL;//(ARM_ZeroCurve* ) myCurve->Clone();
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: CoinZCObjKey
///	Returns: string
///	Action : create the object key for a vol curve
/////////////////////////////////////////////////////////////////
string ARM_MarketData_ManagerImp::CoinVolObjKey( 
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf,
	const string& volMktType,	/// IRG or SWOPT
	const string& volType,		/// ATM or Smile
	const string& source )
{
	return string("IRFWDVOL_") + volType + "_"+ ccy + "_" + indexName + "_" + volMktType + "_"+ cvName + "." + asOf.toString() + "." + source;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: RegisterVolCurve
///	Returns: register a volatility curve
///	Action : 
/////////////////////////////////////////////////////////////////
string ARM_MarketData_ManagerImp::RegisterVolCurve( 
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf,
	const string& volMktType,	/// IRG or SWOPT
	const string& volType,		/// ATM or Smile
	const string& source,
	ARM_VolCurve* volCurve	)
{
	ValidateAsOf( asOf );
	return RegisterData( CoinVolObjKey( indexName, ccy, cvName, asOf, volMktType, volType, source ), NULL);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: GetVolCurve
///	Returns: ARM_VolCurve*
///	Action : get a vol curve
/////////////////////////////////////////////////////////////////
ARM_VolCurve* ARM_MarketData_ManagerImp::GetVolCurve( 
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf,
	const string& volMktType,		/// IRG or SWOPT
	const string& volType,			/// ATM or Smile
	const string& source ) const
{
	ARM_VolCurve* volCurve	= NULL;//static_cast<ARM_VolCurve*>( GetData( CoinVolObjKey( indexName, ccy, cvName, asOf, volMktType, volType, source ) ) );
	if(!volCurve)
	{
		if( ARM_MarketData_ManagerImp::UseMktDataRetriever && ARM_MarketData_ManagerImp::isAuthorizedToChangedMarketData)
		{
			string key	= CoinVolObjKey( indexName, ccy, cvName, asOf, volMktType, volType, "Summit" );

			if( ARM_EventViewerImp::VerboseIsOn )
			{
				/// notice that we will use default curve!
				CC_Ostringstream os;
				os << "\nthe Market Data Manager will the database for the vol curve object with key: " << key;
				ARM_TheEventViewer.Instance()->AddToMessage( os.str() );
			}

			/// Get the default vol from Summit
			volCurve	= NULL;/*ARM_TheMarketDataRetriever.Instance()->GetVolFromSummit( 
				indexName, 
				ccy, 
				cvName, 
				asOf, 
				(ARM_MarketData::VolMktType) ARM_ArgConv_VolMktType.GetNumber( volMktType ), 
				(ARM_MarketData::VolType)    ARM_ArgConv_VolType.GetNumber(    volType    ) );*/

			if( ARM_EventViewerImp::VerboseIsOn )
				ARM_TheEventViewer.Instance()->AddToMessage( "\nloading of the vol curve was successful!\n" );

			/// and obviously register it!
			//RegisterDataNoClone( key, ARM_ObjectPtr( volCurve ) );
		}
		else
			ThrowError_ObjectNoFound( CoinVolObjKey( indexName, ccy, cvName, asOf, volMktType, volType, source ), "VolCurve" );
	}
	return volCurve;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: GetVolCurveAndClone
///	Returns: ARM_VolCurve*
///	Action : get a vol curve
/////////////////////////////////////////////////////////////////
ARM_VolCurve* ARM_MarketData_ManagerImp::GetVolCurveAndClone( 
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf,
	const string& volMktType,	/// IRG or SWOPT
	const string& volType,		/// ATM or Smile
	const string& source ) const
{
	ARM_VolCurve* myVolCurve = GetVolCurve( indexName, ccy, cvName, asOf, volMktType, volType, source );

	if( !myVolCurve )
		ThrowError_ObjectNoFound( CoinVolObjKey( indexName, ccy, cvName, asOf, volMktType, volType, source ), "VolCurve" );

	/// other cases
	return NULL;// (ARM_VolCurve*) myVolCurve->Clone();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: CoinZCObjKey
///	Returns: string
///	Action : create the object key for a vol curve
/////////////////////////////////////////////////////////////////
string ARM_MarketData_ManagerImp::CoinVolMktModelObjKey( 
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf,
	const string& volMktType,	/// IRG or SWOPT
	const string& source )
{
	return string("VOLMKTMODEL_") + ccy + "_" + indexName + "_" + volMktType + "_"+ cvName + "." + asOf.toString() + "." + source;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: RegisterVolMktModel
///	Returns: register a volatility mkt model
///	Action : 
/////////////////////////////////////////////////////////////////
string ARM_MarketData_ManagerImp::RegisterVolMktModel( 
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf,
	const string& volMktType,	/// IRG or SWOPT
	const string& source,
	ARM_Model* volMktModel	)
{
	ValidateAsOf( asOf );
	return RegisterData( CoinVolMktModelObjKey( indexName, ccy, cvName, asOf, volMktType, source ), NULL);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: GetVolMktModel
///	Returns: ARM_Model*
///	Action : return mkt vol model
/////////////////////////////////////////////////////////////////
ARM_Model* ARM_MarketData_ManagerImp::GetVolMktModel( 
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf,
	const string& volMktType,	/// IRG or SWOPT
	const string& source ) const
{
	ARM_Model* volMktModel = NULL;// static_cast<ARM_Model*>( GetData( CoinVolMktModelObjKey( indexName, ccy, cvName, asOf, volMktType, source ) ) );
	if(!volMktModel)
	{
		if( ARM_MarketData_ManagerImp::UseMktDataRetriever && ARM_MarketData_ManagerImp::isAuthorizedToChangedMarketData)
		{
			volMktModel	= ARM_TheMarketDataRetriever.Instance()->GetVolMktModel( indexName,	ccy, 
				cvName, asOf, volMktType, source );
			string key	= CoinVolMktModelObjKey( indexName, ccy, cvName, asOf, volMktType, "Summit" );

			/// and obviously register it!
			//RegisterDataNoClone( key, ARM_ObjectPtr( volMktModel ) );
		}
		else
			ThrowError_ObjectNoFound( CoinVolMktModelObjKey( indexName, ccy, cvName, asOf, volMktType, source ), "VolCurve" );
	}
	return volMktModel;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: GetVolMktModelAndClone
///	Returns: ARM_Model*
///	Action : return mkt vol model
/////////////////////////////////////////////////////////////////
ARM_Model* ARM_MarketData_ManagerImp::GetVolMktModelAndClone( 
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf,
	const string& volMktType,	/// IRG or SWOPT
	const string& source ) const
{
	return NULL;
	/*ARM_Model* myVolMktModel =GetVolMktModel(indexName, ccy, cvName, asOf, volMktType, source );

	if( !myVolMktModel )
		ThrowError_ObjectNoFound( CoinVolMktModelObjKey( indexName, ccy, cvName, asOf, volMktType, source ), "VolCurve" );

	/// other cases
	return (ARM_Model*) myVolMktModel->Clone();*/
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: GetIsAuthorizedTochangeMarketData
///	Returns: check that a user is authorized to change market data
///	Action : 
/////////////////////////////////////////////////////////////////
bool ARM_MarketData_ManagerImp::GetIsAuthorizedTochangeMarketData()
{
	string userName = ARM_USERNAME;
	return UserControl_IsInTheList(GP_User_Level2_Table,GP_User_Level2_Table_Size,userName );
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: toString 
///	Returns: string describing the content of the ARM_MarketData_ManagerImp
///	Action : 
/////////////////////////////////////////////////////////////////
string ARM_MarketData_ManagerImp::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	/*os << "\n\n";
	os << indent << "Market Data Manager\n";
	os << indent << "------------------\n\n";
	os << indent << "As Of Date   : " << itsAsOf.toString() << "\n";
	os << indent << "Detailed Mode: " << (itsDetailViewMode? "on" : "off" ) << "\n";

	stringARM_ObjectPtrMap::const_iterator iter = itsContents.begin(),
		end = itsContents.end();

	if(!itsContents.empty())
	{
		os  << "\n" << indent << CC_NS(std,setw)(50) << CC_NS(std,left) << "key" << "\t" 
			<< "Object type\n";
		os << indent << "-------------------------------------------------------------------\n";

		string name;

		while( iter!= end )
		{
			name = typeid( (*iter).second ).name();

			/// is it in fact a zero curve, a vol curve, a model
			if( dynamic_cast<ARM_ZeroCurve*>( &*(*iter).second ) )
				name = "class ARM_ZeroCurve *";
			else
				if( dynamic_cast<ARM_VolCurve*>( &*(*iter).second ) )
					name = "class ARM_VolCurve *";
				else
					if( dynamic_cast<ARM_Model*>( &*(*iter).second ) )
						name = "class ARM_Model *";
			
			os  << indent << CC_NS(std,setw)(50) << CC_NS(std,left) << (*iter).first << "\t" 
				<< name << "\n";
			++iter;
		}
	}*/
	return os.str();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp 
///	Routine: View
///	Returns: void
///	Action :
/////////////////////////////////////////////////////////////////
void ARM_MarketData_ManagerImp::View(char* id, FILE* ficOut) const
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
    fprintf(fOut, "%s", toString().c_str() );

	if( itsDetailViewMode )
	{
		stringARM_ObjectPtrMap::iterator iter = itsContents.begin(),
			end = itsContents.end();

		while( iter != end )
		{
			((*iter).second)->View(id,fOut);
			++iter;
		}
	}


	/// to allow to have nested view
    if ( NULL == ficOut )
       fclose(fOut);
}



////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: 
///	Returns: ResetAllData
///	Action : removed all data!
/////////////////////////////////////////////////////////////////
void ARM_MarketData_ManagerImp::ResetAllData()
{
	itsContents = stringARM_ObjectPtrMap();
}


////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: CopyClonedData
///	Returns: void
///	Action : copy and clone all data
/////////////////////////////////////////////////////////////////
void ARM_MarketData_ManagerImp::CopyClonedData( const ARM_MarketData_ManagerImp& rhs )
{
	stringARM_ObjectPtrMap::iterator 
		rhsIter = rhs.itsContents.begin(), 
		end		= rhs.itsContents.end();
	
	while( rhsIter != end )
	{
		itsContents.insert( pair< const string, ARM_ObjectPtr>( (*rhsIter).first, ARM_ObjectPtr( ((*rhsIter).second)->Clone() ) ) );
		++rhsIter;
	}
}


////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: CleanData
///	Returns: void
///	Action : just erase everything in itsContents using the fact 
///				that we are using smart pointors
/////////////////////////////////////////////////////////////////
void ARM_MarketData_ManagerImp::CleanData()
{
	/// clean everything as we are using smart pointors!
	itsContents			= stringARM_ObjectPtrMap();
}



////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: ARM_MarketData_ManagerImp copy constructor
///	Returns: 
///	Action : standard copy constructor
/////////////////////////////////////////////////////////////////
ARM_MarketData_ManagerImp::ARM_MarketData_ManagerImp( const ARM_MarketData_ManagerImp& rhs )
:	ARM_RootObject( rhs ), itsAsOf( rhs.itsAsOf ), itsContents(), itsDetailViewMode( rhs.itsDetailViewMode )
{
	CopyClonedData(rhs);
}


////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: assignment operator
///	Returns: 
///	Action : standard operator=
/////////////////////////////////////////////////////////////////
ARM_MarketData_ManagerImp& ARM_MarketData_ManagerImp::operator=(const ARM_MarketData_ManagerImp& rhs)
{
	if( this != &rhs )
	{
		ARM_RootObject::operator=( rhs );
		itsAsOf				= rhs.itsAsOf;
		itsDetailViewMode	= rhs.itsDetailViewMode;
		CleanData();
		CopyClonedData(rhs);
	}
	return *this;
}


////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketData_ManagerImp
///	Routine: 
///	Returns: ResetAllData
///	Action : removed all data!
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_MarketData_ManagerImp::Clone() const
{
	return new ARM_MarketData_ManagerImp(*this); 
}



/// Creation of the market data retriever
ARM_SingletonHolder<ARM_MarketData_ManagerImp> ARM_TheMarketData_Manager;

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

