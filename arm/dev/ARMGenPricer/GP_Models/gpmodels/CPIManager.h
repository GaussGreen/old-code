/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file CPIMAnager.h
 *
 *  \brief prototype model for the generic pricer
 *
 *	\author  A. Schauly
 *	\version 1.0
 *	\date March 2005
 */


#include "gpbase/gpvector.h"
#include "gpbase/typedef.h"

/// includes from stl
#include <map>
#include <utility>


CC_BEGIN_NAMESPACE( ARM )

struct CPIInfo
{
	string itsInfcurveName;
	string itsDCFLag;
	long itsDailyInterp;
	string itsResetLag;

	CPIInfo( const string& InfcurveName, const string& DCFLag, const string& resetLag, long dailyInterp ) : 
	itsInfcurveName(InfcurveName), itsDCFLag(DCFLag), itsResetLag(resetLag), itsDailyInterp(dailyInterp) {}
	CPIInfo( ) : itsInfcurveName(""), itsDCFLag(""), itsResetLag(""), itsDailyInterp(0) {}

	bool operator == (CPIInfo& rhs) 
	{ 
		return ( (itsInfcurveName == rhs.itsInfcurveName) && (itsDCFLag == rhs.itsDCFLag ) && 
			(itsDailyInterp == rhs.itsDailyInterp ) &&(itsResetLag ==rhs.itsResetLag) );
	}

	bool operator != (CPIInfo& rhs) 
	{
		return !((*this) == rhs );
	}
};

class CPIManager;

typedef ARM_CountedPtr<CPIInfo> CPIInfoPtr;
typedef ARM_CountedPtr<CPIManager>	CPIManagerPtr;

/// \class CPIManager
/// Helps to manage CPIs in the past
class CPIManager
{
private: 

/// Maps and pairs used to store information. 
	typedef pair<long,ARM_GP_VectorPtr> CPIandItsRefNb_t;
	typedef map<double,CPIandItsRefNb_t> CPIMap_t;
	typedef pair<double,CPIandItsRefNb_t> CPIMapElt;
	typedef CPIMap_t::iterator CPIMapIterator; 
	typedef map<double,CPIInfoPtr> CPIInfoMap_t;
	typedef pair<double,CPIInfoPtr> CPIInfoMapElt;
	typedef map<double,long> CPIRefsNbMap_t;
	typedef pair<double,long> CPIRefsNbMapElt;
	typedef CPIRefsNbMap_t::iterator CPIRefsIterator; 

/// Information stored :
/// CPIs values during pricing ; How many times CPIs are referenced (useful for AMC) ; Infos on CPIs

	CPIMap_t itsCPIs;
	CPIRefsNbMap_t itsCPIRefsNb;
	CPIInfoMap_t itsCPIInfos;

public: 
	/// constructor/destructor
	CPIManager() {}
	~CPIManager() {}

	CPIManager( const CPIManager& rhs ) { itsCPIs = rhs.itsCPIs; }

	inline CPIInfoPtr getCPIInfo( double CPIDate )
	{
		CPIInfoMap_t::iterator iter = itsCPIInfos.find( CPIDate );

		if ( iter != itsCPIInfos.end() )
			return iter->second;
		else
			return CPIInfoPtr(NULL);
	}

	inline void EnqueueCPIInfo( double CPIDate, const CPIInfoPtr& CpiInfo )
	{
		CPIInfoMap_t::iterator iter = itsCPIInfos.find(CPIDate); 

		if ( iter == itsCPIInfos.end() )
			itsCPIInfos.insert( CPIInfoMapElt( CPIDate,CpiInfo ) );
		else
			if( *(iter->second) != *CpiInfo )
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": CPIDate Already in the map, for a CPI with different parameters. Please Advise. " );

	}

	inline void deleteCPIInfo( double CPIDate )
	{
		CPIInfoMap_t::iterator iter = itsCPIInfos.find(CPIDate); 

		if ( iter != itsCPIInfos.end() )
			itsCPIInfos.erase(iter);
	}

	inline void setCounter( double CPIDate, long counter )
	{ 
		CPIMapIterator iter = itsCPIs.find( CPIDate );  
		CPIRefsIterator iter2 = itsCPIRefsNb.find( CPIDate );
		iter->second.first = counter;
		iter2->second = counter;
	}
	
	inline void increaseCounter( double CPIDate ) 
	{ 
		CPIMapIterator iter = itsCPIs.find( CPIDate ); 
		CPIRefsIterator iter2 = itsCPIRefsNb.find( CPIDate );
		iter->second.first++; 
		iter2->second = ( iter2->second > iter->second.first ) ? iter2->second : iter->second.first;
	}

	inline void increaseCounterIfEqualsOne( double CPIDate ) 
	{ 
		CPIMapIterator iter = itsCPIs.find( CPIDate ); 
		long counter = iter->second.first+1;
		if ( counter == 2 )
			iter->second.first = counter; 
	}

	inline void decreaseCounterAndDeleteAllInformation( double CPIDate ) 
	{ 
		long newcounter;
		CPIMapIterator iter = itsCPIs.find( CPIDate ); 

#if defined(__GP_STRICT_VALIDATION)
		if( iter == itsCPIs.end() )
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": CPIDate Not found in the map " );
#endif

		newcounter = iter->second.first-1; 

		if( newcounter > 0 )
			iter->second.first = newcounter;
		else
		{
			CPIInfoMap_t::iterator cpiInfoIterator = itsCPIInfos.find( CPIDate );
			if(cpiInfoIterator != itsCPIInfos.end())
				itsCPIInfos.erase(cpiInfoIterator);

			itsCPIs.erase( iter );
		}
	}

	inline void decreaseCounter( double CPIDate ) 
	{ 
		long newcounter;
		CPIMapIterator iter = itsCPIs.find( CPIDate ); 

#if defined(__GP_STRICT_VALIDATION)
		if( iter == itsCPIs.end() )
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": CPIDate Not found in the map " );
#endif

		newcounter = iter->second.first-1; 

		if( newcounter > 0 )
			iter->second.first = newcounter;
	}

	inline void decreaseCounterAndDeleteCPI( double CPIDate ) 
	{ 
		long newcounter;
		CPIMapIterator iter = itsCPIs.find( CPIDate ); 

#if defined(__GP_STRICT_VALIDATION)
		if( iter == itsCPIs.end() )
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": CPIDate Not found in the map " );
#endif

		newcounter = iter->second.first-1; 

		if( newcounter > 0 )
			iter->second.first = newcounter;
		else
			itsCPIs.erase( iter );

	}

	inline void EnqueueCPI( double CPIDate, const ARM_GP_VectorPtr& CPIVector ) 
	{ 
		CPIMapIterator iter = itsCPIs.find(CPIDate); 

		if ( iter == itsCPIs.end() ) 
			itsCPIs.insert( CPIMapElt( CPIDate,CPIandItsRefNb_t(1,CPIVector) ) ); 
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": CPIDate Already in the map. Please Advise. " );
	}

	inline bool doesExist( double CPIDate )
	{
		CPIMapIterator iter = itsCPIs.find(CPIDate); 

		if ( iter != itsCPIs.end() ) 
			return true;
		else
			return false;
	}

	inline void setCPI( double CPIDate, const ARM_GP_VectorPtr& CPIVector ) 
	{ 
		CPIMapIterator iter = itsCPIs.find(CPIDate);
		if ( iter != itsCPIs.end() )
				iter->second.second = CPIVector;
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": CPIDate not in the map. Please Advise. " );
	}

	inline void setCPIAndCheckCounter( double CPIDate, const ARM_GP_VectorPtr& CPIVector ) 
	{ 
		CPIMapIterator iter = itsCPIs.find(CPIDate);
		if ( iter != itsCPIs.end() )
		{
			if (iter->second.first > 1)
				iter->second.second = CPIVector;
			else
				itsCPIs.erase( iter );
		}
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": CPIDate not in the map. Please Advise. " );
	}

	inline void AddDate( double CPIDate )
	{
		CPIMapIterator iter = itsCPIs.find(CPIDate);
		CPIRefsIterator iter2 = itsCPIRefsNb.find(CPIDate);

		if( iter2 == itsCPIRefsNb.end() )
			itsCPIRefsNb.insert( CPIRefsNbMapElt( CPIDate, 1 ) );

		iter2 = itsCPIRefsNb.find(CPIDate);

		if ( iter != itsCPIs.end() )
		{
			iter->second.first++;
			iter2->second = (( iter2->second > iter->second.first )? iter2->second : iter->second.first );
		}
		else
			itsCPIs.insert( CPIMapElt( CPIDate,CPIandItsRefNb_t(1,ARM_GP_VectorPtr(NULL)) ) ); 
	}

	inline const ARM_GP_VectorPtr GetCPI( double CPIDate ) 
	{ 
		CPIMapIterator iter = itsCPIs.find(CPIDate);
		if ( iter != itsCPIs.end() )
			return iter->second.second; 
		else
// FIXMEFRED: mig.vc8 (24/05/2007 10:46:29):cast
			return static_cast<ARM_GP_VectorPtr>(NULL);
	}

	inline const ARM_GP_VectorPtr GetCPIAndDecreaseCounter( double CPIDate ) 
	{ 
		CPIMapIterator iter = itsCPIs.find(CPIDate);
		if ( iter != itsCPIs.end() )
		{
			ARM_GP_VectorPtr returnedVector = iter->second.second;
			long newcounter = iter->second.first-1;

			if( newcounter < 1 )
				itsCPIs.erase( iter );

			return returnedVector;
		}
		else
// FIXMEFRED: mig.vc8 (24/05/2007 10:46:29):cast
			return static_cast<ARM_GP_VectorPtr>(NULL);
	}

	inline size_t getCPISize() const
	{
		return itsCPIs.size();
	}

	inline long getOldestCPICounter() const
	{
		return itsCPIs.begin()->second.first;
	}

	inline void RebuildCPIMap()
	{
		itsCPIs = CPIMap_t();
		CPIRefsIterator iter = itsCPIRefsNb.begin(), end = itsCPIRefsNb.end();

		for( ; iter != end ; iter++)
			itsCPIs.insert( CPIMapElt( iter->first, CPIandItsRefNb_t(iter->second,ARM_GP_VectorPtr(NULL)) ) ); 
	}
};


CC_END_NAMESPACE()