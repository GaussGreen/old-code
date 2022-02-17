/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 */

/*! \file cstmanager.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou, R. Guillemot
 *	\version 1.0
 *	\date October 2004
 */

#include "gpinfra/cstmanager.h"
#include "gpbase/checkarg.h"
#include "gpinfra/gramfunction.h"
#include "gpinfra/cstmanager.h"
#include "gpbase/ostringstream.h"
#include "gpbase/gpmatrix.h"

CC_USING_NS(std,pair)

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_CstManager
///	Routine: Constructor
///	Returns: 
///	Action : builds the object
////////////////////////////////////////////////////
	
ARM_CstManager::ARM_CstManager( const vector<string>& names, const CC_NS(std,vector)<double>& csts )
: ARM_RootObject()
{
	CC_NS(ARM_Check,CheckSameArgSize)( names, csts, "names", "csts" );

	for( size_t i=0 ; i<names.size(); ++i )
	{
		if( ARM_GramFunction::IsAFunctionName( names[i] ) )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME + ": name "  + names[i]  + " already taken for function name! " ); 
		double d0 = csts[i];
		ARM_GramFctorArg* d = new ARM_GramFctorArg(d0);
		itsConstMap.insert( pair<const string,ARM_GramFctorArg*>( names[i], d) );
	}
}
////////////////////////////////////////////////////
///	Class  : ARM_CstManager
///	Routine: Constructor
///	Returns: 
///	Action : builds the object unsing GramFunctorArgs
////////////////////////////////////////////////////
ARM_CstManager::ARM_CstManager( const vector<string>& names , const vector <ARM_GramFctorArg>& csts)
: ARM_RootObject()
{
	if(csts.size()!=names.size())
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : The size of namesVector is different from the size of Objects");
	for( size_t i=0 ; i<names.size(); ++i )
	{
		if( ARM_GramFunction::IsAFunctionName( names[i] ) )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME + ": name "  + names[i]  + " already taken for function name! " ); 
	  	itsConstMap.insert( pair<const string,ARM_GramFctorArg*>( names[i], new ARM_GramFctorArg(csts[i]) ));
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_CstManager
///	Routine: Copy Constructor
///	Returns: 
///	Action : builds the object
////////////////////////////////////////////////////
ARM_CstManager::ARM_CstManager( const ARM_CstManager& rhs )
: ARM_RootObject(rhs), itsConstMap( rhs.itsConstMap )
{
	const_iterator iter;
	for(iter = rhs.itsConstMap.begin();iter!=rhs.itsConstMap.end();++iter)
		itsConstMap[iter->first] = new ARM_GramFctorArg((*iter->second));
}


////////////////////////////////////////////////////
///	Class  : ARM_CstManager
///	Routine: Assignment operator
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_CstManager& ARM_CstManager::operator=( const ARM_CstManager& rhs )
{
	if( this!= &rhs )
	{
		ARM_RootObject::operator=( rhs );
		itsConstMap = rhs.itsConstMap;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_CstManager
///	Routine: Destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_CstManager::~ARM_CstManager()
{
	const_iterator iter;
	for(iter = itsConstMap.begin();iter!=itsConstMap.end();++iter)
		delete iter->second;

}

////////////////////////////////////////////////////
///	Class  : ARM_CstManager
///	Routine: DoesCstNameExist
///	Returns: bool
///	Action : tells whether there is a cst name!
////////////////////////////////////////////////////
bool ARM_CstManager::DoesCstNameExist( const string& name ) const
{
	stringDoubleMap::const_iterator iter = itsConstMap.find( name );
	return iter != itsConstMap .end();
}

////////////////////////////////////////////////////
///	Class  : ARM_CstManager
///	Routine: DoesCstNameExist
///	Returns: bool
///	Action : tells whether there is a cst name!
////////////////////////////////////////////////////
void ARM_CstManager::insert( const string& name ,const ARM_GramFctorArg& cst )
{
	stringDoubleMap::iterator pos = itsConstMap.find( name );	
	if(pos != itsConstMap .end())
		pos->second = new ARM_GramFctorArg(cst);
	else
		itsConstMap.insert( pair<const string,ARM_GramFctorArg*>( name, new ARM_GramFctorArg(cst)));
}

////////////////////////////////////////////////////
///	Class  : ARM_CstManager
///	Routine: Clone
///	Returns: ARM_Object
///	Action : clones the object
////////////////////////////////////////////////////
ARM_Object* ARM_CstManager::Clone() const
{
	return new ARM_CstManager(*this);
}


////////////////////////////////////////////////////
///	Class  : ARM_CstManager
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_CstManager::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << "ARM_CstManager\n";
	
	if( ! itsConstMap.empty() )
	{
		stringDoubleMap::const_iterator iter	= itsConstMap.begin();
		stringDoubleMap::const_iterator iterEnd	= itsConstMap.end();

		os << "Names " << "\t Values \n";
		for( ; iter != iterEnd; ++iter )
		{
			os << (*iter).first << "\t" << iter->second->toString() << "\n";
		}
	}
	return os.str();
}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/


