/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *	\file cstmanager.h
 *
 *  \brief 
 *	\author  E. Benhamou, R. Guillemot
 *	\version 1.0
 *	\date October 2004
 */

#ifndef _INGPINFRA_CSTMANAGER_H
#define _INGPINFRA_CSTMANAGER_H

/// this headers has to come first
/// as it contains pragma for stupid warning
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/rootobject.h"
#include "gramfunctorarg.h"
#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gpbase/rootobject.h"

/// STL
#include <map>
CC_USING_NS(std,map)
#include <string>
CC_USING_NS(std,string)
#include <vector>
CC_USING_NS(std,vector)
CC_USING_NS(std,less)


CC_BEGIN_NAMESPACE( ARM )

class ARM_CstManager : public ARM_RootObject
{
private:
	typedef map<string, ARM_GramFctorArg*, less<string> > stringDoubleMap;
	stringDoubleMap itsConstMap;

public:
	ARM_CstManager( const vector<string>& names = CC_NS(std,vector)<string>(), const CC_NS(std,vector)<double>& cst = CC_NS(std,vector)<double>() );
	ARM_CstManager( const vector<string>& names , const vector <ARM_GramFctorArg>& cst );
	ARM_CstManager( const ARM_CstManager& rhs );
	ARM_CstManager& operator=( const ARM_CstManager& rhs );
	virtual ~ARM_CstManager();
	bool DoesCstNameExist( const string& name ) const;

	/// iterator support
	typedef stringDoubleMap::iterator iterator;
	typedef stringDoubleMap::const_iterator const_iterator;

	/// iterator accessor
	iterator begin() { return itsConstMap.begin(); }
	iterator end() { return itsConstMap.end(); }
	const_iterator begin() const { return itsConstMap.begin(); }
	const_iterator end() const { return itsConstMap.end(); }
	iterator find( const string& name ) { return itsConstMap.find( name ); }
	void insert( const string& name,const ARM_GramFctorArg& cst );
	const_iterator find( const string& name ) const { return itsConstMap.find( name ); }

	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="", const string& nextIndent="") const;
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

