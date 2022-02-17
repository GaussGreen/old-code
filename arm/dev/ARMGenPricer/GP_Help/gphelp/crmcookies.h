/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file crmcookies.h
 *
 *  \brief customer relationship management cookies
 *	\author  E Benhamou
 *	\version 1.0
 *	\date FEBRUARY 2004
 */

#ifndef _INGPHELP_CRMCOOKIES_H
#define _INGPHELP_CRMCOOKIES_H

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/port.h"

#include "gpbase/singleton.h"

/// STL
#include <string>
CC_USING_NS(std,string)
#include <set>
CC_USING_NS(std,set)

#include <ctime> /// for clock_t

CC_BEGIN_NAMESPACE( ARM )

/// \struct ARM_CRMCookiesImp
/// is the implementor of the Customer Relationship Management tracer
///	his goal is to monitor for a given user the time spent on a given application
///	and the services used
/// At this stage, this is a very basic CRM cookies

struct ARM_CRMCookiesImp
{
	void RegisterService( const string& userName, const string& serviceName, const string& folderName = "" );
	void EndTracing();

private:
	clock_t itsStartTime;
	string itsFileName;
	set<string> itsServices;

	void InitTracing( const string& userName, const string& folderName );

	/// to forbid client from using it except for the singleton holder
	ARM_CRMCookiesImp();
	friend class ARM_SingletonHolder<ARM_CRMCookiesImp>;
};


extern ARM_SingletonHolder<ARM_CRMCookiesImp> ARM_CRMCookies;


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
