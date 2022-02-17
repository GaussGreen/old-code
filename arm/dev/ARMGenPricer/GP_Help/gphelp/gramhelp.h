/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gramhelp.h
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPHELP_GRAMHELP_H
#define _INGPHELP_GRAMHELP_H

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/port.h"
#include "gpbase/singleton.h"
#include "gpbase/env.h"
#include "gpbase/rootobject.h"

#include <string>
CC_USING_NS( std, string )

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
struct ARM_GramFunctionStruct;

/// \class ARM_GramHelper 
/// very light class to give some help 
///	 on the generic  language!
/// if the FuncName is empty 
/// returns the general help
/// otherwise gives specific
/// help on the function!

struct ARM_GramHelper : ARM_RootObject
{
	void FunctionTableHelpString( ARM_GramFunctionStruct* funcTable, size_t funcTableSize  );
	ARM_GramHelper( const string& FuncName ="ALLFUNC" );

	/// standard ARM Object support
	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="", const string& nextIndent="") const { return itsHelpText; }
private:
	string itsHelpText;
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

