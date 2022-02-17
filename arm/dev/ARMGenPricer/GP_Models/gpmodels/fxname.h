/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 * 
 * \file forexvanilla.h
 *  \brief to create an object for Fx Name
 *
 *	\author  K. Belkheir 
 *  \ reviwer E.Ezzine
 *	\version 1.0
 *	\date January 2007
 */

#ifndef _INGPMODELS_FXNAME_H
#define _INGPMODELS_FXNAME_H

#include "gpbase/assignop.h"
#include "gpbase/rootobject.h"
#include "gpbase/gplinalgtypedef.h"
#include "gpbase/gpvector.h"

CC_BEGIN_NAMESPACE( ARM )

/// ------------------------------
///	--- Forex Base Class ----
/// ------------------------------
class ARM_FXName : public ARM_RootObject
{
	static const string MktNamesTable [];
	static const string InvMktNamesTable [];

private:
	bool itsIsInvMkt;
	string itsMktName;
	ARM_GP_StrVector itsMktNamesVect;
	ARM_GP_StrVector itsInvMktNamesVect;
public: 
	/// standard constructor
	ARM_FXName(const string& mktName);
	//copy constructor
	ARM_FXName( const ARM_FXName& rhs );
	ASSIGN_OPERATOR(ARM_FXName)
	//destructor
	virtual ~ARM_FXName() {}
	
	/// what it is for ...
	bool IsInvMkt( const string& fxName );

	//accesors
	bool GetIsInvMkt() const {return itsIsInvMkt;};
	string GetMktName() const {return itsMktName;};

    /// Standard ARM object support
	virtual ARM_Object* Clone()  const { return new ARM_FXName(*this); };
	virtual string toString(const string& indent="",const string& nextIndent="") const { return string("Forex name ");};
	virtual string ExportShortName() const { return "LFX";}
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/



