/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file InfHybModel.h
 *  \brief file for the inflation hybrid model
 *	\author  E Benhamou
 *	\version 1.0
 *	\date July 2004
 */


#ifndef _INGPINFLATION_INFHYBMODEL_H
#define _INGPINFLATION_INFHYBMODEL_H
/*
#pragma warning(disable : 4786) 

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "mod/model.h"
#include "gpbase/typedef.h"
#include "gpbase\removeidentifiedwarning.h"

#include <map>

CC_USING_NS( std, map )

CC_BEGIN_NAMESPACE( ARM )

class ARM_InfHybridModel: public ARM_Model{

public:

	ARM_InfHybridModel():ARM_Model(){}
	ARM_InfHybridModel( const map< string,ARM_Object*> &);
	ARM_InfHybridModel( const ARM_InfHybridModel &  );

	ASSIGN_OPERATOR	( ARM_InfHybridModel );
	virtual ARM_Object*	Clone	( ) const	{ return new ARM_InfHybridModel (*this);	}
	virtual ~ARM_InfHybridModel( )	{ 
	};

	virtual void CheckMkt();

	virtual ARM_ObjectPtr Get( const string & ) const;

private:
	map< string, ARM_ObjectPtr>		 itsMktMap;
};

CC_END_NAMESPACE()
*/
#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
