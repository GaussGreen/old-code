/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: errviewer.h,v $
 * Revision 1.1  2004/03/01 07:52:17  ebenhamou
 * Initial revision
 *
 */


/*! \file errviewer.h
 *
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */


#ifndef _INGPBASE_ERRVIEWER_H
#define _INGPBASE_ERRVIEWER_H

#include "port.h"
#include "rootobject.h"

#include <string>
CC_USING_NS(std,string)

CC_BEGIN_NAMESPACE( ARM )

struct ARM_ErrViewer : public ARM_RootObject
{
	static const string LogOfAllErrorsFileNameStart;
	static const string LogOfAllErrorsFileNameExt;

	/// constructor,copy constructor,assignment operator and destructor
	ARM_ErrViewer( bool reset = false );
	ARM_ErrViewer( const ARM_ErrViewer& rhs );
	ARM_ErrViewer& operator=( const ARM_ErrViewer& rhs );
	virtual ~ARM_ErrViewer();

	/// standard ARM Object support
	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="",const string& nextIndent="") const;

private:
	void ReadTextFromFile( const string& fileInName, FILE* fOut, const char* contentName ) const;
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

