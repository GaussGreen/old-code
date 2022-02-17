/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file rootobject.h
 *  \brief 
 *	\author  E Benhamou, JM Prié
 *	\version 1.0
 *	\date November 2004
 */


#ifndef _INGPBASE_ROOTOBJECT_H
#define _INGPBASE_ROOTOBJECT_H

#include "port.h"
#include "env.h"

/// kernel
#include <glob/armglob.h>
#include <glob/expt.h>

#include <string>
CC_USING_NS(std,string)


CC_BEGIN_NAMESPACE( ARM )

class ARM_RootObject : public ARM_Object
{
public:
    ARM_RootObject(): ARM_Object() {}
    ARM_RootObject( const ARM_RootObject& rhs )
    : ARM_Object( rhs )
    {}
    ARM_RootObject& operator=( const ARM_RootObject& rhs )
    {
        if( this!= &rhs ) 
            ARM_Object::operator =(rhs);
        return *this;
    }
    virtual ~ARM_RootObject() {}

    /// a default for the view method
    virtual void View(char* id = NULL, FILE* ficOut = NULL) const;
	/// ugly to remove the non const version!
    virtual void View(char* id = NULL, FILE* ficOut = NULL) { ((const ARM_RootObject*) this)->View(id, ficOut ); }

	/// the two methods to reimplement
    virtual string toString(const string& indent="", const string& nextIndent="") const = 0;

	inline string GetExportShortName() const
	{ 
#if defined(__GP_STRICT_VALIDATION)
		string exportShortName = ExportShortName();
		if( exportShortName.size() != 5 )
			ARM_THROW( ERR_INVALID_ARGUMENT, " exportShortName.size() != 5" );
		if( exportShortName[0] != 'L' )
			ARM_THROW( ERR_INVALID_ARGUMENT, " exportShortName should start by 'L'" );
	
		for( size_t i=1; i<5; ++i )
			if(		( exportShortName[i] < 'A' || exportShortName[i] > 'Z' )
				&&  ( exportShortName[i] < '0' || exportShortName[i] > '9' ))
				ARM_THROW( ERR_INVALID_ARGUMENT, " exportShortName following letter can only be capital letters or numbers!" );
#endif
		return ExportShortName();
	}

	virtual string ExportShortName() const { return "LANYC";}
	
    /// the function to implement!
    virtual ARM_Object* Clone() const = 0;
    /// ugly to remove the non const clone
	virtual ARM_Object* Clone() { return ((const ARM_RootObject*) this)->Clone(); }
};


#define CC_ARM_SETNAME(name) SetName(name)


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
