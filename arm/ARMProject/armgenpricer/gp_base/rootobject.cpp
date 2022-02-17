/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file rootobject.cpp
 *  \brief 
 *	\author  E. Benhamou, JM Prié
 *	\version 1.0
 *	\date November 2004
 */

#include "gpbase/rootobject.h"

CC_BEGIN_NAMESPACE( ARM )

/////////////////////////////////////////////////////////////////
///	Class  : ARM_RootObject 
///	Routine: View
///	Returns: void
///	Action :
/////////////////////////////////////////////////////////////////
void ARM_RootObject::View(char* id, FILE* ficOut) const
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
	string allErrors = toString();
    fprintf(fOut, "%s", allErrors.c_str() );

	/// to allow to have nested view
    if ( NULL == ficOut )
       fclose(fOut);
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

