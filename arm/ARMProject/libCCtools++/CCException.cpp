/* -------------------------------------------------------------------------

   File: CCException.cc
   Path: /home/julia/projets/dev/Common/libCCtools++/SCCS/s.CCException.cc
   Description: implantation de la classe CCException 
   Created: 98/06/02
   Author: Charles-Emmanuel MUSY 
   Modified: 99/01/19 17:25:28
   Last maintained by: Charles-Emmanuel MUSY
   Revision: 1.3

   -------------------------------------------------------------------------

   Note:

   ------------------------------------------------------------------------- */

#include <CCcommon.h>
SCCS_ID (CCException_c_SccsId, "@(#)CCException.cc	1.3, modified 99/01/19");

#define CCException_c
#include "CCException.h"

CCException::CCException ()
{
}

CCException::CCException (const char *n_message) : message (n_message)
{
}

CCException::CCException (const CCString& n_message) : message (n_message)
{
}

CCString& CCException::GetMessage ()
{
	return message;
}

CCException::operator const char* ()
{
	return message;
}

// CCException::operator char* ()
// {
// 	return message;
// }

ostream& operator<< (ostream& os, CCException& e)
{
        return os << e.message;
}

// EOF CCException.cc
