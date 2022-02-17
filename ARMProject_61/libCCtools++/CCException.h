/* -------------------------------------------------------------------------

   File: CCException.h
   Path: /home/julia/projets/dev/Common/libCCtools++/SCCS/s.CCException.h
   Description: header CCException 
   Created: 98/06/02
   Author: Charles-Emmanuel MUSY 
   Modified: 99/01/19 17:25:29
   Last maintained by: Charles-Emmanuel MUSY
   Revision: 1.3

   -------------------------------------------------------------------------

   Note:

   ------------------------------------------------------------------------- */

#ifndef CCEXCEPTION_H
#define CCEXCEPTION_H

#include <CCcommon.h>
SCCS_ID (CCException_h_SccsId, "@(#)CCException.h	1.3, modified 99/01/19");

#include "CCString.h"
#include <iosfwd>

using std::ostream;

class CCException
{
friend ostream& operator<< (ostream& os, CCException& e);
public:
	CCException ();
	CCException (const char *n_message);
	CCException (const CCString& n_message);

	CCString& GetMessage ();
	// operator char* ();
	operator const char* ();
		
private:
	CCString message;
};

#ifdef CCException_c

#endif // CCException_c

#endif // CCEXCEPTION_H

// EOF CCException.h
