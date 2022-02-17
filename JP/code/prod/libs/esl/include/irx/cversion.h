/*-----------------------------------------------------------------------
HEADER FILE:  cversion.h

CREATED BY:   Jon Freeman, 1/4/95

PURPOSE:      Version definition.

$Header: /home/alibcvs/alibcvs/.cvsroot-alib/platdep/include/cversion.h,v 1.83 2005/05/26 14:21:32 cbauer Exp $
---------------------------------------------------------------------- 
Proprietary software, whether developed for Morgan by in-house
staff or external contractors, must not be copied or removed from Morgan
premises without the approval of a Senior Vice President and Audit.

This proprietary software has been developed strictly for Morgan's own
internal use.  Any use or misuse, intentional or otherwise which contradicts
or places this policy in jeopardy is strictly forbidden.  Also, any actions or 
inactions resulting in a breach of this goal will be dealt with harshly.

Do Not Distribute this to anyone outside the Analytics Library
Group without authorization from its management. 

Copyright 1995 J.P. Morgan & Co. Incorporated.   All rights reserved.
-------------------------------------------------------------------------  */
#ifndef IRX_CVERSION_H
#define IRX_CVERSION_H

#include "cgeneral.h"

#ifdef __cplusplus
extern "C"
{
#endif

#if 0
#define IRX_VERSION "Analytics Library Version 13.6"
#else
#define IRX_VERSION "Analytics Library Development Version (based on 13.6)"
#endif
#define IRX_VERSION_NUM 1360

/*f
* Returns the IRX_VERSION string - useful to have as a function so that when
* debugging code which uses the Analytics Library you can print the output
* from the this function within the debugger. Note that you cannot print the
* macro IRX_VERSION within a debugger.
*/
char* irxVersion (void);

/*f
* Returns the IRX_VERSION_NUM string - useful to have as a function so that when
* debugging code which uses the Analytics Library you can print the output
* from the this function within the debugger. Note that you cannot print the
* macro IRX_VERSION_NUM within a debugger.
*/
long irxVersionNum (void);

#ifdef __cplusplus
}
#endif

#endif    /* IRX_CVERSION_H */
