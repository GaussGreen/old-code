/*---------------------------------------------------------------------- 
HEADER FILE:    CGeneral.h

CREATED BY:     4/29/92 Ross Boylan

PURPOSE:        General purpose header for ANSI C code.

$Header: /home/alibcvs/alibcvs/.cvsroot-alib/platdep/include/cgeneral.h,v 1.59 2003/10/07 13:33:21 cbauer Exp $
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
#ifndef IRX_CGENERAL_H
#define IRX_CGENERAL_H

#include <stdlib.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C"
{
#endif 

typedef int IrxTBool;

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef SUCCESS
#define SUCCESS 0
#endif

#ifndef FAILURE
#define FAILURE -1
#endif

#ifdef __cplusplus
}
#endif

#endif    /* IRX_CGENERAL_H */
