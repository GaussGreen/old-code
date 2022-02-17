/*
***************************************************************************
** SOURCE FILE: cversion.c
** CREATED BY:  Peter Taylor (21 January 2000)
**
** Some trivial functions returning the version via a function interface.
**
** $Header: /home/alibcvs/alibcvs/.cvsroot-alib/platdep/src/cversion.c,v 1.1 2000/01/21 18:59:05 ptaylor Exp $
***************************************************************************
*/

#include "irx/cversion.h"

/*f
* Returns the IRX_VERSION string - useful to have as a function so that when
* debugging code which uses the Analytics Library you can print the output
* from the this function within the debugger. Note that you cannot print the
* macro IRX_VERSION within a debugger.
*/
const char* IrxVersion (void)
{
    return IRX_VERSION;
}

/*f
* Returns the IRX_VERSION_NUM string - useful to have as a function so that when
* debugging code which uses the Analytics Library you can print the output
* from the this function within the debugger. Note that you cannot print the
* macro IRX_VERSION_NUM within a debugger.
*/
long IrxVersionNum (void)
{
    return IRX_VERSION_NUM;
}

