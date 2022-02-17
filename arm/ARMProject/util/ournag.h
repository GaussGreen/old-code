/*
 * $Log: ournag.h,v $
 * Revision 1.2  2003/08/15 06:55:40  ebenhamou
 * default is now to ignore nwag maxmin
 * use the _SHOW_MAXMIN_WARNING to see warning
 *
 * Revision 1.1  2003/06/23 11:58:23  ebenhamou
 * Initial revision
 *
 *
 */


#ifndef _OURNAG_H
#define _OURNAG_H

/// This is just to remove the 
/// redefinition by nag of various MAX MIN and so on

#ifndef _SHOW_MAXMIN_WARNING
	#ifdef MAX
		#undef MAX
	#endif

	#ifdef MIN
		#undef MIN
	#endif

	#ifdef ROUND
		#undef ROUND
	#endif
#endif

#include "nag.h"
			

#endif
