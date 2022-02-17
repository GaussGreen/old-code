/*
 * $Log: removeournagwarning.h,v $
 * Revision 1.1  2003/06/23 11:58:23  ebenhamou
 * Initial revision
 *
 */


#ifndef _INGPBASE_REMOVEOURNAGWARNING_H
#define _INGPBASE_REMOVEOURNAGWARNING_H

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


#endif 

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
