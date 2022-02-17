 /***************************************************************************
  *	SCCS Keyword Information
  *	------------------------
  *	Module name	:  extern.h
  *	Version		:  2.4
  *	Extracted	:  2/2/94 at 03:24:52
  *	Last Updated	:  12/10/93 at 16:39:32
  ***************************************************************************/

#ifndef EXTERNH

#define EXTERNH

#ifndef extern_h_SCCSINFO
#define extern_h_SCCSINFO
#endif


#ifdef __cplusplus
#define EXTERNC extern "C"
#define ICMO_ExtC extern "C"
#else
#define EXTERNC extern
#define ICMO_ExtC extern
#endif

#endif
