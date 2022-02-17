/***************************************************************
 * Module:	PenGuin
 * Submodule:	
 * File:	kstdinc.h
 * Function:	Standard Definition and Include files.
 * Author:	Christian Daher
 * Revision:	$Header$
 ***************************************************************/
#ifndef	_kstdinc_H
#define	_kstdinc_H

// CMLIB 
#include "General/General.h"
using CM::Object;
using CM::SharedPointer;
using CM::ConstSharedPointer;
using CM::Raw2SharedPointer;
using CM::SharedPointerConvertTo;
using CM::SharedPointerDownCastTo;
using CM::String;
using CM::Date;
using CM::Array;

#include "kexcept.h"
#include "kplatdep.h"
#include "kerrlog.h"		// Logging and error messages


#include <stdio.h>
#include <float.h>
#include <math.h>
#include <limits.h>		/* MIN/MAX */
#include <assert.h>


extern	"C" {
#include "drlstd.h"

#include "ldate.h"
#include "convert.h"
#include "cmemory.h"
#include "macros.h"
#include "bastypes.h"

#include "drlio.h"		/* stdout */
#include "drlproc.h"		/* GlobFlagGet */
};




/*--------------------------------------------------------------
 *
 */


extern	void	DppExecCmd(const char *fmt,...);


typedef KMap(String, double)::value_type String_Double;
typedef KMap(TDate, double)::value_type TDate_Double;
typedef KMap(TDate, TDate)::value_type TDate_TDate;
//typedef KMap(TDate, KVector(double))::value_type TDate_KVectorDouble;


#endif


