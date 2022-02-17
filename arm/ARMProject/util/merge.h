/*
 * $Log: merge.h,v $
 *
 * 2006/06/19 : JLA : Adding const interface
 *
 * Revision 1.3  2003/07/11 12:17:04  emezzine
 *  Add MergeDatesForMc
 *
 * Revision 1.2  2002/11/25 14:22:34  mab
 * Formatting
 *
 */


/*----------------------------------------------------------------------------*

 	merge.h
 
*----------------------------------------------------------------------------*/
#ifndef _MERGE_H
#define _MERGE_H





#include "linalg.h"




extern void MergeSched(ARM_Vector** DateRes, ARM_Vector** CFRes,
                       ARM_Vector* Date1, ARM_Vector* CF1,
                       ARM_Vector* Date2, ARM_Vector* CF2);

extern void MergeDates(ARM_Vector** DateRes,
                       const ARM_Vector* Date1,
                       const ARM_Vector* Date2);

extern void MergeDates(ARM_Vector** DateRes,
                       const ARM_Vector& Date1,
                       const ARM_Vector& Date2);

ARM_Vector* MergeDatesForMc(ARM_Vector* Date1,ARM_Vector* Date2);


#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
