
/*
 * $Log: merge.cpp,v $
 *
 * 2006/06/19 : JLA : Adding const interface
 *
 * Revision 1.7  2003/07/11 12:16:28  emezzine
 * Add MergeDatesForMc
 *
 * Revision 1.6  2002/11/25 14:22:54  mab
 * Formatting
 *
 * Revision 1.5  2001/03/21 15:52:43  abizid
 * Acceleration en passant aux GetElt()
 *
 * Revision 1.4  1999/04/14 15:55:09  mab
 * void MergeDates(ARM_Vector** DateRes, ARM_Vector* Date1, ARM_Vector* Date2)
 * Liberation de tmpdate
 *
 * Revision 1.3  1998/11/13 14:16:23  nicolasm
 * Ajout de l'entete
 *
 */

/*----------------------------------------------------------------------------*
 
    merge.cpp
 
    This file implements merging vectors functions 

*----------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "merge.h"





/*----------------------------------------------------------------------------*/
/*     Merge two schedules (Date1,CF1) and (Date2,CF2) in (DateRes,CFRes)     */
/*----------------------------------------------------------------------------*/

void MergeSched(ARM_Vector** DateRes, ARM_Vector** CFRes,
                ARM_Vector* Date1, ARM_Vector* CF1,
                ARM_Vector* Date2, ARM_Vector* CF2)
{
    int i1 = 0;
    int i2 = 0;
    int ires = 0;

    int size1 = Date1->GetSize();
    int size2 = Date2->GetSize();

    double* tmpdate = new double[size1+size2];
    double* tmpcf   = new double[size1+size2];



    while (( i1 < size1 ) && ( i2 < size2 ))
    {
        if ( Date1->Elt(i1) == Date2->Elt(i2) )
        {
           tmpdate[ires] = Date1->Elt(i1);
           tmpcf[ires]   = CF1->Elt(i1) + CF2->Elt(i2);

           ires++;
           i1++;
           i2++;
        }
        else if ( Date1->Elt(i1) < Date2->Elt(i2) )
        {
           tmpdate[ires] = Date1->Elt(i1);
           tmpcf[ires]  = CF1->Elt(i1);

           ires++;
           i1++;
        }
        else 
        {
           tmpdate[ires] = Date2->Elt(i2);
           tmpcf[ires]  = CF2->Elt(i2);

           ires++;
           i2++;
        }
    }

    while ( i1 < size1 )
    {
        tmpdate[ires] =  Date1->Elt(i1);
        tmpcf[ires]  = CF1->Elt(i1);

        ires++;
        i1++;
    }

    while ( i2 < size2 )
    {
        tmpdate[ires] = Date2->Elt(i2);
        tmpcf[ires]  = CF2->Elt(i2);

        ires++;
        i2++;
    }

    *DateRes = new ARM_Vector(ires, tmpdate);
    *CFRes   = new ARM_Vector(ires, tmpcf);

    if (tmpdate)
       delete tmpdate;

    if (tmpcf)
       delete tmpcf;
}



void MergeDates(ARM_Vector** DateRes,
				const ARM_Vector& Date1,
				const ARM_Vector& Date2)
{
     int i1 = 0;
     int i2 = 0;
     int ires = 0;

     int size1 = Date1.GetSize();
     int size2 = Date2.GetSize();

     double* tmpdate  = new double[size1+size2];
     const double* tmpdate1 = Date1.GetElt();
     const double* tmpdate2 = Date2.GetElt();


     while (( i1 < size1 ) && ( i2 < size2 ))
     {
         if ( tmpdate1[i1] == tmpdate2[i2] )
         {
            tmpdate[ires] = tmpdate1[i1];

            ires++;
            i1++;
            i2++;
         }
         else if ( tmpdate1[i1] < tmpdate2[i2] )
         {
            tmpdate[ires] = tmpdate1[i1];

            ires++;
            i1++;
         }
         else 
         {
            tmpdate[ires] = tmpdate2[i2];

            ires++;
            i2++;
         }
     }

     while ( i1 < size1 )
     {
         tmpdate[ires] = tmpdate1[i1];

         ires++;
         i1++;
     }

     while ( i2 < size2 )
     {
         tmpdate[ires] = tmpdate2[i2];

         ires++;
         i2++;
     }

     *DateRes = new ARM_Vector(ires, tmpdate);


     if (tmpdate)
        delete [] tmpdate;
}
void MergeDates(ARM_Vector** DateRes,
                const ARM_Vector* Date1,
                const ARM_Vector* Date2)
{

	MergeDates(DateRes,*Date1,*Date2); 
}


ARM_Vector* MergeDatesForMc(ARM_Vector* Date1,ARM_Vector* Date2)
{
     int i1 = 0;
     int i2 = 0;
     int ires = 0;

     int size1 = Date1->GetSize();
     int size2 = Date2->GetSize();

     double* tmpdate  = new double[size1+size2];
     double* tmpdate1 = Date1->GetElt();
     double* tmpdate2 = Date2->GetElt();


     while (( i1 < size1 ) && ( i2 < size2 ))
     {
         if ( fabs(tmpdate1[i1] - tmpdate2[i2]) <= FRMVOL_LAG_THRESHOLD)
         {
            tmpdate[ires] = tmpdate1[i1];

            ires++;
            i1++;
            i2++;
         }
         else if ( tmpdate1[i1] < tmpdate2[i2] + FRMVOL_LAG_THRESHOLD)
         {
            tmpdate[ires] = tmpdate1[i1];

            ires++;
            i1++;
         }
         else 
         {
            tmpdate[ires] = tmpdate2[i2];

            ires++;
            i2++;
         }
     }

     while ( i1 < size1 )
     {
         tmpdate[ires] = tmpdate1[i1];

         ires++;
         i1++;
     }

     while ( i2 < size2 )
     {
         tmpdate[ires] = tmpdate2[i2];

         ires++;
         i2++;
     }

     ARM_Vector* DateRes = new ARM_Vector(ires, tmpdate); 

     if (tmpdate)
        delete tmpdate;

     return DateRes;
}




/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
