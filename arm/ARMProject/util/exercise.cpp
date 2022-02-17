/*
 * $Log: exercise.cpp,v $
 * Revision 1.10  2004/06/08 08:10:13  mab
 * Formatting
 *
 * Revision 1.9  2004/01/21 12:04:44  emezzine
 * added #include.
 *
 * Revision 1.8  2001/07/30 09:17:22  smysona
 * correction leak (cleanit(exerciseDates);)
 *
 * Revision 1.7  2001/04/03 12:03:38  nicolasm
 * Construction des exrecices a partir d'une fixleg
 *
 * Revision 1.6  2001/03/16 10:08:52  smysona
 * Correction du isexercise pour les plages americaines
 *
 * Revision 1.5  2001/03/12 19:27:10  smysona
 * Acceleration du isexercise
 *
 * Revision 1.4  2001/02/22 19:43:08  smysona
 * IsExerciseDate travaille sur julainDate
 *
 * Revision 1.3  2000/11/06 09:29:53  sgasquet
 * Ajout include swapleg.h
 *
 * Revision 1.2  1999/04/06 09:58:09  ypilchen
 * Rajout de Log pour RCS et Amelioration des Set(..)
 *
 */


/*----------------------------------------------------------------------------*
    exercise.cpp
 
  Body for the ARM_ExerciseStyle class, a class to deal with option exercice
*----------------------------------------------------------------------------*/



 
#include "exercise.h"
#include "interpol.h"
//#include "swapleg.h"



/*----------------------------------------------------------------------------*
  Default constructor : 1 Year European option from Today
*----------------------------------------------------------------------------*/ 


ARM_ExerciseStyle::ARM_ExerciseStyle(void)
{
    ARM_Date today;


    Init();

    today.AddYears(1);
    
    itsExerciseType = K_EUROPEAN;

    itsExerciseDatesSize = 1;

    itsExerciseStartDates = new ARM_Vector(1, today.GetJulian());

    itsExerciseEndDates = new ARM_Vector(1, today.GetJulian());
}



/*----------------------------------------------------------------------------*
    Constructor for European (exercise at expiry Date)
*----------------------------------------------------------------------------*/ 
void ARM_ExerciseStyle::Set(ARM_Date& xDate)
{
    FreeStyle();

    itsExerciseType = K_EUROPEAN;
    
    itsExerciseDatesSize = 1;

    if (itsExerciseStartDates)
    {
       delete itsExerciseStartDates;
       itsExerciseStartDates = NULL;
    }

    itsExerciseStartDates = new ARM_Vector(1, xDate.GetJulian());

    if (itsExerciseEndDates)
    {
       delete itsExerciseEndDates;
       itsExerciseEndDates = NULL;
    }

    itsExerciseEndDates = new ARM_Vector(1, xDate.GetJulian());
}



ARM_ExerciseStyle::ARM_ExerciseStyle(ARM_Date& xDate)
{
    Init();
    
    Set(xDate);
}



/*----------------------------------------------------------------------------*
    Constructor for American Option 
*----------------------------------------------------------------------------*/ 
void ARM_ExerciseStyle::Set(ARM_Date& xStartDate, ARM_Date& xEndDate)
{
    FreeStyle();

    itsExerciseType = K_AMERICAN;
    
    itsExerciseDatesSize = 1;
    
    if (itsExerciseStartDates)
    {
       delete itsExerciseStartDates;
       itsExerciseStartDates = NULL;
    }

    itsExerciseStartDates = new ARM_Vector(1, xStartDate.GetJulian());


    if (itsExerciseEndDates)
    {
       delete itsExerciseEndDates;
       itsExerciseEndDates = NULL;
    }

    itsExerciseEndDates = new ARM_Vector(1, xEndDate.GetJulian());
}



ARM_ExerciseStyle::ARM_ExerciseStyle(ARM_Date& xStartDate, ARM_Date& xEndDate)
{
    Init();

    Set(xStartDate, xEndDate);
}



/*----------------------------------------------------------------------------*
    Constructor for Bermudan Option 
*----------------------------------------------------------------------------*/ 
void ARM_ExerciseStyle::Set(ARM_Vector* xDates)
{
    itsExerciseType = K_BERMUDAN;
    
    itsExerciseDatesSize = xDates->GetSize();

    SetExerciseStartDates(xDates);

    SetExerciseEndDates(xDates);
}



ARM_ExerciseStyle::ARM_ExerciseStyle(ARM_Vector* xDates)
{
    Init();

    Set(xDates);
}



void ARM_ExerciseStyle::Set(ARM_Date& xStartDate, 
                            ARM_Date& xEndDate, int liborType)
{
    /*FreeStyle();

    Init();

    itsExerciseType = K_BERMUDAN;
    
    ARM_SwapLeg* liborLeg = new ARM_SwapLeg(xStartDate, xEndDate, 
                                            liborType, K_RCV, 0.0);   

    ARM_Vector* xDates = (ARM_Vector *) liborLeg->GetResetDates()->Clone();

    itsExerciseDatesSize = xDates->GetSize();

    SetExerciseStartDates(xDates);

    SetExerciseEndDates(xDates);

    delete xDates;

    delete liborLeg;*/
}



/*----------------------------------------------------------------------------*
   Constructor for Bermudan Style with Libor type reset dates
*-----------------------------------------------------------------------------*/ 
ARM_ExerciseStyle::ARM_ExerciseStyle(ARM_Date& xStartDate, 
                                     ARM_Date& xEndDate, int liborType)
{
    Init();

    Set(xStartDate, xEndDate, liborType);
}



/*----------------------------------------------------------------------------*
    Constructor for Customized Exercise Dates Option 
*-----------------------------------------------------------------------------*/

void ARM_ExerciseStyle::Set(ARM_Vector* xStartDates, ARM_Vector* xEndDates)
{
    FreeStyle();

    if ( xStartDates->GetSize() != xEndDates->GetSize() )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
      "Start Dates and End Dates vectors must have the same size");
    }
    
    itsExerciseType = K_CUSTOMIZED;
    
    itsExerciseDatesSize = xEndDates->GetSize();

    SetExerciseStartDates(xStartDates);

    SetExerciseEndDates(xEndDates);
}



ARM_ExerciseStyle::ARM_ExerciseStyle(ARM_Vector* xStartDates, 
                                     ARM_Vector* xEndDates)
{
    Init();

    Set(xStartDates, xEndDates);
}



/*----------------------------------------------------------------------------*
Constructor for Bermudan Style from swap leg or derived class (eg cap) dates
*----------------------------------------------------------------------------*/
 
ARM_ExerciseStyle::ARM_ExerciseStyle(ARM_SwapLeg* swapLeg, int dateType, 
                                     ARM_Date& firstDate, ARM_Date& lastDate)
{
    Init();
 
    Set(swapLeg, dateType, firstDate, lastDate);
}



void ARM_ExerciseStyle::Set(ARM_SwapLeg* swapLeg, int dateType, 
                            ARM_Date& firstDate, ARM_Date& lastDate)
{
    //FreeStyle();
 
    //Init();
 
    //itsExerciseType = K_BERMUDAN;
 
    //ARM_Vector* xDates;
 
    //// PENDING 
    //// The reset dates, it they have a sense for a fixed leg
    //// should be computed at construction at startFlow - spot
    //// of the currency.
    //// To avoid border effect, we only put it here,
    //// even if it is dirty

    //if (swapLeg->IsFixedLeg())
    //{
    //   xDates = swapLeg->GetFlowStartDates();
    //}
    //else
    //{
    //    switch(dateType)
    //    {
    //        case K_INDEX_RESET:
    //        {
    //            xDates = swapLeg->GetResetDates();
    //        };
    //        break;
 
    //        case K_FLOW_PAYMENT:
    //        {
    //            xDates = swapLeg->GetPaymentDates();
    //        };
    //        break;
 
    //        case K_START_INTEREST:
    //        {
    //            xDates = swapLeg->GetFlowStartDates();
    //        };
    //        break;
 
    //        case K_END_INTEREST:
    //        {
    //            xDates = swapLeg->GetFlowEndDates();
    //        };
    //        break;

    //        default:
    //        {
    //            xDates = swapLeg->GetResetDates();
    //        };
    //    }
    //}

    //int index1 = indexAfterValue(xDates, firstDate.GetJulian());

    //int index2 = indexBeforeValue(xDates, lastDate.GetJulian());

    //ARM_Vector* exerciseDates = new ARM_Vector(xDates, index1, index2); 
 
    //if (swapLeg->IsFixedLeg())
    //{
    //   int i, size = exerciseDates->GetSize();
    //   int   spot  = swapLeg->GetCurrencyUnit()->GetSpotDays();
    //   char* ccy   = swapLeg->GetCurrencyUnit()->GetCcyName();

    //   for (i = 0; i < size; i++)
    //   {
    //       ARM_Date tmp(exerciseDates->Elt(i));
    //       tmp.PreviousBusinessDay(spot, ccy);

    //       exerciseDates->Elt(i) = tmp.GetJulian();
    //   }
    //}

    //itsExerciseDatesSize = exerciseDates->GetSize();
 
    //SetExerciseStartDates(exerciseDates);
 
    //SetExerciseEndDates(exerciseDates);

    //delete exerciseDates;
}



/*----------------------------------------------------------------------------*
    Assignment operator
/*----------------------------------------------------------------------------*/

ARM_ExerciseStyle& ARM_ExerciseStyle::operator = 
                                           (const ARM_ExerciseStyle& xStyle)
{
    (*this).ARM_Object::operator = (xStyle);

    BitwiseCopy(&xStyle);   

    return(*this);
}



/*----------------------------------------------------------------------------*
    Constructor    (copy).
*----------------------------------------------------------------------------*/

ARM_ExerciseStyle::ARM_ExerciseStyle(ARM_ExerciseStyle &xStyle) 
                     : ARM_Object(xStyle)
{
    Init();

    BitwiseCopy(&xStyle);   
}




/*----------------------------------------------------------------------------*
    Check if exercise is allowed at the specified date
*----------------------------------------------------------------------------*/ 

int ARM_ExerciseStyle::IsExerciseDate(long date)
{
    int flag = 0, i = 0;

    double* itsExerciseEndDatesElt = itsExerciseEndDates->GetElt();
    double* itsExerciseStartDatesElt = itsExerciseStartDates->GetElt();


    if (( itsExerciseType == K_EUROPEAN ) || ( itsExerciseType == K_BERMUDAN ))
    {
       while(( i < itsExerciseDatesSize) 
             && ( date != itsExerciseEndDatesElt[i] ))
       {
           i++;
       }

       if ( i < itsExerciseDatesSize ) 
          flag = 1;
    }

    if (( itsExerciseType == K_AMERICAN ) 
        || 
        ( itsExerciseType == K_CUSTOMIZED )
       )
    {
       while (( i < itsExerciseDatesSize ) && 
              ( date <= itsExerciseStartDatesElt[i] ) && 
              ( date >= itsExerciseEndDatesElt[i] ))
       {
           i++;
       }

       if ( i < itsExerciseDatesSize ) 
          flag = 1;
    }

    return(flag);
}



void ARM_ExerciseStyle::View(char* id, FILE* ficOut)
{
    FILE* fOut;
    char fOutName[200];
 
 
    if ( ficOut == NULL )
    {
       ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);

       fOut = fopen(fOutName, "w");
    }
    else
    {
        fOut = ficOut;
    }

    switch(itsExerciseType)
    {
        case K_AMERICAN:
        {
            fprintf(fOut,"\n Exercise Style : AMERICAN\n");
        };
        break;

        case K_CUSTOMIZED:
        {
            fprintf(fOut,"\n Exercise Style : CUSTOMIZED\n");
        };
        break;

        case K_EUROPEAN:
        {
            fprintf(fOut,"\n Exercise Style : EUROPEAN\n");
        };
        break;

        case K_BERMUDAN:
        {
            fprintf(fOut,"\n Exercise Style : BERMUDAN\n");
        };
        break;       
    }

    fprintf(fOut, "\n\n\t =====> Exercise Style Dates \n\n");

 
    int sz = itsExerciseDatesSize;
    char d1[20];
    char d2[20];
    int  i;
 
    fprintf(fOut, "\n\nStart Date\t\tEnd Date \n\n");
   
    for (i = 0; i < sz; i++)
    {
        ARM_Date startDate, endDate;
       

        startDate = ARM_Date (itsExerciseStartDates->Elt(i));

        startDate.JulianToStrDate(d1);
 
        endDate = ARM_Date (itsExerciseEndDates->Elt(i));

        endDate.JulianToStrDate(d2);
 
        fprintf(fOut, "%s\t\t%s\n", d1, d2);
    }

    if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}




/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
