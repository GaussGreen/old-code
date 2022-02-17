/*
 * $Log: exercise.h,v $
 * Revision 1.5  2004/06/08 08:09:58  mab
 * Formatting
 *
 * Revision 1.4  2001/02/22 19:43:37  smysona
 * IsExerciseDate travaille sur julainDate (perf)
 *
 * Revision 1.3  2000/11/06 09:29:30  sgasquet
 * Retrait include swapleg.h
 *
 * Revision 1.2  1999/04/06 09:37:23  ypilchen
 * Rajout de Log pour RCS
 *
 */


/*----------------------------------------------------------------------------*
     exercice.h
 
     Header for the ARM_ExerciseStyle class, a class to deal with 
               option exercice

*----------------------------------------------------------------------------*/
#ifndef _EXERCISE_H
#define _EXERCISE_H




#include "armglob.h"
#include "dates.h"


class ARM_SwapLeg;


class ARM_ExerciseStyle : public ARM_Object 
{
    
    //    variables    

    private:

        int itsExerciseType; // K_EUROPEAN ; K_AMERICAN ; 
                             // K_BERMUDAN; K_CUSTOMIZED

        int itsExerciseDatesSize;

        ARM_Vector* itsExerciseStartDates;

        ARM_Vector* itsExerciseEndDates;



        void Init(void)
        {
            SetName(ARM_EXERCISE_STYLE);
            
			itsExerciseStartDates = NULL;

            itsExerciseEndDates   = NULL;

            itsExerciseType = K_EUROPEAN;

            itsExerciseDatesSize = 0;
        }

        void FreeStyle(void) 
        {
            if (itsExerciseStartDates)
            {
               delete itsExerciseStartDates;
               itsExerciseStartDates = NULL;
            }
 
            if (itsExerciseEndDates)
            {
               delete itsExerciseStartDates;
               itsExerciseEndDates = NULL;
            }

            Init();
        }

    public:

        ARM_ExerciseStyle(void);

        // Constructor for European Style
        ARM_ExerciseStyle(ARM_Date& xDate);  

        // Constructor for American Style
        ARM_ExerciseStyle(ARM_Date& xStartDate, ARM_Date& xEndDate);  

        // Constructor for Bermudan Style
        ARM_ExerciseStyle(ARM_Vector* xDates); 

        // Constructor for Bermudan Style with Libor type reset dates
        ARM_ExerciseStyle(ARM_Date& xStartDate, ARM_Date& xEndDate, 
                          int liborType);  

        // Constructor for Customized Exercise Option 
        ARM_ExerciseStyle(ARM_Vector* xStartDates, ARM_Vector* xEndDates); 

        /* 
           Constructor for Bermudan Style from swap leg
           or derived class (eg cap) dates
                K_FLOW_PAYMENT       1
                K_INDEX_RESET        2
                K_START_INTEREST     3
                K_END_INTEREST       4
        */

        ARM_ExerciseStyle(ARM_SwapLeg* swapLeg, int dateType, 
                          ARM_Date& firstDate = ARM_Date("01/01/1981"), 
                          ARM_Date& lastDate = ARM_Date("31/12/3000"));


        ARM_ExerciseStyle(ARM_ExerciseStyle& xStyle);

        ARM_ExerciseStyle & operator = (const ARM_ExerciseStyle& xStyle);

       ~ARM_ExerciseStyle(void)
        {
            if (itsExerciseStartDates) 
            {
               delete itsExerciseStartDates; 
               itsExerciseStartDates = NULL;
            }

            if (itsExerciseEndDates) 
            {
               delete itsExerciseEndDates; 
               itsExerciseEndDates = NULL;
            }
        }

        void BitwiseCopy(const ARM_Object* srcExe)
        {
            ARM_ExerciseStyle* exe = (ARM_ExerciseStyle *) srcExe;


            itsExerciseType = exe->itsExerciseType;
 
            itsExerciseDatesSize = exe->itsExerciseDatesSize;


            if (itsExerciseStartDates)
            {
               delete itsExerciseStartDates;
               itsExerciseStartDates = NULL;
            }
 
            if (exe->itsExerciseStartDates)
            {
               itsExerciseStartDates = 
                          (ARM_Vector *) exe->itsExerciseStartDates->Clone();
            }

            if (itsExerciseEndDates)
            {
               delete itsExerciseEndDates;
               itsExerciseEndDates = NULL;
            }

            if (exe->itsExerciseEndDates)
            {
               itsExerciseEndDates =
                          (ARM_Vector *) exe->itsExerciseEndDates->Clone();
            }
        }

        void Copy(const ARM_Object* srcExe)
        {
            ARM_Object::Copy(srcExe);
  
            BitwiseCopy(srcExe);
        }
 
        ARM_Object* Clone(void)
        {
            ARM_ExerciseStyle* theClone = new ARM_ExerciseStyle();
 
 
            theClone->Copy(this);
 
            return(theClone);
        }

        ARM_CLASS_NAME GetRootName(void)
        {
            return(ARM_EXERCISE_STYLE);
        }

        void View(char* id = NULL, FILE* fOut = NULL);

        void Set(ARM_Date& xDate);  

        void Set(ARM_Date& xStartDate, ARM_Date& xEndDate);  

        void Set(ARM_Vector* xDates); 

        void Set(ARM_Date& xStartDate, ARM_Date& xEndDate, int liborType);

        void Set(ARM_Vector* xStartDates, ARM_Vector* xEndDates); 

        void Set(ARM_SwapLeg* swapLeg, int dateType, 
                 ARM_Date& firstDate = ARM_Date("01/01/1981"), 
                 ARM_Date& lastDate = ARM_Date("31/12/3000"));


        // Check if exercise is allowed at the specified date
        
        int IsExerciseDate(long date);

        inline int IsExerciseDate(ARM_Date& date)
        {
            return(IsExerciseDate(date.GetJulian()));
        }


        int GetExerciseType(void) const
        {
            return(itsExerciseType); 
        }

        void SetExerciseType(int xStyle)
        {
            itsExerciseType = xStyle;
        }

        int GetExerciseDatesSize(void) const
        {
            return(itsExerciseDatesSize); 
        }

        void SetExerciseDatesSize(int size)
        {
            itsExerciseDatesSize = size;
        }

        ARM_Vector* GetExerciseStartDates(void) 
        {
            return(itsExerciseStartDates); 
        }

        void SetExerciseStartDates(ARM_Vector* xStartDates)
        {
            if ( itsExerciseStartDates == xStartDates )
               return;

            if (itsExerciseStartDates)
            {
               delete itsExerciseStartDates;

               itsExerciseStartDates = NULL;
            }
 
            itsExerciseStartDates = (ARM_Vector *) xStartDates->Clone();
        }

        ARM_Vector* GetExerciseEndDates(void) 
        {
            return(itsExerciseEndDates); 
        }

        void SetExerciseEndDates(ARM_Vector* xEndDates)
        {
            if ( itsExerciseEndDates == xEndDates )
               return;
 
            if (itsExerciseEndDates)
            {
               delete itsExerciseEndDates;

               itsExerciseEndDates = NULL;
            }
 
            itsExerciseEndDates = (ARM_Vector *) xEndDates->Clone();
        }
};


#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
