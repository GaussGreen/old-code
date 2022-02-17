/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : smilecrv.h                                                   */
/*                                                                            */
/* DESCRIPTION : Header for the ARM_SmileCurve class, a class                 */
/*               for dealing with options smiles                              */ 
/*                                                                            */
/* DATE        : Tue Apr  1 1997                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/
#ifndef _SMILE_CRV_H
#define _SMILE_CRV_H



#include <stdio.h>


#include "armglob.h"
#include "dates.h"
#include "volcurv.h"




class ARM_Date;





class ARM_SmileCurve : public ARM_VolCurve
{
    private:

        ARM_Vector* itsStrikes;
        ARM_Vector* itsVolats;


        double VolatilityFunction(double m, double K, double m2 = 0.0);

        void Init(void)
        {
            itsStrikes = NULL;
            itsVolats  = NULL;
        }

    public:

        ARM_SmileCurve(ARM_SmileCurve& volCurve);

        ARM_SmileCurve(ARM_Date& date, ARM_Vector* stk, 
                       ARM_Vector* volats) : ARM_VolCurve(date)
        {
            SetName(ARM_SMILE_CURVE);
 
            SetStrikeType(K_STK_TYPE_MATU_SMILE);

            itsStrikes = NULL;
 
            itsVolats  = NULL;

            if (stk)
               itsStrikes = new ARM_Vector(*stk);
            
            if (volats)
               itsVolats  = new ARM_Vector(*volats);
        }

        ARM_SmileCurve(void)
        {
            SetName(ARM_SMILE_CURVE);
 
            SetStrikeType(K_STK_TYPE_MATU_SMILE);

            itsStrikes = NULL;

            itsVolats  = NULL;
        }

       virtual ~ARM_SmileCurve(void)
        {
            if (itsStrikes)
            {
               delete itsStrikes;
               itsStrikes = NULL;
            }
 
            if (itsVolats)
            {
               delete itsVolats;
               itsVolats  = NULL;
            }
        }

        void Set(ARM_Date& date, ARM_Vector* stk, ARM_Vector* volats)
        {
            if (itsStrikes)
            {
               delete itsStrikes;
               itsStrikes = NULL;
            }
 
            if (itsVolats)
            {
               delete itsVolats;

               itsVolats = NULL;
            }

            SetAsOfDate(date);

            itsStrikes = stk;
            itsVolats  = volats;
        }

        void SetStrikes(ARM_Vector* stk)
        {
            if (itsStrikes)
            {
               delete itsStrikes;
               itsStrikes = NULL;
            }

            if (stk)
               itsStrikes = new ARM_Vector(*stk);
        }

        void SetVolats(ARM_Vector* volats)
        {
            if (itsVolats)
            {
               delete itsVolats;
               itsVolats  = NULL;
            }
 
            if (volats)
               itsVolats  = new ARM_Vector(*volats);
        }

        ARM_SmileCurve& operator = (ARM_SmileCurve& volCurve);
  
        void BitwiseCopy(const ARM_Object* srcVolCurve)
        {
            ARM_SmileCurve* vCurve = (ARM_SmileCurve *) srcVolCurve;
 
            if (itsStrikes)
            {
               delete itsStrikes;
               itsStrikes = NULL;
            }
 
            if (itsVolats)
            {
               delete itsVolats;
               itsVolats  = NULL;
            }
 
            if (vCurve->itsStrikes)
               itsStrikes = new ARM_Vector(*vCurve->itsStrikes);
            
            if (vCurve->itsVolats)
               itsVolats  = new ARM_Vector(*vCurve->itsVolats); 
        }
 
        void Copy(const ARM_Object* vCurve)
        {
            ARM_VolCurve::Copy(vCurve);
 
            BitwiseCopy(vCurve);
        }
 
        virtual ARM_Object* Clone(void)
        {
            ARM_SmileCurve* theClone = new ARM_SmileCurve();
 
 
            theClone->Copy(this);
 
            return(theClone);
        }
};


#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
