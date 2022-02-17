/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : oldvolcurve.h                                                */
/*                                                                            */
/* DESCRIPTION : Header for the ARM_OldVolCurve class,						  */ 
/*               A curve keeping the same volatility from a day to another.	  */
/*                                                                            */
/* DATE        : Tue April  24 2007                                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/

#ifndef _OLD_VOL_CURVE_H
#define _OLD_VOL_CURVE_H

#include "volcurv.h"

class ARM_OldVolCurve : public ARM_VolCurve
{
    private:
		ARM_VolCurve*	itsVolCurve;
		ARM_Date		itsAsOfDate;

    public:

        ARM_OldVolCurve(ARM_Date& asOfDate, ARM_VolCurve* crv);
		~ARM_OldVolCurve() { delete itsVolCurve; }
        
        ARM_OldVolCurve(const ARM_OldVolCurve& rhs);
		ARM_OldVolCurve& operator=(const ARM_OldVolCurve& rhs);

		virtual ARM_Object* Clone(void);
		virtual void View(char* id = NULL, FILE* ficOut = NULL);

		virtual double computeVol(double moneyness, double maturity);
		virtual double VolatilityFunction(double m1, double m2);
		virtual ARM_VolLInterpol* ComputeATMFxVol(double lag=0.0);
		virtual int GetSmileFlag(void) const;
		virtual ARM_Vector* GetExpiryTerms(void) const;
};

#endif