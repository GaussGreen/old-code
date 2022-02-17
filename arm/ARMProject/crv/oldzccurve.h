/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : oldzccurve.h                                                   */
/*                                                                            */
/* DESCRIPTION : Header for the ARM_OldZcCurve class,							  */ 
/*               A curve keeping the discount from a day to another.		  */
/*                                                                            */
/* DATE        : Tue April  24 2007                                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/

#ifndef _OLD_CURVE_H
#define _OLD_CURVE_H

#include "zerocurv.h"


class ARM_OldZcCurve : public ARM_ZeroCurve
{
    private:
		ARM_ZeroCurve*	itsZeroCurve;
		ARM_Date		itsAsOfDate;

    public:

        ARM_OldZcCurve(ARM_Date& asOfDate, ARM_ZeroCurve* crv);
		~ARM_OldZcCurve() { delete itsZeroCurve; }
        
        ARM_OldZcCurve(const ARM_OldZcCurve& rhs);
		ARM_OldZcCurve& operator=(const ARM_OldZcCurve& rhs);

		virtual ARM_Object* Clone(void);
		virtual void View(char* id = NULL, FILE* ficOut = NULL);

		virtual double DiscountFunction(double);
};


#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
