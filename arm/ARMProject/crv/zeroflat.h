/*----------------------------------------------------------------------------*

     zeroflat.h
 
    Header for a flat curve
 
*----------------------------------------------------------------------------*/


#ifndef _ZEROFLAT_H
#define _ZEROFLAT_H

#include "zerocurv.h"



class ARM_ZeroFlat : public ARM_ZeroCurve 
{
    private:

        double       itsYield; //    value of flat rate  


    //  Private methods

    double DiscountFunction(double yearTerm);

    double D1DiscountFunction(double yearTerm)
    {
        return(0.01*itsYield*exp(-yearTerm*itsYield));
    }

    double D2DiscountFunction(double yearTerm)
    {
        return(0.01*itsYield*itsYield*exp(-yearTerm*itsYield));
    }


    public:

        ARM_ZeroFlat(void)
        {
            SetName(ARM_ZERO_FLAT);

            itsYield = 0.0;
        }

        ARM_ZeroFlat(ARM_Date& asOf, double flatRate,
                     ARM_Currency* ccy = ARM_DEFAULT_CURRENCY);

        ARM_ZeroFlat(ARM_ZeroFlat &flatCrv);

       ~ARM_ZeroFlat(void){}
    
virtual double CalcNumericalObjectSignature(void)
        {
             double signature = 0.0;

             signature += GetAsOfDate().GetJulian()*itsYield;

             return(signature);
        }

        void BitwiseCopy(const ARM_Object* srcZf)
        {
            ARM_ZeroFlat* zf = (ARM_ZeroFlat *) srcZf;
 
 
            itsYield = zf->itsYield;
        }
 
        void Copy(const ARM_Object* srcZf)
        {
            ARM_ZeroCurve::Copy(srcZf);
 
            BitwiseCopy(srcZf);
        }
 
        ARM_Object* Clone(void)
        {
            ARM_ZeroFlat* theClone = new ARM_ZeroFlat();
 
 
            theClone->Copy(this);
 
            return(theClone);
        }

        ARM_ZeroFlat & operator = (ARM_ZeroFlat &flatCrv);

        
        double GetFlatRate(void)
        {
            return(itsYield);
        }
 
		void Set(ARM_Date& asOf, double flatRate);

        void SetFlatRate(double flatRate)
        {
            itsYield = flatRate;
        }
};




#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/