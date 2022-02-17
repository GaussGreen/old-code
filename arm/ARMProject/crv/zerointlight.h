/*
 * $Log: zerointlight.h,v $
 * Revision 1.2  2003/09/25 09:17:10  mab
 * "constification"
 *
 */

/*----------------------------------------------------------------------------*
    zerointlight.h

    
    Header for the ARM_ZeroLInterpolLight class, a class for computing 
    a ARM_ZeroLInterpol using linear interpolation of an input zero curve 
    defined by a finite set of yearterms and yields.
    The fonctions are specific in order to speed up the computations
 
*----------------------------------------------------------------------------*/ 
#ifndef _ZEROINTLIGHT_H
#define _ZEROINTLIGHT_H



#include "zeroint.h"


class ARM_ZeroLInterpolLight : public ARM_ZeroLInterpol 
{
    private:

        double* itsZCRates;
        double* itsTerms;

        int     SizeVect;

        // Method
        double DiscountFunction(double);

    public:

        void Init(void)
        {
            SetName(ARM_ZERO_LIN_INTERPOL);

            itsZCRates = NULL;
            itsTerms   = NULL;

            SizeVect   = 0.0;
        }

        ARM_ZeroLInterpolLight(void)
        {
            Init();
        }

        ARM_ZeroLInterpolLight(double* yearTerms, 
                               double* yields, int size);

        ARM_ZeroLInterpolLight(const ARM_ZeroLInterpolLight& zcL);
        
       ~ARM_ZeroLInterpolLight(void)
        {
           if (itsZCRates)
              delete itsZCRates;

           if (itsTerms)
              delete itsTerms;
        }

        ARM_ZeroLInterpolLight& operator = (const ARM_ZeroLInterpolLight &);

        void Copy(const ARM_Object* srczint)
        {
            ARM_ZeroCurve::Copy(srczint);

            BitwiseCopy(srczint);
        }

        ARM_Object* Clone(void)
        {
            ARM_ZeroLInterpolLight* theClone = new ARM_ZeroLInterpolLight();

            theClone->Copy(this);

            return(theClone);
        }

        void Set(double* yearTermsVect, double* ratesVect,
                 int size)
        {
            SizeVect = size;

            SetTerms(yearTermsVect);

            SetRates(ratesVect);
        }

        // Methode sale mais rapide
        void SetTerms(double* yearTermsVect)
        {
            if (itsTerms)
               delete itsTerms;

            itsTerms = yearTermsVect;
        }

        void SetRates(double* zeroRatesVect)
        {
            if (itsZCRates)
               delete itsZCRates;

            itsZCRates = zeroRatesVect;
        }
};


#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
