/*----------------------------------------------------------------------------*

  smileint.h

    
  Header for the ARM_VolLInterpolSmile class, a class for computing 
  a volatility From smiles
 
*----------------------------------------------------------------------------*/ 
#ifndef _SMILEINT_H
#define _SMILEINT_H




#include "dates.h"
#include "volint.h"
#include "linalg.h"
#include "interpol.h"







class ARM_VolLInterpolSmile : public ARM_VolLInterpol
{
    private:


        ARM_VolCurve*** itsVolatAndSmiles; // Matrix of volatilities 
                                           // and smiles



        // Methods

        virtual double VolatilityFunction(double maturity1 , double strike, 
                                          double m2 = 0.0);

        void Init(void)
        {
            itsVolatAndSmiles = NULL;
        }

    public:

        ARM_VolLInterpolSmile(void)
        {
            SetName(ARM_SMILE_LIN_INTERPOL);

            itsStrikes      = NULL;

            itsVolatilities = NULL;
        }

        ARM_VolLInterpolSmile(ARM_Date& asOf, ARM_Vector *yearTerms, 
                              ARM_Vector* strikes, ARM_Matrix* volatilities,
                               int strikeType=k_STK_TYPE_MATU_SMILE);

        ARM_VolLInterpolSmile(ARM_VolLInterpolSmile& src);
    
        virtual ~ARM_VolLInterpolSmile(void)
        {
            if (itsVolatilities)
               delete itsVolatilities;
        }

        ARM_VolLInterpolSmile& operator = (ARM_VolLInterpolSmile &);

        void BitwiseCopy(const ARM_Object* srcVollint)
        {
            ARM_VolLInterpolSmile* vollint = (ARM_VolLInterpolSmile *) srcVollint;
 
 
            if (itsVolatilities)
            {
               delete itsVolatilities;
               itsVolatilities = NULL;
            }
 
            if (vollint->itsVolatilities)
               itsVolatilities = (ARM_Matrix *) vollint->itsVolatilities->Clone();
        }
 
        void Copy(const ARM_Object* vollint)
        {
            ARM_VolCurve::Copy(vollint);
 
            BitwiseCopy(vollint);
        }
 
        virtual ARM_Object* Clone(void)
        {
            ARM_VolLInterpolSmile* theClone = new ARM_VolLInterpolSmile();
 
 
            theClone->Copy(this);
 
            return(theClone);
        }
};




#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
