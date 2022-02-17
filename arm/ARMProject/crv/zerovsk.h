/*----------------------------------------------------------------------------*

  zerovsk.h
 
  Header for the ARM_ZeroVasicek class, a class for computing a ARM_ZeroCurve
               using CDC's research group .

*----------------------------------------------------------------------------*/


#ifndef _ZEROVSK_H
#define _ZEROVSK_H



#include "zerocurv.h"







class ARM_ZeroVasicek : public ARM_ZeroCurve 
{
    private:

        double DiscountFunction(double yt);
   
        double D1DiscountFunction(double yt); 

    public:

        ARM_ZeroVasicek(void);
  
        ARM_ZeroVasicek(ARM_Date& asOf, double* param);

        ARM_ZeroVasicek(ARM_Date& asOf);

        ARM_ZeroVasicek(ARM_ZeroVasicek& );
    
       ~ARM_ZeroVasicek(void)
        {
        }

        ARM_ZeroVasicek& operator = (ARM_ZeroVasicek &);

        void BitwiseCopy(const ARM_Object* srcZvas)
        {
            ARM_ZeroVasicek* zvas = (ARM_ZeroVasicek*) srcZvas;

            if (zvas->GetParameters())
                SetParameters((ARM_Vector *) zvas->GetParameters()->Clone());
        }

        void Copy(const ARM_Object* srcZvas)
        {
            ARM_ZeroCurve::Copy(srcZvas);

            BitwiseCopy(srcZvas);
        }


        ARM_Object* Clone(void)
        {
            ARM_ZeroVasicek* theClone = new ARM_ZeroVasicek();

            theClone->Copy(this);

            return(theClone);
        }

        int GetNbParams(void)
        {
            return(4);
        }
 
        void SetScaleFactor(double a)
        {
            GetParameters()->Elt(0) = a;
        }
 
        void SetParams(double* para)
        {
            ARM_Vector* parameters = new ARM_Vector(4, para);

            SetParameters(parameters);
        }


        void SetParameters(ARM_Vector* param)
        {
            ARM_ZeroCurve::SetParameters(param);

            if (!param)
               return;

            if (param->GetSize() != 4)
            {
               throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                         "Error <SetParameters> method");
            }

        }

 
        void Set(double* coeffs, ARM_Date& asOfDate)
        {
            SetAsOfDate(asOfDate);

            ARM_Vector* parameters = new ARM_Vector(4, &coeffs[0]);

            SetParameters(parameters);

            GenerateFields();
        }

};



double g1(double a, double theta);
 
 
double g2(double a, double theta);


#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
