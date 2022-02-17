/*----------------------------------------------------------------------------*
    zerosplsum.h

    
    Header for the ARM_ZeroSpliSum class, a class for computing 
    a ARM_ZeroCurve using Summit's algorithm 
	to produced a new set of yearterms and cubic splined zero coupon yields
	from an initial set of yearterms and zero coupon yields 
	or a 'Raw' computed zero coupon yields from a Market Data input.
 
*----------------------------------------------------------------------------*/ 
#ifndef _ZEROSPLSUM_H
#define _ZEROSPLSUM_H




#include "dates.h"
#include "zerocurv.h"




class ARM_ZeroSpliSum : public ARM_ZeroCurve 
{
    private:

        ARM_Vector* itsD2Rates;


        // the input year terms and yields vector must have 3 elements
        // in order to compute 1st and 2nd 
        // derivatives -> K_MIN_NUM_YEAR_TERMS=3

        int itsLastBucketInt; // Interpolation method used for
                              // yields beyond the last yearTerm
                              // 0 for flat, 1 if the last line is used  

        // Methods

        virtual double DiscountFunction(double);
        virtual double D1DiscountFunction(double);
        virtual double D2DiscountFunction(double);
        double ComputeSpline(double x);

        void Init(void)
        {
            SetName(ARM_ZERO_SPLSUM);

            itsLastBucketInt = 0;

            itsD2Rates = NULL;
        }

    public:

        ARM_ZeroSpliSum(void)
        {
            Init();
        }

        ARM_ZeroSpliSum(ARM_Date& asOf, 
                        ARM_Vector* yearTerms, 
                        ARM_Vector* yields,
                        int lastBucketInt = 0);

        ARM_ZeroSpliSum(ARM_Date& asOf, 
            ARM_CRV_TERMS& terms, 
            ARM_Vector* mktData, 
            int MMVsFut = 0, 
            int SwapVsFut = 0,
            int raw = K_RAW, 
            int Cont_lin = K_LINEAR,
            int lastBucketInt = 0, 
            ARM_Currency* ccy = ARM_DEFAULT_CURRENCY);

        ARM_ZeroSpliSum(ARM_Date& asOf, 
            ARM_CRV_TERMS& terms, 
            ARM_Vector* mktData, 
            ARM_Container* bonds, 
            ARM_Vector* yields,
            int MMVsFut, 
            int lastBucketInt = 0, 
            ARM_Currency* ccy = ARM_DEFAULT_CURRENCY);


        ARM_ZeroSpliSum(ARM_ZeroSpliSum &);
    
       ~ARM_ZeroSpliSum(void)
        {
           if (itsD2Rates)
              delete itsD2Rates;
        }

        ARM_ZeroSpliSum& operator = (ARM_ZeroSpliSum &);


        void BitwiseCopy(const ARM_Object* srczspl)
        {
            ARM_ZeroSpliSum* zspl = (ARM_ZeroSpliSum*) srczspl;

            itsLastBucketInt = zspl->itsLastBucketInt;

            if (itsD2Rates)
            {
               delete itsD2Rates;
               itsD2Rates = NULL;
            }

            if (zspl->GetD2Rates())
                itsD2Rates = (ARM_Vector *) zspl->GetD2Rates()->Clone();
        }


        void Copy(const ARM_Object* srczspl)
        {
            ARM_ZeroCurve::Copy(srczspl);

            BitwiseCopy(srczspl);
        }


        ARM_Object* Clone(void)
        {
            ARM_ZeroSpliSum* theClone = new ARM_ZeroSpliSum();


            theClone->Copy(this);

            return(theClone);
        }

        void SumGenerateD2Rates();

        
        ARM_Vector *GetD2Rates(void)
        {
            return itsD2Rates;
        }

        void SetD2Rates(ARM_Vector *D2Rates)
        {
            if ( itsD2Rates == D2Rates )
               return;

            if (itsD2Rates)
               delete itsD2Rates;

            itsD2Rates = D2Rates;
        }

      ARM_ZeroCurve* GenerateShiftCurve(ARM_CRV_TERMS& Term, ARM_Vector* epsilon);
        
};




#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
