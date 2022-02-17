/*
 * $Log: zerospl3.h,v $
 * Revision 1.4  2002/11/15 10:16:38  mab
 * correction of the object name : SetName(ARM_ZERO_SPLICUB);
 * instead of SetName(ARM_ZERO_LIN_INTERPOL);
 *
 * Revision 1.3  2002/10/11 08:18:24  mab
 * Improvements
 *
 */


/*----------------------------------------------------------------------------*
    zerospl3.h

    
    Header for the ARM_ZeroSpliCub class, a class for computing 
    a ARM_ZeroCurve using spline cubique interpolation of an input zero curve 
    defined by a finite set of yearterms and yields.
 
*----------------------------------------------------------------------------*/ 
#ifndef _ZEROSPL3_H
#define _ZEROSPL3_H




#include "dates.h"
#include "zerocurv.h"




class ARM_ZeroSpliCub : public ARM_ZeroCurve 
{
    private:

        ARM_Vector* itsD2Rates;


        // the input year terms and yields vector must have 3 elements
        // in order to compute 1st and 2nd 
        // derivatives -> K_MIN_NUM_YEAR_TERMS=3

        int itsCompoundMeth;  // Compounding method of the zero yields


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
            SetName(ARM_ZERO_SPLICUB);

            itsLastBucketInt = 0;

            itsD2Rates = NULL;
        }

    public:

        ARM_ZeroSpliCub(void)
        {
            Init();
        }

        ARM_ZeroSpliCub(ARM_Date& asOf, 
                        ARM_Vector* yearTerms, 
                        ARM_Vector* yields, 
                        int compMeth = K_COMP_CONT,
                        int lastBucketInt = 0);

        ARM_ZeroSpliCub(ARM_Date& asOf, 
            ARM_CRV_TERMS& terms, 
            ARM_Vector* mktData, 
            int MMVsFut = 0, 
            int SwapVsFut = 0,
            int raw = K_PAR, 
            int Cont_lin = K_LINEAR,
            int lastBucketInt = 0, 
            ARM_Currency* ccy = ARM_DEFAULT_CURRENCY);

        ARM_ZeroSpliCub(ARM_Date& asOf, 
            ARM_CRV_TERMS& terms, 
            ARM_Vector* mktData, 
            ARM_Container* bonds, 
            ARM_Vector* yields,
            int MMVsFut, 
            int lastBucketInt = 0, 
            ARM_Currency* ccy = ARM_DEFAULT_CURRENCY);


        ARM_ZeroSpliCub(ARM_ZeroSpliCub &);
    
       ~ARM_ZeroSpliCub(void)
        {
           if (itsD2Rates)
              delete itsD2Rates;
        }

        ARM_ZeroSpliCub& operator = (ARM_ZeroSpliCub &);


        void BitwiseCopy(const ARM_Object* srczspl)
        {
            ARM_ZeroSpliCub* zspl = (ARM_ZeroSpliCub*) srczspl;

            itsCompoundMeth = zspl->itsCompoundMeth;

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
            ARM_ZeroSpliCub* theClone = new ARM_ZeroSpliCub();


            theClone->Copy(this);

            return(theClone);
        }

        void GenerateD2Rates();

        void SetCompoundMeth(int meth)
        {
            itsCompoundMeth = meth;
        }

        virtual int GetCompoundMeth(void) 
        {
            return(itsCompoundMeth);
        }

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
