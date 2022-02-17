/*
 * $Log: zerocbdf.h,v $
 * Revision 1.3  2002/10/11 08:23:31  mab
 * Improvements
 *
 */


/*----------------------------------------------------------------------------*

    zerocbdf.h

*----------------------------------------------------------------------------*/ 
#ifndef _ZEROCBDF_H
#define _ZEROCBDF_H




#include "dates.h"
#include "zerocurv.h"




class ARM_ZeroCubDiff : public ARM_ZeroCurve 
{
    private:

        ARM_Matrix* itsPolyCoef;


        int itsCompoundMeth;  // Compounding method of the zero yields


        // Methods

        virtual double DiscountFunction(double);
        virtual double D1DiscountFunction(double);
        virtual double D2DiscountFunction(double);
        
        double ComputeInterpolValue(double x);
        void ComputePolynomCoef(void);

        void Init(void)
        {
            SetName(ARM_ZERO_CUBDIFF);

            itsPolyCoef = NULL;
        }

    public:

        ARM_ZeroCubDiff(void)
        {
            Init();

            SetName(ARM_ZERO_CUBDIFF);
        }

        ARM_ZeroCubDiff(ARM_Date& asOf, 
                        ARM_Vector* yearTerms, 
                        ARM_Vector* yields, 
                        int compMeth = K_COMP_CONT);

        ARM_ZeroCubDiff(ARM_Date& asOf, 
                        ARM_CRV_TERMS& terms, 
                        ARM_Vector* mktData, 
                        int MMVsFut = 0, 
                        int SwapVsFut = 0,
                        int raw = K_PAR, 
                        int Cont_lin = K_LINEAR,
                        int lastBucketInt = 0, 
                        ARM_Currency* ccy = ARM_DEFAULT_CURRENCY);

        ARM_ZeroCubDiff(ARM_Date& asOf, 
                        ARM_CRV_TERMS& terms, 
                        ARM_Vector* mktData, 
                        ARM_Container* bonds, 
                        ARM_Vector* yields,
                        int MMVsFut, 
                        int lastBucketInt = 0, 
                        ARM_Currency* ccy = ARM_DEFAULT_CURRENCY);


        ARM_ZeroCubDiff(ARM_ZeroCubDiff &);
    
       ~ARM_ZeroCubDiff(void)
        {
           if (itsPolyCoef)
           {
              delete itsPolyCoef;

              itsPolyCoef = NULL;
           }
        }

        ARM_ZeroCubDiff& operator = (ARM_ZeroCubDiff &);

        void BitwiseCopy(const ARM_Object* srczspl)
        {
            ARM_ZeroCubDiff* zspl = (ARM_ZeroCubDiff*) srczspl;

            itsCompoundMeth = zspl->itsCompoundMeth;

            if (itsPolyCoef)
            {
               delete itsPolyCoef;
               itsPolyCoef = NULL;
            }

            if (zspl->GetPolyCoef())
                itsPolyCoef = (ARM_Matrix *) zspl->GetPolyCoef()->Clone();
        }


        void Copy(const ARM_Object* srczspl)
        {
            ARM_ZeroCurve::Copy(srczspl);

            BitwiseCopy(srczspl);
        }


        ARM_Object* Clone(void)
        {
            ARM_ZeroCubDiff* theClone = new ARM_ZeroCubDiff();


            theClone->Copy(this);

            return(theClone);
        }

        void SetCompoundMeth(int meth)
        {
            itsCompoundMeth = meth;
        }

        virtual int GetCompoundMeth(void) 
        {
            return(itsCompoundMeth);
        }

        ARM_Matrix* GetPolyCoef(void)
        {
            return(itsPolyCoef);
        }

        void SetPolyCoef(ARM_Matrix *PolyCoef)
        {
            if ( itsPolyCoef == PolyCoef )
               return;
 
            if (itsPolyCoef)
               delete itsPolyCoef;
 
            itsPolyCoef = PolyCoef;
        }

      ARM_ZeroCurve* GenerateShiftCurve(ARM_CRV_TERMS& Terms, ARM_Vector* epsilon);
};




#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
