/*
 * $Log: zerospli.h,v $
 * Revision 1.2  2002/05/30 13:33:14  mab
 * introducing RCS Mark
 *
 */


/*----------------------------------------------------------------------------*

    zerospli.h
 
    Header for the ARM_ZeroSplines class, a class for computing 
    a ARM_ZeroCurve using FOM's research group splines.

*----------------------------------------------------------------------------*/
#ifndef _ZEROPLI_H
#define _ZEROPLI_H




#include "zerocurv.h"

#include <string.h>




#define    K_MAX_NUM_SPLINE_COEFFS 40






class ARM_ZeroSplines : public ARM_ZeroCurve 
{
    protected:

        int    itsNumSplineCoeffs; // number of splines

        double itsSplineCoeffs[K_MAX_NUM_SPLINE_COEFFS]; // spline coefficients
    
        double itsInt1;
        double itsInt2;
        double itsInt3;

    
        // Private methods

        double DiscountFunction(double yearTerm);
        double D1DiscountFunction(double yearTerm);
        double D2DiscountFunction(double yearTerm);

        void Init(void)
        {
            memset(itsSplineCoeffs, 0, sizeof(itsSplineCoeffs));

            itsInt1 = 0.0;
            itsInt2 = 0.0;
            itsInt3 = 0.0;

            itsNumSplineCoeffs = 0;
        }

    public:

        ARM_ZeroSplines(void)
        {
           Init();

           SetName(ARM_ZERO_SPLINES);

           itsNumSplineCoeffs = 0;
        }


        ARM_ZeroSplines(ARM_Date& asOf, int numSplineCoeffs, 
                        double* splineCoeffs,
                        double Int1=2.0,double Int2=7.0,double Int3=10.0);

        ARM_ZeroSplines(ARM_ZeroSplines &zeroSplines);

       ~ARM_ZeroSplines(void){}
    
        ARM_ZeroSplines& operator = (ARM_ZeroSplines& zeroSplines);

        
        void BitwiseCopy(const ARM_Object* srcZspl)
        {
            ARM_ZeroSplines* zspl = (ARM_ZeroSplines*) srcZspl;


            itsNumSplineCoeffs = zspl->itsNumSplineCoeffs;

            for (int i = 0; i < itsNumSplineCoeffs; i++)
            {
                itsSplineCoeffs[i] = zspl->itsSplineCoeffs[i];
            }

            itsInt1 = zspl->itsInt1;
            itsInt2 = zspl->itsInt2;
            itsInt3 = zspl->itsInt3;
        }

        void EstimeBracket(ARM_Container* Assets, ARM_Matrix* data,int *CST);

        void Copy(const ARM_Object* srcZspl)
        {
            ARM_ZeroCurve::Copy(srcZspl);

            BitwiseCopy(srcZspl);
        }


        ARM_Object* Clone(void)
        {
            ARM_ZeroSplines* theClone = new ARM_ZeroSplines();


            theClone->Copy(this);
 
            return(theClone);
        }


        int GetNbParams(void)
        {
            return(itsNumSplineCoeffs);
        }
 
        void SetParams(double* para)
        {
            for (int i = 0; i < itsNumSplineCoeffs; i++)
            {
                itsSplineCoeffs[i] = para[i];
            }

            ARM_Vector* parameters = new ARM_Vector(itsNumSplineCoeffs,
                                                    para);

            SetParameters(parameters);
        }

        void SetParameters(ARM_Vector* param)
        {
            ARM_ZeroCurve::SetParameters(param);

            if (!param)
               return;

            if (param->GetSize() != itsNumSplineCoeffs)
            {
               throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                         "Error <SetParameters> method");
            }

            for (int i = 0; i < itsNumSplineCoeffs; i++)
            {
                itsSplineCoeffs[i] = (*param)[i];
            }
        }
 
        // a completer si besoin avec Int1,Int2,Int3
        void Set(int sz, double* coeffs, ARM_Date& asOf)
        {
            SetAsOfDate(asOf);

            itsNumSplineCoeffs = sz;

            for (int i = 0; i < itsNumSplineCoeffs; i++)
            {
                itsSplineCoeffs[i] = coeffs[i];
            }

            ARM_Vector* param = new ARM_Vector(itsNumSplineCoeffs, coeffs);

            SetParameters(param);

            GenerateFields();
        }
        
        void SetBracket(double Int1=2.0, double Int2=8.0, double Int3=11.0)
        {
            
            itsInt1 = Int1;
            itsInt2 = Int2;
            itsInt3 = Int3;
        }

        void GetBracket(double *Int1, double *Int2, double *Int3)
        {
            *Int1=itsInt1;
            *Int2=itsInt2;
            *Int3=itsInt3;
        }

        double GetModelFactor(int factorId)
        {
             if ( factorId < GetNbParams() )
             {
                return(itsSplineCoeffs[factorId]);
             }
             else
             {
                printf("\n ??? Factor %d doesnt exist !!!\n", factorId);

                return(0.0);
             }
        }


    virtual int GetNumSplineCoeffs(void);    
    virtual int GetSplineCoeffs(double *);
};




#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
