/*
 * $Log: zerointimp.h,v $
 * Revision 1.7  2002/10/11 08:19:34  mab
 * Improvements
 *
 * Revision 1.6  2002/08/08 10:24:31  mab
 * Formattage
 *
 */


/*----------------------------------------------------------------------------*
    zeroint.h

    
    Header for the ARM_ZeroInterpolation class, a class for computing 
    a ARM_ZeroCurve using implicit interpolation of an input zero curve 
    defined by a finite set of yearterms and yields.
 
*----------------------------------------------------------------------------*/ 
#ifndef _ZEROINTERP_H
#define _ZEROINTERP_H




#include "dates.h"
#include "zerocurv.h"
#include "linalg.h"
#include "zeroint.h"



        // the input year terms and yields vector must have 3 elements
        // in order to compute 1st and 2nd 
        // derivatives -> K_MIN_NUM_YEAR_TERMS=3



#define    K_MIN_NUM_YEAR_TERMS        3



class ARM_ZeroInterpolation : public ARM_ZeroCurve 
{
    private:

        ARM_Matrix* itsPolyCoef; // polynome de splines cubique


        int itsCompoundMeth;  // Compounding method of the zero yields

        int itsLastBucketInt; // Interpolation method used for
                              // yields beyond the last yearTerm
                              // 0 for flat, 1 if the last line is used  

        int itsLambda;   // Interpolation parameter 
        int itsNbPoints; // Interpolation precision

        // Methods

        virtual double DiscountFunction(double);
        virtual double D1DiscountFunction(double);
        virtual double D2DiscountFunction(double);

    /* fonctions de lissage non parametrique */
    void NonParamCurve(ARM_Vector* terms,
                       ARM_Vector* yields,int lambda,
                       int nbPoints,
                       ARM_Vector* newterms,
                       ARM_Vector* newyields);

    double Regularisation(double xkm); // fonction de regularisation

    // fonctions mathematiques pour calculer la hessienne
    void Calc_Hess(ARM_Vector* newterms,ARM_Matrix* Hess);

    void Remp_Theta(ARM_Vector* newterms,ARM_Matrix* theta);

    void factLU(ARM_Matrix* Hess,ARM_Matrix* L,ARM_Matrix* U);
    void resoud(ARM_Matrix* L, ARM_Matrix* U,
                ARM_Vector* Grad,ARM_Vector* Solution);


    /* fonctions d'interpolation splines cubiques */

    double ComputeInterpolValue(double x); 
    void ComputePolynomCoef(void);

    void Init(void)
    {
        itsLastBucketInt = 0;

        if (!itsLambda)
           itsLambda = 500000;

        if (!itsNbPoints)
           itsNbPoints = 300;

        itsPolyCoef = NULL;
    }

    public:

        ARM_ZeroInterpolation(void)
        {
            Init();

            SetName(ARM_ZERO_INTERPOLATION);
        }

        ARM_ZeroInterpolation(ARM_Date& asOf, 
                              ARM_Vector* yearTerms, 
                              ARM_Vector* yields,
                              int compoundMeth = K_COMP_CONT,
                              int lambda = 500000,
                              int nbPoints = 300);

        ARM_ZeroInterpolation(ARM_Date& asOf, 
                              ARM_CRV_TERMS& Terms, 
                              ARM_Vector* mktData, 
                              int MMVsFut = 0, 
                              int SwapVsFut = 0,
                              int raw = K_PAR, 
                              int Cont_lin = K_LINEAR,
                              ARM_Currency* ccy = ARM_DEFAULT_CURRENCY,
                              int lambda = 500000,
                              int nbPoints = 300);

       ARM_ZeroInterpolation(ARM_CRV_TERMS& Terms,
                             ARM_Date& asOf,
                             ARM_Vector* mktData,
                             int MMVsFut = 0,
                             int SwapVsFut = 0,
                             int raw = K_PAR, 
                             int Cont_lin = K_LINEAR,
                             ARM_Currency* ccy = ARM_DEFAULT_CURRENCY,
                             int lambda = 500000,
                             int nbPoints = 300);

        ARM_ZeroInterpolation(ARM_Date& asOf, 
                              ARM_CRV_TERMS& Terms, 
                              ARM_Vector* mktData, 
                              ARM_Container* bonds, 
                              ARM_Vector* Yields,
                              int MMVsFut, 
                              ARM_Currency* ccy = ARM_DEFAULT_CURRENCY,
                              int lambda = 500000,
                              int nbPoints = 300);


        ARM_ZeroInterpolation(ARM_Date& asOf,
                              char* Terms[ARM_NB_TERMS],
                              ARM_Vector* mktData,
                              double mean_rates,
                              int raw,
                              int Cont_Lin,
                              ARM_Currency* ccy,
                              int lambda = 500000,
                              int nbPoints = 300);

        ARM_ZeroInterpolation(const ARM_ZeroInterpolation &);

        ARM_ZeroInterpolation(ARM_ZeroLInterpol* inCurve,
                              int lambda = 500000, 
                              int nbPoints = 300);

    
       ~ARM_ZeroInterpolation(void)
        {
            if (itsPolyCoef)
            {
               delete itsPolyCoef;

               itsPolyCoef = NULL;
            }
        }

        ARM_ZeroInterpolation& operator = (const ARM_ZeroInterpolation &);

        ARM_ZeroCurve* GenerateShiftCurve(ARM_CRV_TERMS& Term, ARM_Vector* epsilon);

        void BitwiseCopy(const ARM_Object* srczint)
        {
            ARM_ZeroInterpolation* zint = (ARM_ZeroInterpolation*) srczint;

            itsCompoundMeth = zint->itsCompoundMeth;

            itsLastBucketInt = zint->itsLastBucketInt;

            itsLambda = zint->itsLambda;

            itsNbPoints = zint->itsNbPoints;

            if (itsPolyCoef) 
            {
               delete itsPolyCoef;

               itsPolyCoef = NULL;
            }

            if (zint->itsPolyCoef)
               itsPolyCoef = (ARM_Matrix *) zint->itsPolyCoef->Clone();
        }


        void Copy(const ARM_Object* srczint)
        {
            ARM_ZeroCurve::Copy(srczint);

            BitwiseCopy(srczint);
        }


        ARM_Object* Clone(void)
        {
            ARM_ZeroInterpolation* theClone = new ARM_ZeroInterpolation();

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
        
        // Ne pas utiliser (ne prend pas en compte le lastBucketInt
        void Set(int sz, double* yt, double* zy, int meth, char* date)
        {
            ARM_Vector* terms = new ARM_Vector(sz, yt);
            
            ARM_Vector* rates = new ARM_Vector(sz, zy);
            
            SetYearTerms(terms);

            SetZeroRates(rates);

            int size = GetYearTerms()->GetSize();
 
            ARM_Vector* BPShifts = new ARM_Vector(size, 0.0);

            SetBPShifts(BPShifts);

            this->SetAsOfDate((ARM_Date) date);

            this->SetCompoundMeth(meth);

            GenerateDateTerms();
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
};


#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
