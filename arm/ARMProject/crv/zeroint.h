/*
 * $Log: zeroint.h,v $
 * Revision 1.18  2003/09/25 09:16:28  mab
 * "Constification"
 *
 * Revision 1.17  2003/09/03 18:07:00  jpriaudel
 * ajout d'un constructeur
 *
 * Revision 1.16  2003/03/19 16:21:49  sgasquet
 * Ajout de la ccy dans le premier constructeur
 *
 * Revision 1.15  2002/10/11 08:20:46  mab
 * Improvements
 *
 * Revision 1.14  2002/09/24 13:18:23  mab
 * Added : TOY management (See TOY section)
 *
 * Revision 1.13  2002/05/30 13:37:32  mab
 * in char terms[ARM_NB_TERMS][6] : 6 replaced by 12
 *
 * Revision 1.12  2001/11/06 11:55:23  mab
 * rajout Commentaire RCS
 *
 */


/*----------------------------------------------------------------------------*
    zeroint.h

    
    Header for the ARM_ZeroLInterpol class, a class for computing 
    a ARM_ZeroCurve using linear interpolation of an input zero curve 
    defined by a finite set of yearterms and yields.
 
*----------------------------------------------------------------------------*/ 
#ifndef _ZEROINT_H
#define _ZEROINT_H


#include "dates.h"
#include "zerocurv.h"
#include "linalg.h"
#include "matlib.h"



#define    K_MIN_NUM_YEAR_TERMS        3



class ARM_ZeroLInterpol : public ARM_ZeroCurve 
{
    private:

        // the input year terms and yields vector must have 3 elements
        // in order to compute 1st and 2nd 
        // derivatives -> K_MIN_NUM_YEAR_TERMS=3

        int itsCompoundMeth;  // Compounding method of the zero yields

        int itsInterpolMeth;  // Interpolation method of the zero yields :
                              // linear = 1 or continuous = 0 (see Summit doc)

        int itsLastBucketInt; // Interpolation method used for
                              // yields beyond the last yearTerm
                              // 0 for flat, 1 if the last line is used  

        // Methods

        virtual double DiscountFunction(double);
        virtual double D1DiscountFunction(double);
        virtual double D2DiscountFunction(double);

    public:

        ARM_ZeroLInterpol(void):ARM_ZeroCurve()
        {
            Init();
        }

        ARM_ZeroLInterpol(ARM_Date asofdate,
                          ARM_ZeroCurve* ZCSpread, ARM_ZeroCurve* zcInit,
                          int MMFreq, 
                          int SwapFreq,
                          ARM_Currency* ccy = ARM_DEFAULT_CURRENCY);

        ARM_ZeroLInterpol(ARM_Date& asOf, 
                          ARM_Vector* yearTerms, 
                          ARM_Vector* yields, 
                          int compMeth = K_COMP_CONT,
                          int lastBucketInt = 0, 
                          int interpolMeth = K_LINEAR, 
                          ARM_Currency* ccy = ARM_DEFAULT_CURRENCY);

        ARM_ZeroLInterpol(ARM_Date& asOf, 
                          ARM_Vector* yields,
                          ARM_CRV_TERMS& terms,
                          int compMeth = K_COMP_CONT,
                          int lastBucketInt = 0, 
                          int interpolMeth = K_LINEAR, 
                          ARM_Currency* ccy = ARM_DEFAULT_CURRENCY);

        ARM_ZeroLInterpol(ARM_Date& asOf, 
                          ARM_CRV_TERMS& terms, 
                          ARM_Vector* mktData, 
                          int MMVsFut = 0, 
                          int SwapVsFut = 0,
                          int raw = K_PAR, 
                          int Cont_lin = K_LINEAR,
                          int lastBucketInt = 0, 
                          ARM_Currency* ccy = ARM_DEFAULT_CURRENCY,
						  int swapFrqId = K_DEF_FREQ,
						  int fixDayCount = KNOBASE);

       ARM_ZeroLInterpol(ARM_CRV_TERMS& terms,
                         ARM_Date& asOf,
                         ARM_Vector* mktData,
                         int MMVsFut = 0,
                         int SwapVsFut = 0,
                         int raw = K_PAR, 
                         int Cont_lin = K_LINEAR,
                         int lastBucketInt = 0,
                         ARM_Currency* ccy = ARM_DEFAULT_CURRENCY);

        ARM_ZeroLInterpol(ARM_Date& asOf, 
                          ARM_CRV_TERMS& terms, 
                          ARM_Vector* mktData, 
                          ARM_Container* bonds, 
                          ARM_Vector* yields,
                          int MMVsFut, 
                          int lastBucketInt = 0, 
                          ARM_Currency* ccy = ARM_DEFAULT_CURRENCY);
 

        ARM_ZeroLInterpol(ARM_Date& asOf,
                          char* Terms[ARM_NB_TERMS],
                          ARM_Vector* mktData,
                          double mean_rates,
                          int raw,
                          int Cont_Lin,
                          int lastBucketInt,
                          ARM_Currency* ccy);


        /*---------------- TOY METHODS ------------------*/

        ARM_ZeroLInterpol(ARM_Date& asOf,
                          ARM_CRV_TERMS& terms,
                          ARM_Vector* mktData,
                          int MMVsFut = 0,
                          int SwapVsFut = 0,
                          int raw = K_PAR,
                          int Cont_lin = K_LINEAR,
                          int lastBucketInt = 0,
						  int Frequency = 0,
                          ARM_Currency* ccy = ARM_DEFAULT_CURRENCY);

        /*------------- END OF TOY METHODS --------------*/


		
		
		ARM_ZeroLInterpol(const ARM_ZeroLInterpol &);

        ARM_ZeroLInterpol(ARM_ZeroCurve &); 
        
       ~ARM_ZeroLInterpol(void)
        {
        }

        ARM_ZeroLInterpol& operator = (const ARM_ZeroLInterpol &);

		ARM_ZeroCurve* GenerateShiftCurve(ARM_CRV_TERMS& Term, ARM_Vector* epsilon);
        
		void BitwiseCopy(const ARM_Object* srczint)
        {
            ARM_ZeroLInterpol* zint = (ARM_ZeroLInterpol*) srczint;

            itsCompoundMeth = zint->itsCompoundMeth;

            itsInterpolMeth = zint->itsInterpolMeth;

            itsLastBucketInt = zint->itsLastBucketInt;
        }


        void Copy(const ARM_Object* srczint)
        {
            ARM_ZeroCurve::Copy(srczint);

            BitwiseCopy(srczint);
        }


        ARM_Object* Clone(void)
        {
            ARM_ZeroLInterpol* theClone = new ARM_ZeroLInterpol();


            theClone->Copy(this);

            return(theClone);
        }
            

        void Init(void)
        {
            SetName(ARM_ZERO_LIN_INTERPOL);

            itsLastBucketInt = 0;

            itsInterpolMeth = K_LINEAR;

            itsCompoundMeth = K_COMP_CONT;
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

            if (GetBPShifts())
                delete GetBPShifts();


            SetBPShifts(BPShifts);


            this->SetAsOfDate((ARM_Date) date);

            this->SetCompoundMeth(meth);

            GenerateDateTerms();
        }


        void Set(ARM_Vector* terms, ARM_Vector* rates, 
                 ARM_Date& date, int meth = 0)
        {
            ARM_ZeroCurve::Set(terms, rates, date);
 
            SetCompoundMeth(meth);

//            GenerateDateTerms(); (enleve, remettre si pb)

// visiblement, ne sert a rien, remettre si pb
//            if (meth == K_CONTINUOUS)
//                GenerateDiscountFactors(meth);

        }

        inline void SetLastBucketInt(int m)
        {
            itsLastBucketInt = m;
        }
 
        // Versions non virtuelles plus rapides
        double nvDiscountPrice(ARM_Date& maturity);

        double nvDiscountPrice(double maturityJul);

		double ComputeSlope(double yearTerm);
};




#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
