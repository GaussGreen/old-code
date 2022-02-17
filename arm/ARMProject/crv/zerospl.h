/*----------------------------------------------------------------------------*
   zeroint.h

    
    Header for the ARM_Zero_Sp_Sm class, a class for computing a ARM_ZeroCurve
    using linear interpolation of an input zero curve defined by a finite 
    set of yearterms and yields.
 
*----------------------------------------------------------------------------*/ 


#ifndef _ZEROSPL_H
#define _ZEROSPL_H




#include "dates.h"
#include "zerocurv.h"
#include "linalg.h"



#define    K_MIN_NUM_YEAR_TERMS 1



class ARM_Zero_Sp_Sm : public ARM_ZeroCurve 
{   
    private:

        // the input year terms and yields vector must have 3 elements
        // in order to compute 1st and 2nd 
        // derivatives -> K_MIN_NUM_YEAR_TERMS=3


		ARM_Matrix* itsCoeffSpline;  	//  Coeff de chaque polynomes
        
        int         itsCompoundMeth;//  compounding method of the zero yields

        //    methods

        double D1DiscountFunction(double);
        double D2DiscountFunction(double);
        double DiscountFunction(double);
		void   InitSansMin(void);
		void   InitAvecMin(void);
    public:

        
		ARM_Zero_Sp_Sm(void)
        {	
			Init();

        }

		void Init(void);

        ARM_Zero_Sp_Sm(ARM_Date& asOf,
						ARM_Vector *yearTerms, 
                        ARM_Vector *yields, 
						int compMeth=0,
						int lastBucketInt=0);

        ARM_Zero_Sp_Sm(ARM_Zero_Sp_Sm &);
    
       ~ARM_Zero_Sp_Sm(void)
        {
           if (itsCoeffSpline)
               delete itsCoeffSpline;

	   }

        ARM_Zero_Sp_Sm& operator = (ARM_Zero_Sp_Sm &);

        void BitwiseCopy(const ARM_Object* srczint)
        {
            ARM_Zero_Sp_Sm* zint = (ARM_Zero_Sp_Sm*) srczint;


            itsCompoundMeth = zint->itsCompoundMeth;

 			
			if (zint->itsCoeffSpline)
               itsCoeffSpline = (ARM_Matrix *) zint->itsCoeffSpline->Clone();
            else
               itsCoeffSpline = NULL;
        }

        void Copy(const ARM_Object* srczint)
        {
            ARM_ZeroCurve::Copy(srczint);

            BitwiseCopy(srczint);
        }


        ARM_Object* Clone(void)
        {
            ARM_Zero_Sp_Sm* theClone = new ARM_Zero_Sp_Sm();


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
        
/*      void Set(int sz,
					double* yt, 
					double* zy, 
					int meth, char* date) 
        {
            if (itsYearTerms)
               delete itsYearTerms;

            itsYearTerms = new ARM_Vector(sz, yt);

            if (itsZeroYields)
               delete itsZeroYields;

            itsZeroYields = new ARM_Vector(sz, zy);

			if (itsCoeffSpline)
				delete itsCoeffSpline;

			itsCoeffSpline = NULL;

            this->SetAsOfDate((ARM_Date) date);

            this->SetCompoundMeth(meth);
        }
 
    virtual int SetYearTerms(ARM_Vector *);
    virtual int SetZeroYields(ARM_Vector *);*/

	double ARM_Zero_Sp_Sm::SpotForwardRate(double yearTerm);

	double ARM_Zero_Sp_Sm::SpotForwardRate(ARM_Date &maturity);

};




#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
