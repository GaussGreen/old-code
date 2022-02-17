#ifndef _ICM_VOLINT_H
#define _ICM_VOLINT_H

/*********************************************************************************/
/*! \class  ICM_VolInterpol icm_volinterpol.h "icm_volinterpol.h"
 *  \author : Damien pouponneau
 *	\version 0.0
 *	\date   14 dec 2005
 *	\file   icm_volinterpol.h
 *	\brief  Vol curve which allows cubic spline interpolation 
 *	\Modified by : Fakher Ben Atig
/**********************************************************************************/


#include "ARMKernel\crv\volint.h"

class ICM_VolInterpol : public ARM_VolLInterpol 
{
    public:

		ICM_VolInterpol() {}
		ICM_VolInterpol(ARM_VolCurve* vol) 
		{Copy(vol);}
        
        // here m2 is not relevant
        double VolatilityFunction(double m1, double k, double m2)
		{
			double result=0.;

			switch (GetInterpType())
			{
			case K_SPLINE:
				{
					int sizeLin = GetVolatilities()->GetNumCols();
					ARM_Vector BaseCorrel(sizeLin);
					
					long i = indexBeforeValue(GetExpiryTerms(), m1);
					
					long iToUpdate;

					if ( i == GetExpiryTerms()->GetSize()-1 )
						iToUpdate = GetExpiryTerms()->GetSize()-1;
					else if (i == -1)
						iToUpdate = 0;
					else if ( fabs(m1 - GetExpiryTerms()->Elt(i))
						< fabs(GetExpiryTerms()->Elt(i+1)-m1)
						)
						iToUpdate = i;
					else
						iToUpdate = i+1;
					for (int j = 0; j < sizeLin; j++)
						BaseCorrel.Elt(j) = GetVolatilities()->Elt(iToUpdate,j);

    				ARM_Vector* Strikes = GetStrikes();
					
					if (Strikes->GetSize() <= 2 || BaseCorrel.GetSize() <= 2)
						throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						"For Spline Interpolation, Base Correlation Matrix should have more than 2 strikes!");
					
					result = SplineInterpolateFunc(Strikes,&BaseCorrel,k,NULL,1);
				}
				break;
			case K_ICM_STEPUP_RIGHT_MATURITY : //values of f(x) with x in [a;b[ are send on f(b) value
				{
					int sizeLin = GetVolatilities()->GetNumCols();
										
					long i = indexBeforeValue(GetExpiryTerms(), m1);
					long j = indexBeforeValue(GetStrikes(), k);
					ARM_Vector* Strikes = GetStrikes();

					if (j<0) j=0;
					
					if (i<GetExpiryTerms()->size()-1)
					{
						if (i<0) 
							{i=0;}
						else if ((m1>GetExpiryTerms()->Elt(i)) || 
								CHECK_EQUAL(GetExpiryTerms()->Elt(i),m1))
							{i++;}
					}
		
					result = GetVolatilities()->Elt(i,j);
				}
				break;
			case K_ICM_STEPUP_LEFT_MATURITY ://values of f(x) with x in ]a;b[ are send on f(b) value
				{							 // a->f(a)	b->f(b)
					int sizeLin = GetVolatilities()->GetNumCols();
										
					long i = indexBeforeValue(GetExpiryTerms(), m1);
					long j = indexBeforeValue(GetStrikes(), k);
					ARM_Vector* Strikes = GetStrikes();

					if (j<0) j=0;
					
					if (i<GetExpiryTerms()->size()-1)
					{
						if (i<0) 
							{i=0;}
						else if ((m1>GetExpiryTerms()->Elt(i)) && 
							!CHECK_EQUAL(GetExpiryTerms()->Elt(i),m1))
							{i++;}
					}
						
					result = GetVolatilities()->Elt(i,j);
				}
				break;
			case K_ICM_STEPUP_LEFT_LINEAR ://values of f(x) with x in ]a;b[ are send on f(b) value
				{							 // a->f(a)	b->f(b)
					int sizeLin = GetVolatilities()->GetNumCols();
										
					long i = indexBeforeValue(GetExpiryTerms(), m1);
					long j = indexBeforeValue(GetStrikes(), k);
					ARM_Vector* Strikes = GetStrikes();

					if (j<0) j=0;
					
					if (i<GetExpiryTerms()->size()-1)
					{
						if (i<0) 
							{i=0;}
						else 
						{return ARM_VolLInterpol::VolatilityFunction(m1, k, m2);}
					}
						
					result = GetVolatilities()->Elt(i,j);
				}
				break;
			default :
				result = ARM_VolLInterpol::VolatilityFunction(m1, k, m2);
			}

			return result;
		}


        void BitwiseCopy(const ARM_Object* srcVolCurve) {}
 
        void Copy(const ARM_Object* vCurve)
        {
            ARM_VolLInterpol::Copy(vCurve);
            BitwiseCopy(vCurve);
        }

        ARM_Object* Clone(void)
        {
            ICM_VolInterpol* theClone = new ICM_VolInterpol();
            theClone->Copy(this);

            return(theClone);
        }

};


#endif
