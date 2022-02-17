/*
 * $Log: volflat.h,v $
 * Revision 1.4  2003/09/19 15:15:00  jpriaudel
 * cumulative bump added
 *
 * Revision 1.3  2003/08/18 09:52:12  mab
 * improvement in the principal method :
 * virtual double VolatilityFunction(double m1, double K, double m2)
 * virtual double VolatilityFunction(double m1, double m2) :
 * No default parameter but 2 methods!
 *
 * Revision 1.2  2002/10/11 08:26:07  mab
 * Added : BumpVolatility
 *
 */


/*----------------------------------------------------------------------------*

     volflat.h
 
*----------------------------------------------------------------------------*/
#ifndef _VOLFLAT_H
#define _VOLFLAT_H



#include "volint.h"



class ARM_VolFlat : public ARM_VolLInterpol 
{
    private:

        double       itsVolatility; //    value of volatility

        // here m2 is not relevant
        double VolatilityFunction(double m1, double k, double m2)
        {
            return(itsVolatility);
        }

        double VolatilityFunction(double m1, double m2)
        {
            // We are in the case of a Matrix
            // the Third param is not relevant

            return(VolatilityFunction(m1, m2, 0.0));
        }

        void Init(void)
        {
            SetName(ARM_VOL_FLAT);

            itsVolatility = 0.0;
        }

    public:

        ARM_VolFlat(void)
        {
            Init();
        }

        virtual ~ARM_VolFlat(void){}


        ARM_VolFlat(ARM_Date& asOf, double vol,
                    ARM_Currency* ccy = ARM_DEFAULT_CURRENCY); 

virtual double CalcNumericalObjectSignature(void)
        {
             double signature = 0.0;

             signature += GetAsOfDate().GetJulian()*itsVolatility;

             return(signature);
        }

        void BitwiseCopy(const ARM_Object* srcZVf)
        {
            ARM_VolFlat* vf = (ARM_VolFlat *) srcZVf;
 
 
            itsVolatility = vf->itsVolatility;
        }
 
        void Copy(const ARM_Object* srcVf)
        {
            ARM_VolLInterpol::Copy(srcVf);
 
            BitwiseCopy(srcVf);
        }
 
        virtual ARM_Object* Clone(void)
        {
            ARM_VolFlat* theClone = new ARM_VolFlat();
 
 
            theClone->Copy(this);
 
            return(theClone);
        }

        void Set(ARM_Date& asOf, double vol); 


        ARM_VolFlat(const ARM_VolFlat& flatCrv) : ARM_VolLInterpol(flatCrv)
        {
            Init();

            BitwiseCopy(&flatCrv);
        }

    
        ARM_VolFlat& operator = (const ARM_VolFlat &flatCrv)
        {
            (*this).ARM_VolLInterpol::operator = (flatCrv);

            BitwiseCopy(&flatCrv);
    
            return(*this);
        }

        void View(char* id = NULL, FILE* fOut = NULL);

        double GetVolatility(void)
        {
            return(itsVolatility);
        }
 
        void SetVolatility(double flatVol)
        {
            itsVolatility = flatVol;
        }

        void BumpVolatility(double value, int nthLine = 0, int nthCol = 0,
                            int isCumul = K_NO, int isAbsolute = K_YES)
        {
            ARM_VolCurve::BumpVolatility(value, nthLine, nthCol, isCumul, isAbsolute);

            if (isAbsolute == K_YES)
				SetVolatility(GetVolatility()+value);
			else
				SetVolatility(GetVolatility()*(1.+value/100.));
        }
};




#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
