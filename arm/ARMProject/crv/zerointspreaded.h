/*
 * $Log: zerointspreaded.h,v $
 * Revision 1.0  2006/08/09 16:33:00 mb
 * "Constification"
 *
 */


/*----------------------------------------------------------------------------*
    zerointspreaded.h

    
    Header for the ARM_ZeroLInterpolSpread class, a class which incorpores the
	Basis Market Data
 
*----------------------------------------------------------------------------*/ 
#ifndef _ZEROINTSPREADED_H
#define _ZEROINTSPREADED_H


#include "dates.h"
#include "zeroint.h"
#include "linalg.h"
#include "matlib.h"




class ARM_BasisCurve : public ARM_ZeroLInterpol 
{
    private:

        ARM_ZeroCurve* itsBasisCurve; 
		ARM_ZeroCurve* itsZeroCurve;

        int            itsMMFreq;

        int            itsSwapFreq;

    public:

        ARM_BasisCurve(void):ARM_ZeroLInterpol()
        {  
            Init();  
        }

        ARM_BasisCurve(ARM_Date		  asofdate,
					   ARM_ZeroCurve* ZCSpread, 
					   ARM_ZeroCurve* ZCInit,
					   int			  MMFreq, 
					   int			  SwapFreq,
					   ARM_Currency*  ccy = ARM_DEFAULT_CURRENCY)
                       :
		               ARM_ZeroLInterpol(asofdate,
						                 ZCSpread, 
						                 ZCInit,
						                 MMFreq, 
						                 SwapFreq,
						                 ccy)
        {	
            Init();
			
            itsBasisCurve = dynamic_cast<ARM_ZeroCurve* > (ZCSpread ->Clone());
											
            itsZeroCurve  = dynamic_cast<ARM_ZeroCurve* > (ZCInit->Clone());
			
			SetCurrencyUnit( itsBasisCurve->GetCurrencyUnit() );

            itsMMFreq     =  MMFreq;

			itsSwapFreq   =  SwapFreq;

			if (ZCSpread->GetMktData())
				SetMktData(static_cast< ARM_MarketData*>(ZCSpread->GetMktData()->Clone()));
		}

		ARM_BasisCurve(const ARM_BasisCurve & basisCurve):	ARM_ZeroLInterpol(basisCurve)
        {
			Init();
				
            BitwiseCopy(& basisCurve);
		}

		// its not a copy constructor. We just need to copy the element of zeroLInterpol
		ARM_BasisCurve(ARM_ZeroLInterpol & zeroLInterpol):ARM_ZeroLInterpol(zeroLInterpol)
		{
			Init();

			BitwiseCopy(& zeroLInterpol);
		} 
        
       ~ARM_BasisCurve(void)	
        { 
			if( itsBasisCurve) { delete itsBasisCurve;	itsBasisCurve	= NULL; }	
			if( itsZeroCurve)  { delete itsZeroCurve;	itsZeroCurve	= NULL; }		
		}

		ARM_BasisCurve &	operator =				( const ARM_BasisCurve & zeroLInterpolSpreaded);
		
		void				BitwiseCopy				( const ARM_Object* srczintsp );

        void				Copy					( const ARM_Object* srczintsp );

        ARM_Object*			Clone					( );

		void				Init					( );

		void				GenerateBasisCurve		( ARM_ZeroCurve* zcurve);

		ARM_ZeroCurve*		GetBasisCurve			( ){	return itsBasisCurve;   }

		void				GenerateZeroCurve		( ARM_ZeroCurve* zcurve);

		ARM_ZeroCurve*		GetZeroCurve			( ){	return itsZeroCurve;   }

		ARM_BasisCurve* GenerateShiftCurve(ARM_CRV_TERMS& Term, ARM_Vector* epsilon, const char* curveToBump);
		
		ARM_BasisCurve* GenerateShiftCurveFwd(ARM_CRV_TERMS& Term, ARM_Vector* epsilon, const char* curveToBump);

  

        void View(char* id = NULL, FILE* ficOut = NULL);
};

#endif

/*----------------------------------------------------------------------------------------*/
/*---- End Of File ----*/
