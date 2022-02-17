/*!
 *
 * Copyright (c) NATIXIS May 2007 Paris
 *
 *	\file infleg.cpp
 *
 *	\author  Francois Poitou
 *	\version 1.0
 *	\date May 2007
 */


#ifndef _INGPINFLATION_INFLEG2_H
#define _INGPINFLATION_INFLEG2_H

/// kernel
#include "inst/swapleg.h"
/// gpinflation
#include "gpinflation/infFloatingCouponVector.h"
#include "gpinflation/instrument.h"
#include "gpinflation/infPayOff.h"

#include "gpbase/gpvector.h"

#include "gpinflation/infComposedLeg.h"


CC_BEGIN_NAMESPACE( ARM )


class InfLeg  : public ARM_SwapLeg , public Instrument<InfFloatingCoupon<YOYPayOff> >, public InfLegBase
{
	public :

		InfLeg()
		{
			Init();
		}

		InfLeg(	const ARM_GP_T_Vector<ARM_CountedPtr<CashFlow> >& vCahFlow,
				//for backward compatibility : to remove
				int interpType,	double firstReset);

		InfLeg( const InfLeg& infLeg) ;
		InfLeg& operator=( const InfLeg& infLeg );

		void Init(void)
		{ 
			itsExpectedFwdRates = NULL;
		    itsDemCPIs          = NULL;
		    itsNumCPIs          = NULL;
		    itsDiscountFactors  = NULL;
		    itsResetTimes       = NULL;
		    itsPaymentTimes     = NULL;
		    itsAccrualTimes     = NULL;
		    itsNumJulianDates   = NULL;
		    itsDemJulianDates   = NULL;
		}

		virtual ~InfLeg()
		{
			if (itsExpectedFwdRates ) 
				delete itsExpectedFwdRates ;
			if (itsDemCPIs ) 
				delete itsDemCPIs ;
			if (itsNumCPIs ) 
				delete itsNumCPIs ;
			if (itsDiscountFactors ) 
				delete itsDiscountFactors ;
			if (itsResetTimes ) 
				delete itsResetTimes ;
			if (itsPaymentTimes ) 
				delete itsPaymentTimes ;
			if (itsAccrualTimes ) 
				delete itsAccrualTimes ;
			if (itsNumJulianDates ) 
				delete itsNumJulianDates ;
			if (itsDemJulianDates ) 
				delete itsDemJulianDates ;
		}

		/* Inspector method : calculate needed*/
		ARM_GP_Vector* GetDiscountFactors()  
		{	
			calculate() ;
			return itsDiscountFactors ;
		}

		ARM_GP_Vector* GetDemCPIs()  
		{
			calculate() ;
			return itsDemCPIs ;
		}

		ARM_GP_Vector* GetNumCPIs() 
		{
			calculate() ;
			return itsNumCPIs ;
		}

		ARM_GP_Vector* GetExpectedFwdRates() 
		{
			calculate() ;
			return itsExpectedFwdRates ;
		}


		/* const Inspector method */
		string							GetIndexName()		const ;
		ARM_CountedPtr<InfPerformance>	GetIndex()			const ; 
		double							GetFirstReset()		const {return itsFirstReset;}
		int								GetInterpType()		const {return itsInterpType;}
		ARM_GP_Vector*					GetNumJulianDates() const ;
		ARM_GP_Vector*					GetDemJulianDates() const ;

		virtual void performCalculation() ;
		virtual double ComputePrice(int mode) ;

	protected :
		
		ARM_GP_Vector* itsExpectedFwdRates ;
		ARM_GP_Vector* itsDemCPIs ;
		ARM_GP_Vector* itsNumCPIs ;
		ARM_GP_Vector* itsDiscountFactors;
		ARM_GP_Vector* itsResetTimes;
		ARM_GP_Vector* itsPaymentTimes;
		ARM_GP_Vector* itsAccrualTimes;
		ARM_GP_Vector* itsNumJulianDates;
		ARM_GP_Vector* itsDemJulianDates ;

		int itsInterpType;
		double itsFirstReset ;
		void CptExpectedFwdRates() ;
		void CptCashFlowValues() ;
		//virtual void SetModel(ARM_CountedPtr<ARM_Model> model) ;



};


CC_END_NAMESPACE()

#endif
/*--------------------------------------------------------------------------*/
/*---- End of file ----*/


