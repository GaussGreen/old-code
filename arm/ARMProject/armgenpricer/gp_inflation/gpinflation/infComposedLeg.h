/*!
 *
 * Copyright (c) NATIXIS May 2007 Paris
 *
 *	\file infComposedLeg.hpp
 *
 *	\author  Francois Poitou
 *	\version 1.0
 *	\date May 2007
 */


#ifndef _INFCOMPOSEDLEG_H
#define _INFCOMPOSEDLEG_H

/// kernel
#include "inst/swapleg.h"
/// gpinflation
#include "gpinflation/infFloatingCouponVector.h"
#include "gpinflation/instrument.h"
#include "gpinflation/infPayOff.h"


CC_BEGIN_NAMESPACE( ARM )

/*
Pure virtual class for interface
*/
class InfLegBase  
{
	public :

		InfLegBase(){}
		virtual ~InfLegBase(){} 

		virtual ARM_GP_Vector* GetDiscountFactors()		= 0 ;
		virtual ARM_GP_Vector* GetDemCPIs()				= 0 ;
		virtual ARM_GP_Vector* GetNumCPIs()				= 0 ;
		virtual ARM_GP_Vector* GetExpectedFwdRates()	= 0 ;

		virtual string							GetIndexName()		const   = 0 ;
		virtual ARM_CountedPtr<InfPerformance>	GetIndex()			const   = 0 ; 
		virtual double							GetFirstReset()		const   = 0 ;
		virtual int								GetInterpType()		const   = 0 ;
		virtual ARM_GP_Vector*					GetNumJulianDates() const   = 0 ;
		virtual ARM_GP_Vector*					GetDemJulianDates() const   = 0 ;
};

class InfComposedLeg  : public Instrument<InfComposedCoupon<YOYPayOff> >, public ARM_SwapLeg, public InfLegBase
{
	public :

		InfComposedLeg()
		{
			Init();
		}

		InfComposedLeg(	const ARM_GP_T_Vector<ARM_CountedPtr<CashFlow> >& vCahFlow,
				//for backward compatibility : to remove
				int interpType,	double firstReset);

		InfComposedLeg( const InfComposedLeg& infLeg) ;
		InfComposedLeg& operator=( const InfComposedLeg& InfComposedLeg );

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

		virtual ~InfComposedLeg(){
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
		ARM_GP_Vector* GetDiscountFactors()  {	
			calculate() ;
			return itsDiscountFactors ;
		}

		ARM_GP_Vector* GetDemCPIs()  {
			calculate() ;
			return itsDemCPIs ;
		}

		ARM_GP_Vector* GetNumCPIs() {
			calculate() ;
			return itsNumCPIs ;
		}
		ARM_GP_Vector* GetExpectedFwdRates() {
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

};



CC_END_NAMESPACE()

#endif
/*--------------------------------------------------------------------------*/
/*---- End of file ----*/


