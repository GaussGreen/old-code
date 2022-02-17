/*!
 *
 * Copyright (c) NATIXIS May 2007 Paris
 *
 *	\file instrument.cpp
 *
 *	\author  Francois Poitou
 *	\version 1.0
 *	\date May 2007
 */


#ifndef _INSTRUMENT_H
#define _INSTRUMENT_H

/// gpbase

/// kernel

/// gpinflation
#include "gpinflation/floatingCoupon.h"


CC_BEGIN_NAMESPACE( ARM )

class SimpleLazyObject /*: public ARM_Object */{
	public :
		SimpleLazyObject() : itsCalculated(false) {} ;
		SimpleLazyObject(const SimpleLazyObject& sLO) : itsCalculated(sLO.itsCalculated) {} ;
		virtual ~SimpleLazyObject(){} ;
		virtual void performCalculation() = 0 ;
		void calculate()
		{
			if (!itsCalculated )
			{
				performCalculation();
				itsCalculated = true ;
			}
		}

	private :
		double itsCalculated ;
} ;


template<class CashFlowType>
class Instrument : public SimpleLazyObject {
	public :
		Instrument(){}
		Instrument(	const ARM_GP_T_Vector<ARM_CountedPtr<CashFlow> >& leg)
			: itsLeg(leg.size())
		{	
			for (int i=0; i<leg.size(); ++i) {
				const CashFlowType& cashflow = dynamic_cast<const CashFlowType&>(*(leg[i]));
				itsLeg[i] = ARM_CountedPtr<CashFlowType>( new CashFlowType(cashflow));
			}
		} 
		Instrument(const Instrument& ins):itsLeg(ins.itsLeg){}
		virtual ~Instrument()
		{
			double z = 0. ;
		}
		/*
		void ComputePrice(const ARM_CountedPtr<ARM_Model>& model)
		{
			SetModel(model) ;
			ComputePrice() ;
		}
		*/
		ARM_GP_T_Vector<ARM_CountedPtr<CashFlow> > GetLeg() const { return itsLeg ; }

	private :
		/*
		virtual void SetModel(ARM_CountedPtr<ARM_Model> model) = 0 ;
		virtual void ComputePrice() 
		{	
			ARMTHROW(ERR_INVALID_ARGUMENT," Not Implemented");
		}
		*/
	protected :
		/*
		void  SetModel(ARM_CountedPtr<ARM_Model> model)			
		{
			itsInstrumentModel = model ;//CreateClonedPtr(&*model);
		}
		*/
		ARM_GP_T_Vector< ARM_CountedPtr<CashFlow> > itsLeg ;
		//ARM_CountedPtr<ARM_Model> itsInstrumentModel ;
} ;


CC_END_NAMESPACE()

#endif
/*--------------------------------------------------------------------------*/
/*---- End of file ----*/


