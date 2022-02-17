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

class SimpleLazyObject {
	public :
		SimpleLazyObject() : itsCalculated(false) {} ;
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




class Instrument : public SimpleLazyObject {
	public :
		Instrument(	ARM_GP_T_Vector<ARM_CountedPtr<CashFlow> > leg):
					 itsLeg(leg){} ;
		Instrument(const Instrument& ins):itsLeg(ins.itsLeg){}
		virtual ~Instrument(){} ;
		void ComputePrice(const ARM_CountedPtr<ARM_Model>& model)
		{
			SetModel(model) ;
			ComputePrice() ;
		}
		ARM_GP_T_Vector<ARM_CountedPtr<CashFlow> > GetLeg() const { return itsLeg ; }

	private :
		virtual void SetModel(ARM_CountedPtr<ARM_Model> model) = 0 ;
		virtual void ComputePrice() = 0 ;
	protected :
			ARM_GP_T_Vector<ARM_CountedPtr<CashFlow> > itsLeg ;
			ARM_CountedPtr<ARM_Model> itsInstrumentModel ;
} ;


CC_END_NAMESPACE()

#endif
/*--------------------------------------------------------------------------*/
/*---- End of file ----*/


