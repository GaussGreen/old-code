/*!
 * Copyright (c) NATIXIS  May 2007 Paris
 *
 *	\file infPerformance.cpp
 *	\author  Francois Poitou
 */

#include "gpinflation/infPerformance.h"
#include "mod\bssmiled.h"

/// gpinflation

/// gpbase

/// kernel

CC_BEGIN_NAMESPACE( ARM )

	double InfPerformance::fixing(const ARM_Date& fixingDate )const {
			double Ii					= itsIndex->FwdCPI( fixingDate );
			ARM_Date firstFixingDate(fixingDate) ;
			firstFixingDate.AddPeriodMult(itsTenor.GetUnit(), -itsTenor.GetLength(), itsIndex->GetCurrencyUnit()->GetCcyName() ) ;
			double Ij = itsIndex->FwdCPI( firstFixingDate );
			if (IS_ZERO(Ij)) ARM_THROW( ERR_INVALID_ARGUMENT, "Null denominator, ratio cannot be calculated");

			double perf = Ii/Ij-1 ;
			return perf ;
	}


CC_END_NAMESPACE()
