/*!
 * Copyright (c) NATIXIS  April 2007 Paris
 *
 *	\file infPerformance.h
 *	\author  Francois Poitou
 */

#ifndef _INFPERF_H
#define _INFPERF_H

/// gpbase
#include <gpbase/countedptr.h> 				//ARM_CountedPtr
#include "gpbase/argconvdefault.h"
/// gpinflation
#include "gpinflation/infidx.h" 			//ARM_InfIdx
#include "gpinflation/performance.h"
/// kernel
#include "glob/dates.h"			//ARM_Date


CC_BEGIN_NAMESPACE( ARM )

	class InfPerformance : public Performance<ARM_InfIdx>{

		public:
			InfPerformance (	ARM_CountedPtr<ARM_InfIdx> index,
								const Period& tenor)
			: Performance<ARM_InfIdx>(index, tenor){
				//for backward compatibility : to remove
				SetDayCount(ARM_ArgConv_LgNameDayCount.GetNumber("30/360")) ;
			}


			InfPerformance (const InfPerformance& perf)
			: Performance<ARM_InfIdx>(perf.itsIndex,perf.itsTenor){}

			virtual ~InfPerformance (){}

			double fixing(const ARM_Date& fixingDate )const ;

		private : 
	};

CC_END_NAMESPACE()

#endif


