/*!
 * Copyright (c) NATIXIS  April 2007 Paris
 *
 *	\file performance.h
 *	\author  Francois Poitou
 */

#ifndef _PERF_H
#define _PERF_H

/// gpbase
#include <gpbase/countedptr.h>		//ARM_CountedPtr
/// kernel
#include "inst/irindex.h"			//ARM_IRIndex
#include "glob/dates.h"				//ARM_Date
///purple
#include "glob/period.h"

CC_BEGIN_NAMESPACE( ARM )


	template<class indexType> 
	class Performance : public ARM_IRIndex{

		public :

			Performance (ARM_CountedPtr<indexType> index,	const Period& tenor)
			: 	ARM_IRIndex( index.operator*() ),itsIndex(index), itsTenor(tenor){}

			Performance (const Performance& perf)
			: itsIndex(perf.itsIndex), itsTenor(perf.itsTenor){}

			virtual ~Performance (){}

			Period GetTenor()const {return itsTenor ;}
			ARM_CountedPtr<indexType> GetIndex()const {return itsIndex ;}

			virtual double fixing(const ARM_Date& fixingDate )const =0;

		protected :
			ARM_CountedPtr<indexType> itsIndex ;
			Period itsTenor ;

	};


CC_END_NAMESPACE()

#endif


