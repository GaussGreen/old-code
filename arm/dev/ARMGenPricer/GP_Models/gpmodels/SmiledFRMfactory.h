/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ModelparamsSFRMFactory.h
 *
 *  \brief class to control the creation of model params of SFRM!
 *	\author  E Benhamou
 *	\version 1.0
 *	\date January 2004
 */


#ifndef _INGPMODELS_SSFRMFACTORY_H
#define _INGPMODELS_SSFRMFACTORY_H

#include "gpbase/port.h"
#include "SmiledMM.h"


CC_BEGIN_NAMESPACE( ARM )


/// forward declaration
template <typename T> class ARM_SingletonHolder;

///-----------------------------------------------------------------------------
/// \class ARM_SmiledFRMfactory
/// \brief Factory class to create model params for SFRM implemented 
/// as a singleton
///-----------------------------------------------------------------------------
struct ARM_SmiledFRMfactoryImp
{
public:
	enum CalibPattern
	{
		LIBOR = 0,
		CMS_OLD,
		CMS,
		VMS_OLD,
		VMS
	};
	
	ARM_SmiledMM* CreateSmiledMarketModel(
		CalibPattern pattern,const ARM_ZeroCurvePtr& zc, const ARM_ModelParams* params=NULL, size_t timeStepsNb=500,size_t gridSize=501,double stdDevNb=6,bool skipPDE=false, bool allowInterpol=false, ARM_ModelParamsSmiled::CalibProxy calibProxy=ARM_ModelParamsSmiled::LocalVolatility ) const;

	ARM_DateStrip* CreateSmiledMarketModelDateStrip(
		CalibPattern pattern,
		const ARM_Date& startDate,			/// startDate of the strip
		const ARM_Date& endDate,			/// end Date
		long resetFreq,			
		long indexFreq,			
		int indexType,			
		long resetTiming,
		long dayCount,
		const char* resetCalendar,
		long fwdRule,
		long intRule,
		long stubRule,
		long resetGap) const;

private:
	/// to forbid client from using it except for the singleton holder
	ARM_SmiledFRMfactoryImp() {};
	friend class ARM_SingletonHolder<ARM_SmiledFRMfactoryImp>;
};

extern ARM_SingletonHolder<ARM_SmiledFRMfactoryImp> ARM_SmiledFRMfactory;

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

