/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file discretisationscheme.h
 *  \brief general class for discretisation method
 *	\author  E Benhamou
 *	\version 1.0
 *	\date January 2004
 */

#ifndef _INGPINFRA_DISCRETISATIONSCHEME_H
#define _INGPINFRA_DISCRETISATIONSCHEME_H

#include "gpbase/port.h"
#include "gpbase/gplinalgtypedef.h"
#include "typedef.h"

const int TIMEINSERTER_LAG_THRESHOLD = 7;


CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_PricingModel;

std::vector<double>& MergeARM_Vector( const std::vector<double>& vec1, const std::vector<double>& vec2 );

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
/// \struct ABSTRACT STRUCT for discretisation
/// discretisation scheme is responsible for creating a model time schedule from timeInfos
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
struct ARM_DiscretisationScheme
{
	virtual std::vector<double>* ModelTimesFromTimesInfo(const ARM_TimeInfoPtrVector& timeInfos, ARM_PricingModel& pricingModel) const = 0;
};

/// event time only inserts eventtime from time infos
struct ARM_EventTime: ARM_DiscretisationScheme
{
	virtual std::vector<double>* ModelTimesFromTimesInfo(const ARM_TimeInfoPtrVector& timeInfos, ARM_PricingModel& pricingModel) const;
};

/// eventTime and model fixing insert event and model time
struct ARM_EventAndModelTime: ARM_DiscretisationScheme
{
	virtual std::vector<double>* ModelTimesFromTimesInfo(const ARM_TimeInfoPtrVector& timeInfos, ARM_PricingModel& pricingModel) const;
};


///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
/// \struct timeinsert: ABSTRACT STRUCT to insert time
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
struct ARM_TimeInserter
{
	virtual ARM_GP_Vector& AddTimes(const ARM_GP_Vector& initialTimes) const = 0;
};

/// \struct timeinsert: abstract struct to insert fix step time between the initial times
class ARM_FixStepTimeInserter
{
public:
	ARM_FixStepTimeInserter(int fixStep,int sensitivity = TIMEINSERTER_LAG_THRESHOLD ):itsFixStep(fixStep), itsSensitivity(sensitivity) {};
	virtual ARM_GP_Vector& AddTimes(const ARM_GP_Vector& initialTimes) const;
private:
	int itsFixStep;
	int itsSensitivity;
};


/// \struct timeinsert: struct to insert fix nb of points between the initial times
struct ARM_FixPointNbTimeInserter
{
	ARM_FixPointNbTimeInserter(int fixPointNb):itsFixPointNb(fixPointNb) {};
	virtual ARM_GP_Vector& AddTimes(const ARM_GP_Vector& initialTimes) const;
private:
	int itsFixPointNb;
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
