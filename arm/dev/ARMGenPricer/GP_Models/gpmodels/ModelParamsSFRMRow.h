/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ModelparamsSFRMRow.h
 *
 *  \brief class for ROW version of model params
 *
 *	\author  E Benhamou
 *	\version 1.0
 *	\date January 2004
 */


#ifndef _INGPMODELS_MODELPARAMSSFRMROW_H
#define _INGPMODELS_MODELPARAMSSFRMROW_H

#include "gpbase/port.h"
#include "ModelParamsSFRM.h"
#include "gpbase/curve.h"
#include "typedef.h"

CC_BEGIN_NAMESPACE( ARM )
//-----------------------------------------------------------------------------
// \class ARM_ModelParamsSFRMRow
// \brief Interface class for the row version of the SFRM model
//-----------------------------------------------------------------------------
class ARM_ModelParamsSFRMRow : public ARM_ModelParamsSFRM
{
private:
	double IntegratedLocalVarianceFromZero(double t) const;
	ARM_CurvePtr itsSquaredIntegratedVol;

public:
	ARM_ModelParamsSFRMRow( const ARM_ModelParamVector& params, ARM_IRIndex* index, size_t factorsNb = 1 );
	ARM_ModelParamsSFRMRow( const ARM_ModelParamsSFRMRow& rhs );
	ARM_ModelParamsSFRMRow& operator=( const ARM_ModelParamsSFRMRow& rhs );
	virtual ~ARM_ModelParamsSFRMRow();
		
	/// pricing function
    virtual double VolatilityFunction(double t, double T) const;
	virtual double IntegratedLocalVariance(double s, double t) const;
	virtual double MaturityTerm(double T) const;
	virtual double MaturityTermSquared(double T) const;
	virtual double VolatilitySpotSquared(double s) const;

    virtual double VarianceToTime(double var,double minTime,double maxTime) const;

    /// because the SFRM Row caches in the squared integrated vol
    /// whenever we change a model param... we update the squaredIntegratedVol
    virtual iterator SetModelParamValue( int paramType, size_t i, double value, double time, double tenor = 0.0 );
    void InitSquaredIntegratedVol();
    void PreComputeSquaredIntegratedVol();

	virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model, int FactorNb = 0 );
	
	/// Standard ARM object support
	virtual ARM_Object* Clone() const;
	virtual string toString(const string& indent="", const string& nextIndent="") const;
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

