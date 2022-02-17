/*
 * Copyright (c) IXIS CIB February 2007 Paris
 *


/*! \file EqModelParams.h
 *
 *  \brief 
 *
 *	\author  M. Bernardo
 *	\version 1.0
 *	\date February 2007
 *
 *	\class ARM_EqModelParams contains the parameters of HW model in the contex of Equity Model with Stochastic Volatility
 */


#ifndef _INGPMODELS_EQMODELPARAMS_H
#define _INGPMODELS_EQMODELPARAMS_H

/// gpbase
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/curve.h"

/// gpinfra
#include "gpinfra/modelparams.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/modelparamutil.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/typedef.h"

#include <math.h>

#define PER_ANNI 365

CC_BEGIN_NAMESPACE( ARM )

/****************************************************

	\class		ARM_EQHWSV_ModelParams
	\author		Mathieu Bernardo
	\version	1.0
	\date		February 2007

*****************************************************/
class ARM_EQHWSV_ModelParams:public ARM_ModelParams{

protected :
	mutable std::vector<double> itsSchedule;

public:

	ARM_EQHWSV_ModelParams(	const	ARM_ModelParamVector & params = ARM_ModelParamVector() );
	ARM_EQHWSV_ModelParams(	const	ARM_EQHWSV_ModelParams& rhs ):ARM_ModelParams(rhs){ itsSchedule= rhs.itsSchedule; }

	ARM_EQHWSV_ModelParams&	operator= (const ARM_EQHWSV_ModelParams& rhs);
	virtual ARM_Object*		Clone()	   const {	return new ARM_EQHWSV_ModelParams(*this);	}

	virtual string ExportShortName() const { return "LEQMP";}
	virtual ~ARM_EQHWSV_ModelParams(){};
 
	double					GetMeanReversion() const;
	double					GetVolMeanReversion() const;

	size_t					FactorCount() const{	return	4;}
	double					G_Function(double t,double T);

	/// Schedule
	const	std::vector<double> & GetSchedule() const								{ return itsSchedule;		}
	void					SetSchedule( const std::vector<double> & schedule )	{ itsSchedule= schedule;	}		
	const	std::vector<double> & CptSchedule( const double& t, const double & T) const;

	void					MergeSchedule();

	/// Calibration
    virtual void PreProcessing	(ARM_ModelFitter& modelFitter,		int factorNb=0) {};
    virtual void PostProcessing	(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,	int factorNb=0 ) {};
};

inline double ARM_EQHWSV_ModelParams::G_Function(double t,double T){ 
	double l = GetMeanReversion();	
	return (1-exp(-l*(T-t) ) )/l; 
};

inline double ARM_EQHWSV_ModelParams::GetMeanReversion() const {
	return GetModelParam( ARM_ModelParamType::MeanReversion).ToCurveModelParam().GetCurve()->GetOrdinates()[0];
};

inline double ARM_EQHWSV_ModelParams::GetVolMeanReversion() const { 
	return GetModelParam( ARM_ModelParamType::VolMeanReversion).ToCurveModelParam().GetCurve()->GetOrdinates()[0];
};

/****************************************************

	\class		ARM_EQHWSV_NumMethod
	\author		Mathieu Bernardo
	\version	1.0
	\date		February 2007

*****************************************************/
class ARM_EQHWSV_NumMethods:public ARM_RootObject{

public:
	ARM_EQHWSV_NumMethods( const double &, const double &, const int & );
	ARM_EQHWSV_NumMethods( const ARM_EQHWSV_NumMethods & );
	~ARM_EQHWSV_NumMethods(){}

	virtual ARM_Object* Clone() const { return new ARM_EQHWSV_NumMethods(*this); }
	ARM_EQHWSV_NumMethods&	operator= (const ARM_EQHWSV_NumMethods& rhs);

	virtual string ExportShortName() const { return "LEQNM";}
	virtual string toString(const string& indent="",const string& nextIndent="") const;

public:
	double			itsMaxDecay;
	double			itsImAxis;
	int				itsIntStep;
};




CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

