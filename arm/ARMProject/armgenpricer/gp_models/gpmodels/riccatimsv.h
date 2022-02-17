/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *  function solution of the Riccati Equation with Piece-Wise functions A(t), B(t) and C(t)
 *  for MSV we have a second equation 
 *	df/dt = A(t) * f(t)² + B(t) * f(t) + C(t)
 *  dg/dt = D(t) * f(t)
 *	\file riccati.h
 *
 *  \brief
 *
 *	\author  A. Triki
 *	\version 1.0
 *	\date October 2005
 */

 
#ifndef _INGPMODELS_RICCATIMSV_H
#define _INGPMODELS_RICCATIMSV_H

#include <glob/firsttoinc.h>

/// gpbase
#include "gpbase/countedptr.h"
#include "gpbase/functor.h"
#include "gpbase/typedef.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/curve.h"

/// gpclosedforms
#include "gpclosedforms/riccati.h"

/// gpinfra
#include "gpinfra/pricingstates.h"

/// gpmodels
#include "gpmodels/typedef.h"
#include "gpmodels/msv.h"

#include <vector>
CC_USING_NS(std,vector)

#include <map>
CC_USING_NS(std,map)



CC_BEGIN_NAMESPACE(ARM)


class ARM_RiccatiMSV:public ARM_Riccati
{
	private:
		//// Swaptions Details for Riccati Calculation

		ARM_SVModels* itsModel;
		string itsCurveName;
		double itsFloatStartTime;
		double itsFloatEndTime;
		std::vector<double> itsFixPayTimes;
		std::vector<double> itsFixPayPeriods;
		ARM_PricingStatesPtr itsStates;
		double itsCstShiftValue;
		double itsTempDt;
		double itsMu;

	public :

	/// Class to store Additional Datas
	struct ARM_FunctDataMSV
    {
		double				itsDt; 
		double				itsPartIntegral;
    };


    /// Internal structures for function description
	struct ARM_Funct0MSV : public ARM_Riccati::ARM_Funct0 
    {
		ARM_Funct0MSV(const ARM_FunctData* data=NULL) : ARM_Funct0(data) {}
        virtual double ComputeB(double t1, double t2) const;
    };

	struct ARM_Funct1MSV : public ARM_Riccati::ARM_Funct1 
    {
		ARM_Funct1MSV(const ARM_FunctData* data=NULL) : ARM_Riccati::ARM_Funct1(data) {}
        /// Discriminant = 0 : one solution
        virtual double ComputeB(double t1, double t2) const;    
	};

	struct ARM_Funct2MSV : public ARM_Riccati::ARM_Funct2
    {
		ARM_Funct2MSV(const ARM_FunctData* data=NULL) : ARM_Riccati::ARM_Funct2(data) {}
        /// Discriminant > 0 : two solutions
        virtual double ComputeB(double t1, double t2) const;
    };

	struct ARM_Funct2OutMSV : public ARM_Riccati::ARM_Funct2Out
    {
		ARM_Funct2OutMSV(const ARM_FunctData* data=NULL) : ARM_Riccati::ARM_Funct2Out(data) {}
        /// Solution outside roots
        virtual double ComputeB(double t1, double t2) const;
    };

	struct ARM_Funct2InMSV : public ARM_Riccati::ARM_Funct2In
    {
		ARM_Funct2InMSV(const ARM_FunctData* data=NULL) : ARM_Riccati::ARM_Funct2In(data) {}
        /// Solution inside roots
        virtual double ComputeB(double t1, double t2) const;
    };

	struct ARM_Funct2SupMSV : public ARM_Riccati::ARM_Funct2Sup
    {
		ARM_Funct2SupMSV(const ARM_FunctData* data=NULL) : ARM_Riccati::ARM_Funct2Sup(data) {}
        /// Solution equal to the greater root
        virtual double ComputeB(double t1, double t2) const;
    };

	struct ARM_Funct2InfMSV : public ARM_Riccati::ARM_Funct2Inf
    {
		ARM_Funct2InfMSV(const ARM_FunctData* data=NULL) : ARM_Riccati::ARM_Funct2Inf(data) {}
        /// Solution equal to the lower root
        virtual double ComputeB(double t1, double t2) const;
    };

	ARM_Curve* itsDtCurve;
	/// Function parameters for all sampling intervals
    vector< ARM_FunctDataMSV > itsFunctionDatasMSV;


	/// Standard ARM object support
	ARM_RiccatiMSV(const ARM_Curve& At,const ARM_Curve& Bt,const ARM_Curve& Ct,const ARM_Curve& Dt);
	ARM_RiccatiMSV(const ARM_Curve& At,const ARM_Curve& Bt,const ARM_Curve& Ct,const ARM_Curve& Dt,
		ARM_SVModels* model,
		string curveName, 
		double floatStartTime,
		double floatEndTime,
		std::vector<double> fixPayTimes,
		std::vector<double> fixPayPeriods,
		ARM_PricingStatesPtr states,
		double cstShiftValue,
		double mu);

	ARM_RiccatiMSV(const ARM_RiccatiMSV& rhs);
	virtual ~ARM_RiccatiMSV();
    ARM_RiccatiMSV& operator = (const ARM_RiccatiMSV& rhs);
	void CopyNoCleanUp(const ARM_RiccatiMSV& rhs);
	virtual ARM_Object* Clone() const;
	
	/// Store Additional Datas
	virtual void ComputeAdditionalDatas();

	/// Accesors
	inline void setMu (double mu) { itsMu = mu;}

	void ComputeRoots();
    /// Collect & save model parameters for function definition
    virtual void InitFunctionDatas();
	void UpdateCstVolatility(double vol);


	virtual void UpdateTempParams(double x) const;
	void derivs(double x, std::vector<double>& yt, std::vector<double>& dyt) const;


	ARM_Riccati::ARM_FunctionsMapIter GenerateFunction(double maturity);
	double ARM_RiccatiMSV::B0T(double T);

};

class ConstantVolatilityFinder : public ARM_GP::UnaryFunc<double,double> 
{
public: 
		ConstantVolatilityFinder(
			ARM_RiccatiMSV* riccatiFunction,
			double time,
			double maturity,
			double target,
			double vol_Zero) :
		itsRiccatiFunction(riccatiFunction),
		itsTime(time),
		itsMaturity(maturity),
		itsTarget(target),
		itsVol_Zero(vol_Zero)
		{
		};

		virtual double operator() (double vol) const 
		{
			//// To be optimized when we know the idex of the ordinate in the Curve
	//		itsRiccatiFunction->itsCtCurve->insert(itsTime,vol);
			//// To be optimized
			itsRiccatiFunction->UpdateCstVolatility(vol);
			double A0T = itsRiccatiFunction->A(0,itsMaturity);
			double B0T = itsRiccatiFunction->B0T(itsMaturity);
			return (exp(A0T - itsVol_Zero * B0T) - itsTarget);

		}
private:
	ARM_RiccatiMSV* itsRiccatiFunction;
	double itsTime;
	double itsMaturity;
	double itsTarget;
	double itsVol_Zero;
};



CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/