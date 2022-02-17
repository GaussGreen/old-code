/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 */

#ifndef _INGPMODELS_PROCESSBUILDERSMFRM_H
#define _INGPMODELS_PROCESSBUILDERSMFRM_H

/// gpbase
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/curve.h"
#include "gpbase/curvetypedef.h"
#include "gpbase/vectormanip.h"

/// gpinfra
#include "gpinfra/pricingstates.h"


CC_BEGIN_NAMESPACE( ARM )

class ARM_VanillaSecurityDensity;

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilder
////////////////////////////////////////////////////

class ARM_ProcessBuilder : public ARM_RootObject
{
protected: 

	inline double	MeanReversion() const			{ return itsMeanReversion; }

	inline double	Horizon() const					{ return itsHorizonDate; }

	inline double	Volatility( double t ) const	{ return itsVolCurve.Interpolate(t);}

	void			ComputeIntegratedVariance();

	double			IntegratedVarianceFromZero(double t) const;

	inline const	ARM_Curve&	GetVolCurve() const		{ return itsVolCurve;}

public:

	ARM_ProcessBuilder ( const ARM_ProcessBuilder& rhs );
	ARM_ProcessBuilder ( double mrs , const ARM_Curve& volCurve );
	virtual ~ARM_ProcessBuilder(){};

	virtual double	GetDelta() const	{ return itsDelta;};

	inline  double	GetFwdRate() const	{ return itsFwd;};

	virtual double	GetEqShift() const = 0;

	virtual double	GetEqVol() const = 0;

	virtual	ARM_GP_VectorPtr Rate (	double t,
									const ARM_PricingStatesPtr& states,
									size_t modelNb,
									size_t k , bool extrapol=false) const = 0;

	virtual	ARM_GP_VectorPtr DxRate (	double t,
										const ARM_PricingStatesPtr& states,
										size_t modelNb,
										size_t k , bool extrapol=false) const = 0;

	virtual	ARM_GP_VectorPtr DFRatio (	double t,
										const ARM_PricingStatesPtr& states,
										size_t modelNb,
										size_t k , bool extrapol=false) const = 0;

	virtual ARM_GP_VectorPtr DxLogDFRatio(	double t,
											const ARM_PricingStatesPtr& states,
											size_t modelNb,
											size_t k , bool extrapol=false) const = 0;

	virtual double	GetDiffusion() const = 0;

	virtual void	ComputeProxyForCalib(size_t calibProxy,const ARM_VanillaSecurityDensity& density,ARM_Curve* volatilityCalib,double strike=0) = 0;

	virtual void	Calibrate(	double horizonDate,
					const ARM_VanillaSecurityDensity& density,
					const std::vector<double>& evalDates,
					const std::vector<double>& numDates,
					bool doPDE=true);

	virtual string toString( const string& indent="",
							 const string& nextIndent="") const;

private:
	double		itsMeanReversion;
	ARM_Curve	itsVolCurve;
	ARM_Curve	itsVarCurve;
	double		itsHorizonDate;
	

protected:
	double		itsFwd;
	double		itsNumFwd;
	double		itsDelta;
};





////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderSLN
////////////////////////////////////////////////////

class ARM_ProcessBuilderSLN : public ARM_ProcessBuilder
{
public:
	ARM_ProcessBuilderSLN ( const ARM_ProcessBuilderSLN& rhs );
	ARM_ProcessBuilderSLN ( double mrs , const ARM_Curve& volCurve );
	virtual ~ARM_ProcessBuilderSLN(){};

	virtual double GetEqShift() const		{return itsShift;};
	virtual double GetEqVol() const			{return itsVol;};

	virtual	ARM_GP_VectorPtr Rate (	double t,
									const ARM_PricingStatesPtr& states,
									size_t modelNb,
									size_t k , bool extrapol=false) const;

	virtual	ARM_GP_VectorPtr DxRate (	double t,
										const ARM_PricingStatesPtr& states,
										size_t modelNb,
										size_t k , bool extrapol=false) const;

	virtual ARM_GP_VectorPtr DFRatio (	double t,
								const ARM_PricingStatesPtr& states,
								size_t modelNb,
								size_t k, bool extrapol=false) const;

	virtual ARM_GP_VectorPtr DxLogDFRatio(	double t,
									const ARM_PricingStatesPtr& states,
									size_t modelNb,
									size_t k , bool extrapol=false) const;

	virtual double	GetDiffusion() const {	return (itsShift + itsFwd) * itsScaling; };
	virtual void	ComputeProxyForCalib(size_t calibProxy,const ARM_VanillaSecurityDensity& density,ARM_Curve* volatilityCalib,double strike=0){};

	virtual void Calibrate(	double horizonDate,
					const ARM_VanillaSecurityDensity& density,
					const std::vector<double>& evalDates,
					const std::vector<double>& numDates,
					bool doPDE=true);

	//standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_ProcessBuilderSLN(*this);	};

private:
	double		itsVol;
	double		itsShift;
	double		itsScaling;

};



class SLN_Approx : public ARM_GP::UnaryFunc<double,double> 
{
public: 
		SLN_Approx(double fwd,double time,double strike1,double target1,double strike2,double target2,const ARM_VanillaSecurityDensity& density) :
		itsFwd(fwd),
		itsTime(time),
		itsTarget1(target1),
		itsStrike1(strike1),
		itsTarget2(target2),
		itsStrike2(strike2),
		itsDensity(density)
		{
		};

		virtual double operator() (double shift) const;
		double GetVol(double shift) const;
private:
	double						itsFwd;
	double						itsTime;
	double						itsStrike1;
	double						itsTarget1;
	double						itsStrike2;
	double						itsTarget2;
	const ARM_VanillaSecurityDensity&	itsDensity;
};

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
////////////////////////////////////////////////////


class ARM_ProcessBuilderPDE : public ARM_ProcessBuilder
{
private:

	inline double NbStdDev() const			{ return itsNbStdDev;}

	inline double GridSize() const			{ return itsGridSize;}

	inline double Theta() const				{ return itsTheta;}

	inline double GlobalVariance() const	{ return itsGlobalVariance; }

	void	BuildGrid();

	void	FreeMatrix();

	bool	IsDateInSchedule( double t) const {	return ExistsInVector( itsDates, t ); };
	
	size_t	NDaysIdxFromDate( const std::vector<double>& dateVector, double date, double daysNb ) const;
	
	size_t	IdxFromDate( const std::vector<double>& dateVector, double date ) const;

	void	StoreDFRatio(size_t i,const ARM_GP_VectorPtr& func);

	void	StoreDxLogDFRatio(size_t i,const ARM_GP_VectorPtr& func);

	void	StoreRate(size_t i,const ARM_GP_VectorPtr& func);

	void	StoreDxRate(size_t i,const ARM_GP_VectorPtr& func);

	void	RescaleRate(size_t index_storage, double scale);

	void	Backward(	double first,
						double last,
							ARM_GP_VectorPtr& func,
						const std::vector<double>& numDates);
	
	ARM_VectorPtr	GetTerminalFunc(const ARM_VanillaSecurityDensity& density);

public:

	ARM_ProcessBuilderPDE ( const ARM_ProcessBuilderPDE& rhs );
	ARM_ProcessBuilderPDE ( double mrs , const ARM_Curve& volCurve, bool withRescalling=false );
	virtual ~ARM_ProcessBuilderPDE(){};

	const std::vector<double>& GetDates() const	{ return itsDates;}

	inline double GetDate(size_t k) const	{ return itsDates[k];} 

	virtual double GetEqShift() const 		{ return itsEqShift; };

	virtual double GetEqVol() const			{return itsEqVol;};

	void	SetPDEParams(	size_t gridSize,
							double nbStdDev,
							double theta);

	virtual	ARM_GP_VectorPtr Rate (	double t,
									const ARM_PricingStatesPtr& states,
									size_t modelNb,
									size_t k , bool extrapol=false) const;

	virtual	ARM_GP_VectorPtr DxRate (	double t,
										const ARM_PricingStatesPtr& states,
										size_t modelNb,
										size_t k , bool extrapol=false) const;

	std::vector<double>& RateAtResetDate (	size_t t,
										const ARM_PricingStatesPtr& states,
										size_t modelNb,
										size_t k ) const;

	std::vector<double>& DxRateAtResetDate (	size_t t,
										const ARM_PricingStatesPtr& states,
										size_t modelNb,
										size_t k ) const;

	std::vector<double>& RateBetweenResetDates (	double ratio, size_t a, size_t b,
											const ARM_PricingStatesPtr& states,
											size_t modelNb,
											size_t k ) const;

	std::vector<double>& DxRateBetweenResetDates ( double ratio, size_t a, size_t b,
											 const ARM_PricingStatesPtr& states,
											 size_t modelNb,
											 size_t k ) const;

	virtual	ARM_GP_VectorPtr DFRatio (	double t,
										const ARM_PricingStatesPtr& states,
										size_t modelNb,
										size_t k , bool extrapol=false) const;

	virtual	ARM_GP_VectorPtr DFRatioAtResetDate (	size_t t,
													const std::vector<double>& states) const;

	virtual	ARM_GP_VectorPtr DFRatioAtResetDate (	size_t t,
													const ARM_PricingStatesPtr& states,
													size_t modelNb,
													size_t k ) const;

	virtual ARM_GP_VectorPtr DxLogDFRatio(	double t,
											const ARM_PricingStatesPtr& states,
											size_t modelNb,
											size_t k , bool extrapol=false) const;

	virtual ARM_GP_VectorPtr DxLogDFRatioAtResetDate(	size_t t,
														const ARM_PricingStatesPtr& states,
														size_t modelNb,
														size_t k ) const;

	virtual double GetDiffusionTerm( double t ) const;
	virtual double GetDLocalVolDRate( double t ) const;

	virtual void	Calibrate(	double horizonDate,
					const ARM_VanillaSecurityDensity& density,
					const std::vector<double>& evalDates,
					const std::vector<double>& numDates,
					bool doPDE=true);

	virtual double	GetDiffusion() const {return itsDiffusionTerm;};
	virtual void	ComputeProxyForCalib(size_t calibProxy,const ARM_VanillaSecurityDensity& density,ARM_Curve* volatilityCalib,double strike=0);

	void	ComputeProxyAtm(const ARM_VanillaSecurityDensity& density,bool black);
	void	ComputeProxyGauss(const ARM_VanillaSecurityDensity& density,bool atm,double strike);
	void	ComputeProxyMoment(const ARM_VanillaSecurityDensity& density);
	void	ComputeProxyLocal(const ARM_VanillaSecurityDensity& density,ARM_Curve* volatilityCalib,bool rescaling,bool black);
	void	ComputeProxyEffective(const ARM_VanillaSecurityDensity& density,ARM_Curve* volatilityCalib,bool rescaling);

	double	RescaleVolatilityCalib(ARM_Curve* volatilityCalib);

	virtual double GetProxyShift(double strike1, double strike2, const ARM_VanillaSecurityDensity& density);
	virtual double SetProxyVol(double strike,const ARM_VanillaSecurityDensity& density);
	
	double GetIthMoment(size_t i) const;
	size_t GetIndex(double x, const std::vector<double>& v) const;
	
	//standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_ProcessBuilderPDE(*this);	};
	virtual string toString( const string& indent="",
							 const string& nextIndent="") const;

private:

	std::vector<double>					itsDates;

	//for PDE
	double							itsGlobalVariance;
	size_t							itsGridSize;
	double							itsNbStdDev;
	double							itsTheta;
	std::vector<double>					itsGrid;

	//Left Tridiag
	std::vector<double>					itsL1;
	std::vector<double>					itsL2;
	std::vector<double>					itsL3;
	//Right Tridiag
	std::vector<double>					itsR1;
	std::vector<double>					itsR2;
	std::vector<double>					itsR3;
	//Right Vector
	std::vector<double>					itsR;
	//For inversion
	std::vector<double>					itsGam;
	std::vector<double>					itsXG;


	//storing
	ARM_VectorPtrVector				itsDFRatio;
	ARM_VectorPtrVector				itsDxLogDFRatio;
	ARM_VectorPtrVector				itsRate;
	ARM_VectorPtrVector				itsDxRate;
	
	//sln approx
	double							itsEqShift;
	double							itsEqVol;	
	double							itsEqScaling;	
	double							itsDiffusionTerm;

	bool							itsWithRescalling;
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

