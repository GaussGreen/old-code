//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : IndependentCopulaModel.hpp
//
//   Description : Declarations for class RFLCopulaModel, a class 
//                 implementing the rfl copula in the context of 
//				   MC CCM
//
//   Date        : Feb 2006
//
//----------------------------------------------------------------------------

#ifndef DR_RFLCOPULAMODEL_HPP
#define DR_RFLCOPULAMODEL_HPP

#include "edginc/CDSParSpreads.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/BaseSimulation.hpp"
#include "edginc/CopulaModelBase.hpp"


DRLIB_BEGIN_NAMESPACE

class RFLCopulaModel	:	public CopulaModelBase
{
public:
	
	RFLCopulaModel() : 
	  CopulaModelBase() 
	{};

	virtual ~RFLCopulaModel() {};

	virtual void getMarket(const CModel* model, const CMarketDataSP market);

	class Simulation :   public BaseSimulation {

		CDoubleMatrix zScoresInv;

	public:

		Simulation(
			const DateTime& valueDate,
			CopulaModelBase &mod, 
			DateTimeArray timeLine,
			int nSamples);

		Simulation() : BaseSimulation() {};

		void setRandoms(const double *rands) {randoms = rands;};

		friend class RFLCopulaModel;

		void computeSample(long idx);

		virtual void UpdateSim();

	};

	DECLARE(Simulation);

	SimulationSP simulation(
		DateTimeArray dates,
		long nSamples)
	{
		return SimulationSP(new Simulation(valueDate, *this, dates, nSamples));
	};

	virtual const CopulaModelBase& getCopulaData() const {return *this;};

	COPULA_MODEL type() const {return RFL;};

private:

	DoubleArray betas;

public:

	void updateBetas(const DoubleArray &inBetas) {betas = inBetas;};

};


DECLARE(RFLCopulaModel);

DRLIB_END_NAMESPACE

#endif




