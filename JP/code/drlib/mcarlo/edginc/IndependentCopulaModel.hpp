//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : IndependentCopulaModel.hpp
//
//   Description : Declarations for class IndependentCopulaModel, a class 
//                 implementing the independence copula in the context of 
//				   MC CCM
//
//   Date        : Feb 2006
//
//----------------------------------------------------------------------------

#ifndef DR_INDEPENDENTCOPULAMODEL_HPP
#define DR_INDEPENDENTCOPULAMODEL_HPP

#include "edginc/BaseSimulation.hpp"
#include "edginc/CopulaModelBase.hpp"

DRLIB_BEGIN_NAMESPACE

/** Implements Independence Copula
*/
class IndependentCopulaModel : public CopulaModelBase
{
public:

	IndependentCopulaModel() 
		: CopulaModelBase()
	{};

	virtual ~IndependentCopulaModel() {}

	virtual COPULA_MODEL type() const;

	virtual const CopulaModelBase& getCopulaData() const {return *this;};

	class Simulation :   public BaseSimulation {

		CDoubleMatrix zScoresInv;

	public:

		Simulation(
			const DateTime& valueDate,
			CopulaModelBase &model, 
			DateTimeArray timeLine,
			long nSamples);

		friend class IndependentCopulaModel;

		void setRandoms(const double *rands) {randoms = rands;};

		void computeSample(long idx);

		virtual void UpdateSim();

		Simulation() 
			: BaseSimulation() 
		{};

	};

	DECLARE(Simulation);

	SimulationSP simulation(
		DateTimeArray dates, 
		long nSamples)
	{
		return SimulationSP(new Simulation(valueDate, *this, dates, nSamples));
	}
};

DECLARE(IndependentCopulaModel);

DRLIB_END_NAMESPACE

#endif




