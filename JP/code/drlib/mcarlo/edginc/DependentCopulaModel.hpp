//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : DependentCopulaModel.hpp
//
//   Description : Declarations for class DependentCopulaModel, a class 
//                 implementing the dependence copula in the context of 
//				   MC CCM
//
//   Date        : Feb 2006
//
//----------------------------------------------------------------------------

#ifndef DR_DEPENDENTCOPULAMODEL_HPP
#define DR_DEPENDENTCOPULAMODEL_HPP

#include "edginc/BaseSimulation.hpp"
#include "edginc/CopulaModelBase.hpp"

DRLIB_BEGIN_NAMESPACE

/** Implements Dependence Copula
*/
class DependentCopulaModel : public CopulaModelBase
{
public:

	DependentCopulaModel() : 
		CopulaModelBase()
	{};

	virtual ~DependentCopulaModel() {}

	virtual COPULA_MODEL type() const;

	virtual const CopulaModelBase& getCopulaData() const {return *this;};

	class Simulation :   public BaseSimulation {

		CDoubleMatrix zScoresInv;

	public:

		Simulation(
			const DateTime& valueDate,
			CopulaModelBase &model, 
			DateTimeArray timeLine,
			long nSamples
		);

		void setRandoms(const double *rands) {randoms = rands;};

		friend class DependentCopulaModel;

		void computeSample(long idx);

		virtual void UpdateSim();

	protected:

		static bool SimulationLoad();

	};

	DECLARE(Simulation);

	SimulationSP simulation(
		DateTimeArray timeLine, 
		long nSamples)
	{
		return SimulationSP(new Simulation(valueDate, *this, timeLine, nSamples));
	}

};

DECLARE(DependentCopulaModel);

DRLIB_END_NAMESPACE
#endif




