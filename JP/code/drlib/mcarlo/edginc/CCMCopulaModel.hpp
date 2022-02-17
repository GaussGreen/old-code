//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : CCMCopulaModel.hpp
//
//   Description : Interface for class CCMCopulaModel, which implements the 
//				   CopulaModel interface as a mixture of an RFL copula, an
//                 independence copula and a dependece copula
//
//   Date        : Feb 2006
//
//----------------------------------------------------------------------------

#ifndef CCMCOPULAMODEL_HPP
#define CCMCOPULAMODEL_HPP

#include "edginc/RFLCopulaModel.hpp"
#include "edginc/DependentCopulaModel.hpp"
#include "edginc/IndependentCopulaModel.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE


/** implementation of Composite Copula Model

	CopulaModelData has the input curves per name (market par spreads)

	Another class should have the CCM specific parameters
*/


class CCMCopulaModel :	public CopulaModelBase
{			

public :

	CCMCopulaModel(
		DependentCopulaModelConstSP depenCop, 
		IndependentCopulaModelConstSP indepCop,
		RFLCopulaModelConstSP	rflCop
	) :
	CopulaModelBase(),
	depenCopula(*depenCop.get()), 
	indepCopula(*indepCop.get()), 
	rflCopula(*rflCop.get())	
	{}

	CCMCopulaModel() : CopulaModelBase(){};

	virtual ~CCMCopulaModel() {}

	virtual COPULA_MODEL type() const {return CCM;}

	virtual void addNames(const CreditAssetWrapperArray& xnames);

	virtual void updateBetas(const vector<double>& betas);

	class Simulation :   public BaseSimulation {

		DependentCopulaModel::SimulationSP dcSim;
		IndependentCopulaModel::SimulationSP indSim;
		RFLCopulaModel::SimulationSP rflSim;

	public:

		friend class CCMCopulaModel;

		// this functions access model and knows how to compute the probs of default for each of the names
		// on each of the copulae
		Simulation(
			const DateTime& valueDate,
			CopulaModelBase &model, 
			DateTimeArray  timeLine, 
			int nSamples);

		Simulation() 
			: BaseSimulation()
		{};

		void computeSample(long idx);

		void UpdateSim();

		void setRandoms(IMCRandomSP rand);

	};

	Simulation simulation(DateTimeArray timeLine, int nSamples)
	{
		return Simulation(valueDate, *this, timeLine, nSamples);
	}

	// this wont be used by the inner simulations, they will go to their own copulas
	virtual const CopulaModelBase& getCopulaData() const {return *this;};

	virtual void populateInnerCopulas();
	
	virtual void ComputeCurves(const DateTimeArray &dates);

	virtual void Update();

	// public?

	DoubleArray           ratios;
	DoubleArray           cataRecFactor;
	DoubleArray			  betaTweak;
	ICDSParSpreadsWrapperArray floorCurves;
	
protected:

	friend class Simulation;

	DependentCopulaModel depenCopula;
	IndependentCopulaModel indepCopula;
	RFLCopulaModel rflCopula;

	bool curvesComputed;

};

DECLARE(CCMCopulaModel)


DRLIB_END_NAMESPACE

#endif




