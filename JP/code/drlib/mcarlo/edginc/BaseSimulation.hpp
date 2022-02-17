//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : BaseSimulation.hpp
//
//   Description : Declares the CopulaModelBase, data and functions common to
//                 all copulas in CCM
//
//   Date        : Feb 2006
//
//----------------------------------------------------------------------------

#ifndef BASESIMULATION_HPP
#define BASESIMULATION_HPP

#include "edginc/DoubleMatrix.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/MCRandom.hpp"
#include "edginc/CDSParSpreads.hpp"
#include "edginc/CleanSpreadCurve.hpp"
#include "edginc/CDSHelper.hpp"
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/IBaseSimulation.hpp"
#include "edginc/ICopulaModelBase.hpp"

DRLIB_BEGIN_NAMESPACE

#define TINY 1e-16
#define ZERO_EVENT 10

class ICopulaModelBase;


class BaseSimulation : public virtual VirtualDestructorBase,
					   public virtual IBaseSimulation
{
public:

	class SimulationState: public virtual IBaseSimulation::ISimulationState
	{

	protected:

		int numNames;
		IntArray defIdx;
		DoubleArray recoveries;

	public:

		SimulationState(int n) :
			numNames(n),
			defIdx(n),
			recoveries(n)
		{};

		~SimulationState();

		virtual int getNumAssets() const {return numNames;};

		virtual const DoubleArray& getRecoveries() const {return recoveries;};

		virtual const IntArray& getDefaultIndices() const {return defIdx;};

		virtual void setDefaultIndex(int i, int value);

		virtual void setRecovery(int i, double value);
	};

protected:

	const ICopulaModelBase *model;

	DateTime valueDate;              // taken from the model

	long numSamples;

	int  randsPerPath;

	// random numbers, directly from generator outside
	const double *randoms;           
	
	SimulationState simState;

	DateTimeArray   timeLine;

	// term structure of default probs per name
	CDoubleMatrix  probabilities;    

public:

	BaseSimulation(
		const DateTime& valueDate,
		const ICopulaModelBase &mod,
		DateTimeArray  dates,
		long  nSamples);

	BaseSimulation() : 
		simState(1),  model(NULL) {};

	virtual ~BaseSimulation() {};


	virtual void computeSample(long idx) 
	{
		throw ModelException("Method not defined");
	};

	virtual const ICopulaModelBase* Model() const {return model;};

	virtual const IBaseSimulation::ISimulationState &simulationState() const 
		{return simState;};

	virtual int getNumAssets() const { return model->getNumAssets();};

	virtual int getNumberRands() const { return randsPerPath;};

	virtual void setRandoms(IMCRandomSP rand)
	{
		randoms = &(rand->getRandomNumbers()[0][0]);
	};

	virtual const DateTimeArray& timeline() const;

	virtual void UpdateSim();

};

DECLARE(BaseSimulation);

DRLIB_END_NAMESPACE

#endif