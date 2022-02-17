//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : IBaseSimulation.hpp
//
//   Description : 
//
//   Date        : Oct 2006
//
//----------------------------------------------------------------------------

#ifndef IBASESIMULATION_HPP
#define IBASESIMULATION_HPP

#include "edginc/DoubleMatrix.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

class ICopulaModelBase;


class IBaseSimulation {

public:

	class ISimulationState
	{

	public:

		virtual int getNumAssets() const = 0;

		virtual const DoubleArray& getRecoveries() const = 0;

		virtual const IntArray& getDefaultIndices() const = 0;

		virtual void setDefaultIndex(int i, int value) = 0;

		virtual void setRecovery(int i, double value) = 0;
	};

	virtual void computeSample(long idx) = 0; 

	virtual const ICopulaModelBase* Model() const = 0;

	virtual const ISimulationState &simulationState() const = 0;

	int virtual getNumAssets() const = 0;

	virtual int getNumberRands() const = 0;

	virtual void setRandoms(IMCRandomSP rand) = 0;

	virtual const DateTimeArray& timeline() const = 0;

	// virtual const IntArray State() const {return simState.getDefaultIndices();};

	virtual void UpdateSim() = 0;

};

DECLARE(IBaseSimulation);

DRLIB_END_NAMESPACE

#endif