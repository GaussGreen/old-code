#ifndef EDR_HazardRegime_HPP
#define EDR_HazardRegime_HPP

#include "edginc/config.hpp"
#include "edginc/RegimeFactor.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/MarketWrapper.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL HazardRegime: public RegimeFactor,
							   virtual public Calibrator::IAdjustable {
public:
	static CClassConstSP const TYPE;
	friend class HazardRegimeAddin;
	friend class MultiRegimeFactor;
	friend class LeastSquareFit;

	void validatePop2Object(); 

	double calcRiskyFactor(double tradYear) const;

	double calcCreditSpread(double tradYear) const;

	DoubleArray calcCreditSpread(DoubleArray tradYears) const;

    double calcEquityCreditDependency(int stateIdx, int whichMethod) const;

	virtual string getName() const;

	DoubleArray	getGeneratorEntryForCalib() const;

    DoubleArray	getJumpEntryForCalib() const;

    double getRecoveryRate() const;

    double	getGenerator(int fromState, int toState) const;

    double	getJumpsOnLogSpot(int fromState, int toState) const;

    double	getState(int iState) const;

    void    setJumpAsymmetry(double value);

    /** Called after adjustments have been made to fields (eg calibrator) */
    virtual void fieldsUpdated(const CFieldArray& fields);

protected:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    HazardRegime();

    static IObject* defaultCtor();

    /** Called after (calibrator) adjustments have been made to fields */
    void update();

	double			recoveryRate;
    double          jumpAsymmetry;              // in terms of lnS
	DoubleArray		generatorEntryForCalib;		// independent entries of generator written in array form
	IntArray		rowIdx;						// row index of generator for each element of generatorEntryForCalib
	IntArray		colIdx;						// col index of generator for each element of generatorEntryForCalib

    DoubleArray		jumpEntryForCalib;		    // independent entries of jumpOnLogSpot written in array form
	IntArray		rowIdxJump;					// row index of jumpOnLogSpot for each element of jumpEntryForCalib
	IntArray		colIdxJump;					// col index of jumpOnLogSpot for each element of jumpEntryForCalib
};
typedef smartPtr<HazardRegime> HazardRegimeSP;
typedef smartConstPtr<HazardRegime> HazardRegimeConstSP;

// support for wrapper class
typedef MarketWrapper<HazardRegime> HazardRegimeWrapper;
#ifndef QLIB_HAZARDREGIME_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<HazardRegime>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<HazardRegime>);
#endif

DRLIB_END_NAMESPACE

#endif
