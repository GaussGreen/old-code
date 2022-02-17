//----------------------------------------------------------------------------
#ifndef QR_MultiRegimeFactor_HPP
#define QR_MultiRegimeFactor_HPP

#include "edginc/DateTime.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/GetMarket.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/RateRegime.hpp"
#include "edginc/VolRegime.hpp"
#include "edginc/HazardRegime.hpp"
#include "edginc/CreditSpreadRhoParallel.hpp"
#include "edginc/VolParallel.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL MultiRegimeFactor: public MarketObject,
				                    virtual public IGetMarket,
                                    virtual public CreditSpreadRhoParallel::RestorableShift,
                                    virtual public IRestorableWithRespectTo<VolParallel>{
public:
    static CClassConstSP const TYPE;
    friend class MultiRegimeFactorHelper;
	friend class FD1DRegimes;

	virtual void validatePop2Object();

	virtual string getName() const;

    virtual void getMarket(const IModel* model, const MarketData* market, const string& name);

	virtual void calcTransitionProb(DateTime date);

	virtual void calcTransitionProb(double tradYear);

    HazardRegimeSP  getHazardRegime() const;

    void update();

    virtual ~MultiRegimeFactor();

    /** CreditSpreadRhoParallel sensitivity support */
    string sensName(CreditSpreadRhoParallel* shift) const;
    bool sensShift(CreditSpreadRhoParallel* shift);
    void sensRestore(CreditSpreadRhoParallel* shift);

    /** Returns name identifying vol for vega parallel */
    virtual string sensName(const VolParallel*) const;
    /** Shifts the object using given shift */
    virtual TweakOutcome sensShift(const PropertyTweak<VolParallel>& tweak);
    /** Restores the object to its original form */
    virtual void sensRestore(const PropertyTweak<VolParallel>& tweak);

private:
	MultiRegimeFactor& operator=(const MultiRegimeFactor& rhs);

	/** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    MultiRegimeFactor();

    static IObject* defaultCtor();

	int getRateIdx(int iState) const;	
	int getVolIdx(int iState)const;
	int getHazardIdx(int iState)const;

protected:
    MultiRegimeFactor(const CClassConstSP& clazz);

    void buildCache();

	// fields ...
	string				name;
	RateRegimeSP		rate;
	VolRegimeSP			vol;
	HazardRegimeSP		hazard;

	// transient fields ...
	DoubleArray			rateStates;
	DoubleArray			volStates;
	DoubleArray			hazardStates;

	DoubleMatrix		generator;
	DoubleMatrix		jumpsOnLogSpot;
	DoubleMatrix		jumpsAtDrift;
	DoubleMatrix		transitionProb;

    IntArray            rateIdx;
    IntArray            volIdx;
    IntArray            hazardIdx;

	int					initRateStateIdx;
	int					initVolStateIdx;
	int					initHazardStateIdx;
	int					initState;

	int					nbRateStates;
	int					nbVolStates;
	int					nbHazardStates;
	int					nbStates;
	int					nbFactors;
};

typedef smartConstPtr<MultiRegimeFactor> MultiRegimeFactorConstSP;
typedef smartPtr<MultiRegimeFactor> MultiRegimeFactorSP;

// support for wrapper class
typedef MarketWrapper<MultiRegimeFactor> MultiRegimeFactorWrapper;
#ifndef QLIB_MULTIREGIMEFACTOR_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<MultiRegimeFactor>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<MultiRegimeFactor>);
#endif

DRLIB_END_NAMESPACE

#endif
