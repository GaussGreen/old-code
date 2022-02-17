//----------------------------------------------------------------------------
#ifndef QR_RegimeFactor_HPP
#define QR_RegimeFactor_HPP

#include "edginc/DateTime.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/GetMarket.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/AtomicArray.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL RegimeFactor: public MarketObject,
				  virtual public IGetMarket {
public:
    static CClassConstSP const TYPE;
    friend class RegimeFactorHelper;

	virtual void validate();

	virtual string getName() const;

    virtual void getMarket(const IModel* model, const MarketData* market, const string& name);

	virtual void calcTransitionProb(double tradYear);

	double expect(double tradYear);

    int getInitStateIdx() const;
    
    int getNbStates() const;

    DoubleArray getStates() const;

    void setInitState(double level);

    void setStates(DoubleArray levels);

    virtual ~RegimeFactor();

private:
	RegimeFactor& operator=(const RegimeFactor& rhs);

	/** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    RegimeFactor();

    static IObject* defaultCtor();

protected:
    RegimeFactor(const CClassConstSP& clazz);

    void buildCache();

	// fields ...
	string			name;
	DoubleArray		states;
	DoubleMatrix	generator;
	DoubleMatrix	jumpsOnLogSpot;
	int				initStateIdx;

	// transient fields ...
	int				nbStates;
	DoubleMatrix	transitionProb;
};

typedef smartConstPtr<RegimeFactor> RegimeFactorConstSP;
typedef smartPtr<RegimeFactor> RegimeFactorSP;
typedef array<RegimeFactor, RegimeFactor> RegimeFactorArray;

DRLIB_END_NAMESPACE

#endif