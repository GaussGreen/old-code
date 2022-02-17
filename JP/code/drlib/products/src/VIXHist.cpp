//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VIXHist.cpp
//
//   Description : Construct Historical VIX instrument
//
//   Author      : Xiaolan zhang
//
//   Date        : 6 February 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/UtilFuncs.hpp"
#include "edginc/VolRequestTime.hpp"
#include "edginc/Format.hpp"
#include "edginc/Results.hpp"
#include "edginc/Control.hpp"

#include "edginc/Instrument.hpp"
#include "edginc/ClosedForm.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/Theta.hpp"
#include "edginc/YieldCurve.hpp"

#include <float.h>

//#include "edginc/OutputRequestUtil.hpp"

DRLIB_BEGIN_NAMESPACE

/** VIXHist instrument - constructing VIX from Call/Put prices*/

class VIXHist: public CInstrument, 
				public virtual ClosedForm::IIntoProduct,
				public virtual LastSensDate,
				public virtual Theta::Shift {
public:
    static CClassConstSP const TYPE;
    friend class VIXHistHelper;
    friend class VIXHistClosedForm;
             
    /** copy market data relevant to the instrument */
    virtual void GetMarket(const IModel*, const CMarketDataSP);

    virtual void Validate();

    /** Implementation of ClosedForm::IntoProduct interface */
    virtual ClosedForm::IProduct* createProduct(ClosedForm* model) const;

    /** what's today ? */
    virtual DateTime getValueDate() const;

    /** when to stop tweaking */
	virtual DateTime endDate(const Sensitivity* sensControl) const;

    virtual bool sensShift(Theta* shift);

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;

private:
	VIXHist();
	VIXHist(CClassConstSP clazz);
	VIXHist(CClassConstSP clazz, const string&  discountName) ;

    VIXHist(const VIXHist& rhs);
    VIXHist& operator=(const VIXHist& rhs);

	void price(Control* control, CResults* results) const ;
    void requests(Control* control, CResults* results) const;

	/** calc fwd Index Level for a given maturity 
		return atmKIndex */
	double calcFwdIndexLevel(const DateTime maturity, const DoubleArray& Ks, const DoubleArray& calls, const DoubleArray& puts) const;
	int calcK0Index(double fwd, const DoubleArray& Ks) const;
	int selectCalls(const DateTime maturity, int KZeroIndex, const DoubleArray& calls)const;
	int selectPuts(const DateTime maturity, int KZeroIndex,  const DoubleArray& puts) const;

	/** using the midpoint of the bid-ask spread for each option with strike K_i */
	DoubleArraySP prepareOptions(int bot, int top, int KZeroIndex, 
						const DoubleArray& callsBid, const DoubleArray& callsAsk, 
						const DoubleArray& putsBid, const DoubleArray& putsAsk) const;

	/** deltaKs_j = 0.5 *(K_(j+1) - K_j)
		deltaKs_bot = (K_(bot+1) - K_bot)
		deltaKs_top = (K_(top) - K_(top-1)) */
	void prepareDeltaKs(int bot, int top, const DoubleArray& Ks, DoubleArray& deltaKs)const;

	/** calc vol^2 for a given maturity */
	double calcVar(DateTime maturity, double timeToMaturity, 
						const DoubleArray& Ks, 
						const DoubleArray& callsBid, const DoubleArray& callsAsk, 
						const DoubleArray& putsBid, const DoubleArray& putsAsk)const;


	/** calc vol^2 for a given maturity */
	double calcVIXHist() const;

	/** inputs */
	/** near term */
    DateTime		nearTerm;				//near term options' maturity
	DoubleArray		nearTermKs;				//near term strikes 
	DoubleArray		nearTermCallsBid;		//near term Call Price (Bid)
	DoubleArray		nearTermCallsAsk;		//near term Call Price (Ask)
	DoubleArray		nearTermPutsBid;		//near term Put Price (Bid)
	DoubleArray		nearTermPutsAsk;		//near term Put Price (Ask)

	/** next term */
    DateTime		nextTerm;				//next term options' maturity
	DoubleArray		nextTermKs;				//next term strikes 
	DoubleArray		nextTermCallsBid;		//next term Call Price (Bid)
	DoubleArray		nextTermCallsAsk;		//next term Call Price (Ask)
	DoubleArray		nextTermPutsBid;		//next term Put Price (Bid)
	DoubleArray		nextTermPutsAsk;		//next term Put Price (Ask)

	DateTime		VIXTerm;				//valueDate + 30	
    DateTime		valueDate;
	YieldCurveWrapper discount;
};

void VIXHist::Validate() {
    static const string method = "VIXHist::Validate";
    try {
		int sizeNearTerm = nearTermKs.size();
		/** all vectors should be the same size */
		if (sizeNearTerm != nearTermCallsBid.size() ||
			sizeNearTerm != nearTermCallsAsk.size() ||
			sizeNearTerm != nearTermPutsBid.size() ||
			sizeNearTerm != nearTermPutsAsk.size()
			 ){
            throw ModelException(method,
                                    "The number of near term strikes (size = " + Format::toString(nearTermKs.size()) +
                                    "), Bid Call prices (size =" + Format::toString(nearTermCallsBid.size()) +
                                    "), Ask Call prices (size =" + Format::toString(nearTermCallsAsk.size()) + 
                                    "), \n Bid Put prices (size =" + Format::toString(nearTermPutsBid.size()) + 
                                    "), Ask Put prices  (size =" + Format::toString(nearTermPutsAsk.size()) + 
                                    ") must be equal.  ");
		}		
        int sizeNextTerm = nextTermKs.size();
		if (sizeNextTerm != nextTermKs.size() ||			
			sizeNextTerm != nextTermCallsBid.size() ||
			sizeNextTerm != nextTermCallsAsk.size() ||
			sizeNextTerm != nextTermPutsBid.size() ||
			sizeNextTerm != nextTermPutsAsk.size() ){
            throw ModelException(method,
                                    "The number of next term strikes (size = " + Format::toString(nextTermKs.size()) +
                                    "), Bid Call prices (size =" + Format::toString(nextTermCallsBid.size()) +
                                    "), Ask Call prices (size =" + Format::toString(nextTermCallsAsk.size()) + 
                                    "), \n Bid Put prices (size =" + Format::toString(nextTermPutsBid.size()) + 
                                    "), Ask Put prices  (size =" + Format::toString(nextTermPutsAsk.size()) + 
                                    ") must be equal.  ");
		}

		/** nearTermKs should be increasing */
		int i = 0;
		while (i < sizeNearTerm-1 && nearTermKs[i] < nearTermKs[i+1]){
			i++;
		}
		if (i != sizeNearTerm -1){
			throw ModelException(method, "The near term strikes aren't increasing at strikes [" + 
				Format::toString(i) + "] = " +
				Format::toString(nearTermKs[i]) + " and "  + 
				"at strikes[" + Format::toString(i+1) + "] = " + 
				Format::toString(nearTermKs[i+1])+ ".");
		}

		/** nextTermKs should be increasing */
		i = 0;
		while (i < sizeNextTerm-1 && nextTermKs[i] < nextTermKs[i+1]){
			i++;
		}
		if (i != sizeNextTerm -1){
			throw ModelException(method, "The next term strikes aren't increasing at " + 
				Format::toString(i) + "with strikes[" + Format::toString(i) + "] = " +
				Format::toString(nextTermKs[i]) + " and "  + 
				Format::toString(i+1) + "with strikes[" + Format::toString(i+1) + "] = " + 
				Format::toString(nextTermKs[i+1])+ ".");
		}

		/** nearTerm < VIXTerm < nextTerm */
//		if (nearTerm.isGreaterOrEqual(VIXTerm) || VIXTerm.isGreaterOrEqual(nextTerm)){
		if (nearTerm.isGreaterOrEqual(nextTerm) ){
			throw ModelException(method, " Near term  should be smaller than VIXTerm, VIXTerm should be smaller than Next Term. But, near term = " +
											nearTerm.toString() + ", VIXTerm = " + VIXTerm.toString() +
											" and Next Term = " +
											nextTerm.toString() +".");
		}

		/** VIXTerm should be 30 days */
		int vixTermDays = VIXTerm.daysDiff(valueDate);
		if (vixTermDays != 30){
			throw ModelException(method, "VIX Tenor should be 30 days instead of " + Format::toString(vixTermDays));
		}

		/** if nearTerm < VIXTerm, then, VixTerm - nearTerm < 8 days, should switch options */
		int t = nearTerm.daysDiff(valueDate);
		if (t < 8 && t > 0 ){
			throw ModelException(method, "The near term's time to maturity " + Format::toString(t) + 
								" is smaller than 8 days. Should switch to next set of options!");
		}
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// copy market data relevant to the instrument
void VIXHist::GetMarket(const IModel* model, const CMarketDataSP market) {
    market->GetReferenceDate(valueDate);
    discount.getData(model, market);
}

void VIXHist::price(Control* control, CResults* results) const {
    static const string method = "VIXHist::price";
    try {        
		double value = calcVIXHist();

        results->storePrice(value, discount->getCcy());

        if (control && control->isPricing() ) {
            requests(control, results);
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

void VIXHist::requests(Control* control, CResults* results) const {
    static const string method = "VIXHist::requests";
    try {
/*
        OutputRequest* request =
                control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request) {
            DateTimeArraySP dates(new DateTimeArray(cfl->getLength()));
            for (int i = 0; i < cfl->getLength(); i++) {
                (*dates)[i] = (*cfl)[i].date;
            }
            OutputRequestUtil::recordPaymentDates(control,results,dates.get()); 
        }
        
        request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request) {
            OutputRequestUtil::recordKnownCashflows(control,
                                                    results,
                                                    discount->getCcy(),
                                                    cfl.get()); 
        }
*/
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** private class */
class VIXHistClosedForm: public ClosedForm::IProduct{
private:
    const VIXHist* vix; // a reference

public:
    VIXHistClosedForm(const VIXHist* vix): vix(vix){}

    void price(ClosedForm* model,
               Control*    control, 
               CResults*   results) const {
        vix->price(control, results);
    }
};
    
/** Implementation of ClosedForm::IntoProduct interface */
ClosedForm::IProduct* VIXHist::createProduct(
    ClosedForm* model) const{
    return new VIXHistClosedForm(this);
}

/** what's today ? */
DateTime VIXHist::getValueDate() const {
    return valueDate;
}

/** when to stop tweaking */
/** not used , to review */
DateTime VIXHist::endDate(const Sensitivity* sensControl) const {
    return (nextTerm);
}

bool VIXHist::sensShift(Theta* shift) {
    try {
        valueDate = shift->rollDate(valueDate);
    }
    catch (exception& e) {
        throw ModelException(e, "VIXHist::sensShift (theta)");
    }    
    return true; // our components have theta type sensitivity
}

/** Returns the name of the instrument's discount currency */
string VIXHist::discountYieldCurveName() const {
    return discount.getName();
}

VIXHist::VIXHist(CClassConstSP clazz,
				const string&  discountName) : CInstrument(clazz), discount(discountName)
{
}

// for reflection
VIXHist::VIXHist(): CInstrument(TYPE) {}
VIXHist::VIXHist(CClassConstSP clazz): CInstrument(clazz) {}

class VIXHistHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VIXHist, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(ClosedForm::IIntoProduct);
        IMPLEMENTS(LastSensDate);
        IMPLEMENTS(Theta::Shift);
        EMPTY_SHELL_METHOD(defaultVIXHist);
		
		FIELD(discount, "identifies discount curve");
        FIELD(valueDate, "valuation date");
        //FIELD_MAKE_OPTIONAL(valueDate);

        FIELD(VIXTerm, "VIX Term = value Date + 30 days");

		//near term
        FIELD(nearTerm,			"maturity of options for the near term");
        FIELD(nearTermKs,		"strikes of options for the near term");
		FIELD(nearTermCallsBid,	"near term Call Price (Bid)");
        FIELD(nearTermCallsAsk,	"near term Call Price (Ask)");
		FIELD(nearTermPutsBid,	"near term Put Price (Bid)");
		FIELD(nearTermPutsAsk,	"near term Put Price (Ask)");

		//next term
        FIELD(nextTerm,			"maturity of options for the next term");
        FIELD(nextTermKs,		"strikes of options for the next term");
		FIELD(nextTermCallsBid,	"next term Call Price (Bid)");
        FIELD(nextTermCallsAsk,	"next term Call Price (Ask)");
		FIELD(nextTermPutsBid,	"next term Put Price (Bid)");
		FIELD(nextTermPutsAsk,	"next term Put Price (Ask)");
    }

    static IObject* defaultVIXHist(){
        return new VIXHist();
    }
};

CClassConstSP const VIXHist::TYPE = CClass::registerClassLoadMethod(
    "VIXHist", typeid(VIXHist), VIXHistHelper::load);
   
/** calc fwd Index Level for a given maturity */
double VIXHist::calcFwdIndexLevel(const DateTime maturity, const DoubleArray& Ks, const DoubleArray& calls, const DoubleArray& puts) const{
    static const string method("VIXHist::calcFwdIndexLevel");
    try {
		double diffP = DBL_MAX;  // infinity : big positive number
		int atmKi  = -1;

		/** determine the ATM strike: 
			ATM strike is the strike price at which the difference between the call and put prices is smallest.*/
		for (int i = 0; i < Ks.size(); ++i){
			double tempDiff = fabs(calls[i] - puts[i]);
			if (tempDiff < diffP ){
				diffP = tempDiff;
				atmKi = i;
			}
		}
		if (atmKi == -1){
            throw ModelException("Unable to determinate the ATM strike. The difference between Call and Put prices are all infinity.");
		}

		/** determine the fwd Index level
			F = Ks[atmKi]+ exp(RT) * (diffP) */
		double capi = 1.0/ discount->pv(valueDate, maturity);  //to check if should use fwd rate and with basis = continuous, to review
                                                                //and also, there maybe 1 day off due to exchange rules
		return Ks[atmKi] + capi * diffP;
	}
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** calc the index corresponding to K_0 */
int VIXHist::calcK0Index(double fwd, const DoubleArray& Ks) const{
    static const string method("VIXHist::calcK0Index");
    try {
		/* mode = -1,  position of the largest item smaller than or equal to target*/
		return  Neighbour(fwd, &*Ks.begin(), 0, Ks.size()-1,  -1);
	}
    catch (exception& e) {
        throw ModelException(e, method);
    }	
}

/** select Option used to construct the vol^2:
	1). Select call options that have strike prices > K_0 and a non_zero bid price.
		After encountering two consecutive calls with a bid price of zero, don't select any other calls
	2). Select put options that have strike prices < K_0 and a non_zero bid price.
		After encountering two consecutive puts with a bid price of zero, don't select any other puts 
	3). Averaging the quoted bid-ask prices for each option
	4). Select both the put and call with strikes price K_0. and averaging them to get a single value */
	
/** return number of calls we're going to use */
int VIXHist::selectCalls(const DateTime maturity, int K0Index, const DoubleArray& calls) const{
	int i;
	/** select calls */
	i = K0Index;
	while ( i < calls.size() && calls[i] != 0.0){
		i++;
	}
	if (i > 0 )  i--;
	return i;
}
/** return number of puts we're going to use */
int VIXHist::selectPuts(const DateTime maturity, int K0Index,  const DoubleArray& puts)const{
	int i;

	/** select puts */
	i = K0Index;
	while ( i >= 0 && puts[i] != 0.0){
		i--;
	}
	if (i < puts.size()-1) i++;
	return i ;
}

/** using the midpoint of the bid-ask spread for each option with strike K_i */
DoubleArraySP VIXHist::prepareOptions(int bot, int top, int K0Index, 
							 const DoubleArray& callsBid, const DoubleArray& callsAsk, 
							 const DoubleArray& putsBid, const DoubleArray& putsAsk) const {
	int i ;	 
	/** Option Prices used to construct the vol^2 */
	DoubleArraySP optionPrices = DoubleArraySP(new DoubleArray(callsBid.size()));

	/** puts */
	for (i = bot; i < K0Index; ++i){
		(*optionPrices)[i] = 0.5 * ( putsBid[i] + putsAsk[i]);
	}
	/** call */
	for (i = K0Index + 1; i <= top; ++i){
		(*optionPrices)[i] = 0.5 * ( callsBid[i] + callsAsk[i]);
	}

	/** at K_0 */
	(*optionPrices)[K0Index] = 0.25*(putsBid[K0Index] + putsAsk[K0Index]
								+callsBid[K0Index] + callsAsk[K0Index]);
	return optionPrices;
}

/** deltaKs_j = 0.5 *(K_(j+1) - K_j)
	deltaKs_bot = (K_(bot+1) - K_bot)
	deltaKs_top = (K_(top) - K_(top-1)) */
void VIXHist::prepareDeltaKs(int bot, int top, const DoubleArray& Ks,  DoubleArray& deltaKs) const{
	for (int i = bot + 1; i < top; ++i){
		deltaKs[i] = 0.5 * (Ks[i+1] - Ks[i-1]);
	}
	deltaKs[bot] = Ks[bot+1] - Ks[bot];
	deltaKs[top] = Ks[top] - Ks[top-1];
}

/** calc vol^2 for a given maturity */
double VIXHist::calcVar(DateTime maturity, double timeToMaturity, 
						const DoubleArray& Ks, 
						const  DoubleArray& callsBid, const DoubleArray& callsAsk, 
						const DoubleArray& putsBid, const DoubleArray& putsAsk
					  )const{
						  						 
	/** calc fwd Index Level for a given maturity */
	double fwd = calcFwdIndexLevel(maturity, Ks, callsBid, putsBid);  //to review, should use bid or ask??

	/** calc the index corresponding to K_0 */
	int K0Index = calcK0Index(fwd, Ks);

	/** return number of calls we're going to use */
	int top = selectCalls(maturity, K0Index, callsBid);  // Ks[bot]: lowest strike to use
	int bot = selectPuts(maturity, K0Index, putsBid);	// Ks[top]: highest strike to use

	/** using the midpoint of the bid-ask spread for each option with strike K_i */
	/** calc optionPrices */
	/** Option Prices used to construct the vol^2 */
	DoubleArraySP optionPrices = prepareOptions(bot, top, K0Index, 
								callsBid, callsAsk, 
								putsBid, putsAsk);

	/** calc delta Ks*/
	DoubleArray		deltaKs;			// deltaKs_j = 0.5 *(K_(j+1) - K_j)
										// deltaKs_bot = (K_(bot+1) - K_bot)
										// deltaKs_top = (K_(top) - K_(top-1))

	deltaKs.resize(Ks.size());
	prepareDeltaKs(bot, top, Ks, deltaKs);

	/** calc vol^2 */
	double capi = 1.0/ discount->pv(valueDate, maturity);  //to check if should use fwd rate and with basis = continuous, to review
	double sum = 0;
	for (int i = bot; i <= top; ++i){
		sum += deltaKs[i] * (*optionPrices)[i] / (Ks[i]*Ks[i]);
	}
	sum *= 2.0 * capi;
	sum -= (fwd / Ks[K0Index] - 1.0)*(fwd / Ks[K0Index] - 1.0);
	sum /= timeToMaturity;
	return sum;
}

double VIXHist::calcVIXHist() const{
	
	/** calc time to maturity */
    double timeToMaturity1 = nearTerm.daysDiff(valueDate);
    double timeToMaturityVIX = VIXTerm.daysDiff(valueDate);
    double timeToMaturity2 = nextTerm.daysDiff(valueDate);
    
    timeToMaturity1 = (timeToMaturity1 - 1.0 )/365.0;
    timeToMaturity2 = (timeToMaturity2 - 1.0 )/365.0;
    timeToMaturityVIX = timeToMaturityVIX/365.0;

	//near term
	double varNearTerm = calcVar(nearTerm, timeToMaturity1, 
                            nearTermKs, 
                            nearTermCallsBid, nearTermCallsAsk, 
                            nearTermPutsBid, nearTermPutsAsk);
	//next term
	double varNextTerm = calcVar(nextTerm, timeToMaturity2, 
							nextTermKs, 
							nextTermCallsBid, nextTermCallsAsk, 
							nextTermPutsBid, nextTermPutsAsk);

	double dist = timeToMaturity2 - timeToMaturity1;
	double ratio1 = (timeToMaturityVIX- timeToMaturity1 )* timeToMaturity2/dist;
	double ratio2 = (timeToMaturity2 - timeToMaturityVIX)* timeToMaturity1/dist;

	double vol = sqrt((ratio2 * varNearTerm + ratio1 * varNextTerm)/timeToMaturityVIX)*100.0;
	return vol;
}

bool VIXHistLoad()
{
    return (VIXHist::TYPE != 0);
}
   
DRLIB_END_NAMESPACE
