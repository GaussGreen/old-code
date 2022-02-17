//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : COPPER.cpp
//
//   Description : 
//
//   Date        : Dec 2003
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/Maths.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/CliquetVolRequest.hpp"
#include <algorithm>
#include "edginc/HandlePaymentEvents.hpp"

DRLIB_BEGIN_NAMESPACE

/*****************************************************************************/

class COPPER: public GenericNFBase, 
              virtual public IMCIntoProduct{
private: 
    /* fields */
    DoubleArray                  weights;            // weights for initial cliquet
    DoubleArray                  cliquetWeights;     // weights for futur cliquet
    DateTimeArray                averageOutDates;    // average out dates
    DateTimeArray                cliquetDates;       // maturity dates of cliquets
    string                       perfType;           // performance type for the payoff of each cliquet
    double                       overallStrike;      // overall strike
    double                       cliquetStrike;      // cliquet strike
    bool                         payAtCliquetMat;    // if there is a payment at each cliquet date or at the option maturity
    
public:
    static CClassConstSP const TYPE;
    friend class COPPERMC;
    
    // validation
    void validatePop2Object(){
	static const string routine = "COPPER::validatePop2Object";
	GenericNFBase::validatePop2Object();
	/* we need some average out dates */
	if (averageOutDates.empty()) {
	    throw ModelException(routine, "No average out date given!");
	}
	
	/* we need some cliquet dates */
	if (cliquetDates.empty()) {
	    throw ModelException(routine, "No cliquet date given!");
	}

	/* Require some averaging out dates for each cliquet period */
	int iAvgOut(0);
	int nbrAverageOut;
	
	for(int iCliquet = 0; iCliquet < cliquetDates.size(); iCliquet++){
	    nbrAverageOut = 0;
            bool flag = true;
	    while(flag) {
		nbrAverageOut ++;
		iAvgOut ++;
                if (iAvgOut >= averageOutDates.size()) {
                    flag = false;}
                else {
                    flag = cliquetDates[iCliquet].isGreaterOrEqual(averageOutDates[iAvgOut]);
                }
	    }
	    if (nbrAverageOut == 0) {
		throw ModelException(routine, "No average out date in Cliquet " + Format::toString(iCliquet));
	    }
	}
        
	/* Makes no sense to have averaging after final cliquet date */
	const DateTime& lastAvg = averageOutDates[averageOutDates.size()-1];
	const DateTime& lastCpn = cliquetDates[cliquetDates.size()-1];
	
	if (lastAvg > lastCpn) {
	    throw ModelException(routine, "Cannot average on " + lastAvg.toString() + 
				 " since after final cliquet date " + lastCpn.toString());
	}
    
	/* check that we've got as many weights as there are assets 
	   and if % weights that they sum to 100% */
	AssetUtil::checkWeights(weights, assets->NbAssets());
	AssetUtil::checkWeights(cliquetWeights, assets->NbAssets());
	
	/* check that if the payments are at each cliquet maturity, the overall strike is 0 */
	if (payAtCliquetMat && !(Maths::isZero(overallStrike))) {
	    throw ModelException(routine, "The payment is at each cliquet maturity but the overall strike is not 0");
	}
        
	/* check the performance type */
	if (!(perfType == "C" || perfType == "Call" || 
	      perfType == "F" || perfType == "Forward" || 
	      perfType == "P" || perfType == "Put")) {
	    throw ModelException(routine, "The performance type is "+perfType+" but it should be C/Call/F/Forward/P/Put");
	}
    }
    
    // returns the array of all sampling dates the instrument will ever need
    // excluding (possibly) the ref level dates (handled in GenericNFBase)
    const DateTimeArray samplingDates() const {
        return DateTime::merge(averageOutDates,cliquetDates);
    }

    /* Implementation of MonteCarlo::IntoProduct interface - the
       implementation of this is below */
    virtual IMCProduct* createProduct(
	const MonteCarlo* model) const; // see below
    
private:
    COPPER(): GenericNFBase(TYPE) {} // for reflection
    COPPER(const COPPER& rhs); // not implemented
    COPPER& operator=(const COPPER& rhs); // not implemented
    
    static IObject* defaultCOPPER(){
	return new COPPER();
    }
    
    /* Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
	REGISTER(COPPER, clazz);
	SUPERCLASS(GenericNFBase);
	IMPLEMENTS(IMCIntoProduct);
	EMPTY_SHELL_METHOD(defaultCOPPER);
	FIELD(weights,         "Weights defining the basket for the first cliquet");
	FIELD(cliquetWeights,  "Weights defining the basket for the latter cliquets");
	FIELD(averageOutDates, "Simulation dates");
	FIELD(cliquetDates,    "Cliquet maturity dates");
	FIELD(perfType,        "Performance type of cliquet options: use F/Forward/C/Call/P/Put");
	FIELD(overallStrike,   "Overall Strike in Pct., relates to Notional");
	FIELD(cliquetStrike,   "Strike for each cliquet in Pct., relates to Notional");
	FIELD(payAtCliquetMat, "Payment at each cliquet maturity if True, payment at option maturity if False");
	clazz->setPublic();           // make visible to EAS/spreadsheet
	clazz->setDescription("Call on cliquets of basket in which the weight associated to each component is based on the performances in the previous cliquet");
    }
};

/* MC product class for super COPPER */
class COPPERMC: public IMCProduct,
		virtual public IMCProductLN,
                virtual public IHandlePaymentEvents {
private:
    const COPPER*            inst;     // reference to original instrument
    int                      nbAssets; // nicer
    
    /* historical values */
    double                   sumCliquetPerfSoFar;   // Cliquet performance
    DoubleArray              sumAverageOutSoFar;   // [nbAssets]
    DoubleArray              sumAverageInSoFar;    // [nbAssets]
    DoubleArray              preStockPerfSoFar;    // [nbAssets] - stock performances in the previous cliquet
    int                      iCliquetSoFar;        // number of the actual cliquet (begining from 0)
    
    /* averageOut dates merged with cliquet dates */
    DateTimeArray            allDates;            // the array of all the useful dates
    
    /* identifies averageOut dates and cliquet dates */
    IntArray                 averageOutMap;          // [nbAllDates] - convenient way to track averageOut dates
    IntArray                 cliquetMap;             // [nbAllDates] - convenient way to track cliquet dates
    IntArray                 nbAvgOutPerCliquetDate; // [nbCliquets] 
    
    /* payment date for each cliquet, discount factor curve for cliquet payment date */
    DateTimeArray            cliquetPayDates;
    DoubleArray              discountFactor;
 
    /* already known coupons in the case payAtCliquetMat=True */
    CashFlowArray            knownCoupons;

    /* used to sort performance */
    IntArray                 permutation;           // [nbAssets]
    
    /* function used to sort the weights with respect to an array of performances
       we define a comparison function over indices:
       i1 is strictly greater than i2 if 
       the performance of indice i1 is greater than the performance of indiice i2 */
    class DoubleSort{
	DoubleArray dbles;
    public:
	DoubleSort(DoubleArray dbles): dbles(dbles){}
	bool operator()(int i1, int i2){
	    return (dbles[i1] > dbles[i2]);
	}
    };

public: 
    COPPERMC(const COPPER*       inst,
	     const SimSeriesSP&  simSeries):
	IMCProduct(inst->assets.get(),
                   inst->valueDate,
                   inst->discount.get(),
                   inst->refLevel,
                   simSeries,
                   inst->pastValues,
                   inst->instSettle.get(),
                   simSeries->getLastDate()),
	inst(inst),
	nbAssets(getNumAssets()),
	sumCliquetPerfSoFar(0.0),
	sumAverageOutSoFar(nbAssets, 0.0),
	sumAverageInSoFar(nbAssets, 0.0),
	preStockPerfSoFar(nbAssets,0.0),
	iCliquetSoFar(0),
	nbAvgOutPerCliquetDate(inst->cliquetDates.size(), 0),
        cliquetPayDates(inst->cliquetDates.size(), DateTime()),
	discountFactor(inst->cliquetDates.size(), 0.0),
	knownCoupons(0, CashFlow()),
	permutation(nbAssets,0) {      
        bool isTrivial;
        
        /* identifies the averageOut dates */
	
        allDates = DateTime::merge(inst->averageOutDates,inst->cliquetDates);
	
        averageOutMap = DateTime::createMapping(allDates,
                                                inst->averageOutDates,
                                                isTrivial);
        
        /* identifies the cliquet dates */
        cliquetMap    = DateTime::createMapping(allDates,
                                                inst->cliquetDates,
                                                isTrivial);
        
        /* number of AvgOut per period */
        int iCliquet = 0; 
	
        for(int iStep = 0; iStep < inst->averageOutDates.size(); iStep++) {    
            if( (inst->averageOutDates[iStep]).isGreater(inst->cliquetDates[iCliquet])) {
                iCliquet++;
            }
            nbAvgOutPerCliquetDate[iCliquet]++;
        }
	
        /* in case of payment at each cliquet maturity, we compute the discount factors */
        if (inst->payAtCliquetMat) {
            for(iCliquet = 0; iCliquet < inst->cliquetDates.size(); iCliquet++) {
                cliquetPayDates[iCliquet] = inst->instSettle->settles(inst->cliquetDates[iCliquet], NULL/*asset*/);
                discountFactor[iCliquet] = inst->discount->pv(Today, cliquetPayDates[iCliquet]);
            }
        }
	
        /* initialise the permutation array */
        for(int iAsset=0; iAsset<nbAssets; iAsset++) {
            permutation[iAsset] = iAsset;
        }
    }
  
    /* Use this opportunity to do any LogNormal driven initialisation
       of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
	// empty
    }
    
    /* override the initial function in case of payment at each cliquet date implying no final discount */
    double pvFromPaymentDate() const {
	if (inst->payAtCliquetMat) {
	    return(1.0);
	}
	else {
	    return(IMCProduct::pvFromPaymentDate());
	}
    }
    
    void payoff(const IPathGenerator*  pathGen,
		IMCPrices&          prices) {
	int    beginIdx = pathGen->begin(0); // 0 <- same for all assets
	int    endIdx   = pathGen->end(0);   // 0 <- same for all assets
	int    iAsset;
	int    iStep;
	double        sumCliquetPerf = sumCliquetPerfSoFar;        
	DoubleArray   sumAverageOut = sumAverageOutSoFar;     
	DoubleArray   sumAverageIn  = sumAverageInSoFar;  /* contains only one date but keep the symetry betwen 
							     averageOutDates and averageInDates for futur improvment */  
	DoubleArray   preStockPerf  = preStockPerfSoFar;
	int           iCliquet      = iCliquetSoFar;  
        
	/* Special treatment for average in for first cliquet */
	if (iCliquet==0) {
	    for (iAsset = 0; iAsset < nbAssets; iAsset ++) {
		sumAverageIn[iAsset] = pathGen->refLevel(iAsset,0);
	    }
	}
	else {
	    sumAverageIn = sumAverageInSoFar;         
	}

	bool isAverageOutDate;   // true iff an averageOut date
	bool isCliquetDate;      // true iff a cliquet date
	
	for (iStep=beginIdx; iStep<endIdx; iStep++) {
	    isAverageOutDate = (averageOutMap[iStep]==0);   // true iff an averageOut date
	    isCliquetDate    = (cliquetMap[iStep]   ==0);   // true iff a cliquet date
      
	    /* if we are at an averageOut date, we update the sum of the averageOut values */
	    if (isAverageOutDate) {
		for(iAsset=0; iAsset<nbAssets; iAsset++) {
		    sumAverageOut[iAsset] += pathGen->Path(iAsset, 0)[iStep];
		}
	    }
          
	    /* if we are at a cliquet maturity date */
	    if (isCliquetDate) {
		DoubleArray averageOut(nbAssets,0.0);
		DoubleArray stockPerf(nbAssets,0.0);
		double basketPerf(0.0);
		double cliquetPerf(0.0);
	
		/* we compute the average of the values at the averageOut dates for each asset */
		for(iAsset=0; iAsset<nbAssets; iAsset++) {
		    averageOut[iAsset] = sumAverageOut[iAsset] / nbAvgOutPerCliquetDate[iCliquet];
		}
	
		// we compute the stock performance for each asset */
		for(iAsset=0; iAsset<nbAssets; iAsset++) {
		    stockPerf[iAsset] = averageOut[iAsset] / sumAverageIn[iAsset];    
		}
	
	
		/* for the futur cliquet, we compute the permutation of stockPerf with respect to the previous stockPerf */
		if (iCliquet > 0) {
		    DoubleSort sortHelper(preStockPerf);
		    sort(permutation.begin(),permutation.end(), sortHelper);
		}
    
		/* we compute basketPerf for each asset : 
		   Special treatment for the first cliquet */
		basketPerf = 0.0;
		
		if (iCliquet==0) {
		    for (iAsset = 0; iAsset < nbAssets; iAsset ++) {
			basketPerf += inst->weights[iAsset] * stockPerf[iAsset];
		    }
		}
		else {
		    for (iAsset = 0; iAsset < nbAssets; iAsset ++) {
			basketPerf += inst->cliquetWeights[iAsset] * stockPerf[permutation[iAsset]];
		    }
		}
		
		/* we compute cliquetPerf */
		if ( inst->perfType == "C" || inst->perfType=="Call" ) {
		    cliquetPerf = max(0.0, basketPerf - inst->cliquetStrike);
		}
		else {
		    if ( inst->perfType == "P" || inst->perfType=="Put" ) {
			cliquetPerf = max(0.0, inst->cliquetStrike - basketPerf);
		    }
		    else {
			cliquetPerf = basketPerf - inst->cliquetStrike;
		    }
		}
		
		// in case of payment at each cliquet maturity, we discount the payment
		if (inst->payAtCliquetMat) {
		    sumCliquetPerf += cliquetPerf * discountFactor[iCliquet]; 
                    if (pathGen->doingPast()) { 
                        knownCoupons.resize(iCliquet+1);
                        knownCoupons[iCliquet] = CashFlow(cliquetPayDates[iCliquet],inst->notional * cliquetPerf);
                    }
		}
		else {
		    sumCliquetPerf += cliquetPerf; 
		}
	
		/* Reset average out and in level for next cliquet*/
		/* remember the stock performance for the next cliquet */
		preStockPerf = stockPerf;       
	
		/* the sum for the average out dates is 0 */
		for(iAsset=0; iAsset<nbAssets; iAsset++) {
		    sumAverageOut[iAsset] = 0.0;
		}
	
		/* the end of this cliquet is the only average in date for the next cliquet */
		for(iAsset=0; iAsset<nbAssets; iAsset++) {
		    sumAverageIn[iAsset] = pathGen->Path(iAsset, 0)[iStep];
		}
		
		/* we move to the next cliquet */ 
		iCliquet ++;
	    }
	}
        
	// preserve values for past - before we modify/sort the cliquets.
	if (pathGen->doingPast()) { 
	    sumAverageOutSoFar = sumAverageOut; 
	    sumAverageInSoFar  = sumAverageIn; 
	    preStockPerfSoFar  = preStockPerf; 
      	    iCliquetSoFar      = iCliquet;
      
	    // in case of payment at each cliquet maturity, the past value is 0
	    if (inst->payAtCliquetMat) 
		sumCliquetPerfSoFar = 0.; 
	    else  
		sumCliquetPerfSoFar = sumCliquetPerf;
	}
	
	// Compute a payoff, but only when we have a "complete" situation : either 
	// doingPast() and all is past, or !doingPast().
	if (!pathGen->doingPast() || !hasFuture())
	    prices.add(inst->notional * max(0.0,sumCliquetPerf-inst->overallStrike)); 
    }
	
    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
				    int                     iAsset) const {
	static const string routine = "COPPERMC::getVolInterp";
	CVolRequestLNArray reqarr(1);
	const DateTime&    startDate = getRefLevel()->getAllDates().front();
	const DateTime&    today = getToday();
	const DateTime&    lastSimDate = getSimSeries()->getLastDate();
	bool               fwdStarting = startDate.isGreater(today);
	double             interpLevel = inst->cliquetStrike;
	
	// get hold of the future strike dates
	int numLiveCliqs = inst->cliquetDates.size() - iCliquetSoFar;
	if (numLiveCliqs<=0) {
	    throw ModelException(routine, "No future cliquets!?");
	}
	
	DateTimeArray liveCliqStartDates(numLiveCliqs);
	for (int iLiveCliquet = 0; iLiveCliquet < numLiveCliqs; iLiveCliquet++){
	    int iCliquet = iCliquetSoFar+iLiveCliquet-1;
	    liveCliqStartDates[iLiveCliquet] = iCliquet<0?startDate:inst->cliquetDates[iCliquet];
	}
    
	// same strike levels per cliquet (but may need to adjust first one)
	DoubleArray  strikes(numLiveCliqs, interpLevel);
	if (!fwdStarting){
	    // need to set first level to absolute strike - adjusted
	    // additionally for any average out samples for this cliquet
	    int iStep;
	    const DateTime& thisCliquetStart = iCliquetSoFar==0?
		startDate:inst->cliquetDates[iCliquetSoFar-1];
	    
	    // find first avg date of this cliquet
	    for(iStep = 0; iStep < allDates.size() && 
		    allDates[iStep] <= thisCliquetStart; iStep++)
		; // empty
	    
	    // then walk through counting past avg dates in this cliq
	    int numRemaining = nbAvgOutPerCliquetDate[iCliquetSoFar];
	    for(; iStep < allDates.size() && 
		    allDates[iStep] <= today; iStep++)
		numRemaining--;
	    
	    // something wrong!
	    if (numRemaining<0)
		throw ModelException(routine, "INTERNAL ERROR : numRemaining is " + Format::toString(numRemaining));
	    
	    // Can't set up refLevel earlier, 'cos need PathGen. First cliq has standard ref level
	    double refLevel =  iCliquetSoFar==0 ? pathGen->refLevel(iAsset, 0) : sumAverageInSoFar[iAsset];
	    // if numRemaining is 0 then the next date is a cliquet Dates, the value does not matter
	    if (numRemaining == 0) {
		strikes[0]=0.;
	    }
	    else {
		strikes[0] = ((nbAvgOutPerCliquetDate[iCliquetSoFar]) * refLevel * interpLevel
			      - sumAverageOutSoFar[iAsset])/ numRemaining;
	    }
	}
	reqarr[0] =  CVolRequestLNSP(new CliquetVolRequest(fwdStarting, 
							   liveCliqStartDates, 
							   lastSimDate,
							   strikes));
	return reqarr;
    }

    // Satisfy IHandlePaymentEvents interface
    void recordEvents(Control* control,
                      Results* results) {
        static const string method("BoostedNFBMC::recordEvents");
        try {
            // PAYMENT_DATES is a list of all dates on which payments may occur
            // including past and potential future dates.
           
            OutputRequest* request =
                control->requestsOutput(OutputRequest::PAYMENT_DATES);
            if (request && !request->getHasFinished()) {
                if (inst->payAtCliquetMat) {
                    OutputRequestUtil::recordPaymentDates(control,results,&cliquetPayDates);
                }
                else {
                    DateTime matDate = inst->instSettle->settles(inst->cliquetDates[inst->cliquetDates.size()-1], NULL /* asset */);
                    OutputRequestUtil::recordPaymentDates(control,results,&DateTimeArray(1,matDate));
                }
            }
        
            // KNOWN_CASHFLOWS should have dates a subset of PAYMENT_DATES
            // and be supplied for all past cash flows, and any future ones
            // that are determined.
            
            request =control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
            if (request && !request->getHasFinished()) {
                if (inst->payAtCliquetMat) {
                    if (knownCoupons.size()>0) {
                        OutputRequestUtil::recordKnownCashflows(control,
                                                                results,
                                                                discount->getCcy(),
                                                                &knownCoupons); 
                    }
                }
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

};
    

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* COPPER::createProduct(const MonteCarlo* model) const {
    
    // we create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); 
    /* create empty one */
    simSeries->addDates(samplingDates());
    return new COPPERMC(this, simSeries);
}

CClassConstSP const COPPER::TYPE = 
    CClass::registerClassLoadMethod("COPPER", typeid(COPPER), COPPER::load);

// * for class loading (avoid having header file) */
bool COPPERLoad() {
    return (COPPER::TYPE != 0);
}
    
DRLIB_END_NAMESPACE





