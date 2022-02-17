//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RollingAveAmer.cpp
//
//   Description : MC model for rolling average american exercise
//                 There is no optimisation performed by this model !
//                 The exercise strategy is taken as input !
//                 The two exercise criteria are:
//                   alpha = Spot/S(rolling ave)
//                   beta = S(rolling ave)/Strike
//
//                 Note that current implementation only supports exercise on rolling average.
//
//   Date        : 14 March 2003
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/Generic1Factor.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/Vanilla.hpp"

DRLIB_BEGIN_NAMESPACE

/** Rolling average american */
//////// instrument class //////////
class RollingAveAmer: public Generic1Factor, 
                      virtual public LastSensDate,
                      virtual public IMCIntoProduct{
protected:
    /// fields ////////
    SampleListSP            spotSamples;    // all historic spots
    bool                    isCall;
    DateTime                matDate;
    double                  strike;
    int                     rollNum;
    DateTimeArray           alphaBetaDates;
    DoubleArray             alpha;
    DoubleArray             beta;
    double                  shortStockCharge;

    // not to appear in IMS
    bool                    DEBUG_useCtrlVariate;
    string                  DEBUG_interp;

public:
    static CClassConstSP const TYPE;
    friend class RollingAveAmerProd;

    virtual void Validate(){
        static const string routine("RollingAveAmer::Validate");
        // NB don't call parent's validate - issues with fwdStart
        if (fwdStarting){
            throw ModelException(routine, "Fwd starting flag on "
                                 "Generic1Factor is not used");
        }

        

        if (oneContract){
            throw ModelException(routine, "oneContract flag on "
                                 "Generic1Factor is not used");
        }

        if (Maths::isNegative(shortStockCharge) || shortStockCharge >= 1.0) {
            throw ModelException(routine, "shortStockCharge must be non negative, and less than 1.");
        }

        int i;
        for (i = 0; i < alpha.size(); i++)
        {
            if (Maths::isNegative(alpha[i]) || alpha[i] > 1.0)
                throw ModelException(routine, "alpha must be non negative, and less than 1.");
            if (Maths::isNegative(beta[i]) || beta[i] < 1.0)
                throw ModelException(routine, "beta must be non negative, and greater than 1.");
        }

        AssetUtil::assetCrossValidate(asset.get(),
                                      false, //fwdStarting,
                                      startDate,
                                      valueDate,
                                      discount,
                                      this);

    }

    void validatePop2Object(){
        static const string routine("RollingAveAmer::validatePop2Object");

        int num = alphaBetaDates.size();
        if (num == 0){
            throw ModelException(routine, "no exercise strategy dates supplied");
        }
        if (alpha.size() != num){
            throw ModelException(routine, "alpha array must be the same size as alphaBetaDates");
        }
        if (beta.size() != num){
            throw ModelException(routine, "beta array must be the same size as alphaBetaDates");
        }
        if (!isCall){
            throw ModelException("only call is implemented now");
        }
    }
   
    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const;

    /** when to stop tweaking - default implementation assumes product can
        be priced with a MC */
    DateTime endDate(const Sensitivity* sensControl) const{
        DateTime instEnd  = instSettle->settles(matDate, asset.get());
        DateTime assetEnd = asset->settleDate(matDate);
        DateTime end = instEnd.isGreater(assetEnd) ? instEnd : assetEnd;
        return end;
    }

private:
    RollingAveAmer(): Generic1Factor(TYPE) {
        DEBUG_useCtrlVariate = true;
        DEBUG_interp = "L";
    }

    // for reflection
    RollingAveAmer(const RollingAveAmer& rhs); // not implemented
    RollingAveAmer& operator=(const RollingAveAmer& rhs); // not implemented

    static IObject* defaultRollingAveAmer(){
        return new RollingAveAmer();
    }

    //// roll through time (setting historic values)
    bool sensShift(Theta* theta){
        // use valueDate before it changes
        spotSamples->roll(theta->getUtil(valueDate), 0 /* iAsset */,
                         asset.get()); // then roll our past values
        Generic1Factor::sensShift(theta); // and then call parent's method
        return true; // continue to tweak components which implement Theta
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(RollingAveAmer, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(IMCIntoProduct);
        IMPLEMENTS(LastSensDate);   
        IMPLEMENTS(Theta::Shift);
        EMPTY_SHELL_METHOD(defaultRollingAveAmer);
        FIELD(spotSamples, "historical values");
        FIELD(isCall,    "is it a call option");
        FIELD(matDate,    "maturity date");
        FIELD(strike,    "Strike");
        FIELD(rollNum, "num of of points for rolling average");
        FIELD(alphaBetaDates, "date list corresponding to alpha and beta params");
        FIELD(alpha, "Spot over rolling average");
        FIELD(beta, "average over strike");
        FIELD(shortStockCharge, "a percentage charge on exercise");
        FIELD(DEBUG_useCtrlVariate, "true(default) = plain vanilla closed used, false=straight MC pricing only");
        FIELD_MAKE_OPTIONAL(DEBUG_useCtrlVariate);
        FIELD(DEBUG_interp, "L(default)=linear, S=stairs, N=none");
        FIELD_MAKE_OPTIONAL(DEBUG_interp);

        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

/* MC product class for super rainbow */
//////// product class //////////
class RollingAveAmerProd : public IMCProduct, virtual public IMCProductLN{
private:
    const RollingAveAmer*   inst; // reference to original instrument
    double                  callPut;
    DateTimeArray           sampleDates;
    DoubleArray             alphaArr;
    DoubleArray             betaArr;

public:
    
    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        // empty
    }

    /** equivalent to InstIntoMCProduct. */
    RollingAveAmerProd(const RollingAveAmer*         inst,
                IRefLevelSP refLevel,
                const SimSeriesSP&             simSeries) :
                IMCProduct(inst->asset.get(),
                          inst->valueDate,
                          inst->discount.get(),
                          refLevel,
                          simSeries,
                          inst->spotSamples,
                          inst->instSettle.get(),
                          simSeries->getLastDate()),
                          inst(inst){
                    
                    sampleDates = simSeries->getAllDates();
                    alphaArr.resize(sampleDates.size());
                    betaArr.resize(sampleDates.size());
                    // use schedule for interpolation
                    Schedule alpha(inst->alphaBetaDates, inst->alpha, inst->DEBUG_interp);
                    Schedule beta(inst->alphaBetaDates, inst->beta, inst->DEBUG_interp);
                    for (int i=0; i<sampleDates.size(); i++)
                    {
                        alphaArr[i] = alpha.interpolate(sampleDates[i]);
                        betaArr[i] = beta.interpolate(sampleDates[i]);
                    }
                }
        
    /** Called within the simulation loop */
    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) {

        const double* path = pathGen->Path(0,0); // access path
        
        // compute rolling average, size of rollAvg is pathGen->end() - pathGen->begin()
        int rollNum = inst->rollNum;
        int i_start = pathGen->begin(0);        
        int i_end = pathGen->end(0);
        i_start = Maths::max(i_start, rollNum-1);

        int idx = i_end - i_start-1;
        idx = Maths::max(idx, 0);

//        if (i_end < rollNum)
//            throw ModelException("RollingAveAmer::payoff", "not enough points for rolling average");

        //nothing happens before rollNum -1, pure init average period.
        //start to check the Am exericise rule at rollNum-1
        if (i_end >= rollNum ){

            DoubleArraySP rollAvg(new DoubleArray(i_end - i_start, 0.0));
            DoubleArray& rollingAvg= *rollAvg; // for convenience

            // do last point first
            for (int j=1; j<= rollNum; j++)
                rollingAvg[idx] += path[i_end - j];

            rollingAvg[idx] /= rollNum;
            // progessively go back
            idx --;
            while (idx >= rollNum-i_start-1 && idx >=0)
            {
                rollingAvg[idx] = rollingAvg[idx+1] + (path[idx+i_start-rollNum+1] - path[idx+i_start+1])/rollNum;
                idx --;
            }
        
            // get exercise value
            double pathValue;
            double strike = inst->strike;
            bool isExercised = false;
            for (int i = i_start; i<i_end; i++)
            {
                if (rollingAvg[i-i_start]>0.0
                    && path[i]/rollingAvg[i-i_start] < alphaArr[i]
                    && rollingAvg[i-i_start]/strike > betaArr[i])
                {
                    isExercised = true;
                    pathValue = path[i]*(1.0-inst->shortStockCharge)*(1.0 - strike/rollingAvg[i-i_start]);

                    // need to remove pv factor to mat + settle.
                    pathValue /= inst->instSettle->pv(sampleDates[i],   
                        inst->matDate,
                        inst->discount.get(),
                        inst->asset.get());
                    // and add pv factor from exerDate to exerDate + settle.
                    pathValue *= discount->pv(sampleDates[i],
                        inst->instSettle->settles(sampleDates[i], inst->asset.get()));

                    break;
                }
            }

            // maturity value if not exercised
            if (!isExercised)
            {
                pathValue = Maths::max(path[i_end-1]*(1.0-inst->shortStockCharge)-strike, 
                                path[i_end-1]*(1.0-inst->shortStockCharge)*(1.0 - strike/rollingAvg[i_end-1-i_start]));
                if (pathValue < 0.0)
                    pathValue = 0.0;
            }

            // control variate is done as a spread, can explicitly compute both amer and euro prices if needed
            if (inst->DEBUG_useCtrlVariate)
                pathValue -= Maths::max(0.0, path[i_end-1]*(1.0-inst->shortStockCharge)-strike);

            if (!pathGen->doingPast())
                prices.add(pathValue);
        }
    }

    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        DateTime imntStartDate = inst->fwdStarting? 
            inst->startDate: inst->valueDate;
        CVolRequestLNArray   reqarr(1); // one interp level/path per asset here

        reqarr[0] = CVolRequestLNSP(new LinearStrikeVolRequest(
            inst->strike,               // should really be strike/(1-shortStockCharge)
            imntStartDate, 
            inst->matDate,
            inst->fwdStarting));
        
        return reqarr;
    }

    // control variate is done here
    virtual void recordExtraOutput(Control* control, Results* results, const IMCPrices& ) const
    {
        if (inst->DEBUG_useCtrlVariate)
        {
            double euroPrice = CVanilla::priceBS(inst->getValueDate(),
                                                inst->startDate,
                                                inst->matDate,
                                                inst->isCall,
                                                false, // inst->fwdStarting, ??
                                                true, // oneContract=true ???
                                                inst->notional,
                                                inst->initialSpot,
                                                inst->strike/(1.0-inst->shortStockCharge),
                                                inst->instSettle.get(),
                                                inst->asset.get(),
                                                inst->discount.get());
            
            euroPrice *= (1.0-inst->shortStockCharge);

            results->storePrice(euroPrice + results->retrievePrice(),
                                results->getCcyName());

            if (control && control->isPricing())
            {
                OutputNameConstSP van1(new OutputName("vanilla-brrw-adj"));
                OutputNameConstSP van2(new OutputName("vanilla-unadj"));
                results->storeGreek(CDoubleSP(CDouble::create(euroPrice)), Results::DEBUG_PACKET, van1);
                
                // compute unadjusted vanilla price
                double euroPriceUnAdj = CVanilla::priceBS(inst->getValueDate(),
                    inst->startDate,
                    inst->matDate,
                    inst->isCall,
                    false, // inst->fwdStarting, ??
                    true, // oneContract=true ???
                    inst->notional,
                    inst->initialSpot,
                    inst->strike,
                    inst->instSettle.get(),
                    inst->asset.get(),
                    inst->discount.get());

                results->storeGreek(CDoubleSP(CDouble::create(euroPriceUnAdj)), Results::DEBUG_PACKET, van2);
                
            }
        }
    }
}; // end of class RollingAveAmerProd


/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* RollingAveAmer::createProduct(const MonteCarlo* model) const {
    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(1)); /* create empty one */
    simSeries->addDates(spotSamples->getAllDates());

    IRefLevelSP refLevel(IRefLevel::Util::makeTrivialAverage(spotSamples->getDates()[0]));
    return new RollingAveAmerProd(this, refLevel, simSeries);
}

CClassConstSP const RollingAveAmer::TYPE = CClass::registerClassLoadMethod(
    "RollingAveAmer", typeid(RollingAveAmer), RollingAveAmer::load);

// force linker to include this file (avoid having header file) */
bool RollingAveAmerLoad() {
    return (RollingAveAmer::TYPE != 0);
}

DRLIB_END_NAMESPACE
