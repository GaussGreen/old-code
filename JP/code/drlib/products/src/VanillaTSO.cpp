//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VanillaTSO.cpp
//
//   Description : VanillaTSO option with flat value div adjustment.
//                 Flat value adjustment is achieved using array of prices in strike dimension. Interpolation
//                 is used at the ex-div date to determine flat value strike adjustment.
//                 Allows gamma scaling of local vol to reflect vol hedging cost(vol reduced or increased according to 
//                 local gamma).
//
//                 Payoff takes in array of options (strikes, maturities, size(+/-), call/put, and div-adj flag. isAmer)
//
//                 price[] array usage: 
//                          [0] is the aggregate price
//                          [1] to [strikeGrid] for 1st tso option
//                          [strikeGrid+1] to [2*strikeGrid] for 2nd tso option
//                          [2*strikeGrid+1] to [3*strikeGrid] for 3nd tso option
//                          ... to NumTSO...
//                          [strikeGrid*NumTSO+1] to [NumOfPrice-1] are hedge options. no div adj.
//
//                 *** Not meant for LN use but if requested, vol is interpolated for the longest maturity option ***
//                 *** and if there are more than one, the last found will be used. So use flat vol in LN for testing. ***
//
//                 *** if dollar div flag is on it is assumed to be call type treatment !!!***
//
//                 *** tree exercise step set up also assumes call option !!! ***
//
//                 *** fwd start and averaging not supported  ***
//
//   Author      : Ning Shen
//
//   Date        : 21 October 2003
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Nrfns.hpp"
#include "edginc/UtilFuncs.hpp"
#include "edginc/Black.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/Tree1fLVGD.hpp"
#include "edginc/AggregateModel.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/Equity.hpp" // for dollar div

#include "edginc/Vanilla.hpp"
#include "edginc/DividendAdjusted.hpp"

#include "edginc/IndexSpecEQ.hpp"
#include "edginc/LatticeProdEDR.hpp"

DRLIB_BEGIN_NAMESPACE

/**********  local helpers **********/
// interp function for strike adjustment
inline void interpStrikeAdj(double* strikes, vector<double>& adjP, vector<double>& price)
{
// interp/extrap

    int i;
    int num = adjP.size();
    vector<double> adjK(num);
    vector<double> v = price;

    unsigned long j;
    // linear interp
    // first interp for adjusted strikes given adjusted price
    for (i=0; i< num; i++)
    {
        locate(&*adjP.begin()-1, num, price[i], &j);
        if ((int)j == num) j--; // outside upper, shift down
        if (j > 0) j--; // shift for 0 offset

        // interp strike adjusted
        adjK[i] = strikes[j]+(price[i]-adjP[j])*(strikes[j+1]-strikes[j])/(adjP[j+1]-adjP[j]);
    }
    // interp prices for original strikes
    double* k = &*adjK.begin()-1;
    for (i=0; i< num; i++)
    {
        locate(k, num, strikes[i], &j);
        if ((int)j == num) j--; // outside upper, shift down
        if (j > 0) j--; // shift for 0 offset

        price[i] = v[j]+(strikes[i]-adjK[j])*(v[j+1]-v[j])/(adjK[j+1]-adjK[j]);
    }

/*
//  spline interp, does not work as too many extrapolations?
    const double y1 = 2e30; // default to natural spline
    const double yn = 2e30;
    double result;

    vector<double> y2(num);

    // first interp for adjusted strikes given adjusted price
    double* p = &*adjP.begin()-1;
    spline(p, strikes-1, num, y1, yn, &*y2.begin()-1);
    for (i=0; i< num; i++)
    {
        splint(p, strikes-1, &*y2.begin()-1, num, (p[i]), &result);


if (fabs(result/strikes[i]-1.0)>1.0 || result<strikes[i])
{
    double what = strikes[i];
}

        adjK[i] = result;
    }
    // interp prices for original strikes
    double* k = &*adjK.begin()-1;
    p = &*v.begin()-1;
    spline(k, p, num, y1, yn, &*y2.begin()-1);
    for (i=0; i< num; i++)
    {
        splint(k, p, &*y2.begin()-1, num, strikes[i], &result);

if (fabs(result/(p[i])-1.0)>1.0 || result<(p[i]))
{
    double what = (p[i]);
}

        (p[i]) = result;
    }
*/

}
// ********** end of local helper ********

/** input for options. This is private data populated by createAggregateInstrument */
class TSOInput : public CObject
{
public:

    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);

    static IObject* defaultTSOInput(){
        return new TSOInput();
    }

    /** object validation */
    void validatePop2Object() {}

    TSOInput():  CObject(TYPE){}
    friend class VanillaTSO;

    bool        isCall;
    double      strike;
    DateTime    maturity;
    bool        isAmerican;
    double      weight;
};

typedef array<TSOInput, TSOInput> TSOInputArray;

///////////////////////////////////////////////////////
/** input for options */
CClassConstSP const TSOInput::TYPE = CClass::registerClassLoadMethod(
    "TSOInput", typeid(TSOInput), load);

// array has to have its own type, can we get rid of this ?
DEFINE_TEMPLATE_TYPE(TSOInputArray);

/** specialisations of arrayObjectCast */
template <> class arrayObjectCast<TSOInput>{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const TSOInput& value){
        IObjectConstSP objValue(IObjectConstSP::attachToRef(&value));
        return objValue;
    }
    /** Casts array element to an IObject */
    static IObjectSP toIObject(TSOInput& value){
        return IObjectSP::attachToRef(&value);
    }
    /** Turns the IObjectSP into an object */
    static TSOInput fromIObject(IObjectSP& value){
        TSOInput *ptr = dynamic_cast<TSOInput *>(value.get());
        if (!ptr){
            throw ModelException("arrayObjectCast::fromIObject", "Object is not"
                                 " a TSOInput");
        }
        return *ptr;
    }
};

void TSOInput::load(CClassSP& clazz)
{
//    clazz->setPublic(); // private data
    REGISTER(TSOInput, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultTSOInput);
    FIELD(isCall, "CALL=true, PUT=false");
    FIELD(strike, "option strike");
    FIELD(maturity, "option maturity date");
    FIELD(isAmerican, "true(default)=american, false=european");
    FIELD(weight, "position size of the option (+ve for long and -ve for short)");
}

///*********************************//
/** VanillaTSO instrument */
///*********************************//
class VanillaTSO: public CInstrument, //public Generic1Factor,
                  public virtual LastSensDate,
                  public virtual ISensitiveStrikes,
                  public FDModel::IIntoProduct, 
                  public Theta::IShift{
public:
    static CClassConstSP const TYPE;

    /** instrument validation */
    virtual void Validate();

    /** retrieve market data needed by VanillaTSO - just valueDate, asset and
        discount yield curve */
    void GetMarket(const IModel*          model, 
                   const CMarketDataSP    market);

    /** create a fd payoff tree - new for all fd/tree state variable interface */
    virtual FDProductSP createProduct(FDModel* model) const;

    /** indicates whether VEGA_MATRIX is sensible for this instrument */
    virtual bool avoidVegaMatrix(const IModel* model){
        return CTree1fLV::TYPE->isInstance(model);
    }

    /** returns all strikes on the vol surface to which 
        this instrument is sensitive */
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*      model);

    /** Returns rolls value date and sets initial spot for Theta,
        return true if sub objects need to be tweaked */
    virtual bool sensShift(Theta* shift);
  
    virtual DateTime getValueDate() const;
 
    virtual DateTime endDate(const Sensitivity* sensControl) const;

    /** price a dead instrument until settlement - exercised, expired, knocked out etc.
    returns true if it is dead (and priced), false if it is not dead */
    virtual bool priceDeadInstrument(CControl* control, CResults* results) const{
        // do nothing for now
        return false;
    }

    // method that create an aggregate instrument
    CInstrument* createAggregateInstrument(
            const CInstrumentArray& instruments,
            const DoubleArray&      weights);

    /** copy results array [1] to [n] */
    void getResults(CResultsArray& results, const CTree1f* model);

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;

private:
    friend class VanillaTSOFDProd;

    VanillaTSO();
    VanillaTSO(const VanillaTSO& rhs);
    VanillaTSO& operator=(const VanillaTSO& rhs);
    static void load(CClassSP& clazz);
    static IObject* defaultVanillaTSO(){
        return new VanillaTSO();
    }

    DateTimeArray getAllMatDates() const;
    LinearStrikeVolRequestSP GetLNRequest() const;

    void addOutputRequests(Control* control,
                           Results* results,
                           const double& fairValue,
                           const double& indVol) const;

    // option spec's
    int             strikeGrid; // size of strike dimension for div-adjustment grid
    bool            tsoAmerican;

    // transient
    DateTime                valueDate;

    InstrumentSettlementSP  instSettle;       /* Instrument settlement details */
    InstrumentSettlementSP  premiumSettle;    /* When premium is paid */

    CAssetWrapper           asset;            /* The underlying */
    YieldCurveWrapper       discount;         /* Ccy to discount payoff */
    
    CashFlowArray           assumedDivs;
    TSOInputArray           instArr;
    IntArray                instIndex; // keeps instrument index (order) after possible re-arraging
    DateTime                lastMatDate;
    int                     NumTSO;
};

// end of intrument class declaration 

void VanillaTSO::Validate()
{
    static const string method = "VanillaTSO::Validate";
    // just check the things that aren't/cannot be checked in 
    if (!asset){
        throw ModelException(method, "Asset is null");
    }
    if (!discount){
        throw ModelException(method, "Discount YC is null");
    }

    AssetUtil::assetCrossValidate(asset.get(),
                                  false,
                                  valueDate,
                                  valueDate,
                                  discount,
                                  this);  

    // remove past dates
    for (int i=0; i<assumedDivs.size(); i++)
    {
        if (assumedDivs[i].date <= valueDate)
        {
            assumedDivs.erase(assumedDivs.begin()+i);
            i--;
        }
    }
    if (strikeGrid !=1)
        strikeGrid =1; // !!! *** use 1 only for now
}

void VanillaTSO::GetMarket(const IModel*          model, 
                         const CMarketDataSP    market)
{
    market->GetReferenceDate(valueDate);

    CAsset::getAssetMarketData(model, market.get(), "V", 
                               discount, asset);

    discount.getData(model, market);
    instSettle->getMarket(model, market.get());
    if (premiumSettle.get())
    {
        premiumSettle->getMarket(model, market.get());
    }
}

DateTime VanillaTSO::getValueDate() const
{
  return valueDate;
}

/** when to stop tweaking */
DateTime VanillaTSO::endDate(const Sensitivity* sensControl) const {
    
    DateTime instEnd  = instSettle->settles(lastMatDate, asset.get());
    DateTime assetEnd = asset->settleDate(lastMatDate);
    DateTime end = instEnd.isGreater(assetEnd) ? instEnd : assetEnd;
    return end;
}

void VanillaTSO::addOutputRequests(Control* control,
                                 Results* results,
                                 const double& fairValue,
                                 const double& indVol) const
{
    // take care of additional outputs
    if ( control && control->isPricing() )
    {
        // DELAY_PRICE
        InstrumentUtil::delayPriceHelper(control,
                                         results,
                                         fairValue,
                                         valueDate,
                                         discount.get(),
                                         asset.get(),
                                         premiumSettle.get());
        
        // FWD_AT_MAT, this is to the lat maturity date
        InstrumentUtil::recordFwdAtMat(control,
                                       results,
                                       lastMatDate,
                                       valueDate,
                                       asset.get());
    }
}

DateTimeArray VanillaTSO::getAllMatDates() const
{
    DateTimeArray dates(instArr.size());
    
    for (int i=1; i<instArr.size(); i++)
        dates[i] = instArr[i].maturity;

    return dates;
}

/** returns a vol request for log-normal vol */
LinearStrikeVolRequestSP VanillaTSO::GetLNRequest() const
{
    // get strike for longest maturity
    DateTime matDate = valueDate;
    double volStrike =0.0;

    for (int i=1; i<instArr.size(); i++)
    {
        if (matDate < instArr[i].maturity)
        {
            matDate = instArr[i].maturity;
            volStrike  = instArr[i].strike;
        }
    }

    DateTime imntStartDate = valueDate;

    LinearStrikeVolRequestSP volRequest(
        new LinearStrikeVolRequest(volStrike, imntStartDate, 
                                   matDate, false));
    return volRequest;
}

// method that create an aggregate instrument
CInstrument* VanillaTSO::createAggregateInstrument(
        const CInstrumentArray& instruments,
        const DoubleArray&      weights)
{
    string method = "VanillaTSO::createAggregateInstrument";

    DividendAdjBase* divAdj = 0;
    CVanilla* van = 0;

    if ( strikeGrid < 1 || strikeGrid > 100)
    {
        throw ModelException(method, "strikeGrid must be between 1 and 100");
    }

    // for re-arrange order, all tso's first
    TSOInputArray  instTSO, instVan;
    IntArray       indexTSO, indexVan;

    NumTSO = 0;

    lastMatDate = valueDate;

    instArr.clear();
    instIndex.clear();

    for (int k=0; k<instruments.size(); k++)
    {
        // for now, only support vanilla and div adjusted types
        if ((divAdj = dynamic_cast<DividendAdjBase*>(instruments[k].get())))
        {
            if (divAdj->ccyTreatment != "N" && divAdj->ccyTreatment != "V")
                throw ModelException(method, "ccy treatment not supported");

            instTSO.resize(NumTSO+1);
            instTSO[NumTSO].isCall = divAdj->isCall;
            instTSO[NumTSO].strike  = divAdj->strike;
            instTSO[NumTSO].maturity  = divAdj->matDate;
            instTSO[NumTSO].isAmerican  = tsoAmerican;
            instTSO[NumTSO].weight  = weights[k];
            NumTSO++; // tso's are div adjusted for now
            indexTSO.push_back(k);
            // copy assumed divs
            if (assumedDivs.size() == 0)
                assumedDivs = divAdj->assumedDivs;
            else if (assumedDivs[assumedDivs.size()-1].date < divAdj->assumedDivs[divAdj->assumedDivs.size()-1].date)
            {// copy more assumed divs
                for (int i=0; i<divAdj->assumedDivs.size(); i++)
                {
                    if (assumedDivs[assumedDivs.size()-1].date < divAdj->assumedDivs[i].date)
                    {
                        assumedDivs.insert(assumedDivs.end()+i, divAdj->assumedDivs.begin()+i, divAdj->assumedDivs.end());
                        break;
                    }
                }
            }
            // get longest maturity
            if (lastMatDate < divAdj->matDate)
                lastMatDate = divAdj->matDate;
        }
        else if ((van = dynamic_cast<CVanilla*>(instruments[k].get())))
        {
            if (van->fwdStarting)
                throw ModelException(method, "fwd start not supported");
            if (van->ccyTreatment != "N" && van->ccyTreatment != "V")
                throw ModelException(method, "ccy treatment not not supported");

            int num = instVan.size();
            instVan.resize(num+1);
            instVan[num].isCall = van->isCall;
            instVan[num].strike  = van->exerciseSchedule->lastValue();
            instVan[num].maturity  = van->exerciseSchedule->lastDate();
            instVan[num].isAmerican  = van->canExerciseEarly;
            instVan[num].weight  = weights[k];
            indexVan.push_back(k);
            // get longest maturity
            if (lastMatDate < instVan[num].maturity)
                lastMatDate = instVan[num].maturity;
        }
        else
        {
            throw ModelException(method, "component instrument not supported by VanillaTSO");
        }
        // fill market data
        if (k==0)
        {
            if (divAdj)
            {
                instSettle = InstrumentSettlementSP(dynamic_cast<InstrumentSettlement*>(divAdj->instSettle->clone()));
                if (divAdj->premiumSettle.get())
                    premiumSettle = InstrumentSettlementSP(dynamic_cast<InstrumentSettlement*>(divAdj->premiumSettle->clone()));
                asset = divAdj->asset;
                discount = divAdj->discount;
            }
            else
            {
                instSettle = InstrumentSettlementSP(dynamic_cast<InstrumentSettlement*>(van->instSettle->clone()));
                if (van->premiumSettle.get())
                    premiumSettle = InstrumentSettlementSP(dynamic_cast<InstrumentSettlement*>(van->premiumSettle->clone()));
                asset = van->asset;
                discount = van->discount;
            }
        }
    }
    // put (possibly re-arranged) data together
    instArr.insert(instArr.end(), instTSO.begin(), instTSO.end());
    instArr.insert(instArr.end(), instVan.begin(), instVan.end());
    instIndex.insert(instIndex.end(), indexTSO.begin(), indexTSO.end());
    instIndex.insert(instIndex.end(), indexVan.begin(), indexVan.end());

    // we create one more inst for the target
    // just copy a dummy there
    instArr.insert(instArr.begin(), *instArr.begin()); 
    instIndex.insert(instIndex.begin(), 0); 
    instArr[0].weight  = 0.0; // but set its weight to 0

    return dynamic_cast<CInstrument*>(clone());
}


/** copy results array [1] to [n] */
void VanillaTSO::getResults(CResultsArray& results, const CTree1f* tree1f)
{
    int i;
    // copy TSO results
    vector<double> v(results.size(),0.0);

    // copy prices to the correct index, in case re-ordered
    for (i=0; i<results.size(); i++)
    {
        v[instIndex[i+1]] = results[i]->retrievePrice();
    }
    // store results
    for (i=0; i<results.size(); i++)
        results[i]->storePrice(v[i], discount->getCcy());
}

// SmoothMax routine
void VanillaTSOSMax(int bot, int top, double strike, int callPutMult, 
                 double settlePV, double* option, const double* under)
{
    const double h = 0.5; // can discuss all day what's best for this one
    double diffUp = 0;
    double valueDown, diffDown, hMaxDiff;
    double valueMid = option[top] - under[top]*callPutMult;
    for (int idx=top; idx > bot; idx--)          
    {
        valueDown = option[idx-1] - under[idx-1]*callPutMult;
        diffDown = fabs(valueMid - valueDown);
        hMaxDiff = h*Maths::max(diffUp, diffDown);
        option[idx] = SMax((under[idx]-strike)*callPutMult, option[idx], hMaxDiff);
        diffUp = diffDown;
        valueMid = valueDown;
    }
    hMaxDiff = diffUp;
    option[bot] = settlePV*SMax((under[bot]-strike)*callPutMult, option[bot], hMaxDiff);
}


/******************************************************************************************************************************/

class VanillaTSOFDProd : public LatticeProdEDR
{
public:
    VanillaTSOFDProd(const VanillaTSO* VanillaTSO, FDModel * model) :
        LatticeProdEDR(model), inst(VanillaTSO)
    {
        if( ! tree1f )
        {
            throw ModelException( "VanillaTSOFDProd::VanillaTSOFDProd", "Instrument of type "+
                                 inst->getClass()->getName() +
                                 " can be priced by CTree1f only" );
        }

        // first: set discount curve
        if( tree1f )
            tree1f->setDiscountCurve( inst->discount.getSP() );

        // second: create spot payoff
        payoffIndex = model->createProduct( IProdCreatorSP( new
            IndexSpecEQ( inst->asset.getName(), inst->asset, "V" ) ) );
    }

    /** initialise model1f - allow product customisation */
    virtual void init(CControl* control) const;

    virtual DateTime getStartDate() const 
    {
        return inst->valueDate;
    }

    /** premium scaling */ 
    double scalePremium(const double& fairValue, 
                              YieldCurveConstSP disc)                              
    {    
        return fairValue;
    }

    /** initialising and setting product variables */
    // this is called per pricing call before tree sweep call (after InitTree)
    virtual void initProd();

    void recordOutput(Control* control, 
                      YieldCurveConstSP disc, 
                      Results* results);

    /** product payoff method at maturity */
    void prod_BWD_T(const TreeSlice & spot,
                          int step,
                          int bot,
                          int top,
                          int pStart,
                          int pEnd,
                          const vector< TreeSliceSP > & price);

    // this is payoff boundary condition, for KO, early exercise etc.
    void prod_BWD(const TreeSlice & spot,
                        int step,
                        int bot,
                        int top,
                        int pStart,
                        int pEnd,
                        const vector< TreeSliceSP > & price);    

    virtual void update(int& step, FDProduct::UpdateType type);

    // temp here 
    void AdjustDeltaShift(CControl* control) const 
    {/* does nothing*/}

    virtual bool GetFwdStartLV() 
    {
        return false;
    }

    virtual DateTime GetFwdStartDateLV() 
    {
        return inst->valueDate;
    }

    // aggregate component prices
    void AggregateValue(const double*s, 
                              int step, 
                              int bot, 
                              int top, 
                              int pStart, 
                              int pEnd, 
                              const vector< double * > & p);

    // maturity payoff helper
    inline void payoffMat(const double* s, 
                                int iNum, 
                                int bot, 
                                int top, 
                                double df, 
                                const vector< double * > & p);


private:
    const VanillaTSO*   inst;
   
    // can exercise flag for each step, for isAmerican=true
    vector<bool>    stepCanExercise;
    vector<double>  strikes;
    vector<int>     matStep;
    DateTimeArray   adjustmentDates;
    vector<bool>     adjustStep;
    vector<double>   divDiff;
};

FDProductSP VanillaTSO::createProduct(FDModel* model) const
{
    return FDProductSP( new VanillaTSOFDProd(this, model) );
}

/** initialise tree1f - allow product customisation */
void VanillaTSOFDProd::init(CControl* control) const
{
    static const string method = "VanillaTSOFDProd::Init";
    try 
    {
        if( tree1f )
        {
            if (tree1f->GetSmoothMethod() == CTree1f::DEFAULT) 
                tree1f->SetSmoothMethod(CTree1f::NODE_INSERTION);
            
            tree1f->NumOfPrice = inst->NumTSO * inst->strikeGrid + inst->instArr.size() - inst->NumTSO;
            tree1f->NumOfInsertNode = 0;                
        }

        /** customize tree parameters here and set up the tree */
        DateTimeArray segDates;
        segDates.resize(2);

        segDates[0] = inst->valueDate; // no fwd start supported
        segDates[1] = inst->lastMatDate;

        IntArray density( 1, 1 );

        // all exercise dates are copied to critical dates
        DateTimeArray critDates = inst->getAllMatDates();
        
        // remove last exercise date from crit date, as it's there
        critDates.erase(critDates.end()-1);
        // add div event dates if needed
        EventAssetMove divEvent;
        DateTimeArraySP divCritDates;

        DateTime end(segDates[1]);
        if (tree1f->DivsTreatedAsAbsolute()) 
        {
            end = max(end, Equity::calcDivTransPeriodEndDate(inst->valueDate));
        }
        //int numDivs = (int)(4*start.yearFrac(end))+1; // 4 divs per year selected as critical dates
        if (AssetUtil::getJumpEvents(inst->asset.get(), 
                                     inst->valueDate, 
                                     inst->lastMatDate,
                                     1000000, // get them all
                                     divEvent)) 
        {

            // calculate critical dates
            divCritDates = divEvent.getCritDate(0, true/*inst->isCall*/);
            /* If dol divs, filter past critical dividend dates out as they are used as bench 
               dates by pseudo asset. Have to do this here (as opposed to in Equity), since
               if there are no critical div dates left after filtering I revert to 
               DivAmountTreatment == false immediately below */
            if (tree1f->DivsTreatedAsAbsolute() && divCritDates->size() > 0)
            {
                DateTimeArray temp(segDates[0].getFutureDates(*divCritDates));
                divCritDates = DateTimeArraySP(copy(&temp));
            }
        }

        /* If no dividend critical dates, no point to bother with dol divs */
        if (tree1f->DivsTreatedAsAbsolute() && (!divCritDates || divCritDates->size() == 0))
        {
            tree1f->SetDivAmountTreatment(false);
        }

        // use simple delta size adjustment if needed, note that TreeDeltaShift is stored in base tree only.
        if (!tree1f->DEBUG_SameGridDelta && control->isPricing()) 
        {
            AdjustDeltaShift(control);
        }
        
        // add critical dates
        model->addCritDates( critDates );
        if( divCritDates.get() )
            tree1f->addDivCritDates( *divCritDates );

        // prepare timeline set up
        model->initSegments( segDates, density );
    }
    catch (exception& e) 
    {
        throw ModelException(e, method);
    }
}

/** initialising and setting product variables
    this is called per pricing call before tree sweep call (after InitTree)  */
void VanillaTSOFDProd::initProd() 
{
    static const string method = "VanillaTSOFDProd::InitProd";
    
    try 
    {
        int i,j;
        double fwdAtStart = 1.0;

        initSlices(inst->NumTSO * inst->strikeGrid + inst->instArr.size() - inst->NumTSO );
 
        stepCanExercise.resize(model->getLastStep() + 1);
        // ask tree to decide first about steps that can exercise
        ScheduleSP exerciseSchedule(new Schedule(DateTimeArray(1,inst->lastMatDate),
                                                DoubleArray(1, 0.0),
                                                "L"));
       
        AssetUtil::setStepExercise(stepCanExercise,
                                 model->getDates(),
                                 exerciseSchedule,
                                 true,
                                 inst->asset.getSP());

        // signal ex-date - 1 step for strike adjustment
        adjustStep.resize(model->getLastStep() + 1, false);
        divDiff.resize(model->getLastStep() + 1, 0.0);
        // set up strikes
        int numOfPriceArray = inst->NumTSO * inst->strikeGrid + inst->instArr.size() - inst->NumTSO; 
        strikes.resize(numOfPriceArray);
        if (inst->NumTSO >0)
        {
            DividendListSP divs = AssetUtil::getAllDivsBetweenDates(inst->asset,
                                                                    inst->valueDate,
                                                                    inst->lastMatDate);
            DoubleArrayConstSP amounts = divs->getDivAmounts();

            EventAssetMove divEvent;
            if (AssetUtil::getJumpEvents(inst->asset.get(), 
                                         inst->valueDate, 
                                         inst->lastMatDate,
                                         1000000, // get them all
                                         divEvent)) 
            {
                // copy critical dates
                adjustmentDates = (*divEvent.getCritDate(0, true/*inst->isCall*/));
            }
            if (adjustmentDates.size() != amounts->size())
            {
                throw ModelException(method, "number of divs = " + Format::toString(amounts->size()) +
                                " does not equal to num of adjustment dates = " + Format::toString(adjustmentDates.size()));
            }

            vector<double> div_diff(amounts->size());
            double minDiff = 0.0;
            double maxDiff = 0.0;
            double netDiff = 0.0;
            for (j=0; j<amounts->size(); j++)
            {
                div_diff[j] = (*amounts)[j] - inst->assumedDivs[j].amount;
                netDiff += div_diff[j];
                minDiff = Maths::min(Maths::min(netDiff, div_diff[j]), minDiff);
                maxDiff = Maths::max(Maths::max(netDiff, div_diff[j]), maxDiff);
            }

            // signal ex-date - 1 step for strike adjustment
            if (adjustmentDates.size() != 0) {    
                for (i=0, j=0; i <= model->getLastStep(); i++)
                {
                    if (adjustmentDates[j] == model->getDate(i))
                    {
                        adjustStep[i] = true;
                        divDiff[i] = div_diff[j];
                        j++;
                        if (j == adjustmentDates.size())
                            break;
                    }
                }
                if (adjustmentDates.size() != j)
                {
                    throw ModelException(method, "number of adjustment dates = " + Format::toString(adjustmentDates.size()) +
                                    " does not equal to num of adjustment steps = " + Format::toString(j));
                }
            }
            // tso strikes
            for (i=1; i<=inst->NumTSO; i++)
            {
                double adjAmount = 0.0;
                if (inst->strikeGrid > 1)
                    adjAmount  =(maxDiff - minDiff)/(inst->strikeGrid - 1.0);

                if (adjAmount < FP_MIN)
                    adjAmount = FP_MIN;

                for (j=0; j<inst->strikeGrid; j++)
                {// allow sqrt(n) for grid range
                    if (inst->strikeGrid > 1)
                        strikes[(i-1)*inst->strikeGrid+j+1] = inst->instArr[i].strike 
                                                - sqrt((double)(inst->strikeGrid))*(minDiff + j*adjAmount); 
                    else
                        strikes[(i-1)*inst->strikeGrid+j+1]  = inst->instArr[i].strike;
                }
            }
        }

        // hedge option strikes
        j = inst->NumTSO*inst->strikeGrid - inst->NumTSO;
        for (i=inst->NumTSO+1; i<inst->instArr.size(); i++)
        {
            strikes[i+j] = fwdAtStart*inst->instArr[i].strike;
        }

        // signal maturity step
        matStep.resize(inst->instArr.size(), 0);
        for (j=1; j< (int)matStep.size(); j++) // skip the first aggregate inst
        {
            for (i=0; i <=  model->getLastStep(); i++)
            {
                if (inst->instArr[j].maturity == model->getDate(i))
                {
                    matStep[j] = i;
                    i = 0; // re-start again
                    break;
                }
            }
        }
    }
    catch (exception& e) 
    {
        throw ModelException(e, method);
    }
}

void VanillaTSOFDProd::prod_BWD_T(const TreeSlice & spot,
                                        int step, 
                                        int bot, 
                                        int top, 
                                        int pStart, 
                                        int pEnd,
                                        const vector< TreeSliceSP > & price) 
{
    double * s = spot.getValues();
    const vector< double * > & p = getValues( price );

    double settlementPV = inst->instSettle->pvAdjust(inst->lastMatDate,
                                                     inst->discount.get(), 
                                                     inst->asset.get());
    int i,j;
    // init to 0's first
    for (i=pStart; i<=pEnd; i++)
    {
        for (j=bot; j<=top; j++)
            (p[i])[j] = 0.0;
    }
    // skip the 1st component which is used as target here
    // compute mat value
    for (i=1; i<inst->instArr.size(); i++)
    {
        if (matStep[i] == step)
            payoffMat(s, i, bot, top, settlementPV, p);
    }
    // add all together to [0] array
    AggregateValue(s, step, bot, top, pStart, pEnd, p);
}

void VanillaTSOFDProd::prod_BWD(const TreeSlice & spot,
                                        int step, 
                                        int bot, 
                                        int top, 
                                        int pStart, 
                                        int pEnd,
                                        const vector< TreeSliceSP > & price)
{
    static const string method = "VanillaTSOFDProd::prod_BWD";
    try {
        double * s = spot.getValues();
        const vector< double * > & p = getValues( price );

        int i,j, k;
        int numGrid, gridStart;
        bool fastRoll;
        vector <double> vol_arr, drift_arr;

        double settlementPV = inst->instSettle->pvAdjust(model->getDate(step),
                                                         inst->discount.get(), 
                                                         inst->asset.get());
        // loop all options
        // skip the 1st in aggregate which is the target
        for (i=1; i<inst->instArr.size(); i++)
        {
            if (matStep[i] >= step)
            {// only if maturity is in play
                if (i <= inst->NumTSO)
                {
                    numGrid = inst->strikeGrid;
                    gridStart = (i-1)*inst->strikeGrid +1;
                }
                else
                {
                    numGrid = 1;
                    gridStart = inst->NumTSO*inst->strikeGrid - inst->NumTSO + i;
                }
                if (matStep[i] == step)
                {
                    payoffMat(s, i, bot, top, settlementPV, p);
                }
                // do penultimate smoothing
                if (tree1f && step == matStep[i]-1)
                {
                    if (vol_arr.size() == 0)
                    {
                        fastRoll = tree1f->CalcStepDriftAndVar(s, bot, top, vol_arr, &drift_arr); // just to get drift
                        // use GetStepVol() to get vol, do not use GetStepVar() which may not have simple BS conversion
                        tree1f->GetStepVol(step, vol_arr, s, bot, top);
                    }

                    double dt = tree1f->getTrdYrFrac(step+1);
                    double variance = vol_arr[0]*vol_arr[0]*dt;
                    double drift = drift_arr[0];
                    double df = inst->discount->pv(model->getDate(step),
                                                   model->getDate(step+1));

                    /* If dollar dividend treatment, we must adjust fwd and strike by stock floor.
                       If not dollar dividend treatment, stock floor == 0.0 */
                    double StockFloor0 = tree1f->GetStockFloor(step);
                    double StockFloor1 = tree1f->GetStockFloor(step + 1);

                    // loop each grid
                    for (k=0; k<numGrid; k++)
                    {
                        for (j=bot; j<=top; j++)
                        {
                            if (!fastRoll)
                            {
                                if (vol_arr.size() > 1)
                                { 
                                    variance = vol_arr[j - bot]*vol_arr[j - bot]*dt; // needs one vol per node - local vol
                                }
                                if (drift_arr.size() > 1)
                                {
                                    drift = drift_arr[j - bot]; // needs drift per node
                                }
                            }

                            if (s[j]>0.0 && fabs(log(s[j]/strikes[gridStart+k])) < tree1f->TruncationStd*sqrt(variance))
                            {
                                /* If dollar div treatment (in which case stock floor is <> 0.0), the price of the option
                                   is that of an option with same maturity written on the pseudo asset (i.e., S - StockFloor)
                                   and with strike adjusted by the stock floor at maturity (i.e., K - StockFloor).
                                   NB Black::price takes care of potentially negative fwd and strike values. */
                                (p[gridStart+k])[j] = Black::price(inst->instArr[i].isCall, 
                                                             drift * (s[j] - StockFloor0),  // pseudo asset's (step+1)-maturity fwd
                                                             // as viewed from t = step
                                                             strikes[gridStart+k] - StockFloor1,  // adjusted strike
                                                             df,
                                                             variance);
                            }
                        }
                    }
                }
                // takes care of exercise
                if (inst->instArr[i].isAmerican)
                {
                    double intrinsic;
                    double callput = (inst->instArr[i].isCall? 1.0:-1.0);
                    // loop each grid
                    for (k=0; k<numGrid; k++)
                    {
                        if (tree1f->GetSmoothMethod() == CTree1f::SMAX)
                        {
                            VanillaTSOSMax(bot, top, strikes[gridStart+k], (int)callput,
                                            settlementPV, (p[gridStart+k]), s);
                        }
                        else
                        {
                            for (j=bot; j<=top; j++)
                            {
                                intrinsic = settlementPV*callput*(s[j] - strikes[gridStart+k]);
                                if ((p[gridStart+k])[j] < intrinsic)
                                {
                                    (p[gridStart+k])[j] = intrinsic; // American
                                }
                            }
                        }
                    }
                }
            }
        }
        // add all together to [0] array
        AggregateValue(s, step, bot, top, pStart, pEnd, p);
    }
    catch (exception& e) 
    {
        throw ModelException(e, method);
    }
}

// perform adjustment for div-protection
// aggregate component prices at each node for [0] array for use with gamma scaling
void VanillaTSOFDProd::AggregateValue(const double*s, 
                                            int step, 
                                            int bot, 
                                            int top, 
                                            int pStart, 
                                            int pEnd, 
                                            const vector< double * > & p)
{
    int j, k, i, n;
    bool interp;

    if (adjustStep[step])
    {// strike adjustment step for TSO
        int grid = inst->strikeGrid;
        vector<double> adjPrice(grid);
        vector<double> v(grid);
        for (k=1; k<=inst->NumTSO; k++)
        {
            for (j=bot+1; j<=top; j++) // no adjustment for bottom node (s=0)
            {
                n = (k-1)*grid + 1;
                // copy prices for interp
                for (i = 0; i<grid; i++)
                {
                    v[i] = (p[n+i])[j];
                    double delta_s = ((p[n+i])[j] - (p[n+i])[j-1])/(s[j] - s[j-1]);
                    interp = (fabs(delta_s) > FP_MIN);
                    if (!interp)
                        break;

                    adjPrice[i]= (p[n+i])[j] + delta_s * divDiff[step];
                }
                // re-interp strike and price if needed
                if (interp)
                {
                    if (grid == 1)
                        (p[n])[j] = Maths::max(0.0, adjPrice[0]);
                    else
                    {// interpolate adjusted strikes crossing div date, and then re-interp for values of original strikes
                        interpStrikeAdj(&*strikes.begin()+n, adjPrice, v);
                        // copy results
                        for (i = 0; i<grid; i++)
                        {
                            if(v[i] <= 0.0) // option payoff assumed !!!
                                (p[n+i])[j] = 0.0;
                            else
                                (p[n+i])[j] = v[i];
                        }
                    }
                }
            }
        }
    }
    
    if (step >0)
    {// sum all values to [0] for gamma-vol scaling
        for (j=bot; j<=top; j++)
        {
            (p[0])[j] = 0.0;
            // process TSO's first
            for (k=pStart+1; k<=inst->NumTSO; k++)
            {
                int n = (k-1)*inst->strikeGrid + 1 + inst->strikeGrid/2; // choose the middle grid point
                (p[0])[j] += inst->instArr[k].weight*(p[n])[j];
            }
            //for vanilla's
            for (k=inst->NumTSO*inst->strikeGrid+1; k<=pEnd; k++)
            {
               (p[0])[j] += inst->instArr[k].weight*(p[k])[j];
            }
        }
    }
    else if (inst->strikeGrid>1) // copy results into 1st numOfPrice array for getResults() use
    {// this is step 0
        // process TSO's first
        for (k=1; k<=inst->NumTSO; k++)
        {
            unsigned long loc;
            int n = (k-1)*inst->strikeGrid; // remember we skip 1st target and NR has offset 1
            locate(&*strikes.begin()+n, inst->strikeGrid, inst->instArr[k].strike, &loc);
            if ((int)loc == inst->strikeGrid) loc--; // outside upper, shift down
            // linear interp, n is offset 1 but this is ok as we skip first target component
            for (j=bot; j<=top; j++)
            {
                (p[k])[j] = (p[n+loc])[j]+(inst->instArr[k].strike-strikes[n+loc])
                                *((p[n+loc+1])[j]-(p[n+loc])[j])/(strikes[n+loc+1]-strikes[n+loc]);
            }
        }
        //for vanilla's
        for (k=1; k<=pEnd-inst->NumTSO*inst->strikeGrid; k++)
        {
            for (j=bot; j<=top; j++)
            {
                (p[k])[j] = (p[k+inst->NumTSO*inst->strikeGrid])[j];
            }
        }
    }
}

inline void VanillaTSOFDProd::payoffMat(const double* s, 
                                              int iNum, 
                                              int bot, 
                                              int top, 
                                              double df, 
                                              const vector< double * > & p)
{
    int j, k;
    int numGrid, gridStart;
    if (iNum <= inst->NumTSO)
    {
        numGrid = inst->strikeGrid;
        gridStart = (iNum-1)*numGrid+1;
    }
    else
    {
        numGrid = 1;
        gridStart = inst->NumTSO*numGrid - inst->NumTSO + iNum;
    }

    for (k=0; k<numGrid; k++)
    {
        for (j=bot; j<=top; j++)
        {
            (p[gridStart+k])[j] = df* GetIntrinsic(s[j], strikes[gridStart+k], inst->instArr[iNum].isCall, true/* this allow fwd */);
        }
    }
}

void VanillaTSOFDProd::update(int& step, FDProduct::UpdateType type)
{
    // we assume just need one und level for spot here
    const TreeSlice & s = payoffIndex->getValue( step );
    int bot, top;
    s.getCalcRange( bot, top );

    const vector< TreeSliceSP > & price = slices;
    int pStart = 0, pEnd = price.size() - 1;

    if (type == FDProduct::BWD_T)
    {
        prod_BWD_T(   s,
                      step,
                      bot,
                      top,
                      pStart, 
                      pEnd,
                      price);       
    }
    else if(type == FDProduct::BWD)
    {
        prod_BWD( s,
                  step,
                  bot,
                  top,
                  pStart, 
                  pEnd,
                  price);      
    }
}

/** returns all strikes on the vol surface to which 
    this instrument is sensitive */
DoubleArraySP VanillaTSO::getSensitiveStrikes(OutputNameConstSP outputName,
                                            const IModel*      model)
{
    DoubleArraySP sensStrikes = DoubleArraySP(new DoubleArray(0));

    if (avoidVegaMatrix(model)) {
        throw ModelException("VanillaTSO::getSensitiveStrikes", 
                             "VEGA_MATRIX is not valid for this instrument");
    }

    // create a vol request object to be passed on
    LinearStrikeVolRequestSP volRequest = GetLNRequest();

    SensitiveStrikeDescriptor sensStrikeDesc;
    sensStrikeDesc.forwardOnly = false;

    asset->getSensitiveStrikes(volRequest.get(), outputName, 
                               sensStrikeDesc, sensStrikes);

    return sensStrikes;
}

/** Rolls the value date and sets initial spot if rolling over start date */
bool VanillaTSO::sensShift(Theta* shift)
{    
    // roll today 
    valueDate = shift->rollDate(valueDate);
    
    return true;
};


/** Returns the name of the instrument's discount currency */
string VanillaTSO::discountYieldCurveName() const {
    return discount.getName();
}


// for reflection
VanillaTSO::VanillaTSO(): CInstrument(TYPE), strikeGrid(1), tsoAmerican(false), NumTSO(0){};

void VanillaTSO::load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VanillaTSO, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(CTree1f::IIntoProduct);
        IMPLEMENTS(FDModel::IIntoProduct);
        IMPLEMENTS(ISensitiveStrikes);
        IMPLEMENTS(Theta::Shift);
        IMPLEMENTS(LastSensDate);
        EMPTY_SHELL_METHOD(defaultVanillaTSO);
        FIELD(strikeGrid,"strike grid size, between 1 and 10");
        FIELD_MAKE_OPTIONAL(strikeGrid);
        FIELD(tsoAmerican,"true=TSO options are american exercise ");
        FIELD_MAKE_OPTIONAL(tsoAmerican);

        // transient
        FIELD(valueDate,"valuation Date");
        FIELD_MAKE_TRANSIENT(valueDate);
        FIELD(instSettle, "Instrument settlement at maturity");
        FIELD_MAKE_TRANSIENT(instSettle);
        FIELD(premiumSettle, "Premiumsettlement");
        FIELD_MAKE_TRANSIENT(premiumSettle);
        FIELD(asset,"Underlying of option");
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(asset);
        FIELD(discount,"Discount curve");
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(discount);

        FIELD(assumedDivs,"assumed dividends");
        FIELD_MAKE_TRANSIENT(assumedDivs);
        FIELD(instArr,"array of options(of TSOInput type)");
        FIELD_MAKE_TRANSIENT(instArr)
        FIELD(instIndex,"");
        FIELD_MAKE_TRANSIENT(instIndex)
        FIELD(lastMatDate,"");
        FIELD_MAKE_TRANSIENT(lastMatDate);
        FIELD(NumTSO,"");
        FIELD_MAKE_TRANSIENT(NumTSO);
    }

CClassConstSP const VanillaTSO::TYPE = CClass::registerClassLoadMethod(
    "VanillaTSO", typeid(VanillaTSO), load);

/** Combine VanillaTSO aggregate type with Tree1fLVGD **
    If EAS can take separate params it would be much better not to do this. 
    e.g. VanillaTSO cannot be used with log-normal tree and each product needs 
    to be added into the combined model for support */

// put two together for now

typedef smartPtr<VanillaTSO>  VanillaTSODSP;

class ModelTSO : public CTree1fLVGD, public IAggregateModel
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);

    static IObject* defaultModelTSO(){
        return new ModelTSO();
    }

    ModelTSO() : CTree1fLVGD(TYPE) {
        target= VanillaTSODSP(   );
    }

    virtual CInstrument* createAggregateInstrument(
        const CInstrumentArray& instruments,
        const DoubleArray&      weights)
    {
        if (target.get() == 0)
            throw ModelException("ModelTSO", "target not supplied.");

        return target->createAggregateInstrument(instruments, weights);
    }
    
    virtual void price(CInstrument*         tso, 
                       Control*             control, 
                       CResultsArray&       results)
    {
        priceEnd = CResultsArraySP::attachToRef(&results);
        // call base method
        CTree1fLVGD::Price(tso, control, results[0].get());
        dynamic_cast<VanillaTSO*>(tso)->getResults(results, this);
    }

    VanillaTSODSP   target;
    CResultsArraySP priceEnd; // recordOutput() to store results 
};

CClassConstSP const ModelTSO::TYPE = CClass::registerClassLoadMethod(
    "ModelTSO", typeid(ModelTSO), load);

void ModelTSO::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(ModelTSO, clazz);
    SUPERCLASS(CTree1fLVGD);
    IMPLEMENTS(IAggregateModel);
    EMPTY_SHELL_METHOD(defaultModelTSO);
    FIELD(target, "aggregate target type");
    FIELD(priceEnd, "");
    FIELD_MAKE_TRANSIENT(priceEnd);
}

/** record output */
void VanillaTSOFDProd::recordOutput(Control* control, 
                                    YieldCurveConstSP disc, 
                                    Results* results)
{
    // get prices at t=0
    double price0 = model->getPrice0( *slices[0] );
    // save price
    double price = scalePremium(price0, disc);
    results->storePrice(price, disc->getCcy());

    ModelTSO* m = dynamic_cast<ModelTSO*>(tree1f);
    if(m){// store component prices, skip the first one which is the aggregate
        double p;
        for (int i=1; i<m->NumOfPrice; i++){
            p = scalePremium(model->getPrice0(*slices[i]) , disc);
            if (m->priceEnd->size() < m->NumOfPrice-1)
                double x = 1.0;
            (*(m->priceEnd))[i-1]->storePrice(p, disc->getCcy());
        }
    }

    // take care of additional outputs
    if ( control && control->isPricing() )
    {
        DateTime       matDate = inst->lastMatDate;
        double         indVol;
        // calculate indicative vol
        if ( matDate.isGreater(inst->valueDate) )
        {
            DateTime imntStartDate = inst->valueDate;

            // get vol request
            CVolRequestConstSP lnVolRequest = inst->GetLNRequest();

            try{
                // interpolate the vol
                CVolProcessedSP  vol(inst->asset->getProcessedVol(lnVolRequest.get()));
                // cast to the BS vol we're expecting
                CVolProcessedBSSP volBS = CVolProcessedBSSP::dynamicCast(vol);

                // this should never happen if our get market data has worked properly
                if (!vol){
                    throw ModelException("VanillaTSOFDProd::recordOutput", 
                                         "No Black Scholes Vol");
                }
                // calculate the indicative vol
                indVol = volBS->CalcVol(imntStartDate, matDate);
            }
            catch (exception& ) {
                indVol = 0.0;
            }
        }
        else
        {
            indVol = 0.0;
        }

        inst->addOutputRequests(control,
                                results,
                                price,
                                indVol);
    }
}

extern bool VanillaTSOLoad()
{
    return VanillaTSO::TYPE && ModelTSO::TYPE;
}

DRLIB_END_NAMESPACE

