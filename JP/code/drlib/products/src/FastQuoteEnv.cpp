//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FastQuoteEnv.cpp
//
//   Description : pricing environment for market making and exchange linked trading.
//                 performance is the key requirement to be close to real time..
//                 the pricing environment performs optimisation and contains cached data.
//
//   Author      : Ning Shen
//
//   Date        : 24 June 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Nrfns.hpp"
#include "edginc/FastQuoteEnv.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/ImpliedVol.hpp"
#include "edginc/VegaParallel2Sided.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/VegaParallel.hpp"

DRLIB_BEGIN_NAMESPACE

/******* local helper  *************/
//** set vol  */
static void setVol(CAsset* asset, double vol)
{
    EquityCache* assetCache = dynamic_cast<EquityCache*>(asset);
    if (assetCache) // this is american case which has cache for now
        assetCache->setVol(vol);
    else
    {
        // temp solution for euro env
        EquityBase* eq = dynamic_cast<EquityBase*>(asset);
        VolLevel shift(vol);
        VolLevel::Shift* tweak = dynamic_cast<VolLevel::Shift*>(const_cast<CVolBase*>(eq->getVol().get()));
        tweak->sensShift(&shift);
    }
}

static void setSpot(CAsset* asset, double spot)
{
    EquityCache* assetCache = dynamic_cast<EquityCache*>(asset);
    if (assetCache) // this is american case which has cache for now
        assetCache->setSpot(spot);
    else
    {
        // temp solution for euro env
        EquityBase* eq = dynamic_cast<EquityBase*>(asset);
        SpotLevel newSpot(spot);
        eq->getEquity()->sensShift(&newSpot);
    }
}

static double fPrice(double vol, CInstrument* inst, CAsset* asset, IModel* model)
{
    CResults result;
    setVol(asset, vol);
    model->Price(inst, 0, &result);
    return result.retrievePrice();
}


static DateTime lastDivDate(const EquityBase* eq)
{
    // get the 5yr transition date
    DateTime divEnd = eq->getEquity()->getDivTransPeriodEndDate();
    const DividendArray& divs = eq->getDivList()->getArray();
    DateTime maxDollarDiv;

    // find last $div ex-date < divEnd
    for (int i = 0; 
         i < divs.size() && divs[i].getExDate() < divEnd; 
         i++) {
        if (divs[i].getExDate() > maxDollarDiv &&
            divs[i].getDivType() == Dividend::AMOUNT) {
            maxDollarDiv = divs[i].getExDate();
        }
    }
    
    return (divEnd < maxDollarDiv ? divEnd : maxDollarDiv);
}

//******* end of local helper  *************/

///// FastImpVol class members ///////////////////
///////////////////////////////////////////////////
// FastImpVol  class
///////////////////////////////////////////////////
/** compute implied vol. returns true if result converged, false if not */
bool FastImpVol::FastImplied(double spot, double price, double volGuess, 
                             CInstrument* inst, IModel* model,
                             CAsset* asset, CResults& out)
{
    const double tol = 0.0001;
    const double vega_tiny = 1.0e-6*asset->getSpot();
    const double vol_max = 5.0;
    const double volGuess_default = 0.3;
    const int max_iter = 5; // max num of iterations

    int i;
    double price2;
    double price1 = oldPrice;
    double vol1 = oldVol;
    double vol2 = 0.0;
    double vega = oldVega;
    bool  pricedOnce = false;

    try{
        // first estimate using old vega and price
        if (price1 > 0.0 && vega > 0.0)
            vol2 =  Maths::max(1e-8, vol1+(price-price1)/vega);
        else
        {
	    // 0% is not a good guess then use the default guess
            vol1 = Maths::isZero(volGuess) ? volGuess_default : volGuess;
	    price1 = fPrice(vol1, inst, asset, model);
            vol2 = vol1*price/price1;
            if (vol2 > vol_max)
                vol2 = vol_max;
            pricedOnce = true;
        }

        // calc if outside tolerance
        if (fabs(vol2 - vol1) > tol ||
            (fabs(oldSpot/spot -1.0) > tol && !pricedOnce))
        {
            // loop up to 5 times, simple Newton
            int n_iter = 0;
            for (i=0; i<max_iter && fabs(vol2 - vol1) > tol; i++)
            {
                price2 = fPrice(vol2, inst, asset, model);
                vega = (price2 - price1)/(vol2-vol1);
                if (vega < vega_tiny)
                    return false; // does not handle tiny vega case
                vol2 = (vol1=vol2) + (price - (price1=price2))/vega;
                vol2 = Maths::min(Maths::max(0.5*vol1, vol2), 1.5*vol1); // restrict size of change
                n_iter ++;
            }
            // store num of iterations used in debug
            out.storeScalarGreek(n_iter, Results::DEBUG_PACKET,
                                 OutputNameSP(new OutputName("ITERATIONS_FAST_IV")));
        }
    }
    catch (exception& ){
        return false;
    }

    if (fabs(vol2 - vol1) > tol) // not converged
        return false;

    // cache results
    oldVol = vol2;
    oldPrice = price;
    oldVega = vega;
    oldSpot = spot;

    OutputNameSP valueOutName(new OutputName("VALUE"));
    out.storeScalarGreek(vol2, "IMPLIED_VOL", valueOutName);
    return true;
}

///////////////////////////////////////////////////////
// FastQuoteSpline class 
////////////////////////////
bool FastQuoteSpline::canInterpolate(const string& req, double spot, double vol)
{
    return (spot >= spotBot && spot <= spotTop
            && fabs(volCached - vol) < 0.00001);
}

/** local helper function. returns
    -1 not to include this result
    0 this is price
    1 this is a greek  */
inline int checkResult(string& packetName, OutputNameSP& outName)
{
    int check = -1;
    
    if (packetName == Results::INSTRUMENT_PACKET && outName->equals(Results::VALUE))
        check = 0;
    else if (packetName == Results::FWD_AT_MAT_PACKET
             || packetName == Delta::NAME
             || packetName == Delta::SECOND_ORDER_NAME
             || packetName == RhoParallel::NAME
             || packetName == RhoBorrowParallel::NAME
             || packetName == VegaParallel::NAME
             || packetName == VegaSkewParallel::NAME
             || packetName == VegaParallel2Sided::NAME
             || packetName == VegaParallel2Sided::SECOND_ORDER_NAME
             || packetName == DDeltaDVol::NAME
             || packetName == Results::INSTRUMENT_PACKET && outName->equals(OutputRequest::IND_VOL)
             || packetName == Results::INSTRUMENT_PACKET && outName->equals(Theta::NAME)
             || packetName == Theta::NAME
             || packetName == MuParallel::NAME)
        check = 1;

    return check;;
}

/** store and prepare spline interpolation */
void FastQuoteSpline::prepareSpline(const string& inFlags, double vol, const CResultsArray& in)
{
    // for spline
    const double y1 = 2e30; // default to natural spline
    const double yn = 2e30;

    int i,j;
    int resultType;

    // store these for checking if can interp
    resultFlag = inFlags;
    resultCcy = in[0]->getCcyName();
    volCached = vol;
    
    DoubleArray     results;
    // retrieve all results of the first grid point
    int numResults = in[0]->retrieveAllScalarGreeks(packetNames, outNames, results);
    int numGrid = in.size();
    // prepare cubic spline coeff
    for (i=0; i<numResults; i++)
    {
        if((resultType=checkResult(packetNames[i], outNames[i]))>=0)
        {
            y2.resize(i+1);
            valueGrid.resize(i+1);
            y2[i].resize(numGrid);
            valueGrid[i].resize(numGrid);
            for (j=0; j<numGrid; j++)
            {
                if (resultType == 0)
                    valueGrid[i][j] = in[j]->retrievePrice();
                else
                    valueGrid[i][j] = in[j]->retrieveScalarGreek(packetNames[i], outNames[i]);
            }
            spline(&*spotGrid.begin()-1, &*valueGrid[i].begin()-1, 
                   numGrid, y1, yn, &*y2[i].begin()-1);
        }
        else
        {// remove unused names
            packetNames.erase(packetNames.begin()+i);
            outNames.erase(outNames.begin()+i);
            i --; 
            numResults --;
        }
    }

    // set interp limits
    spotBot = spotGrid[0];
    spotTop = spotGrid[numGrid-1];
}

/** interpolate a price  */
void FastQuoteSpline::interpolate(double spot, CResults& out)
{
    double result;
    for (int i = 0; i < outNames.size(); i++)
    {
        splint(&*spotGrid.begin()-1, &*valueGrid[i].begin()-1,
               &*y2[i].begin()-1, 
               spotGrid.size(), spot, &result);

        // store result
        if(checkResult(packetNames[i], outNames[i]) == 0)
        {
            out.storePrice(result, resultCcy);
        }
        else
            out.storeScalarGreek(result, packetNames[i], outNames[i]);
    }
}

///////////////////////////////////////////////////////
/** input for FastQuote */
CClassConstSP const FastQuoteInput::TYPE = CClass::registerClassLoadMethod(
    "FastQuoteInput", typeid(FastQuoteInput), load);

void FastQuoteInput::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(FastQuoteInput, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultFastQuoteInput);
    FIELD(optionType, "CALL , PUT or STRADDLE");
    FIELD(strike, "option strike");
    FIELD(maturity, "option maturity date");
    FIELD(volatility, "volatility");
    FIELD(marketPrice, "market price for implied vol");
    FIELD(resultFlag, "[I][P][D][V][R][B][W][E][S][T][M] for imp vol, price, delta&gamma, vega, rho, rho borrow, vega2sided&volgamma, ddeltadvol, vega skew, theta, mu"
                 " Calculation skipped for this option if no valid flag found.");
    FIELD(isAmerican, "true(default)=american, false=european, this is to replace the flag in FastQuoteEnv, use this one ! ");
    FIELD_MAKE_OPTIONAL(isAmerican);
    Addin::registerConstructor("FAST_QUOTE_INPUT",
                               Addin::RISK,
                               "create input for FastQuote",
                               TYPE);
}

string FastQuoteInput::needCalc()
{
    string flagUsed;

    // skip unsupported flags
    //if(resultFlag.find("I") != string::npos)
    // flagUsed += "I"; // implied vol is done first separately to allow its use for price etc.
    if(resultFlag.find("P") != string::npos)
        flagUsed += "P";
    if(resultFlag.find("D") != string::npos)
        flagUsed += "D";
    if (resultFlag.find("V") != string::npos)
        flagUsed += "V";
    if (resultFlag.find("R") != string::npos)
        flagUsed += "R";
    if (resultFlag.find("B") != string::npos)
        flagUsed += "B";
    if (resultFlag.find("W") != string::npos)
        flagUsed += "W"; // wisoo (actually Vega2Sided)
    if (resultFlag.find("E") != string::npos)
        flagUsed += "E"; // ddeltadvol
    if (resultFlag.find("S") != string::npos)
        flagUsed += "S"; // VegaSkew

    // always add indicative vol
    flagUsed += "K";  

    // if only implied vol is requested, don't return indicative vol 
    if (flagUsed.size() == 1 && resultFlag.find("I") != string::npos)
    {
        return string("");
    }
    else
    {
        return flagUsed;
    }
}

void FastQuoteInput::validatePop2Object()
{
    static const string method = "FastQuoteInput::validatePop2Object";
/*    if (resultFlag.find("T") != string::npos)
    {
        throw ModelException(method, "Theta not yet supported");
    } */
}

/////////////////////// FastQuoteEnv /////////////////////////

/** validate and set up */
void FastQuoteEnv::validatePop2Object()
{
    static const string method = "FastQuoteEnv::validatePop2Object";
    try {
        if (!market){
            throw ModelException(method, "market not provided");
        }
        if (!model){
            throw ModelException(method, "model not provided");
        }
        if (!baseInst){
            throw ModelException(method, "baseInstument not provided");
        }
    } 
    catch (exception& e){
        throw ModelException(e, method);
    }
}

// initialise pricing environment
void FastQuoteEnv::init(double spot)
{
    static const string method = "FastQuoteEnv::init";

    if (!needInit)
    {// already initialised
        setSpot(baseInstCopy->asset.get(), spot);
        return;
    }

    // initialise pricing env
    try
    {
        baseInstCopy = CVanillaSP::dynamicCast((IObjectSP)baseInst->clone());
        model->getInstrumentAndModelMarket(market.get(), baseInstCopy.get());
        baseInstCopy->Validate();

        // some validation
        if (baseInstCopy->fwdStarting){
            throw ModelException(method, "forward starting not supported.");
        }
        if (baseInstCopy->ccyTreatment !="N"
            && baseInstCopy->ccyTreatment !="V"){
            throw ModelException(method, "quanto or ccy struck not supported.");
        }

        if (spotRef <=0) // uninitialised
            spotRef = baseInstCopy->asset->getSpot();
        else
        {
            // override spot in baseInstCopy
            SpotLevel shift(spotRef);
            dynamic_cast<EquityBase*>(baseInstCopy->asset.get())->getEquity()->sensShift(&shift);
        }

        // set smoothing to default
        CTree1f* tree = dynamic_cast<CTree1f*> (model.get());
        if (tree !=0)
            tree->SetSmoothString("DEFAULT");

        // set last maturity using settlement info
        DateTime instEnd  = baseInstCopy->instSettle->settles(lastMaturity, baseInstCopy->asset.get());
        DateTime assetEnd = baseInstCopy->asset->settleDate(lastMaturity);
        lastMaturity = instEnd.isGreater(assetEnd) ? instEnd : assetEnd;
       
        if (lastMaturity < baseInstCopy->endDate(0))
            throw ModelException(method, "last maturity date is before base instrument maturity settlement date.");

        EquityBase* asset = dynamic_cast<EquityBase*>(baseInstCopy->asset.get());

        // get time metric for creating a vol
        ATMVolRequest atmReq;
        TimeMetricConstSP metric = CVolProcessedSP(asset->getVol()->getProcessedVol(&atmReq,asset))
            ->GetTimeMetric();
        // create vol
        if(isFlatVol)
        {  // flat vol
            FlatVol  vol(asset->getVol()->getName(),
                         baseInstCopy->getValueDate(),
                         metric.get(),
                         1.0); // 100% will not be used - otherwise there is a bug
            // create an equity cache
            assetCache = EquityCacheSP(new EquityCache(asset->getEquity().get(),
                                                       &vol));
        }
        else
        {// just use input vol 
            // create an equity cache
            assetCache = EquityCacheSP(new EquityCache(asset->getEquity().get(),
                                                       asset->getVol().get()));
        }
        // pre-calc forwrads
        DateTime divEnd = lastDivDate(dynamic_cast<EquityBase*>(baseInstCopy->asset.get()));
        assetCache->preCalcFwd(asset, baseInstCopy->getValueDate(), lastMaturity > divEnd ? lastMaturity : divEnd);
        // replace asset in instrument so that equityCache is used
        baseInstCopy->asset = CAssetWrapper(assetCache);

        // do this again after baseInstCopy has been modified
        model->getMarket(market.get(), IInstrumentCollection::singleton( baseInstCopy.get() ));

        needInit = false;
        // just set the spot 
        init(spot); 
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

//** set input to the required array element */
void FastQuoteEnv::setInput(const FastQuoteInput& input)
{
    static const string method = "FastQuoteEnv::setInput";

    // call or put
    if (input.optionType == "CALL")
        baseInstCopy->isCall = true;
    else if (input.optionType == "PUT")
        baseInstCopy->isCall = false;
    else
        throw ModelException(method, "option type not supported");

    // exercise schedule
    DateTimeArray mat(1, input.maturity);
    DoubleArray   strike(1, input.strike);;
    ScheduleSP exercise(new Schedule(mat, 
                                     strike, 
                                     input.isAmerican ? Schedule::INTERP_LINEAR : Schedule::INTERP_NONE));

    baseInstCopy->exerciseSchedule = exercise;
    baseInstCopy->canExerciseEarly = input.isAmerican;
    // vol
    if (isFlatVol)
        setVol(baseInstCopy->asset.get(), input.volatility);

}

/** call from high up, examine performance  */
void FastQuoteEnv::Price(FastQuoteInput& in, CResults& out)
{
    static const string method = "FastQuoteEnv::Price";
    bool iv_requested = false;
    bool theta_requested = false;
    bool mu_requested = false;
    bool other_requested = false;
    double iv = -1.0;

    double spot = baseInstCopy->asset->getSpot();

    // need to reset fwd arrays if cache used to get correct fwd at mat
    EquityCache* eq = dynamic_cast<EquityCache*>(baseInstCopy->asset.get());
    if (eq)
        eq->chooseFwdCache(EquityCache::BASE);

    string key = in.maturity.toString() + Format::toString(in.strike)
        + in.optionType + Format::toString(in.isAmerican);

    // calc fwd first here directly if requested, this avoids rho tweaked fwdArray being used for this
    if(in.resultFlag.find("F") != string::npos)
    {
        out.storeScalarGreek(baseInstCopy->asset->fwdValue(in.maturity),
                             CResults::FWD_AT_MAT_PACKET,
                             OutputNameSP(new OutputName(baseInstCopy->asset->getName())));
    }

    // choose tghe model to use
    IModelSP modelToUse = (in.isAmerican? model : closedForm);

    // implied vol is done before pricing etc. to allow vol override using implied
    if(in.resultFlag.find("I") != string::npos)
    {
        // set input override values
        setInput(in);
        iv_requested = true;
        bool success = false;
        if (useFastImplied)
        {
            success = impVolCache[key].FastImplied(spot, in.marketPrice, in.volatility,
                                                   baseInstCopy.get(), modelToUse.get(), 
                                                   baseInstCopy->asset.get(), out);
        }
        if (!success)
        {
            // control for implied vol
            SensitivityArraySP   sens(new SensitivityArray(0));
            OutputRequestArraySP dummyReq(new OutputRequestArray(0));
            sens->push_back(SensitivitySP(new ImpliedVol(in.marketPrice, in.volatility, 0.0001,"VALUE")));  
            Control ctrl(sens, dummyReq, false, "");

            // calc implied vol
            ctrl.calculate(modelToUse.get(), baseInstCopy.get(), &out);
        }

        OutputNameSP valueOutName(new OutputName("VALUE"));
        // if implied vol has an error (due to price < intrinsic) set value to -1
        // label is hard coded
        if (!out.isValidScalarGreek("IMPLIED_VOL", valueOutName))
        {
            out.storeScalarGreek(-1.0, "IMPLIED_VOL", valueOutName);
        }
        else
        {
            iv = out.retrieveScalarGreek("IMPLIED_VOL", valueOutName);
            // cache values
            impVolCache[key] = FastImpVol(iv, in.marketPrice, -1.0, spot);
        }
    }

    // do other requests
    string flag = in.needCalc();
    if (flag.size() > 0)
        other_requested = true;
    if (in.resultFlag.find("T") != string::npos)
        theta_requested = true;
    if (in.resultFlag.find("M") != string::npos)
        mu_requested = true;

    
    if (other_requested || theta_requested || mu_requested) {
        // set input override values if not done yet
        if(!iv_requested) // already done otherwise
            setInput(in);
        else
        {// reset vol
            if (IVpricing && iv > 0.0)
                setVol(baseInstCopy->asset.get(), iv);
            else
                setVol(baseInstCopy->asset.get(), in.volatility);
        }
    }
    
    double volForInterp = -1.0;

    if (other_requested)
    {
        // create control according to flag
        // "[I][P][D][V][R][W][E][S]" for  I=imp vol, P=price, D=delat/gamma, V=vega, R=rho, W = vega2sided/volgamma, 
        // E = DDeltaDVol, S = VegaSkew
        CControlSP ctrl(Control::makeFromFlags(flag, in.marketPrice));
        // use interp if possible
        if(!iv_requested && useInterp) 
        {
            string keyOther = flag + key; // can only deal if the output requested are the same

            volForInterp = in.volatility;
            if (!isFlatVol)
            {
                // get vol from asset
                LinearStrikeVolRequest volRequest(in.strike, 
                                                  baseInstCopy->getValueDate(), 
                                                  in.maturity,
                                                  false);
                // interpolate the vol using our LN request
                CVolProcessedBSSP volBS(baseInstCopy->asset->getProcessedVol(&volRequest));
                volForInterp = volBS->CalcVol(baseInstCopy->getValueDate(), in.maturity);
            }

            if (resultCache.find(keyOther) != resultCache.end()
                && resultCache[keyOther].canInterpolate(flag, spot, volForInterp))
            {// inpterp from spline
                resultCache[keyOther].interpolate(spot, out);
            }
            else
            {// calculate price and store spline coeffs
                computeSpline(modelToUse, ctrl, spot, volForInterp, flag, out, resultCache[keyOther]);
            }
        }
        else
        {   // call pricing from control level
            ctrl->calculate(modelToUse.get(), baseInstCopy.get(), &out);
        }
    }

    // do theta if needed. must be done separatelly because it messes up pointers to forwardCache
    if(theta_requested) {
        // create control according to flag
        flag = "T";
        CControlSP ctrl(Control::makeFromFlags(flag, in.marketPrice));
        // cannot use interp
        ctrl->calculate(modelToUse.get(), baseInstCopy.get(), &out);

        // move forward cache array back
        eq->chooseFwdCache(EquityCache::BASE);
    }

    // do mu if needed. must be done separatelly because it messes up pointers to forwardCache
    if(mu_requested) {
        // create control according to flag
        flag = "M";
        CControlSP ctrl(Control::makeFromFlags(flag, in.marketPrice));
        // cannot use interp
        ctrl->calculate(modelToUse.get(), baseInstCopy.get(), &out);
        
        // move forward cache array back
        eq->chooseFwdCache(EquityCache::BASE);
    }
}

/** compute spline and result */
void FastQuoteEnv::computeSpline(IModelSP& model, CControlSP& ctrl, double spot, double vol, const string& flag, 
                                 CResults& out, FastQuoteSpline& store)
{
    const int NUM_PTS = 11;
    const int MIDDLE_PT = 5;
    const double grid[NUM_PTS] = {0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.03, 1.04, 1.05};

    // calc initial result
    CResultsArray results(NUM_PTS);
    for (int i=0; i<NUM_PTS; i++)
        results[i] = CResultsSP(new CResults());

    // use middle point
    results[MIDDLE_PT] = CResultsSP::attachToRef(&out); 

    // set spot grids
    store.spotGrid.resize(NUM_PTS);
    for (int j=0; j<NUM_PTS; j++)
    {
        store.spotGrid[j] = grid[j]*spot;
        setSpot(baseInstCopy->asset.get(), grid[j]*spot);                  
        ctrl->calculate(model.get(), baseInstCopy.get(), results[j].get());         
    }

    //restore spot
    setSpot(baseInstCopy->asset.get(), spot);

    store.prepareSpline(flag, vol, results);
}

/** Invoked when Class is 'loaded' */
void FastQuoteEnv::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(FastQuoteEnv, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultFastQuoteEnv);
    FIELD(market, "market");
    FIELD(model, "model");
    FIELD(baseInst, "base instrument for getting settlment, growth curve data.");
    FIELD(lastMaturity, "longest maturity for pricing environment");
    FIELD(isAmerican, "no longer used, please use the flag in FastQuoteInput");
    FIELD_MAKE_OPTIONAL(isAmerican);
    FIELD(isFlatVol, "true=flat vol treatment, flase=vol term structure(offset) treatment");
    FIELD_MAKE_OPTIONAL(isFlatVol);
    FIELD(spotRef, "spot level the pricing environment is set up for - use current spot");
    FIELD_MAKE_OPTIONAL(spotRef);
    FIELD(IVpricing, "false=use input vol for price and sensitivities, true(default)=use implied vol if requested, input vol if IV not requested");
    FIELD_MAKE_OPTIONAL(IVpricing);
    FIELD(useInterp, "false(default)=do not use interpolation, true = calc grid and use interpolation when possible. Interpolation is not used if implied vol requested.");
    FIELD_MAKE_OPTIONAL(useInterp);
    FIELD(useFastImplied, "false(default)=use standard implied method, true = use fast implied method.");
    FIELD_MAKE_OPTIONAL(useFastImplied);
    Addin::registerConstructor("FAST_QUOTE_ENV",
                               Addin::RISK,
                               "create a pricing environment for real time quote apps",
                               TYPE);
}

CClassConstSP const FastQuoteEnv::TYPE = CClass::registerClassLoadMethod(
    "FastQuoteEnv", typeid(FastQuoteEnv), load);

// array has to have its own type, can we get rid of this ?
DEFINE_TEMPLATE_TYPE(FastQuoteInputArray);

DRLIB_END_NAMESPACE

