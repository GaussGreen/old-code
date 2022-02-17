//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VanillaGrid.cpp
//
//   Description : Grid of European Vanillas
//
//   Author      : Regis Guichard
//
//   Date        : 27 May 02
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
//#include "edginc/UtilFuncs.hpp"
#include "edginc/VanillaGrid.hpp"
//#include "edginc/mathlib.hpp"
//#include "edginc/Maths.hpp"
#include "edginc/Black.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
//#include "edginc/LinearStrikeSpreadVolRequest.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/InstrumentUtil.hpp"
//#include "edginc/XCB.hpp"
//#include "edginc/ShiftSizeCollector.hpp"
//#include "edginc/Events.hpp"
//#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/CashSettlePeriod.hpp"
#include "edginc/Enum2StringListHelper.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/VolRequestTime.hpp"
#include "edginc/AtMaturity.hpp"
#include "edginc/DeltaToStrike.hpp"
#include "edginc/Delta.hpp"
#include "edginc/Scenario.hpp"
#include "edginc/BenchmarkDate.hpp"
//#include "edginc/YieldCurve.hpp"

#include "edginc/IndexSpecEQ.hpp"
#include "edginc/LatticeProdEDR.hpp"
#include "edginc/TreeSliceOperExt.hpp"

#include "edginc/VolSpline.hpp"

DRLIB_BEGIN_NAMESPACE

// VANILLA GRID OUTPUT

const int VanillaGrid::Output::ErrorCode::IMP_VOL_ERROR = -6;           // excel #NUM! error
const int VanillaGrid::Output::ErrorCode::IMP_VOL_NOT_REQUESTED = -7;   // excel #N/A error
const int VanillaGrid::Output::ErrorCode::PRICE_NOT_REQUESTED = -7;     // excel #N/A error

/** Invoked when Class is 'loaded' */
void VanillaGrid::Output::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(VanillaGrid::Output, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultCtor);
    FIELD(values, "values");
    FIELD(errors, "errors");
}

IObject* VanillaGrid::Output::defaultCtor(){
    return new Output();
}

VanillaGrid::Output::Output():
CObject(TYPE){}

VanillaGrid::Output::Output(int numCols, int numRows):
CObject(TYPE),
values(numCols, numRows),
errors(numCols){
    int iCol = 0;
    for(; iCol < numCols; ++iCol){
        errors[iCol] = IntArraySP(new IntArray(numRows));
    }
}

int VanillaGrid::Output::numRows() const{
    return values.numRows();
}

int VanillaGrid::Output::numCols() const{
    return values.numCols();
}

void VanillaGrid::Output::setValue(int iCol, int iRow, double value){
    values[iCol][iRow] = value;
    (*errors[iCol])[iRow] = 0;
}

double VanillaGrid::Output::getValue(int iCol, int iRow) const{
    if (isError(iCol, iRow)){
        throw ModelException("VanillaGrid::Output::getValue",
                             "Cannot access value for elt (" +
                             Format::toString(iCol) +
                             ", " +
                             Format::toString(iRow) +
                             "), as this elt is a (" +
                             Format::toString((*errors[iCol])[iRow]) +
                             " code) error");
    }
    return values[iCol][iRow];
}

void VanillaGrid::Output::setError(int iCol, int iRow, int code){
    (*errors[iCol])[iRow] = code;
}

int VanillaGrid::Output::getError(int iCol, int iRow) const{
    if (!isError(iCol, iRow)){
        throw ModelException("VanillaGrid::Output::getValue",
                             "Cannot access error type for elt (" +
                             Format::toString(iCol) +
                             ", " +
                             Format::toString(iRow) +
                             "), as this elt is not an error");
    }
    return (*errors[iCol])[iRow];
}

bool VanillaGrid::Output::isError(int iCol, int iRow) const{
    return !!(*errors[iCol])[iRow];
}

CClassConstSP const VanillaGrid::Output::TYPE = CClass::registerClassLoadMethod(
    "VanillaGrid::Output", typeid(VanillaGrid::Output), load);

/* VanillaGrid class */
void VanillaGrid::Validate() {
    static const string method = "VanillaGrid::Validate";
    // just check the things that aren't/cannot be checked in
    // validatePop2Object
    try{
        if (!asset){
            throw ModelException(method, "Asset is null");
        }
        if (!discount){
            throw ModelException(method, "Discount YC is null");
        }

        AssetUtil::assetCrossValidate(asset.get(),
                                      fwdStarting,
                                      startDate,
                                      valueDate,
                                      discount,
                                      this);

        // Removed this validation as FX should work for the next release
        //if (FXAsset::TYPE->isInstance(asset.get())){
        //    throw ModelException(method,
        //                         "Options on FX assets are not allowed yet");
        //}

        if (!valueDate.isLess((*maturities)[0])){
            throw ModelException(method, "Dead instrument not supported");
        }
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

void VanillaGrid::validatePop2Object() {
    static const string method("VanillaGrid::validatePop2Object");
    try{
        // can't get instrument settlement from Market - fail if it is NULL
        if(!instSettle ){
            throw ModelException(method, "Instrument settlement is NULL");
        }

        // if fwd starting, can't be one contract
        if (fwdStarting && oneContract){
            throw ModelException(method, "Can't be forward starting and one contract");
        }

        if (oneContract){
            /* one contract is pricing for a fixed notional where
               notional/initial spot = 1.0 */
            notional    = 1.0;
            initialSpot = 1.0;
        }

        // Make sure at least weightMatrix is passed or (strikes, weights, maturities)
        if (!weightMatrix) {
            if (!strikes || !maturities || !weights) {
                throw ModelException(method, "At least a weightMatrix or (strikes, maturities, weights) must be passed");
            }
        }

        if (!((!weightMatrix.getMO())
                ||(!asset.getMO())
                ||(!discount.getMO()))){

            // compute absolute strikes and types
            weightMatrix->computeStrikesAndTypes(asset,discount,instSettle);

            strikeUnits = weightMatrix->getStrikeUnits();
            instType = weightMatrix->getInstType();
            strikes = weightMatrix->getStrikes();
            maturities = weightMatrix->getDateTimeMaturities();
            weights = weightMatrix->getWeights();
            instStrikesUsed = weightMatrix->getInstStrikesUsed();
            instTypesUsed = weightMatrix->getInstTypesUsed();
        }
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

void VanillaGrid::GetMarket(const IModel*          model,
                            const CMarketDataSP    market){
    static const string method("VanillaGrid::GetMarket");
    try{
        market->GetReferenceDate(valueDate);
        CAsset::getAssetMarketData(model,
                                   market.get(),
                                   ccyTreatment,
                                   discount,
                                   asset);

        discount.getData(model, market);
        instSettle->getMarket(model, market.get());
        if (premiumSettle.get()){
            premiumSettle->getMarket(model, market.get());
        }

        if (weightMatrix.isEmpty()){
            // if no weightmatrix is available, create one
            string instType;
            instType = isCall ? WeightMatrix::CALL : WeightMatrix::PUT;
            WeightMatrixSP wM(new WeightMatrix("",
                                                   maturities,
                                                   strikes,
                                                   weights,
                                                   WeightMatrix::ABSOLUTE,
                                                   instType,
                                                   valueDate));
            weightMatrix = WeightMatrixWrapper(wM);
        }

        weightMatrix.getData(model, market.get());

        //compute absolute strikes and types
        weightMatrix->computeStrikesAndTypes(asset,discount,instSettle);

        //copy fields
        strikeUnits = weightMatrix->getStrikeUnits();
        instType = weightMatrix->getInstType();
        strikes = weightMatrix->getStrikes();
        maturities = weightMatrix->getDateTimeMaturities();
        weights = weightMatrix->getWeights();
        instStrikesUsed = weightMatrix->getInstStrikesUsed();
        instTypesUsed = weightMatrix->getInstTypesUsed();

        //Check that strikeUnits is absolute if forward starting
        if (fwdStarting&&!CString::equalsIgnoreCase(strikeUnits,WeightMatrix::ABSOLUTE)
            &&!CString::equalsIgnoreCase(strikeUnits,WeightMatrix::DEFAULT)) {
            throw ModelException(method,"strikeUnits " + strikeUnits + " is not supported when forward starting");
        }
    }
    catch(exception& e){
        throw ModelException(e, "VanillaGrid::GetMarket");
    }
}

DateTime VanillaGrid::getValueDate() const{
  return valueDate;
}

/** Returns the name of the instrument's discount currency. */
string VanillaGrid::discountYieldCurveName() const {
    return discount.getName();
}

/** The weight matrix */
CDoubleMatrixConstSP VanillaGrid::getWeights() const {
    return weights;
}

/** The notional */
double VanillaGrid::getNotional() const {
    return notional;
}

/** when to stop tweaking */
DateTime VanillaGrid::endDate(const Sensitivity* sensControl) const {
    DateTime matDate = maturities->back();
    DateTime instEnd  = instSettle->settles(matDate, asset.get());
    DateTime assetEnd = asset->settleDate(matDate);
    DateTime end = instEnd.isGreater(assetEnd) ? instEnd : assetEnd;
    return end;
}

#if 0
void VanillaGrid::addOutputRequests(Control* control,
                                    Results* results,
                                    const double& fairValue,
                                    const double& indVol) const
{
    // take care of additional outputs
    if ( control && control->isPricing() )
    {
        DateTime matDate = exerciseSchedule->lastDate();
        // DELAY_PRICE
        InstrumentUtil::delayPriceHelper(control,
                                         results,
                                         fairValue,
                                         valueDate,
                                         discount.get(),
                                         asset.get(),
                                         premiumSettle.get());
        // IND_VOL
        OutputRequest* request = NULL;
        if (control->requestsOutput(OutputRequest::IND_VOL, request))
        {
             if (isCall)
             {
                 double indVolThreshold = 0.01;
                 double strike = exerciseSchedule->lastValue();
                 if (!fwdStarting)
                 {
                     double spot = asset->getSpot();
                     indVolThreshold *= spot;
                 }
                 // If strike is less than 1% of the spot store NotApplicable.
                 // This is required by the result combining code later
                 if (strike < (indVolThreshold))
                 {
                     results->storeNotApplicable(request);
                 }
                 else
                 {
                     results->storeRequestResult(request, indVol);
                 }
             }
             else
             {
                 results->storeRequestResult(request, indVol);
             }
        }

        // FWD_AT_MAT
        InstrumentUtil::recordFwdAtMat(control,
                                       results,
                                       matDate,
                                       valueDate,
                                       asset.get());
    }
}

/** price a dead instrument until settlement - exercised, expired, knocked out etc.
returns true if it is dead (and priced), false if it is not dead */
bool VanillaGrid::priceDeadInstrument(CControl* control,
                                      CResults* results) const{
    static string method = "VanillaGrid::priceDeadInstrument";

    return false;
}
#endif

/** private class */
class VanillaGridClosedForm: public CClosedFormLN::IProduct{
private:
    const VanillaGrid* inst; // a reference

public:
    VanillaGridClosedForm(const VanillaGrid* inst):
    inst(inst){}

    void price(CClosedFormLN*   model,
               Control*         control,
               CResults*        results) const{
        static const string method = "VanillaGridClosedForm::price";

        try {
            if (!control){
                return;
            }

            OutputRequest* indVolRequest = 0;
            control->requestsOutput(OutputRequest::IND_VOL, indVolRequest);
            OutputRequest* optionPriceRequest = 0;
            control->requestsOutput(OutputRequest::OPTION_PRICE, optionPriceRequest);
            OutputRequest* optionVegaRequest = 0;
            control->requestsOutput(OutputRequest::OPTION_VEGA, optionVegaRequest);
            bool doSThing = indVolRequest || optionPriceRequest || optionVegaRequest;

            if (!doSThing){
                return; // nothing to do
            }

            DateTime startDate = inst->fwdStarting ? inst->startDate: inst->valueDate;

            /* Calculate variances */
            const DateTimeArray& maturities = *inst->maturities;
            const DoubleArray& strikes = *inst->strikes;
            int nbStrikes = strikes.size();
            int nbMats = maturities.size();
            const DoubleMatrix& instStrikesUsed = *inst->instStrikesUsed;
            const IntArrayArray& instTypesUsed = *inst->instTypesUsed;
            vector<DoubleArraySP> variances(nbStrikes);
            DateTime lastMatDate = maturities.back();
            CVolProcessedBSSP volBS;
            int iStrike, iMat;

            VolRequestLNStrikeSP volRequest(new LinearStrikeTSVolRequest(1.0,
                                                                         startDate,
                                                                         lastMatDate, // used ?
                                                                         inst->fwdStarting));
            volRequest->allowNegativeFwdVar(true);

            for (iStrike = 0; iStrike < nbStrikes; iStrike++) {
                variances[iStrike] = DoubleArraySP(new DoubleArray(nbMats));
                for (iMat = 0; iMat < nbMats; iMat++) {
                    volRequest->setStrike(instStrikesUsed[iMat][iStrike]);
                    volBS = CVolProcessedBSSP(inst->asset->getProcessedVol(volRequest.get()));
                    (*(variances[iStrike]))[iMat] = volBS->CalcVar(startDate, maturities[iMat]);
                }
            }

            /* Calculate vols and time fractions and export what is requested */
            DoubleArraySP timeFracs;
            VanillaGrid::OutputSP vols;
            if (indVolRequest || optionVegaRequest){
                timeFracs = CDoubleArraySP(new DoubleArray(nbMats));
                vols = VanillaGrid::OutputSP(new VanillaGrid::Output(nbMats, nbStrikes));
                iMat = 0;
                for (; iMat < nbMats; ++iMat){
                    (*timeFracs)[iMat] = volBS->calcTradingTime(startDate,
                                                                maturities[iMat]);
                    for (iStrike = 0; iStrike < nbStrikes; ++iStrike){
                        vols->setValue(iMat, iStrike,
                                       sqrt((*variances[iStrike])[iMat] / (*timeFracs)[iMat]));
                    }
                }
                if (indVolRequest){
                    results->storeRequestResult(indVolRequest, vols);
                }
            }

            /* If neither option price nor option vega is requested, we're done */
            if (!optionPriceRequest && !optionVegaRequest){
                return; // nothing left to do
            }

            /* If it's forward starting, convert percentage strikes
               to absolute values based on spot at start date */
            double fwdAtStart;
            CDoubleMatrixSP absStrikesUsedSP(copy(&instStrikesUsed));
            if (inst->fwdStarting){
                fwdAtStart = inst->asset->fwdValue(inst->startDate);
                absStrikesUsedSP->scale(fwdAtStart);
            }

            /* Calculate fwds */
            DoubleArray fwdPrices(nbMats);
            inst->asset->fwdValue(maturities,
                                  fwdPrices);

            double scalingFactor = InstrumentUtil::scalePremium(inst->oneContract,
                                                                inst->fwdStarting,
                                                                inst->notional,
                                                                fwdAtStart,
                                                                inst->initialSpot);

            /* Calculate and export premiums and vegas */
            VanillaGrid::OutputSP premiums;
            VanillaGrid::OutputSP vegas;
            if (optionPriceRequest){
                premiums = VanillaGrid::OutputSP(new VanillaGrid::Output(nbMats, nbStrikes));
            }
            if (optionVegaRequest){
                vegas = VanillaGrid::OutputSP(new VanillaGrid::Output(nbMats, nbStrikes));
            }

            iMat = 0;
            for (; iMat < nbMats; ++iMat){
                double discFactor = inst->instSettle->pv(inst->valueDate,
                                                         maturities[iMat],
                                                         inst->discount.get(),
                                                         inst->asset.get());

                for (iStrike = 0; iStrike < nbStrikes; ++iStrike){
                    bool isCall = (*instTypesUsed[iMat])[iStrike] ? true : false;
                    if (optionPriceRequest){
                        premiums->setValue(iMat, iStrike,
                                           scalingFactor * Black::price(isCall,
                                                                        fwdPrices[iMat],
                                                                        (*absStrikesUsedSP)[iMat][iStrike],
                                                                        discFactor,
                                                                        (*variances[iStrike])[iMat]));
                    }

                    if (optionVegaRequest){
                        vegas->setValue(iMat, iStrike,
                                        scalingFactor * Black::vega(isCall,
                                                                    fwdPrices[iMat],
                                                                    (*absStrikesUsedSP)[iMat][iStrike],
                                                                    discFactor,
                                                                    (*timeFracs)[iMat],
                                                                    vols->getValue(iMat, iStrike)));
                    }
                }
            }

            if (optionPriceRequest){
                results->storeRequestResult(optionPriceRequest, premiums);
            }
            if (optionVegaRequest){
                results->storeRequestResult(optionVegaRequest, vegas);
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }
};

/// ------------  start of FD product -----------
/////////////////////////////////////////////////////////
//           private class for FD/tree product - new state variable interface
/////////////////////////////////////////////////////////
/** vanilla product payoff for a FD */
class VanillaGridProd: public LatticeProdEDR
{
public:
    VanillaGridProd(const VanillaGrid* vanilla, FDModel* m) :
        LatticeProdEDR(m),
        instrVan(vanilla)
    {
        // first: set discount curve
        if( tree1f )
            tree1f->setDiscountCurve( instrVan->discount.getSP() );

        // second: create spot payoff
        payoffIndex = model->createProduct( IProdCreatorSP( new
            IndexSpecEQ( instrVan->asset.getName(), instrVan->asset, instrVan->ccyTreatment ) ) );

        try{
            if (!CString::equalsIgnoreCase(vanilla->strikeUnits,WeightMatrix::ABSOLUTE)){
                throw ModelException("VanillaGridProd::VanillaGridProd", "strikeUnits:"
                    + vanilla->strikeUnits + "is not supported by VanillaGridProd");
            }
            if (!CString::equalsIgnoreCase(vanilla->instType,WeightMatrix::CALL)
                &&!CString::equalsIgnoreCase(vanilla->instType,WeightMatrix::PUT)){
                throw ModelException("VanillaGridProd::VanillaGridProd", "instType:"
                    + vanilla->instType + "is not supported by VanillaGridProd");
            }
        }
        catch (exception& e) {
            throw ModelException(e, "VanillaGridProd::VanillaGridProd");
        }

    }


    /** isInitValue == true, payoff at T for backward or value at t=0 for fwd induction
        isInitValue == false, payoff boundary condition, for KO, early exercise etc. */
    virtual void update(int& step, FDProduct::UpdateType type)
    {
        // we assume just need the first und level for spot here for now
        const TreeSlice & s = payoffIndex->getValue( step );

        const vector< TreeSliceSP > & price = slices;
        int pStart = 0, pEnd = price.size() - 1;

        // backward induction, at maturity date t = T
        if (type == FDProduct::BWD_T)
        {
            prod_BWD_T( s,
                        step,
                        pStart,
                        pEnd,
                        price);
        }
        // backward induction, at date t < T
        else if (type == FDProduct::BWD)
        {
            prod_BWD(   s,
                        step,
                        pStart,
                        pEnd,
                        price);
        }
        // forward induction, at initial date t = 0
        else if (type == FDProduct::FWD_0)
        {
            prod_FWD_0( s,
                        step,
                        pStart,
                        pEnd,
                        price);
        }
        // forward induction, at date t > 0
        else
        {
            prod_FWD(  s,
                        step,
                        pStart,
                        pEnd,
                        price);
        }
    };

    /** initialise model1f - allow product customisation */
    virtual void init(CControl* control) const;

    /** initialising and setting product variables */
    // this is called per pricing call before tree sweep call
    virtual void initProd();

    /** override output method */
    virtual void recordOutput(Control* control, YieldCurveConstSP disc, Results* results);

    /** product payoff method at maturity date for backward induction */
    virtual void prod_BWD_T(
        const TreeSlice & s,
        int step,
        int pStart,
        int pEnd,
        const vector< TreeSliceSP > & price);

    // product payoff method before maturity date for backward induction
    virtual void prod_BWD(
        const TreeSlice & s,
        int step,
        int pStart,
        int pEnd,
        const vector< TreeSliceSP > & price);

    virtual DateTime getStartDate() const{
        return instrVan->fwdStarting ? instrVan->startDate : instrVan->valueDate;
    }

    /** product payoff method at initial date for forward induction */
    virtual void prod_FWD_0(
        const TreeSlice & s,
        int step,
        int pStart,
        int pEnd,
        const vector< TreeSliceSP > & price);

    // product payoff method after initial date for forward induction
    virtual void prod_FWD(
        const TreeSlice & s,
        int step,
        int pStart,
        int pEnd,
        const vector< TreeSliceSP > & price);

    // to be removed when new mdf is done
    virtual DateTime getValueDate(){ return instrVan->getValueDate();}

private:
    const VanillaGrid* instrVan;
    double scalingFactor; // scaling factor
    // stepIsMat indicates if the current step is one VanillaGrid maturities
    vector<bool> stepIsMat;
    // stepMatIndex indicates the index of the maturity
    // for steps corresponding to VanillaGrid maturities
    vector<int> stepMatIndex; // index of input maturities
    DoubleArrayArray pricesFwdBwd; // pricesFwdInduction stores the prices for forward induction
};

/** initialise tree1f - allow product customisation
    must not initialise product variables here - use initProd() */
void VanillaGridProd::init(CControl* control) const{

    static const string method = "VanillaGridProd::Init";
    try {
        // add maturity dates to critical dates
        model->addCritDates( *instrVan->maturities );

        /** customize tree parameters here and set up the tree */
        DateTimeArray segDates;
        segDates.resize(2);

        // start date if forward starting and value date otherwise
        if (instrVan->fwdStarting && instrVan->startDate>instrVan->valueDate)
        {
            segDates[0] = instrVan->startDate;
        }
        else
        {
            segDates[0] = instrVan->valueDate;
        }

        // longest maturity in VanillaGrid instrument
        segDates[1] = (*instrVan->maturities)[(*instrVan->maturities).size()-1];

        IntArray density(1,1);

        // prepare timeline set up
        model->initSegments( segDates, density );
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** initialising and setting product variables
    this is called per pricing call before tree sweep call (after InitTree)  */
void VanillaGridProd::initProd() {

    const DateTimeArray& maturities = *instrVan->maturities;
    int nbMats = maturities.size();

    const DoubleArray& strikes = *instrVan->strikes;
    int nbStrks = strikes.size();

    // init slice
    initSlices(getModel()->IsFwdInduction() ? 1 : nbStrks * nbMats);

    // pricesFwdInduction stores the prices for forward induction
    pricesFwdBwd.resize(nbMats);
    int iMat;
    for (iMat = 0 ; iMat < nbMats ; iMat++)
    {
        pricesFwdBwd[iMat].resize(nbStrks);
    }

    // stepIsMat indicates if the current step is one VanillaGrid maturities
    stepIsMat.resize(getModel()->getLastStep() + 1);
    // stepMatIndex indicates the index of the maturity
    // for steps corresonding to VanillaGrid maturities
    stepMatIndex.resize(getModel()->getLastStep() + 1);

    int it;
    for (it = 0 ; it <= getModel()->getLastStep() ; it++)
    {
        stepIsMat[it] = false; // initialize
        stepMatIndex[it] = 0; // initialize
    }

    it = 0, iMat = 0;
    while (iMat < nbMats)
    {
        while ((getModel()->getDate(it) < maturities[iMat])
               && (it <= getModel()->getLastStep()))
        {
            it++;
        }
        stepIsMat[it] = true;
        stepMatIndex[it] = iMat;
        iMat++;
    }

    // scaling factor
    if (instrVan->fwdStarting)
    {
        // forward starting => notional based, and % strike
        scalingFactor = instrVan->notional;
    }
    else
    {
        scalingFactor = (instrVan->oneContract) ? 1.0 : (instrVan->notional / instrVan->initialSpot);
    }
}

/* output results as a matrix */
void VanillaGridProd::recordOutput(Control* control, YieldCurveConstSP disc, Results* results)
{
    static const string method = "VanillaGridProd::recordOutput";
    try {
        OutputRequest* indVolRequest = 0;
        control->requestsOutput(OutputRequest::IND_VOL, indVolRequest);

        OutputRequest* optionPriceRequest = 0;
        control->requestsOutput(OutputRequest::OPTION_PRICE, optionPriceRequest);

        bool doSThing = indVolRequest || optionPriceRequest;

        if (!doSThing)
        {
            return; // nothing to do
        }

        // strikes and maturities
        const DateTimeArray& maturities = *instrVan->maturities;
        int nbMats = maturities.size(); // nber of maturities
        const DoubleArray& strikes = *instrVan->strikes;
        int nbStrks = strikes.size(); // nber of strikes
        const DoubleMatrix& instStrikesUsed = *instrVan->instStrikesUsed;

        // forward prices
        DoubleArray fwds(nbMats);
        instrVan->asset->fwdValue(maturities, fwds);

        VanillaGrid::OutputSP prices;
        VanillaGrid::OutputSP indVols;
        if (optionPriceRequest)
        {
            prices = VanillaGrid::OutputSP(new VanillaGrid::Output(nbMats, nbStrks));
        }
        if (indVolRequest)
        {
            indVols = VanillaGrid::OutputSP(new VanillaGrid::Output(nbMats, nbStrks));
        }

        // time metric (time fractions from value date to maturities)
        TimeMetricConstSP timeMetric = getModel()->getTimeMetric();
        // time fractions
        DoubleArray timeFracs = DoubleArray(nbMats);
        int iMat;
        for(iMat = 0 ; iMat < nbMats ; iMat++ )
        {
            timeFracs[iMat] = timeMetric->yearFrac(instrVan->valueDate, maturities[iMat]);
        }

        // pv factors
        DoubleArray pvFactors = DoubleArray(nbMats);
        for(iMat = 0 ; iMat < nbMats ; iMat++ )
        {
            pvFactors[iMat] = instrVan->instSettle->pv(instrVan->valueDate,
                                                       maturities[iMat],
                                                       instrVan->discount.get(),
                                                       instrVan->asset.get()); // discount factor
        }

        // if backward induction, need to set pricesFwdBwd
        if (!getModel()->IsFwdInduction())
        {
            // get prices at t = 0
            // save a value to top level instrument price
            results->storePrice(model->getPrice0( *slices[0] ), disc->getCcy());

            int iMat;
            for (iMat = 0 ; iMat < nbMats ; ++iMat)
            {
                int iStrk;
                for (iStrk = 0 ; iStrk < nbStrks ; ++iStrk)
                {
                    int iPrice = iStrk + (nbMats - 1 - iMat) * nbStrks;
                    pricesFwdBwd[iMat][iStrk] = model->getPrice0( *slices[iPrice] );
                }
            }
        }
        // if forward induction, need to discount pricesFwdBwd
        if (getModel()->IsFwdInduction())
        {
            int iMat;
            for (iMat = 0 ; iMat < nbMats ; ++iMat)
            {
                int iStrk;
                for (iStrk = 0 ; iStrk < nbStrks ; ++iStrk)
                {
                    pricesFwdBwd[iMat][iStrk] *= pvFactors[iMat];
                }
            }
        }

        // loop on maturities
        for (iMat = 0 ; iMat < nbMats ; ++iMat)
        {
            // weights for the this maturity
            const double* weights = (*instrVan->weights)[iMat];

            // loop on strikes
            int iStrk;
            for (iStrk = 0 ; iStrk < nbStrks ; ++iStrk)
            {
                if (!weights || Maths::isPositive(weights[iStrk]))
                {
                    double price = pricesFwdBwd[iMat][iStrk];
                    double variance = 0.0;

                    if (indVolRequest){

                        bool isCall = (*((*instrVan->instTypesUsed)[iMat]))[iStrk] ? true : false;
                        // call version of impliedVariance that returns status
                        if (Black::impliedVariance(isCall, // option call or put
                                                   fwds[iMat], // forward price
                                                   instStrikesUsed[iMat][iStrk], // strike
                                                   pvFactors[iMat], // pv
                                                   0.3 * 0.3 * timeFracs[iMat], // initial var guess
                                                   price, // option price
                                                   2.0 * 0.3 * 1.0e-5 * timeFracs[iMat], // var accuracy
                                                   variance)) // variance
                        {
                            indVols->setValue(iMat, iStrk, sqrt(variance / timeFracs[iMat]));
                        }
                        else
                        {
                            indVols->setError(iMat, iStrk, VanillaGrid::Output::ErrorCode::IMP_VOL_ERROR);
                        }
                    }
                    if (optionPriceRequest)
                    {
                        prices->setValue(iMat, iStrk, scalingFactor * price);
                    }
                }
                else
                {
                    if (indVolRequest)
                    {
                        indVols->setError(iMat, iStrk, VanillaGrid::Output::ErrorCode::IMP_VOL_NOT_REQUESTED);
                    }
                    if (optionPriceRequest)
                    {
                        prices->setError(iMat, iStrk, VanillaGrid::Output::ErrorCode::PRICE_NOT_REQUESTED);
                    }
                }
            }
        }

        if (optionPriceRequest)
        {
            results->storeRequestResult(optionPriceRequest, prices);
        }
        if (indVolRequest)
        {
            results->storeRequestResult(indVolRequest, indVols);
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

/** product payoff method at maturity date for backward induction */
void VanillaGridProd::prod_BWD_T(
    const TreeSlice & s,
    int step,
    int pStart,
    int pEnd,
    const vector< TreeSliceSP > & price)
{
    static const string method = "VanillaGridProd::prod_BWD_T";

    double settlePV = instrVan->instSettle->pvAdjust((*instrVan->maturities.get())[instrVan->maturities->size()-1],
                                                     instrVan->discount.get(),
                                                     instrVan->asset.get());
    int nbStrikes = instrVan->strikes->size();
    int nbMats = instrVan->maturities->size();
    int nbOfP = nbStrikes * nbMats;

    // check nbOfP is consistent
    if (nbOfP != pEnd - pStart + 1)
    {
        throw ModelException(method,
                             "The number of prices is not consistent.");
    }

    // fill the first nK payoff at the longest maturity
    // and also set the correct pStart and pEnd;
    pStart = 0;
    pEnd = nbStrikes -1;

    int iStrk;
    for (iStrk = 0 ; iStrk < nbStrikes ; iStrk++)
    {
// Using template slice operators
        (*price[iStrk]) =
            smax( 0., ( s - (*instrVan->strikes)[iStrk] ) * ( instrVan->isCall ? settlePV : -settlePV ) );
    }

    // fill other price by zero
    int iT;
    for (iT = 1 ; iT < nbMats ; iT++)
    {
        int iStrk;
        for (iStrk = 0 ; iStrk < nbStrikes ; iStrk++)
        {
// Using template slice operators
            (*price[iStrk + iT * nbStrikes]) = 0.;
        }
    }
}

/** product payoff method at initial date for forward induction */
void VanillaGridProd::prod_FWD_0(
    const TreeSlice & s,
    int step,
    int pStart,
    int pEnd,
    const vector< TreeSliceSP > & price)
{
    static const string method = "VanillaGridProd::prod_FWD_0";

    ASSERT( dynamic_cast< const TreeSliceEQ * >( price[0].get() ) );

    //for test
    int nbOfP = 1;

    // check nbOfP is consistent
    if (nbOfP != pEnd - pStart + 1)
    {
        throw ModelException(method,
                             "The number of prices is not consistent.");
    }

    IntArray initialIndex; // array storing the index of initial conditions
    DoubleArray initialValue; // array storing the values of initial conditions
    getModel()->getInitialConditions(initialIndex, initialValue);

    // currently, the forward induction is only supported fo 2D models
    if (initialIndex.size() == 2) //2D
    {
// Using template slice operators
        (*price[0]) = 0.;
        static_cast< TreeSliceEQ & >( (*price[0]) )( initialIndex[0], initialIndex[1] ) = 1.;
    }
    else if (initialIndex.size() == 1) // 1D
    {
// Using template slice operators
        (*price[0]) = 0.;
        static_cast< TreeSliceEQ & >( (*price[0]) )( initialIndex[0] ) = 1.;
    }
    else
    {
        throw ModelException(method,
                             "The forward induction is not supported for 1D or 3D models.");
    }
}

/** product payoff before maturity date for backward induction */
void VanillaGridProd::prod_BWD(
    const TreeSlice & s,
    int step,
    int pStart,
    int pEnd,
    const vector< TreeSliceSP > & price)
{
    static const string method = "VanillaGridProd::prod_BWD";
    try {
        //test if it's a maturity of some product
        if (stepIsMat[step])
        {
            const DateTimeArray& maturities = *instrVan->maturities;
            int nbMats = maturities.size();
            int iMat = stepMatIndex[step]; // index of the maturity
            if ((iMat < 0) || (iMat >= nbMats))
            {
                throw ModelException(method,
                                     "VanillaGridProd::prodBCDS: The maturities index is inconsistent.");
            }

            //set up the last price we need to calculate
            const DoubleArray& strikes = *instrVan->strikes;
            int nbStrikes = strikes.size();

            pEnd = (nbMats - iMat) * nbStrikes - 1;

            double settlePV = instrVan->instSettle->pvAdjust((*instrVan->maturities.get())[iMat],
                                                             instrVan->discount.get(),
                                                             instrVan->asset.get());

            int iStrk;
            for (iStrk = 0 ; iStrk < nbStrikes ; iStrk++)
            {
// Using template slice operators
                (*price[iStrk + (nbMats - 1 - iMat) * nbStrikes]) =
                    smax( 0., ( s - (*instrVan->strikes)[iStrk] ) * ( instrVan->isCall ? settlePV : -settlePV ) );
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** product payoff after initial date for forward induction */
void VanillaGridProd::prod_FWD(
    const TreeSlice & s,
    int step,
    int pStart,
    int pEnd,
    const vector< TreeSliceSP > & price)
{
    static const string method = "VanillaGridProd::prod_FWD";

    try {
        //test if it's a maturity of some product
        if (stepIsMat[step])
        {
            // maturities
            const DateTimeArray& maturities = *instrVan->maturities;
            int nbMats = maturities.size();
            int iMat = stepMatIndex[step]; // index of the maturity
            if ((iMat < 0) || (iMat >= nbMats))
            {
                throw ModelException(method,
                                     "VanillaGridProd::prodBCDS: The maturities index is inconsistent.");
            }

            // strikes
            const DoubleArray& strikes = *instrVan->strikes;
            const DoubleMatrix& instStrikesUsed = *instrVan->instStrikesUsed;
            int nbStrks = strikes.size();
            int iStrk;

            // weights for the this maturity and these strikes
            const double* weights = (*instrVan->weights)[iMat];

            // compute the total probability (sum should be 1.0 by construction
// Using template slice operators
            double totalProba = esum( (*price[0]), (*price[0]) );

            for (iStrk = 0 ; iStrk < nbStrks ; iStrk++)
            {
                // checks that weights are not 0.0
                if (!weights || Maths::isPositive(weights[iStrk]))
                {
// Using template slice operators
                    pricesFwdBwd[iMat][iStrk] = esum(
                            (*price[0]),
                            (*price[0]) / totalProba * 
                            smax( 0., ( s - instStrikesUsed[iMat][iStrk] ) * ( instrVan->isCall ? 1. : -1. ) ) );
                }
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** create a fd payoff product */
FDProductSP VanillaGrid::createProduct(FDModel* model) const
{
    return FDProductSP( new VanillaGridProd(this, model) );
}

/// ------------  end of FD product -----------


/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* VanillaGrid::createProduct(CClosedFormLN* model) const{
    return new VanillaGridClosedForm(this);
}

#if 0
/** create a tree payoff product */
CTree1f::IProduct* VanillaGrid::createProduct(CTree1f* model) const
{
    if (model->DivsTreatedAsAbsolute() && fwdStarting){
        model->SetDivAmountTreatment(false);
    }

    Vanilla1fProd* treeProd = new Vanilla1fProd(this);
    treeProd->tree1f = model;
    treeProd->model1F = model;

    return treeProd;
}

/** create a fd payoff product */
FD1F::IProduct* VanillaGrid::createProduct(FD1F* model) const
{

    Vanilla1fProd* fdProd = new Vanilla1fProd(this);
    fdProd->fdModel = model;
    fdProd->model1F = model;

    return fdProd;
}

/** indicates whether VEGA_MATRIX is sensible for this instrument */
bool VanillaGrid::avoidVegaMatrix(const IModel* model)
{
    /* this should possibly be false for local vol models etc. */
    return false;
}

/** returns all strikes on the vol surface to which
    this instrument is sensitive */
DoubleArraySP VanillaGrid::getSensitiveStrikes(OutputNameConstSP outputName,
                                               const IModel*      model)
{
    DoubleArraySP sensStrikes = DoubleArraySP(new DoubleArray(0));

    if (avoidVegaMatrix(model)) {
        throw ModelException("VanillaGrid::getSensitiveStrikes",
                             "VEGA_MATRIX is not valid for this instrument");
    }

    // get start date for vol interpolation
    DateTime imntStartDate = fwdStarting ? startDate:valueDate;
    // get last exercise date in exercise schedule
    DateTime maturityDate = exerciseSchedule->lastDate();
    // get last strike in exercise schedule
    double   strike       = exerciseSchedule->lastValue();

    // create a vol request object to be passed on
    LinearStrikeVolRequestSP volRequest(new LinearStrikeVolRequest(strike,
                                                                   imntStartDate,
                                                                   maturityDate,
                                                                   fwdStarting));

    SensitiveStrikeDescriptor sensStrikeDesc;
    sensStrikeDesc.forwardOnly = false;

    asset->getSensitiveStrikes(volRequest.get(), outputName,
                               sensStrikeDesc, sensStrikes);

    return sensStrikes;
}

#endif

/** Rolls the value date and sets initial spot if rolling over start date */
bool VanillaGrid::sensShift(Theta* shift)
{
    throw ModelException("VanillaGrid::sensShift()", "FIXME need to review this first");
    DateTime newDate = shift->rollDate(valueDate);

    if (fwdStarting && newDate.isGreaterOrEqual(startDate) &&
        startDate.isGreaterOrEqual(valueDate))
    {
        fwdStarting = false;
        initialSpot = asset->getThetaSpotOnDate(shift, startDate);
    }

    // roll today
    valueDate = newDate;

    return true;
};

#if 0

double VanillaGrid::priceSpread(const DateTime& valueDate,
                             const DateTime& startDate,
                             const DateTime& matDate,
                             bool isCall,
                             bool fwdStarting,
                             bool oneContract,
                             double notional,
                             double initialSpot,
                             double lowStrike,
                             double highStrike,
                             const InstrumentSettlement* instSettle,
                             const Asset* asset,
                             const YieldCurve* discount)
{
    static const string method = "VanillaGrid::priceSpread";

    /* if it's forward starting, convert payoffToProcFreqBoundWeight strikes
       to absolute values based on spot at start date */
    DateTime imntStartDate = fwdStarting? startDate: valueDate;

    double fwdAtStart;
    double lowAbsStrike = lowStrike;
    double highAbsStrike = highStrike;
    if (fwdStarting)
    {
        fwdAtStart = asset->fwdValue(startDate);
        lowAbsStrike *= fwdAtStart;
        highAbsStrike *= fwdAtStart;
    }

    // get the fwd price - only need to do this once
    double fwdPrice = asset->fwdValue(matDate);

    // choose how to interpolate the vol - go for traditional route for now
    LinearStrikeVolRequestSP lowVolRequest(new LinearStrikeVolRequest(
                                                lowStrike,
                                                imntStartDate,
                                                matDate,
                                                fwdStarting));

    LinearStrikeVolRequestSP highVolRequest(new LinearStrikeVolRequest(
                                                highStrike,
                                                imntStartDate,
                                                matDate,
                                                fwdStarting));

    // interpolate the vol for each strike
    CVolProcessedSP lowVol(asset->getProcessedVol(lowVolRequest.get()));
    CVolProcessedSP highVol(asset->getProcessedVol(highVolRequest.get()));

    // cast to the type of vol we're expecting
    CVolProcessedBSSP lowVolBS = CVolProcessedBSSP::dynamicCast(lowVol);
    CVolProcessedBSSP highVolBS = CVolProcessedBSSP::dynamicCast(highVol);

    // this should never happen if our get market data has worked properly
    if (!lowVol || !highVol)
    {
        throw ModelException(method, "No Black Scholes Vol");
    }

    // calculate the variance
    double lowVariance = lowVolBS->CalcVar(imntStartDate, matDate);
    double highVariance = highVolBS->CalcVar(imntStartDate, matDate);

    // calculate the discount factor back to today
    double discFactor = instSettle->pv(valueDate,
                                       matDate,
                                       discount,
                                       asset);

    // finally call Black model for each strike
    double lowStrikePremium = Black::price(isCall, fwdPrice,
                                           lowAbsStrike, discFactor, lowVariance);

    double highStrikePremium = Black::price(isCall, fwdPrice,
                                            highAbsStrike, discFactor, highVariance);

    double scalingFactor = InstrumentUtil::scalePremium(oneContract,
                                                        fwdStarting,
                                                        notional,
                                                        fwdAtStart,
                                                        initialSpot);


    lowStrikePremium *= scalingFactor;
    highStrikePremium *= scalingFactor;


    // Now combine the premiums
    double spreadPrice = isCall?
        lowStrikePremium - highStrikePremium : /* call spread */
        highStrikePremium - lowStrikePremium;  /* put spread */

    return spreadPrice;
}

CSensControl* VanillaGrid::AlterControl(const IModel*       modelParams,
                                        const CSensControl* sensControl) const
{
    CSensControl* alteredControl = NULL;
    if (Delta::TYPE->isInstance(sensControl)         &&
        CClosedFormLN::TYPE->isInstance(modelParams) )
    {
       const Delta* delta =
            dynamic_cast<const Delta*>((IObject*)sensControl);
        double  strike  = exerciseSchedule->lastValue();
        ShiftSizeCollectorSP shiftSizeVisitor(new ShiftSizeCollector(
            delta,
            strike,
            fwdStarting?
            ShiftSizeCollector::FWD_START_ADJUSTMENT:
            ShiftSizeCollector::SPOT_START_ADJUSTMENT));

        asset->accept(shiftSizeVisitor.get());

        if ( Maths::isPositive(shiftSizeVisitor->getShiftSize()) )
        {
            alteredControl = new Delta(shiftSizeVisitor->getShiftSize());
            alteredControl->
                setMarketDataName(sensControl->getMarketDataName());
        }
    }

    return alteredControl;
}
#endif

// for reflection
VanillaGrid::VanillaGrid():
CInstrument(TYPE),
isCall(false),
fwdStarting(false),
oneContract(true),
notional(0.0),
initialSpot(0.0),
strikeUnits(WeightMatrix::DEFAULT),
instType(WeightMatrix::DEFAULT){}//,
//spotAtMaturity(0.0){}


static const int MC_CALL_FLAG = 0;
static const int MC_PUT_FLAG = 1;
static const int MC_OTM_FLAG =2;
static const int MC_ITM_FLAG =3;


/////////////////////////////////////////////////////////
//           private class for Monte-Carlo product
/////////////////////////////////////////////////////////
/** vanilla grid product payoff for a Monte-Carlo */
class VanillaGridMC : public IMCProduct, virtual public IMCProductLN {
private:
    const VanillaGrid* inst;
    int                nbStrikes;    // nber of strikes
    int                nbMaturities; // nber of maturities
    DoubleArraySP      timeFrac;     // time fractions from value date to maturities
    DoubleArraySP      forwards;     // forwards
    DoubleArraySP      pv;           // pv factors
    double             mult;
    bool               useRatio;

    int                instFlag;


    class VanillaGridPrices;
    typedef refCountPtr<VanillaGridPrices> VanillaGridPricesSP;

    class VanillaGridPrices: public IMCPricesSimple {
    public:

        /** adds supplied price to this set of IMCPrices */
        void add(double price,
                 int    idxRow,
                 int    idxCol){
            simplePrices[idxRow][idxCol]->add(price);
        }

        /** Returns the averaged result, together with the standard error */
        virtual void getResult(double* result,
                               double* resultStdErr) const {
            simplePrices[0][0]->getResult(result, resultStdErr);
        }

        /** Returns the averaged result, together with the standard error */
        virtual void getResult(int     idxRow,
                               int     idxCol,
                               double* result,
                               double* resultStdErr) const {
            simplePrices[idxRow][idxCol]->getResult(result, resultStdErr);
        }

        /** Returns the grid (of size nbMaturities, nbStrikes) of averaged results */
        DoubleArrayArraySP getPricesGrid(int nberMats, int nberStrks) const {
            DoubleArrayArraySP pricesGrid = DoubleArrayArraySP(new DoubleArrayArray(nberMats));
            double dummy; // used to store the std error

            int iMat, iStrk;
            // loop on maturities
            for (iMat = 0 ; iMat < nberMats ; iMat++) {
                (*pricesGrid)[iMat].resize(nberStrks);
                // loop on strikes
                for (iStrk = 0 ; iStrk < nberStrks ; iStrk++) {
                    getResult(iMat, iStrk, &(*pricesGrid)[iMat][iStrk], &dummy);
                }
            }

            return pricesGrid;
        }

        /** Returns the grid (of size nbMaturities, nbStrikes) of standard errors */
        DoubleArrayArraySP getStdErrorsGrid(int nberMats, int nberStrks) const {
            DoubleArrayArraySP stdErrorsGrid = DoubleArrayArraySP(new DoubleArrayArray(nberMats));
            double dummy; // used to store the averaged result

            int iMat, iStrk;
            // loop on maturities
            for (iMat = 0 ; iMat < nberMats ; iMat++) {
                (*stdErrorsGrid)[iMat].resize(nberStrks);
                // loop on strikes
                for (iStrk = 0 ; iStrk < nberStrks ; iStrk++) {
                    getResult(iMat, iStrk, &dummy, &(*stdErrorsGrid)[iMat][iStrk]);
                }
            }
            return stdErrorsGrid;
        }

        /** Returns true if the path, identified by pathIdx, should be
            repriced. If false is returned then there will be no add()
            method called for this path and the IMCPrices object must
            take any appropriate action */
        virtual bool repriceForGreek(int pathIdx){
            return true;
        }

        /** Support product-level caching across tweaks */
        virtual int storagePerPath(IMCProduct* product) const{
            return 0;
        }

        virtual void configureCache(const IntArray& changedAssets){}

        /** Returns a deep copy of this object */
        IMCPrices* clone() const {
            VanillaGridPricesSP copy(new VanillaGridPrices());
            copy->simplePrices.resize(simplePrices.size());

            unsigned int iMat, iStrk;
            // loop on maturities
            for (iMat = 0; iMat < copy->simplePrices.size() ; ++iMat) {
                copy->simplePrices[iMat].resize(simplePrices[iMat].size());
                // loop on strikes
                for (iStrk = 0 ; iStrk < simplePrices[iMat].size() ; ++iStrk) {
                    IMCPricesSP temp(this->simplePrices[iMat][iStrk]->clone());
                    copy->simplePrices[iMat][iStrk] = DYNAMIC_POINTER_CAST<IMCPricesSimple>(temp);
                }
            }
            return copy.get();
        }

        VanillaGridPrices(int nberMats,
                          int nberStrks,
                          int NbIter,
                          int NbSubSamples):
        simplePrices(nberMats) {
            int iMat, iStrk;
            // loop on maturities
            for (iMat = 0 ; iMat < nberMats ; ++iMat) {
                simplePrices[iMat].resize(nberStrks);
                // loop on strikes
                for (iStrk = 0 ; iStrk < nberStrks ; ++iStrk){
                    simplePrices[iMat][iStrk] = IMCPricesSimpleSP(new MCPricesSimple(NbIter, NbSubSamples));
                }
            }
        }

    private:
        VanillaGridPrices() {}

        /** adds supplied price to this set of IMCPrices */
        virtual void add(double price) {
            throw ModelException("IMCPrices::add(double)",
                                 "internal error");
        }

        /** Returns the last price stored. Undefined behaviour if no
            prices yet stored */
        virtual double lastPrice() const {
            throw ModelException("IMCPrices::lastPrice",
                                 "internal error");
        }

        /** On pricing run returns MAX(x, 0.0). It should be used for applying
            the 'final point of optionality' within the product. This allows
            QuickGreeks type of IMCPrices not to apply the max when doing first
            order greeks (this may sound strange but see the docs for why) */
        virtual double maxWithZero(double x) const {
            throw ModelException("IMCPrices::maxWithZero",
                                 "internal error");
        }

        /** Reset this object so that it can be used for the same operation
            again. Normally, a new IMCPrices object is created for each
            pricing run. However, for quick x gamma, it is important to use
            the same one for each pair of assets */
        virtual void reset() {
            throw ModelException("IMCPrices::reset",
                                 "internal error");
        }

        /** Ease cloning */
        virtual IMCPrices* emptyConstructor() const {
            throw ModelException("IMCPrices::emptyConstructor",
                                 "internal error");
        }

        vector<vector<IMCPricesSimpleSP> > simplePrices;
    };

public:

    virtual IMCPrices* createOrigPrices(int nbIter,
                                     int nbSubSamples,
                                     int mode) {
        return new VanillaGridPrices(nbMaturities, nbStrikes, nbIter, nbSubSamples);
    }

    /** invoked after final simulated path is run. Default does nothing.
        Allows derived classes to store debug information for example */
    virtual void recordExtraOutput(Control*      control,
                                   Results*      results,
                                   const IMCPrices& prices) const {
        const VanillaGridPrices& myprices = static_cast<const VanillaGridPrices&>(prices);
        if (control) {

            // options prices
            OutputRequest* request = control->requestsOutput(OutputRequest::OPTION_PRICE);
            if (request) {
                VanillaGrid::OutputSP pricesGrid;
                pricesGrid = VanillaGrid::OutputSP(new VanillaGrid::Output(nbMaturities, nbStrikes));

                DoubleArrayArraySP pricesGridTemp = myprices.getPricesGrid(nbMaturities, nbStrikes);

                int iMat, iStrk;
                // loop on maturities
                for (iMat = 0 ; iMat < nbMaturities ; iMat++) {
                    // loop on strikes
                    for (iStrk = 0 ; iStrk < nbStrikes ; iStrk++) {
                        if (!inst->weights || Maths::isPositive((*inst->weights)[iMat][iStrk])) {
                            pricesGrid->setValue(iMat, iStrk, (*pricesGridTemp)[iMat][iStrk]);
                        }
                        else {
                            pricesGrid->setError(iMat, iStrk, VanillaGrid::Output::ErrorCode::PRICE_NOT_REQUESTED);
                        }
                    }
                }
                results->storeRequestResult(request, pricesGrid);
            }
            // options std errors
            request = control->requestsOutput(OutputRequest::VALUE_STE);
            if (request) {
                VanillaGrid::OutputSP stdErrorsGrid;
                stdErrorsGrid = VanillaGrid::OutputSP(new VanillaGrid::Output(nbMaturities, nbStrikes));

                DoubleArrayArraySP stdErrorsGridTemp = myprices.getStdErrorsGrid(nbMaturities, nbStrikes);

                int iMat, iStrk;
                // loop on maturities
                for (iMat = 0 ; iMat < nbMaturities ; iMat++) {
                    // loop on strikes
                    for (iStrk = 0 ; iStrk < nbStrikes ; iStrk++) {
                        if (!inst->weights || Maths::isPositive((*inst->weights)[iMat][iStrk])) {
                            stdErrorsGrid->setValue(iMat, iStrk, (*stdErrorsGridTemp)[iMat][iStrk]);
                        }
                        else {
                            stdErrorsGrid->setError(iMat, iStrk, VanillaGrid::Output::ErrorCode::PRICE_NOT_REQUESTED);
                        }
                    }
                }
                results->storeRequestResult(request, stdErrorsGrid);
            }
            // options implied volatilities
            request = control->requestsOutput(OutputRequest::IND_VOL);
            if (request) {
                VanillaGrid::OutputSP impVolsGrid;
                impVolsGrid = VanillaGrid::OutputSP(new VanillaGrid::Output(nbMaturities, nbStrikes));

                // forward option prices
                DoubleArrayArraySP pricesGridTemp = myprices.getPricesGrid(nbMaturities, nbStrikes);

                int iMat, iStrk;
                // loop on maturities
                for (iMat = 0 ; iMat < nbMaturities ; iMat++) {
                    // loop on strikes
                    for (iStrk = 0 ; iStrk < nbStrikes ; iStrk++) {
                        if (!inst->weights || Maths::isPositive((*inst->weights)[iMat][iStrk])) {
                            // call version of impliedVariance that returns status
                            double variance;
                            bool isCall = (*(*inst->instTypesUsed)[iMat])[iStrk] ? true : false;
                            if (Black::impliedVariance(isCall,
                                                  (*forwards)[iMat],
                                                  (*inst->instStrikesUsed)[iMat][iStrk],
                                                  (*pv)[iMat],                            // pv
                                                  0.3 * 0.3 * (*timeFrac)[iMat],          // initial var guess
                                                  (*pricesGridTemp)[iMat][iStrk] / mult,
                                                  2.0 * 0.3 * 1.0e-5 * (*timeFrac)[iMat], // var accuracy
                                                  variance)) {
                                impVolsGrid->setValue(iMat, iStrk, sqrt(variance / (*timeFrac)[iMat]));
                            }
                            else {
                                impVolsGrid->setError(iMat, iStrk, VanillaGrid::Output::ErrorCode::IMP_VOL_ERROR);
                            }
                        }
                        else {
                            impVolsGrid->setError(iMat, iStrk, VanillaGrid::Output::ErrorCode::IMP_VOL_NOT_REQUESTED);
                        }
                    }
                }
                results->storeRequestResult(request, impVolsGrid);
            }
            // option vegas: don't want to support vega request (only useful for closed form)
            request = control->requestsOutput(OutputRequest::OPTION_VEGA);
            if (request){
                results->storeRequestResult(request, IObjectSP(new NotApplicable()));
            }
        }
    }

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen) const {
        // empty
    }

    // equivalent to InstIntoMCProduct
    VanillaGridMC(const VanillaGrid*        inst,
                  const IRefLevelConstSP&   refLevel,      // how to 'avg in'
                  const SimSeriesConstSP&   simSeries,     // simulation dates
                  const IPastValuesConstSP& mcPastValues): // historic values
    IMCProduct(inst->asset.get(),
        inst->valueDate,
        inst->discount.get(),
        refLevel,
        simSeries,
        mcPastValues,
        inst->instSettle.get(),
        simSeries->getLastDate()),
    inst(inst),
    nbStrikes(inst->strikes->size()),
    nbMaturities(inst->maturities->size()) {

        static const string method = "VanillaGridMC::VanillaGridMC";
            try{

            if(CString::equalsIgnoreCase(inst->instType, WeightMatrix::CALL)) {
                instFlag = MC_CALL_FLAG;
            } else if(CString::equalsIgnoreCase(inst->instType, WeightMatrix::PUT)) {
                instFlag = MC_PUT_FLAG;
            } else if(CString::equalsIgnoreCase(inst->instType, WeightMatrix::OTM)) {
                instFlag = MC_OTM_FLAG;
            } else if(CString::equalsIgnoreCase(inst->instType, WeightMatrix::ITM)) {
                instFlag = MC_ITM_FLAG;
            } else {
                 throw ModelException("instType is not recognized");
            }

            if (inst->fwdStarting) {
                // forward starting => notional based, and % strike
                useRatio = true;
                mult = inst->notional;
            } else {
                useRatio = false;
                mult = inst->oneContract ? 1.0 : inst->notional / inst->initialSpot;
            }

            // time metric (time fractions from value date to maturities)
            VolRequestTimeSP volTimeReq(new VolRequestTime());
            CVolProcessedSP volTimeProcessed(getMultiFactors()->factorGetProcessedVol(0, volTimeReq.get()));
            TimeMetricConstSP timeMetric = volTimeProcessed->GetTimeMetric();

            /*
            ATMVolRequestSP requestATM(new ATMVolRequest());
            requestATM->allowNegativeFwdVar(false);
            CVolProcessedSP volTemp(getMultiFactors()->factorGetProcessedVol(0, requestATM.get()));
            CVolProcessedBSSP volTempBS(CVolProcessedBSSP::dynamicCast(volTemp));

            TimeMetricConstSP timeMetric = TimeMetricConstSP(volTempBS->GetTimeMetric().get());*/

            timeFrac = DoubleArraySP(new DoubleArray(nbMaturities));
            // loop on maturities
            int iMat;
            for(iMat = 0; iMat < nbMaturities ; iMat++ ) {
                (*timeFrac)[iMat] = timeMetric->yearFrac(inst->valueDate, (*inst->maturities)[iMat]);
            }

            // pv factors
            pv = DoubleArraySP(new DoubleArray(nbMaturities));
            // loop on maturities
            for (iMat = 0 ; iMat < nbMaturities ; iMat++) {
                (*pv)[iMat] = inst->instSettle->pv(inst->valueDate,
                                              (*inst->maturities)[iMat],
                                              inst->discount.get(),
                                              inst->asset.get());
            }

            // forwards
            forwards = DoubleArraySP(new DoubleArray(nbMaturities));
            if (inst->fwdStarting) {
                DateTimeArray dates(1 + nbMaturities);
                dates[0] = inst->startDate;
                int iMat;
                int iDate = 1;
                // loop on maturities
                for (iMat = 0 ; iMat < nbMaturities ; iMat = iDate++) {
                    dates[iDate] = (*inst->maturities)[iMat];
                }
                CDoubleArray myFwds(dates.size());
                mAsset->assetFwdValue(0, dates, myFwds);
                // loop on maturities
                for (iMat = 0, iDate = 1 ; iMat < nbMaturities ; iMat = iDate++) {
                    (*forwards)[iMat] = myFwds[iDate] / myFwds[0];
                }
            }
            else {
                mAsset->assetFwdValue(0, (*inst->maturities), *forwards);
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&     prices) {
        static const string routine("VanillaGridMC::payoff");
        try {
            VanillaGridPrices& myprices = static_cast<VanillaGridPrices&>(prices);
            // average in, average out value and price
            double inValue = pathGen->refLevel(0,0);
            double outValue, price;
            // index of the first maturity and the last maturity in the path generator
            int beginIdx = pathGen->begin(0); // only one asset
            int endIdx   = pathGen->end(0);   // only one asset

            const double* path = pathGen->Path(0,0);

            int iMat, iStrk;
            bool isZero;
            const DoubleMatrix& instStrikesUsed = *inst->instStrikesUsed;
            const IntArrayArray& instTypesUsed  = *inst->instTypesUsed;

            // loop on maturities
            for (iMat = beginIdx ; iMat < endIdx ; iMat++) {
                outValue = path[iMat];

                const double* weights = (*inst->weights)[iMat];
                const double* strikesUsed = &(*inst->instStrikesUsed)[iMat][0];

                isZero = false;

                // loop on strikes for calls
                if (instFlag == MC_CALL_FLAG)
                {
                    // loop for increasing strike for calls
                    for (iStrk = 0 ; iStrk < nbStrikes ; iStrk++) {
                        if (!weights || Maths::isPositive(weights[iStrk])) {
                            if (!isZero)
                            {
                                price = useRatio ?
                                    (outValue / inValue - strikesUsed[iStrk]) :
                                        (outValue - strikesUsed[iStrk]);

                                // update isZero flag
                                if (price < 0.0) {
                                    isZero = true;
                                }

                                price = price < 0.0 ? 0.0 : (price * mult);
                            }
                            // if the price is zero, all the calls with greater strikes are worth zero
                            else {
                                price = 0.0;
                            }
                            // line below to be reviewed
                            // here for consistency with the SV compliant MC product
                            if (!doingPast()) {
                                price *= (*pv)[iMat];
                            }
                            myprices.add(price, iMat - beginIdx, iStrk);
                        }
                    }
                } else if (instFlag == MC_PUT_FLAG) {
                    // loop for increasing strike for calls
                    for (iStrk = nbStrikes-1 ; iStrk >= 0 ; iStrk--) {
                        if (!weights || Maths::isPositive(weights[iStrk])) {
                            if (!isZero){
                                price = useRatio ?
                                        (strikesUsed[iStrk] - outValue / inValue) :
                                            (strikesUsed[iStrk] - outValue);

                                // update isZero flag
                                if (price < 0.0) {
                                    isZero = true;
                                }

                                price = price < 0.0 ? 0.0 : (price * mult);
                            }
                            // if the price is zero, all the puts with greater strikes are worth zero
                            else {
                                price = 0.0;
                             }
                             // line below to be reviewed
                             // here for consistency with the SV compliant MC product
                             if (!doingPast()) {
                                price *= (*pv)[iMat];
                             }
                             myprices.add(price, iMat - beginIdx, iStrk);
                        }
                    }
                } else if ((instFlag == MC_OTM_FLAG)
                            ||(instFlag == MC_ITM_FLAG)) {
                    // loop for increasing strike for other
                    for (iStrk = 0 ; iStrk < nbStrikes ; iStrk++) {
                        if (!weights || Maths::isPositive(weights[iStrk])) {
                            price = useRatio ?
                                (outValue / inValue - strikesUsed[iStrk]) :
                                    (outValue - strikesUsed[iStrk]);
                            if(!(*instTypesUsed[iMat])[iStrk]) {
                                price *= -1.0;
                            }
                            price = price < 0.0 ? 0.0 : (price * mult);

                            // line below to be reviewed
                            // here for consistency with the SV compliant MC product
                            if (!doingPast()) {
                                price *= (*pv)[iMat];
                            }
                            myprices.add(price, iMat - beginIdx, iStrk);
                        }
                    }
                } else {
                    throw ModelException("insType not recognized");
                }
            }
        }
        catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGenerator,
                                    int                     iAsset) const {
        CVolRequestLNArray reqarr(1);
        double interpLevel;
        DateTime imntStartDate;

        // atm vol request
        if (inst->fwdStarting) {
            interpLevel = 1.0;
            imntStartDate = inst->startDate;
        }
        else {
            interpLevel = inst->initialSpot;
            imntStartDate = inst->valueDate;
        }

        reqarr[0] = CVolRequestLNSP(new LinearStrikeVolRequest(
            interpLevel,              // atm
            imntStartDate,            // start date
            inst->maturities->back(),
            inst->fwdStarting));

        return reqarr;
    }

};


///////////////////////////////////////////////////////////////////////////////////
//           private class for state var compliant Monte-Carlo product
///////////////////////////////////////////////////////////////////////////////////
/** vanilla grid product payoff for a Monte-Carlo */
class VanillaGridMCSV : public MCProductClient, virtual public IMCProductLN {
private:
    const VanillaGrid* inst;
    int                nbStrikes;    // nber of strikes
    int                nbMaturities; // nber of maturities
    DoubleArraySP      timeFrac;     // time fractions from value date to maturities
    DoubleArraySP      forwards;     // forwards
    DoubleArraySP      pv;           // pv factors
    double             mult;
    bool               useRatio;

    int                instFlag;

    // state vars
    SVGenSpotSP                  spotGen;     // generator for spot
    SVGenSpot::IStateVarSP       spotSV;      // spot state variable
    IRefLevel::IStateVarGenSP refLevelGen; // generator for ref level
    IRefLevel::IStateVarSP    refLevelSV;  // ref level state variable
    SVGenDiscFactorSP            dfGen;       // generator for discount factors
    SVDiscFactorSP dfSV;        // df state variable

    /** Override default method on IMCProduct. This method is called every time
        the path generator is changed (which is, at the moment, when the
        past path generator is created, and then when the future path
        generator is created  */
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen) {
        static const string routine = "VanillaGridMCSV::pathGenUpdated";

        try {
            spotSV = spotGen->getSpotSV(newPathGen);
            refLevelSV = refLevelGen->getRefLevelSV(refLevelSV, newPathGen);
            dfSV = dfGen->getSVDiscFactor(newPathGen);
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    };

    class VanillaGridPrices;
    typedef refCountPtr<VanillaGridPrices> VanillaGridPricesSP;

    class VanillaGridPrices: public IMCPricesSimple {
    public:

        /** adds supplied price to this set of IMCPrices */
        void add(double price,
                 int    idxRow,
                 int    idxCol){
            simplePrices[idxRow][idxCol]->add(price);
        }

        /** Returns the averaged result, together with the standard error */
        virtual void getResult(double* result,
                               double* resultStdErr) const {
            simplePrices[0][0]->getResult(result, resultStdErr);
        }

        /** Returns the averaged result, together with the standard error */
        virtual void getResult(int     idxRow,
                               int     idxCol,
                               double* result,
                               double* resultStdErr) const {
            simplePrices[idxRow][idxCol]->getResult(result, resultStdErr);
        }

        /** Returns the grid (of size nbMaturities, nbStrikes) of averaged results */
        DoubleArrayArraySP getPricesGrid(int nberMats, int nberStrks) const {
            DoubleArrayArraySP pricesGrid = DoubleArrayArraySP(new DoubleArrayArray(nberMats));
            double dummy; // used to store the std error

            int iMat, iStrk;
            // loop on maturities
            for (iMat = 0 ; iMat < nberMats ; iMat++) {
                (*pricesGrid)[iMat].resize(nberStrks);
                // loop on strikes
                for (iStrk = 0 ; iStrk < nberStrks ; iStrk++) {
                    getResult(iMat, iStrk, &(*pricesGrid)[iMat][iStrk], &dummy);
                }
            }

            return pricesGrid;
        }

        /** Returns the grid (of size nbMaturities, nbStrikes) of standard errors */
        DoubleArrayArraySP getStdErrorsGrid(int nberMats, int nberStrks) const {
            DoubleArrayArraySP stdErrorsGrid = DoubleArrayArraySP(new DoubleArrayArray(nberMats));
            double dummy; // used to store the averaged result

            int iMat, iStrk;
            // loop on maturities
            for (iMat = 0 ; iMat < nberMats ; iMat++) {
                (*stdErrorsGrid)[iMat].resize(nberStrks);
                // loop on strikes
                for (iStrk = 0 ; iStrk < nberStrks ; iStrk++) {
                    getResult(iMat, iStrk, &dummy, &(*stdErrorsGrid)[iMat][iStrk]);
                }
            }
            return stdErrorsGrid;
        }

        /** Returns true if the path, identified by pathIdx, should be
            repriced. If false is returned then there will be no add()
            method called for this path and the IMCPrices object must
            take any appropriate action */
        virtual bool repriceForGreek(int pathIdx){
            return true;
        }

        /** Support product-level caching across tweaks */
        virtual int storagePerPath(IMCProduct* product) const{
            return 0;
        }

        virtual void configureCache(const IntArray& changedAssets){}

        /** Returns a deep copy of this object */
        IMCPrices* clone() const {
            VanillaGridPricesSP copy(new VanillaGridPrices());
            copy->simplePrices.resize(simplePrices.size());

            unsigned int iMat, iStrk;
            // loop on maturities
            for (iMat = 0; iMat < copy->simplePrices.size() ; ++iMat) {
                copy->simplePrices[iMat].resize(simplePrices[iMat].size());
                // loop on strikes
                for (iStrk = 0 ; iStrk < simplePrices[iMat].size() ; ++iStrk) {
                    IMCPricesSP temp(this->simplePrices[iMat][iStrk]->clone());
                    copy->simplePrices[iMat][iStrk] = DYNAMIC_POINTER_CAST<IMCPricesSimple>(temp);
                }
            }
            return copy.get();
        }

        VanillaGridPrices(int nberMats,
                          int nberStrks,
                          int NbIter,
                          int NbSubSamples):
        simplePrices(nberMats) {
            int iMat, iStrk;
            // loop on maturities
            for (iMat = 0 ; iMat < nberMats ; ++iMat) {
                simplePrices[iMat].resize(nberStrks);
                // loop on strikes
                for (iStrk = 0 ; iStrk < nberStrks ; ++iStrk){
                    simplePrices[iMat][iStrk] = IMCPricesSimpleSP(new MCPricesSimple(NbIter, NbSubSamples));
                }
            }
        }

    private:
        VanillaGridPrices() {}

        /** adds supplied price to this set of IMCPrices */
        virtual void add(double price) {
            throw ModelException("IMCPrices::add(double)",
                                 "internal error");
        }

        /** Returns the last price stored. Undefined behaviour if no
            prices yet stored */
        virtual double lastPrice() const {
            throw ModelException("IMCPrices::lastPrice",
                                 "internal error");
        }

        /** On pricing run returns MAX(x, 0.0). It should be used for applying
            the 'final point of optionality' within the product. This allows
            QuickGreeks type of IMCPrices not to apply the max when doing first
            order greeks (this may sound strange but see the docs for why) */
        virtual double maxWithZero(double x) const {
            throw ModelException("IMCPrices::maxWithZero",
                                 "internal error");
        }

        /** Reset this object so that it can be used for the same operation
            again. Normally, a new IMCPrices object is created for each
            pricing run. However, for quick x gamma, it is important to use
            the same one for each pair of assets */
        virtual void reset() {
            throw ModelException("IMCPrices::reset",
                                 "internal error");
        }

        /** Ease cloning */
        virtual IMCPrices* emptyConstructor() const {
            throw ModelException("IMCPrices::emptyConstructor",
                                 "internal error");
        }

        vector<vector<IMCPricesSimpleSP> > simplePrices;
    };

public:

    virtual IMCPrices* createOrigPrices(int nbIter,
                                     int nbSubSamples,
                                     int mode) {
        return new VanillaGridPrices(nbMaturities, nbStrikes, nbIter, nbSubSamples);
    }

    /** invoked after final simulated path is run. Default does nothing.
        Allows derived classes to store debug information for example */
    virtual void recordExtraOutput(Control*      control,
                                   Results*      results,
                                   const IMCPrices& prices) const {
        const VanillaGridPrices& myprices = static_cast<const VanillaGridPrices&>(prices);
        if (control) {

            // options prices
            OutputRequest* request = control->requestsOutput(OutputRequest::OPTION_PRICE);
            if (request) {
                VanillaGrid::OutputSP pricesGrid;
                pricesGrid = VanillaGrid::OutputSP(new VanillaGrid::Output(nbMaturities, nbStrikes));

                DoubleArrayArraySP pricesGridTemp = myprices.getPricesGrid(nbMaturities, nbStrikes);

                int iMat, iStrk;
                // loop on maturities
                for (iMat = 0 ; iMat < nbMaturities ; iMat++) {
                    // loop on strikes
                    for (iStrk = 0 ; iStrk < nbStrikes ; iStrk++) {
                        if (!inst->weights || Maths::isPositive((*inst->weights)[iMat][iStrk])) {
                            pricesGrid->setValue(iMat, iStrk, (*pricesGridTemp)[iMat][iStrk]);
                        }
                        else {
                            pricesGrid->setError(iMat, iStrk, VanillaGrid::Output::ErrorCode::PRICE_NOT_REQUESTED);
                        }
                    }
                }
                results->storeRequestResult(request, pricesGrid);
            }
            // options std errors
            request = control->requestsOutput(OutputRequest::VALUE_STE);
            if (request) {
                VanillaGrid::OutputSP stdErrorsGrid;
                stdErrorsGrid = VanillaGrid::OutputSP(new VanillaGrid::Output(nbMaturities, nbStrikes));

                DoubleArrayArraySP stdErrorsGridTemp = myprices.getStdErrorsGrid(nbMaturities, nbStrikes);

                int iMat, iStrk;
                // loop on maturities
                for (iMat = 0 ; iMat < nbMaturities ; iMat++) {
                    // loop on strikes
                    for (iStrk = 0 ; iStrk < nbStrikes ; iStrk++) {
                        if (!inst->weights || Maths::isPositive((*inst->weights)[iMat][iStrk])) {
                            stdErrorsGrid->setValue(iMat, iStrk, (*stdErrorsGridTemp)[iMat][iStrk]);
                        }
                        else {
                            stdErrorsGrid->setError(iMat, iStrk, VanillaGrid::Output::ErrorCode::PRICE_NOT_REQUESTED);
                        }
                    }
                }
                results->storeRequestResult(request, stdErrorsGrid);
            }
            // options implied volatilities
            request = control->requestsOutput(OutputRequest::IND_VOL);
            if (request) {
                VanillaGrid::OutputSP impVolsGrid;
                impVolsGrid = VanillaGrid::OutputSP(new VanillaGrid::Output(nbMaturities, nbStrikes));

                // forward option prices
                DoubleArrayArraySP pricesGridTemp = myprices.getPricesGrid(nbMaturities, nbStrikes);

                int iMat, iStrk;
                // loop on maturities
                for (iMat = 0 ; iMat < nbMaturities ; iMat++) {
                    // loop on strikes
                    for (iStrk = 0 ; iStrk < nbStrikes ; iStrk++) {
                        if (!inst->weights || Maths::isPositive((*inst->weights)[iMat][iStrk])) {
                            // call version of impliedVariance that returns status
                            bool isCall = (*(*inst->instTypesUsed)[iMat])[iStrk] ? true : false;
                            double variance;
                            if (Black::impliedVariance(isCall,
                                                  (*forwards)[iMat],
                                                  (*inst->instStrikesUsed)[iMat][iStrk],
                                                  (*pv)[iMat],                            // pv
                                                  0.3 * 0.3 * (*timeFrac)[iMat],          // initial var guess
                                                  (*pricesGridTemp)[iMat][iStrk] / mult,
                                                  2.0 * 0.3 * 1.0e-5 * (*timeFrac)[iMat], // var accuracy
                                                  variance)) {
                                impVolsGrid->setValue(iMat, iStrk, sqrt(variance / (*timeFrac)[iMat]));
                            }
                            else {
                                impVolsGrid->setError(iMat, iStrk, VanillaGrid::Output::ErrorCode::IMP_VOL_ERROR);
                            }
                        }
                        else {
                            impVolsGrid->setError(iMat, iStrk, VanillaGrid::Output::ErrorCode::IMP_VOL_NOT_REQUESTED);
                        }
                    }
                }
                results->storeRequestResult(request, impVolsGrid);
            }
            // option vegas: don't want to support vega request (only useful for closed form)
            request = control->requestsOutput(OutputRequest::OPTION_VEGA);
            if (request){
                results->storeRequestResult(request, IObjectSP(new NotApplicable()));
            }
        }
    }

    /** Appends 'true' (ie non derived) state variable generators
        required to the supplied collector.*/
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const {
        // ask for a reference level State Variable
        svCollector->append(refLevelGen.get());
        svCollector->append(spotGen.get());
        svCollector->append(dfGen.get());
    }

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen) const {
        // empty
    }

    // equivalent to InstIntoMCProduct
    VanillaGridMCSV(const VanillaGrid*        inst,
                    const IRefLevelConstSP&   refLevel,      // how to 'avg in'
                    const SimSeriesConstSP&   simSeries,     // simulation dates
                    const IPastValuesConstSP& mcPastValues): // historic values
    MCProductClient(inst->asset.get(),
                    inst->valueDate,
                    inst->discount.get(),
                    refLevel,
                    simSeries,
                    mcPastValues,
                    inst->instSettle.get(),
                    simSeries->getLastDate()),
    inst(inst),
    nbStrikes(inst->strikes->size()),
    nbMaturities(inst->maturities->size()),
    spotGen(new SVGenSpot(simSeries)),
    refLevelGen(refLevel->createStateVarGen(getMultiFactors(), inst->valueDate)),
    dfGen(new SVGenDiscFactor(inst->valueDate,
                           inst->discount.getSP(),
                           inst->instSettle,
                           *inst->maturities)) {

        static const string method = "VanillaGridMCSV::VanillaGridMCSV";
        try{

            if(CString::equalsIgnoreCase(inst->instType, WeightMatrix::CALL)) {
                instFlag = MC_CALL_FLAG;
            } else if(CString::equalsIgnoreCase(inst->instType, WeightMatrix::PUT)) {
                instFlag = MC_PUT_FLAG;
            } else if(CString::equalsIgnoreCase(inst->instType, WeightMatrix::OTM)) {
                instFlag = MC_OTM_FLAG;
            } else if(CString::equalsIgnoreCase(inst->instType, WeightMatrix::ITM)) {
                instFlag = MC_ITM_FLAG;
            } else {
                throw ModelException("instType is not recognized");
            }

            if (inst->fwdStarting) {
                // forward starting => notional based, and % strike
                useRatio = true;
                mult = inst->notional;
            } else {
                useRatio = false;
                mult = inst->oneContract ? 1.0 : inst->notional / inst->initialSpot;
            }

            // time metric (time fractions from value date to maturities)
            VolRequestTimeSP volTimeReq(new VolRequestTime());
            CVolProcessedSP volTimeProcessed(getMultiFactors()->factorGetProcessedVol(0, volTimeReq.get()));
            TimeMetricConstSP timeMetric = volTimeProcessed->GetTimeMetric();

            /*
            ATMVolRequestSP requestATM(new ATMVolRequest());
            requestATM->allowNegativeFwdVar(false);
            CVolProcessedSP volTemp(getMultiFactors()->factorGetProcessedVol(0, requestATM.get()));
            CVolProcessedBSSP volTempBS(CVolProcessedBSSP::dynamicCast(volTemp));

            TimeMetricConstSP timeMetric = TimeMetricConstSP(volTempBS->GetTimeMetric().get());*/

            timeFrac = DoubleArraySP(new DoubleArray(nbMaturities));
            // loop on maturities
            int iMat;
            for (iMat = 0; iMat < nbMaturities ; iMat++ ) {
                (*timeFrac)[iMat] = timeMetric->yearFrac(inst->valueDate, (*inst->maturities)[iMat]);
            }

            // pv factors
            pv = DoubleArraySP(new DoubleArray(nbMaturities));
            // loop on maturities
            for (iMat = 0 ; iMat < nbMaturities ; iMat++) {
                (*pv)[iMat] = inst->instSettle->pv(inst->valueDate,
                                                (*inst->maturities)[iMat],
                                              inst->discount.get(),
                                              inst->asset.get());
            }

            // forwards
            forwards = DoubleArraySP(new DoubleArray(nbMaturities));
            if (inst->fwdStarting) {
                DateTimeArray dates(1 + nbMaturities);
                dates[0] = inst->startDate;
                int iMat;
                int iDate = 1;
                // loop on maturities
                for (iMat = 0 ; iMat < nbMaturities ; iMat = iDate++) {
                    dates[iDate] = (*inst->maturities)[iMat];
                }
                CDoubleArray myFwds(dates.size());
                mAsset->assetFwdValue(0, dates, myFwds);
                // loop on maturities
                for (iMat = 0, iDate = 1 ; iMat < nbMaturities ; iMat = iDate++) {
                    (*forwards)[iMat] = myFwds[iDate] / myFwds[0];
                }
            }
            else {
                mAsset->assetFwdValue(0, (*inst->maturities), *forwards);
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&     prices) {
        static const string routine("VanillaGridMCSV::payoff");
        try {
            VanillaGridPrices& myprices = static_cast<VanillaGridPrices&>(prices);
            // average in, average out value and price
            double inValue = refLevelSV->refLevel(0);
            double outValue, price;
            // index of the first maturity and the last maturity in the path generator
            int beginIdx = spotSV->path(0).begin(); // only one asset
            int endIdx   = spotSV->path(0).end();   // only one asset

            const SVPath& path = spotSV->path(0);

            int iMat, iStrk;
            bool isZero;
            const DoubleMatrix& instStrikesUsed = *inst->instStrikesUsed;
            const IntArrayArray& instTypesUsed  = *inst->instTypesUsed;

            // loop on maturities
            for (iMat = beginIdx ; iMat < endIdx ; iMat++) {
                outValue = path[iMat];

                const double* weights = (*inst->weights)[iMat];
                const double* strikesUsed = &(*inst->instStrikesUsed)[iMat][0];

                isZero = false;

                // loop on strikes for calls
                if (instFlag == MC_CALL_FLAG)
                {
                    // loop for increasing strike for calls
                    for (iStrk = 0 ; iStrk < nbStrikes ; iStrk++) {
                        if (!weights || Maths::isPositive(weights[iStrk])) {
                            if (!isZero)
                            {
                                price = useRatio ?
                                    (outValue / inValue - strikesUsed[iStrk]) :
                                        (outValue - strikesUsed[iStrk]);

                                // update isZero flag
                                if (price < 0.0) {
                                    isZero = true;
                                }

                                price = price < 0.0 ? 0.0 : (price * mult);
                            }
                            // if the price is zero, all the calls with greater strikes are worth zero
                            else {
                                price = 0.0;
                            }
                            // line below to be reviewed. The issue being that
                            // (currently) we can't access values in the future when
                            // doing the past. NB apply df before vanillaReprice->store()
                            if (!doingPast()) {
                                price *= dfSV->path()[iMat];
                            }
                            myprices.add(price, iMat - beginIdx, iStrk);
                        }
                    }
                } else
                // loop on strikes for puts
                if (instFlag == MC_PUT_FLAG)
                {
                    // loop for increasing strike for calls
                    for (iStrk = nbStrikes-1 ; iStrk >= 0 ; iStrk--) {
                        if (!weights || Maths::isPositive(weights[iStrk])) {
                            if (!isZero)
                            {
                                price = useRatio ?
                                    (strikesUsed[iStrk] - outValue / inValue) :
                                    (strikesUsed[iStrk] - outValue);

                                // update isZero flag
                                if (price < 0.0) {
                                    isZero = true;
                                }
                                    price = price < 0.0 ? 0.0 : (price * mult);
                            }
                                // if the price is zero, all the puts with greater strikes are worth zero
                            else {
                                price = 0.0;
                            }
                            // line below to be reviewed. The issue being that (currently)
                            // we can't access values in the future when doing the past
                            if (!doingPast()) {
                                price *= dfSV->path()[iMat];
                            }
                            myprices.add(price, iMat - beginIdx, iStrk);
                        }
                    }
                } else if ((instFlag == MC_OTM_FLAG)
                            ||(instFlag == MC_ITM_FLAG)){
                    // loop for increasing strike for calls
                    for (iStrk = 0 ; iStrk < nbStrikes ; iStrk++) {
                        if (!weights || Maths::isPositive(weights[iStrk])) {
                            price = useRatio ?
                                (outValue / inValue - strikesUsed[iStrk]) :
                                (outValue - strikesUsed[iStrk]);

                            if(!(*instTypesUsed[iMat])[iStrk]) {
                                price *= -1.0;
                            }

                            price = price < 0.0 ? 0.0 : (price * mult);

                            // line below to be reviewed. The issue being that
                            // (currently) we can't access values in the future when
                            // doing the past. NB apply df before vanillaReprice->store()
                            if (!doingPast()) {
                                price *= dfSV->path()[iMat];
                            }
                            myprices.add(price, iMat - beginIdx, iStrk);
                        }
                    }
                } else {
                    throw ModelException("instType is not recognized");
                }
            }
        }
        catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGenerator,
                                    int                     iAsset) const {
        CVolRequestLNArray reqarr(1);
        double interpLevel;
        DateTime imntStartDate;

        // atm vol request
        if (inst->fwdStarting) {
            interpLevel = 1.0;
            imntStartDate = inst->startDate;
        }
        else {
            interpLevel = inst->initialSpot;
            imntStartDate = inst->valueDate;
        }

        reqarr[0] = CVolRequestLNSP(new LinearStrikeVolRequest(
            interpLevel,              // atm
            imntStartDate,            // start date
            inst->maturities->back(),
            inst->fwdStarting));

        return reqarr;
    }

};

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* VanillaGrid::createProduct(const MonteCarlo* model) const {

    static const string routine("VanillaGrid::createProduct");

    // v simple simSeries
    SimSeriesSP simSeries(new SimSeries(1));
    simSeries->addDates(*maturities);

    // v simple RefLevel
    IRefLevelSP refLevel(IRefLevel::Util::makeFwdStart(startDate));

    // v simple PastValues
    IPastValuesSP pastValues(IPastValues::Util::makeTrivial(valueDate, initialSpot));

    // if state vars requested
    if (model->stateVarUsed()) {
        return new VanillaGridMCSV(this, refLevel, simSeries, pastValues);
    }
    // otherwise, use old methodology
    return new VanillaGridMC(this, refLevel, simSeries, pastValues);
}

class VanillaGridISAPHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(VanillaGridISAP, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(FourierEngine::ISAP);
        EMPTY_SHELL_METHOD(defaultVanillaGridISAP);
        FIELD(useOneIntegral, "Use 1-integral approach if true; use 2-integral approach otherwise");
        FIELD_MAKE_OPTIONAL(useOneIntegral);
        FIELD(payoffToProcFreqBoundWeight, "Determines the imaginary line of integration. "
                                                  "The closer to 1.0, the closer to the payoff's relevant frequency boundary");
        FIELD_MAKE_OPTIONAL(payoffToProcFreqBoundWeight);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultVanillaGridISAP(){
        return new VanillaGridISAP();
    }
};

CClassConstSP const VanillaGridISAP::TYPE = CClass::registerClassLoadMethod(
    "VanillaGridISAP", typeid(VanillaGridISAP), VanillaGridISAPHelper::load);

/** Single integrand - Default integrator case */
template <class Process, class Product>
class VanillaGridSingleIntegrand: public Function1DDouble {
public:
    VanillaGridSingleIntegrand(const Process&   process,
                    const Product&   product,
                    double           strike,
                    double           fwd,
                    const DateTime&  maturity,
                    double           omega):
    Function1DDouble(Range(OpenBoundary(0.0), Infinity(Infinity::Plus))),
    logMoneyness(log(strike/fwd)),
    matdate(maturity),
    omega(omega),
    process(process),
    product(product){}

    virtual double operator()(double  u) const {    // u == frequency
        Complex  z(omega, u);
        Complex  Laplace = exp(process.scalelessCumulant(product, z, matdate) - z * logMoneyness) / (z * (z - 1.0));
        return Laplace.real();
    }

private:
    double logMoneyness;
    const DateTime& matdate;
    double omega;
    const Process& process;
    const Product& product;
};

/** Single integrand - FFT case */
template <class Process, class Product>
class VanillaGridSingleIntegrandFFT: public Function1DComplex {
public:
    VanillaGridSingleIntegrandFFT(const Process&   process,
                       const Product&   product,
                       const DateTime&  maturity,
                       double           omega):
    // Infinite range by default
    matdate(maturity),
    omega(omega),
    process(process),
    product(product){}

    const Range& getInterval() const {return this->interval;}

    virtual Complex operator()(double  u) const {    // u == frequency
        Complex  z(omega, u);
        Complex  Laplace = exp(process.scalelessCumulant(product, z, matdate)) / (z * (z - 1.0));
        return Laplace;
    }

private:
    const DateTime& matdate;
    double omega;
    const Process& process;
    const Product& product;
};

/** Double integrand - Default integrator case */
template <class Process, class Product>
class VanillaGridDoubleIntegrand: public Function1DDouble {
public:
    VanillaGridDoubleIntegrand(const Process&   process,
                    const Product&   product,
                    double           strike,
                    double           fwd,
                    const DateTime&  maturity,
                    int              j):    // 0 or 1
    Function1DDouble(Range(OpenBoundary(0.0), Infinity(Infinity::Plus))),
    logMoneyness(log(strike/fwd)),
    matdate(maturity),
    j(j),
    process(process),
    product(product){}

    virtual double operator()(double  u) const {    // u == frequency
        Complex  z(0.0, u); // = i * u
        Complex  Laplace = exp(process.scalelessCumulant(product, z + j, matdate) - z * logMoneyness) / z;
        return Laplace.real();
    }

private:
    double logMoneyness;
    const DateTime& matdate;
    double j;
    const Process& process;
    const Product& product;
};

/** Double integrand - FFT integrator case */
template <class Process, class Product>
class VanillaGridDoubleIntegrandFFT: public Function1DComplex {
public:
    VanillaGridDoubleIntegrandFFT(const Process&   process,
                       const Product&   product,
                       const DateTime&  maturity,
                       int              j):    // 0 or 1
    // Infinite range by default
    matdate(maturity),
    j(j),
    process(process),
    product(product){}

    virtual Complex operator()(double  u) const {    // u == frequency
        if (Maths::isZero(u)){
            throw ModelException("VanillaGridDoubleIntegrandFFT::operator(double)",
                                 "Zero frequency is not supported yet");
        }
        Complex  z(0.0, u); // = i * u
        Complex  Laplace = exp(process.scalelessCumulant(product, z + j, matdate)) / z;
        return Laplace;
    }

private:
    const DateTime& matdate;
    double j;
    const Process& process;
    const Product& product;
};

/** Fourier Product */
class VanillaGridFP: public FourierProduct,
                     public FourierProductIntegrator1D,     // supports default 1D integrators
                     public FourierProductFFTIntegrator1D,  // supports FFT 1D integrators
                     public StFourierProductLogRtn,         // requires a StFourierProcessLogRtn
                     public FwdStFourierProductLogRtn{      // requires a FwdStFourierProcessLogRtn
public:
    // equivalent to InstIntoFourierProduct
    VanillaGridFP(const VanillaGrid* inst):    // maturity date
    FourierProduct(inst->asset.get(),
                   inst->valueDate,
                   inst->discount.get(),
                   inst->instSettle.get()),
    inst(inst),
    indVolRequest(0),
    optionPriceRequest(0),
    optionVegaRequest(0){
        if (inst->fwdStarting) {
            // forward starting => notional based, and % strike
            mult = inst->notional;
        }
        else {
            mult = inst->oneContract ? 1.0 : inst->notional / inst->initialSpot;
        }
    }

    // loop on integration
    virtual void price(const FourierEngine* model,
                       Control*             control,
                       Results*             results){
        static const string method = "VanillaGridFP::price";

        try{
            if (!control){
                return;
            }

            /* Don't want to support vega request (only useful for closed form) */
            control->requestsOutput(OutputRequest::OPTION_VEGA, optionVegaRequest);
            if (optionVegaRequest){
                results->storeRequestResult(optionVegaRequest, IObjectSP(new NotApplicable()));
            }
            control->requestsOutput(OutputRequest::IND_VOL, indVolRequest);
            control->requestsOutput(OutputRequest::OPTION_PRICE, optionPriceRequest);

            if (!(indVolRequest || optionPriceRequest)){
                return; // nothing to do
            }

            FourierProduct::price(model,
                                  control,
                                  results);
        }
        catch (exception& e){
            throw ModelException(e, method);
        }
    }

    /** Constructs default 1D integrands for VanillaGrid */
    virtual Function1DDoubleArrayConstSP Integrand(const FourierEngine* model,
                                                   const Integrator1D*  integrator){
        static const string method = "VanillaGridFP::Integrand";
        try{
            const VanillaGridISAP* isap = dynamic_cast<const VanillaGridISAP*>(&model->getISAP());
            if(!isap) {
                throw ModelException(method, "FourierEngine needs to be provided with a VanillaGridISAP" );
            }
            useOneIntegral = isap->useOneIntegral;

            const DateTimeArray& maturities = *inst->maturities;
            int nbMats = maturities.size();

            fwds = DoubleArraySP(new DoubleArray(nbMats));
            if (inst->fwdStarting) {
                DateTimeArray dates(1 + nbMats);
                dates[0] =  getStartDate();
                int iMat = 0;
                int iDate = 1;
                for (; iMat < nbMats; iMat = iDate++){
                    dates[iDate] = maturities[iMat];
                }
                CDoubleArray myFwds(dates.size());
                mAsset->assetFwdValue(0, dates, myFwds);
                for (iMat = 0, iDate = 1; iMat < nbMats; iMat = iDate++){
                    (*fwds)[iMat] = myFwds[iDate] / myFwds[0];
                }
            }
            else {
                mAsset->assetFwdValue(0, maturities, *fwds);
            }

            const DoubleArray& strikes = *inst->strikes;
            int nbStrikes = strikes.size();
            const DoubleMatrix& instStrikesUsed = *inst->instStrikesUsed;

            int nbIntegrals = isap->useOneIntegral ? 1 : 2;
            nbIntegrals *= nbMats * nbStrikes;
            Function1DDoubleArraySP functions(new Function1DDoubleArray(0));
            functions->reserve(nbIntegrals);

            if (useOneIntegral) {
                if (inst->fwdStarting) {
                    const FwdStFourierProductLogRtn& thisProd = *this;
                    const FwdStFourierProcessLogRtn* thisProc = dynamic_cast<const FwdStFourierProcessLogRtn*>(&model->getProcess());
                    if(!thisProc) {
                        throw ModelException(method, "Process does not support FwdStFourierProcessLogRtn interface");
                    }

                    int iMat = 0;
                    for (; iMat < nbMats; ++iMat){
                        const DateTime& maturity = maturities[iMat];
                        double fwd = (*fwds)[iMat];

                        int iStrike = 0;
                        for(; iStrike < nbStrikes; ++iStrike){
                            bool isCall = (*((*inst->instTypesUsed)[iMat]))[iStrike] || !useOneIntegral;
                            double omega = FourierProduct::intersectRanges(isCall ? callLowerBound : putLowerBound,
                                isCall ? callUpperBound : putUpperBound,
                                thisProc->lowerRealBound(thisProd, maturity),
                                thisProc->upperRealBound(thisProd, maturity),
                                isCall ? isap->payoffToProcFreqBoundWeight : 1.0 - isap->payoffToProcFreqBoundWeight);

                            if (!inst->weights || Maths::isPositive((*inst->weights)[iMat][iStrike])){
                                double strike = instStrikesUsed[iMat][iStrike];
                                functions->push_back(Function1DDoubleSP(new VanillaGridSingleIntegrand<FwdStFourierProcessLogRtn, FwdStFourierProductLogRtn>
                                                                                        (*thisProc,
                                                                                         thisProd,
                                                                                         strike,
                                                                                         fwd,
                                                                                         maturity,
                                                                                         omega)));
                                payoffStrikes.push_back(strike);
                                payoffMaturities.push_back(maturity);
                            }
                        }
                    }
                }
                else {  // started
                    const StFourierProductLogRtn& thisProd = *this;
                    const StFourierProcessLogRtn* thisProc = dynamic_cast<const StFourierProcessLogRtn*>(&model->getProcess());
                    if(!thisProc) {
                        throw ModelException(method, "Process does not support StFourierProcessLogRtn interface");
                    }

                    int iMat = 0;
                    for (; iMat < nbMats; ++iMat){
                        const DateTime& maturity = maturities[iMat];
                        double fwd = (*fwds)[iMat];

                        int iStrike = 0;
                        for(; iStrike < nbStrikes; ++iStrike){
                            bool isCall = (*((*inst->instTypesUsed)[iMat]))[iStrike] || !useOneIntegral;
                            double omega = FourierProduct::intersectRanges(isCall ? callLowerBound : putLowerBound,
                                isCall ? callUpperBound : putUpperBound,
                                thisProc->lowerRealBound(thisProd, maturity),
                                thisProc->upperRealBound(thisProd, maturity),
                                isCall ? isap->payoffToProcFreqBoundWeight : 1.0 - isap->payoffToProcFreqBoundWeight);

                            if (!inst->weights || Maths::isPositive((*inst->weights)[iMat][iStrike])){
                                double strike = instStrikesUsed[iMat][iStrike];
                                functions->push_back(Function1DDoubleSP(new VanillaGridSingleIntegrand<StFourierProcessLogRtn, StFourierProductLogRtn>
                                                                                        (*thisProc,
                                                                                         thisProd,
                                                                                         strike,
                                                                                         fwd,
                                                                                         maturity,
                                                                                         omega)));
                                payoffStrikes.push_back(strike);
                                payoffMaturities.push_back(maturity);
                            }
                        }
                    }
                }
            }
            else {  // !useOneIntegral
                if (inst->fwdStarting) {
                    const FwdStFourierProductLogRtn& thisProd = *this;
                    const FwdStFourierProcessLogRtn* thisProc = dynamic_cast<const FwdStFourierProcessLogRtn*>(&model->getProcess());
                    if(!thisProc) {
                        throw ModelException(method, "Process does not support FwdStFourierProcessLogRtn interface");
                    }

                    int iMat = 0;
                    for (; iMat < nbMats; ++iMat){
                        const DateTime& maturity = maturities[iMat];
                        double fwd = (*fwds)[iMat];

                        int iStrike = 0;
                        for(; iStrike < nbStrikes; ++iStrike){
                            if (!inst->weights || Maths::isPositive((*inst->weights)[iMat][iStrike])){
                                double strike = instStrikesUsed[iMat][iStrike];
                                int j = 0;
                                for (; j < 2; ++j){
                                    functions->push_back(Function1DDoubleSP(new VanillaGridDoubleIntegrand<FwdStFourierProcessLogRtn, FwdStFourierProductLogRtn>
                                                                                            (*thisProc,
                                                                                             thisProd,
                                                                                             strike,
                                                                                             fwd,
                                                                                             maturity,
                                                                                             j)));
                                    payoffStrikes.push_back(strike);
                                    payoffMaturities.push_back(maturity);
                                }
                            }
                        }
                    }
                }
                else {  // started
                    const StFourierProductLogRtn& thisProd = *this;
                    const StFourierProcessLogRtn* thisProc = dynamic_cast<const StFourierProcessLogRtn*>(&model->getProcess());
                    if(!thisProc) {
                        throw ModelException(method, "Process does not support StFourierProcessLogRtn interface.");
                    }

                    int iMat = 0;
                    for (; iMat < nbMats; ++iMat){
                        const DateTime& maturity = maturities[iMat];
                        double fwd = (*fwds)[iMat];

                        int iStrike = 0;
                        for(; iStrike < nbStrikes; ++iStrike){
                            if (!inst->weights || Maths::isPositive((*inst->weights)[iMat][iStrike])){
                                double strike = instStrikesUsed[iMat][iStrike];
                                int j = 0;
                                for (; j < 2; ++j){
                                    functions->push_back(Function1DDoubleSP(new VanillaGridDoubleIntegrand<StFourierProcessLogRtn, StFourierProductLogRtn>
                                                                                            (*thisProc,
                                                                                             thisProd,
                                                                                             strike,
                                                                                             fwd,
                                                                                             maturity,
                                                                                             j)));
                                    payoffStrikes.push_back(strike);
                                    payoffMaturities.push_back(maturity);
                                }
                            }
                        }
                    }
                }
            }
            return functions;
        }
        catch (exception& e){
            throw ModelException(e, method, "Failed to construct integrand(s)");
        }
    }

    /** Post process method for default integrator */
    virtual void postResults(
        const FourierEngine* model,
        const Integrator1D*  integrator,
        const FourierProductIntegrator1D::IntegralArray& integrals,
        CControl*            control,
        CResults*            results) {
        static const string method = "VanillaGridFP::postResults";

        try {
            const DateTimeArray& maturities = *inst->maturities;
            int nbMats = maturities.size();
            const DoubleArray& strikes = *inst->strikes;
            int nbStrikes = strikes.size();
            const DoubleMatrix& instStrikesUsed = *inst->instStrikesUsed;

            VanillaGrid::OutputSP prices;
            VanillaGrid::OutputSP indVols;
            if (optionPriceRequest){
                prices = VanillaGrid::OutputSP(new VanillaGrid::Output(nbMats, nbStrikes));
            }
            if (indVolRequest){
                indVols = VanillaGrid::OutputSP(new VanillaGrid::Output(nbMats, nbStrikes));
            }
            const TimeMetric& timeMetric = model->getProcess().getTimeMetric();

            int iMat = 0; int iIntegral = 0;
            for (; iMat < nbMats; ++iMat){
                const DateTime& maturity = maturities[iMat];
                double pv = inst->instSettle->pv(inst->valueDate,
                                                 maturity,
                                                 inst->discount.get(),
                                                 inst->asset.get());
                double fwd = (*fwds)[iMat];
                double timeFrac = 0.0;
                if (indVolRequest){
                    timeFrac = timeMetric.yearFrac(getStartDate(), maturity);
                }

                int iStrike = 0;
                for (; iStrike < nbStrikes; ++iStrike){
                    if (!inst->weights ||
                        Maths::isPositive((*inst->weights)[iMat][iStrike])){
                        double strike = instStrikesUsed[iMat][iStrike];
                        double price;

                        if (useOneIntegral){
                            price = integrals[iIntegral++] * strike/Maths::PI;
                        }
                        else{
                            double integral0 = integrals[iIntegral++];
                            double integral1 = integrals[iIntegral++];
                            price = (fwd * integral1 - strike * integral0) /
                                Maths::PI;
                            if ((*((*inst->instTypesUsed)[iMat]))[iStrike]) {
                                price += 0.5 * (fwd - strike);
                            }
                            else {  // put
                                price -= 0.5 * (fwd - strike);
                            }
                        }

                        if (indVolRequest){
                            bool isCall = (*(*inst->instTypesUsed)[iMat])[iStrike] ? true : false;

                            double variance;
                            /* call version of impliedVariance that
                               returns status */
                            if (Black::impliedVariance(
                                    isCall,
                                    fwd,
                                    strike,
                                    1.0,       // pv
                                    0.3 * 0.3 * timeFrac,  // initial var guess
                                    price,
                                    2.0 * 0.3 * 1.0e-5 * timeFrac,// var accuracy
                                    variance)){
                                indVols->setValue(iMat, iStrike,
                                                  sqrt(variance / timeFrac));
                            }
                            else {
                                indVols->setError(iMat, iStrike,
                                                  VanillaGrid::Output::ErrorCode::IMP_VOL_ERROR);
                            }
                        }
                        if (optionPriceRequest){
                            prices->setValue(iMat, iStrike,
                                             price * mult * pv);
                        }
                    }
                    else{
                        if (indVolRequest){
                            indVols->setError(iMat, iStrike,
                                              VanillaGrid::Output::ErrorCode::IMP_VOL_NOT_REQUESTED);
                        }
                        if (optionPriceRequest){
                            prices->setError(iMat, iStrike,
                                             VanillaGrid::Output::ErrorCode::PRICE_NOT_REQUESTED);
                        }
                    }
                }
            }

            if (optionPriceRequest){
                results->storeRequestResult(optionPriceRequest, prices);
            }
            if (indVolRequest){
                results->storeRequestResult(indVolRequest, indVols);
            }
        }
        catch (exception& e){
            throw ModelException(e, method);
        }
    }

    /** Gives error reporting when failing to price the i-th payoff */
    virtual string errorHandling(int                  iIntegrand,
                                 const Integrator1D*  integrator) const {
        string errormsg = Format::toString("Failed to price vanilla with strike %f and maturity ",
                                           payoffStrikes[iIntegrand]) +
                          payoffMaturities[iIntegrand].toString() + ".";
        return errormsg;
    }


    /** Constructs FFT 1D integrands for VanillaGrid */
    virtual Function1DComplexArrayConstSP Integrand(const FourierEngine*    model,
                                                    const FFTIntegrator1D*  integrator){
        static const string method = "VanillaGridFP::Integrand";
        try{
            if(!CString::equalsIgnoreCase(inst->instType,WeightMatrix::CALL)
             && !CString::equalsIgnoreCase(inst->instType,WeightMatrix::PUT)) {
                throw ModelException(method, "instType is "+ inst->instType +" while it should be " +
                    WeightMatrix::CALL + " or " + WeightMatrix::PUT);
            }

            const VanillaGridISAP* isap = dynamic_cast<const VanillaGridISAP*>(&model->getISAP());
            if(!isap) {
                throw ModelException(method, "FourierEngine needs to be provided with a VanillaGridISAP." );
            }

            useOneIntegral = isap->useOneIntegral;
            bool isCall = inst->isCall || !useOneIntegral;

            const DateTimeArray& maturities = *inst->maturities;
            int nbMats = maturities.size();

            fwds = DoubleArraySP(new DoubleArray(nbMats));
            if (inst->fwdStarting) {
                DateTimeArray dates(1 + nbMats);
                dates[0] =  getStartDate();
                int iMat = 0;
                int iDate = 1;
                for (; iMat < nbMats; iMat = iDate++){
                    dates[iDate] = maturities[iMat];
                }
                CDoubleArray myFwds(dates.size());
                mAsset->assetFwdValue(0, dates, myFwds);
                for (iMat = 0, iDate = 1; iMat < nbMats; iMat = iDate++){
                    (*fwds)[iMat] = myFwds[iDate] / myFwds[0];
                }
            }
            else {
                mAsset->assetFwdValue(0, maturities, *fwds);
            }

            int nbIntegrals = isap->useOneIntegral ? 1 : 2;
            nbIntegrals *= nbMats;
            Function1DComplexArraySP functions(new Function1DComplexArray(nbIntegrals));

            if (useOneIntegral) {
                omegas = DoubleArraySP(new DoubleArray(nbMats));
                if (inst->fwdStarting) {
                    const FwdStFourierProductLogRtn& thisProd = *this;
                    const FwdStFourierProcessLogRtn* thisProc = dynamic_cast<const FwdStFourierProcessLogRtn*>(&model->getProcess());
                    if(!thisProc) {
                        throw ModelException(method, "Process does not support FwdStFourierProcessLogRtn interface");
                    }

                    int iMat = 0; int iFunc = 0;
                    for (; iMat < nbMats; ++iMat){
                        const DateTime& maturity = maturities[iMat];
                        (*omegas)[iMat] = FourierProduct::intersectRanges(isCall ? callLowerBound : putLowerBound,
                                                                          isCall ? callUpperBound : putUpperBound,
                                                                          thisProc->lowerRealBound(thisProd, maturity),
                                                                          thisProc->upperRealBound(thisProd, maturity),
                                                                          isCall ? isap->payoffToProcFreqBoundWeight : 1.0 - isap->payoffToProcFreqBoundWeight);

                        (*functions)[iFunc++] = Function1DComplexSP(new VanillaGridSingleIntegrandFFT<FwdStFourierProcessLogRtn, FwdStFourierProductLogRtn>
                                                                                 (*thisProc,
                                                                                 thisProd,
                                                                                 maturity,
                                                                                 (*omegas)[iMat]));
                         payoffMaturities.push_back(maturity);
                    }
                }
                else {  // started
                    const StFourierProductLogRtn& thisProd = *this;
                    const StFourierProcessLogRtn* thisProc = dynamic_cast<const StFourierProcessLogRtn*>(&model->getProcess());
                    if(!thisProc) {
                        throw ModelException(method, "Process does not support StFourierProcessLogRtn interface");
                    }

                    int iMat = 0; int iFunc = 0;
                    for (; iMat < nbMats; ++iMat){
                        const DateTime& maturity = maturities[iMat];
                        (*omegas)[iMat] = FourierProduct::intersectRanges(isCall ? callLowerBound : putLowerBound,
                                                                          isCall ? callUpperBound : putUpperBound,
                                                                          thisProc->lowerRealBound(thisProd, maturity),
                                                                          thisProc->upperRealBound(thisProd, maturity),
                                                                          isCall ? isap->payoffToProcFreqBoundWeight : 1.0 - isap->payoffToProcFreqBoundWeight);

                        (*functions)[iFunc++] = Function1DComplexSP(new VanillaGridSingleIntegrandFFT<StFourierProcessLogRtn, StFourierProductLogRtn>
                                                                                 (*thisProc,
                                                                                 thisProd,
                                                                                 maturity,
                                                                                 (*omegas)[iMat]));
                        payoffMaturities.push_back(maturity);
                    }
                }
            }
            else {  // !useOneIntegral
                if (inst->fwdStarting) {
                    const FwdStFourierProductLogRtn& thisProd = *this;
                    const FwdStFourierProcessLogRtn* thisProc = dynamic_cast<const FwdStFourierProcessLogRtn*>(&model->getProcess());
                    if(!thisProc) {
                        throw ModelException(method, "Process does not support FwdStFourierProcessLogRtn interface");
                    }

                    int iMat = 0; int iFunc = 0;
                    for (; iMat < nbMats; ++iMat){
                        const DateTime& maturity = maturities[iMat];

                        int j = 0;
                        for (; j < 2; ++j){
                            (*functions)[iFunc++] = Function1DComplexSP(new VanillaGridDoubleIntegrandFFT<FwdStFourierProcessLogRtn, FwdStFourierProductLogRtn>
                                                                                     (*thisProc,
                                                                                     thisProd,
                                                                                     maturity,
                                                                                     j));
                            payoffMaturities.push_back(maturity);
                        }
                    }
                }
                else {  // started
                    const StFourierProductLogRtn& thisProd = *this;
                    const StFourierProcessLogRtn* thisProc = dynamic_cast<const StFourierProcessLogRtn*>(&model->getProcess());
                    if(!thisProc) {
                        throw ModelException(method, "Process does not support StFourierProcessLogRtn interface.");
                    }

                    int iMat = 0; int iFunc = 0;
                    for (; iMat < nbMats; ++iMat){
                        const DateTime& maturity = maturities[iMat];

                        int j = 0;
                        for (; j < 2; ++j){
                            (*functions)[iFunc++] = Function1DComplexSP(new VanillaGridDoubleIntegrandFFT<StFourierProcessLogRtn, StFourierProductLogRtn>
                                                                                    (*thisProc,
                                                                                     thisProd,
                                                                                     maturity,
                                                                                     j));
                            payoffMaturities.push_back(maturity);
                        }
                    }
                }
            }
            return functions;
        }
        catch (exception& e){
            throw ModelException(e, method, "Failed to construct integrand(s)");
        }
    }

    /** Post process method for FFT integrator */
    virtual void postResults(
        const FourierEngine* model,
        const FFTIntegrator1D*  integrator,
        const FourierProductFFTIntegrator1D::IntegralArray& integrals,
        CControl*            control,
        CResults*            results) {
        static const string method = "VanillaGridFP::postResults";

        try {
            const DateTimeArray& maturities = *inst->maturities;
            int nbMats = maturities.size();
            const DoubleArray& strikes = *inst->strikes;
            int nbStrikes = strikes.size();
            const DoubleMatrix& instStrikesUsed = *inst->instStrikesUsed;

            VanillaGrid::OutputSP prices;
            VanillaGrid::OutputSP indVols;
            if (optionPriceRequest){
                prices = VanillaGrid::OutputSP(new VanillaGrid::Output(nbMats, nbStrikes));
            }
            if (indVolRequest){
                indVols = VanillaGrid::OutputSP(new VanillaGrid::Output(nbMats, nbStrikes));
            }
            const TimeMetric& timeMetric = model->getProcess().getTimeMetric();

            int iMat = 0; int iIntegral = 0;
            for (; iMat < nbMats; ++iMat){
                const DateTime& maturity = maturities[iMat];
                double pv = inst->instSettle->pv(inst->valueDate,
                                                 maturity,
                                                 inst->discount.get(),
                                                 inst->asset.get());
                double fwd = (*fwds)[iMat];
                double timeFrac = 0.0;
                if (indVolRequest){
                    timeFrac = timeMetric.yearFrac(getStartDate(), maturity);
                }

                double omega = 0.0;
                FourierProductFFTIntegrator1D::Integral integral0;
                FourierProductFFTIntegrator1D::Integral integral1;
                if (useOneIntegral){
                    omega = (*omegas)[iMat];
                    integral0 = integrals[iIntegral++];
                }
                else{
                    integral0 = integrals[iIntegral++];
                    integral1 = integrals[iIntegral++];
                }

                int iStrike = 0;
                for (; iStrike < nbStrikes; ++iStrike){
                    double strike = instStrikesUsed[iMat][iStrike];
                    double logMoneyness = log(strike / fwd);

                    double price;
                    if (useOneIntegral){
                        price = exp(-omega * logMoneyness) *
                            integral0->getValue(logMoneyness) * strike;
                    } else{
                        price = (fwd * integral1->getValue(logMoneyness) -
                                 strike * integral0->getValue(logMoneyness)) /
                            Maths::PI;
                        if ((*((*inst->instTypesUsed)[iMat]))[iStrike]) {
                            price += 0.5 * (fwd - strike);
                        } else {  // put
                            price -= 0.5 * (fwd - strike);
                        }
                    }

                    if (indVolRequest){
                        bool isCall = (*(*inst->instTypesUsed)[iMat])[iStrike] ? true : false;
                        double variance;
                        // call version of impliedVariance that returns status
                        if (Black::impliedVariance(
                                isCall,
                                fwd,
                                strike,
                                1.0,       // pv
                                0.3 * 0.3 * timeFrac,  // initial var guess
                                price,
                                2.0 * 0.3 * 1.0e-5 * timeFrac,
                                variance)){
                            indVols->setValue(iMat, iStrike,
                                              sqrt(variance / timeFrac));
                        } else {
                            indVols->setError(iMat, iStrike,
                                             VanillaGrid::Output::ErrorCode::IMP_VOL_ERROR);
                        }
                    }
                    if (optionPriceRequest){
                        prices->setValue(iMat, iStrike,
                                         price * mult * pv);
                    }
                }
            }

            if (optionPriceRequest){
                results->storeRequestResult(optionPriceRequest, prices);
            }
            if (indVolRequest){
                results->storeRequestResult(indVolRequest, indVols);
            }
        }
        catch (exception& e){
            throw ModelException(e, method);
        }
    }

    const DateTime& getStartDate() const {
        return (inst->fwdStarting ? inst->startDate: inst->valueDate);
    }

    /** Gives error reporting when failing to price the i-th payoff */
    virtual string errorHandling(int                    iIntegrand,
                                 const FFTIntegrator1D* integrator) const {
        string errormsg = Format::toString("Failed to price vanillas for maturity ") +
                          payoffMaturities[iIntegrand].toString() + ".";
        return errormsg;
    }

private:
    CDoubleArray  payoffStrikes;
    DateTimeArray payoffMaturities;

    const VanillaGrid* inst;
    double             mult;
    //bool               isCall;
    bool               useOneIntegral;
    CDoubleArraySP     fwds;
    CDoubleArraySP     omegas;
    OutputRequest*     indVolRequest;
    OutputRequest*     optionPriceRequest;
    OutputRequest*     optionVegaRequest;

    static const double callLowerBound;
    static const double callUpperBound;
    static const double putLowerBound;
    static const double putUpperBound;
};

const double VanillaGridFP::callLowerBound = 1.0;
const double VanillaGridFP::callUpperBound = 100.0;     // infinity
const double VanillaGridFP::putLowerBound = -100.0;     // infinity
const double VanillaGridFP::putUpperBound = 0.0;

/** Implementation of FourierEngine::IntoProduct interface */
FourierProduct* VanillaGrid::createProduct(const FourierEngine* model) const {
    return new VanillaGridFP(this);
}

#if 0
/** make a simple started vanilla grid ready for pricing */
VanillaGrid* VanillaGrid::make(const DateTime&             valueDate,
                         bool                        isCall,
                         bool                        american,
                         const Schedule*             exerciseSchedule,
                         const CAsset*               asset,
                         const YieldCurve*           discount,
                         const InstrumentSettlement* settle,
                         int                         noExerciseWindow) {
    static const string routine = "VanillaGrid::make";
    try {
        VanillaGridSP vanilla(new VanillaGrid());

        vanilla->valueDate        = valueDate;
        vanilla->isCall           = isCall;
        vanilla->canExerciseEarly = american;
        vanilla->exerciseSchedule = ScheduleSP(copy(exerciseSchedule));

        vanilla->oneContract      = true;
        vanilla->notional         = 1.0;
        vanilla->initialSpot      = 1.0;
        vanilla->asset            = CAssetWrapper(copy(asset));
        vanilla->ccyTreatment     = CAsset::CCY_TREATMENT_NONE;
        vanilla->discount         = YieldCurveWrapper(copy(discount));
        vanilla->instSettle       = InstrumentSettlementSP(copy(settle));
        vanilla->noExerciseWindow = noExerciseWindow;

        return vanilla.release();
    }
    catch (exception& e){
        throw ModelException(e, routine);
    }
}
#endif

InstrumentSettlementSP VanillaGrid::getDefaultSettlement() {
    return (InstrumentSettlementSP(new CashSettlePeriod(0)));
}

VanillaGrid::VanillaGrid(const CAssetWrapper&          asset,
                         const YieldCurveWrapper&      discount,
                         const WeightMatrixWrapper&     weightMatrix,
                         const InstrumentSettlementSP& instSettle,
                         const string&                 ccyTreatment):
CInstrument(TYPE),
fwdStarting(false),
oneContract(true),
notional(0.0),  // not used
initialSpot(0.0),   // not used
asset(asset),
ccyTreatment(ccyTreatment),
discount(discount),
weightMatrix(weightMatrix),
instSettle(instSettle) {
    validatePop2Object();
}

VanillaGrid::VanillaGrid(const DateTime&               valueDate,
                         const CAssetWrapper&          asset,
                         const YieldCurveWrapper&      discount,
                         const WeightMatrixWrapper&    weightMatrix,
                         const InstrumentSettlementSP& instSettle):
CInstrument(TYPE),
startDate(valueDate),
fwdStarting(false),
oneContract(true),
notional(0.0),  // not used
initialSpot(0.0),   // not used
valueDate(valueDate),
asset(asset),
ccyTreatment(CAsset::CCY_TREATMENT_NONE),
discount(discount),
weightMatrix(weightMatrix),
instSettle(instSettle) {
    validatePop2Object();
}

VanillaGrid::VanillaGrid(bool                          fwdStarting,
                         DateTime                      startDate,
                         bool                          oneContract,
                         double                        notional,
                         double                        initialSpot,
                         const CAssetWrapper&          asset,
                         const YieldCurveWrapper&      discount,
                         const WeightMatrixWrapper&    weightMatrix,
                         const InstrumentSettlementSP& instSettle):
CInstrument(TYPE),
startDate(startDate),
fwdStarting(fwdStarting),
oneContract(oneContract),
notional(notional),
initialSpot(initialSpot),
asset(asset),
ccyTreatment(CAsset::CCY_TREATMENT_NONE),
discount(discount),
weightMatrix(weightMatrix),
instSettle(instSettle) {
    validatePop2Object();
}

VanillaGrid::VanillaGrid(bool                          fwdStarting,
                         DateTime                      startDate,
                         bool                          oneContract,
                         double                        notional,
                         double                        initialSpot,
                         const DateTime&               valueDate,
                         const CAssetWrapper&          asset,
                         const string&                 ccyTreatment,
                         const YieldCurveWrapper&      discount,
                         const WeightMatrixWrapper&    weightMatrix,
                         const InstrumentSettlementSP& instSettle,
                         const InstrumentSettlementSP& premiumSettle):
CInstrument(TYPE),
startDate(startDate),
fwdStarting(fwdStarting),
oneContract(oneContract),
notional(notional),
initialSpot(initialSpot),
valueDate(valueDate),
asset(asset),
ccyTreatment(ccyTreatment),
discount(discount),
weightMatrix(weightMatrix),
instSettle(instSettle),
premiumSettle(premiumSettle) {
    validatePop2Object();
}

class VanillaGridHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VanillaGrid, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(CClosedFormLN::IIntoProduct);
        IMPLEMENTS(IMCIntoProduct);
        IMPLEMENTS(FourierEngine::IIntoProduct);
        IMPLEMENTS(FDModel::IIntoProduct);
//        IMPLEMENTS(ISensitiveStrikes);
//        IMPLEMENTS(Theta::Shift);
//        IMPLEMENTS(LastSensDate);
        EMPTY_SHELL_METHOD(defaultVanillaGrid);
        FIELD(startDate, "Option start date");
        FIELD_MAKE_OPTIONAL(startDate);
        FIELD(isCall, "Is it a call option");
        FIELD(fwdStarting, "Is it a fwd starting option");
        FIELD(maturities, "Array of maturities - ignored if weightMatrix passed");
        FIELD(strikes, "Array of strikes - ignored if weightMatrix passed");
        FIELD(weights, "Matrix of weights - ignored if weightMatrix passed");
        FIELD_MAKE_OPTIONAL(maturities);
        FIELD_MAKE_OPTIONAL(strikes);
        FIELD_MAKE_OPTIONAL(weights);
        FIELD(oneContract, "Calc price for 1 contract");
        FIELD(notional, "Option notional");
        FIELD_MAKE_OPTIONAL(notional);
        FIELD(initialSpot,      "Initial spot price");
        FIELD_MAKE_OPTIONAL(initialSpot);

        FIELD(valueDate, "Value date");
        FIELD_MAKE_OPTIONAL(valueDate);
        FIELD(asset, "Asset");
        FIELD(discount, "Discount curve");
        FIELD(ccyTreatment, "Currency treatment");
        FIELD(instSettle, "Instrument settlement at maturity");
        FIELD(premiumSettle, "Settlement of option premium");
        FIELD_MAKE_OPTIONAL(premiumSettle);
        FIELD(weightMatrix, "weightMatrix");
        FIELD_MAKE_OPTIONAL(weightMatrix);
        FIELD(strikeUnits, "strikeUnits");
        FIELD_MAKE_TRANSIENT(strikeUnits);
        FIELD(instType, "instType");
        FIELD_MAKE_TRANSIENT(instType);
        FIELD(instStrikesUsed, "instStrikesUsed");
        FIELD_MAKE_TRANSIENT(instStrikesUsed);
        FIELD(instTypesUsed, "instTypesUsed");
        FIELD_MAKE_TRANSIENT(instTypesUsed);
    }

    static IObject* defaultVanillaGrid(){
        return new VanillaGrid();
    }
};

CClassConstSP const VanillaGrid::TYPE = CClass::registerClassLoadMethod(
    "VanillaGrid", typeid(VanillaGrid), VanillaGridHelper::load);

// LEAST SQUARE FIT
typedef VanillaGrid::LeastSquareFit::FitType VGLSFitType;
template<> string nameForType<VGLSFitType>(VGLSFitType*){
    return "VanillaGrid::LeastSquareFit::FitType";
}

typedef Enum2StringListHelper<VGLSFitType> VGLSFitTypeHelper;
template<> string VGLSFitTypeHelper::names[VGLSFitTypeHelper::EnumList::NB_ENUMS] = {
    "PRICE",
    "VEGA_NORMALIZED_PRICE",
    "IMPLIED_VOL"
};

void VanillaGrid::LeastSquareFit::validatePop2Object(){
    static const string method("VanillaGrid::LeastSquareFit::validatePop2Object");
    try{
        fitWhat = VGLSFitTypeHelper::getIndex(fitType);

        // what request
        switch(fitWhat){
        case VGLSFitTypeHelper::EnumList::IMPLIED_VOL:
            outReqName = OutputRequest::IND_VOL;
            break;
        default:
            outReqName = OutputRequest::OPTION_PRICE;
        }

        vanillaGridOriginal = VanillaGridSP(copy(vanillaGrid.get()));

        /* Create control */
        SensitivityArraySP sens(new SensitivityArray(0));
        OutputRequestArraySP outReqs(new OutputRequestArray(1));
        (*outReqs)[0] = OutputRequestSP(new OutputRequest(outReqName));
        control = CControlSP(new Control(sens,
                                         outReqs,
                                         false,
                                         ""));
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

void VanillaGrid::LeastSquareFit::getMarket(MarketData* market){
    static const string method = "VanillaGrid::LeastSquareFit::getMarket";

    try{
        /* Get market for closed form model */
        IModelSP mktmodel(new CClosedFormLN(volType));
        SensitivityArraySP sens(new SensitivityArray(0));
        OutputRequestArraySP outReqs(new OutputRequestArray(1));
        (*outReqs)[0] = OutputRequestSP(new OutputRequest(outReqName));
        if (fitWhat == VGLSFitTypeHelper::EnumList::VEGA_NORMALIZED_PRICE){
            outReqs->push_back(OutputRequestSP(new OutputRequest(OutputRequest::OPTION_VEGA)));
        }
        CControlSP ctrl(new Control(sens,
                                      outReqs,
                                      false,
                                      ""));
        Calibrator::ObjFunc::Helper::getMarket(vanillaGrid.get(),
                                               mktmodel.get(),
                                               ctrl.get(),
                                               market);

        /* Compute mkt vols or prices + vegas */
        CResultsSP results(mktmodel->Run(vanillaGrid.get(), ctrl.get()));
        IObjectConstSP mktvalsobj(results->retrieveRequestResult(outReqName));
        mktVals = VanillaGrid::OutputSP::constCast(VanillaGrid::OutputConstSP::dynamicCast(mktvalsobj));
        if (fitWhat == VGLSFitTypeHelper::EnumList::VEGA_NORMALIZED_PRICE){
            IObjectConstSP mktvegasobj(results->retrieveRequestResult(OutputRequest::OPTION_VEGA));
            mktVegas = VanillaGrid::OutputSP::constCast(VanillaGrid::OutputConstSP::dynamicCast(mktvegasobj));
        }

        /* Get market for inputted model */
        Calibrator::ObjFunc::Helper::getMarket(vanillaGrid.get(),
                                               model.get(),
                                               control.get(),
                                               market);
        Calibrator::ObjFunc::Helper::getMarket(vanillaGridOriginal.get(),
                                               model.get(),
                                               control.get(),
                                               market);
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}


void VanillaGrid::LeastSquareFit::validate(){
    static const string method = "VanillaGrid::LeastSquareFit::validate";
    try{
        /* Return an error message if using Monte-Carlo engine and settlement type 4 ("at maturity")
           indeed the implied volatility are badly calculated when using stochastic discount factors.
           This validation should be updated once it will be possible to know if the model is using
           stochastic discount factors or not */
        if (MonteCarlo::TYPE->isInstance(model) && AtMaturity::TYPE->isInstance(vanillaGrid->instSettle))
        {
            throw ModelException(method, "the settlement type can't be At Maturity "
                                 "when using Monte-Carlo engine");
        }

        /* If no weights provided, create default ones (a la Denis Mustafa's !)*/
        if (!vanillaGrid->weights){
            double spot = vanillaGrid->asset.get()->getSpot();

            
            CDoubleMatrixSP weightsSP = DefaultWeightMatrix::create(vanillaGrid->valueDate,
                                                                    spot,
                                                                    *vanillaGrid->maturities,
                                                                    *vanillaGrid->strikes,
                                                                    vanillaGrid->weightMatrix->getStrikeUnits(),
                                                                    vanillaGrid->weightMatrix->getInstType(),
                                                                    DefaultWeightMatrix::OLDMODE);
            vanillaGrid->weightMatrix->setWeights(weightsSP);
            vanillaGrid->weights = weightsSP;

            vanillaGrid->weights->checkNonNegative();
        }

        /* Calculate total nb of funcs (that's the number of non-zero weights)
           and compute square root of weights */
        CDoubleMatrix& weights = *vanillaGrid->weights;
        int nbRows = weights.numRows();
        int nbCols = weights.numCols();
        nbFuncs = 0;

        // to normalize weights
        DoubleArray sumOfWeights(nbCols);
        int idxCol = 0,
            idxRow = 0;
        for(idxCol = 0 ; idxCol < nbCols ; idxCol++){
            double sum = 0.;
            for(idxRow = 0 ; idxRow < nbRows ; idxRow++){
                double weight = weights[idxCol][idxRow];
                sum += weight;
            }
            sumOfWeights[idxCol] = sum;
        }

        CDoubleMatrix sumOfWeightsMatrix(nbCols, nbRows);
        for(idxCol = 0 ; idxCol < nbCols ; idxCol++){
            for(idxRow = 0 ; idxRow < nbRows ; idxRow++){
                sumOfWeightsMatrix[idxCol][idxRow] = sumOfWeights[idxCol];
            }
        }


        int iRow = 0;
        for (; iRow < nbRows; ++iRow){
            int iCol = 0;
            for (; iCol < nbCols; ++iCol){
                double weight = weights[iCol][iRow];
                double _sumOfWeights = sumOfWeightsMatrix[iCol][iRow];
                if (Maths::isPositive(weight)){
                    if (CString::equalsIgnoreCase(normalizeWeights,VanillaGrid::LeastSquareFit::NORMALIZE)){
                        // normalize weights so that they sum to 1
                        if(Maths::isPositive(_sumOfWeights)){
                            weights[iCol][iRow] = sqrt(weight/_sumOfWeights);
                         }
                         if(!Maths::isPositive(_sumOfWeights)){
                             weights[iCol][iRow] = sqrt(weight);
                          }
                    }
                    else{
                        weights[iCol][iRow] = sqrt(weight);
                    }
                    ++nbFuncs;
                    // if vega weighted, must ensure we don't divide by zero vega
                    if (fitWhat == VanillaGrid::LeastSquareFit::FitType::VEGA_NORMALIZED_PRICE){
                        if (!Maths::isPositive(mktVegas->getValue(iCol, iRow))){
                            throw ModelException(method,
                                                 "the weight ("
                                                 + Format::toString(weight)
                                                 + ") of the option with strike "
                                                 + Format::toString((*vanillaGrid->instStrikesUsed)[iCol][iRow])
                                                 + "\nand maturity "
                                                 + (*vanillaGrid->maturities)[iCol].toString()
                                                 + " is non-zero, yet its vega is zero");
                        }
                        else{
                            weights[iCol][iRow] /= mktVegas->getValue(iCol, iRow);
                        }
                    }
                }
            }
        }
        vanillaGridOriginal = VanillaGridSP(copy(vanillaGrid.get()));
        if (nbFuncs == 0){
            throw ModelException(method,
                                 "The weight matrix is identically zero");
        }
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

void VanillaGrid::LeastSquareFit::useNormalizedWeights(string normalize){
    normalizeWeights = normalize;
}


const string VanillaGrid::LeastSquareFit::NORMALIZE  = "normalize"; 

IObjectSP VanillaGrid::LeastSquareFit::getAdjustableGroup(){
    return IObjectSP::attachToRef(this);
}

int VanillaGrid::LeastSquareFit::getNbFuncs() const{
    return nbFuncs;
}

void VanillaGrid::LeastSquareFit::calcValue(CDoubleArray& funcvals) const{
    try{
        CResultsSP results(model->Run(vanillaGrid.get(), control.get()));
        IObjectConstSP resobj(results->retrieveRequestResult(outReqName));
        VanillaGrid::OutputConstSP res(VanillaGrid::OutputConstSP::dynamicCast(resobj));
        const DoubleMatrix& weights = *vanillaGrid->weights;
        int nbRows = res->numRows();
        int nbCols = res->numCols();
        int iFunc = 0;
        int iRow = 0;
        for (; iRow < nbRows; ++iRow){
            int iCol = 0;
            for (; iCol < nbCols; ++iCol){
                double weight = weights[iCol][iRow];
                if (Maths::isPositive(weight)){
                    funcvals[iFunc++] = weight * (res->getValue(iCol, iRow) - mktVals->getValue(iCol, iRow));
                }
            }
        }
    }
    catch(exception& e){
        throw ModelException(e, "VanillaGrid::LeastSquareFit::calcValue");
    }
}


void VanillaGrid::LeastSquareFit::getSensitivities(CDoubleMatrix& sensitivities,
                                                   CDoubleMatrix& valGuess,
                                                   DoubleArray shift,
                                                   int index) const{
    try{
        CResultsSP results(model->Run(vanillaGrid.get(), control.get()));
        IObjectConstSP resobj(results->retrieveRequestResult(outReqName));
        VanillaGrid::OutputConstSP res(VanillaGrid::OutputConstSP::dynamicCast(resobj));
        const DoubleMatrix& weights = *vanillaGrid->weights;
        int nbCols = res->numCols();
        int nbRows = res->numRows();
        int idx= 0;
        int iCol = 0;

        for (; iCol < nbCols; ++iCol){
            int iRow = 0;
            double weightCenter = weights[iCol][index];
            for (; iRow < nbRows; ++iRow){
                double weight = weights[iCol][iRow];
                if (Maths::isPositive(weight)){
                        double shiftVal = weight * (res->getValue(iCol,iRow) - mktVals->getValue(iCol,iRow));
                        sensitivities[iCol-idx][iRow] = (shiftVal - valGuess[iCol-idx][iRow])/shift[iCol-idx];
                    }
            }
            if(!Maths::isPositive(weightCenter)){
                ++ idx;
            }
        }

      
    }catch(exception& e){
        throw ModelException(e,"VanillaGrid::LeastSquareFit::getSensitivities");
    }
}

void VanillaGrid::LeastSquareFit::calcInitialVals(CDoubleMatrix& valGuess,
                                                  int index) const{
    try{
        CResultsSP results(model->Run(vanillaGrid.get(), control.get()));
        IObjectConstSP resobj(results->retrieveRequestResult(outReqName));
        VanillaGrid::OutputConstSP res(VanillaGrid::OutputConstSP::dynamicCast(resobj));
        const DoubleMatrix& weights = *vanillaGrid->weights;
        int nbCols = res->numCols();
        int nbRows = res->numRows();
        int idx= 0;
        int iCol = 0;
        
        for (; iCol < nbCols; ++iCol){
            int iRow = 0;
            double weightCenter = weights[iCol][index];
            for (; iRow < nbRows; ++iRow){
                double weight = weights[iCol][iRow];
                if (Maths::isPositive(weight)){
                    valGuess[iCol-idx][iRow] = weight * (res->getValue(iCol,iRow) - mktVals->getValue(iCol,iRow));
                }
            }
            if(!Maths::isPositive(weightCenter)){
                ++ idx;
            }
        }

    }catch(exception& e){
        throw ModelException(e,"VanillaGrid::LeastSquareFit::calcInitialVals");
    }
}

/** Makes additional validations for calibration with bootstrapping */
void VanillaGrid::LeastSquareFit::validate(const Calibrator::InstanceID::IBootstrappableArray& ids) const{
    static const string method("VanillaGrid::LeastSquareFit::validate(const InstanceID::IBootstrappableArray&)");

    int nberMat = (*vanillaGridOriginal->maturities).size(); // nber of maturities in VanillaGrid
    int nberStrk = vanillaGridOriginal->strikes->size(); // nber of strikes in VanillaGrid

    // Weak since there seems to be an assumption that startDate is referenced only if !fwdStarting
    // but the code to support that assumption is spread everywhere (though originally it wasn't done here...)
    const DateTime& startDate = vanillaGridOriginal->fwdStarting? vanillaGridOriginal->startDate :
        vanillaGridOriginal->valueDate;

    // check that the weight matrix has at least a line with non-zero weights
    BoolArray NonZeroLine(nberMat, false);
    bool hasNonZeroLine = false;

    int iMat;
    for (iMat = 0 ; iMat < nberMat ; iMat++)
    {
        const double* weights = (*vanillaGridOriginal->weights)[iMat]; // weights for the matutity iMat
        // loop on the strikes
        int iStrk;
        for (iStrk = 0 ; iStrk < nberStrk ; iStrk++)
        {
            if (Maths::isPositive(weights[iStrk]))
            {
                NonZeroLine[iMat] = true; // the line has a non-zero weight
                if (!hasNonZeroLine)
                {
                    hasNonZeroLine = true;
                }
                break;
            }
        }
    }

    if (!hasNonZeroLine){
        throw ModelException(method, "there are no non-zero weights in the weights matrix.");
    }

    // check that the weight matrix doesn't have a line with only zero weights
    // separating two blocks with non-zero weights
    // loop on the maturities
    int FirstTrue = nberMat - 1; // first line with non zero weights
    int LastTrue = 0; // last line with non zero weights

    // first line with non zero weights or nberMat - 1 if only zero lines
    for (iMat = 0 ; iMat < nberMat ; iMat++)
    {
        if (NonZeroLine[iMat])
        {
            FirstTrue = iMat;
            break;
        }
    }

    // last line with non zero weights or 0 if only zero lines
    for (iMat = nberMat - 1 ; iMat >= 0 ; iMat--)
    {
        if (NonZeroLine[iMat])
        {
            LastTrue = iMat;
            break;
        }
    }

    if (nberMat > 2 && (FirstTrue + 1 < LastTrue))
        // there could be problems only in the case where the number of maturities is stricly greater than 2
        // there are at least two non-zero lines which are not one just after the other
    {
        for (iMat = FirstTrue + 1 ; iMat < LastTrue ; iMat++)
        {
            if (!NonZeroLine[iMat])
                throw ModelException(method, "there is a zero line in the weights matrix "
                                     "separating two blocks with non zero weights.");
        }
    }

    // validate against vanilla grid dates....
    int nbIds = ids.size();
    if (!nbIds){
        throw ModelException(method, "the ids were not found.");
    }

    IObjectConstSP me(IObjectConstSP::attachToRef(this));
    ExpiryArraySP expiries = ids[0]->getExpiries(me);
    int iIds;
    for (iIds = 0 ; iIds < nbIds ; ++iIds)
    {
        ExpiryArraySP newExpiries = ids[iIds]->getExpiries(me);
        if (!Expiry::equals(startDate, newExpiries.get(), expiries.get())){
            throw ModelException(method, "all the ids don't have the same expiries.");
        }
    }

    // check that the expiries correspond to the maturities with non-zero weights
    int nbExpiries = expiries->size();
    int nberNonZeroMats = LastTrue - FirstTrue + 1; // number of maturity dates with non-zero weights

    if (nbExpiries != nberNonZeroMats){
        throw ModelException(method, "the number of expiries (" + Format::toString(nbExpiries) +
                             ") is different from the number of maturities with non zero weights (" +
                             Format::toString(nberNonZeroMats) + ")");
    }

    int iExpiry;
    for (iExpiry = 0 ; iExpiry < nbExpiries ; iExpiry++){
        DateTime date = (*expiries)[iExpiry]->toDate(startDate);
        if (date != (*vanillaGridOriginal->maturities)[iExpiry + FirstTrue]){
            throw ModelException(method, "the expiries are different from the maturities:"
                                 "Expiry #" + Format::toString(iExpiry+1) + " is " +
                                 (*expiries)[iExpiry]->toString() + " giving a date of " +
                                 date.toString() + " compared to a maturity of " +
                                 (*vanillaGridOriginal->maturities)[iExpiry + FirstTrue].toString());
        }
    }
}

/** Gets the first maturity with non zero weight */
const int VanillaGrid::LeastSquareFit::getIdxFirstMat(int nberMat, int nberStrk) const{
    static const string method("VanillaGrid::LeastSquareFit::getIdxFirstMat");

    int idxFirstMat = 0;
    bool foundNonZero = false;

    // loop on the maturities
    int iMat;
    for (iMat = 0 ; iMat < nberMat ; iMat++)
    {
        const double* weights = (*vanillaGridOriginal->weights)[iMat]; // weights for the matutity iMat
        // loop on the strikes
        int iStrk;
        for (iStrk = 0 ; iStrk < nberStrk ; iStrk++)
        {
            if (! foundNonZero && Maths::isPositive(weights[iStrk]))
            {
                idxFirstMat = iMat;
                foundNonZero = true;
            }
        }
    }

    // all the weights are equal to 0.0
    if (!foundNonZero){
        throw ModelException(method, "the weights in the weights matrix are all <= 0.0");
    }

    return idxFirstMat;

}

/** Updates the instrument for calibration with bootstrapping
    before each run of the calibrator */
void VanillaGrid::LeastSquareFit::update(int idxSmileMat) {

    int nberMat = vanillaGridOriginal->maturities->size(); // nber of maturities in VanillaGrid
    int nberStrk = vanillaGridOriginal->strikes->size(); // nber of strikes in VanillaGrid

    // the index idxSmileMat is related to the equity smile parameters
    // the correponding maturity has an index idxMat which can be different
    // the assumption is that the matrix weight doesn't have a line with only zero weights
    // separating 2 blocks with non-zero weights
    int idxFirstMat = getIdxFirstMat(nberMat, nberStrk);

    // retrieve the whole VanillaGrid
    vanillaGrid = VanillaGridSP(copy(vanillaGridOriginal.get()));

    // if useSameRandonNumber flag is set to true
    // we only need to set all the weights to 0.0
    // except on the maturity we are calibrating on
    if (useSameRandomNumbers)
    {
        // loop on the maturities
        int iMat;
        for (iMat = 0 ; iMat < nberMat ; iMat++)
        {
            if (iMat != (idxFirstMat + idxSmileMat))
            {
                // loop on the strikes
                int iStrk;
                for (iStrk = 0 ; iStrk < nberStrk ; iStrk++)
                {
                    (*vanillaGrid->weights)[iMat][iStrk] = 0.0;
                }
            }
        }
    }
    // if useSameRandonNumber flag is set to true
    // we only need to set all the weights to 0.0
    // except on the maturity we are calibrating on
    else
    {
        if (idxFirstMat + idxSmileMat + 1 < nberMat -1)
        {
            // erase maturities after maturity we are calibrating on (for field maturities)
            vanillaGrid->maturities->erase(vanillaGrid->maturities->begin() + idxFirstMat + idxSmileMat + 1,
                                           vanillaGrid->maturities->begin() + nberMat - 1);

            // erase maturities after maturity we are calibrating on (for field weights)
            (*vanillaGrid->weights).transpose();

            int iMat;
            for (iMat = 0 ; iMat < idxFirstMat + idxSmileMat ; iMat++)
            {
                int iStrk;
                for (iStrk = 0 ; iStrk < nberStrk ; iStrk++)
                {
                    (*vanillaGrid->weights)[iMat][iStrk] = 0.0;
                }
            }

            for (iMat = idxFirstMat + idxSmileMat + 1 ; iMat < nberMat ; iMat++)
            {
                (*vanillaGrid->weights).removeLastCol();
            }

            (*vanillaGrid->weights).transpose();
        }
    }
}


/** Reset the instrument for calibration with bootstrapping
    after each run of the calibrator */
void VanillaGrid::LeastSquareFit::reset() {

    // retrieve the whole VanillaGrid
    vanillaGrid = VanillaGridSP(copy(vanillaGridOriginal.get()));
}


// for reflection
VanillaGrid::LeastSquareFit::LeastSquareFit():
Calibrator::ObjFuncLeastSquare(TYPE),
useSameRandomNumbers(true),
fitType(VGLSFitTypeHelper::getDefaultName()),
nbFuncs(0), volType(VanillaGrid::LeastSquareFit::VOLSURFACE){}

VanillaGrid::LeastSquareFit::LeastSquareFit(const IModel&       model,
                                            const VanillaGrid& vanillaGrid,
                                            const string&      fitType,
                                            const string&      volType):
Calibrator::ObjFuncLeastSquare(TYPE),
model(copy(&model)),
vanillaGrid(copy(&vanillaGrid)),
useSameRandomNumbers(true),
fitType(fitType),
nbFuncs(0), volType(volType)
{
    validatePop2Object();

}

class LeastSquareFitHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VanillaGrid::LeastSquareFit, clazz);
        SUPERCLASS(Calibrator::ObjFuncLeastSquare);
        IMPLEMENTS(Calibrator::ObjFunc::IBootstrappable);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(model, "Model");
        FIELD(vanillaGrid, "VanillaGrid instrument");
        FIELD(vanillaGridOriginal, "VanillaGrid instrument for bootstrapping routine");
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(vanillaGridOriginal);
        FIELD(useSameRandomNumbers, "use of random numbers caching for bootstrapping routine");
        FIELD_MAKE_OPTIONAL(useSameRandomNumbers);
        FIELD(fitType, "Fit Type ("+ VGLSFitTypeHelper::getNameList() +")");
        FIELD(mktVals, "Market volatilities/prices");
        FIELD_MAKE_OPTIONAL(mktVals);
        FIELD(mktVegas, "Market vegas");
        FIELD_MAKE_OPTIONAL(mktVegas);
        FIELD(control, "");
        FIELD_MAKE_TRANSIENT(control);
        FIELD(nbFuncs, "");
        FIELD_MAKE_TRANSIENT(nbFuncs);
        FIELD(fitWhat, "");
        FIELD_MAKE_TRANSIENT(fitWhat);
        FIELD(outReqName, "");
        FIELD_MAKE_TRANSIENT(outReqName);
        FIELD(volType, "volType");
        FIELD_MAKE_OPTIONAL(volType);
    }

    static IObject* defaultCtor(){
        return new VanillaGrid::LeastSquareFit();
    }
};


const string VanillaGrid::LeastSquareFit::VOLSURFACE = "VolSurface";

CClassConstSP const VanillaGrid::LeastSquareFit::TYPE = CClass::registerClassLoadMethod(
    "VanillaGrid::LeastSquareFit", typeid(VanillaGrid::LeastSquareFit), LeastSquareFitHelper::load);

/*  LeastSquareFit Ctor / Addin class */
class LeastSquareFitCtor: public CObject {
public:
    static CClassConstSP const TYPE;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(LeastSquareFitCtor, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(maturities, "Array of maturities");
        FIELD(strikes, "Array of strikes");
        FIELD(weights, "Matrix of weights");
        FIELD(asset, "Asset");
        FIELD(discount, "Discount curve");
        FIELD(fitType, "Fit Type ("
                              + VGLSFitTypeHelper::getNameList() +")");
        FIELD(model, "Model. Will de defaulted to Closed Form if not specified");
        FIELD_MAKE_OPTIONAL(model);
        FIELD(volType, "Volatility Type. Needed if Model is not specified");
        FIELD_MAKE_OPTIONAL(volType);

        // registration for addin function
        Addin::registerClassObjectMethod("LEAST_SQUARE_FIT_CREATE",
                                         Addin::RISK,
                                         "Returns a least-square objective function",
                                         LeastSquareFitCtor::TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)create);
    }

private:
    LeastSquareFitCtor(const LeastSquareFitCtor& rhs);
    LeastSquareFitCtor& operator=(const LeastSquareFitCtor& rhs);

    static IObjectSP create(LeastSquareFitCtor* params){
        if (!params->model){
            if (params->volType.empty()){
                throw ModelException("LeastSquareFitCtor::create",
                                     "Need to specify a volatility type when the model is not provided");
            }
            params->model = IModelSP(new CClosedFormLN(params->volType));
        }

        // create weightmatrix for vanillagrid
        DateTime refDate = params->discount->getToday();
        WeightMatrixSP weightMatrix(new WeightMatrix("",
                                                     params->maturities,
                                                     params->strikes,
                                                     params->weights,
                                                     WeightMatrix::ABSOLUTE,
                                                     WeightMatrix::CALL,
                                                     refDate));
        WeightMatrixWrapper weightMatrixWrapper(weightMatrix);

        VanillaGrid vanillaGrid(params->asset,
                                params->discount,
                                weightMatrixWrapper);
        return VanillaGrid::LeastSquareFitSP(new VanillaGrid::LeastSquareFit(*params->model,
                                                                             vanillaGrid,
                                                                             params->fitType,
                                                                             params->volType));
    }

    LeastSquareFitCtor():
    CObject(TYPE),
    fitType(VGLSFitTypeHelper::getDefaultName()){}

    static IObject* defaultCtor(){
        return new LeastSquareFitCtor();
    }

    DateTimeArraySP     maturities;
    DoubleArraySP       strikes;
    CDoubleMatrixSP     weights;

    CAssetWrapper       asset;
    YieldCurveWrapper   discount;
    IModelSP            model;
    string              volType;
    string              fitType;
};

CClassConstSP const LeastSquareFitCtor::TYPE = CClass::registerClassLoadMethod(
    "LeastSquareFitCtor", typeid(LeastSquareFitCtor), LeastSquareFitCtor::load);

//----------------------------------------------------------------------------------
// LeastSquareSimple. A simplified interface onto LeastSquareFit above.

typedef VanillaGrid::LeastSquareSimple::FitType VGLSSFitType;
template<> string nameForType<VGLSSFitType>(VGLSSFitType*){
    return "VanillaGrid::LeastSquareSimple::FitType";
}

typedef Enum2StringListHelper<VGLSSFitType> VGLSSFitTypeHelper;
template<> string VGLSSFitTypeHelper::names[VGLSSFitTypeHelper::EnumList::NB_ENUMS] = {
    "PRICE",
    "VEGA_NORMALIZED_PRICE",
    "IMPLIED_VOL"
};

void VanillaGrid::LeastSquareSimple::validatePop2Object(){
    static const string method("VanillaGrid::LeastSquareSimple::validatePop2Object");
    try{
        int    fitWhat = VGLSSFitTypeHelper::getIndex(fitType);
        string outReqName;

        // what request
        switch(fitWhat){
        case VGLSSFitTypeHelper::EnumList::IMPLIED_VOL:
            outReqName = OutputRequest::IND_VOL;
            break;
        default:
            outReqName = OutputRequest::OPTION_PRICE;
        }

        /* Create control */
        SensitivityArraySP sens(new SensitivityArray(0));
        OutputRequestArraySP outReqs(new OutputRequestArray(1));
        (*outReqs)[0] = OutputRequestSP(new OutputRequest(outReqName));
        control = CControlSP(new Control(sens,
                                         outReqs,
                                         false,
                                         ""));
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

void VanillaGrid::LeastSquareSimple::getMarket(MarketData* market){
    static const string method = "VanillaGrid::LeastSquareSimple::getMarket";
    try{
        CAsset::getAssetMarketData(model.get(),
                                   market,
                                   CAsset::CCY_TREATMENT_NONE,
                                   discount,
                                   asset);
        discount.getData(model.get(), market);
        instSettle->getMarket(model.get(), market);
        weightMatrix.getData(model.get(), market);

        // XXX This construction needs to be not before getMarket, and
        // XXX not after, so though it's not an obvious place for this
        // XXX construction there is little alternative.
        DateTime valueDate;
        market->GetReferenceDate(valueDate);
        // XXX Seem to be problems using the same constructor as LeastSquareFit - startDate
        // XXX and oneContract are not set as expected ...
        vanillaGrid = VanillaGridSP(new VanillaGrid(false, // fwdStarting
                                                    valueDate, // startDate
                                                    false, // oneContract
                                                    100., // notional
                                                    asset->getSpot(), // initialSpot
                                                    //*(weightMatrix->getMaturities(valueDate)),
                                                    //*(weightMatrix->getStrikes()),
                                                    //*(weightMatrix->getWeights()),
                                                    asset,
                                                    discount,
                                                    weightMatrix,
                                                    instSettle));

        leastSquareFit = LeastSquareFitSP(new LeastSquareFit(*model.get(),
                                                             *vanillaGrid.get(),
                                                             fitType,
                                                             volType));
        leastSquareFit->getMarket(market);
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

void VanillaGrid::LeastSquareSimple::validate(){
    static const string method = "VanillaGrid::LeastSquareSimple::validate";
    try{
        leastSquareFit->useNormalizedWeights(getNormalizedWeights());
        leastSquareFit->validate();
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

IObjectSP VanillaGrid::LeastSquareSimple::getAdjustableGroup(){
    return IObjectSP::attachToRef(this);
}

int VanillaGrid::LeastSquareSimple::getNbFuncs() const{
    return leastSquareFit->getNbFuncs();
}

void VanillaGrid::LeastSquareSimple::calcValue(CDoubleArray& funcvals) const{
    try{
        leastSquareFit->calcValue(funcvals);
    }
    catch(exception& e){
        throw ModelException(e, "VanillaGrid::LeastSquareSimple::calcValue");
    }
}

void VanillaGrid::LeastSquareSimple::getSensitivities(CDoubleMatrix& sensitivities,
                                                      CDoubleMatrix& valGuess,
                                                      DoubleArray shift,
                                                      int index) const{
    try{
        leastSquareFit->getSensitivities(sensitivities,valGuess,shift,index);
    }
    catch(exception& e){
        throw ModelException(e, "VanillaGrid::LeastSquareSimple::getSensitivities");
    }
}

void VanillaGrid::LeastSquareSimple::calcInitialVals(CDoubleMatrix& valGuess,
                                                     int index) const{
    try{
        leastSquareFit->calcInitialVals(valGuess,index);
    }
    catch(exception& e){
        throw ModelException(e, "VanillaGrid::LeastSquareSimple::calcInitialVals");
    }
}


void VanillaGrid::LeastSquareSimple::useNormalizedWeights(string normalize){
    normalizeWeightsLSS = normalize;
}

string VanillaGrid::LeastSquareSimple::getNormalizedWeights(){
   return normalizeWeightsLSS;
}

/** Makes additional validations for calibration with bootstrapping */
void VanillaGrid::LeastSquareSimple::validate(const Calibrator::InstanceID::IBootstrappableArray& ids) const{
    try {
        leastSquareFit->validate(ids);
    }
    catch(exception& e){
        throw ModelException(e, "VanillaGrid::LeastSquareSimple::validate(const InstanceID::IBootstrappableArray&)");
    }
}

/** Gets the first maturity with non zero weight */
const int VanillaGrid::LeastSquareSimple::getIdxFirstMat(int nberMat, int nberStrk) const{
    return leastSquareFit->getIdxFirstMat(nberMat, nberStrk);
}

IModelSP VanillaGrid::LeastSquareSimple::getModel() {
    return model;
}

VanillaGridSP VanillaGrid::LeastSquareSimple::getVanillaGrid() {
    return vanillaGrid;
}

WeightMatrixWrapper VanillaGrid::LeastSquareSimple::getWeightMatrix() {
    return weightMatrix;
}

CAssetWrapper VanillaGrid::LeastSquareSimple::getAsset() {
    return asset;
}


/** Updates the instrument for calibration with bootstrapping 
    before each run of the calibrator */
void VanillaGrid::LeastSquareSimple::update(int idxSmileMat) {
    leastSquareFit->update(idxSmileMat);
}

/** Reset the instrument for calibration with bootstrapping
    after each run of the calibrator */
void VanillaGrid::LeastSquareSimple::reset() {
    leastSquareFit->reset();
}

// for reflection
VanillaGrid::LeastSquareSimple::LeastSquareSimple():
Calibrator::ObjFuncLeastSquare(TYPE),
fitType(VGLSSFitTypeHelper::getDefaultName()),
useSameRandomNumbers(true), volType(VanillaGrid::LeastSquareFit::VOLSURFACE){}

VanillaGrid::LeastSquareSimple::LeastSquareSimple(const IModel&                   model,
                                                  string                         fitType,
                                                  bool                           useSameRandomNumbers,
                                                  const CAssetWrapper&           asset,
                                                  const YieldCurveWrapper&       discount,
                                                  const InstrumentSettlementSP&  instSettle,
                                                  const WeightMatrixWrapper&     weightMatrix,
                                                  const string&                  volType):
    Calibrator::ObjFuncLeastSquare(TYPE),
                model(copy(&model)),
                fitType(fitType),
                useSameRandomNumbers(true),
                asset(asset),
                discount(discount),
                instSettle(instSettle),
                weightMatrix(weightMatrix),
                volType(volType){
    validatePop2Object();
}

class LeastSquareSimpleHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VanillaGrid::LeastSquareSimple, clazz);
        SUPERCLASS(Calibrator::ObjFuncLeastSquare);
        IMPLEMENTS(Calibrator::ObjFunc::IBootstrappable);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(model, "Model");
        FIELD(fitType, "Fit Type ("+ VGLSSFitTypeHelper::getNameList() +")");
        FIELD(useSameRandomNumbers, "use of random numbers caching for bootstrapping routine");
        FIELD_MAKE_OPTIONAL(useSameRandomNumbers);
        FIELD(asset, "Asset");
        FIELD(discount, "Discount yield curve");
        FIELD(instSettle, "Instrument settlement at maturity");
        FIELD(weightMatrix, "Calibration weight matrix");
        FIELD(mktVals, "Market volatilities/prices");
        FIELD_MAKE_OPTIONAL(mktVals);
        FIELD(mktVegas, "Market vegas");
        FIELD_MAKE_OPTIONAL(mktVegas);
        FIELD(vanillaGrid, "VanillaGrid instrument used for bootstrapping routine");
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(vanillaGrid);
        FIELD_NO_DESC(leastSquareFit);
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(leastSquareFit); // Should this be the TRANSIENT_BUT_TWEAKABLE?
        FIELD(control, "");
        FIELD_MAKE_TRANSIENT(control);
        FIELD(volType, "volType");
        FIELD_MAKE_OPTIONAL(volType);
    }

    static IObject* defaultCtor(){
        return new VanillaGrid::LeastSquareSimple();
    }
};

CClassConstSP const VanillaGrid::LeastSquareSimple::TYPE = CClass::registerClassLoadMethod(
    "VanillaGrid::LeastSquareSimple", typeid(VanillaGrid::LeastSquareSimple), LeastSquareSimpleHelper::load);

//----------------------------------------------------------------------------------

bool VanillaGridLoad() {
    return (VanillaGrid::TYPE != 0
            && VanillaGrid::LeastSquareFit::TYPE != 0
            && LeastSquareFitCtor::TYPE != 0);
}

/* Addin. VanillaGridAddin is already defined in VolCheck */
class VanillaGridAddin2: public CObject{
public:
    static CClassConstSP const TYPE;

    struct RequestType{
        enum {
            OPTION_PRICE = 0,
            IND_VOL,
            OPTION_VEGA,
            NB_ENUMS
        };
    };
    typedef Enum2StringListHelper<RequestType> RequestTypeHelper;

    virtual void validatePop2Object();

private:

    CMarketDataSP       market;
    CAssetWrapper       assetWrapper;
    YieldCurveWrapper   yieldWrapper;
    bool                isFwdStarting;
    DateTime            startDate;
    bool                isCall;
    bool                oneContract;
    double              notional;
    double              initialSpot;
    DateTimeArray       maturities;
    CDoubleArray        strikes;
    CDoubleMatrixSP     weights;
    IModelSP            model;
    string              volType;
    StringArray         requests;


    static IObjectSP prices(VanillaGridAddin2* params) {
        static const string routine = "VanillaGridAddin2::prices";
        try {
            // create weight matrix for vanillagrid
            DateTimeArraySP maturitiesSP(copy(&(params->maturities)));
            DoubleArraySP strikesSP(copy(&(params->strikes)));
            string instType;
            if (params->isCall){
                instType = WeightMatrix::CALL;
            } else {
                instType = WeightMatrix::PUT;
            }
            WeightMatrixSP weightMatrix(new WeightMatrix("",
                                                         maturitiesSP,
                                                         strikesSP,
                                                         params->weights,
                                                         WeightMatrix::ABSOLUTE,
                                                         instType,
                                                         params->startDate));
            WeightMatrixWrapper weightMatrixWrapper(weightMatrix);

            VanillaGrid inst(params->isFwdStarting,
                             params->startDate,
                             params->oneContract,
                             params->notional,
                             params->initialSpot,
                             params->assetWrapper,
                             params->yieldWrapper,
                             weightMatrixWrapper);

            // control
            StringArray& requests = params->requests;
            int nbReqs = requests.size();
            OutputRequestArraySP outputArray(new OutputRequestArray(nbReqs));
            int iReq = 0;
            for (; iReq < nbReqs; ++iReq){
                (*outputArray)[iReq] = OutputRequestSP(new OutputRequest(requests[iReq]));
            }
            Control control(SensitivityArrayConstSP(new SensitivityArray(0)),
                            outputArray,
                            false,
                            string(""));

            // get market and price
            CResultsSP results(
                params->model->go(CInstrumentSP::attachToRef(&inst),
                                  ScenarioSP(),
                                  CControlSP::attachToRef(&control),
                                  params->market));

            ObjectArraySP rtn(new ObjectArray(nbReqs));
            for (iReq = 0; iReq < nbReqs; ++iReq){
                IObjectConstSP obj(results->retrieveRequestResult(requests[iReq]));
                (*rtn)[iReq] = IObjectSP(copy(obj.get()));
            }
            return rtn;
        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }


    /** for reflection */
    VanillaGridAddin2():
    CObject(TYPE),
    isFwdStarting(false),
    isCall(true),
    volType("VolSurface"){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(VanillaGridAddin2, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultVanillaGridAddin2);
        FIELD(market, "market");
        FIELD(assetWrapper, "Asset wrapper");
        FIELD(yieldWrapper, "Yield curve wrapper");
        FIELD(isCall, "Call if true; put, otherwise");
        FIELD(isFwdStarting, "Fwd starting if true; Started, otherwise");
        FIELD(startDate, "Start Date");
        FIELD_MAKE_OPTIONAL(startDate);
        FIELD(oneContract, "oneContract");
        FIELD(notional, "notional");
        FIELD_MAKE_OPTIONAL(notional);
        FIELD(initialSpot, "initialSpot");
        FIELD_MAKE_OPTIONAL(initialSpot);
        FIELD(maturities, "Array of maturities");
        FIELD(strikes, "Array of strikes");
        FIELD(weights, "weights");
        FIELD_MAKE_OPTIONAL(weights);
        FIELD(model, "Model");
        FIELD_MAKE_OPTIONAL(model);
        FIELD(volType, "Type of vol to use");
        FIELD_MAKE_OPTIONAL(volType);
        FIELD(requests, "Requests");

        Addin::registerClassObjectMethod("VANILLAGRID_PRICE",
                                         Addin::RISK,
                                         "IMCPrices a Vanilla Grid (no ccy treatment)",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)prices);

    }

    static IObject* defaultVanillaGridAddin2(){
        return new VanillaGridAddin2();
    }
};

template<> string nameForType(VanillaGridAddin2::RequestType*){
    return "VanillaGridAddin2::RequestType";
}

template<> string VanillaGridAddin2::RequestTypeHelper::names[
    VanillaGridAddin2::RequestTypeHelper::EnumList::NB_ENUMS] = {
        "OPTION_PRICE",
        "IND_VOL",
        "OPTION_VEGA"
};

// need to define / specialize nameForType prior to actually
// instantiating Enum2StringListHelper (the instantiation is
// implicitly done when Enum2StringListHelper::getIndex()
// gets called below
void VanillaGridAddin2::validatePop2Object(){
    static const string method("VanillaGridAddin2::ValidatePop2Object");
    try{
        model = IModelSP(new CClosedFormLN(volType));
        // check request types
        int nbReqs = requests.size();
        if (requests.size() == 0){
            throw ModelException(method,
                                 "At least one price request must be provided");
        }
        int iReq = 0;
        for (; iReq < nbReqs; ++iReq){
            // will throw exception if request type is not recognized
            RequestTypeHelper::getIndex(requests[iReq]);
        }
        // if no weights provided,
        if (!weights){
            int nbStrikes = strikes.size();
            int nbSteps   = maturities.size();
            weights = CDoubleMatrixSP(new CDoubleMatrix(nbSteps, nbStrikes));
            int iStep = 0;
            for (; iStep < nbSteps; ++iStep){
                int iStrike = 0;
                for (; iStrike < nbStrikes; ++iStrike){
                    (*weights)[iStep][iStrike] = 1.0;
                }
            }
        }
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

CClassConstSP const VanillaGridAddin2::TYPE = CClass::registerClassLoadMethod(
    "VanillaGridAddin2", typeid(VanillaGridAddin2), load);


////////////////////////////////////////////////////////////////////////////


/** Creates an Implied VolSurface corresponding to Model Parameters e.g.
    VolSurface corresponding to VolSV */
class ModelVolSurface: public CObject{
public:
    /** Returns output of addin function */
    class Output: public CObject {
    public:
        static CClassConstSP const TYPE;

        /** Default constructor */
        Output(): CObject(TYPE) {}

        /** Full constructor */
        Output(VolSurfaceSP volSurface,
               DoubleArrayArraySP fittingError,
               IntArrayArraySP interpFlags): 
        CObject(TYPE), volSurface(volSurface), fittingError(fittingError), interpFlags(interpFlags) {}

    private:
        /** Invoked when Class is 'loaded' */
        static void load(CClassSP& clazz){
            REGISTER(ModelVolSurface::Output, clazz);
            SUPERCLASS(CObject);
            EMPTY_SHELL_METHOD(defaultModelVolSurfaceOutput);
            FIELD(volSurface, "Vol Surface");
            FIELD(fittingError, "Fitting error");
            FIELD(interpFlags, "Interpolation flags");
        }

        /** Default constructor */
        static IObject* defaultModelVolSurfaceOutput(){
            return new Output();
        }

        VolSurfaceSP        volSurface;         //!< Model vol surface
        DoubleArrayArraySP  fittingError;       //!< Model - Market implied vol
        IntArrayArraySP     interpFlags;        //!< Whether model vol is evaluated or interpolated
    };

    typedef smartPtr<Output> OutputSP;
    
    static CClassConstSP const TYPE;

    ModelVolSurface::OutputSP createVolsurface() {
        static const string method = "ModelVolSurface::createVolsurface";
        try {
            // Vol surface i.e. original market vol
            CAssetWrapper assetVolPreferred(asset.getName());
            IModelSP modelVolPreferred(new CClosedFormLN("VolPreferred"));
            CAsset::getAssetMarketData(modelVolPreferred.get(), 
                                       market.get(), 
                                       CAsset::CCY_TREATMENT_NONE,
                                       yieldCurve,
                                       assetVolPreferred);
            
            CAssetWrapper assetVolSpline(asset.getName());
            
            // 1) Create components of Vanilla Grid
            DateTime valueDate;
            market->GetReferenceDate(valueDate);

            // Extract VolSurface from market cache
            IModelSP modelLN(new CClosedFormLN("VolSurface"));
            vol.getData(modelLN.get(), market);
            VolSurfaceSP volSurface = VolSurfaceSP::dynamicCast(vol.getSP());
            TimeMetricSP metric = volSurface->getTimeMetric();
            const DateTimeArray& dates = volSurface->getDates();
            const DoubleArray& strikes = volSurface->getStrikes();

            CAsset::getAssetMarketData(model.get(), 
                                       market.get(), 
                                       CAsset::CCY_TREATMENT_NONE,
                                       yieldCurve,
                                       asset);
            yieldCurve.getData(model.get(), market);
            
            InstrumentSettlementSP instSettle(new CashSettlePeriod(0));
            instSettle->getMarket(model.get(), market.get());
            
            int numStrikes = strikes.size();
            int numMaturities = dates.size();
            CDoubleMatrixSP weights(new CDoubleMatrix(numMaturities, numStrikes));
            weights->scalarAdd(1.0);
            
            DateTimeArraySP datesSP(new DateTimeArray(dates));

            DoubleArraySP strikesSP(new DoubleArray(strikes));
            WeightMatrixSP weightMatrix(new WeightMatrix(
                 asset.getName(), 
                 datesSP,
                 strikesSP,
                 weights,
                 WeightMatrix::ABSOLUTE,
                 WeightMatrix::CALL,
                 valueDate));
            WeightMatrixWrapper weightMatrixWrapper(weightMatrix);
            
            // 2) Create Vanilla Grid
            VanillaGridSP vanillaGridNew(new VanillaGrid(false,              // fwdStarting
                                                         valueDate,          // startDate
                                                         false,              // oneContract
                                                         100.,               // notional
                                                         asset->getSpot(),   // initialSpot
                                                         asset,
                                                         yieldCurve,
                                                         weightMatrixWrapper,
                                                         instSettle));
            
            // 3) Price Vanilla Grid
            OutputRequestSP request(new OutputRequest(OutputRequest::IND_VOL));
            OutputRequestArrayConstSP outReqs(new OutputRequestArray(1, request));
            SensitivityArrayConstSP sens(new SensitivityArray(0));
            CControlSP control(new Control(sens, outReqs, false, ""));
            
            CResultsSP results = model->go(CInstrumentSP(vanillaGridNew),
                                           ScenarioSP(0),
                                           control,
                                           market);
            
            // 4) Extract Implied Vols
            IObjectConstSP impVolsObj = results->retrieveRequestResult(OutputRequest::IND_VOL);
            VanillaGrid::OutputConstSP impVols = VanillaGrid::OutputConstSP::dynamicCast(impVolsObj);
            
            // 5) Spline Vols, interpolate missing vols, create new VolSurface
            DoubleMatrix newVols(impVols->numRows(), impVols->numCols());
            IntArrayArraySP interpFlags(new IntArrayArray(impVols->numCols()));
            DoubleArrayArraySP fittingErrors(new DoubleArrayArray(impVols->numCols()));
            IModelSP modelVolSpline(new CClosedFormLN("VolSpline"));
            for(int iDate = 0; iDate< impVols->numCols(); iDate++) {
                try {
                    (*interpFlags)[iDate] = IntArraySP(new IntArray(impVols->numRows()));
                    (*fittingErrors)[iDate] = DoubleArray(impVols->numRows());
                    
                    // Get strikes & vols
                    DoubleArray tmpStrikes;
                    DoubleArray tmpVols;
                    int iStrike;
                    for(iStrike = 0; iStrike < impVols->numRows(); iStrike++) {
                        bool isError = impVols->isError(iDate, iStrike);
                        if(!isError) {
                            tmpStrikes.push_back(strikes[iStrike]);
                            tmpVols.push_back(impVols->getValue(iDate, iStrike));
                            (*(*interpFlags)[iDate])[iStrike] = 0;
                        } else {
                            (*(*interpFlags)[iDate])[iStrike] = 1;
                        }
                    }

                    DoubleMatrix tmpVolMatrix(tmpVols);
                    tmpVolMatrix.transpose();
                    
                    // Spline and interpolate
                    ExpiryArray tmpExpiries(1, ExpirySP(new BenchmarkDate(dates[iDate])));
                    VolSurfaceSP tmpVolSurface(new VolSurface(
                        vol.getName(), metric.get(), tmpStrikes, tmpVolMatrix, &tmpExpiries, valueDate));
                    
                    smartPtr<VolSpline> spline(new VolSpline(*tmpVolSurface, asset->getSpot()));
                    
                    IObjectSP tmpMarketObj(market->clone());
                    CMarketDataSP tmpMarket = CMarketDataSP::dynamicCast(tmpMarketObj);
                    tmpMarket->AddData(spline);
                    CAsset::getAssetMarketData(modelVolSpline.get(), 
                                               tmpMarket.get(), 
                                               CAsset::CCY_TREATMENT_NONE,
                                               yieldCurve,
                                               assetVolSpline);

                    VolRequestLNStrikeSP volRequest(new LinearStrikeTSVolRequest(
                        1.0,
                        valueDate, 
                        dates[iDate], // used ?
                        false));

                    CDoubleMatrixSP matrix;
                    for(iStrike = 0; iStrike < impVols->numRows(); iStrike++) {
                        // New splined vol
                        volRequest->setStrike(strikes[iStrike]);
                        CVolProcessedBSSP volProcessedVolSpline(assetVolSpline->getProcessedVol(volRequest.get()));
                        double vol = volProcessedVolSpline->CalcVol(valueDate, dates[iDate]);
                        newVols[iStrike][iDate] = vol;

                        // Market vol
                        CVolProcessedBSSP volProcessedVolPreferred(assetVolPreferred->getProcessedVol(volRequest.get()));
                        double marketVol = volProcessedVolPreferred->CalcVol(valueDate, dates[iDate]);
                        (*fittingErrors)[iDate][iStrike] = vol - marketVol;
                    }
                } catch(exception& e) {
                    string message = "Failed when Splining model Vols at expiry " + dates[iDate].toString();
                    throw ModelException::addTextToException(e, message);
                }
            }

            // 6) Create Results
            VolSurfaceSP newVolSurface(new VolSurface(vol.getName(), metric.get(), 
                strikes, newVols, volSurface->getExpiries().get(), valueDate));

            OutputSP output(new Output(newVolSurface, fittingErrors, interpFlags));
            return output;
        } catch(exception& e) {
            throw ModelException(e, method);
        }
    }

private:
    // Mandatory fields
    CMarketDataSP       market;         //!< Market Data
    IModelSP            model;          //!< Model
    CAssetWrapper       asset;          //!< Asset
    CVolBaseWrapper     vol;            //!< Vol
    YieldCurveWrapper   yieldCurve;     //!< Yield curve
    
    /** for reflection */
    ModelVolSurface(): CObject(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ModelVolSurface, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultModelVolSurface);
        FIELD(market, "Market cache");
        FIELD(model, "Model");
        FIELD(asset, "Asset");
        FIELD(vol, "Vol");
        FIELD(yieldCurve, "Yield curve");

        Addin::registerObjectMethod("MODEL_VOLSURFACE",
                                    Addin::MARKET,
                                    "Creates a VolSurface from pricing a VanillaGrid",
                                    true,
                                    Addin::returnHandle,
                                    &ModelVolSurface::createVolsurface);
    }

    /** Default constructor */
    static IObject* defaultModelVolSurface(){
        return new ModelVolSurface();
    }
};


CClassConstSP const ModelVolSurface::TYPE = CClass::registerClassLoadMethod(
    "ModelVolSurface", typeid(ModelVolSurface), load);


CClassConstSP const ModelVolSurface::Output::TYPE = CClass::registerClassLoadMethod(
    "ModelVolSurface::Output", typeid(ModelVolSurface::Output), load);


////////////////////////////////////////////////////////////////////////


/** Creates an Implied VolSurface corresponding to Model Parameters e.g.
    VolSurface corresponding to VolSV */
class ModelScenario: public CObject{
public:
    static CClassConstSP const TYPE;

    /** Returns output of addin function */
    class Output: public CObject {
    public:
        static CClassConstSP const TYPE;

        /** Default constructor */
        Output(): CObject(TYPE) {}

        /** Full constructor */
        Output(CResultsArraySP  paramModelResults,
               CResultsArraySP  convModelAtParamVolSurfaceResults,
               CResultsArraySP  convModelResults): 
        CObject(TYPE), paramModelResults(paramModelResults), 
            convModelAtParamVolSurfaceResults(convModelAtParamVolSurfaceResults), convModelResults(convModelResults) {}

    private:
        /** Invoked when Class is 'loaded' */
        static void load(CClassSP& clazz){
            REGISTER(ModelScenario::Output, clazz);
            SUPERCLASS(CObject);
            EMPTY_SHELL_METHOD(defaultModelScenarioOutput);
            FIELD(paramModelResults, "Parametric model results");
            FIELD(convModelAtParamVolSurfaceResults, "Conventional model at parametric vol surface results");
            FIELD(convModelResults, "Conventional model results");
        }

        /** Default constructor */
        static IObject* defaultModelScenarioOutput(){
            return new Output();
        }

        CResultsArraySP  paramModelResults;                     //!< Parametric model results
        CResultsArraySP  convModelAtParamVolSurfaceResults;     //!< Conventional model at parametric vol surface results
        CResultsArraySP  convModelResults;                      //!< Conventional model results
    };

    typedef smartPtr<Output> OutputSP;


    ModelScenario::OutputSP runModelScenario() {
        static const string method = "ModelScenario::runModelScenario";
        try {
            // 1) Price using both models and identical market data
            CResultsArraySP conventionalResults = conventionalModel->go(inst, ScenarioSP(0), control, market);
            CResultsArraySP parametricResults = parametricModel->go(inst, ScenarioSP(0), control, market);

            // Create new market by adding market data
            IObjectSP newMarketObj(market->clone());
            CMarketDataSP newMarket = CMarketDataSP::dynamicCast(newMarketObj);
            int iObj;
            for(iObj = 0; iObj < marketObjects->size(); iObj++) {
                newMarket->AddData((*marketObjects)[iObj]);
            }

            CResultsArraySP conventionalNewResults = conventionalModel->go(inst, ScenarioSP(0), control, newMarket);

            // Return results
            OutputSP output(new Output(parametricResults, conventionalNewResults, conventionalResults));
            return output;
        } catch(exception& e) {
            throw ModelException(e, method);
        }
    }

private:
    // Mandatory fields
    CMarketDataSP           market;                 //!< Market Data
    IModelSP                parametricModel;        //!< Parametric model e.g. FourierEngine VolSV
    IModelSP                conventionalModel;      //!< Conventional model e.g. MCLocalVol
    IInstrumentCollectionSP inst;                   //!< Instrument
    CControlSP              control;                //!< Control
    MarketObjectArraySP     marketObjects;          //!< New market objects for conventional model
    
    /** for reflection */
    ModelScenario(): CObject(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ModelScenario, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultModelScenario);
        FIELD(market, "Market cache");
        FIELD(parametricModel, "Parametric model e.g. StochVol");
        FIELD(conventionalModel, "Conventional model e.g. Local Vol");
        FIELD(inst, "Instrument");
        FIELD(control, "Control");
        FIELD(marketObjects, "Market Objects for conventional model e.g. VolSurface for StochVol");

        Addin::registerObjectMethod("MODEL_SCENARIO",
                                    Addin::RISK,
                                    "Model scenario on instrument using 2 models",
                                    true,
                                    Addin::returnHandle,
                                    &ModelScenario::runModelScenario);
    }

    /** Default constructor */
    static IObject* defaultModelScenario(){
        return new ModelScenario();
    }
};


CClassConstSP const ModelScenario::TYPE = CClass::registerClassLoadMethod(
    "ModelScenario", typeid(ModelScenario), load);


CClassConstSP const ModelScenario::Output::TYPE = CClass::registerClassLoadMethod(
    "ModelScenario::Output", typeid(ModelScenario::Output), load);


DRLIB_END_NAMESPACE

