//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolCheck.cpp
//
//   Description : Vol Check
//
//   Author      : Regis Guichard
//
//   Date        : 03 Jan 03
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VanillaGrid.hpp"
#include "edginc/AtMaturity.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/Scenario.hpp"
#include "edginc/Weightmatrix.hpp"
#include "edginc/BenchmarkDate.hpp"
#include "edginc/ClientRunnable.hpp"

DRLIB_BEGIN_NAMESPACE

//////////////////////////////////////////
// Information on VolSurface validation
//////////////////////////////////////////
class ArbitrageChecks: public CObject {
public:
    static CClassConstSP const TYPE;

    // For reflection
    ArbitrageChecks(): CObject(TYPE) {}

    // constructor from calls
    ArbitrageChecks(const VanillaGrid::Output&  callsobj,
                    const CDoubleArraySP&       strikes,
                    const CAssetWrapper&        assetWrapper,
                    const DateTime              valueDate,
                    const DateTimeArraySP&      maturities,
                    const double                accuracy):
    CObject(TYPE), strikes(strikes), ok(true) {
        static string routine = "ArbitrageChecks::ArbitrageChecks";

        try {
            // copy values over and transpose
            // warning: may throw an exception if elt (iRow, iCol) of callsobj is an error
            calls = CDoubleMatrixSP(new DoubleMatrix(callsobj.numRows(), callsobj.numCols()));
            int iCol = 0;
            for (; iCol < callsobj.numCols(); ++iCol){
                int iRow = 0;
                for (; iRow < callsobj.numRows(); ++iRow){
                    (*calls)[iRow][iCol] = callsobj.getValue(iCol, iRow);
                }
            }

            CDoubleArray& strikes = *this->strikes;
            int nbStrikes = strikes.size();
            int nbSteps   = calls->numRows();

            int iStrike;
            for(iStrike = 1; iStrike < nbStrikes; iStrike++) {
                if(!Maths::isPositive(strikes[iStrike] - strikes[iStrike - 1])) {
                    throw ModelException("Strikes must be in strictly increasing order");
                }
            }

            if(calls->numCols() != nbStrikes) {
                throw ModelException("Number of strikes must equal number of columns in call prices matrix");
            }

            // create the string messages
            messages = StringArrayArraySP(new StringArrayArray(nbStrikes));
            for(iStrike = 0; iStrike < nbStrikes; iStrike++) {
                (*messages)[iStrike] = StringArraySP(new StringArray(nbSteps));
            }

            int iStep;

            if(nbStrikes >= 1) {
                // create callSpreads
                callSpreads = CDoubleMatrixSP(new DoubleMatrix(nbStrikes - 1, nbSteps));
                DoubleMatrix& callSpreadsMatrix = *callSpreads;
                DoubleMatrix& callMatrix        = *calls;
                for(iStrike = 1; iStrike < nbStrikes; iStrike++) {
                    double dStrike = strikes[iStrike] - strikes[iStrike - 1];
                    StringArray& stringArray = *(*messages)[iStrike];
                    for(iStep = 0; iStep < nbSteps; iStep++) {
                        callSpreadsMatrix[iStrike - 1][iStep] =
                            (callMatrix[iStrike - 1][iStep] - callMatrix[iStrike][iStep]) / dStrike;

                        // check callSpreads
                        if(Maths::isNegative(callSpreadsMatrix[iStrike - 1][iStep] + accuracy)) {
                            stringArray[iStep] = ArbitrageChecks::negCallSpreads;
                            ok = false;
                        } else if(Maths::isNegative(1.0 + accuracy - callSpreadsMatrix[iStrike - 1][iStep])) {
                            stringArray[iStep] = ArbitrageChecks::callSpreadsExceedOne;
                            ok = false;
                        }
                    }
                }
                // check that implied variances are an increasing function of time
                impliedVariances = CDoubleMatrixSP(new DoubleMatrix(nbStrikes, nbSteps));
                DoubleArray fwdValues(nbSteps);
                assetWrapper->fwdValue(*maturities, fwdValues);
                double spot = assetWrapper->fwdValue(valueDate);
                LinearStrikeVolRequestSP volRequest(new LinearStrikeVolRequest(spot,
                                                                               valueDate,
                                                                               (*maturities)[nbSteps-1],
                                                                               false)); // isfwdstarting
                // allowNegativeFwdVar = true, since we want to detect sthg like this
                volRequest->allowNegativeFwdVar(true);
                for (iStrike = 0; iStrike < nbStrikes; iStrike++) {
                    StringArray& stringArray = *(*messages)[iStrike];
                    // calculate precentage strike
                    double percentageStrike = strikes[iStrike] / spot;
                    for (iStep = 0; iStep < nbSteps; iStep++) {
                        // calculate current strike (use percentage strike and fwd (rather than spot))
                        double currentStrike = fwdValues[iStep] * percentageStrike;
                        volRequest->setStrike(currentStrike); // absolute strike
                        CVolProcessedBSSP volBS(assetWrapper->getProcessedVol(volRequest.get()));
                        (*impliedVariances)[iStrike][iStep] = volBS->CalcVar(valueDate, (*maturities)[iStep]);
                        if(iStep>0) {
                            double diff =  (*impliedVariances)[iStrike][iStep] - (*impliedVariances)[iStrike][iStep-1];
                            if (Maths::isNegative(diff + accuracy)){
                                stringArray[iStep] += ArbitrageChecks::negFwdVar;
                                ok = false;
                                }
                        }
                    }
                }

                // create bflySpreads
                if(nbStrikes >= 2) {
                    bflySpreads = CDoubleMatrixSP(new DoubleMatrix(nbStrikes - 2, nbSteps));
                    DoubleMatrix& bflySpreadsMatrix = *bflySpreads;
                    for(iStrike = 2; iStrike < nbStrikes; iStrike++) {
                        StringArray& stringArray = *(*messages)[iStrike];
                        double dStrike = strikes[iStrike] - strikes[iStrike - 2];
                        for(iStep = 0; iStep < nbSteps; iStep++) {
                            bflySpreadsMatrix[iStrike - 2][iStep] =
                                (callSpreadsMatrix[iStrike - 2][iStep] - callSpreadsMatrix[iStrike - 1][iStep]) / dStrike;

                            // check bflySpreads
                            if(Maths::isNegative(bflySpreadsMatrix[iStrike - 2][iStep] + accuracy)) {
                                stringArray[iStep] += ArbitrageChecks::negBflySpreads;
                                ok = false;
                            }
                        }
                    }
                }
            }

            // create OK messages for all the rest
            for(iStrike = 0; iStrike < nbStrikes; iStrike++) {
                StringArray& stringArray = *(*messages)[iStrike];
                for(iStep = 0; iStep < nbSteps; iStep++) {
                    if(stringArray[iStep].empty()) {
                        stringArray[iStep] = ArbitrageChecks::okMessage;
                    }
                }
            }
        } catch(exception& e){
            throw ModelException(e, routine);
        }
    }

    static string nullMessage;
    static string callSpreadsExceedOne;
    static string negCallSpreads;
    static string negBflySpreads;
    static string negFwdVar;
    static string okMessage;

private:
    static void load(CClassSP& clazz){
        REGISTER(ArbitrageChecks, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultArbitrageChecks);
        FIELD(strikes, "Strikes");
        FIELD(calls, "Call prices");
        FIELD(messages, "Error messages");
        FIELD(callSpreads, "Call spread prices");
        FIELD(bflySpreads, "Butterfly spread prices");
        // FIELD(impliedVariances, "Implied Variances");
        FIELD(ok, "Is VolSurface arbitrage free");
    }

    static IObject* defaultArbitrageChecks(){
        return new ArbitrageChecks();
    }

    CDoubleArraySP      strikes;
    CDoubleMatrixSP     calls;
    StringArrayArraySP  messages;
    CDoubleMatrixSP     callSpreads;
    CDoubleMatrixSP     bflySpreads;
    CDoubleMatrixSP     impliedVariances; // $unregistered
    bool                ok;
};

string ArbitrageChecks::nullMessage = "-";
string ArbitrageChecks::callSpreadsExceedOne = "CS > 100%. ";
string ArbitrageChecks::negCallSpreads = "CS < 0%. ";
string ArbitrageChecks::negBflySpreads = "BS < 0. ";
string ArbitrageChecks::negFwdVar = "FwdVar < 0";
string ArbitrageChecks::okMessage = "OK";

typedef smartPtr<ArbitrageChecks> ArbitrageChecksSP;
typedef smartConstPtr<ArbitrageChecks> ArbitrageChecksConstSP;

CClassConstSP const ArbitrageChecks::TYPE = CClass::registerClassLoadMethod(
    "ArbitrageChecks", typeid(ArbitrageChecks), load);


/////////////////////////////////
// VANILLA GRID ADDIN
/////////////////////////////////
/** ADDIN method for interpolating with some interpolant */
class VanillaGridAddin: public CObject{
public:
    static CClassConstSP const TYPE;
    friend class VolGridChecker;

private:

    CMarketDataSP       market;
    // CControlSP          control;
    DateTimeArraySP     maturities;
    CDoubleArraySP      strikes;

    CAssetWrapper       assetWrapper;
    YieldCurveWrapper   yieldWrapper;

    // string              asset;
    // string              discount;
    string              volType;
    bool                allowNegativeFwdVar;
    string              ccyTreatment;
    double              accuracy;

    static IObjectSP prices(VanillaGridAddin* params) {
        static const string routine = "VanillaGridAddin::prices";
        try {
            if(!Maths::equals((*params->strikes)[0], 0.0)) {
                // add zero strike
                vector<double>::iterator begin = params->strikes->begin();
                params->strikes->insert(begin, 0.0);
            }

            int nbStrikes = params->strikes->size();
            int nbSteps   = params->maturities->size();

            // model
            CClosedFormLN model(params->volType,
                                params->allowNegativeFwdVar);

            // instrument
            CDoubleMatrix weights(nbSteps, nbStrikes);
            int iStep = 0;
            for (; iStep < nbSteps; ++iStep){
                int iStrike = 0;
                for (; iStrike < nbStrikes; ++iStrike){
                    weights[iStep][iStrike] = 1.0;
                }
            }
            // CAssetWrapper assetWrapper(params->asset);
            // YieldCurveWrapper yieldWrapper(params->discount);
            InstrumentSettlementSP instSettle(new AtMaturity());

            // create a weightmatrix for the vanillagrid
            DateTime refDate;
            params->market->GetReferenceDate(refDate);
            CDoubleMatrixSP weightsSP(copy(&weights));
            WeightMatrixSP weightMatrix(new WeightMatrix("", 
                                                        params->maturities,
                                                        params->strikes,
                                                        weightsSP,
                                                        WeightMatrix::ABSOLUTE,
                                                        WeightMatrix::CALL,
                                                        refDate));
            WeightMatrixWrapper weightMatrixWrapper(weightMatrix);
            
            VanillaGrid inst(params->assetWrapper,
                             params->yieldWrapper,
                             weightMatrixWrapper,
                             instSettle,
                             params->ccyTreatment);

            // control
            OutputRequestArraySP outputArray(new OutputRequestArray(1));
            (*outputArray)[0] = OutputRequestSP(new OutputRequest(OutputRequest::OPTION_PRICE));

            Control control(SensitivityArrayConstSP(new SensitivityArray(0)),
                            outputArray,
                            false,
                            string(""));

            CResultsSP results(model.go(CInstrumentSP::attachToRef(&inst),
                                        ScenarioSP(),
                                        CControlSP::attachToRef(&control),
                                        //params->control.get(),
                                        params->market));

        IObjectConstSP callsObj(results->retrieveRequestResult(OutputRequest::OPTION_PRICE));
        VanillaGrid::OutputConstSP calls(VanillaGrid::OutputConstSP::dynamicCast(callsObj));

        // still a null pointer, in parts
        params->assetWrapper.getData(&model, params->market);
        DateTime valueDate = params->market->GetReferenceDate();

        // Initialize the output packet
        ArbitrageChecksSP checks(new ArbitrageChecks(*calls,
                                                     params->strikes,
                                                     params->assetWrapper,
                                                     valueDate,
                                                     params->maturities,
                                                     params->accuracy));

        return IObjectSP(checks.clone());

        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }


    /** for reflection */
    VanillaGridAddin():  CObject(TYPE), 
        allowNegativeFwdVar(false), 
        ccyTreatment(CAsset::CCY_TREATMENT_NONE),
        accuracy(1.0e-10) {}

    /** for VolGridChecker */
    VanillaGridAddin( CMarketDataSP     market,
                     DateTimeArraySP   maturities,
                     CDoubleArraySP    strikes,
                     CAssetWrapper     assetWrapper,
                     YieldCurveWrapper yieldWrapper,
                     string            volType ) :
        CObject(TYPE),
        market(market),
        maturities(maturities), 
        strikes(strikes),
        assetWrapper(assetWrapper),
        yieldWrapper(yieldWrapper),
        volType(volType),
        allowNegativeFwdVar(false), 
        ccyTreatment(CAsset::CCY_TREATMENT_NONE),
        accuracy(1.0e-10) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(VanillaGridAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultVanillaGridAddin);
        FIELD(market, "market");
        //FIELD(control, "control");
        FIELD(maturities, "Array of maturities");
        FIELD(strikes, "Array of strikes");
        FIELD(assetWrapper, "Asset wrapper");
        FIELD(yieldWrapper, "Yield curve wrapper");
        // FIELD(asset, "Asset");
        // FIELD(discount, "Discount curve");
        FIELD(volType, "Type of vol to use");
        FIELD(allowNegativeFwdVar, "allow negative fwd var");
        FIELD_MAKE_OPTIONAL(allowNegativeFwdVar);
        FIELD(ccyTreatment, "Currency treatment");
        FIELD_MAKE_OPTIONAL(ccyTreatment);
        FIELD(accuracy, "accuracy of arbitrage checks");
        FIELD_MAKE_OPTIONAL(accuracy);

        Addin::registerClassObjectMethod("VANILLAGRID_CLOSEDFORM_LN",
                                         Addin::RISK,
                                         "Prices a started Vanilla Grid",
                                         TYPE,
                                         false,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)prices);

    }

    static IObject* defaultVanillaGridAddin(){
        return new VanillaGridAddin();
    }
};

typedef smartPtr<VanillaGridAddin> VanillaGridAddinSP;
typedef smartConstPtr<VanillaGridAddin> VanillaGridAddinConstSP;

CClassConstSP const VanillaGridAddin::TYPE = CClass::registerClassLoadMethod(
    "VanillaGridAddin", typeid(VanillaGridAddin), load);


/////////////////////////////////
// VOL GRID CHECKER
/////////////////////////////////

/** ADDIN method for checking whether vol grids are arbitrage free */
class VolGridChecker: public CObject, public ClientRunnable {
public:
    static CClassConstSP const TYPE;

    // Override the <<ClientRunnable>> interface
    virtual IObjectSP run() {
        return check(this);
    }

private:

    CMarketDataSP       market;

    ExpiryArraySP       expiries;
    CDoubleArraySP      strikes;
    bool                isAbsolute;

    string              asset;
    string              yield;

    string              volType;

    static IObjectSP check(VolGridChecker* params) {
        static const string routine = "VolGridChecker::prices";
        try {
            // Construct the MarketWrappers
            CAssetWrapper       assetWrapper(params->asset);
            YieldCurveWrapper   yieldWrapper(params->yield);

            // Construct the maturity dates from the expiry array
            DateTime refDate( params->market->GetReferenceDate() );
            ExpiryArraySP expiries = params->expiries;
            DateTimeArraySP maturities( new DateTimeArray(expiries->size()) );
            for( int i=0; i<expiries->size(); i++ ) {
                (*maturities)[i] = (*expiries)[i]->toDate(refDate);
            }

            // The arbitrage checking mechanism expects absolute strikes.
            if( !params->isAbsolute ) {
                // If necessary, scale the percentage strikes by the spot.
                // Construct tempory model and asset to retrieve spot details
                CAssetWrapper tmpAsset(params->asset);
                CClosedFormLN model(params->volType,
                                    false);
                tmpAsset.getData(&model, params->market);

                double spot = tmpAsset->getSpot();
                for(int iStrike = 0; iStrike<params->strikes->size(); iStrike++) {
                    (*params->strikes)[iStrike] *= spot;
                }
            }

            // Construct a new VanillaGridAddin
            VanillaGridAddinSP vgridParams( new VanillaGridAddin(
                params->market,
                maturities,
                params->strikes,
                assetWrapper,
                yieldWrapper,
                params->volType) );

            IObjectSP arbitrageObj = VanillaGridAddin::prices( vgridParams.get() );

            // Return the ArbitrageChecks object
            return arbitrageObj;
        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }


    /** for reflection */
    VolGridChecker():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(VolGridChecker, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultVolGridChecker);
        FIELD(market, "market");
        FIELD(expiries, "Array of expiries");
        FIELD(strikes, "Array of strikes");
        FIELD(isAbsolute,"Strikes are absolute, not percentages");
        FIELD(asset, "Name of Asset");
        FIELD(yield, "Name of Yield curve");
        FIELD(volType, "Type of vol to use");
    }

    static IObject* defaultVolGridChecker(){
        return new VolGridChecker();
    }

};

CClassConstSP const VolGridChecker::TYPE = CClass::registerClassLoadMethod(
    "VolGridChecker", typeid(VolGridChecker), load);

bool VolCheckLoad() {
    return (VanillaGridAddin::TYPE != 0
            && VolGridChecker::TYPE != 0
            && ArbitrageChecks::TYPE != 0);
}

DRLIB_END_NAMESPACE

