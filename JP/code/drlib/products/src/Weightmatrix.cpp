//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : WeightMatrix.cpp
//
//   Description : Weight Matrix
//
//   Author      : Regis Guichard
//
//   Date        : 03 Jan 03
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/NonPricingModel.hpp"
#include "edginc/CashSettlePeriod.hpp"
#include "edginc/VanillaGrid.hpp"
#include "edginc/DeltaToStrike.hpp"
#include "edginc/Scenario.hpp"
#include "edginc/Delta.hpp"
#include "edginc/BenchmarkDate.hpp"
#include "edginc/Maths.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/VolProcessedBS.hpp"

DRLIB_BEGIN_NAMESPACE

// DEFAULT REF EXPIRIES

/** Invoked when Class is 'loaded' */
void DefaultWeightMatrix::DefaultRefExpiries::load(CClassSP& clazz){
    REGISTER(DefaultRefExpiries, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultCtor);

    // registration for addin function
    Addin::registerClassObjectMethod("WEIGHT_MATRIX_REF_EXPIRIES_GET_DEFAULT",
                                     Addin::RISK,
                                     "Returns the weight matrix default reference expiries",
                                     DefaultRefExpiries::TYPE,
                                     false,
                                     Addin::expandMulti,
                                     (Addin::ObjMethod*)run);
}

ExpiryArray DefaultWeightMatrix::DefaultRefExpiries::create(){
    ExpiryArray refExpiries(4);
    refExpiries[0] = ExpirySP(new MaturityPeriod("1D"));
    refExpiries[1] = ExpirySP(new MaturityPeriod("3M"));
    refExpiries[2] = ExpirySP(new MaturityPeriod("1Y"));
    refExpiries[3] = ExpirySP(new MaturityPeriod("7Y"));
    return refExpiries;
}

IObjectSP DefaultWeightMatrix::DefaultRefExpiries::run(DefaultRefExpiries* params){
    ExpiryArray expiries(create());
    int n = expiries.size();
    StringArray rtn(n);
    int i = 0;
    for (; i < n; ++i){
        rtn[i] = expiries[i]->toString();
    }
    return IObjectSP(rtn.clone());
}

DefaultWeightMatrix::DefaultRefExpiries::DefaultRefExpiries():
CObject(TYPE){}

IObject* DefaultWeightMatrix::DefaultRefExpiries::defaultCtor(){
    return new DefaultRefExpiries();
}

CClassConstSP const DefaultWeightMatrix::DefaultRefExpiries::TYPE = CClass::registerClassLoadMethod(
    "DefaultWeightMatrix::DefaultRefExpiries", typeid(DefaultWeightMatrix::DefaultRefExpiries), load);

// DEFAULT REF LOWER STRIKES

/** Invoked when Class is 'loaded' */
void DefaultWeightMatrix::DefaultRefLowerStrikes::load(CClassSP& clazz){
    REGISTER(DefaultRefLowerStrikes, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultCtor);
    FIELD(strikeUnits, "strikeUnits");
    FIELD_MAKE_OPTIONAL(strikeUnits);
    FIELD(instType, "instType");
    FIELD_MAKE_OPTIONAL(instType);
    FIELD(userMode, "userMode");
    FIELD_MAKE_OPTIONAL(userMode);

    // registration for addin function
    Addin::registerClassObjectMethod("WEIGHT_MATRIX_REF_LOWER_STRIKES_GET_DEFAULT",
                                     Addin::RISK,
                                     "Returns the weight matrix default reference lower strikes",
                                     DefaultRefLowerStrikes::TYPE,
                                     false,
                                     Addin::expandMulti,
                                     (Addin::ObjMethod*)run);
}

DoubleArray DefaultWeightMatrix::DefaultRefLowerStrikes::create(const string& strikeUnits,
                                                                const string& instType,
                                                                const string& userMode){
    DoubleArray refStrikes(4);
    // backward compatibility
    if (CString::equalsIgnoreCase(userMode,DefaultWeightMatrix::OLDMODE)){
        refStrikes[0] = 0.95;
        refStrikes[1] = 0.80;
        refStrikes[2] = 0.50;
        refStrikes[3] = 0.95;
    } else if (CString::equalsIgnoreCase(userMode,DefaultWeightMatrix::NEWMODE)) {

        // new fonctinnality
        if (CString::equalsIgnoreCase(strikeUnits,WeightMatrix::ABSOLUTE)){
            refStrikes[0] = 0.0;
            refStrikes[1] = 0.0;
            refStrikes[2] = 0.0;
            refStrikes[3] = 0.0;
        }
        if (CString::equalsIgnoreCase(strikeUnits,WeightMatrix::SPOTMONEYNESS)
            ||CString::equalsIgnoreCase(strikeUnits,WeightMatrix::FWDMONEYNESS)){
            refStrikes[0] = 0.95;
            refStrikes[1] = 0.80;
            refStrikes[2] = 0.50;
            refStrikes[3] = 0.95;
        }

        if (CString::equalsIgnoreCase(strikeUnits,WeightMatrix::DELTA)){
           if (CString::equalsIgnoreCase(instType,WeightMatrix::CALL)){
               refStrikes = DoubleArray(4,0.95);
            } 
           if (CString::equalsIgnoreCase(instType,WeightMatrix::PUT)){
                refStrikes = DoubleArray(4,-0.05);
            }
            if (CString::equalsIgnoreCase(instType,WeightMatrix::OTM)){
                refStrikes = DoubleArray(4,-0.05);
            } 
           if (CString::equalsIgnoreCase(instType,WeightMatrix::ITM)){
                refStrikes = DoubleArray(4,0.95);
            }
        }
    } else {
         throw ModelException("DefaultWeightMatrix::DefaultRefLowerStrikes::create: userMode"
             + userMode + "is not recognized");

    }
    return refStrikes;
}

IObjectSP DefaultWeightMatrix::DefaultRefLowerStrikes::run(DefaultRefLowerStrikes* params){
    return IObjectSP(create(params->strikeUnits,params->instType,params->userMode).clone());
}

DefaultWeightMatrix::DefaultRefLowerStrikes::DefaultRefLowerStrikes():
CObject(TYPE),
strikeUnits(WeightMatrix::DEFAULT),
instType(WeightMatrix::DEFAULT),
userMode(DefaultWeightMatrix::OLDMODE){}

IObject* DefaultWeightMatrix::DefaultRefLowerStrikes::defaultCtor(){
    return new DefaultRefLowerStrikes();
}

CClassConstSP const DefaultWeightMatrix::DefaultRefLowerStrikes::TYPE = CClass::registerClassLoadMethod(
    "DefaultWeightMatrix::DefaultRefLowerStrikes", typeid(DefaultWeightMatrix::DefaultRefLowerStrikes), load);

// DEFAULT REF UPPER STRIKES

/** Invoked when Class is 'loaded' */
void DefaultWeightMatrix::DefaultRefUpperStrikes::load(CClassSP& clazz){
    REGISTER(DefaultRefUpperStrikes, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultCtor);
    FIELD(strikeUnits, "strikeUnits");
    FIELD_MAKE_OPTIONAL(strikeUnits);
    FIELD(instType, "instType");
    FIELD_MAKE_OPTIONAL(instType);
    FIELD(userMode, "userMode");
    FIELD_MAKE_OPTIONAL(userMode);

    // registration for addin function
    Addin::registerClassObjectMethod("WEIGHT_MATRIX_REF_UPPER_STRIKES_GET_DEFAULT",
                                     Addin::RISK,
                                     "Returns the weight matrix default reference upper strikes",
                                     DefaultRefUpperStrikes::TYPE,
                                     false,
                                     Addin::expandMulti,
                                     (Addin::ObjMethod*)run);
}

DoubleArray DefaultWeightMatrix::DefaultRefUpperStrikes::create(const string& strikeUnits,
                                                                const string& instType,
                                                                const string& userMode){
    DoubleArray refStrikes(4);
    // backward compatibility
    if (CString::equalsIgnoreCase(userMode,DefaultWeightMatrix::OLDMODE)){
        refStrikes[0] = 1.05;
        refStrikes[1] = 1.20;
        refStrikes[2] = 1.50;
        refStrikes[3] = 1.05;
    } else if ((CString::equalsIgnoreCase(userMode,DefaultWeightMatrix::NEWMODE))) {
        // new fonctionnality
        if (CString::equalsIgnoreCase(strikeUnits,WeightMatrix::ABSOLUTE)){
            refStrikes[0] = 0.0;
            refStrikes[1] = 0.0;
            refStrikes[2] = 0.0;
            refStrikes[3] = 0.0;
        }
        if (CString::equalsIgnoreCase(strikeUnits,WeightMatrix::SPOTMONEYNESS)
            ||CString::equalsIgnoreCase(strikeUnits,WeightMatrix::FWDMONEYNESS)){
            refStrikes[0] = 1.05;
            refStrikes[1] = 1.20;
            refStrikes[2] = 1.50;
            refStrikes[3] = 1.05;
        }
        if (CString::equalsIgnoreCase(strikeUnits,WeightMatrix::DELTA)){
           if (CString::equalsIgnoreCase(instType,WeightMatrix::CALL)){
                refStrikes = DoubleArray(4,0.05);
            } 
           if (CString::equalsIgnoreCase(instType,WeightMatrix::PUT)){
               refStrikes = DoubleArray(4,-0.95);
            }
            if (CString::equalsIgnoreCase(instType,WeightMatrix::OTM)){
                refStrikes = DoubleArray(4,0.05);
            } 
           if (CString::equalsIgnoreCase(instType,WeightMatrix::ITM)){
               refStrikes = DoubleArray(4,0.95);
            }
        }
    } else {
         throw ModelException("DefaultWeightMatrix::DefaultRefLowerStrikes::create: userMode"
             + userMode + "is not recognized");
    }
    return refStrikes;
}

IObjectSP DefaultWeightMatrix::DefaultRefUpperStrikes::run(DefaultRefUpperStrikes* params){
    return IObjectSP(create(params->strikeUnits,params->instType,params->userMode).clone());
}

DefaultWeightMatrix::DefaultRefUpperStrikes::DefaultRefUpperStrikes():
CObject(TYPE),
strikeUnits(WeightMatrix::DEFAULT),
instType(WeightMatrix::DEFAULT),
userMode(DefaultWeightMatrix::OLDMODE){}

IObject* DefaultWeightMatrix::DefaultRefUpperStrikes::defaultCtor(){
    return new DefaultRefUpperStrikes();
}

CClassConstSP const DefaultWeightMatrix::DefaultRefUpperStrikes::TYPE = CClass::registerClassLoadMethod(
    "DefaultWeightMatrix::DefaultRefUpperStrikes", typeid(DefaultWeightMatrix::DefaultRefUpperStrikes), load);

// WEIGHT MATRIX
CDoubleMatrixSP DefaultWeightMatrix::create(const DateTime&      baseDate,
                                            double               spot,
                                            const DateTimeArray& maturities,
                                            const DoubleArray&   strikes,
                                            const string&        strikeUnits,
                                            const string&        instType,
                                            const string&        userMode){
    return create(baseDate,
                  spot,
                  maturities,
                  strikes,
                  DefaultRefExpiries::create(),                  
                  DefaultRefLowerStrikes::create(strikeUnits,instType,userMode),
                  DefaultRefUpperStrikes::create(strikeUnits,instType,userMode),
                  false,
                  DoubleArray(0),
                  strikeUnits,
                  instType,
                  userMode);
}

CDoubleMatrixSP DefaultWeightMatrix::create(const DateTime&      baseDate,
                                            double               spot,
                                            const DateTimeArray& maturities,
                                            const DoubleArray&   strikes,
                                            const ExpiryArray&   refExpiries,
                                            const DoubleArray&   refLowerStrikes,
                                            const DoubleArray&   refUpperStrikes,
                                            bool                 useFwd,
                                            const DoubleArray    fwds,
                                            const string&        strikeUnits,
                                            const string&        instType,
                                            const string&        userMode){
    const string method = "DefaultWeightMatrix::create";
    try{
        int nbMats = maturities.size();
        if (useFwd && (fwds.size() != nbMats)){
        throw ModelException(method,
                             "Nber of maturities and nber of fwds must be the same; got "
                             + Format::toString(nbMats) + " and "
                             + Format::toString(fwds.size()) + ", respectively");
        }

        int nbStrikes = strikes.size();
        CDoubleMatrixSP weights(new DoubleMatrix(nbMats, nbStrikes));

        // create ref date array and check is increasing
        int nbRefDates = refExpiries.size();
        DateTimeArray refDates(nbRefDates);
        int iRefDate = 0;
        DateTime lastDate = baseDate;
        for (; iRefDate < nbRefDates; ++iRefDate){
            refDates[iRefDate] = refExpiries[iRefDate]->toDate(baseDate);
            if (!refDates[iRefDate].isGreater(lastDate)){
                throw ModelException(method,
                                     "refExpiries must be increasing");
            }
            lastDate = refDates[iRefDate];
        }

        // loop over maturities
        iRefDate = 0;
        double lastYearFrac = 0.0;
        double yearFracDiff = baseDate.yearFrac(refDates[0]) - lastYearFrac;
        double lastLowerStrike = refLowerStrikes[0];
        double lowerStrikeDiff = refLowerStrikes[0] - lastLowerStrike;
        double lastUpperStrike = refUpperStrikes[0];
        double upperStrikeDiff = refUpperStrikes[0] - lastUpperStrike;
        int iMat = 0;
        for (; iMat < nbMats; ++iMat){
            // locate which maturity period we are in
            while (iRefDate < nbRefDates && refDates[iRefDate].isLess(maturities[iMat])){
                lastYearFrac = baseDate.yearFrac(refDates[iRefDate]);
                lastLowerStrike = refLowerStrikes[iRefDate];
                lastUpperStrike = refUpperStrikes[iRefDate];
                ++iRefDate;
                if (iRefDate < nbRefDates){
                    yearFracDiff = baseDate.yearFrac(refDates[iRefDate]) - lastYearFrac;
                    lowerStrikeDiff = refLowerStrikes[iRefDate] - lastLowerStrike;
                    upperStrikeDiff = refUpperStrikes[iRefDate] - lastUpperStrike;
                }
            }
            // calculate current boundary
            double yearFrac = baseDate.yearFrac(maturities[iMat]);
            double lowerStrike;
            double upperStrike;
            if (iRefDate < nbRefDates){
                lowerStrike = lastLowerStrike + lowerStrikeDiff * (yearFrac - lastYearFrac) / yearFracDiff;
                upperStrike = lastUpperStrike + upperStrikeDiff * (yearFrac - lastYearFrac) / yearFracDiff;
            }
            else{
                lowerStrike = lastLowerStrike;
                upperStrike = lastUpperStrike;
            }
            // fill in weights at current mat
            double ref = (useFwd ? fwds[iMat] : spot);
            int iStrike = 0;
            for (; iStrike < nbStrikes; ++iStrike){
                // backward compatibility
                if (CString::equalsIgnoreCase(userMode,DefaultWeightMatrix::OLDMODE)){
                    double moneyness = strikes[iStrike] / ref;
                    if (Maths::isNegative(moneyness - lowerStrike)){
                        (*weights)[iMat][iStrike] = 0.0;
                    }
                    else if (Maths::isPositive(moneyness - upperStrike)){
                        (*weights)[iMat][iStrike] = 0.0;
                    }
                    else{
                        (*weights)[iMat][iStrike] = 1.0;
                    }
                } else if (CString::equalsIgnoreCase(userMode,DefaultWeightMatrix::NEWMODE)){
                    // new interface
                    if (!CString::equalsIgnoreCase(strikeUnits,WeightMatrix::DELTA)){
                        // if strikeunits is not Delta
                        if (Maths::isNegative(strikes[iStrike] - lowerStrike)){
                            (*weights)[iMat][iStrike] = 0.0;
                        }
                        else if (Maths::isPositive(strikes[iStrike] - upperStrike)){
                            (*weights)[iMat][iStrike] = 0.0;
                        }
                        else{
                            (*weights)[iMat][iStrike] = 1.0;
                        }
                    } else {
                        // if strikeUnits is Delta
                        if (CString::equalsIgnoreCase(instType,WeightMatrix::CALL)
                            ||CString::equalsIgnoreCase(instType,WeightMatrix::PUT)){
                            if (Maths::isPositive(strikes[iStrike] - lowerStrike)){
                                (*weights)[iMat][iStrike] = 0.0;
                            }
                            else if (Maths::isNegative(strikes[iStrike] - upperStrike)){
                                (*weights)[iMat][iStrike] = 0.0;
                            }
                            else{
                                (*weights)[iMat][iStrike] = 1.0;
                            }
                        }
                        if (CString::equalsIgnoreCase(instType,WeightMatrix::OTM)){
                            if (Maths::isPositive(strikes[iStrike] - lowerStrike)
                                &&Maths::isNegative(strikes[iStrike])){
                                (*weights)[iMat][iStrike] = 0.0;
                            }
                            else if (Maths::isNegative(strikes[iStrike] - upperStrike)
                                    &&Maths::isPositive(strikes[iStrike])){
                                (*weights)[iMat][iStrike] = 0.0;
                            }
                            else{
                                (*weights)[iMat][iStrike] = 1.0;
                            }
                        } 
                        if (CString::equalsIgnoreCase(instType,WeightMatrix::ITM)){
                            if (Maths::isPositive(strikes[iStrike] - lowerStrike)
                                &&Maths::isPositive(strikes[iStrike])){
                                (*weights)[iMat][iStrike] = 0.0;
                            }
                            else if (Maths::isNegative(strikes[iStrike] - upperStrike)
                                    &&Maths::isNegative(strikes[iStrike])){
                                (*weights)[iMat][iStrike] = 0.0;
                            }
                            else{
                                (*weights)[iMat][iStrike] = 1.0;
                            }
                        }
                    }
                } else {
                    throw ModelException("DefaultWeightMatrix::create: userMode"
                        + userMode + "is not recognized");
                }
            }
        }
        return weights;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

void DefaultWeightMatrix::validatePop2Object(){
    const string method("DefaultWeightMatrix::validatePop2Object");
    try{
        // validate sizes
        if ((refExpiries.size() != refLowerStrikes.size())
            || (refExpiries.size() != refUpperStrikes.size())){
            throw ModelException(method,
                                 "Nber of refExpiries must be the same as nber of refLowerStrikes "
                                 "and nber of refUpperStrikes; got "
                                 + Format::toString(refExpiries.size()) + ", "
                                 + Format::toString(refLowerStrikes.size()) + " and "
                                 + Format::toString(refUpperStrikes.size()) + ", respectively");
        }
        // get date
        market->GetReferenceDate(baseDate);
        // create a non pricing model in order to get the market data
        MarketDataFetcherSP mdf(new MarketDataFetcherLN("VolSurface"));
        NonPricingModel dummyModel(mdf);
        asset.getData(&dummyModel, market.get());
        // calc spot or fwd depending on what's requested
        if (useFwd){
            fwds.resize(maturities.size());
            asset->fwdValue(maturities,
                            fwds);
        }
        else{
            spot = asset->getSpot();
        }

        // validate userMode
        if(!CString::equalsIgnoreCase(userMode,DefaultWeightMatrix::OLDMODE)
            &&!CString::equalsIgnoreCase(userMode,DefaultWeightMatrix::NEWMODE)){
            throw ModelException(method,"userMode should be either"
                + DefaultWeightMatrix::OLDMODE
                + DefaultWeightMatrix::NEWMODE);
        }
        
        // construct a weight matrix for validation
        DateTimeArraySP maturitiesSP(copy(&maturities));
        DoubleArraySP strikesSP(copy(&strikes));
        CDoubleMatrixSP weightsSP(new CDoubleMatrix(maturities.size(),strikes.size()));
        
        WeightMatrix wm("", 
                        maturitiesSP,
                        strikesSP,
                        weightsSP,
                        strikeUnits,
                        instType,
                        baseDate);
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/** Invoked when Class is 'loaded' */
void DefaultWeightMatrix::load(CClassSP& clazz){
    REGISTER(DefaultWeightMatrix, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultCtor);
    FIELD(asset, "asset");
    FIELD(market, "market");        
    FIELD(maturities, "Maturities");
    FIELD(strikes, "Strikes");
    FIELD(refExpiries, "refExpiries");
    FIELD_MAKE_OPTIONAL(refExpiries);
    FIELD(refLowerStrikes, "refLowerStrikes");
    FIELD_MAKE_OPTIONAL(refLowerStrikes);
    FIELD(refUpperStrikes, "refUpperStrikes");
    FIELD_MAKE_OPTIONAL(refUpperStrikes);
    FIELD(useFwd, "If true, centered at fwd; centerd a spot, otherwise");
    FIELD_MAKE_OPTIONAL(useFwd);
    FIELD(isVegaWeighted, "If true, non-zero weight are vega-weighted");
    FIELD_MAKE_OPTIONAL(isVegaWeighted);
    FIELD(yield, "yield");
    FIELD_MAKE_OPTIONAL(yield);
    FIELD(strikeUnits, "strikeUnits");
    FIELD_MAKE_OPTIONAL(strikeUnits);
    FIELD(instType, "instType");
    FIELD_MAKE_OPTIONAL(instType);
    FIELD(userMode, "userMode");
    FIELD_MAKE_OPTIONAL(userMode);
    
    // transient
    FIELD(baseDate, "");
    FIELD_MAKE_TRANSIENT(baseDate)
    FIELD(spot, "");
    FIELD_MAKE_TRANSIENT(spot)
    FIELD(fwds, "");
    FIELD_MAKE_TRANSIENT(fwds);


    // registration for addin function
    Addin::registerClassObjectMethod("WEIGHT_MATRIX_GET_DEFAULT",
                                     Addin::RISK,
                                     "Returns a matrix of default weights",
                                     DefaultWeightMatrix::TYPE,
                                     false,
                                     Addin::expandMulti,
                                     (Addin::ObjMethod*)run);
}

IObjectSP DefaultWeightMatrix::run(DefaultWeightMatrix* params){
    const string method("DefaultWeightMatrix::run");
    try{
        CDoubleMatrixSP weights(create(params->baseDate,
                                       params->spot,
                                       params->maturities,
                                       params->strikes,
                                       params->refExpiries,
                                       params->refLowerStrikes,
                                       params->refUpperStrikes,
                                       params->useFwd,
                                       params->fwds,
                                       params->strikeUnits,
                                       params->instType,
                                       params->userMode));
        if (params->isVegaWeighted
            &&(params->userMode==DefaultWeightMatrix::OLDMODE)){
            // model
            CClosedFormLN model("VolSurface");
            // settlement
            InstrumentSettlementSP instSettle(new CashSettlePeriod(0));
            // instrument

            // create a weightmatrix for the vanillagrid            
            DateTimeArraySP maturitiesSP(copy(&(params->maturities)));
            DoubleArraySP strikesSP(copy(&(params->strikes)));
            CDoubleMatrixSP weightsSP(copy(&(*weights)));
            WeightMatrixSP weightMatrix(new WeightMatrix("", 
                                                        maturitiesSP,
                                                        strikesSP,
                                                        weightsSP,
                                                        WeightMatrix::ABSOLUTE,
                                                        WeightMatrix::CALL,
                                                        params->baseDate));
            WeightMatrixWrapper weightMatrixWrapper(weightMatrix);

            VanillaGrid inst(params->asset,
                             params->yield,
                             weightMatrixWrapper,
                             instSettle);
            // control
            OutputRequestArraySP outputArray(new OutputRequestArray(1));
            (*outputArray)[0] = OutputRequestSP(new OutputRequest(OutputRequest::OPTION_VEGA));
            Control control(SensitivityArrayConstSP(new SensitivityArray(0)),
                            outputArray,
                            false,
                            string(""));
            // run model
            CResultsSP results(model.go(CInstrumentSP::attachToRef(&inst),
                                        ScenarioSP(),
                                        CControlSP::attachToRef(&control),
                                        params->market));
            // read off result
            IObjectConstSP vegasobj(results->retrieveRequestResult(OutputRequest::OPTION_VEGA));
            VanillaGrid::OutputConstSP vegas(VanillaGrid::OutputConstSP::dynamicCast(vegasobj));
            // calculate maximum vega
            double maxvega = 0.0;
            int iCol = 0;
            for (; iCol < vegas->numCols(); ++iCol){
                int iRow = 0;
                for (; iRow < vegas->numRows(); ++iRow){
                    if (Maths::isPositive((*weights)[iCol][iRow])){
                        maxvega = Maths::max(vegas->getValue(iCol, iRow), maxvega);
                    }
                }
            }
            if (!Maths::isPositive(maxvega)){
                maxvega = 1.0;
            }
            // multiply weights by vegas                
            for (iCol = 0; iCol < weights->numCols(); ++iCol){
                int iRow = 0;
                for (; iRow < weights->numRows(); ++iRow){
                    if (Maths::isPositive((*weights)[iCol][iRow])){
                        (*weights)[iCol][iRow] *= (vegas->getValue(iCol, iRow) / maxvega);
                    }
                }
            }
        }
        return weights;
    }
    catch(ModelException& e){
        throw ModelException(e, method);
    }
}

DefaultWeightMatrix::DefaultWeightMatrix():
CObject(TYPE),
refExpiries(DefaultRefExpiries::create()),
refLowerStrikes(DefaultRefLowerStrikes::create(WeightMatrix::DEFAULT,WeightMatrix::DEFAULT,DefaultWeightMatrix::OLDMODE)),
refUpperStrikes(DefaultRefUpperStrikes::create(WeightMatrix::DEFAULT,WeightMatrix::DEFAULT,DefaultWeightMatrix::OLDMODE)),
useFwd(false),
isVegaWeighted(false),
spot(0.0),
strikeUnits(WeightMatrix::DEFAULT),
instType(WeightMatrix::DEFAULT),
userMode(DefaultWeightMatrix::OLDMODE){}

IObject* DefaultWeightMatrix::defaultCtor(){
    return new DefaultWeightMatrix();
}
    
CClassConstSP const DefaultWeightMatrix::TYPE = CClass::registerClassLoadMethod(
    "DefaultWeightMatrix", typeid(DefaultWeightMatrix), load);


const string DefaultWeightMatrix::OLDMODE = "OldMode";
const string DefaultWeightMatrix::NEWMODE = "NewMode";



// --------------------------------------------------------------------
// WeightMatrix

WeightMatrix::WeightMatrix(const string&            name, 
                           DateTimeArraySP          dateTimeMaturities,
                           DoubleArraySP            strikes,
                           CDoubleMatrixSP          weights,
                           const string&            strikeUnits,
                           const string&            instType,
                           const DateTime&          refDate) : 
    MarketObject(TYPE),
    name(name), dateTimeMaturities(dateTimeMaturities),
    strikes(strikes), weights(weights),
    strikeUnits(strikeUnits), instType(instType), refDate(refDate) {
    
    static const string method = "WeightMatrix::WeightMatrix";
    try{

        if (!dateTimeMaturities){
            throw ModelException("Internal error: dateTimeMaturities is NULL");
        }
    
        maturities = ExpiryArraySP(new ExpiryArray((*dateTimeMaturities).size()));
        int iMat;
        for (iMat=0; iMat<(*dateTimeMaturities).size(); iMat++){
            ExpirySP singleExp(new BenchmarkDate((*dateTimeMaturities)[iMat]));
            (*maturities)[iMat] = singleExp;
        }
        validatePop2Object();
    }
    catch(exception& e){
        throw ModelException(e, method);
    }   
}

WeightMatrix::~WeightMatrix(){}

string WeightMatrix::getName() const {
    return name;
}
// Assuming this is not used in a performance-critical place...
DateTimeArraySP WeightMatrix::getDateTimeMaturities() const {
    return dateTimeMaturities;
}

DoubleArraySP WeightMatrix::getStrikes() const {
    return strikes;
}

CDoubleMatrixSP WeightMatrix::getWeights() const {
    return weights;
}

const string& WeightMatrix::getStrikeUnits() const{
    return strikeUnits;
}

const string& WeightMatrix::getInstType() const{
    return instType;
}

CDoubleMatrixSP WeightMatrix::getInstStrikesUsed() const{
    static const string method = "WeightMatrix::getInstTypesUsed";
    try{
        if(!instStrikesUsed) {
            throw ModelException("Internal error: instStrikesUsed is NULL. Client must call computeStrikesAndTypes method");
        }
        return instStrikesUsed;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }   
}

IntArrayArraySP WeightMatrix::getInstTypesUsed() const{
    static const string method = "WeightMatrix::getInstTypesUsed";
    try{
        if(!instTypesUsed) {
            throw ModelException("Internal error: instTypesUsed is NULL. Client must call computeStrikesAndTypes method");
        }
        return instTypesUsed;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }   
}

void WeightMatrix::setWeights(CDoubleMatrixSP weightsSP) {
    weights = weightsSP;
}

void WeightMatrix::validatePop2Object() {
    
    static const string method = "WeightMatrix::validatePop2Object";
    try{
      
        // check that we have at least one maturity
        int nbMats = maturities->size();
        if (nbMats < 1){
            throw ModelException(method, "maturities is empty");

        }
        
        // check that we have at least one strike
        int nbStrikes = strikes->size();
        if (nbStrikes < 1){
            throw ModelException(method, "strikes is empty");
        }

        // check strikes are increasing
        if (!CString::equalsIgnoreCase(strikeUnits,WeightMatrix::DELTA)){
            int iStrike = 1;
            for (; iStrike < nbStrikes; ++iStrike){
                if (!Maths::isPositive((*strikes)[iStrike] - (*strikes)[iStrike-1])){
                    throw ModelException(method, 
                                         Format::toString("%d-th strike must be greater than %d-th strike",
                                                          iStrike+1,
                                                          iStrike));
                }
            }
        }
        
        // check that there are weights
        if (!weights) 
        {
            throw ModelException(method, "weights is empty");

        } else {
            int nbCols = weights->numCols();
            if (nbCols != nbMats){
                throw ModelException(method,
                                     Format::toString("Nb of columns (%d) in weights should be equal to nb of maturities (%d)",
                                                      nbCols,
                                                      nbMats));
            }

            int nbRows = weights->numRows();
            if (nbRows != nbStrikes){
                throw ModelException(method,
                                     Format::toString("Nb of rows (%d) in weights should be equal to nb of strikes (%d)",
                                                      nbRows,
                                                      nbStrikes));
            }
            weights->checkNonNegative();
        }

        // validate strikeUnits input
        if (!CString::equalsIgnoreCase(strikeUnits,WeightMatrix::ABSOLUTE)
            &&!CString::equalsIgnoreCase(strikeUnits,WeightMatrix::SPOTMONEYNESS)
            &&!CString::equalsIgnoreCase(strikeUnits,WeightMatrix::FWDMONEYNESS)
            &&!CString::equalsIgnoreCase(strikeUnits,WeightMatrix::DELTA)
			&&!CString::equalsIgnoreCase(strikeUnits,WeightMatrix::ATMVARIANCE)
            &&!CString::equalsIgnoreCase(strikeUnits,WeightMatrix::DEFAULT)) {
            throw ModelException(method,"strikeUnits " + strikeUnits + " is not recognized");
        } else {
            if (CString::equalsIgnoreCase(strikeUnits,WeightMatrix::DEFAULT)) {
                strikeUnits = WeightMatrix::ABSOLUTE;
            }
        }
    
        // validate instType input
        if (!CString::equalsIgnoreCase(instType,WeightMatrix::CALL)
            &&!CString::equalsIgnoreCase(instType,WeightMatrix::PUT)
            &&!CString::equalsIgnoreCase(instType,WeightMatrix::OTM)
            &&!CString::equalsIgnoreCase(instType,WeightMatrix::ITM)
            &&!CString::equalsIgnoreCase(instType,WeightMatrix::DEFAULT)
            &&!CString::equalsIgnoreCase(instType,WeightMatrix::EMPTY)) {
            throw ModelException(method,"instType " + instType + " is not recognized");
        } else {
            if (CString::equalsIgnoreCase(instType,WeightMatrix::DEFAULT)) {
                instType = WeightMatrix::EMPTY;
            }
            if (CString::equalsIgnoreCase(instType,WeightMatrix::EMPTY)) {
                instType = WeightMatrix::CALL;
            }
        }
            
        // strikes validation when strikeUnits is delta
        if (CString::equalsIgnoreCase(strikeUnits,WeightMatrix::DELTA)) {
            
            // strikes validation when instType is Call or Put
            if ((CString::equalsIgnoreCase(instType,WeightMatrix::CALL))||(CString::equalsIgnoreCase(instType,WeightMatrix::PUT))){
                int iStrike;
                for (iStrike = 0; iStrike < nbStrikes; iStrike++){
                    if ((CString::equalsIgnoreCase(instType,WeightMatrix::CALL))&&(!Maths::isPositive((*strikes)[iStrike]))) {
                        throw ModelException(method, "Strike" + Format::toString((*strikes)[iStrike]) +"should be positive");
                    }
                    if ((CString::equalsIgnoreCase(instType,WeightMatrix::PUT))&&(!Maths::isNegative((*strikes)[iStrike]))) {
                        throw ModelException(method, "Strike" + Format::toString((*strikes)[iStrike]) +"should be negative");
                    }
                }
                for (iStrike = 1; iStrike < nbStrikes; iStrike++) {
                    if (Maths::isPositive((*strikes)[iStrike] - (*strikes)[iStrike-1])) {
                        throw ModelException(method, Format::toString("%d-th strike must be lower than %d-th strike",iStrike+1,iStrike));
                    }
                }
            }

            // strikes validation when instType is OTM
            if (CString::equalsIgnoreCase(instType,WeightMatrix::OTM)) {
                int lStrike = -1;
                int iStrike = 1;
                while ((lStrike==-1)&&(iStrike < nbStrikes)) {
                    if ((Maths::isPositive((*strikes)[iStrike])) 
                        && (Maths::isNegative((*strikes)[iStrike-1]))) {
                        lStrike = iStrike;
                    }
                    iStrike++;
                }
                if (lStrike==-1){
                   throw ModelException(method, "strikes is not consistent with instType");
                } else {
                    for (iStrike = 0; iStrike < nbStrikes; iStrike++) {
                        if (iStrike < lStrike){
                            if ((!Maths::isNegative((*strikes)[iStrike]))
                                ||(Maths::isNegative((*strikes)[iStrike]+0.5))) {
                                throw ModelException(method,"strikes is not consistent with instType");
                            }
                        }
                        if (iStrike>=lStrike) {
                            if ((!Maths::isPositive((*strikes)[iStrike]))
                                ||(!Maths::isNegative((*strikes)[iStrike]-0.5))) {
                                throw ModelException(method,"strikes is not consistent with instType");
                            }
                        }
                        // check that strikes are decreasing
                        if (iStrike>0) {
                            if ((Maths::isPositive((*strikes)[iStrike]-(*strikes)[iStrike-1]))
                                &&(Maths::isNegative((*strikes)[iStrike]))
                                &&(Maths::isNegative((*strikes)[iStrike-1]))) {
                                throw ModelException(method,"strikes is not consistent with instType");
                            }
                            if ((Maths::isPositive((*strikes)[iStrike]-(*strikes)[iStrike-1]))
                                &&(Maths::isPositive((*strikes)[iStrike]))
                                &&(Maths::isPositive((*strikes)[iStrike-1]))) {
                                throw ModelException(method,"strikes is not consistent with instType");
                            }
                        }
                    }
                }
            }

            // strikes validation when instType is ITM
            if (CString::equalsIgnoreCase(instType,WeightMatrix::ITM)) {
                int lStrike = -1;
                int iStrike = 1;
                while ((lStrike==-1)&&(iStrike < nbStrikes)) {
                    if ((Maths::isNegative((*strikes)[iStrike]))
                        && (Maths::isPositive((*strikes)[iStrike-1]))) {
                        lStrike = iStrike;
                    }
                    iStrike++;
                }
                if (lStrike==-1){
                   throw ModelException(method, "strikes is not consistent with instType");
                } else {
                    for (iStrike = 0; iStrike < nbStrikes; iStrike++) {
                        if (iStrike < lStrike) {
                            if ((!Maths::isPositive((*strikes)[iStrike]))
                                ||(Maths::isNegative((*strikes)[iStrike]-0.5))) {
                                throw ModelException(method,"strikes is not consistent with instType");
                            }
                        }
                        if (iStrike>=lStrike) {
                            if ((!Maths::isNegative((*strikes)[iStrike]))
                                ||(!Maths::isNegative((*strikes)[iStrike]+0.5))) {
                                throw ModelException(method,"strikes is not consistent with instType");
                            }
                        }
                        // check that strikes are decreasing
                        if (iStrike>0) {
                            if ((Maths::isPositive((*strikes)[iStrike]-(*strikes)[iStrike-1]))
                                &&(Maths::isNegative((*strikes)[iStrike]))
                                &&(Maths::isNegative((*strikes)[iStrike-1]))) {
                                throw ModelException(method,"strikes is not consistent with instType");
                            }
                            if ((Maths::isPositive((*strikes)[iStrike]-(*strikes)[iStrike-1]))
                                &&(Maths::isPositive((*strikes)[iStrike]))
                                &&(Maths::isPositive((*strikes)[iStrike-1]))) {
                                throw ModelException(method,"strikes is not consistent with instType");
                            }
                        }
                    }
                }
            }
        }
    }
    catch(exception& e){
        throw ModelException(e, method);
    }         
}

void WeightMatrix::getMarket(const IModel* model, const MarketData* market) {
    static const string method("WeightMatrix::GetMarket");
    try{
        // get reference date from the market
        market->GetReferenceDate(refDate);

        // compute maturities from expiries
        dateTimeMaturities = DateTimeArraySP(new DateTimeArray);
        for(int iMat=0; iMat<maturities->size(); iMat++) {
        dateTimeMaturities->push_back((*maturities)[iMat]->toDate(refDate));
        if (iMat>0 &&
            (*dateTimeMaturities)[iMat]<=(*dateTimeMaturities)[iMat-1]) {
            throw ModelException("WeightMatrix::getMaturities",
                                 "Maturities are not strictly increasing : " +
                                 (*maturities)[iMat]->toString() + " is not after " +
                                 (*maturities)[iMat-1]->toString());
            }
        }
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

void WeightMatrix::computeStrikesAndTypes(const CAssetWrapper&            asset,
                                          const YieldCurveWrapper&        discount,
                                          InstrumentSettlementSP          instSettle) {
    
    static const string method("WeightMatrix::ComputeStrikesAndTypes");
    try
    {    
        int iMat,iStrike;
		double spot;
		if ((!CString::equalsIgnoreCase(strikeUnits,WeightMatrix::ABSOLUTE))
			||(!CString::equalsIgnoreCase(instType,WeightMatrix::CALL))){
			spot = asset->getSpot();
		}

        instStrikesUsed = CDoubleMatrixSP(new CDoubleMatrix(dateTimeMaturities->size(),strikes->size()));
        instTypesUsed = IntArrayArraySP(new IntArrayArray(dateTimeMaturities->size()));
        
        for (iMat=0; iMat<dateTimeMaturities->size();iMat++){
            (*instTypesUsed)[iMat] = IntArraySP(new IntArray(strikes->size()));
            double fwdSpot;
			if (CString::equalsIgnoreCase(strikeUnits,WeightMatrix::FWDMONEYNESS)){			
				fwdSpot = asset->fwdValue((*dateTimeMaturities)[iMat]);
			}
            for (iStrike=0; iStrike<strikes->size();iStrike++){
                // computation of strikes used
                if (CString::equalsIgnoreCase(strikeUnits,WeightMatrix::ABSOLUTE)) {
                    (*instStrikesUsed)[iMat][iStrike] = (*strikes)[iStrike];
                } else if (CString::equalsIgnoreCase(strikeUnits,WeightMatrix::SPOTMONEYNESS)) {
                    (*instStrikesUsed)[iMat][iStrike] = (*strikes)[iStrike]*spot;
                } else if (CString::equalsIgnoreCase(strikeUnits,WeightMatrix::FWDMONEYNESS)) {
                    (*instStrikesUsed)[iMat][iStrike] = (*strikes)[iStrike]*fwdSpot;
                } else if (CString::equalsIgnoreCase(strikeUnits,WeightMatrix::DELTA)){
                    try {
                        // computation of the strike corresponding to the inputed delta
                        string className("CVanilla::DeltaImpliedStrikeMakerLN");
                        IDeltaToStrikeMakerSP maker(dynamic_cast<IDeltaToStrikeMaker*>(CClass::forName(className)->newInstance()));
                        if(!maker) {
                            throw ModelException(method, "Failed to create CVanilla::DeltaImpliedStrikeMakerLN");
                        }

                        double targetDelta = (*strikes)[iStrike];
                        
                        IDeltaToStrikeMaker::IDeltaToStrikeConstSP calc(maker->make(
                            refDate,
                            (*dateTimeMaturities)[iMat],
                            targetDelta < 0.0 ? false : true,
                            asset.get(),
                            discount.get(),
                            instSettle.get(),
                            Delta::DEFAULT_SHIFT,
                            targetDelta,
                            "VolPreferred",
                            true));
            
                        double modelDelta = calc->calcStrike(0.1*spot, 10.0*spot, 0.001*spot);
                        (*instStrikesUsed)[iMat][iStrike] = modelDelta;
                    } catch(exception& e){
                        throw ModelException(e, "Failed to compute delta for strike " +
                            Format::toString((*strikes)[iStrike]) + 
                            " and maturity " + 
                            (*dateTimeMaturities)[iMat].toString());
                    }
                } else if (CString::equalsIgnoreCase(strikeUnits,WeightMatrix::ATMVARIANCE)){
					ATMVolRequest req;
					CVolProcessedBSSP proc(asset->getProcessedVol(&req));
					double var = proc->CalcVar(refDate,(*dateTimeMaturities)[iMat]);
					(*instStrikesUsed)[iMat][iStrike] = spot*exp((*strikes)[iStrike]*sqrt(var));
				}

                // definition of types (Call or Put)
                if (CString::equalsIgnoreCase(instType,WeightMatrix::CALL)) {
                    (*(*instTypesUsed)[iMat])[iStrike] = 1;
                }
                if (CString::equalsIgnoreCase(instType, WeightMatrix::PUT)) {
                    (*(*instTypesUsed)[iMat])[iStrike] = 0;
                }
                if (CString::equalsIgnoreCase(instType,WeightMatrix::OTM)) {
                    if ((*instStrikesUsed)[iMat][iStrike]>=spot) {
                        (*(*instTypesUsed)[iMat])[iStrike] = 1;
                    } else {
                        (*(*instTypesUsed)[iMat])[iStrike] = 0;
                    }
                }
                if (CString::equalsIgnoreCase(instType,WeightMatrix::ITM)) {
                    if ((*instStrikesUsed)[iMat][iStrike]>=spot) {
                        (*(*instTypesUsed)[iMat])[iStrike] = 0;
                    } else {
                        (*(*instTypesUsed)[iMat])[iStrike] = 1;
                    }
                }
            }
        }    
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/** for reflection */
WeightMatrix::WeightMatrix():MarketObject(TYPE),
                                strikeUnits(DEFAULT),
                                instType(DEFAULT){}

class WeightMatrixHelper{
public:
/** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(WeightMatrix, clazz);
        SUPERCLASS(MarketObject);
        EMPTY_SHELL_METHOD(defaultWeightMatrix);
        FIELD(name, "Name identifying 'WeightMatrix's");
        FIELD(maturities, "Maturities (Expiries)");
        FIELD(strikes, "Strikes");
        FIELD(weights, "Weights");
        FIELD(strikeUnits, "strikeUnits");
        FIELD_MAKE_OPTIONAL(strikeUnits);
        FIELD(instType, "instType");
        FIELD_MAKE_OPTIONAL(instType);
        FIELD(instStrikesUsed, "instStrikesUsed");
        FIELD_MAKE_TRANSIENT(instStrikesUsed);
        FIELD(instTypesUsed, "instTypesUsed");
        FIELD_MAKE_TRANSIENT(instTypesUsed);
        FIELD(dateTimeMaturities, "dateTimeMaturities");
        FIELD_MAKE_TRANSIENT(dateTimeMaturities);
        FIELD(refDate, "RefDate");
        FIELD_MAKE_TRANSIENT(refDate);

//        ClassSetAcceptMethod(WeightMatrix::acceptWrapperNameCollector);

        Addin::registerConstructor("WEIGHT_MATRIX",
                                   Addin::MARKET,
                                   "Creates a WeightMatrix object",
                                   WeightMatrix::TYPE);
    }
    
    static IObject* defaultWeightMatrix(){
        return new WeightMatrix();
    }
};

CClassConstSP const WeightMatrix::TYPE = CClass::registerClassLoadMethod(
    "WeightMatrix", typeid(WeightMatrix), WeightMatrixHelper::load);

// definition of TYPE for MarketWrapper template class
DEFINE_TEMPLATE_TYPE(WeightMatrixWrapper);

const string WeightMatrix::DEFAULT = "Default";
const string WeightMatrix::EMPTY = "";
const string WeightMatrix::CALL = "Call";
const string WeightMatrix::PUT = "Put";
const string WeightMatrix::OTM = "OTM";
const string WeightMatrix::ITM = "ITM";
const string WeightMatrix::ABSOLUTE = "Absolute";
const string WeightMatrix::SPOTMONEYNESS = "SpotMoneyness";
const string WeightMatrix::FWDMONEYNESS = "FwdMoneyness";
const string WeightMatrix::DELTA = "Delta";
const string WeightMatrix::ATMVARIANCE = "AtmVariance";


// --------------------------------------------------------------------

DRLIB_END_NAMESPACE

