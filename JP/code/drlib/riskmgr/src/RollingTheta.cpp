//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RollingTheta.cpp
//
//   Description : Rolling Theta shift
//
//   Author      : Andre Segger
//
//   Date        : 05 Oct 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/RollingTheta.hpp"
#include "edginc/Delta.hpp"
#include "edginc/ThetaFwdSpot.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/Format.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include ext_hash_map

DRLIB_BEGIN_NAMESPACE


/** Sens Control for RollingTheta */
const string RollingTheta::NAME           = "ROLLING_VALUE";
const int    RollingTheta::DEFAULT_SHIFT  = 1;
const string RollingDeltaNAME             = "ROLLING_DELTA";


/** Is this sensitivity made using a discrete shift (ie a jump) or a
    an approximately continuous one. Returns true */
bool RollingTheta::discreteShift() const{
    return true;
}

/* silly class needed by hash_map template */
struct OutputNameHashUtil{
    bool operator()(const OutputNameConstSP& name1, 
                    const OutputNameConstSP& name2) const{
        return name1->equals(name2.get());
    }
    size_t operator()(const OutputNameConstSP& name) const{ 
        return name->hashCode();
    }
};

typedef hash_map<OutputNameConstSP, CashFlowArraySP, OutputNameHashUtil,
    OutputNameHashUtil> SensitivityData;


RollingTheta::RollingTheta(const CClassConstSP& clazz,
             const string&        outputName): 
    Sensitivity(clazz), offset(0), useAssetFwds(false){}

/** identifies the name used for storing associated results in the output*/
const string& RollingTheta::getSensOutputName() const{
    return NAME;
}

/** identifies the packet in which the results are stored. RollingTheta
    results are stored in the instrument packet */
const string& RollingTheta::getPacketName() const{
    return Results::INSTRUMENT_PACKET;
}

/** for reflection */
RollingTheta::RollingTheta(): Sensitivity(TYPE), 
                              offset(0), useAssetFwds(false) {}

void RollingTheta::validatePop2Object() {
   static const string method = "RollingTheta::validatePop2Object";
    if ( thetaInterval->size() == 0 ) {
        throw ModelException(method, 
                             "At least one pricing interval must be provided");
    }

    for (int i=1 ; i<thetaInterval->size() ; ++i) {
        const DateTime& thisStartDate = (*thetaInterval)[i].getStartDate();
        const DateTime& previousEndDate = (*thetaInterval)[i-1].getEndDate();
        if ( thisStartDate < previousEndDate ) {
            string m("Pricing intervals cannot be overlapping.\n"
                     "Start date " + Format::toString(i) + " (" + 
                     thisStartDate.toString() + 
                     ") is before end date " + Format::toString(i-1) + " (" + 
                     previousEndDate.toString() + ")");
            throw ModelException(method, m);
        }
    }
}

/** calculates given sensitivity - invoked by calculateSens */
void RollingTheta::calculate(TweakGroup*  tweakGroup,
                             CResults*    results) {

    OutputNameConstSP outputName(new OutputName(getSensOutputName()));
    bool untweakableDetected = false;

    if (!results->exists(Results::INSTRUMENT_PACKET, outputName)) {
        CashFlowArraySP   thetaArray(new CashFlowArray(0));
        // create theta shift
        ThetaSP     thetaSens;
        if ( !useAssetFwds ) {
            thetaSens = ThetaSP(new Theta(offset, hols.getSP()));
        } else {
            thetaSens = ThetaSP(new ThetaFwdSpot(offset, hols.getSP()));
        }
        thetaSens->control   = control;
        thetaSens->algorithm = algorithm;
        /* create a copy of the instrument/model - better than trying
           to restore after each roll */
        TweakGroupSP  shifted(copy(tweakGroup));

        // create a fresh control for pricing
        DeltaSP              delta(new Delta(deltaShift));
        SensitivityArraySP   sens(
            new SensitivityArray(1,  SensitivitySP::
                                 dynamicCast((IObjectSP)delta)));
        OutputRequestArraySP req(new OutputRequestArray(0));
        CControlSP ctrl(new Control(sens, req, false, ""));

        OutputNameArrayConstSP deltaNames(delta->names(shifted.get()));
        SensitivityData deltaHash; // create one on stack

        int i;
        for (i=0 ; i<deltaNames->size() ; ++i ) {
            CashFlowArraySP   deltaArray(new CashFlowArray(0));
            deltaHash[(*deltaNames)[i]] = deltaArray;
        }

        for (i=0 ; i<thetaInterval->size() ; ++i ) {
            // roll forward to the start date
            while ( shifted->getInstrument()->getValueDate() < 
                    (*thetaInterval)[i].getStartDate() ) {
                thetaSens->applyScenario(shifted);
            }
            try {
                while (shifted->getInstrument()->getValueDate() <= 
                        (*thetaInterval)[i].getEndDate() ) {
                    Results       interimResult;
                    DateTime      thetaDate;
                    try {
                        // shift and calculate new price
                        thetaSens->applyScenario(shifted);
                        ctrl->calculate(shifted->getModel(), 
                                        shifted->getInstrument(),
                                        &interimResult);
                        // then store price against date
                        double shiftedPrice = interimResult.retrievePrice();
                        thetaDate = shifted->getInstrument()->getValueDate();
                        CashFlow dailyTheta(thetaDate, shiftedPrice);
                        thetaArray->push_back(dailyTheta);
                    } catch (exception& e) {
                        results->storeGreek(IObjectSP(new Untweakable(e)),
                                            NAME, 
                                            OutputNameSP(new OutputName("")));
                        return; //exit this routine (but not quite as an error)
                    }

                    // copy the interimResults to the main result
                    vector<pair<OutputNameConstSP, IObjectConstSP> > outputs;
                    outputs = interimResult.listPacketResults(Delta::NAME);

                    for (unsigned int j=0; j < outputs.size(); j++) {
                        double delta = interimResult.retrieveScalarGreek(
                            Delta::NAME,
                            outputs[j].first);

                        CashFlowArraySP   deltaArray = 
                            deltaHash[outputs[j].first];
                        CashFlow dailyDelta(thetaDate, delta);
                        deltaArray->push_back(dailyDelta);
                    }
                }
            } 
            catch (exception& e) {
                results->storeGreek(IObjectSP(new Untweakable(e)),
                                    Results::INSTRUMENT_PACKET, 
                                    outputName);
                untweakableDetected = true;
            }
        }

        if (!untweakableDetected)
        {
            results->storeGreek(thetaArray,
                                Results::INSTRUMENT_PACKET,
                                outputName);
        }
        
        for ( i=0 ; i<deltaNames->size() ; ++i ) {
            CashFlowArraySP   deltaArray = deltaHash[(*deltaNames)[i]];
            results->storeGreek(deltaArray,
                                RollingDeltaNAME,
                                (*deltaNames)[i]);
        }
    }
}

/** populate from market cache */
void RollingTheta::getMarket(const IModel* model, const MarketData* market) {
    hols.getData(model, market);
}

class RollingThetaHelper{
public:
    /** Invoked when SmallClass is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(RollingTheta, clazz);
        SUPERCLASS(Sensitivity);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultTheta);
        FIELD(thetaInterval,       "interval for rolling theta calculation");
        FIELD(deltaShift,   "shift size for delta shift");
        FIELD(offset,       "number of days to roll");
        FIELD(hols,         "holidays");
        FIELD(useAssetFwds, "use fwd prices when rolling to the "
        "next date?");
        FIELD_MAKE_OPTIONAL(useAssetFwds);
    }

    static IObject* defaultTheta(){
        return new RollingTheta();
    }
};

CClassConstSP const RollingTheta::TYPE = CClass::registerClassLoadMethod(
    "RollingTheta", typeid(RollingTheta), RollingThetaHelper::load);


DRLIB_END_NAMESPACE
