//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : QPVanilla.cpp
//
//   Description : Simple vanilla pricer suitable for QuickPricer spreadsheet
//
//   Author      : Andrew J Swain
//
//   Date        : 28 March 2002
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/Vanilla.hpp"
#include "edginc/SimpleEquity.hpp"
#include "edginc/Tree1fLN.hpp"
#include "edginc/CashSettlePeriod.hpp"
#include "edginc/ImpliedScalarShift.hpp"
#include "edginc/ATMVolRequest.hpp"

DRLIB_BEGIN_NAMESPACE

class QPVanilla: public CObject {
public:
    static CClassConstSP const TYPE;

private:
    string                 callPut;
    string                 amerEuro;
    double                 strike;
    DateTime               maturity;
    EquitySP               equity;
    CVolBaseSP             vol;
    InstrumentSettlementSP settle;
    double                 userSpot;
    double                 userVol;
    double                 target;
    string                 control;

    static void quickCredit(
        const CAsset*   asset,
        const DateTime& today,
        const DateTime& maturity,
        bool            isCall,
        double          strike,
        Results*        results) {
        static const string routine = "QPVanilla::quickCredit";
        try {
            ATMVolRequestSP atm(new ATMVolRequest());

            CVolProcessedBSSP volBS(asset->getProcessedVol(atm.get()));
           
            double variance = volBS->CalcVar(today, maturity);
            double fwd = asset->fwdValue(maturity);
            double sign = isCall ? 1.0 : -1.0;

            double peak = Maths::max(sign*(fwd*exp(-variance/2.0 + 
                                                   sign*1.96*sqrt(variance))
                                           -strike), 
                                     0.0);

            OutputNameSP name(new OutputName("PEAK"));
            results->storeScalarGreek(peak, "CREDIT", name);
        } 
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }  

    /** the 'addin function' */
    static IObjectSP run(QPVanilla* params){
        static const string routine = "QPVanilla::run";
        try {
            bool isCall = CString::equalsIgnoreCase(params->callPut,"C");
            bool american = CString::equalsIgnoreCase(params->amerEuro,"A",1);

            DateTimeArray dates;
            DoubleArray   strikes;

            dates.push_back(params->maturity);
            strikes.push_back(params->strike);

            ScheduleSP exercise(new Schedule(dates, 
                                             strikes, 
                                             american ? Schedule::INTERP_LINEAR : 
                                             Schedule::INTERP_NONE));
            
            // if american, see if exercise window is passed in (e.g. A2 means
            // american with 2 day window
            int exerWindow = 0;

            if (american && params->amerEuro.length() > 1) {
                sscanf(params->amerEuro.c_str(), "%*c%d", &exerWindow);
            }

            // allow user to skip settlement - get T+0 then
            InstrumentSettlementSP settle;
            if (!params->settle) {
                HolidaySP hols(Holiday::noHolidays());

                settle=InstrumentSettlementSP(new CashSettlePeriod(0, 
                                                                   hols.get()));
            }
            else {
                settle = params->settle;
            }

            // spot/vol override
            CVolBaseSP vol;
            EquitySP   equity;

            if (params->userSpot > 0.0) {
                equity = EquitySP(copy(params->equity.get()));
                SpotLevel shift(params->userSpot);
                equity->sensShift(&shift);
            }
            else {
                equity = params->equity;
            }

            // bit tough on trading time
            if (params->userVol > 0.0) {
                HolidaySP  hols(Holiday::noHolidays());
                TimeMetric metric(1.0, hols.get());

                vol = CVolBaseSP(new FlatVol(params->vol->getName(),
                                             params->equity->getValueDate(),
                                             &metric,
                                             params->userVol));
            }
            else {
                vol = params->vol;
            }

            SimpleEquitySP asset(new SimpleEquity(equity.get(), vol.get())); 

            CVanillaSP vanilla(CVanilla::make(params->equity->getValueDate(),
                                              isCall,
                                              american,
                                              exercise.get(),
                                              asset.get(),
                                              params->equity->getYC().get(),
                                              settle.get(),
                                              exerWindow));

            IModelSP model;
            string volType = params->vol->getClass()->getName();

            if (american) {
                model = IModelSP(Tree1fLN::make(volType));
            }
            else {
                model = IModelSP(new CClosedFormLN(volType));
            }

            MarketDataConstSP market( new MarketData( params->equity->getValueDate() ) );
            model->getInstrumentAndModelMarket( market.get(), vanilla.get() );

            CControlSP ctrl(Control::makeFromFlags(params->control, 
                                                   params->target));

            // can't realistically add implied vol shift on control construction
            // unless we want a model in the interface
            bool credit = false;
            for (int i = 0; i < (int)params->control.length(); i++) {
                char iflag = toupper(params->control[i]);
                if (iflag == 'J') {
                    smartPtr<ScalarPerturbation> volScenario(
                        new VolParallelShift(0.16));
                    OutputNameSP  impliedName(new OutputName("VALUE"));
                    RootFinder1D::TwoInitValNoDerivSP rootFinder(
                        new ImpliedScalarShift::DefaultRootFinder(0.0001));

                    SensitivitySP impliedShift(
                        new ImpliedScalarShift(volScenario,
                                               params->target,
                                               rootFinder,
                                               impliedName));
                    ctrl->addSensitivity(impliedShift);
                }
                else if (iflag == 'E') {
                    credit = true;
                }
            }

            CResultsSP results(model->Run(vanilla.get(), ctrl.get()));

            if (credit) {
                quickCredit(asset.get(),
                            params->equity->getValueDate(),
                            params->maturity,
                            isCall,
                            params->strike,
                            results.get());
            }

            return results;
        } 
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    QPVanilla(): CObject(TYPE), userSpot(0.0), userVol(0.0), target(0.0){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(QPVanilla, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultQPVanilla);
        FIELD(callPut, "callPut");
        FIELD(amerEuro, "amerEuro");
        FIELD(strike, "strike");
        FIELD(maturity, "maturity");
        FIELD(equity, "equity");
        FIELD(vol, "vol");
        FIELD(settle, "settle");
        FIELD_MAKE_OPTIONAL(settle);
        FIELD(userSpot, "userSpot");
        FIELD_MAKE_OPTIONAL(userSpot);
        FIELD(userVol, "userVol");
        FIELD_MAKE_OPTIONAL(userVol);
        FIELD(target, "target");
        FIELD_MAKE_OPTIONAL(target);
        FIELD(control, "control");
        FIELD_MAKE_OPTIONAL(control);

        Addin::registerClassObjectMethod("QP_VANILLA",
                                         Addin::RISK,
                                         "Simple vanilla pricer",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)run);
    }

    static IObject* defaultQPVanilla(){
        return new QPVanilla();
    }
    
};

CClassConstSP const QPVanilla::TYPE = CClass::registerClassLoadMethod(
    "QPVanilla", typeid(QPVanilla), load);

/* for class loading */
bool QPVanillaLoad() {
    return (QPVanilla::TYPE != 0);
}

DRLIB_END_NAMESPACE
