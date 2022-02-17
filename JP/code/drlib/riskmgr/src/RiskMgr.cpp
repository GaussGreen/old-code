
#include "edginc/config.hpp"
#include "edginc/RiskMgr.hpp"
#include "edginc/MultiRiskMgrInterface.hpp"
#include "edginc/RiskMgrInterface.hpp"
#include "edginc/ArrayInstrumentCollection.hpp"
#include "edginc/OutputFile.hpp"
#include "edginc/SpreadSheetMode.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/ScaleOutputs.hpp"
#include "edginc/DerivativeAsset.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/IntrinsicMTM.hpp"

DRLIB_BEGIN_NAMESPACE

/** Calculates price and greeks etc for possibly two ways of viewing the
    world. One is totally theoretical and the other is using mark to market
    for instruments which could be priced theortically but are also traded
    and therefore have a mark to market (mtm) value. */
CResultsArraySP RiskMgr::calculateTheoAndMtm(IModelSP                mdl,
                                             IInstrumentCollectionSP insts,
                                             CControlSP              control){
    static const string method("CRiskMgr::calculateTheoAndMtm");
    try{
        // deal with theo/mtm request
        bool wantTheoPrice = false;
        bool wantMTMPrice = false;
        const string& assetPriceSource = control->getAssetPriceSource();
        if (CString::equalsIgnoreCase(assetPriceSource, Control::ASSET_THEO)){
            wantTheoPrice = true;
        } else if (CString::equalsIgnoreCase(assetPriceSource,
                                             Control::ASSET_MTM)){
            wantMTMPrice = true;
        } else if (CString::equalsIgnoreCase(assetPriceSource,
                                             Control::ASSET_THEO_MTM)){
            wantTheoPrice = true;
            wantMTMPrice = true;
        } else {
            throw ModelException(method, "Unrecognised asset price source: "+
                                 assetPriceSource);
        }

        IInstrumentCollectionSP firstInstToUse(insts);
        // find out if we have any derivative assets
        bool derivAssetsExist = DerivativeAsset::derivativeAssetsExist(insts);
        if (derivAssetsExist){
            // swap the model over to make sure we always update the
            // control before each pricing
            mdl = IModelSP(DerivativeAsset::createModel(mdl, insts));
            // additionally clone the instrument in case roll to now is
            // different for theo and mtm
            if (wantTheoPrice && wantMTMPrice){
                firstInstToUse = IInstrumentCollectionSP(insts.clone());
            }
        } else {
            wantTheoPrice = true; // mtm irrelevant so better price something
        }
        // calculation call
        CResultsArraySP results;
        if (wantTheoPrice){
            // Dummy try/catch to try to avoid vc6.opt crashes
			try {
                results = CResultsArraySP(calculateTheoOrMtm(true,
                                                             derivAssetsExist,
                                                             mdl, firstInstToUse,
                                                             control));
            }
            catch (...) {
                throw;
            }
        } else {
            // mtm results set is always returned inside top level theo even
            // if it is empty
            results = insts->emptyResults();
        }
        firstInstToUse.reset(); // don't use again
        // price again for mtm if ASSET_THEO_MTM or first time if ASSET_MTM

        if (derivAssetsExist && wantMTMPrice){
            // price using original instrument
            CResultsArraySP mtmResults(calculateTheoOrMtm(false, true, mdl,
                                                          insts, control));

            {for (int i = 0; i < insts->size(); ++i) {
                // FIXME this is just horrible, is purely for an OCB hack
                // update with MTM if required
                IInstrument* inst = (*insts)[i].get();
                if (IIntrinsicMTM::TYPE->isInstance(inst)) {
                    IIntrinsicMTM* mtmImnt = dynamic_cast<IIntrinsicMTM*>(inst);
                    mtmImnt->calculateIntrinsic(control, (*mtmResults)[i]);

                    // perform scaling if needed
                    if (control->scaleOutputs() &&
                        IScaleOutputs::TYPE->isInstance(inst)) {
                        IObject* obj   = inst;
                        IScaleOutputs* scaleableInst = dynamic_cast<IScaleOutputs*>(obj);
                        scaleableInst->scaleOutputs(control, (*mtmResults)[i]);
                    }

                }
            }}

            // store sub results in the main object
            OutputNameConstSP mtmName(
                new OutputName(DerivativeAsset::ASSET_MTM_ID));
            for (int i = 0; i < results->size(); ++i) {
                (*results)[i]->storeGreek((*mtmResults)[i],
                                          Results::INSTRUMENT_PACKET, mtmName);
            }
        }

        // perform scaling if needed
        if (wantTheoPrice && control->scaleOutputs()) {
            insts->scaleOutputs(control, results);
        }

        return results;
    } catch(exception& e){
        throw ModelException(e, method);
    }
}
// See below --- this is a workaround for obscure VC6 crashes
//@{

static void maybeWriteToFile(CControl* control, ModelException& e) {
    if (control->getWriteToFile()) {
        string fn = OutputFile::createOutputFileName(
            control->getFileName(false));
        OutputFile n(fn);
        n.write(e);
    }
}

static string theMessage;

void setTheMessage(IInstrumentCollection* instruments, IModel* model) {
    CClassConstSP clazz =
        dynamic_cast<ArrayInstrumentCollection *>(instruments) &&
        instruments->size() > 0 ?
            (*instruments)[0]->getClass() : instruments->getClass();

    theMessage =
        "CRiskMgr::calculateTheoOrMtm: Failed for " + clazz->getName() +
        " using " + model->getClass()->getName() + " model";
}

//@}

/** similar to above, but requires an instrument who has
    already called GetMarket - also writes outputs to file. Note that the
    supplied instrument is modified via a "roll to now" */
CResultsArraySP CRiskMgr::calculateTheoOrMtm(
    bool                    useTheoAssetPrice,
    bool                    derivAssetsExist,
    IModelSP                model,
    IInstrumentCollectionSP instruments,
    CControlSP              control)
{
    static const string method("CRiskMgr::calculateTheoOrMtm");
    CResultsArraySP results;

    CControl* controlp = control.get();

    // Ridiculuous, but putting this message construction logic inside the
    // "catch" caused no end of obscure crashes in VC6.opt.  The crashes
    // are in IInstrumentCollectionSP::~IInstrumentCollectionSP, but the
    // key thing that seems to help is making theMessage a static
    // variable (I kid you not).

    setTheMessage(instruments.get(), model.get());

    try {
        // Inform derivative assets what's expected of them. Need to do
        // this before roll to now in case asset spots/fwds needed there
        control->setUseTheoAssetPrice(useTheoAssetPrice);
        if (derivAssetsExist){
            DerivativeAsset::setUseTheoAssetPrice(instruments,
                                                  useTheoAssetPrice);
        }

        // create TweakGroup (holds what we want to tweak)
        MultiTweakGroupSP tweakGroup(new MultiTweakGroup(instruments, model));
        /* roll value date forward by 0 days - this is to populate
           any samples which should be set now, but are not
           populated yet, for example when running overnight grids
           for instruments which have a SOD sample */
        Theta   thetaShift(0, HolidaySP(Holiday::noHolidays()));
        thetaShift.applyScenario(tweakGroup);

        // and then price
        // Dummy try/catch to try to avoid vc6.opt crashes
		try {
            results = model->RunMulti(instruments, control.get());
		}
        catch (...) {
            throw;
        }
    }
    catch (ModelException& e) {
        // See above --- putting logic directly in here causes obscure VC60.opt
        // crashes
        e.addMsg(theMessage);
        maybeWriteToFile(controlp, e);
        throw e;
    }
    return results;
}

void CRiskMgr::validate(const string& method,
                        IModelSP model,
                        IInstrumentCollectionSP instruments,
                        CControlSP control,
                        MarketDataSP market) {
    if (!model){
        throw ModelException(method, "NULL model");
    }
    if (!instruments){
        throw ModelException(method, "NULL instruments");
    }
    if (!control){
        throw ModelException(method, "NULL control");
    }
    if (!market){
        throw ModelException(method, "NULL market");
    }
}

/** Constructs a MultiRiskMgrInterface object and writes it to file. Such an
    object can be used directly with the regression tester */
void CRiskMgr::writeInputs(IModelSP                model,
                           IInstrumentCollectionSP instruments,
                           CControlSP              control,
                           MarketDataSP            market)
{
    MultiRiskMgrInterface rmi(model, instruments, control, market);
    XMLWriter xml(control->getFileName(true));
    rmi.write("RISK-MGR", &xml);
}

/** Constructs a RiskMgrInterface object and writes it to file. Such an
    object can be used directly with the regression tester */
void CRiskMgr::writeInputs(IModelSP                model,
                           CInstrumentSP           instrument,
                           CControlSP              control,
                           MarketDataSP            market)
{
    RiskMgrInterface rmi(model, instrument, control, market);
    XMLWriter xml(control->getFileName(true));
    rmi.write("RISK-MGR", &xml);
}

string CRiskMgr::getOutputFileName(Control* control){
    return OutputFile::createOutputFileName(control->getFileName(false));
}

/** Writes an output file containing supplied object. Name is based on fileName
    in control */
void CRiskMgr::writeOutputs(Control*      control,
                            IObjectSP     output) {
    OutputFile fileOut(getOutputFileName(control));
    fileOut.write(output.get());
}
DRLIB_END_NAMESPACE
