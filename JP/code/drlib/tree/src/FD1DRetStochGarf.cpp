//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FD1DRetStochGarf.cpp
//
//   Description : mixing of one factor trinomial tree for local vol process.
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/VolProcessedStochGarf.hpp"

#include "edginc/FD1DRetStochGarfSolver.hpp"

//#include "edginc/Format.hpp"

DRLIB_BEGIN_NAMESPACE

////////////////////////// FD1DRetStochGarf /////////////////////////
void FD1DRetStochGarf::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FD1DRetStochGarf, clazz);
    SUPERCLASS(FD1DRetLV);
    EMPTY_SHELL_METHOD(defaultFD1DRetStochGarf);

/*
    // transient for tweeking 
    FIELD(und,"");
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(und);
//    FIELD(vol, "");
//    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(vol);
    FIELD(discYC, "");
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(discYC);
    FIELD(FD1DRetLVs, "");
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(FD1DRetLVs);

    FIELD(volType, "Type of log-normal vol to use");
    FIELD_MAKE_OPTIONAL(volType);
*/

}

//----------------------------------------------------------------
// helpers
CClassConstSP const FD1DRetStochGarf::TYPE = CClass::registerClassLoadMethod(
    "FD1DRetStochGarf", typeid(FD1DRetStochGarf), load);

//----------------------------------------------------------------

FD1DRetStochGarf::FD1DRetStochGarf():FD1DRetLV(TYPE){
    numOfVols = 3;
    FD1DRetLVs.resize(3);
    for (int i=0;i<numOfVols;i++){
        FD1DRetLVs[i] = FD1DRetLVSP(new FD1DRetLV());
    }
};

//----------------------------------------------------------------

FD1DRetStochGarf::~FD1DRetStochGarf(){
}


//----------------------------------------------------------------

//return mdf->fetch(market, name, type); 


/** Override default createMDF in order to set the right MDF */

MarketDataFetcherSP FD1DRetStochGarf::createMDF() const {
    for (int i=0;i<numOfVols;i++)
        FD1DRetLVs[i]->createMDF();
    return MarketDataFetcherSP(new MarketDataFetcherLN(volType));
}

//----------------------------------------------------------------

/** get processed vol */
// Preparing for each LV models, here.
void FD1DRetStochGarf::initModel()
{
    static const string method("FD1DRetStochGarf::initModel");
    try {
        // first, init MV model by using FD1DRetLV....
        FD1DRetLV::initModel();

        int i;
        IIntoProduct* intoProd = dynamic_cast<IIntoProduct*>(myInst);
        // create the product
        for (i=0;i<numOfVols;i++) {
            modelsToUse[i]->prod = modelsToUse[i]->createProduct(IProdCreatorSP::attachToRef(intoProd));
        }

        VolRequestRaw request;
        VolProcessedStochGarfSP rawVol(dynamic_cast<VolProcessedStochGarf*>(underlying->getProcessedVol(&request)));
        // customize fd set up by product, rebuild flag allows same fd grid tweak
        for (i=0;i<numOfVols;i++){
            modelsToUse[i]->prod->init(myControl);
            modelsToUse[i]->FD1DRetLV::initModel();
            VolProcessedStochGarfSP vol(copy(rawVol.get()));
            vol->setScalar(-1.0+double(i));
            modelsToUse[i]->setVolLV(vol);
            
//                CVolProcessedDVFSP(rawVol.get()->VolStochGarfgetProcessedVol(-1.0+double(i))));

            // need to overwrite VolLV with different vol level.
//            CVolProcessedSP rawVol(und->getProcessedVol(&request));
//            VolStochGarf test(dynamic_cast<VolStochGarf*>(rawVol.get()), 0.0);
//            VolStochGarfSP myGarf( new VolStochGarf(dynamic_cast<VolStochGarf*>(rawVol.get()), 0.1));
            //VolStochGarfSP myGarf( new VolStochGarf(dynamic_cast< const VolStochGarf& >(*rawVol), 0.1));
//            modelsToUse[i]->FD1DRetLV::setVolLV(
//                    CVolProcessedDVFSP(new VolProcessedStochGarf(myGarf)));
        }
       
    } 
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// create FD1DRetLVs at here.
void FD1DRetStochGarf::finaliseModel(CControl*    control)
{
    static const string method("FD1DRetStochGarf::finaliseModel");
    try {
        int i;    
        // initialise product variables 
        for (i=0;i<numOfVols;i++)
            modelsToUse[i]->prod->initProd();

        // model param validation
        for (i=0;i<numOfVols;i++)
            modelsToUse[i]->FD1DRetLV::finaliseModel(control);

        // copy to member classes.
        for (i=0;i<numOfVols;i++)
            FD1DRetLVs[i] = modelsToUse[i];

    }catch (exception& e) {
        throw ModelException(e, method);
    }
}

//----------------------------------------------------------------

/*virtual*/ 
FDModel::IFDSolverSP FD1DRetStochGarf::createSolver()
{
    return IFDSolverSP(new FD1DRetStochGarfSolver(this)); 
}

/** main model entry point */
//----------------------------------------------------------------
void FD1DRetStochGarf::Price(CInstrument* instrument, 
                    CControl*    control, 
                    CResults*    results){

    static const string method = "FDModel::Price";

    if (!IIntoProduct::TYPE->isInstance(instrument)){
        throw ModelException("FDModel::Price", "Instrument of type "+
                             instrument->getClass()->getName() +
                             " does not support FDModel::IIntoProduct");
    }

    try {

        // keep control inside of this model.
        myControl = control;
        myInst = instrument;

        // need to prepare the LV models before call MV's InitProd
        // oterwise fails to copy InitData.
        int i;
        modelsToUse.resize(numOfVols);
        for (i=0;i<numOfVols;i++){
            modelsToUse[i] = FD1DRetLVSP(copy(this));
        }            


        FDModel::Price(instrument,control,results);

/*
        int i;
        FD1DRetLVSP modelToUse;
        vector<FD1DRetLVSP > modelsToUse;
        modelsToUse.resize(numOfVols);
        
        modelToUse = FD1DRetLVSP::attachToRef(this);

        for (i=0;i<numOfVols;i++){
            modelsToUse[i] = FD1DRetLVSP(copy(this));
            //modelsToUse[i] = FD1DRetLVSP::attachToRef(this);
        }            

        bool rebuild = true;
        // rebuild fd nodes when required
//        sameGridTweak = true;
//        bool rebuild = isRebuilt(control, sameGridTweak);
//        if(!control || control->isPricing() || !rebuild){
//            for (i=0;i<numOfVols;i++){
//                modelsToUse[i] = FDModelSP::attachToRef(this);
//            }            
//        } else {
//            for (i=0;i<numOfVols;i++){
//                modelsToUse[i] = FDModelSP(copy(this));
//            }            
//        }

        // prod creator
        IIntoProduct& intoProd = dynamic_cast<IIntoProduct&>(*instrument);
        // create the product
        modelToUse->prod = intoProd.createProduct(this);
        for (i=0;i<numOfVols;i++){
            modelsToUse[i]->prod = intoProd.createProduct(modelsToUse[i].get());
            //modelsToUse[i]->prod = intoProd.createProduct(this);
            //modelsToUse[i]->prod = modelToUse->prod;
        }

        // customize fd set up by product, rebuild flag allows same fd grid tweak
        if (rebuild) {
            for (i=0;i<numOfVols;i++){
                modelsToUse[i]->prod->init(control);
                modelsToUse[i]->initModel();    
            }
        }

        // initialise product variables 
        for (i=0;i<numOfVols;i++)
            modelsToUse[i]->prod->initProd();

        // model param validation
        for (i=0;i<numOfVols;i++)
            modelsToUse[i]->finaliseModel();

        // copy to member classes.
        for (i=0;i<numOfVols;i++)
            FD1DRetLVs[i] = modelsToUse[i];
        // create FD solver
        IFDSolverSP solver(modelToUse->createSolver());
        solver->roll();

        if (control && control->isPricing()){
            double deltaSize = control->getDeltaShiftSize();
            // default tweak size (0.5%) will give 5% of grid spacing
            //sameGridDeltaShiftSize = 10.0*deltaSize*(Stock[CurrIdx][1] - Stock[CurrIdx][-1])/Stock[CurrIdx][0];

            // store some extra results
            results->storeScalarGreek(modelsToUse[0]->timeLine->NumOfStep, Results::DEBUG_PACKET, 
                                      OutputNameSP(new OutputName("STEPS_USED")));
        }

        // record price and additional outputs 
        // need to review
        FD1DRetLVs[1]->prod->recordOutput(control, results);
*/
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

//----------------------------------------------------------------


bool FD1DRetStochGarfLoad(){
    return (FD1DRetStochGarf::TYPE && FD1DRetLVArray::TYPE);
}

DRLIB_END_NAMESPACE
