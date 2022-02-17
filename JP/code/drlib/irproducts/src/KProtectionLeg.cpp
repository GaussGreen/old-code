//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KProtectionLeg.cpp
//
//   Description : Protection
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/ICDSParSpreads.hpp"
#include "edginc/KProtectionLeg.hpp"
#include "edginc/CreditTree.hpp"

DRLIB_BEGIN_NAMESPACE

/////////// KProtectionLegTree class /////////////
class KProtectionLegTree : public FDProduct {
public:

    /******************** variables ********************/
    KProtectionLegConstSP inst;
    ProtectionLegSP mainProtection;  // for protection from start to end
    ProtectionLegSP frontProtection; // if unconditional settlement
    CreditTree*      crTree;

    /******************** methods ********************/
    // constructor
    KProtectionLegTree(const KProtectionLegConstSP &inst, FDModel* model);

    /** initialisation, called ONCE only after InitState() for each new model instance */
    virtual void init(Control*) const;

    /** initialising and setting product variables */
    // this is called per pricing call before each pricing
    virtual void initProd(void);

    /** isInitValue == true, payoff at T for backward or value at t=0 for fwd induction
        isInitValue == false, payoff boundary condition, for KO, early exercise etc. */
    virtual void update(int & step, FDProduct::UpdateType);
    virtual const TreeSlice & getValue(int step, DateTime eventDate) const;

    virtual DateTime getStartDate(void) const { return model->getDate(0);}

    virtual string getOutputName() const;
    void printInfo(ostream& outputStream) const;

    virtual void recordOutput(Control* ctrl, YieldCurveConstSP, Results* results) {
        string ccy;
        if (inst->discount.get()) ccy = inst->discount->getCcy();
        recordSliceToOutputName(ctrl, results, model, 
            inst->isTopComponent(), inst->outputName, ccy, getValue(0, model->getDate(0)));
        inst->recordExtraOutput(ctrl, results);
    }
private:
    TreeSliceSP mainSlice, getValueSlice;
    // to help debugging and know which product instance this is, take copy of the instrument's
    // outputName field as the debuggers are unable to get from instrumentSP in the watch windows
    string instOutputName;  
};

/******************************* KProtectionLegTree ******************************/

KProtectionLegTree::KProtectionLegTree(const KProtectionLegConstSP &inst, FDModel* model) :
    FDProduct(model), inst(inst), 
        mainProtection(new ProtectionLeg(inst->discountYieldCurveName(),
            inst->creditCurveName(),inst->protectionStartDate, inst->protectionEndDate,
            false /*pay defaults as they occur*/,
            inst->kLo,inst->kHi)), frontProtection(0) {

    const static char* method = "KProtectionLegTree::KProtectionLegTree";
    crTree = dynamic_cast<CreditTree*>(model);
    if (!crTree) {
        throw ModelException(method,
        "Only CreditTree models are supported!");
    }

    const string& discYCName = inst->discountYieldCurveName();
    if (discYCName.empty())
        throw ModelException(method, 
                             "Discount YieldCurve name of component must be set.");
    const string& credCrvName = inst->creditCurveName();
    if (credCrvName.empty()) {
        throw ModelException(method, 
                             "Credit Curve name of component must be set.");
    }


    /* Add protection until start date payable at start date, if unconditional */
    if (inst->settlementIsUnconditional) {
        frontProtection = ProtectionLegSP(new ProtectionLeg(inst->discountYieldCurveName(),
            inst->creditCurveName(),min(model->getValueDate(),inst->protectionStartDate), inst->protectionStartDate,
            true /*Pay front protection at start date*/,
            inst->kLo,inst->kHi));
    }
}

void KProtectionLegTree::init(Control*) const{
    /* ADD CRITICAL DATES: PROTECTION START AND END DATES */
    DateTimeArray cdt;
    cdt.push_back(inst->protectionStartDate);
    cdt.push_back(inst->protectionEndDate);
    model->addCritDates(cdt);
}

void KProtectionLegTree::initProd(void){

    string discCurve = inst->discount->getName();

    mainSlice = model->createSlice(discCurve, discCurve);
    mainSlice->name = inst->outputName + "_mainSlice";
    *mainSlice = 0.;

    getValueSlice = model->createSlice("", discCurve);
    getValueSlice->name = inst->outputName + "_getValueSlice";

    instOutputName = inst->outputName;  // to help debugging, store local copy of instrument name
}

/* Extract value of instrument at given time step */
const TreeSlice & KProtectionLegTree::getValue(int step, DateTime eventDate) const {
    DateTime stepDate = model->getDate(step);

    if (eventDate != stepDate) {
        throw ModelException(__FUNCTION__, "Cannot be valued at a date != currentDate");
    }

    *getValueSlice = *mainSlice;
    if (stepDate>inst->protectionStartDate && stepDate<=inst->protectionEndDate) {
        // during main protection period: value needs to be added
        // since not in mainSlice yet
        *getValueSlice += inst->notional * crTree->getProtectionLegValue(*mainProtection, step);;
    }
    if (inst->settlementIsUnconditional && stepDate<inst->protectionStartDate) {
        // during front protection period: value needs to be added
        // since not in mainSlice yet
        *getValueSlice += inst->notional * crTree->getProtectionLegValue(*frontProtection, step);
    }
    return *getValueSlice; 
}


void KProtectionLegTree::update(int & step, FDProduct::UpdateType update) {
    try {
        DateTime stepDate = model->getDate(step);

        if (stepDate==inst->protectionStartDate) {
            // add main protection value when known at protection start date
            *mainSlice += crTree->getProtectionLegValue(*mainProtection,step);
            startDEV(mainSlice);
        }
        
        // NB never add front protection in update - just add in getValue
    }
    catch (exception& e) {
        throw ModelException(e, "KFloatLegTree::update");
    }
}

string KProtectionLegTree::getOutputName() const{
    return inst->outputName;
}

void KProtectionLegTree::printInfo(ostream& outputStream) const {
    string prodType = "KProtectionLegTree";
    string name = getOutputName();

    outputStream << "Product type: " << prodType << " name: " << name << endl << endl;

    /* WRITE SECTION PRINTING OUT PRODUCT DETAILS */

    // internal slices and slice dimensions info - change slice type, though!
    outputStream << endl;
    outputStream << "Internal product slice details:" << endl;
    outputStream << "instanceName   dimension   devCurveName   descriptiveName" << endl;
    outputStream << "mainSlice   " << DYNAMIC_POINTER_CAST<TreeSliceRates>(mainSlice)->getDim() << "   " <<
                     mainSlice->getCurveToDEV() << endl; // "   " << mainSlice->getSliceName() << endl;
    outputStream << "getValueSlice   " << DYNAMIC_POINTER_CAST<TreeSliceRates>(getValueSlice)->getDim() << "   " <<
                     getValueSlice->getCurveToDEV() << endl; // "   " << getValueSlice->getSliceName() << endl;

}

/*=================================================================================================
 * KProtectionLeg methods
 *===============================================================================================*/
void KProtectionLeg::validatePop2Object(void) {
    const static char* method = "KProtectionLeg::validatePop2Object";
    if (protectionStartDate>=protectionEndDate) {
        throw ModelException(method, " Protection start date must be strictly before protection end date.");
    } else if (kLo<0) {
        throw ModelException(method, " Lower attachment point must be >=0.");
    } else if (kHi>1) {
        throw ModelException(method, " Upper attachment point must be <=1.");
    } else if (kLo>=kHi) {
        throw ModelException(method, " Lower attachment point must be lower than upper attachment point.");
    }
}

/* NOT SURE WHAT TO DO WITH THIS YET - THERE USUALLY AREN'T ANY! */
void KProtectionLeg::reportEvents(const KnownCashflows*, IModel* model,
    const DateTime& eDate, EventResults* events) const {
    
}

void KProtectionLeg::setup(const IModel* model, const MarketData* market) {
    // TODO: should maybe repeat checking - otherwise no more to do.
}

DateTime KProtectionLeg::getLastDate() const{
    return protectionEndDate;
}

FDProductSP KProtectionLeg::createProduct( FDModel * model ) const{
    return FDProductSP(new KProtectionLegTree(KProtectionLegConstSP(this), model));
}

void KProtectionLeg::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDescription("Protection leg component");
    REGISTER(KProtectionLeg, clazz);
    SUPERCLASS(KRiskyComponent);
    IMPLEMENTS(FDModel::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(notional, "Notional amount: positive is long protection.");
    FIELD(protectionStartDate, "Date at which protection starts");
    FIELD(protectionEndDate, "Date at which protection ends");
    FIELD(isMarketRecovery, "True (default) if losses are determined by the curve's market recovery rate; if false, digital.");
    FIELD_MAKE_OPTIONAL(isMarketRecovery);
    FIELD(digitalRecoveryRate, "Recovery rate to use is marketRecovery is false.");
    FIELD_MAKE_OPTIONAL(digitalRecoveryRate);
    FIELD(kLo, "Losses are recovered only after they have passed this level (default = 0)");
    FIELD_MAKE_OPTIONAL(kLo);
    FIELD(kHi, "Losses are recovered only until they reach this level (default = 1)");
    FIELD_MAKE_OPTIONAL(kHi);
    FIELD(settlementIsUnconditional, "If true (default is false), losses before protectionStartDate are paid at protectionStartDate; otherwise, no payments losses before then.");
    FIELD_MAKE_OPTIONAL(settlementIsUnconditional);

    Addin::registerConstructor(Addin::UTILITIES, KProtectionLeg::TYPE);
}

CClassConstSP const KProtectionLeg::TYPE = CClass::registerClassLoadMethod(
    "KProtectionLeg", typeid(KProtectionLeg), KProtectionLeg::load);

bool KProtectionLegLoad(){ return (KProtectionLeg::TYPE != 0); }

DRLIB_END_NAMESPACE
