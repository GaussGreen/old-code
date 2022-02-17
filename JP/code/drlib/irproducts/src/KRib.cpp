//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KRib.cpp
//
//   Description : Rib component
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/KRib.hpp"
#include "edginc/KKnockOut.hpp"
#include "edginc/RatesUtils.hpp"

DRLIB_BEGIN_NAMESPACE

class KRibTree : public FDProduct {
public:
    /************************ variables ************************/
    const KRib*     inst;    
   
    FDProductSP     obs;     /* observation underlying is a composite KKnockOut */
    FDProductSP     cpn;     /* coupon composed of weighted underlyings  */

    int             nbObs;   /* number of observation dates per accrual period */

    /************************ methods ************************/

    KRibTree(const KRib* inst, FDModel* model);

    virtual DateTime getStartDate(void) const {return model->getValueDate();}
    virtual string getOutputName(void) const { return inst->outputName; }
    virtual void update(int& step, UpdateType type);
    virtual const TreeSlice & getValue(int step, DateTime eventDate) const { return *cashFlow; }
    virtual void recordOutput(Control* ctrl, YieldCurveConstSP, Results* results) 
    {
        /* It does not make sense to record the rib component value at t=0
        ecordSliceToOutputName(ctrl, results, model, 
            false, inst->outputName, "", getValue(0));
        */
    }

    /** initialisation, called ONCE only for each new model instance */
    virtual void init(Control*) const;

    /** initialising and setting product variables */
    virtual void initProd();

protected:
    TreeSliceSP cashFlow;
    TreeSliceSP ribCpn;
    TreeSliceSP ribIndex;
};

/******************************** KRibTree ********************************/

KRibTree::KRibTree(const KRib* inst, FDModel* model) :
    FDProduct(model), inst(inst)
{
    obs = model->createProduct(inst->obsUnd);
    obs->addModelResetDates(inst->obsDates);

    cpn = model->createProduct(inst->cpnUnd); 
    const DateTimeArray& vDates = (inst->resetInArrears ? inst->accEndDates : inst->accStartDates);
    cpn->addModelResetDates(vDates);
}

void KRibTree::init(Control*) const
{
    /* add observation, coupon reset, accrual period start dates to the list of critical dates */
    if (inst->obsEffDates.size() != inst->obsDates.size())
        throw ModelException("KRibTree::init", 
                             "Rates tree model requires RIB observation effective/spot dates as "
                             "well as actual observation dates");

    model->addCritDates(inst->obsEffDates);
    model->addCritDates(inst->accStartDates);
    model->addCritDates(inst->accEndDates);
}

/** initialising and setting product variables */
void KRibTree::initProd(void)
{
    try {
        nbObs = 0;

        // CREATE SLICES: 
        // !!! review ccy names!!!
        string curve = inst->discountYieldCurveName();

        cashFlow = model->createSlice(curve);
        cashFlow->name = inst->outputName + "_cashFlow";

        ribCpn = model->createSlice(curve); 
        ribCpn->name = inst->outputName + "_ribCpn";

        ribIndex = model->createSlice(curve);
        ribIndex->name = inst->outputName + "_ribIndex";

        //RIB index (different in case of reset in-advance or in-arrear)
    }
    catch (exception& e) {
        throw ModelException(e, "KRibTree::initProd()");
    }
}

/** update products slice at each step */
void KRibTree::update(int& step, UpdateType type)
{
    try 
    {    
        DateTime stepDate = model->getDate(step);
        int accrualPos;

        /* at the observation date the observation profile is created, and RIB coupon is updated */

        if (RatesUtils::happensNow(inst->obsEffDates, stepDate))
        {  
            const TreeSlice& obsVal = obs->getValue(step, stepDate);
            *ribCpn += *ribIndex * obsVal;
            //ribCpn += ribIndex * obs->getValue(step);
            nbObs++;
            startDEV(ribCpn);
        }

        /* at the accrual start: 
        if reset-in-advance -- multiply [ W_true*N_true + W_false*N_false ] (accumulated in ribCpn), 
        and weighted sum of coupon underlyings: w_0*und_0 + w_1*und_1 + ... (accumulated in cpnTmp)
        if reset-in-arrears -- ribCpn already contains all contributions */

        if (RatesUtils::happensNow(inst->accStartDates, stepDate, &accrualPos))
        {
            if ( inst->resetInArrears == false)
            {                        
                // The spread should be added here, not in KFloatLeg's update().
                // This means that during the knock-out period, nothing is accrued.
                const TreeSlice& cpnVal = cpn->getValue(step, stepDate);
                double spread = 0.;
                if (inst->spreads.get())
                    spread = (*inst->spreads)[accrualPos];
                *ribCpn *= ( cpnVal + spread );
                //ribCpn *= ( cpn->getValue(step) + inst->spreads[event.accStart] );           
                //startDEV(ribCpn);
            }

            /* divide the coupon by the number of observation dates in the current accrual period */ 
            *ribCpn *= 1./nbObs;
            // The parent of the KRib, which is a KFloatLeg, always collects ribCpn on the accrual start date.
            // We save ribCpn value in cashFlow 
            *cashFlow = *ribCpn;
            // startDEV(ribCpn);
        }

        /* at coupon reset date:
        if reset-in-advance -- set ribIndex to 1, ribCpn to 0;
        if reset-in-arrears -- set ribIndex to the weighted sum of the coupon underlyings w_0*und_0 + w_1*und_1 + ...
        evaluated at the current date, and set ribCpn to 0 */

        if (RatesUtils::happensNow(inst->accEndDates, stepDate, &accrualPos))
        {
            if (inst->resetInArrears == true)
            {                           
                // The spread should be added here, not in KFloatLeg's update().
                // This means that during the knock-out period, nothing is accrued.
                const TreeSlice &cpnVal = cpn->getValue(step, stepDate);
                double spread = 0.;
                if (inst->spreads.get())
                    spread = (*inst->spreads)[accrualPos];

                *ribIndex = cpnVal + spread;
                //ribIndex = cpn->getValue(step) + inst->spreads[event.accEnd];                    
            }
            else 
                *ribIndex = 1.;

            startDEV(ribIndex);
            *ribCpn = 0.;
            nbObs = 0;
            stopDEV(ribCpn);
        }
    } 
    catch (exception& e) 
    {
        throw ModelException(e, "KRibTree::update()");
    }
}

/********************************** KRib **********************************/

/** implement KComponent virtual function */
void KRib::setup(const IModel* model, const MarketData* market) {
    try {
        KComponent::setup(model, market);
        if (obsUnd.get()) {
            // ??? for now, this must be a KKnockOut object
            KKnockOut* ko = dynamic_cast<KKnockOut*>(obsUnd.get());
            if (!ko)
                throw ModelException("obsUnd field must be of type KKnockOut component");
    	    obsUnd->addResetDates(obsDates);
            obsUnd->setup(model, market);

        }
        if (cpnUnd.get()) {
            const DateTimeArray& vDates = (resetInArrears ? accEndDates : accStartDates);
            cpnUnd->addResetDates(vDates);
            cpnUnd->setup(model, market);
        }

        if (obsDates.size() < 1) 
            throw ModelException("At least one  observation date is required");
        if (accStartDates.size() < 1) 
            throw ModelException("At least one  accrual start date is required");
        
        if (accStartDates.size() != accEndDates.size())
            throw ModelException("Number of accrual start dates " + Format::toString(accStartDates.size()) +
                                    " must equal number of accrual end dates" +
                                    Format::toString(accEndDates.size()));
        if (spreads.get() && accStartDates.size() != spreads->size())
        throw ModelException("Number of spread dates " + Format::toString(spreads->size()) +
                                " must equal number of accrual start/end dates " + 
                                Format::toString(accStartDates.size()));

        // if obsEffDates supplied, must be same size as obsDates
        if (obsEffDates.size()) {
            int i;
            if (obsEffDates.size() != obsDates.size())
                throw ModelException("If supplied, number of observation effective dates " +
                                     Format::toString(obsEffDates.size()) + 
                                     " must equal number of observation dates " +
                                     Format::toString(obsDates.size()));

            for (i = 0; i < obsDates.size(); i++) {
                if (obsDates[i] > obsEffDates[i])
                    throw ModelException("obsDates[" + Format::toString(i) + "] " + 
                                         obsDates[i].toString() + " must be <= obsEffDate " +
                                         obsEffDates[i].toString());
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, "KRib::setup, "+outputName);
    }

}


/** implement FDModel product interface */
FDProductSP KRib::createProduct(FDModel * model) const {
    if (!setupCalled) throw ModelException(getClass()->getName()+" \""+outputName+"\" not setup");
    return FDProductSP(new KRibTree(this, model));
}

double KRib::getValue(DateTime date, CashflowInfo &cfi) const
{
    try {
        int accrualPos;
        CashflowInfo cfiLoc;
        if (!RatesUtils::happensNow(accStartDates, date, &accrualPos)) {
            throw ModelException("Can be called only on accrual start dates");
        }
        DateTime accStart = accStartDates[accrualPos];
        DateTime accEnd = accEndDates[accrualPos];
        int nbObs = 0;
        double coupon = 0.;
        for (int i=0; i<obsDates.size(); ++i) {
            if (accStart <= obsDates[i] && obsDates[i] < accEnd) {
                coupon += 1.;
                ++nbObs;
            }
        }
        coupon /= nbObs;
        if (cpnUnd.get()) {
            coupon *= cpnUnd->getValue(resetInArrears ? accEnd : accStart, cfiLoc);
        }
        cfiLoc.push("nbObs", (double)nbObs);
        cfi.merge(cfiLoc, outputName);
        return coupon;
    }
    catch (exception& e) {
        throw ModelException(e, "KRib::getPastValue");
    }
}

void KRib::load(CClassSP& clazz)
{
    clazz->setPublic();
    clazz->setDescription("Rib component");
    REGISTER(KRib, clazz)
    SUPERCLASS(KComponent);
    EMPTY_SHELL_METHOD(KRib::defaultConstructor);
    IMPLEMENTS(FDModel::IIntoProduct);
    //IMPLEMENTS(LastSensDate);
    FIELD(obsUnd, "RIB events observation composite index");
    FIELD_MAKE_OPTIONAL(obsUnd);
    FIELD(cpnUnd, "payment coupon rate");
    FIELD_MAKE_OPTIONAL(cpnUnd);
    FIELD(obsDates, "RIB actual observation date");
    FIELD(obsEffDates, "ie. spot date for observation");
    FIELD_MAKE_OPTIONAL(obsEffDates);
    FIELD(accEndDates, "");
    FIELD_MAKE_OPTIONAL(accEndDates);
    FIELD(accStartDates, "");
    FIELD_MAKE_OPTIONAL(accStartDates);
    FIELD(resetInArrears, "");
    FIELD_MAKE_OPTIONAL(resetInArrears);
    FIELD(spreads,"Spreads on rib payment index");
    FIELD_MAKE_OPTIONAL(spreads);
    Addin::registerConstructor(Addin::UTILITIES, KRib::TYPE);
}

CClassConstSP const KRib::TYPE = CClass::registerClassLoadMethod(
    "KRib", typeid(KRib), KRib::load);

/**********************************/
bool KRibLoad()
{
    return KRib::TYPE !=0;
}

DRLIB_END_NAMESPACE
