//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KTarn.cpp
//
//   Description : Tarn component
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/KTarn.hpp"
#include "edginc/RatesUtils.hpp"
#include "edginc/TreeSliceOperExt.hpp"
#include "edginc/RateTree.hpp"

DRLIB_BEGIN_NAMESPACE

class KTarnTree : public FDProduct,
                  public TreeSliceLayer::StateSupport
{
    FDModel * model;
public:
    KTarnTree(KTarnConstSP inst, FDModel* model);

    /* TreeSliceLayer::StateSupport */
    virtual bool gridChanges(int step) const { return gridChanged; }
    virtual void setGrid(int step, bool init = false);
    virtual void transition(
        const vector< const TreeSlice * > & srcSlices,
        const vector< double > & gridIn,
        vector< double > & gridOut ) const;

    /* FDProduct */
    virtual DateTime getStartDate(void) const { return model->getDate(0);}
    virtual string getOutputName(void) const { return inst->outputName; }
    virtual void init(Control*) const;
    virtual void initProd(void);
    virtual void update(int& step, UpdateType type);
    virtual const TreeSlice & getValue(int step, DateTime eventDate) const { 
        return *value = *cplxValue - *fundValue; 
    }
    virtual void recordOutput(Control* ctrl, YieldCurveConstSP, Results* results);
    

protected:
    bool gridChanged;
    KTarnConstSP inst; // this overrides base
    FDProductSP complexLegProd;
    FDProductSP fundingLegProd;

    TreeSliceSP cashFlow;
    TreeSliceSP cplxValue;
    TreeSliceSP fundValue;
    TreeSliceSP value;
    TreeSliceArray cplxLegCF; // to keep the complex leg cashFlows until we use them
    TreeSliceArray fundLegCF; // to keep the complex leg cashFlows until we use them
    DateTimeArray cplxAddDates; // contains min of cf dates for cplx and fund coupons.
    DateTimeArray fundAddDates; // contains min of cf dates for cplx and fund coupons.
    DateTimeArray cplxCfDates;
    DateTimeArray fundCfDates;
};

/****************************** KTarnTree ********************************/

KTarnTree::KTarnTree(KTarnConstSP inst, FDModel* model) :
    FDProduct(model),
    TreeSliceLayer::StateSupport( "TarnAccuCF", inst->pastRealizedAmount, INTERP_LINEAR ),
    model(model),
    gridChanged(false),
    inst(inst)
{
    try {
        complexLegProd = model->createProduct(inst->complexLeg);
        fundingLegProd = model->createProduct(inst->fundingLeg);

        DYNAMIC_POINTER_CAST<IGetCfDate>(complexLegProd)->getCfDate(cplxCfDates);
        DYNAMIC_POINTER_CAST<IGetCfDate>(fundingLegProd)->getCfDate(fundCfDates);
        cplxLegCF.resize(cplxCfDates.size());
        fundLegCF.resize(fundCfDates.size());
        cplxAddDates = cplxCfDates;
        fundAddDates = fundCfDates;
        
        int fundCpn, cplxCpn;
        for (cplxCpn=0; cplxCpn<cplxCfDates.size(); ++cplxCpn) {
            DateTime accEnd = inst->cplxAccEndDates[cplxCpn];
            for (fundCpn=0; fundCpn<fundCfDates.size(); ++fundCpn) {
                if (inst->fundAccEndDates[fundCpn] > accEnd) {
                    cplxAddDates[cplxCpn] = min(cplxAddDates[cplxCpn], fundCfDates[fundCpn]);
                }
            }
        }

        for (fundCpn=0; fundCpn<fundCfDates.size(); ++fundCpn) {
            DateTime accEnd = inst->fundAccEndDates[fundCpn];
            for (cplxCpn=0; cplxCpn<cplxCfDates.size(); ++cplxCpn) {
                if (inst->cplxAccEndDates[cplxCpn] >= accEnd) {
                    fundAddDates[fundCpn] = min(fundAddDates[fundCpn], cplxCfDates[cplxCpn]);
                }
            }
        }
    }
    catch (exception& e){ 
        throw ModelException(e,"KTarnTree::KTarnTree"); 
    }
}

void KTarnTree::recordOutput(Control* ctrl, YieldCurveConstSP, Results* results) {
    if (inst->recordOutputName) {
        recordSliceToOutputName(ctrl, results, model, false, 
            "ComplexTarn", inst->discount->getCcy(), *cplxValue);
        recordSliceToOutputName(ctrl, results, model, false, 
            "FundingTarn", inst->discount->getCcy(), *fundValue);
        recordSliceToOutputName(ctrl, results, model, 
            inst->isTopComponent(), inst->outputName, inst->discount->getCcy(), getValue(0, model->getDate(0)));
    }
    inst->recordExtraOutput(ctrl, results);
}

void KTarnTree::init(Control*) const {
    ASSERT(cplxAddDates.size()==inst->cplxPay.size());
    for (int i=0; i<cplxAddDates.size(); ++i) {
        model->registerZero(cplxAddDates[i], inst->cplxPay[i], inst->discountYieldCurveName());
    }
}

/** initializing and setting product variables */
void KTarnTree::initProd(void){
    const string curve = inst->discount->getName();

    model->registerStateSupport(this);
    cashFlow = model->createSlice();
    cashFlow->name = inst->outputName+"_cashFlow";

    cplxValue = model->createSlice(curve);
    cplxValue->name = inst->outputName+"_cplxValue";
    *cplxValue = 0.0;
    startDEV(cplxValue);
    
    fundValue = model->createSlice(curve);
    fundValue->name = inst->outputName+"_fundValue";
    *fundValue = 0.0;
    startDEV(fundValue);

    value = model->createSlice(curve);
    value->name = inst->outputName+"_value";

    setGrid(model->getLastStep(), true);
}

void KTarnTree::setGrid(int step, bool init) {
    int nbGrid = inst->nbStateVar + (step%2) /* just to randomize it a bit for testing */;
    gridChanged = true;
    currGrid.resize(nbGrid);
    double loVal = inst->pastRealizedAmount;
    double hiVal = inst->targetLevel;
    if (Maths::isZero(hiVal-loVal)) 
        hiVal += 1;

    double slope = (hiVal - loVal) / (nbGrid - 1.0);
    for (int i=0; i<nbGrid; i++){ 
        currGrid[i] = slope * i + loVal;
    }

    if( init )
        prevGrid = currGrid;

    populateGrid( step );
}

void KTarnTree::transition(
    const vector< const TreeSlice * > & srcSlices,
    const vector< double > & gridIn,
    vector< double > & gridOut ) const
{
    int nbGrid = gridOut.size();
    ASSERT( nbGrid == (int)gridIn.size() && nbGrid == (int)srcSlices.size() );

    for( int i = 0; i < nbGrid; ++i )
        gridOut[ i ] = gridIn[ i ] + srcSlices[ i ]->calc();
}

void KTarnTree::update(int& step, UpdateType type) 
{
    DateTime stepDate = model->getDate(step);
    
    int cplxCpn=0, fundCpn=0;

    if (RatesUtils::happensNow(cplxCfDates, stepDate,  &cplxCpn)) {
        cplxLegCF[cplxCpn] = complexLegProd->getCashFlow(step).clone();
        // !!! model->cloneSlice( complexLegProd->getCashFlow(step), inst->discount->getName() );
        startDEV(cplxLegCF[cplxCpn]);
    }
    if (RatesUtils::happensNow(fundCfDates, stepDate,  &fundCpn)) {
        if (inst->knockInOut==KTarn::KnockType::KNOCK_OUT) {
            fundLegCF[fundCpn] = fundingLegProd->getCashFlow(step).clone();
            // !!! model->cloneSlice( fundingLegProd->getCashFlow(step), inst->discount->getName() );
        }
        else {
            fundLegCF[fundCpn] = fundingLegProd->getValue(step, stepDate).clone();
            // !!! model->cloneSlice( fundingLegProd->getValue(step, stepDate), inst->discount->getName() );
        }
        startDEV(fundLegCF[fundCpn]);
    }
    while (cplxCpn>=0 || fundCpn >=0) {

        for (fundCpn=fundAddDates.size()-1; fundCpn >=0 ; --fundCpn) {
            if (fundAddDates[fundCpn] == stepDate
            && fundLegCF[fundCpn].get())
                break;
        }

        for (cplxCpn=cplxAddDates.size()-1; cplxCpn >=0; --cplxCpn) {
            if (cplxAddDates[cplxCpn] == stepDate
            && cplxLegCF[cplxCpn].get())
                break;
        }
        if (cplxCpn>=0 && (fundCpn<0 || inst->fundAccEndDates[fundCpn]<=inst->cplxAccEndDates[cplxCpn]))
        {
            TreeSliceSP undCF = cplxLegCF[cplxCpn];
            
            if (inst->cplxPay[cplxCpn]!=stepDate) 
            {
                TreeSliceSP tmp;
                model->getZero(stepDate, inst->cplxPay[cplxCpn], inst->discountYieldCurveName(), tmp);
                *tmp = *undCF / *tmp;
                undCF = tmp;
            } 

            setGrid(step);

            const TreeSlice & stateVar = getGridSlice();
            if (inst->excessCoupon) 
            {
                *cashFlow = cond(stateVar < inst->targetLevel*(1-0.00001), *undCF, 0.);
            } 
            else {
                *cashFlow = smin(smax(inst->targetLevel - stateVar,0.), *undCF);
            }  
            if (cplxCpn == cplxAddDates.size()-1 && inst->finalPayment) 
            {
               *cashFlow += smax(inst->targetLevel - stateVar - *cashFlow, 0.);
            }

            const vector< TreeSliceSP > & slices =
                static_cast< const TreeSliceLayer & >( *cashFlow ).getSlices();
            int nbSlices = slices.size();
            vector< const TreeSlice * > srcSlices( nbSlices );
            for( int i = 0; i < nbSlices; ++i )
                srcSlices[ i ] = slices[ i ].get();
            interpolate( step, srcSlices );

            if (inst->knockInOut==KTarn::KnockType::KNOCK_IN) {
                *cashFlow = *undCF - *cashFlow;
            }
            {
                TreeSliceSP zero;
                model->getZero(stepDate, inst->cplxPay[cplxCpn], inst->discountYieldCurveName(), zero);
                *cplxValue += *cashFlow * *zero;
            }
            stopDEV(cplxLegCF[cplxCpn]);
            cplxLegCF[cplxCpn]=TreeSliceSP();
            continue;
        }
        if (fundCpn>=0) {
            // for efficiency don't need to populate slice again - it's already done
            const TreeSlice & stateVar = getGridSlice();
            //!!! REVIEW THIS: stateVar might have been populated at a different step
            const_cast< TreeSlice & >( stateVar ).treeStep = step;
            const vector< TreeSliceSP > slices =
                static_cast< const TreeSliceLayer & >( stateVar ).getSlices();
            int nbSlices = slices.size();
            for( int i = 0; i < nbSlices; ++i )
                slices[ i ]->treeStep = step;

            if (inst->knockInOut==KTarn::KnockType::KNOCK_OUT) {
                *fundValue += *fundLegCF[fundCpn];
                *fundValue = cond(stateVar < inst->targetLevel*(1-0.00001), *fundValue, 0.);
            }
            else *fundValue = cond(stateVar < inst->targetLevel*(1-0.00001), *fundValue, *fundLegCF[fundCpn]);

            stopDEV(fundLegCF[fundCpn]);
            fundLegCF[fundCpn]=TreeSliceSP();
        }
    }
}

/********************************** KTarn **********************************/
/** implement FDModel product interface */

void KTarn::setup(const IModel* model, const MarketData* market) {
    try {
        KComponent::setup(model, market);
        complexLeg->setup(model, market);
        fundingLeg->setup(model, market);
        
        {
            KComponentSP comp = KComponentSP::dynamicCast(complexLeg);
            comp->getAccEndDates(cplxAccEndDates);
            comp->getCfDates(cplxPay);
        }
        {
            KComponentSP comp = KComponentSP::dynamicCast(fundingLeg);
            comp->getAccEndDates(fundAccEndDates);
        }
    }
    catch (exception& e) {
        makeException(e,"setup");
    }
}

FDProductSP KTarn::createProduct(FDModel * model) const {
    if (!setupCalled) 
        throw ModelException(getClass()->getName()+" \""+outputName+"\" not setup");
        
    return FDProductSP(new KTarnTree(KTarnConstSP(this), model));
}

void KTarn::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(KTarn, clazz)
    SUPERCLASS(KComponent);
    IMPLEMENTS(FDModel::IIntoProduct);
    EMPTY_SHELL_METHOD(KTarn::defaultConstructor);
    FIELD(targetLevel, "");
    FIELD(knockInOut, ""); FIELD_MAKE_OPTIONAL(knockInOut);
    FIELD(finalPayment, "");       FIELD_MAKE_OPTIONAL(finalPayment);
    FIELD(payoffSmoothing, "");    FIELD_MAKE_OPTIONAL(payoffSmoothing);
    FIELD(pastRealizedAmount, ""); FIELD_MAKE_OPTIONAL(pastRealizedAmount);
    FIELD(excessCoupon, "");       FIELD_MAKE_OPTIONAL(excessCoupon);
    FIELD(nbStateVar, "");
    FIELD(stateVarNbStdDeviations, "");
    FIELD(complexLeg,"");
    FIELD(fundingLeg,"");

    /* transient fields */
    FIELD(cplxPay,"");         FIELD_MAKE_TRANSIENT(cplxPay);
    FIELD(fundAccEndDates,""); FIELD_MAKE_TRANSIENT(fundAccEndDates);
    FIELD(cplxAccEndDates,""); FIELD_MAKE_TRANSIENT(cplxAccEndDates);
    Addin::registerConstructor(Addin::UTILITIES, KTarn::TYPE);
}


CClassConstSP const KTarn::TYPE = CClass::registerClassLoadMethod(
    "KTarn", typeid(KTarn), KTarn::load);

/**********************************/
bool KTarnLoad()
{
    return KTarn::TYPE !=0;
}


START_PUBLIC_ENUM_DEFINITION(KTarn::KnockType::Enum, "");
ENUM_VALUE_AND_NAME(KTarn::KnockType::KNOCK_IN, "KNOCK_IN", "");
ENUM_VALUE_AND_NAME(KTarn::KnockType::KNOCK_OUT, "KNOCK_OUT", "");
END_ENUM_DEFINITION(KTarn::KnockType::Enum);

DRLIB_END_NAMESPACE
