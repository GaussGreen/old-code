//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : FDModel.cpp
//
//   Description : FD/tree model base class.
//
//   Author      : Ning Shen
//
//   Date        : February 8, 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FDModel.hpp"

#include "edginc/Instrument.hpp"
#include "edginc/BasketDelta.hpp"
#include "edginc/DeltaNextDay.hpp"
#include "edginc/StruckEquity.hpp"
#include "edginc/CrossGamma.hpp"
#include "edginc/FXCrossGamma.hpp"
#include "edginc/RollingTheta.hpp"
#include "edginc/Maths.hpp"
#include "edginc/IndexSpec.hpp"
#include "edginc/CorrelationCategory.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include <fstream>

DRLIB_BEGIN_NAMESPACE

/************* FDModel::IProdModifier *************/ 

void FDModel::IProdModifier::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER_INTERFACE(FDModel::IProdModifier, clazz);
    EXTENDS(IObject);
}

CClassConstSP const FDModel::IProdModifier::TYPE = CClass::registerInterfaceLoadMethod(
    "FDModel::IProdModifier", typeid(FDModel::IProdModifier), FDModel::IProdModifier::load);

typedef FDModel::IProdModifierArray FDModelIProdModifierArray;
DEFINE_TEMPLATE_TYPE_WITH_NAME("FDModel::IProdModifierArray", FDModelIProdModifierArray);

///********  FDModel ***********
// registration codes
// FDModel registration
CClassConstSP const FDModel::TYPE = CClass::registerClassLoadMethod(
    "FDModel", typeid(FDModel), load);

void FDModel::load(CClassSP& clazz){
    REGISTER(FDModel, clazz);
    SUPERCLASS(CModel);
    IMPLEMENTS(LastProductSensDate);

    FIELD(isFwdInduction, "Specifies whether to try using Forward instead of Backward PDE");
    FIELD_MAKE_OPTIONAL(isFwdInduction);
    FIELD(inductionType, "Specifies whether to try using Forward, Backward or Express Backward");
    FIELD_MAKE_OPTIONAL(inductionType);
    FIELD(productInfoFileName, "Name of file to print product debug/additional info to");
    FIELD_MAKE_OPTIONAL(productInfoFileName);
    FIELD(prodModifiers, "List of component modifiers");
    FIELD_MAKE_OPTIONAL(prodModifiers);
    FIELD(segDates,"");
    FIELD_MAKE_TRANSIENT(segDates);
    FIELD(density,"");
    FIELD_MAKE_TRANSIENT(density);
    FIELD(critDates,"");
    FIELD_MAKE_TRANSIENT(critDates);
    FIELD(stepsPerYearFD,"");
    FIELD_MAKE_TRANSIENT(stepsPerYearFD);
    FIELD(sameGridTweak,"");
    FIELD_MAKE_TRANSIENT(sameGridTweak);
    FIELD(sameGridDeltaShiftSize,"");
    FIELD_MAKE_TRANSIENT(sameGridDeltaShiftSize);


    FIELD(discYC,"")
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(discYC)
    FIELD(factors,"")
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(factors)
    FIELD(correls,"")
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(correls)
}

//QENUM_REGISTER_EMBEDDED(FDModel, InductionType);
START_PUBLIC_ENUM_DEFINITION(FDModel::InductionType::Enum, "");
ENUM_VALUE_AND_NAME(FDModel::InductionType::FWD, "FWD", "");
ENUM_VALUE_AND_NAME(FDModel::InductionType::BACK, "BACK", "");
ENUM_VALUE_AND_NAME(FDModel::InductionType::EXPRESSBACK, "EXPRESSBACK", "");
END_ENUM_DEFINITION(FDModel::InductionType::Enum);

// product interface registration
CClassConstSP const FDModel::IIntoProduct::TYPE =
CClass::registerInterfaceLoadMethod("FDModel::IIntoProduct",
                                    typeid(FDModel::IIntoProduct), load);


void FDModel::IIntoProduct::load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER_INTERFACE(FDModel::IIntoProduct, clazz);
        EXTENDS(Model::IModelIntoProduct);
        EXTENDS(IProdCreator);
}

FDModel::FDModel(const CClassConstSP &type) :
    CModel(type),
    stepsPerYearFD(50),
    sameGridTweak(false),
    isFwdInduction(false),
    inductionType(FDModel::InductionType::BACK),
    sameGridDeltaShiftSize(0.0),
    sliceCreationMode(true)
{}

//----------------------------------------------------------------
/** local helpers */
bool rebuildRequest(const SensitivitySP sens)
{
    bool rebuild = Theta::TYPE->isInstance(sens.get())
                || MuParallel::TYPE->isInstance(sens.get())
                || MuPointwise::TYPE->isInstance(sens.get())
                || MuSpecial::TYPE->isInstance(sens.get())
                || FXDelta::TYPE->isInstance(sens.get())
                || DeltaNextDay::TYPE->isInstance(sens.get())
                || DeltaSurface::TYPE->isInstance(sens.get())
                || RollingTheta::TYPE->isInstance(sens.get());

    return rebuild;
}

////////////////////////// FDModel /////////////////////////
//----------------------------------------------------------------

const FDProductArray & FDModel::getProducts() const
{
    return products;
}

/** add critical dates to initData for timeline setup */
void FDModel::addCritDates(const DateTimeArray& critDates)
{
    this->critDates.insert(this->critDates.end(), critDates.begin(), critDates.end() );
#ifdef DEBUG
    static int expectedI = 19980128; // set to desired value with the debugger (use "QuickWatch")
    DateTime expectedD = DateTime::fromIrDate(expectedI);
    for (int i=0; i<critDates.size(); ++i) if (critDates[i]==expectedD) {
        int dummy=0;
        ++dummy; // put breakpoint here
    }
#endif
}

void FDModel::addCritDate(const DateTime critDate)
{
    DateTimeArray tmp(1);
    tmp[0] = critDate;
    addCritDates(tmp);
}

/** add data for fd initialisation, timeline setup */
void FDModel::initSegments(const DateTimeArray& seg, 
                           const IntArray&      dens)
{
    segDates = seg;
    density = dens;
}

/**  helper function, aggregate segDates, critDates */
void FDModel::arrangeDates()
{
   // sort and remove duplicates
   TimeMetric::SortDate(true, segDates);
   TimeMetric::SortDate(true, critDates);
}



/**  loop through all the products in the product list and decide              *
 *   if forward induction is possible or if we need to store StatePrice points */
void FDModel::decideDirection() {
    FDProductArray::iterator it;
    DateTimeArray  statePriceDates;

    // initialise with neutral setting
    FDProduct::FDDirection portfolioDirection = FDProduct::FWD_BACK_INDUCTION; 
    inductionType = FDModel::InductionType::BACK;

    // loop through the products and check which FD induction they support
    for (it=products.begin(); it!=products.end(); ++it)
    {           
        FDProduct::FDDirection prodDirection = (*it)->getDirectionPref( );
    
        // we only need to do something if the products have different preferences
        if ( portfolioDirection != prodDirection)
        {
            if (portfolioDirection == FDProduct::FWD_BACK_INDUCTION)
            {
                portfolioDirection = prodDirection;
            }
            else if (prodDirection != FDProduct::FWD_BACK_INDUCTION)
            {
                inductionType = FDModel::InductionType::EXPRESSBACK;
            }
        }

        // collect the information for dates where state prices are required 
        DateTimeArray prodSPDates = (*it)->getStatePriceDates();

        statePriceDates = DateTime::merge( prodSPDates , statePriceDates );
    }

    // do forward induction if we can, otherwise use what has been set (BACK or EXPRESSBACK)
    if (portfolioDirection == FDProduct::FWD_INDUCTION)
        inductionType = FDModel::InductionType::FWD; // EXPRESSBACK ; // test.....should be FWD

    // store a list of possible dates to compute the state prices
    statePriceDates = DateTime::doSortUniq(statePriceDates); // 

    this->storeStatePriceDates( statePriceDates );
}




bool FDModel::isRebuilt(const CControl* control)
{
    bool buildTree = (!control || control->isPricing());
    if (!buildTree)
    {
        SensitivitySP sens = control->getCurrentSensitivity();
        if (Delta::TYPE->isInstance(sens.get())       ||
            BasketDelta::TYPE->isInstance(sens.get()) ||
            DeltaDDE::TYPE->isInstance(sens.get())    ||
            CrossGamma::TYPE->isInstance(sens.get())  ||
            FXCrossGamma::TYPE->isInstance(sens.get()))
        {
            sameGridTweak = true;
        }
        buildTree = rebuildRequest(sens);
    }
    return buildTree;
}

/** call init() on all elementary products first and then on all products */
void FDModel::initCall(CControl* control)
{
    critDates.clear(); 
    for( int i = 0; i < (int)products.size(); ++i )
        products[i]->init(control);
}

/** call initProd() on all elementary products first and then on all products */
void FDModel::initProdCall()
{
    for( int i = 0; i < (int)products.size(); ++i )
        products[i]->initProd();
}

/** print to stream any debug/other information each FDProduct chooses to do */
void FDModel::printProductInfo(ostream& outputStream)
{
    for( int i = 0; i < (int)products.size(); ++i )
    {
        outputStream << endl << "PRODUCT NUMBER " << i << endl;
        products[i]->printInfo(outputStream);
    }
}

/** call recordOutput() on all products, decorating results with instrument/component/indexSpec name:: */
void FDModel::recordOutputCall(Control* control, YieldCurveConstSP disc, Results* results)
{
    if (control && control->isPricing())
    {
        for( int i = 0; i < (int)products.size()-1; ++i ){ // top product must be the last one, not done here
//            products[i]->recordOutput(control, disc, results, true);
            products[i]->recordOutput(control, disc, results);
        }
    }
    // top product does not have decorator
    //prod->recordOutput(control, discYC, results, false);
    prod->recordOutput(control, discYC, results);
}

/** collect model initialisation data, set up timeline */
void FDModel::initModel(){
    try{
        arrangeDates();
        // create time line
        timeLine = TimeLineSP(new CTimeLine());
        timeLine->CreateTimeLineSimple(segDates,
                                     timeMetric,
                                     stepsPerYearFD,
                                     critDates);
    }
    catch(exception& e){
        throw ModelException(e, "FDModel::initModel()");
    }
}

void FDModel::price(const IProdCreatorSP & creator,
                    CControl* control, CResults* results, bool rebuild)
{
    try
    {
        // empty products and creators cache
        products.clear();
        creators.clear();

        // empty StateSupport cache
        stateSupportCache.clear();
        sliceCreationMode = true;

        // create the product
        prod = createProduct(creator);

        decideDirection(); 

        // further selecting market data
        retrieveFactor();

        // customize set up by product, rebuild flag allows same fd grid tweak
        if (rebuild) {
            initCall(control);
            initModel();    
        }
        // initialise product variables 
        initProdCall();

        // print out product hierarchy and any other info products choose to print
        if (!productInfoFileName.empty()) {
            ofstream outFile(productInfoFileName.c_str());
            if (!outFile.good())
                throw ModelException("Unable to open productInfo file name " + productInfoFileName);
            printProductInfo(outFile);
            outFile.close();
        }

        // finalise model
        finaliseModel(control);

        sliceCreationMode = false;

        // create FD solver
        IFDSolverSP solver(createSolver());
        solver->roll();

        if (control && control->isPricing()){
            // just an example to store some extra results
            results->storeScalarGreek(getLastStep(), Results::DEBUG_PACKET, 
                                      OutputNameSP(new OutputName("STEPS_USED")));

            //store any results from the model not captured through the product
            recordOutput(control, results);
        }

        // record price and additional outputs
        recordOutputCall(control, discYC, results);
    }
    catch (exception& e) {
        throw ModelException(e, "FDModel::price");
    }
}

/** main model entry point */
//----------------------------------------------------------------
void FDModel::Price(CInstrument* instrument, 
                    CControl*    control, 
                    CResults*    results){

    static const string method = "FDModel::Price";

    IIntoProduct * creator = dynamic_cast< IIntoProduct * >( instrument );
    if( ! creator )
    {
        throw ModelException(method,
            "Instrument of type "+ instrument->getClass()->getName() +
            " does not support FDModel::IIntoProduct");
    }

    try 
    {
        // check and price a dead instrument
        if (instrument->priceDeadInstrument(control, results))
            return; // done for a dead instrument

        // rebuild fd nodes when required
        sameGridTweak = true;
        bool rebuild = isRebuilt(control);
        if(!control || control->isPricing() || !rebuild)
            price(IProdCreatorSP::attachToRef(creator), control, results, rebuild);
        else
        {
            FDModelSP modelToUse( copy(this) ); 
            modelToUse->price(IProdCreatorSP::attachToRef(creator), control, results, rebuild);
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


// Pretty expensive way of determine last step date (used by endDate).
DateTime FDModel::lastStepDate(IIntoProduct &creator,
                               const Sensitivity *sensControl)
{
    try
    {
        // empty products and creators cache
        products.clear();
        creators.clear();

        // create the product
        prod = createProduct(IProdCreatorSP::attachToRef(&creator));

        // further selecting market data
        retrieveFactor();

        initCall(sensControl->getControl());
        initModel();

        return getDate(getLastStep());
    }
    catch (exception& e)
    {
        throw ModelException(e, "FDModel::lastStepDate failed.");
    }
}

/** Product provided end date for when to stop tweaking */
DateTime FDModel::endDate(const CInstrument *instrument,
                          const Sensitivity *sensControl) const
{
    static const string method = "FDModel::endDate";

    // Need to work with copies to preserve const
    FDModelSP modelToUse(copy(this));
    CInstrumentSP instToUse(copy(instrument));

    IIntoProduct *creator = dynamic_cast<IIntoProduct*>(instToUse.get());
    if (!creator)
    {
        throw ModelException(method,
            "Instrument of type "+ instrument->getClass()->getName() +
            " does not support FDModel::IIntoProduct");
    }

    return modelToUse->lastStepDate(*creator, sensControl);
}


/** override a control shift (eg for delta on trees)
    returns new control if constructed else returns 0 */
SensControl* FDModel::AlterControl(const SensControl* currSensControl) const
{
    if (!Delta::TYPE->isInstance(currSensControl) 
        || Maths::isZero(sameGridDeltaShiftSize /* this is not the best test */)){
        return 0;
    }

    //double deltaSize = control->getDeltaShiftSize();
    // default tweak size (0.5%) will give 5% of grid spacing
    //sameGridDeltaShiftSize = 10.0*deltaSize*(Stock[CurrIdx][1] - Stock[CurrIdx][-1])/Stock[CurrIdx][0];

    SensControlPerName* sensControl = new Delta(sameGridDeltaShiftSize);
    sensControl->setMarketDataName(currSensControl->getMarketDataName());

    return sensControl;
}

/** ask the model whether it is doing forward or backward induction */

bool FDModel::IsFwdInduction() const {
    return isFwdInduction;
}

/** ask the model whether it is doing forward, backward or express induction */

FDModel::InductionType::Enum FDModel::getInductionType() const {
    return inductionType;
}

void FDModel::getMarket(const MarketData* market, IInstrumentCollectionSP instruments)
{
    // recurse the instruments to retrieve domestic yield curve
    class RetrieveYC : public ObjectIteration::IAction
    {
        const string & name;
        YieldCurveConstSP & yc;
    public:
        RetrieveYC( const string & name, YieldCurveConstSP & yc ) :
            name( name ),
            yc( yc )
        {}

        virtual bool invoke(ObjectIteration::State& state, IObjectConstSP obj)
        {
            const YieldCurveConstSP & ryc = YieldCurveConstSP::dynamicCast( obj );
            // set yc only if it's empty yet
            if( ! yc && name == ryc->getName() )
                yc = ryc;

            // don't recurse inside the YieldCurve
            return false;
        }
    };
    RetrieveYC retrieveYC( getDomesticYCName(), discYC );
    ObjectIteration YCIter( YieldCurve::TYPE );
    YCIter.setSkipTransient( true );
    YCIter.recurse( retrieveYC, instruments );

    // recurse the instruments to retrieve factors
    class RetrieveFactors : public ObjectIteration::IAction
    {
        FDModel * model;
        const MarketData * market;
        IMarketFactorArray & factors;

    public:
        RetrieveFactors( FDModel * model, const MarketData * market, IMarketFactorArray & factors ) :
            model( model ), market( market ), factors( factors ) { factors.clear(); }

        virtual bool invoke(ObjectIteration::State& state, IObjectConstSP obj)
        {
            IMarketFactor * factor = dynamic_cast< IMarketFactor * >( state.getObject().get() );

            // don't collect and recurse into duplicates
            string name = factor->getName();
            string type = factor->getClass()->getName();
            for( int i = 0; i < factors.size(); ++i )
            {
                if( name == factors[ i ]->getName() && type == factors[ i ]->getClass()->getName() )
                     return false;
            }

            // if factor is accepted by derived model
            if( model->acceptFactor( factor ) )
            {
                factors.push_back( IMarketFactorSP::attachToRef( factor ) );
                // don't recurse inside
                return false;
            }
            else
                // continue recursing
                return true;
        }
    };
    RetrieveFactors retrieveFactors( this, market, factors );
    ObjectIteration factorsIter( IMarketFactor::TYPE );
    factorsIter.setSkipTransient( true );
    factorsIter.recurse( retrieveFactors, instruments );

    // retrieve correlations
    int numFactors = factors.size();
    correls = CDoubleMatrixSP( new DoubleMatrix( numFactors, numFactors ) );
    for( int i = 0; i < numFactors; ++i )
    {
        (*correls)[i][i] = 1.; // set diagonal element
        for( int j = i + 1; j < numFactors; ++j )
        {
            string corrName;
            if ( market->hasCorrelationData( factors[i]->getName(), factors[j]->getName() ) )
            {
                // get hold of the actual object
                CorrelationBaseSP corrBase = getCorrelation(
                    market->getCorrelationName(
                        factors[i]->getName(),
                        factors[j]->getName() ),
                    factors[i]->getClass(),
                    factors[j]->getClass(),
                    Correlation::TYPE,
                    market );
                (*correls)[i][j] = (*correls)[j][i] =
                    CorrelationSP::dynamicCast( corrBase )->getCorrelation();
            } // The following would take up too much memory for sampras (couple of gigs).
            else {
                (*correls)[i][j] = (*correls)[j][i] = 0.;
            }
        }
    }
}

void FDModel::registerStateSupport( TreeSliceLayer::StateSupport * owner )
{
    stateSupportCache.push_back( owner );
    owner->init( this, stateSupportCache.size() - 1 );
}

// helper function for createSlice
TreeSliceSP FDModel::createLayer( const TreeSliceSP & slice ) const
{
    TreeSliceSP s = slice;

    if (!sliceCreationMode 
    && stateSupportCache.size() /* for backward compatibility */) {
        throw ModelException("FDModel::createLayer","You can create slices only from initProd");
        /* otherwise you might create a TreeSliceLayer in a product that does not
         * handle TreeSliceLayers */
    }

    for( int i = 0; i < (int)stateSupportCache.size(); ++i )
        s = TreeSliceSP( new TreeSliceLayer( stateSupportCache[ i ], *s, true ) );
    return s;
}

/** product creation request is always done through this function
    it takes topmost base class (IProdCreator not FDModel::IIntoProduct)
    to avoid down-casting before each call */
FDProductSP FDModel::createProduct( const IProdCreatorSP & creator )
{
    // creators are collected in creation order
    // and checked for duplicates
    vector< IProdCreatorSP >::iterator iter = find( creators.begin(), creators.end(), creator );
    if( iter == creators.end() )
    {
        creators.push_back( creator );

        int idx = prodIdxs.size();
        prodIdxs.push_back( -1 );

        // create the product
        FDProductSP product = makeProduct( creator );

        // products are collected in dependency order
        // and checked for duplicates
        FDProductArray::iterator iter = find( products.begin(), products.end(), product );
        if( iter == products.end() )
        {
            prodIdxs[ idx ] = products.size();
            products.push_back( product );
        }
        else
            prodIdxs[ idx ] = iter - products.begin();

        return product;
    }
    else
        return products[ prodIdxs[ iter - creators.begin() ] ];
}

/** creates actual products
    can be overridden by inherited models */
FDProductSP FDModel::makeProduct(const IProdCreatorSP & creator)
{
    if (!creator)
        throw ModelException("FDModel::makeProduct", "Product creator is a NULL object pointer"
            " - has not been defined.  Parent component must point to valid object ");

    const FDModel::IIntoProduct * intoProduct = dynamic_cast< const FDModel::IIntoProduct * >( creator.get() );
    if( ! intoProduct )
    {
        throw ModelException( "FDModel::makeProduct",
            creator->getClass()->getName() + " does not support FDModel::IIntoProduct" );
        
    }

    return intoProduct->createProduct( this );
}

FDModel::IProdModifierSP FDModel::getProdModifier(CClassConstSP const &type) {
    int i,  nb = prodModifiers.size();
    for (i=0; i<nb; ++i) {
        IProdModifierSP pm = prodModifiers[i];
        if (type->isInstance(pm))
            return pm;
    }
    return IProdModifierSP();
}

// type loading
bool FDModelLoad(){
    return (FDModel::TYPE && FDModel::IIntoProduct::TYPE);
}

DRLIB_END_NAMESPACE
