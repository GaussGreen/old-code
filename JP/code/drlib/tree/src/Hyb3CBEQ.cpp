//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : Hyb3CBEQ.hpp
//
//   Description : hyb3 domestic equity mode FD model class
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Hyb3CBEQ.hpp"
#include "edginc/TreeSlice.hpp"
#include "edginc/IRVolPair.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/IrConverter.hpp"
#include "edginc/IndexSpecEQ.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SimpleEquity.hpp"
#include "edginc/VolSurface.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/SRMEQVol.hpp"
#include "edginc/FixedSettlement.hpp"
#include "edginc/RollingSettlement.hpp"
#include "edginc/IndexSpecFX.hpp"
#include "esl_log.h"
#include "esl_market.h"


DRLIB_BEGIN_NAMESPACE


static int toDivType(int t) {
    switch (t) {
        case 0: return 'D'; // AMOUNT
        case 1: return 'Y'; // PERCENT
        case 2: return 'C'; // CONTINUOUS
        default: ASSERT(0); return -1;
    }
}

 // ??? reuse as much as existing for now - rewrite later   
 static void PopulateEQData(
    EQ_DATA    *eq_data,
    long       today,
    double     *eqSpotVolOverride,
    const SimpleEquity *eq,
    double     *EQCalibrationStrike,
    DividendListConstSP &dividendList)
{
    char *routine = "PopulateEQData";
    int i;
    
    if (eqSpotVolOverride==NULL)    /* 1.Disallow cut off      */
    {
        eq_data->EqCutOffFlag = FALSE;
        eq_data->EqCutOffLast = FALSE;
        eq_data->EqCutOffLevel   = EQSPOT_CUTOFF_R; /* Ratio to fwd vol level*/
        eq_data->EqBootStrapMode = EQ_NO_FAILURE_ALLOWED;
    }
    else                                         /* 3.Cut off at user level */
    {
        /* Cut off is allowed and there */
        /* is volatility cutoff data    */
        eq_data->EqCutOffFlag = TRUE;
        eq_data->EqCutOffLast = FALSE;
        eq_data->EqCutOffLevel = *eqSpotVolOverride;

        if (eq_data->EqCutOffLevel < TINY) /* The FAILURE return */
            throw ModelException(routine, "Unable to convert EQ cut off value!\n");

        eq_data->EqBootStrapMode = EQ_CONSTANT_SPOT_VOL;     
    }

    if (eq)
    {
        IrConverter::to_EQ_DATA(*eq_data, *eq);
        eq_data->ValueDate = today;

        // if no smile parameters in the environment, throw an error for now - sort out
        // properly when smile families are implemented
        if (eq_data->NbSmilePt == 0)
            throw ModelException("At least one equity smile date and a1,a2,a3 data set must "
                                 "be provided.  To turn smile off set a1=a2=a3 = 0.0");

        eq_data->Spot = eq->getSpot();

        // ??? Need to use a different vol representation, but not sure what
        const CVolBaseWrapper& eqVol = eq->getVol();
        const SRMEQVol *vol = dynamic_cast<const SRMEQVol*>(eqVol.get());
        if (!vol) 
            throw ModelException(routine, "Equity vol market structure must be of type SRM3EQVol");

        if (EQCalibrationStrike)
            throw ModelException(routine, "EQCalibrationStrike field not yet supported");
        /*
        LinearStrikeVolRequest volRequest(EQCalibrationStrike, 
                                            baseVolDates[0],
                                            baseVolDates.back(),
                                            false);

        DoubleArray vols(baseVolDates.size());
        CVolProcessedBS::calcVol(vol,
                                    &volRequest,
                                    eq,
                                    baseVolDates[0],
                                    baseVolDates,
                                    vols);
        */

        DateTimeArray compVolDates = vol->getCompVolDate();
        DoubleArray compVols = vol->getCompVol();

        eq_data->NbVol = compVolDates.size();
        for (i = 1; i <= eq_data->NbVol; i++)
        {
            eq_data->VolDate[i] = compVolDates[i-1].toIrDate();
            eq_data->Vol[i] = compVols[i-1];;
        }
        if (eq_data->NbVol > 0)
        {
            eq_data->VolDate[0] = today;
            eq_data->Vol[0] = eq_data->Vol[1];
        }

        DateTimeArray spotVolDates = vol->getSpotVolDate();
        DoubleArray spotVols = vol->getSpotVol();

        eq_data->NbInpSpotVol = spotVolDates.size();
        for (i = 0; i < eq_data->NbInpSpotVol; i++)
        {
            eq_data->InpSpotVolDate[i] = spotVolDates[i].toIrDate();
            eq_data->InpSpotVol[i] = spotVols[i];
        }

        // populate "static" equity data
        EquitySP equity = eq->getEquity();

        BorrowCurveSP borrow = equity->getBorrow();
        RhoBorrowPointwise rhoBorrow(0.0);  // ??? dummy sens
        ExpiryArrayConstSP borrowDates = equity->sensExpiries(&rhoBorrow);

        dividendList = equity->getDivList();
        const DividendArray& divs = dividendList->getArray();

        eq_data->NbBorrow = (*borrowDates).size();
        if (eq_data->NbBorrow > MAXNBDATE)
            throw ModelException(routine, "Nb of borrow curve points " +
                                Format::toString(eq_data->NbBorrow) + 
                                " exceeds maximum internal tree limit of " +
                                Format::toString(MAXNBDATE));

        // ??? have to check this properly
        for (i = 0; i < eq_data->NbBorrow; i++)
        {
            DateTime baseDate = equity->getValueDate();  // ??? no idea if this is correct date
            eq_data->BorrowDate[i] = (*borrowDates)[i]->toDate(baseDate).toIrDate();
            eq_data->Borrow[i] = borrow->interpAtDate(DateTime::fromIrDate(eq_data->BorrowDate[i]));
        }

        eq_data->NbFwd = divs.size();  // ??? need to check what dates in equity.sta are exactly
        if (eq_data->NbFwd > MAXNBEQDIV)
            throw ModelException(routine, "Number of dividends " +
                                Format::toString(eq_data->NbFwd) +
                                " exceeds maximum internal tree limit of " +
                                Format::toString(MAXNBEQDIV));
      
        // ??? first attempt at getting something to work
        for (i = 0; i < eq_data->NbFwd; i++)
        {
            eq_data->FwdDate[i] = divs[i].getExDate().toIrDate();
            eq_data->Fwd[i] = divs[i].getDivAmount();
            eq_data->FwdType[i] = toDivType(divs[i].getDivType()); 
        }

        const Settlement& settlement = equity->getSettlement();
        eq_data->SettleType = '0';
        {
            const FixedSettlement *fs = dynamic_cast<const FixedSettlement*>(&settlement);
            if (fs) {
                eq_data->SettleType = 'F';
                eq_data->NbSettle = fs->getFirstTradingDates().size();
                for (int i=0; i<eq_data->NbSettle; ++i) {
                    eq_data->LastTrading[i] = fs->getFirstTradingDates()[i].toIrDate();
                    eq_data->SettleDate[i] = fs->getSettlementDates()[i].toIrDate();
                }
            }
        }
        {
            const RollingSettlement *rs = dynamic_cast<const RollingSettlement*>(&settlement);
            if (rs) {
                eq_data->SettleType = 'R';
                eq_data->NbSettle = rs->getPeriod();
            }
        }
        if (eq_data->SettleType == '0') 
            throw ModelException("Unsupported settlement type "+settlement.getClass()->getName());
    }
    else
    {
        // if equity not used in product, set default values
        long dummyDate = Nxtday(today, 100);  // dummy first vol date

        eq_data->Spot = 0.0;
        eq_data->NbVol = 1;
        eq_data->ValueDate = today;
        eq_data->VolDate[0] = today;
        eq_data->VolDate[1] = dummyDate;
        eq_data->Vol[0] = 0.1;
        eq_data->Vol[1] = 0.1;
        eq_data->NbInpSpotVol = 0;
        eq_data->InpSpotVolDate[0] = today;
        eq_data->InpSpotVol[0] = 0.0;

        eq_data->NbBorrow = 1;
        eq_data->BorrowDate[0] = dummyDate;
        eq_data->Borrow[0] = 0.0;

        eq_data->NbFwd = 1;
        eq_data->FwdDate[0] = dummyDate;
        eq_data->Fwd[0] = 0.0;
        eq_data->FwdType[0] = 'D';

        eq_data->NbSettle = 0;
        eq_data->SettleType = 'R';  // ??? other option is fixed 'F'

        eq_data->NbSmilePt = 1;
        eq_data->SmileDate[0] = dummyDate;
        eq_data->a1[0] = 0.0;
        eq_data->a2[0] = 0.0;
        eq_data->a3[0] = 0.0;
    }

    // Check validity of input
    if (Hyb3_Eq_Check_W(eq_data) != SUCCESS)
        throw ModelException(routine, IrConverter::slog.pop());

    if (Hyb3_Eq_Check_Sta_W(eq_data) != SUCCESS)
        throw ModelException(routine, IrConverter::slog.pop());
}


 /***********************************
 *
 * Hyb3 Domestic Equity tree mode
 *
************************************/


void Hyb3CBEQ::load(CClassSP& clazz) {

    clazz->setPublic();
    REGISTER(Hyb3CBEQ, clazz);
    SUPERCLASS(Hyb3CB);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(ppy, "minimum number of tree nodes(points)/year");
    FIELD(maxStdDeviations, "number of standard deviations to trim the tree at");
    FIELD(IRParams, "Collection of interest rate model parameters");
    FIELD(EQCalibrationStrike, "Strike used to determine volatility matrix column "
        "for tree EQ vol calibration");
    FIELD_MAKE_OPTIONAL(EQCalibrationStrike);
    FIELD(lastDividendDate,"");
    FIELD_MAKE_OPTIONAL(lastDividendDate);
    FIELD(eqSpotVolOverride,"Overrides market comp/spot vols for equity - sets constant "
        "spot vol for EQ diffusion (units = decimal NOT percentage)");
    FIELD_MAKE_OPTIONAL(eqSpotVolOverride);
    FIELD(IREQCorrelationOverride, "Overrides the market defined IR/Equtiy correlation "
        "value");
    FIELD_MAKE_OPTIONAL(IREQCorrelationOverride);
    
    // ensure fields are copied when tree is copied
    FIELD(corrIREQ,"");
    FIELD_MAKE_TRANSIENT(corrIREQ);
    FIELD(eqName,"");
    FIELD_MAKE_TRANSIENT(eqName);
}

CClassConstSP const Hyb3CBEQ::TYPE = CClass::registerClassLoadMethod(
    "Hyb3CBEQ", typeid(Hyb3CBEQ), Hyb3CBEQ::load);

// for class loading
bool Hyb3CBEQLoad(void) { return (Hyb3CBEQ::TYPE != 0); }

Hyb3CBEQ::Hyb3CBEQ() : Hyb3CB(TYPE), EQCalibrationStrike(0), eqSpotVolOverride(0),  
    IREQCorrelationOverride(0), corrIREQ(0.0) 
{}

Hyb3CBEQ::~Hyb3CBEQ() {
    clear();
}

// reset the underlying structures to NULL etc.
void Hyb3CBEQ::reset() {

    // initialise structures
    MktVol_Init(&mktVolData[0]);  // ??? AK: DOMESTIC is defined as "1"....needs change in RateTree

    Hyb3_Tree_Init(&treeData);
    Hyb3_Dev_Init(&devData);

    mTpIdx = -1;
    //    mSmoothFlag = 0;  ???

    // initialiase the number of factors of the rates-modes
    mktVolData[0].NbFactor = IRParams->nbFactors;

    treeBuilt = false;
}



void Hyb3CBEQ::clear() {

    if (!treeBuilt)
        return;
    treeBuilt = false;

    int i;

    map<string, NamedClaimBank>::iterator cb;
    for (cb = claimBank.begin(); cb != claimBank.end(); cb++) {
        Hyb3_CbkFree(&cb->second.zeroBank, &treeData, 1);
        Hyb3_CbkInit(&cb->second.zeroBank);
    }

    Hyb3_Dev_Free(&devData, &treeData);
    Hyb3_Tree_Free(&treeData);

    for (i=0; i<NBCRV; ++i)
        today2ValDateZeros[0][i] = -1.0;

    reset(); // set all the data structures to the initial state

    mTpIdx = -1;
}

// Implementation of IRVegaPointwise::ISensitivePoints for smart vega tweaking
// and volatility exposure reporting.
IRGridPointAbsArraySP Hyb3CBEQ::getSensitiveIRVolPoints(OutputNameConstSP  outputName, 
                                                      const CInstrument* inst) const {
    try {
        IRGridPointAbsArraySP volExposuresSP(new IRGridPointAbsArray);

        T_CURVE diffusionCurve;
        IrConverter::to_T_CURVE(diffusionCurve, IRParams->curveToDiffuse.get(),
                                false);

        IRVolSelector volSelector(getIRVolRaw(IRParams), diffusionCurve,
                                IRParams->volCalibIndex,
                                cetSmoothing, IRParams);
        volSelector.getVolExposures(volExposuresSP, outputName);
     
        return volExposuresSP;
    }
    catch (exception& e){
        throw ModelException(e, __FUNCTION__);
    }
}


void Hyb3CBEQ::getMarket(const MarketData*  market, IInstrumentCollectionSP instruments) {
    try {
        today = market->GetReferenceDate();
        // ??? need final decision by business what to do here
        // set valueDate = today for now, but today may be SOD or EOD but we
        // want valueDate to be EOD no matter what time today is - this means
        // payments are always dropped anytime today
        valueDate = DateTime(today.getDate(), DateTime::END_OF_DAY_TIME);

        IRParams->getData(this, market);

        // Get Engine
        if (!engineTable.isEmpty())
            engineTable.getData(this, market);

        // recurse the instruments to retrieve domestic yield curve
        class RetrieveYC : public ObjectIteration::IAction {
            const string & name;
            YieldCurveConstSP & yc;
        public:
            RetrieveYC( const string & name, YieldCurveConstSP & yc ) :
                name( name ),
                yc( yc )
            {}

            virtual bool invoke(ObjectIteration::State& state, IObjectConstSP obj) {
                const YieldCurveConstSP & ryc = YieldCurveConstSP::dynamicCast( obj );
                if( ! yc && name == ryc->getName() )
                    yc = ryc;

                // don't recurse inside the YieldCurve
                return false;
            }
        };
        RetrieveYC retrieveYC( getDomesticYCName(), discYC );
        ObjectIteration retrieveYCIter( YieldCurve::TYPE );
        retrieveYCIter.recurse( retrieveYC, instruments );

        // retrieve domestic YieldCurve
        if(!discYC)
            throw ModelException("Internal error - Product must populate FDModel::discYC field");

        if(discYC->getName().empty())
            throw ModelException("Internal error - Product must populate FDModel::discYC field");
        
        // in EQ mode, must provide curve to discount
        if (!IRParams->curveToDiscount.isEmpty())
            if (IRParams->curveToDiscount->getName() != discYC->getName())
                throw ModelException("Supplied model curveToDiscount " + 
                                    IRParams->curveToDiscount->getName() +
                                    " must be same as instruments discount curve " +
                                    discYC->getName());


        // recurse the components to let IndexSpec::getMarket
        class RetrieveFactors : public ObjectIteration::IAction {
            FDModel * model;
            const MarketData * market;
            IMarketFactorArray & factors;

        public:
            RetrieveFactors( FDModel * model, const MarketData * market, IMarketFactorArray & factors ) :
                model( model ), market( market ), factors( factors ) { factors.clear(); }

            virtual bool invoke(ObjectIteration::State& state, IObjectConstSP obj) {

                IMarketFactor * factor = dynamic_cast< IMarketFactor * >( state.getObject().get() );

                string name = factor->getName();
                string type = factor->getClass()->getName();
                int i;
                for( i = 0; i < factors.size(); ++i ) {
                    if( name == factors[ i ]->getName() && type == factors[ i ]->getClass()->getName() )
                        break;
                }
                if( i >= factors.size() && model->acceptFactor(factor)) {
                    factors.push_back( IMarketFactorSP::attachToRef( factor ) );
                    return false;
                }
                else
                    return true;
            }
        };
        RetrieveFactors retrieveFactors( this, market, factors );
        ObjectIteration iteration( IMarketFactor::TYPE );
        iteration.recurse( retrieveFactors, instruments );

        // check equity factor was supplied
        int i;
        
        int eqIndex = -1;
        for (i = 0; i < factors.size(); i++) {

            IMarketFactorConstSP &mf = factors[i];
            const SimpleEquity* eq = dynamic_cast<const SimpleEquity*>(mf.get());
            const IndexSpecEQ* eqSpec = dynamic_cast<const IndexSpecEQ*>(mf.get());

            if (eq) {
                if (eqIndex == -1) {
                    eqIndex = i;
                    eqName = eq->getName();
                }
                else {
                    if (eq->getName() != eqName)
                        throw ModelException("Hyb3CBEQ only supports ONE equity - at least two equities "
                                            "supplied (" + eq->getName() + ", " +
                                            factors[eqIndex]->getName() + ")");
                }
            }
            // ??? serious issues here - indexSpecEQ should be registering market factors but doesn't seem to,
            // so I will do it here for now - indexSpecEQ getting complicated and don't want to mess change anything
            // now in case it breaks something else
            if (eqSpec) {
                if (eqIndex == -1) {
                    const IMarketFactorConstSP &mfEQ = eqSpec->getFactor();
                    IMarketFactorSP nonConstmfEQ(const_cast<IMarketFactor*>(mfEQ.get()));
                    string name = mfEQ->getName();
                    string type = mfEQ->getClass()->getName();

                    const SimpleEquity* eqEQ = dynamic_cast<const SimpleEquity*>(mfEQ.get());
                    if (!eqEQ)
                        throw ModelException("Hyb3CBEQ currently only supports SimpleEquity equity type");
                    eqName = eqEQ->getName();
                    // for now, add this simple equity as a marketFactor so it can be more easily retrieved
                    // in the retrieve factors function
                    factors.push_back(nonConstmfEQ);
                }
            }
        }

        // recurse the components to store specific index specs 
        class RetrieveIndexSpecs : public ObjectIteration::IAction {
            IndexSpecArray & indexSpecs;

        public:
            RetrieveIndexSpecs(IndexSpecArray & indexSpecs) :
                indexSpecs( indexSpecs ) {}

            virtual bool invoke(ObjectIteration::State& state, IObjectConstSP obj) {
                int i;
                IndexSpec * indexSpec = dynamic_cast< IndexSpec * >( state.getObject().get() );

                string name = indexSpec->getName();
                string type = indexSpec->getClass()->getName();

                for( i = 0; i < indexSpecs.size(); ++i ) {
                    if(name == indexSpecs[ i ]->getName() && 
                    type == indexSpecs[ i ]->getClass()->getName() )
                    break;
                }
                if(i >= indexSpecs.size()) {
                    indexSpecs.push_back( IndexSpecSP::attachToRef( indexSpec ) );
                    return false;
                }
                else
                    return true;
            }
        };
        RetrieveIndexSpecs retrieveIS ( indexSpecs );
        ObjectIteration retrieveISIter( IndexSpec::TYPE );
        retrieveISIter.recurse( retrieveIS, instruments );

        // retrieve correlation for equity and IR
        bool corrFound = false;

        // check correlation of equity for all IR curves provided
        if (market->hasCorrelationData(IRParams->curveToDiffuse->getName(), eqName)) {

            string corrName = market->getCorrelationName(IRParams->curveToDiffuse->getName(), eqName);
            CorrelationBaseSP corrBase = getCorrelation(corrName,
                                                        IRParams->curveToDiffuse->getClass(),
                                                        SimpleEquity::TYPE,
                                                        Correlation::TYPE,
                                                        market);
            corrFound = true;
            corrIREQ = CorrelationSP::dynamicCast(corrBase)->getCorrelation();
        }

        // check if eq/(ir discount) correlation supplied
        if (market->hasCorrelationData(discYC->getName(), eqName)) {

            string corrName = market->getCorrelationName(discYC->getName(), eqName);
            CorrelationBaseSP corrBase = getCorrelation(corrName,
                                                        discYC->getClass(),
                                                        SimpleEquity::TYPE,
                                                        Correlation::TYPE,
                                                        market);
            double result = CorrelationSP::dynamicCast(corrBase)->getCorrelation();

            if (corrFound) {  //  check same value
                if (!Maths::equals(corrIREQ, result))
                    throw ModelException("Correlation found for discount yieldCurve " +
                                        discYC->getName() + " and equity " + eqName + "(value = " + 
                                        Format::toString(result) +
                                        ") does not equal previously found correlation for "
                                        "curveToDiffuse " + IRParams->curveToDiffuse->getName() +
                                        " and equity " + Format::toString(corrIREQ));
            }
            else {
                corrFound = true;
                corrIREQ = CorrelationSP::dynamicCast(corrBase)->getCorrelation();
            }
        }

        // ??? maybe check 3rd curve/Eq correlation if supplied but don't know if that curve
        // has been supplied until retrieve factor - think about refactoring the whole thing later
        if (eqName.empty())
            corrIREQ = 0.0;
        else
            if (corrFound == false)
                throw ModelException("Correlation not found in market for IR and equity "
                                    "searched for the following 2 correlations in the market: " + 
                                    market->getCorrelationName(IRParams->curveToDiffuse->getName(), eqName) +
                                    " and " +
                                    market->getCorrelationName(discYC->getName(), eqName));
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void Hyb3CBEQ::initModel() {
    try {
        treeBuilt = true; 

        CRIT_DATE* critDate=NULL;
        int nbCritDate = 0;
        IrConverter::AutoDeleteDRArray<CRIT_DATE>
            critDateDelete(&critDate, CRITDATE, &nbCritDate);

        // ??? note:  ZeroInterpTypeFlag is global variable defined in both hyb3.lib
        // and fix3.lib - tidy up but have to check all hyb3 products
        EslSetZeroInterpolation(zeroInterpStyle == ZeroInterpStyle::LINEAR ?
            ESL_INTERP_LINEAR : ESL_INTERP_FLATFWD);

        int j;

        map<string, NamedClaimBank>::iterator cb;

        // needed to setup timeline
        critDate = (CRIT_DATE *) DR_Array (CRITDATE, 0, 0);
        if (critDate == NULL)
            throw ModelException("Unable to allocate memory for critDate array");

        // today always critical date
        if (Add_To_DateList(&nbCritDate, 
                            &critDate,
                            today.toIrDate(),
                            0, 0, 0, 0, 0, 0, 0, 0, 0) != SUCCESS)
            throw ModelException("Add_To_DateList function failed : " + IrConverter::slog.pop());

        // valueDate always critical date
        if (Add_To_DateList(&nbCritDate, 
                            &critDate,
                            valueDate.toIrDate(),
                            0, 0, 0, 0, 0, 0, 0, 0, 0) != SUCCESS)
            throw ModelException("Add_To_DateList function failed : " + IrConverter::slog.pop());

        // add critical dates registered directly by products
        for (j = 0; j < critDates.size(); j++) {

            DateTime critDateL = critDates[j];
            if (critDateL > valueDate) {
                if (Add_To_DateList(&nbCritDate, 
                                    &critDate,
                                    critDateL.toIrDate(),
                                    0, 0, 0, 0, 0, 0, 0, 0, 0) != SUCCESS)
                    throw ModelException("Add_To_DateList function failed : " 
                                         + IrConverter::slog.pop());
            }
        }

        // store zero/discount factor between today and the value date to
        // in order to forward the price back to the value date
        for (j=0; j<NBCRV; j++) {
            today2ValDateZeros[0][j] = ::pow(1.0 + treeCurves[0][j].Zero[0], 
                -Daysact(today.toIrDate(), currencyValueDate[0].toIrDate())/365.0);

            // adjust the zero dates from the value to today by scaling all of 
            // the zeroRates in the TCurves - domestic eq mode only uses index 0 of TCurve
            RateTree::ExtendTreeCurve(&treeCurves[0][j], today);
        }

        // configure and add claim banks
        configureClaimBank(critDate, nbCritDate);

        // add claim bank critical dates to list
        for (cb = claimBank.begin(); cb != claimBank.end(); cb++) {

            // add combined critical zero bank dates
            int nbCrit = cb->second.critZeroUseDates.size();

            for (j=0; j<nbCrit; ++j) {
                long matDate = cb->second.critZeroMatDates[j];
                long useDate = cb->second.critZeroUseDates[j];
                if (Add_To_DateList(&nbCritDate,
                                    &critDate, 
                                    matDate,
                                    0, 0, 0, 0, 
                                    0, 0, 0, 0, 0) != SUCCESS)
                    throw ModelException("Add_To_DateList function failed : " 
                                         + IrConverter::slog.pop());

                if (Add_To_DateList(&nbCritDate, 
                                    &critDate,
                                    useDate,
                                    0,0, 0, 0, 
                                    0, 0, 0, 0, 0) != SUCCESS)
                    throw ModelException("Add_To_DateList function failed : " 
                                         + IrConverter::slog.pop());
            }
        }
    
        // ??? review
        if (dividendList.get()) {

            const DividendArray& divs = dividendList->getArray();
            for (j = 0; j < divs.size(); j++) {

                DateTime exDate = divs[j].getExDate();
                if (exDate.getDate() >= getValueDate().getDate() /* exDate occurs @time XBS, not SOD */
                    &&  (lastDividendDate.empty() || exDate <= lastDividendDate)) {

                    if (Add_To_DateList (&nbCritDate,
                                        &critDate,
                                        exDate.toIrDate(),
                                        12,
                                        divs[j].getDivAmount(),
                                        0, 0, 0, 0,
                                        0, 0, 0) != SUCCESS) {
                        throw ModelException("Add_To_DateList function failed : " 
                                             + IrConverter::slog.pop());
                    }
                }
            }
        }

        // remove duplicates, sort critical dates and store them to
        // help debugging by printing them out in the debug data file
        if (Sort_CritDate(nbCritDate, critDate) != SUCCESS)
            throw ModelException("Sort_CritDate failed: " + IrConverter::slog.pop());

        sortedCritDates.resize(1);
        sortedCritDates[0] = critDate[0].CritDate;
        int i;
        for (i = 1; i < nbCritDate; i++) {
            if (sortedCritDates.back() == critDate[i].CritDate)
                continue;
            sortedCritDates.push_back(critDate[i].CritDate);
        }
        if (sortedCritDates.size()<=1) {
            range.reset(new TreeSliceRates::Range(0,0,0,0,0,0));
            treeData.NbTP = -1;
            return;
        }
        // construct the timeline 
        if (Hyb3_Time_Line(today.toIrDate(), nbCritDate, critDate, 'I', &treeData) 
            != SUCCESS)
            throw ModelException("Hyb3_Time_Line failed: " + IrConverter::slog.pop());

        // if printing out a debug file, turn on printing the default CET debug file also
        int printCET = FALSE;
        if (!treeDataDebugFile.empty())
            printCET = TRUE;
        if (Hyb3_Build_Tree(printCET, treeCurves, mktVolData, &fxData, &eqData, 
                            &treeData) != SUCCESS)
            throw ModelException("Hyb3_Build_Tree failed: " + IrConverter::slog.pop());

        // Equity 1D mode is always 2D tree
        range.reset(new TreeSliceRates::Range(
                    -treeData.HalfWidth[0],
                    treeData.HalfWidth[0],
                    -treeData.HalfWidth[1],
                    treeData.HalfWidth[1]));

        mTpIdx = treeData.NbTP;

        // allocate claimBank internal memory
        for (cb = claimBank.begin(); cb != claimBank.end(); cb++) {
            Hyb3_CbkInit(&cb->second.zeroBank);
            int nbCrit = cb->second.critZeroMatDates.size();

            if (nbCrit > 0) {
                // claim bank is always 1 dimension as it is always an interest rate
                // and interest rate is always the first dimension in the tree
                if (Hyb3_CbkAlloc(&cb->second.zeroBank, nbCrit, &treeData, 1) != SUCCESS)
                    throw ModelException("Failed to allocate claim bank memory for curve " +
                                         cb->second.curveName);
            }
        }

        if (Hyb3_Dev_Alloc(&devData, &treeData) != SUCCESS)
            throw ModelException("Hyb3_Dev_Alooc failed: " + IrConverter::slog.pop());

        print();
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


void Hyb3CBEQ::retrieveFactor() {
    try{
        int i;

        if (!discYC)
            throw ModelException("FDModel::discYC field is not populated - must contain the top product's "
                                 "discount curve. Internal error in FDModel/product code");

        currency[0] = discYC->getCcy();

        if (IRParams->curveToDiffuse->getCcy() != currency[0])
            throw ModelException("Currency of supplied domestic curve to diffuse " + 
                                 IRParams->curveToDiffuse->getCcy() +
                                 " must be the same as currency defined by the product " +
                                 currency[0]);

        reset();  // initialise hyb3 internal structures

        treeData.Ppy = ppy;
        treeData.NbSigmaMax = maxStdDeviations;
        treeData.TreeType = TTYPE_EQ1IR;

        treeData.CvDisc[0] = 1;
        treeData.CvDiff[0] = 0;
        treeData.CvIdx1[0] = 1;
        treeData.CvIdx2[0] = 2;

        // by convention defined above, set 0=diff curve, 1=discount curve, 2=other
        curveName[0][0] = IRParams->curveToDiffuse->getName();
        IrConverter::to_T_CURVE(treeCurves[0][0], IRParams->curveToDiffuse.get(), false);
        curveName[0][1] = discYC->getName();
        IrConverter::to_T_CURVE(treeCurves[0][1], discYC.get(), false);

        bool thirdCurve;
        thirdCurve = false;
        int eqIndex = -1;

        for (i = 0; i < factors.size(); i++) {
            IMarketFactorConstSP &mf = factors[i];
            const YieldCurve* yc = dynamic_cast<const YieldCurve*>(mf.get());
            const SimpleEquity* eq = dynamic_cast<const SimpleEquity*>(mf.get());
            const IndexSpecEQ* eqSpec = dynamic_cast<const IndexSpecEQ*>(mf.get());

            if (yc) {
                // for now must be either Foreign or Domestic
                if (yc->getCcy() != currency[0])
                    throw ModelException("Yield curve " + yc->getName() +
                                         " defined in instrument is currency " +
                                         yc->getCcy() +
                                         ". HybEq is in single currency (" + currency[0] + 
                                         ") equity mode ");

                if (yc->getName() == curveName[0][0] ||
                    yc->getName() == curveName[0][1]) {
                    continue;
                }
                else {

                    // define third currency or error if already defined
                    if (thirdCurve)
                        throw ModelException("Hyb3 engine supports maximum of three yield curves/currency "
                                             "and supplied yield curve/factor " + yc->getName() +
                                             "is a fourth for currency " + yc->getCcy());
                    curveName[0][2] = yc->getName();
                    IrConverter::to_T_CURVE(treeCurves[0][2], yc, false);
                    thirdCurve = true;
                }
            }
            else if (eq) {
                eqIndex = i;  // store for simple lookup later on
            }
            else if (eqSpec) {

                // ??? should eq index Spec be a factor?
            }
            else {
                throw ModelException("Hyb3CBEQ only supports yield curve and equity market "
                                     "factors - not type " + mf->getClass()->getName());
            }
        }
        // if not supplied, set the third curve to equal the second/discount
        if (thirdCurve == false) {
            curveName[0][2] = curveName[0][1];
            treeCurves[0][2] = treeCurves[0][1];
        }

        // check that the curves have the same value date, and set this to the
        // spot date of the IR curves
        currencyValueDate[0] = DateTime::fromIrDate(treeCurves[0][0].ValueDate);
        for (i = 1; i < 3; i++) {
            if (treeCurves[0][i].ValueDate != currencyValueDate[0].toIrDate())
                throw ModelException("ValueDate " + 
                                     DateTime::fromIrDate(treeCurves[0][i].ValueDate).toString() +
                                     " defined in zero curve index [" + Format::toString(i) + 
                                     "] is not the same "
                                     "as the valueDate defined in the other zero curve(s) " +
                                     currencyValueDate[0].toString());
        }

        // only valid IR and EQ indexSpecs supported - double check
        // if indexSpecFX found, ensure it is domestic/domestic - ie value = 1.0
        for (i = 0; i < indexSpecs.size(); i++) {

            string type = indexSpecs[i]->getClass()->getName();
            if (type != "IndexSpecIR" && 
                type != "IndexSpecEQ") {

                // check indexSpecFX is dom-dom and thus guaranteed a value = 1.0
                IndexSpecFX* spec = dynamic_cast<IndexSpecFX*>(indexSpecs[i].get());
                if (spec) {
                    if (spec->sameCurve == false)
                        throw ModelException("IndexSpecFX " + spec->getName() +
                            "must be a single currency currency FX (ie. rate = 1.0, "
                            "same foreign and domestic curve) for this model");
                }
                else
                    throw ModelException("Index spec type " + type +
                                         " is not supported by Hyb3CBEQ model "
                                         "- must be either EQ or IR");
            }
        }
        populateTreeIRParams();
        strcpy(treeData.Index[0], IRParams->volCalibIndex.c_str());

        // populate equity vols and other equity data into tree
        const SimpleEquity* eq;
        if (eqIndex == -1)
            eq = 0;
        else
            eq = dynamic_cast<SimpleEquity*>(factors[eqIndex].get());

        // correlation - apply override value if supplied
        if (IREQCorrelationOverride)
            eqData.Rho[0][0] = *IREQCorrelationOverride;
        else
            eqData.Rho[0][0] = corrIREQ;

        if ((eqData.Rho[0][0] < -MAXCORRELATION) || 
            (eqData.Rho[0][0] >  MAXCORRELATION))
            throw ModelException("IR-EQ correlation " + 
                                 Format::toString(eqData.Rho[0][0]) + 
                                 " is outside the acceptable range of between +/-" +
                                 Format::toString(MAXCORRELATION));

        PopulateEQData(&eqData,
                       today.toIrDate(),
                       eqSpotVolOverride,
                       eq,
                       EQCalibrationStrike,
                       dividendList);
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

// ??? how do I know for sure if this is an equity value or a rate value - if 
// rate set to dim 1, if equity set to dim 2
void Hyb3CBEQ::sliceExpand(TreeSlice& treeSlice, const string& curveName) const {
    try {
        TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);

        if (slice.getDim()==0)
            slice.expand(1);
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void Hyb3CBEQ::sliceDev(TreeSlice& treeSlice, int curveIdx) const
{
    TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);
    try {
        // nothing to do at the back of the tree 
        if (mTpIdx == treeData.NbTP)
            return;

        // curveIdx == crvIdx * 10 + currIdx, see getCurveIdx
        //int currIdx = curveIdx % 10; // always = 0 in this mode
        int crvIdx = curveIdx / 10;

        int devMode;
        int sliceDim = slice.getDim();

        if (sliceDim == 1)
            devMode = DISC_1D_NOCUPS;
        else if (sliceDim == 2)
            devMode = DISC_2D_NOCUPS;
        else {
            throw ModelException("Slice \""+slice.name+"\" of dimension "+Format::toString(sliceDim)
                +" not supported by Hyb3CBEQ in EQ1IR mode");
        }

        if (Hyb3_Dev(slice.getValuePtr(),
                    mTpIdx, 
                    treeData.NbTP, 
                    crvIdx, 
                    devMode, 
                    const_cast<HYB3_DEV_DATA*>(&devData),
                    const_cast<HYB3_TREE_DATA*>(&treeData)) != SUCCESS) {
            throw ModelException("Hyb3_Dev failed ," + IrConverter::slog.pop());
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__, " slice name = " + slice.name);
    }
}

// ??? needs a rethink when we unify dates - for now implement current logic
DateTime Hyb3CBEQ::getCurveValueDate(string curveName) const {
    try {
        if (curveName.empty()) {
            throw ModelException("curveName string function argument is empty");
        }
        else {  // return value date of one of the interest rate(s)
            int currIdx = getCurrIdx(curveName);
            int index = getCrvIdx(curveName, currIdx);
            // set time to TIME_MIN because SOD if after dividend time XBS
            return DateTime::fromIrDate(treeCurves[currIdx][index].ValueDate);
        }
        return DateTime();
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


void Hyb3CBEQ::populateTreeIRParams()
{
    IRVolRawSP volSP;
    volSP = getIRVolRaw(IRParams);

    IRVolSelector volSelector(volSP, treeCurves[0][0],
                              IRParams->volCalibIndex, cetSmoothing, IRParams);
    volSelector.getMktVolData(mktVolData[0]);

    // ??? temp "intermediate" structure to reuse logic from hyb3 for now
    RateTree::IRModelParams mp;

    // retrieve model parameters from market cache
    IRExoticParamTable* pSmileTable = IRParams->smileTable.get();
    IRExoticParamTable* pModelTable = IRParams->modelTable.get();
    if (pSmileTable || pModelTable)
    {
        // Use new MarketTable object
        IrConverter::to_MODELPARAMETERS_DATA(mp, 
                                            IRParams->smileSet, 
                                            IRParams->modelSet,
                                            engineSet,
                                            *(engineTable.get()),
                                            *pSmileTable,
                                            *pModelTable);
    }
    else
    {
        // Use deprecated IRCalib object 
        IrConverter::to_MODELPARAMETERS_DATA(mp,
                                            IRParams->smileSet, 
                                            IRParams->modelSet,
                                            *(IRParams->irCalib.get()));
    }

    // populate tree structures from the model overwrite parameters or
    // underlying model parameters in the market data */
    populateIRParams(&mktVolData[0], &treeData, &mp);
}


// ??? review fully
void Hyb3CBEQ::print()
{
    int i;

    map<string, NamedClaimBank>::iterator cb;

    if (treeDataDebugFile.empty())
        return;

    FILE *stream = fopen(treeDataDebugFile.c_str(), "w");
    if (stream == NULL) {
        throw ModelException(__FUNCTION__, "Unable to open file " +
                             treeDataDebugFile + " for writing");
    }
    int count=0;

    fprintf(stream, "Model Settings:\n");
    fprintf(stream, "TodayDate: %s\n", today.toString().c_str());
    fprintf(stream, "valueDate: %s\n", valueDate.toString().c_str());
    fprintf(stream, "currency IR spot/value date: %s\n", 
            currencyValueDate[0].toString().c_str());

    fprintf(stream, "\nFollowing critical dates registered with the engine:\n\n");
    for (i = 0; i < (int)sortedCritDates.size(); i++) {
        fprintf(stream, "%d   %ld\n", i, sortedCritDates[i]);
    }

    fprintf(stream, "\nFollowing claim banks registered with the engine:\n\n");

    for (cb = claimBank.begin(); cb != claimBank.end(); cb++) {
        fprintf(stream, "%d, name: %s, curveName: %s, critDatesIdx %d\n",
                count,
                cb->first.c_str(),
                cb->second.curveName.c_str(),
                cb->second.critDatesIdx);

        fprintf(stream, "\n\tCritical Zero Dates used by tree\n");
        fprintf(stream, "\t[idx]  (zero Use/Obs Date)   (zero Maturity Date)\n");

        for (i = 0; i < (int)cb->second.critZeroMatDates.size(); i++) {
            fprintf(stream, "\t%d   %ld     %ld\n",
                    i,
                    cb->second.critZeroUseDates[i],
                    cb->second.critZeroMatDates[i]);
        }

        fprintf(stream, "\n\tOptional Zero Dates\n");
        fprintf(stream, "\t[idx]  (zero Use/Obs Date)   (zero Maturity Date)\n");

        for (i = 0; i < (int)cb->second.optZeroMatDates.size(); i++) {
            fprintf(stream, "\t%d   %ld     %ld\n",
                    i,
                    cb->second.optZeroUseDates[i],
                    cb->second.optZeroMatDates[i]);
        }
        fprintf(stream, "\n\tUnsorted Critical Zero Dates (as registered by products "
                        "before internal sorting and duplicate removal)\n");
        fprintf(stream, "\t[idx]  (zero Use/Obs Date)   (zero Maturity Date)    (label)\n");

        for (i = 0; i < (int)cb->second.critZeroMatDatesUnsorted.size(); i++) {
            fprintf(stream, "\t%d   %ld     %ld    %s\n",
                    i,
                    cb->second.critZeroUseDatesUnsorted[i],
                    cb->second.critZeroMatDatesUnsorted[i],
                    cb->second.critZeroLabel[i].c_str());
        }
    }
    
    fprintf(stream, "\nDomestic Currency %s\n\n", currency[0].c_str());
    fprintf(stream, "Nb  DomesticCurves\n");
    for (i = 0; i < 3; i++) {

        fprintf(stream, "%d  %s\n", 
                i,
                curveName[0][i].c_str());
    }

    fprintf (stream, "\n\n");

    // Slice sizes used for allocation
    fprintf(stream, "  W1   HW1    W2   HW2    W3   HW3\n");
    fprintf(stream, "%5d %5d %5d %5d %5d %5d\n\n\n",
                    treeData.Width[0], treeData.HalfWidth[0],
                    treeData.Width[1], treeData.HalfWidth[1],
                    treeData.Width[2], treeData.HalfWidth[2]);

    // Schedule of the underlying swap 
    fprintf (stream,"\nTIMELINE INFORMATION (IR)\n");

    fprintf (stream,
            "Node      Date     Days  Max  Forward   "
            "Zero0   Discount0   Zero1   Discount1   "
            "Zero2   Discount2    IrMidNode    SpotVol \n");

    for (i = 0; i <= treeData.NbTP; i++) {

        double Forward = 100. * treeData.FwdRate[0][0][i]/treeData.Length[i];
        int daysL = Daysact (treeData.TPDate[0], treeData.TPDate[i]);

        fprintf(stream,
                "[%3d] \t%8ld  %3d  %3d  %7.4f  %7.4f   %8.6f  "
                "%7.4f   %8.6f  %7.4f   %8.6f    %9.6f    %6.2f \n",
                i,
                treeData.TPDate[i],
                daysL,
                treeData.Top1[i],
                Forward,
                treeData.ZeroRate[0][0][i] * 100.,
                treeData.ZeroCoupon[0][0][i],
                treeData.ZeroRate[0][1][i] * 100.,
                treeData.ZeroCoupon[0][1][i],
                treeData.ZeroRate[0][2][i] * 100.,
                treeData.ZeroCoupon[0][2][i],
                exp(treeData.IrZCenter[0][i]),
                treeData.SpotVol[0][i] * 100.);
        fflush(stream);
    }

    fprintf (stream, "\n\n");
    fprintf (stream,"\nTIMELINE INFORMATION (Stock)\n");
    fprintf (stream, "Node \t  Date \t\t   Forward      "
                     "Center      Eq vol  Spotvol   Rho's (Irf-Eq, Ird-Eq, FX-Eq)   "
                     "Dividend   DivAmt   SettlDate\n");

    for (i = 0; i <= treeData.NbTP; i++) {

        fprintf (stream, "[%3d] \t%8ld \t%10.4f \t%10.4f  \t"
                         "%5.2f \t%5.2f \t     %4.2f       %4.2f      %4.2f     "
                         "(%d)     %7.3f    %8ld\n",
                 i,
                 treeData.TPDate[i],
                 treeData.FwdEq[i],
                 treeData.EqMidNode[i],
                 treeData.EqVol[i] * 100.,
                 treeData.SpotEqVol[i] * 100.,
                 treeData.Rho[3][i],
                 treeData.Rho[4][i],
                 treeData.Rho[5][i],
                 treeData.TPType[12][i],
                 treeData.CritDate[12][i].Value[0],
                 treeData.NodeSettleDate[i]);
    }

    // Print out input zero curves
    int k;
    fprintf (stream, "\n\n");
    fprintf (stream, "%s Currency Curves (Index, COF, Risk)\n\n", 
             currency[0].c_str());
    for (k = 0; k < 3; k++) {

        fprintf (stream, "Maturity      Zero       Discount \n");

        for (i = 0; i < (treeCurves[0][k]).NbZero; i++) {

            double days = Daysact(valueDate.toIrDate(), 
                                (treeCurves[0][k]).ZeroDate[i]);
            // Discount factor up to the current date
            double discount = ::pow(1. + (treeCurves[0][k]).Zero[i], -days/365.);

            fprintf (stream,
                    "%ld   %9.6f     %8.6f \n",
                    (treeCurves[0][k]).ZeroDate[i],
                    (treeCurves[0][k]).Zero[i] * 100.,
                    discount);
        }
        fprintf (stream, "\n");
    }

    fclose(stream);
}

// ??? figure out how to turn this on/off but preferably not defining a full
// outputRequest for each one - maybe have a "debug info" output request
// that switches on all of this debug packet info
void Hyb3CBEQ::recordOutput(Control* ctrl, Results* results) const {
    try {
        int i, j;

        if (!ctrl->requestsOutput(OutputRequest::DBG))
            return; // no debug info requested

        // record today and valueDate
        DateTimeSP todaySP(new DateTime(today));
        results->storeGreek(todaySP,
                            Results::DEBUG_PACKET,
                            OutputNameConstSP(new OutputName("TODAY")));

        DateTimeSP valueDateSP(new DateTime(valueDate));
        results->storeGreek(valueDateSP,
                            Results::DEBUG_PACKET,
                            OutputNameConstSP(new OutputName("VALUE_DATE")));

        DateTimeSP curveValueDateSP(new DateTime(currencyValueDate[0]));
        results->storeGreek(curveValueDateSP,
                            Results::DEBUG_PACKET,
                            OutputNameConstSP(new OutputName("IR_SPOT_DATE")));

        // ??? record equity spot date

        // store tree critical dates
        DateTimeArraySP treeDates(new DateTimeArray(treeData.NbTP+1));
        for (i = 0; i <= treeData.NbTP; i++)
            (*treeDates)[i] = DateTime::fromIrDate(treeData.TPDate[i]);

        results->storeGreek(treeDates,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("TREE_DATES")));

        // store number of days offsets from start of tree for each critical date
        IntArraySP daysOffsets(new IntArray(treeData.NbTP+1));
        for (i = 0; i <= treeData.NbTP; i++)
        {
            double days = Daysact(treeData.TPDate[0], treeData.TPDate[i]);
            (*daysOffsets)[i] = (int)days;
        }

        results->storeGreek(daysOffsets,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("TREE_DATE_OFFSETS")));

        // store DOMESTIC ir rates
        CDoubleMatrixSP domZeros(new CDoubleMatrix(3, treeData.NbTP+1));
        for (i = 0; i < 3; i++)
            for (j = 0; j <= treeData.NbTP; j++)
                (*domZeros)[i][j] = treeData.ZeroRate[0][i][j];

        results->storeGreek(domZeros,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("DOMESTIC_IR_ZEROS")));

        // store domestic IR spotVols
        DoubleArraySP domIRSpotVol(new DoubleArray(treeData.NbTP+1));
        for (i = 0; i <= treeData.NbTP; i++)
            (*domIRSpotVol)[i] = treeData.SpotVol[0][i];

        results->storeGreek(domIRSpotVol,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("DOMESTIC_IR_SPOT_VOLS")));

        // store EQ spot vols
        DoubleArraySP eqSpotVols(new DoubleArray(treeData.NbTP+1));
        for (i = 0; i <= treeData.NbTP; i++)
            (*eqSpotVols)[i] = treeData.SpotEqVol[i];

        results->storeGreek(eqSpotVols,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("EQ_SPOT_VOLS")));

        // store EQ composite vols
        DoubleArraySP eqCompVols(new DoubleArray(treeData.NbTP+1));
        for (i = 0; i <= treeData.NbTP; i++)
            (*eqCompVols)[i] = treeData.EqVol[i];

        results->storeGreek(eqCompVols,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("EQ_COMP_VOLS")));

        // store forward EQ
        DoubleArraySP eqFwds(new DoubleArray(treeData.NbTP+1));
        for (i = 0; i <= treeData.NbTP; i++)
            (*eqFwds)[i] = treeData.FwdEq[i];

        results->storeGreek(eqFwds,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("EQ_FWDS")));

        /*
        // ??? implement properly for equity smile 
        // Store local volatility function for various strikes
        CDoubleMatrixSP smileVols(new CDoubleMatrix(5, treeData.NbTP+1));
        DoubleArraySP smileDeltas(new DoubleArray(5));
        (*smileDeltas)[0] = 0.1;
        (*smileDeltas)[1] = 0.25;
        (*smileDeltas)[2] = 0.5;
        (*smileDeltas)[3] = 0.75;
        (*smileDeltas)[4] = 0.9;

        for (i = 0; i <= treeData.NbTP;i++)
        {
            double time = (*daysOffsets)[i] / 365.0;
            double vol = sqrt(time)* treeData.FxVol[i]; // use comp Vol for estiamting delta
            double fxSpotVol = treeData.SpotFxVol[i];   // spot Vol for local vol function

            for (j = 0; j < 5; j++)
            {
                double delta = (*smileDeltas)[j];
                double strikeFwdFx = deltaToStrike(delta, vol, treeData.FwdFx[i]);
                double moneyness = strikeFwdFx/treeData.FwdFx[i];
                double smile;

                if (Hyb3_Gfunc(&smile,
                            moneyness,
                            treeData.A1C[i],
                            treeData.A2C[i],
                            treeData.A3C[i]) != SUCCESS)
                    throw ModelException("Hyb3_Gfunc failed: "+ IrConverter::slog.pop());

                smile *= fxSpotVol/moneyness;

                (*smileVols)[j][i] = smile;
            }
        }

        results->storeGreek(smileVols,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("FX_SMILE_VOLS")));

        results->storeGreek(smileDeltas,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("FX_SMILE_DELTAS")));
        */
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


void Hyb3CBEQ::getEQIndex(TreeSlice& treeSlice)
{
    const char* routine = "Hyb3CBEQ::getEQIndex";
    TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);

    if (slice.getDim() != 2)
        slice.allocDim(2);

    // for now simply copy the EQ Spot slice into the provided slice
    if (Hyb3_CopySlice(slice.getValuePtr(),
                       devData.EqSpot,
                       2,       // EQ always 2D slice in this mode
                       mTpIdx,  // current step
                       &treeData) != SUCCESS) {
        throw ModelException(routine, IrConverter::slog.pop());
    }
    slice.treeStep = range->treeStep;
}

// ??? todo - date adjustment, and checking if index spec has an overridden value
// that should be used instead of returning the tree estimated equity value
void Hyb3CBEQ::getEQIndex(const IndexSpecEQ& eqSpec, DateTime currentDate, DateTime resetDate,
                          TreeSlice& treeSlice) {
    try {
        if (currentDate != resetDate)
            throw ModelException("Currently resetDate " + resetDate.toString() +
                                 " must equal currentDate " + currentDate.toString() +
                                 " as date offset logic is not yet implemented ");

        TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);

        if (slice.getDim() != 2)
            slice.allocDim(2);

        // for now simply copy the EQ Spot slice into the provided slice
        if (Hyb3_CopySlice(slice.getValuePtr(),
                           devData.EqSpot,
                           2,       // EQ always 2D slice in this mode
                           mTpIdx,  // current step
                           &treeData) != SUCCESS) {
            throw ModelException("Hyb3_CopySlice failed: " + IrConverter::slog.pop());
        }
        slice.treeStep = range->treeStep;
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

// hyb3 eq mode always uses index [0] for domestic yield curve
int Hyb3CBEQ::getCurrIdx(const string& /*curveName*/) const {
    
    return 0;
}

// interest rate curves are always one/first dimension
int Hyb3CBEQ::getCrvDim(const string& curveName) const {

    return 1;
}

// forward the price back to the valueDate as the convention is to report values to
// the valueDate,
double Hyb3CBEQ::getPrice0(const TreeSlice& price) const {
    try {
        double vdPrice;
        const TreeSliceLayer * layer = dynamic_cast< const TreeSliceLayer * >( &price );
        if (layer)
            vdPrice = layer->getPrice0();
        else
            vdPrice = price.getCentre();

        // ??? old method:
        // now forward to valueDate using precalculate discount factors between today/valueDate
        // tree is setup so curveToDiffuse = [1], which is what we are reporting the price in
        // so we use that curve's discountFactor to forward back to valueDate
        // vdPrice /= today2ValDateZeros[0][1];

        // ??? new method:
        // as today=valueDate in the current implementation, no need to forward value
        // from today to valueDate
        return vdPrice;
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

DRLIB_END_NAMESPACE
