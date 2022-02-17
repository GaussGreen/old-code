//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : Hyb3Eq.hpp
//
//   Description : hyb3 domestic equity mode FD model class
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Hyb3EQ.hpp"
#include "edginc/TreeSlice.hpp"
#include "edginc/IRVolPair.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/IrConverter.hpp"
#include "edginc/IndexSpecEQ.hpp"
#include "edginc/IndexSpecFX.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SimpleEquity.hpp"
#include "edginc/VolSurface.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/SRMEQVol.hpp"
#include "edginc/FixedSettlement.hpp"
#include "edginc/RollingSettlement.hpp"
#include "esl_log.h"
#include "esl_market.h"


DRLIB_BEGIN_NAMESPACE

/***********************************
 *
 *
 *
 ***********************************/

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
    long       valueDate,
    double     *eqSpotVolOverride,
    const SimpleEquity *eq,
    double     EQCalibrationStrike,
    DividendListConstSP &dividendList,
    DateTime today)
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
        eq_data->ValueDate = today.toIrDate();

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
            eq_data->VolDate[0] = valueDate;
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
        if (eq_data->SettleType == '0') throw ModelException("Unsupported settlement type "+settlement.getClass()->getName());
    }
    else
    {
        // if equity not used in product, set default values
        long dummyDate = Nxtday(valueDate, 100);  // dummy first vol date

        eq_data->Spot = 0.0;
        eq_data->NbVol = 1;
        eq_data->ValueDate = today.toIrDate();
        eq_data->VolDate[0] = valueDate;
        eq_data->VolDate[1] = dummyDate;
        eq_data->Vol[0] = 0.1;
        eq_data->Vol[1] = 0.1;
        eq_data->NbInpSpotVol = 0;
        eq_data->InpSpotVolDate[0] = valueDate;
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
    if (Hyb3_Eq_Check_Dyn_WType4(eq_data,
                                 valueDate) != SUCCESS)
        throw ModelException(routine, IrConverter::slog.pop());

    if (Hyb3_Eq_Check_Sta_W(eq_data) != SUCCESS)
        throw ModelException(routine, IrConverter::slog.pop());
}


 /***********************************
 *
 * Hyb3 Domestic Equity tree mode
 *
************************************/


void Hyb3EQ::load(CClassSP& clazz){

    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(Hyb3EQ, clazz);
    SUPERCLASS(Hyb3);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(ppy, "minimum number of tree nodes(points)/year");
    FIELD(maxStdDeviations, "number of standard deviations to trim the tree at");
    FIELD(IRParams, "Collection of interest rate model parameters");
    FIELD(EQSmileParams, "Collection of named Equity smile parameters - assumes no Equity smile if not supplied");
    FIELD_MAKE_OPTIONAL(EQSmileParams);
    FIELD(IREQCorrelation, "User defined correlation override")
    FIELD_MAKE_OPTIONAL(IREQCorrelation);
    FIELD(EQCalibrationStrike, "Strike used to determine volatility matrix column for tree EQ vol calibration");
    FIELD_MAKE_OPTIONAL(EQCalibrationStrike);
    FIELD(treeDataDebugFile, "Print tree debug data to named output file");
    FIELD_MAKE_OPTIONAL(treeDataDebugFile);
    FIELD(zeroInterpStyle, "zero curve interpolation style - LINEAR or FLAT_FWD.  Defaults to "
                                  "LINEAR if not supplied");
    FIELD_MAKE_OPTIONAL(zeroInterpStyle);
    FIELD(lastDividendDate,"");
    FIELD_MAKE_OPTIONAL(lastDividendDate);
    FIELD(eqSpotVolOverride,"Overrides market comp/spot vols for equity - sets constant spot vol for EQ "
                            "diffusion (units = decimal NOT percentage)");
    FIELD_MAKE_OPTIONAL(eqSpotVolOverride);

    FIELD(corrIREQ,"");
    FIELD_MAKE_TRANSIENT(corrIREQ);
    FIELD(mEQName,"");
    FIELD_MAKE_TRANSIENT(mEQName);
}

CClassConstSP const Hyb3EQ::TYPE = CClass::registerClassLoadMethod(
    "Hyb3EQ", typeid(Hyb3EQ), Hyb3EQ::load);

// for class loading
bool Hyb3EQLoad(void) { return (Hyb3EQ::TYPE != 0); }


void Hyb3EQ::getMarket(const MarketData*  market, IInstrumentCollectionSP instruments)
{
    static const string routine = "Hyb3EQ::getMarket";
    try 
    {
    mToday = market->GetReferenceDate();
    IRParams->getData(this, market);

    // Get Engine
    if (!engineTable.isEmpty())
        engineTable.getData(this, market);

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

            string name = factor->getName();
            string type = factor->getClass()->getName();
            int i;
            for( i = 0; i < factors.size(); ++i )
            {
                if( name == factors[ i ]->getName() && type == factors[ i ]->getClass()->getName() )
                     break;
            }
            if( i >= factors.size() && model->acceptFactor( factor ) )
            {
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
    for (i = 0; i < factors.size(); i++)
    {
        IMarketFactorConstSP &mf = factors[i];
        const SimpleEquity* eq = dynamic_cast<const SimpleEquity*>(mf.get());
        const IndexSpecEQ* eqSpec = dynamic_cast<const IndexSpecEQ*>(mf.get());

        if (eq)
        {
            if (eqIndex == -1)
            {
                eqIndex = i;
                mEQName = eq->getName();
            }
            else
            {
                if (eq->getName() != mEQName)
                    throw ModelException("Hyb3EQ only supports ONE equity - at least two equities "
                                         "supplied (" + eq->getName() + ", " +
                                         factors[eqIndex]->getName() + ")");
            }
        }
        // ??? serious issues here - indexSpecEQ should be registering market factors but doesn't seem to,
        // so I will do it here for now - indexSpecEQ getting complicated and don't want to mess change anything
        // now in case it breaks something else
        if (eqSpec)
        {
            if (eqIndex == -1)
            {
                const IMarketFactorConstSP &mfEQ = eqSpec->getFactor();
                IMarketFactorSP nonConstmfEQ(const_cast<IMarketFactor*>(mfEQ.get()));
                string name = mfEQ->getName();
                string type = mfEQ->getClass()->getName();

                const SimpleEquity* eqEQ = dynamic_cast<const SimpleEquity*>(mfEQ.get());
                if (!eqEQ)
                    throw ModelException("Hyb3EQ currently only supports SimpleEquity equity type");
                mEQName = eqEQ->getName();
                // for now, add this simple equity as a marketFactor so it can be more easily retrieved
                // in the retrieve factors function
                factors.push_back(nonConstmfEQ);
            }
        }
    }
//    if (mEQName.empty())
//        throw ModelException("One equity must be supplied in instrument for Hyb3EQ model");

    // recurse the components to store specific index specs 
    class RetrieveIndexSpecs : public ObjectIteration::IAction
    {
        IndexSpecArray & indexSpecs;

    public:
        RetrieveIndexSpecs(IndexSpecArray & indexSpecs) :
              indexSpecs( indexSpecs ) {}

        virtual bool invoke(ObjectIteration::State& state, IObjectConstSP obj)
        {
            int i;
            IndexSpec * indexSpec = dynamic_cast< IndexSpec * >( state.getObject().get() );

            string name = indexSpec->getName();
            string type = indexSpec->getClass()->getName();

            for( i = 0; i < indexSpecs.size(); ++i )
            {
                if( name == indexSpecs[ i ]->getName() && type == indexSpecs[ i ]->getClass()->getName() )
                    break;
            }
            if( i >= indexSpecs.size() )
            {
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
    if (market->hasCorrelationData(IRParams->curveToDiffuse->getName(), mEQName))
    {
        string corrName = market->getCorrelationName(IRParams->curveToDiffuse->getName(), mEQName);

        CorrelationBaseSP corrBase = getCorrelation(corrName,
                                                    IRParams->curveToDiffuse->getClass(),
                                                    SimpleEquity::TYPE,
                                                    Correlation::TYPE,
                                                    market);
        corrFound = true;
        corrIREQ = CorrelationSP::dynamicCast(corrBase)->getCorrelation();
    }

    // check if eq/(ir discount) correlation supplied
    if (market->hasCorrelationData(discYC->getName(), mEQName))
    {
        string corrName = market->getCorrelationName(discYC->getName(), mEQName);

        CorrelationBaseSP corrBase = getCorrelation(corrName,
                                                    discYC->getClass(),
                                                    SimpleEquity::TYPE,
                                                    Correlation::TYPE,
                                                    market);
        double result = CorrelationSP::dynamicCast(corrBase)->getCorrelation();
        if (corrFound)  //  check same value
        {
            if (!Maths::equals(corrIREQ, result))
                throw ModelException("Correlation found for discount yieldCurve " +
                                     discYC->getName() + " and equity " + mEQName + "(value = " + 
                                     Format::toString(result) +
                                     ") does not equal previously found correlation for "
                                     "curveToDiffuse " + IRParams->curveToDiffuse->getName() +
                                     " and equity " + Format::toString(corrIREQ));
        }
        else
        {
            corrFound = true;
            corrIREQ = CorrelationSP::dynamicCast(corrBase)->getCorrelation();
        }
    }

    // ??? maybe check 3rd curve/Eq correlation if supplied but don't know if that curve
    // has been supplied until retrieve factor - think about refactoring the whole thing later
    if (mEQName.empty())
        corrIREQ = 0.0;
    else
        if (corrFound == false)
            throw ModelException("Correlation not found in market for IR and equity "
                                 "searched for the following 2 correlations in the market: " + 
                                 market->getCorrelationName(IRParams->curveToDiffuse->getName(), mEQName) +
                                 " and " +
                                 market->getCorrelationName(discYC->getName(), mEQName));
    }
    catch (exception& e) {
        throw ModelException(e, routine);
    }
}

void Hyb3EQ::initModel()
{
try 
{
    // ??? note:  ZeroInterpTypeFlag is global variable defined in both hyb3.lib
    // and fix3.lib - tidy up but have to check all hyb3 products

    mTreeBuilt = true; 
    ESL_INTERP_TYPE interpType = ESL_INTERP_LINEAR;

    if (zeroInterpStyle.empty())
        interpType = ESL_INTERP_LINEAR;
    else if (zeroInterpStyle == "LINEAR")  // ??? make case insensitive
        interpType = ESL_INTERP_LINEAR;
    else if (zeroInterpStyle == "FLAT_FWD")
        interpType = ESL_INTERP_FLATFWD;
    else
        throw ModelException("Supplied zeroInterpStyle flag (" +
                             zeroInterpStyle +
                             ") is not valid - must be either LINEAR or FLAT_FWD");

    EslSetZeroInterpolation(interpType);
    //EslSetZeroInterpolationStub(ESL_INTERP_LINEAR);

    mToday = DateTime::fromIrDate(mTCurves[0][0].ValueDate);

    CRIT_DATE* critDate=NULL;
    int nbCritDate = 0;
    int i,j;

    map<string, ZeroBond>::iterator zero;
    map<string, NamedZeroBank>::iterator zb;

    // needed to setup timeline
    critDate = (CRIT_DATE *) DR_Array (CRITDATE, 0, 0);
    if (critDate == NULL)
        throw ModelException(IrConverter::slog.pop());

    // value date always critical date
    if (Add_To_DateList(&nbCritDate, 
                        &critDate,
                        getValueDate().toIrDate(),
                        0, 0, 0, 0, 0, 0, 0, 0, 0) != SUCCESS)
        throw ModelException(IrConverter::slog.pop());

    arrangeDates();

    for (j = 0; j < critDates.size(); j++)
    {

        DateTime critDateL = critDates[j];
        if (critDateL > getValueDate())
        {
            if (Add_To_DateList(&nbCritDate, 
                                &critDate,
                                critDateL.toIrDate(),
                                0, 0, 0, 0, 0, 0, 0, 0, 0) != SUCCESS)
                throw ModelException(IrConverter::slog.pop());
        }
    }

    // sort critical dates, so we can stop adding dividends after last exDiv date
    if (Sort_CritDate(nbCritDate, critDate) != SUCCESS )
    {
        throw ModelException(IrConverter::slog.pop());
    }

    // ??? review
    if (dividendList.get())
    {
        const DividendArray& divs = dividendList->getArray();
        IRDate lastTreeDate = critDate[MAX(nbCritDate-1,0)].CritDate;

        for (j = 0; j < divs.size(); j++)
        {
            DateTime exDate = divs[j].getExDate();
            if (exDate.getDate() >= getValueDate().getDate() /* exDate occurs @time XBS, not SOD */
                && (exDate.toIrDate() <= lastTreeDate ) /* make sure we do not build beyond infinity */
                && (lastDividendDate.empty() || exDate <= lastDividendDate))
            {
                if (Add_To_DateList (&nbCritDate,
                                     &critDate,
                                     exDate.toIrDate(),
                                     12,
                                     divs[j].getDivAmount(),
                                     0, 0, 0, 0,
                                     0, 0, 0) != SUCCESS)
                {
                    throw ModelException(IrConverter::slog.pop());
                }
            }
        }
    }

    // sort critical dates for debug output
    if (Sort_CritDate(nbCritDate, critDate) != SUCCESS)
    {
        throw ModelException(IrConverter::slog.pop());
    }

    sortedCritDates.resize(1);
    sortedCritDates[0] = critDate[0].CritDate;
    for (i = 1; i < nbCritDate; i++)
    {
        if (sortedCritDates.back() == critDate[i].CritDate)
            continue;
        sortedCritDates.push_back(critDate[i].CritDate);
    }

     // if in CLAIMBANK mode, insert critical dates associated with the claimbank
    if (zeroBankMode == ZeroBankMode::CLAIMBANK)
    {
        throw ModelException("Claim bank not yet suppported");
    }
/*
    if (!treeDataDebugFile.empty())
    {
        FILE* stream = fopen("c:/hyb3EQCritDates.dat", "w");
        // sort dates before printing them
        if (Sort_CritDate(nbCritDate,
                          critDate) != SUCCESS)
        {
            throw ModelException(IrConverter::slog.pop());
        }

        PrintCriticalDates(stream, nbCritDate, critDate);
        fclose(stream);
    }
*/
    // make the timeline 
    if (Hyb3_Time_Line (mFxData.ValueDate, nbCritDate, critDate, 'I', &mTreeData) != SUCCESS)
        throw ModelException(IrConverter::slog.pop());


    if (Hyb3_Build_Tree(1, mTCurves, mMktVolData, &mFxData, &mEqData, 
                        &mTreeData) != SUCCESS)
    {
        throw ModelException(IrConverter::slog.pop());
    }

    range.reset(new TreeSliceRates::Range(
                -mTreeData.HalfWidth[0],
                 mTreeData.HalfWidth[0],
                -mTreeData.HalfWidth[1],
                 mTreeData.HalfWidth[1]));

    mTpIdx = mTreeData.NbTP;

    if (zeroBankMode == ZeroBankMode::ZEROBANK)
    {   
        for (zb = mZeroBank.begin(); zb != mZeroBank.end(); zb++)
        {
            string crvName = zb->second.curveName;  // ??? to help debugging
            int crvIdx = getCrvIdx(crvName, 0);

            // zero bank always domestic - 1D
            zb->second.zeroBankDim = 1;
            zb->second.currencyEnum = DOMESTIC;
            zb->second.discCurve = crvIdx;
            zb->second.devMode = DISC_1D_NOCUPS;

            int nbZeros = zb->second.nbZeros;
            zb->second.zeroCritDates.resize(nbZeros);
            zb->second.zeroBank.resize(nbZeros);

            for (i = 0; i < nbZeros; i++)
            {
                zb->second.zeroBank[i] = Hyb3_Alloc_Slice(&mTreeData, zb->second.zeroBankDim);
                if (zb->second.zeroBank[i] == NULL)
                {
                    throw ModelException("Unable to allocate memory for zero banks");
                }
            }
        }

        // loop through each stored named zero bond/discounter
        for (zero = mZeroBond.begin(); zero != mZeroBond.end(); zero++)
        {
            zero->second.currencyEnum = DOMESTIC;   // ?? remove
            zero->second.sliceDim = 1;
            zero->second.currency = mCurrency[0];
            zero->second.devMode = DISC_1D_NOCUPS;

            zero->second.zeroSlice = Hyb3_Alloc_Slice(&mTreeData, zero->second.sliceDim);
            if (zero->second.zeroSlice == NULL)
            {
                throw ModelException("Unable to allocate memory for named zero bond "
                    +zero->second.getName());
            }
        }
    }
    
    if (Hyb3_Dev_Alloc(&mDevData, &mTreeData) != SUCCESS)
    {
        throw ModelException(IrConverter::slog.pop()
        +" Failed to allocate storage for DEV data.");
    }

    print();
}
catch (exception& e) {
    throw ModelException(e, __FUNCTION__);  // catch anything else
}
}



void Hyb3EQ::retrieveFactor(){
    try{
        int i;

        if (!discYC)
            throw ModelException("FDModel::discYC field is not populated - must contain the top product's "
                                 "discount curve. Internal error in FDModel/product code");

        mCurrency[0] = discYC->getCcy();

        if (IRParams->curveToDiffuse->getCcy() != mCurrency[0])
            throw ModelException("Currency of supplied domestic curve to diffuse " + 
                                 IRParams->curveToDiffuse->getCcy() +
                                 " must be the same as currency defined by the product " +
                                 mCurrency[0]);

        Hyb3::init();
        mTreeData.Ppy = ppy;
        mTreeData.NbSigmaMax = maxStdDeviations;
        mTreeData.TreeType = TTYPE_EQ1IR;

        mTreeData.CvDisc[0] = 1;
        mTreeData.CvDiff[0] = 0;
        mTreeData.CvIdx1[0] = 1;
        mTreeData.CvIdx2[0] = 2;

        // by convention defined above, set 0=diff curve, 1=discount curve, 2=other
        mCurveName[0][0] = IRParams->curveToDiffuse->getName();
        IrConverter::to_T_CURVE(mTCurves[0][0], IRParams->curveToDiffuse.get(), false);
        mCurveName[0][1] = discYC->getName();
        IrConverter::to_T_CURVE(mTCurves[0][1], discYC.get(), false);

        bool thirdCurve;
        thirdCurve = false;
        int eqIndex = -1;
        for (i = 0; i < factors.size(); i++)
        {
            IMarketFactorConstSP &mf = factors[i];
            const YieldCurve* yc = dynamic_cast<const YieldCurve*>(mf.get());
            const SimpleEquity* eq = dynamic_cast<const SimpleEquity*>(mf.get());
            const IndexSpecEQ* eqSpec = dynamic_cast<const IndexSpecEQ*>(mf.get());

            if (yc)
            {
                // for now must be either Foreign or Domestic
                if (yc->getCcy() != mCurrency[0])
                    throw ModelException("Yield curve " + yc->getName() +
                                         " defined in instrument is currency " +
                                         yc->getCcy() +
                                         ". HybEq is in single currency (" + mCurrency[0] + 
                                         ") equity mode ");

                if (yc->getName() == mCurveName[0][0] ||
                    yc->getName() == mCurveName[0][1])
                    continue;
                else
                {
                    // define third currency or error if already defined
                    if (thirdCurve)
                        throw ModelException("Hyb3 engine supports maximum of three yield curves/currency "
                                             "and supplied yield curve/factor " + yc->getName() +
                                             "is a fourth for currency " + yc->getCcy());
                    mCurveName[0][2] = yc->getName();
                    IrConverter::to_T_CURVE(mTCurves[0][2], yc, false);
                    thirdCurve = true;
                }
            }
            else if (eq)
            {
                eqIndex = i;  // store for simple lookup later on
            }
            else if (eqSpec)
            {
                // ??? should eq index Spec be a factor?
            }
            else
            {
                throw ModelException("Hyb3EQ only supports yield curve and equity market factors - not type " +
                                     mf->getClass()->getName());
            }
        }
        // if not supplied, set the third curve to equal the second/discount
        if (thirdCurve == false)
        {
            mCurveName[0][2] = mCurveName[0][1];
            mTCurves[0][2] = mTCurves[0][1];
        }

        // only IR and EQ indexSpecs supported - double check
        for (i = 0; i < indexSpecs.size(); i++)
        {
            string type = indexSpecs[i]->getClass()->getName();
            if (type != "IndexSpecIR" && 
                type != "IndexSpecEQ")
            {
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
        strcpy(mTreeData.Index[0], IRParams->volCalibIndex.c_str());

        // populate equity vols and other equity data into tree
        const SimpleEquity* eq;
        if (eqIndex == -1)
            eq = 0;
        else
            eq = dynamic_cast<SimpleEquity*>(factors[eqIndex].get());

        if (!IREQCorrelation.empty())
            throw ModelException("Model correlation family overrides not yet supported, but will be "
                                 "a requirement to be able to easily override correlations");
        else
        {
            mEqData.Rho[0][0] = corrIREQ;
            if ((mEqData.Rho[0][0] < -MAXCORRELATION) || 
                (mEqData.Rho[0][0] >  MAXCORRELATION))
                throw ModelException("IR-EQ correlation " + 
                                     Format::toString(mEqData.Rho[0][0]) + 
                                     " is outside the acceptable range of between +/-" +
                                     Format::toString(MAXCORRELATION));
        }

        
        if (!EQSmileParams.empty())
            throw ModelException("Model defined EQSmile Parameter Family name is not yet "
                                 "supported - model reads FXSmile params directly from "
                                 "the SRMEQVol market structure");

        PopulateEQData(&mEqData,
                       mTCurves[0][0].ValueDate,  // ??? equity mode uses this to define value date
                       eqSpotVolOverride,   // ??? hardwire spot cutoff parameter for now
                       eq,
                       EQCalibrationStrike,
                       dividendList,
                       mToday);
    }
    catch (exception& e){
        throw ModelException(e, "Hyb3EQ::retrieveFactor");
    }
}

void Hyb3EQ::sliceExpand(TreeSlice& treeSlice, const string& curveName) const {
    try {
        TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);

        if (slice.getDim()==0)
            slice.expand(2);
    }
    catch (exception& e) {
        throw ModelException(e, "Hyb3EQ::sliceExpand");
    }
}

void Hyb3EQ::sliceDev(TreeSlice& treeSlice, int curveIdx) const
{
    TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);
    try {

        // nothing to do at the back of the tree 
        if (mTpIdx == mTreeData.NbTP)
            return;

        // curveIdx == crvIdx * 10 + currIdx, see getCurveIdx
        int crvIdx = curveIdx / 10;

        int devMode;
        int sliceDim = slice.getDim();

        if (sliceDim == 1)
            devMode = DISC_1D_NOCUPS;
        else if (sliceDim == 2)
            devMode = DISC_2D_NOCUPS;
        else
        {
            throw ModelException("Slice \""+slice.name+"\" of dimension "+Format::toString(sliceDim)
                +" not supported by Hyb3EQ in EQ1IR mode");
        }

        if (Hyb3_Dev(slice.getValuePtr(),
                    mTpIdx, 
                    mTreeData.NbTP, 
                    crvIdx, 
                    devMode, 
                    (/*const_cast*/HYB3_DEV_DATA*)&mDevData, 
                    (/*const_cast*/HYB3_TREE_DATA*)&mTreeData) != SUCCESS)
        {
            throw ModelException(IrConverter::slog.pop());
        }
    }
    catch (exception& e) {
        throw ModelException(e, "Hyb3EQ::sliceDev, slice.name=\""+slice.name+"\"");
    }
}

// ??? needs a rethink when we unify dates - for now implement current logic
DateTime Hyb3EQ::getCurveValueDate(string curveName) const
{
    if (curveName.empty())
    {
        if (mTreeData.TreeType == TTYPE_EQ1IR)
            return DateTime::fromIrDate(mTCurves[0][0].ValueDate);
    }
    else  // return value date of one of the interest rate(s)
    {
        int currIdx = getCurrIdx(curveName);
        int index = getCrvIdx(curveName, currIdx);
        // set time to TIME_MIN because SOD if after dividend time XBS
        return DateTime::fromIrDate(mTCurves[currIdx][index].ValueDate);
    }
    return DateTime();
}


void Hyb3EQ::populateTreeIRParams()
{
    IRVolRawSP volSP;
    volSP = getIRVolRaw(IRParams);

    IRVolSelector volSelector(volSP, mTCurves[0][0],
                              IRParams->volCalibIndex,
                              cetSmoothing, IRParams);
    volSelector.getMktVolData(mMktVolData[0]);

    // ??? temp "intermediate" structure to reuse logic from hyb3 for now
    RateTree::IRModelParams mp;

    // retrieve smile parameters object from smile market table
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
    populateIRParams(&mMktVolData[0], &mTreeData, &mp);
}

// ??? tidy up later, but for now essentially replicate term.prn from calleqrib_t
void Hyb3EQ::print()
{
    static char const* routine  = "Hyb3EQ::Print";
    int i,k;

    map<string, ZeroBond>::iterator zero;
    map<string, NamedZeroBank>::iterator zb;    

    if (treeDataDebugFile.empty())
        return;

    FILE *stream = fopen(treeDataDebugFile.c_str(), "w");
    if (stream == NULL)
    {
        throw ModelException(routine, "Unable to open tree debug info file  (" + 
                             treeDataDebugFile + ") for writing");
    }
    int count=0;

    fprintf(stream, "Following critical dates registered with the engine:\n\n");
    for (i = 0; i < (int)sortedCritDates.size(); i++)
    {
        fprintf(stream, "%d   %ld\n", i, sortedCritDates[i]);
    }

    fprintf(stream, "Domestic Currency %s\n\n", mCurrency[0].c_str());
    fprintf(stream, "Nb  DomesticCurves\n");
    for (i = 0; i < 3; i++)
    {
        fprintf(stream, "%d  %s\n", 
                i,
                mCurveName[0][i].c_str());
    }

    fprintf(stream, "\nFollowing zero banks registered with the engine:\n\n");
    for (zb = mZeroBank.begin(); zb != mZeroBank.end(); zb++)
    {
        fprintf(stream, "%d, name: %s, internalName: %s, nbZeroInBank %d, zeroBankDim %d, currency %s\n",
                count,
                zb->first.c_str(),
                zb->second.name.c_str(),
                zb->second.nbZeros,
                zb->second.zeroBankDim,
                zb->second.currency.c_str());
        count++;
    }

    count = 0;
    fprintf(stream, "\n\nFollowing zero bond/discounting slices registered with the engine:\n\n");
    for (zero = mZeroBond.begin(); zero != mZeroBond.end(); zero++)
    {
        fprintf(stream, "%d, name: %s, nbZeros(ie. non-overlapping) %d, sliceDim %d, currency %s\n",
                count,
                zero->first.c_str(),
                zero->second.startDates.size(),
                zero->second.sliceDim,
                zero->second.currency.c_str());

        fprintf(stream, "   Non-overlapping zero start/end dates stored in zeroIndexSpec\n");
        for (i = 0; i < zero->second.startDates.size(); i++)
        {
            fprintf(stream, "   %d  %ld  %ld\n", i, 
                    zero->second.startDates[i].toIrDate(),
                    zero->second.matDates[i].toIrDate());
        }
        count++;
    }

    fprintf (stream, "\n\n");

    // Slice sizes used for allocation
    fprintf(stream, "  W1   HW1    W2   HW2    W3   HW3\n");
    fprintf(stream, "%5d %5d %5d %5d %5d %5d\n\n\n",
                    mTreeData.Width[0], mTreeData.HalfWidth[0],
                    mTreeData.Width[1], mTreeData.HalfWidth[1],
                    mTreeData.Width[2], mTreeData.HalfWidth[2]);

    // Schedule of the underlying swap 
    fprintf (stream,"\nTIMELINE INFORMATION (IR)\n");

    fprintf (stream,
            "Node      Date     Days  Max  Forward   "
            "Zero0   Discount0   Zero1   Discount1   "
            "Zero2   Discount2    IrMidNode    SpotVol \n");

    for (i = 0; i <= mTreeData.NbTP; i++)
    {
        double Forward = 100. * mTreeData.FwdRate[0][0][i]/mTreeData.Length[i];

        int daysL = Daysact (mTreeData.TPDate[0], mTreeData.TPDate[i]);

        fprintf(stream,
                "[%3d] \t%8ld  %3d  %3d  %7.4f  %7.4f   %8.6f  "
                "%7.4f   %8.6f  %7.4f   %8.6f    %9.6f    %6.2f \n",
                i,
                mTreeData.TPDate[i],
                daysL,
                mTreeData.Top1[i],
                Forward,
                mTreeData.ZeroRate[0][0][i] * 100.,
                mTreeData.ZeroCoupon[0][0][i],
                mTreeData.ZeroRate[0][1][i] * 100.,
                mTreeData.ZeroCoupon[0][1][i],
                mTreeData.ZeroRate[0][2][i] * 100.,
                mTreeData.ZeroCoupon[0][2][i],
                exp(mTreeData.IrZCenter[0][i]),
                mTreeData.SpotVol[0][i] * 100.);
        fflush(stream);
    }

    fprintf (stream, "\n\n");
    fprintf (stream,"\nTIMELINE INFORMATION (Stock)\n");
    fprintf (stream, "Node \t  Date \t\t   Forward      "
                     "Center      Eq vol  Spotvol   Rho's (Irf-Eq, Ird-Eq, FX-Eq)   "
                     "Dividend   DivAmt   SettlDate\n");

    for (i = 0; i <= mTreeData.NbTP; i++)
    {
        fprintf (stream, "[%3d] \t%8ld \t%10.4f \t%10.4f  \t"
                         "%5.2f \t%5.2f \t     %4.2f       %4.2f      %4.2f     "
                         "(%d)     %7.3f    %8ld\n",
                 i,
                 mTreeData.TPDate[i],
                 mTreeData.FwdEq[i],
                 mTreeData.EqMidNode[i],
                 mTreeData.EqVol[i] * 100.,
                 mTreeData.SpotEqVol[i] * 100.,
                 mTreeData.Rho[3][i],
                 mTreeData.Rho[4][i],
                 mTreeData.Rho[5][i],
                 mTreeData.TPType[12][i],
                 mTreeData.CritDate[12][i].Value[0],
                 mTreeData.NodeSettleDate[i]);
    }

    // Print out input zero curves
    fprintf (stream, "\n\n");
    fprintf (stream, 
            "%s Currency Curves (Index, COF, Risk)\n\n",mCurrency[0].c_str());
    for (k = 0; k < 3; k++)
    {
        fprintf (stream, "Maturity      Zero       Discount \n");

        for (i = 0; i < (mTCurves[0][k]).NbZero; i++)
        {

            double days = Daysact((mTCurves[0][k]).ValueDate, 
                                (mTCurves[0][k]).ZeroDate[i]);
            // Discount factor up to the current date
            double discount = ::pow(1. + (mTCurves[0][k]).Zero[i], -days/365.);

            fprintf (stream,
                    "%ld   %9.6f     %8.6f \n",
                    (mTCurves[0][k]).ZeroDate[i],
                    (mTCurves[0][k]).Zero[i] * 100.,
                    discount);
        }
        fprintf (stream, "\n");
    }

    // print out raw hyb3 structures
    fprintf(stream, "\nMKTVOL_DATA structure for IR\n\n");
    MktVol_PrintStructure(stream, &mMktVolData[0]);

    fprintf(stream, "\nTREE_DATA structure\n\n");
    PrintTree(stream, mTreeData, 0);

    fclose(stream);
}

void Hyb3EQ::getEQIndex(TreeSlice& treeSlice)
{
    const char* routine = "Hyb3EQ::getEQIndex";
    TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);

    if (slice.getDim() != 2)
        slice.allocDim(2);

    // for now simply copy the FX Spot slice into the provided slice
    if (Hyb3_CopySlice(slice.getValuePtr(),
                       mDevData.EqSpot,
                       2,       // FX always 3D slice
                       mTpIdx,  // current step
                       &mTreeData) != SUCCESS)
    {
        throw ModelException(routine, IrConverter::slog.pop());
    }
    slice.treeStep = range->treeStep;
}

Hyb3EQ::Hyb3EQ() : Hyb3(TYPE), eqSpotVolOverride(0) {}
Hyb3EQ::~Hyb3EQ() {}


DRLIB_END_NAMESPACE
