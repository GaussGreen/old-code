//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : Fix3.cpp
//
//   Description : Fix3 tree
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Fix3.hpp"
#include "edginc/IRVolPair.hpp"
#include "edginc/IrConverter.hpp"
#include "edginc/IndexSpecIR.hpp"
#include "edginc/IndexSpecFX.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/Format.hpp"
#include "edginc/RatesUtils.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/MDFUtil.hpp"
#include "esl_log.h"

#include <stdio.h>


DRLIB_BEGIN_NAMESPACE


/************************************** Fix3 *************************************/

void Fix3::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(Fix3, clazz);
    SUPERCLASS(RateTree);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(nbFactors,"Number of interest rate factors for tree");
    FIELD(IRParams, "Collection of interest rate model parameters")
    FIELD(treeDataDebugFile, "Print tree debug data to named output file");
    FIELD_MAKE_OPTIONAL(treeDataDebugFile);
    FIELD(zeroInterpStyle, "zero curve interpolation style.  Defaults to "
                                  "FLAT_FWD if not supplied");
    FIELD_MAKE_OPTIONAL(zeroInterpStyle);
    FIELD(zeroBankMode, "internal tree zero bank mode - defaults to"
                               "to ZEROBANK mode if not supplied");
    FIELD_MAKE_OPTIONAL(zeroBankMode);
    FIELD(cetSmoothing,"");
    FIELD_MAKE_OPTIONAL(cetSmoothing);
}

CClassConstSP const Fix3::TYPE = CClass::registerClassLoadMethod(
    "Fix3", typeid(Fix3), Fix3::load);


void Fix3::getMarket(const MarketData*  market, IInstrumentCollectionSP instruments)
{
    static const string routine = "Fix3::getMarket";
    try 
    {
        mToday = market->GetReferenceDate().toIrDate();
        IRParams->getData(this, market);  // required for market wrapper

        // Get Engine
        if (!engineTable.isEmpty())
            engineTable.getData(this, market);

/*        // top level instrument must supply discount factor curve in order to define the
        // pricing/PnL currency of the trade.  This is retrieved through a complicated tangle
        // of code that eventually populates a marketDataFetcher by KComponent implementing the
        // IInstrumentCollection interface function:
        // virtual string discountYieldCurveName() const = 0
        // which returns the name of the discount curve defined in the instrument
        string discYCName = getDomesticYCName();
        if (discYCName.empty())
            throw ModelException("instrument must supply a discount Curve");

        if (!!IRParams->curveToDiscount)  // is optional parameter
        {
        string modelDiscYCName = IRParams->curveToDiscount->getName();
        if (discYCName != modelDiscYCName)
            throw ModelException("Name of curveToDiscount field set in model " + modelDiscYCName +
                                " must be same as defined on the instrument " + discYCName);

        discYC = IRParams->curveToDiscount.getSP();
*/
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
        if (!discYC)
            throw ModelException("Instrument does not contain discount YieldCurve type - "
                                 "unable to populate FDModel::discYC field");

        if (!!IRParams->curveToDiscount)  // if supplied check same as instrument discYC
        {
            string discYCName = discYC->getName();
            string modelDiscYCName = IRParams->curveToDiscount->getName();
            if (discYCName != modelDiscYCName)
                throw ModelException("Name of curveToDiscount field set in Fix3 model " + modelDiscYCName +
                                    " must be same as defined on the instrument " + discYCName);
        }

        // recurse the instrument components to find and store all the IMarketFactor types found
        // in the exported fields.
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

        // recurse the instrument components to find and store all the IndexSpec types found
        // in the exported fields.
        class RetrieveIndexSpecs : public ObjectIteration::IAction
        {
            IndexSpecArray & indexSpecs;

        public:
            RetrieveIndexSpecs(IndexSpecArray & indexSpecs) :
                indexSpecs( indexSpecs ) { indexSpecs.clear(); }

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

    }
    catch (exception& e) {
        throw ModelException(e, routine);
    }
}

void Fix3::initTreeData(void)
{
    try
    {
        // ??? for now reuse this OverWriteString mechanism until final model/market rewrite
        char OverWriteString[6][MAXBUFF];
        RateTree::IRModelParams mp;

        mTreeData.NbFactor = nbFactors;

        // retrieve model parameters from market cache
        IRExoticParamTable* pSmileTable = IRParams->smileTable.get();
        IRExoticParamTable* pModelTable = IRParams->modelTable.get();
        if (pSmileTable || pModelTable)
        {
            // Use new method using MarketTable objects
            IrConverter::to_MODELPARAMETERS_DATA(mp, 
                                                IRParams->smileSet, 
                                                IRParams->modelSet,
                                                engineSet,
                                                *(engineTable.get()),
                                                *pSmileTable,
                                                *pModelTable);
            if (nbFactors != mp.nbFactors)
                throw ModelException(" - mismatch in nbFactors on IR Model object and Fix3 object!");
        }
        else
        {
            if (nbFactors != 1)
                throw ModelException("Currently only support fix3 in 1 factor mode");

            // Use deprecated IRCalib object 
            IrConverter::to_MODELPARAMETERS_DATA(mp, 
                                                IRParams->smileSet, 
                                                IRParams->modelSet,
                                                *(IRParams->irCalib.get()));
        }
        mTreeData.NbSigmaMax = mp.nbStdDevs;

        // ppy
        sprintf(OverWriteString[0], "%d", mp.PPY);
        // smile
        sprintf(OverWriteString[1], "%lf %lf %lf %d",mp.QLeft, mp.QRight, mp.FwdShift, mp.CetNbIter);

        switch (nbFactors)
        {
            case 1:
            // factor weights
            sprintf(OverWriteString[2], "%lf", mp.FactorVol[0]);
            // mean reversion
            sprintf(OverWriteString[3], "%lf", mp.FactorMR[0]);
            // correlation override
            sprintf(OverWriteString[4], "nil");  // ??? N/A for 1F mode, but need to add for 2/3 factor support
            break;

            case 2:
            // ..as above
            sprintf(OverWriteString[2], "%lf %lf", mp.FactorVol[0], mp.FactorVol[1]);
            sprintf(OverWriteString[3], "%lf %lf", mp.FactorMR[0], mp.FactorMR[1]);
            sprintf(OverWriteString[4], "%lf", mp.FactorCorr[0]);
            break;

            case 3:
            // ..as above
            sprintf(OverWriteString[2], "%lf %lf %lf", mp.FactorVol[0], mp.FactorVol[1], mp.FactorVol[2]);
            sprintf(OverWriteString[3], "%lf %lf %lf", mp.FactorMR[0], mp.FactorMR[1], mp.FactorMR[2]);
            sprintf(OverWriteString[4], "%lf %lf %lf", mp.FactorCorr[0], mp.FactorCorr[1], mp.FactorCorr[2]);
            break;
        }

        // backbone 
        sprintf(OverWriteString[5], "%lf", mp.backBone);

        // this follows the convention of the allocation of curves in retrieveFactors
        mTreeData.CvDiff = 0;
        mTreeData.CvIdx1 = 1;
        mTreeData.CvIdx2 = 2;
        mTreeData.CvDisc = 1;

        if (Fix3_Param_Input (&mMktVolData, &mTreeData, mTreeData.NbFactor, OverWriteString,"") != SUCCESS)
        {
            throw ModelException("Fix3_Param_Input falied: "+IrConverter::slog.pop());
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

//------------------------------------------------------------------------------
// Need mapping between literal curve name and numeric id
//------------------------------------------------------------------------------
int Fix3::getCrvIdx(const string& curveName) const
{
    if(curveName.empty())
        throw ModelException(__FUNCTION__, "Curve id not supplied");
    
    int i, index=-1;
    for (i = 0; i < 3; i++)
    {
        if (curveName == mCurveName[i]) {
            index=i;
            break;
        }
    }
    if (index == -1)
        throw ModelException(__FUNCTION__, "Unable to find curveName " +
                             curveName + 
                             " in list of fix3 registered curve names");
    
    return index;
}

//------------------------------------------------------------------------------
// Build tree using market and engine params
//------------------------------------------------------------------------------
Fix3::Fix3(const CClassConstSP &type) :
    RateTree(type), 
    mToday(0), mValueDate(0), 
    zeroBankMode(ZeroBankMode::ZEROBANK), 
    mTreeBuilt(false), 
    nbFactors(0),
    zeroInterpStyle(ZeroInterpStyle::FLAT_FWD),
    cetSmoothing(false)
{
    memset(&mMktVolData, 0, sizeof(mMktVolData));
    memset(&mIRSim, 0, sizeof(mIRSim));
    memset(&mTreeData, 0, sizeof(mTreeData));
    memset(&mDevData, 0, sizeof(mDevData));
    memset(&mTCurves, 0, sizeof(mTCurves));

    Fix3_Tree_Init( &mTreeData );
    Fix3_Dev_Init( &mDevData );
}

//------------------------------------------------------------------------------
// Destroy the tree
//------------------------------------------------------------------------------
Fix3::~Fix3() {
    clear();
}

//------------------------------------------------------------------------------
// Register required ParYield on given date with the tree
// if in CLAIMBANK mode, update the claim bank, otherwise simply
// adds the reset date to the critical date list and expects the associated
// zerobank to be inserted with a separate function call
//------------------------------------------------------------------------------
void Fix3::insertIRIndex(const IndexSpecIR& rate, DateTime date)
{
    int status = FAILURE;
    int i;

    // ignore rate resets in the past
    long resDate = date.toIrDate();
    if (resDate < getValueDate().toIrDate())
        return;

    // nothing more to do if in zeroBank mode, as zeroDates are inserted
    // through the namedZeroBank function
    if (zeroBankMode == ZeroBankMode::ZEROBANK)
        return;

    // all that follows here is optimized/unoptimized claimBank logic
    int nbIdxZMatDate = 0;
    int nbIdxZUseDate = 0;

    long* IdxZMatDL = NULL;
    long* IdxZUseDL = NULL;
    bool isCrit = true;

    string curveName = rate.getFactor()->getName();
    int tenorInMonths = rate.tenor->toMonths();
    int rateFreq  = rate.frequency->annualFrequency();

    long resetDate[1];
    resetDate[0] = date.toIrDate();
    long resetEffDate[1]; 
    resetEffDate[0] = resetDate[0];

    map<string, NamedClaimBank>::iterator cb;

    // Generate all cashflow dates
    if (ZbkDLFromIdx(1,
                     resetDate,
                     resetEffDate,
                     tenorInMonths,
                     rateFreq,
                     &nbIdxZMatDate,
                     &IdxZMatDL,
                     &nbIdxZUseDate,
                     &IdxZUseDL) != SUCCESS)
        goto RETURN;

    // Add to zero bank dates
    
    if (isCrit)
    {
        for (i=0; i<nbIdxZUseDate; i++)
        {
            cb->second.critZeroUseDates.push_back(IdxZUseDL[i]);
            cb->second.critZeroMatDates.push_back(IdxZMatDL[i]);
        }
    }
    else
    {
        for (i=0; i<nbIdxZUseDate; i++)
        {
            cb->second.optZeroUseDates.push_back(IdxZUseDL[i]);
            cb->second.optZeroMatDates.push_back(IdxZMatDL[i]);
        }
    }

    status = SUCCESS;

    RETURN:
    if (status != SUCCESS)
        throw ModelException(__FUNCTION__, IrConverter::slog.pop());
}


//------------------------------------------------------------------------------
// Insert directly into the tree the zero bank dates and support fields 
// that have been calculated by the product when the tree is running in
// zeroBank mode
//------------------------------------------------------------------------------
void Fix3::insertNamedZeroBank(NamedZeroBank& zeroBank)
{
    const char* routine = "Fix3::insertNamedZeroBank";
    if (zeroBankMode == ZeroBankMode::CLAIMBANK)
        throw ModelException(routine, "Tree must be in ZEROBANK mode to insert " 
                             "client calculated zero bank");

    // simply insert the zero -map ensures a unique name
    if (mZeroBank.find(zeroBank.name)!=mZeroBank.end()) {
        return; //!!! temporary
        throw ModelException(routine,
        "Trying to insert two ZeroBanks for the name IndexSpecIR="+zeroBank.name);
    }
    mZeroBank[zeroBank.name] = zeroBank; 
}


//------------------------------------------------------
// Register the zero discounting required by the pricer
//------------------------------------------------------
void Fix3::insertZero(const ZeroBond& zero)
{
    const char* routine = "Fix3::insertZero";
    int i;

    if (zero.startDates.size() != zero.matDates.size())
        throw ModelException(routine, "Number of start dates " + 
                             Format::toString(zero.startDates.size()) +
                             " must equal number of maturity dates " +
                             Format::toString(zero.matDates.size()));

    // check all the zero dates are in ascending order
    // don't strictly have to be, but makes things simpler
    for (i = 1; i < zero.matDates.size(); i++)
    {
        if (zero.matDates[i-1] >= zero.matDates[i])
            throw ModelException(routine, "zero mat dates for named zero " +
                                 zero.getName() +
                                 "must be in ascending order");

        if (zero.startDates[i-1] >= zero.startDates[i])
            throw ModelException(routine, "zero start dates for named zero " +
                                 zero.getName() +
                                 "must be in ascending order");
    }

    for (i = 0; i < zero.matDates.size(); i++)
    {
        if (zero.startDates[i] >= zero.matDates[i])
            throw ModelException(routine, "zero start date " + 
                                 zero.startDates[i].toString() + 
                                 " must be before the associated zero maturity date " +
                                 zero.matDates[i].toString());
    }

    // simply insert the zero -map ensures a unique name
    mZeroBond[zero.getName()] = zero;
}

//---------------------------------------------------------------
// Get the zero slice - currently assume that the current date
// must map to an existing zero startDates - otherwise it's a bit
// difficult to figure out which zero slice we want
// ?? this doesn't work generally as we could have 2 zeros starting
// on the same day and maturing at different dates, but rethink
// this logic later
//---------------------------------------------------------------

void Fix3::getZero(const ZeroBond& zero, int step, TreeSlice& treeSlice)
{
    TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);

    map<string, ZeroBond>::iterator zeroIter = mZeroBond.find(zero.getName());
    if (zeroIter == mZeroBond.end())
        throw ModelException(__FUNCTION__, "Unable to find named zero bond " + zero.getName());

    // ?? put logic to find zero slice associated with current zero.  For the 
    // meantime, we only have 1 zero slice so it is assumed it is correctly
    // reset, dev'd etc. - must assumes contiguous coupons for now
    if (slice.getDim() < nbFactors)
        slice.allocDim(nbFactors);

    double* slicePtr = slice.getValuePtr();

    // ?? use this for now until sorting out use of TreeSlice internally
    if (Fix3_Copy_Slice(slicePtr,
                        zeroIter->second.zeroSlice,
                        step,
                        &mTreeData) != SUCCESS)
    {
        throw ModelException(__FUNCTION__, IrConverter::slog.pop());
    }
    slice.treeStep = range->treeStep;
}

void Fix3::registerZero(
        DateTime obsDate, 
        DateTime matDate, 
        string curveName)
{
    try {
        if (matDate < obsDate) {
            throw ModelException("matDate < obsDate");
        }
        DateTimeArray startList, endList;
        startList.push_back(obsDate);
        endList.push_back(matDate);
        ZeroBondProdSP zeroProd = DYNAMIC_POINTER_CAST<ZeroBondProd>(createProduct(IProdCreatorSP(new ZeroBond(
            startList, endList, curveName))));

        string key = DateTime::dateFormat(obsDate.getDate())+"_"+DateTime::dateFormat(matDate.getDate());
        ZeroProdMap::iterator it = zeroProdList.find(key);
        if (it == zeroProdList.end())
            zeroProdList.insert(pair<string, ZeroBondProdSP>(key, zeroProd));
    }
    catch (exception& e) { 
        throw ModelException(e, __FUNCTION__); 
    }
}


void Fix3::getZero(
        DateTime obsDate, 
        DateTime matDate, 
        string curveName, 
        TreeSliceSP &slice)
{
    try {
        if (!dynamic_cast<TreeSliceRates*>(slice.get())) {
            slice = RateTree::createSimpleSlice(curveName);
            slice->name = "Zero "+obsDate.toString()+" to "+matDate.toString();
        }

        if (obsDate==matDate) {
            *slice = 1.;
            return;
        }
        string key = DateTime::dateFormat(obsDate.getDate())+"_"+DateTime::dateFormat(matDate.getDate());
        ZeroProdMap::iterator it = zeroProdList.find(key);
        if (it == zeroProdList.end())
            throw ModelException("ZeroProd "+key+" not found");

        *slice = it->second->getValue(range->treeStep, matDate);
    }
    catch (exception& e) { 
        throw ModelException(e, __FUNCTION__); 
    }
}


void Fix3::getDiffusionTCurve(T_CURVE& diffusionCurve) const
{
    IrConverter::to_T_CURVE(diffusionCurve, IRParams->curveToDiffuse.get(),
                            false);
}


/** retrieving market data */
void Fix3::retrieveFactor(){
    static const string method = "Fix3::retrieveFactor";
    try{
        if (!discYC)
            throw ModelException("FDModel::discYC field is not populated - must contain the product's "
                                 "discount curve. Internal error in FDModel/product code");

        string mCurrency = discYC->getCcy();

        if (!IRParams->curveToDiffuse)
            throw ModelException("Internal error - Model curve to diffuse not populated by market object");

        // check curves are same currency
        if (IRParams->curveToDiffuse->getCcy() != mCurrency)
            throw ModelException("Currency of model curve to diffuse " + 
                                 IRParams->curveToDiffuse->getCcy() +
                                 " must be the same as currency of instrument " +
                                 mCurrency);

        // for consistency, check model curve to discount if supplied is that same as
        // the instrument curve to discount
        if (!IRParams->curveToDiscount.isEmpty())
        {
            if (IRParams->curveToDiscount->getName() != discYC->getName())
                throw ModelException("model curveToDiscount supplied (name = " +
                                     IRParams->curveToDiscount->getName() +
                                     ") is not same name as discount curve defined by the product " +
                                     discYC->getName());
        }
        // by convention set 0=diff curve, 1=discount curve, 2=other
        // assume vols are associated with diffusion curve
        mCurveName[0] = IRParams->curveToDiffuse->getName();
        getDiffusionTCurve(mTCurves[0]);

        mCurveName[1] = discYC->getName();
        IrConverter::to_T_CURVE(mTCurves[1], discYC.get(), false);

        // all factors supplied must be yield curves and in same currency for fix3.
        // fix3 only supports max 3 curves, so check if this limit is exceeded
        bool thirdCurve=false;
        int i;
        for (i = 0; i < factors.size(); i++)
        {
            IMarketFactorConstSP &mf = factors[i];
            const YieldCurve* yc = dynamic_cast<const YieldCurve*>(mf.get());\

            if (!yc)
                throw ModelException("Fix3 only supports yield curve type market objects - "
                                     "type supplied = " +
                                     mf->getClass()->getName());
            // now check currency
            if (yc->getCcy() != mCurrency)
                throw ModelException("Yield curve defined in instrument of currency " +
                                     yc->getCcy() +
                                     ". Fix3 is single currency model in currency - " +
                                     mCurrency);
            // check if third curve (ie. not discount or diffusion)
            if (yc->getName() == mCurveName[0] ||
                yc->getName() == mCurveName[1])
                continue;
            else
            {
                // error if thirdCurve is already defined
                if (thirdCurve)
                    throw ModelException("Fix3 engine supports maximum of three yield curves "
                                         "and supplied yield curve/factor " + yc->getName() +
                                         "is a fourth.  Current 3 curves registered are " +
                                         mCurveName[0] + ", " + mCurveName[1] + ", " +
                                         mCurveName[2]);
                mCurveName[2] = yc->getName();
                IrConverter::to_T_CURVE(mTCurves[2], yc, false);
                thirdCurve = true;
            }
        }
        if (thirdCurve == false)
        {
            // assign 3rd curve same as second curve - rates default behaviour
            mCurveName[2] = mCurveName[1];
            mTCurves[2] = mTCurves[1];
        }

        // now check that the curves have the same value date, and set this accordingly in the tree
        mValueDate = mTCurves[0].ValueDate;
        for (i = 1; i < 3; i++)
        {
            if (mTCurves[i].ValueDate != mValueDate)
                throw ModelException("ValueDate " + 
                                     DateTime::fromIrDate(mTCurves[i].ValueDate).toString() +
                                     " defined in zero curve index [" + Format::toString(i) + "] is not the same "
                                     "as the valueDate defined in the other zero curve(s) " +
                                     DateTime::fromIrDate(mValueDate).toString());
        }

        IRVolSelector volSelector(getIRVolRaw(IRParams), mTCurves[0],
                                  IRParams->volCalibIndex,
                                  cetSmoothing, IRParams);
        volSelector.getMktVolData(mMktVolData);
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

// Implementation of IRVegaPointwise::ISensitivePoints for smart vega tweaking
// and volatility exposure reporting.
IRGridPointAbsArraySP Fix3::getSensitiveIRVolPoints(
    OutputNameConstSP  outputName,
    const CInstrument* inst) const
{
    IRGridPointAbsArraySP volExposuresSP(new IRGridPointAbsArray);

    T_CURVE diffusionCurve;
    getDiffusionTCurve(diffusionCurve);

    IRVolSelector volSelector(getIRVolRaw(IRParams), diffusionCurve,
                              IRParams->volCalibIndex,
                              cetSmoothing, IRParams);
    volSelector.getVolExposures(volExposuresSP, outputName);
 
    return volExposuresSP;
}


/**  collect model initialisation data, set up timeline  */
void Fix3::initModel(void)
{
    int status = FAILURE;
    mTreeBuilt=true;

    EslSetZeroInterpolation(zeroInterpStyle == ZeroInterpStyle::LINEAR ?
        ESL_INTERP_LINEAR : ESL_INTERP_FLATFWD);
    
    CRIT_DATE* critDate=NULL;
    int nbCritDate = 0;
    IrConverter::AutoDeleteDRArray<CRIT_DATE>
        critDateDelete(&critDate, CRITDATE, &nbCritDate);

    int i, j;

    map<string, ZeroBond>::iterator zero;
    map<string, NamedZeroBank>::iterator zb;

    initTreeData();

    // needed to setup fix3 timeline structures - have to register critical dates 
    // that are > valuedate, as dates are dropped on <= value date
    critDate = (CRIT_DATE *) DR_Array (CRITDATE, 0, 0);
    if (critDate == NULL)
        goto RETURN;

    long treeStartDate;
    if (zeroBankMode == ZeroBankMode::ZEROBANK)
        treeStartDate = getValueDate().toIrDate();
    else
        treeStartDate = mToday;

    // add tree startDate as critical date
    if (Add_To_DateList(&nbCritDate, 
                        &critDate,
                        getValueDate().toIrDate(),
                        0, 0, 0, 0, 0, 0, 0, 0, 0) != SUCCESS)
        goto RETURN;

    arrangeDates();

    // add product specific critical dates
    for (j = 0; j < critDates.size(); j++)
    {
        DateTime critDateL = critDates[j];
        if (critDateL > getValueDate())
        {
            if (Add_To_DateList(&nbCritDate, 
                                &critDate,
                                critDateL.toIrDate(),
                                0, 0, 0, 0, 0, 0, 0, 0, 0) != SUCCESS)
                goto RETURN;
        }
    }

    // if in CLAIMBANK mode, construct and insert critical zero dates
    if (zeroBankMode == ZeroBankMode::CLAIMBANK)
    {
        DR_Error("Claim bank not yet suppported");
        goto RETURN;
    }

    // construct the tree timeline
    if (Fix3_Time_Line (getValueDate().toIrDate(),
                        nbCritDate,
                        critDate,
                        'I',                            
                        &mTreeData) != SUCCESS)
    {
        goto RETURN;
    }

    // numerical volatility calibration routine
    if( Fix3_Cet_Main(FALSE,
                      mTCurves,
                      &mMktVolData,
                      &mTreeData ) != SUCCESS )
    {
        DR_Error("Fix3_Cet_Main failed");
        goto RETURN;
    }

    if( Fix3_Build_Tree(mTCurves,
                        &mMktVolData,
                        &mTreeData ) != SUCCESS )
    {
        DR_Error("Fix3_Build_Tree falied");
        goto RETURN;
    }

    if( mTreeData.NbFactor == 1 )
    {
        range.reset(new TreeSliceRates::Range(
            -mTreeData.HalfWidth[0], 
            mTreeData.HalfWidth[0]));
    }
    else if( mTreeData.NbFactor == 2 )
    {
        range.reset(new TreeSliceRates::Range(
            -mTreeData.HalfWidth[0],
            mTreeData.HalfWidth[0],
            -mTreeData.HalfWidth[1], 
            mTreeData.HalfWidth[1]));
    }
    else if( mTreeData.NbFactor == 3 )
    {
        range.reset(new TreeSliceRates::Range(
            -mTreeData.HalfWidth[0],
            mTreeData.HalfWidth[0],
            -mTreeData.HalfWidth[1], 
            mTreeData.HalfWidth[1],
            -mTreeData.HalfWidth[2], 
            mTreeData.HalfWidth[2]));
    }

    ASSERT(zeroBankMode == ZeroBankMode::ZEROBANK);

    // ?? allocate memory for zero banks
    // loop through each registered zero bank
    for (zb = mZeroBank.begin(); zb != mZeroBank.end(); zb++)
    {
        int nbZeros = zb->second.nbZeros;
        zb->second.zeroCritDates.resize(nbZeros, 0);
        zb->second.zeroBank.resize(nbZeros);
        zb->second.zeroBankDim = mTreeData.NbFactor;

        for (i = 0; i < nbZeros; i++)
        {
            zb->second.zeroBank[i] = Fix3_Alloc_Slice(&mTreeData);
            if (zb->second.zeroBank[i] == NULL)
            {
                DR_Error("Unable to allocate memory for zero banks");
                goto RETURN;
            }
        }
    }

    // loop through each stored named zero bank
    for (zero = mZeroBond.begin(); zero != mZeroBond.end(); zero++)
    {
        zero->second.zeroSlice = Fix3_Alloc_Slice(&mTreeData);
        zero->second.sliceDim = mTreeData.NbFactor;
        if (zero->second.zeroSlice == NULL)
        {
            DR_Error("Unable to allocate memory for named zero bond %s",
                        zero->second.getName().c_str());
            goto RETURN;
        }
    }
    print();

    status = SUCCESS;

    RETURN:

    if (status != SUCCESS)
        throw ModelException(__FUNCTION__, IrConverter::slog.pop());
}



/**  finalise model initialisation   */
void Fix3::finaliseModel(CControl* control)
{
    if( Fix3_Dev_Alloc( &mDevData, &mTreeData ) != SUCCESS )
        throw ModelException( "Fix3::finaliseModel", "Fix3_Dev_Alloc falied" );

    mTpIdx = mTreeData.NbTP;
}

//------------------------------------------------------------------------------
// Get today's date
//------------------------------------------------------------------------------
DateTime Fix3::getToday() const
{
    DateTime today = RateTree::getValueDate();
    return today;
}

int Fix3::getCurveIdx( const string & curveName ) const
{
    return curveName.empty() ? -1 : getCrvIdx( curveName );
}

//------------------------------------------------------------------------------
// Get date on the timeline pointed by the index
//------------------------------------------------------------------------------
DateTime Fix3::getDate(int dateIdx) const
{
    if (dateIdx < 0 || dateIdx > mTreeData.NbTP) {
        throw ModelException(__FUNCTION__, 
        "Invalid time point index (" + Format::toString(dateIdx) 
        + "). Must be in the range 0..." + Format::toString(mTreeData.NbTP));
    }
    return DateTime::fromIrDate(mTreeData.TPDate[dateIdx]);
}

DateTimeArray Fix3::getDates() const
{
    DateTimeArray dates( mTreeData.NbTP + 1 );
    for( int i = 0; i <= mTreeData.NbTP; ++i )
        dates[ i ] = DateTime::fromIrDate( mTreeData.TPDate[ i ] );
    return dates;
}

/** get last step (total num of steps) on time line */
int Fix3::getLastStep() const
{
    return mTreeData.NbTP;
}

//------------------------------------------------------------------------------
// Get value date for given curve
//------------------------------------------------------------------------------
DateTime Fix3::getCurveValueDate(string curveName) const
{
    int idx;

    // if curve is not supplied, default value date is the 1st T_CURVE date
    if (curveName.empty())
        idx = 0;
    else
        idx = getCrvIdx(curveName);

    long valDate = mTCurves[idx].ValueDate;
    return DateTime::fromIrDate(valDate);
}

//------------------------------------------------------------------------------
// Print
//------------------------------------------------------------------------------
void Fix3::print()
{
    static char const* routine  = "Fix3::print";

    int i;

    map<string, ZeroBond>::iterator zero;
    map<string, NamedZeroBank>::iterator zb;

    char   ErrorString[2][MAXBUFF];
    strcpy(ErrorString[0], "foreign");
    strcpy(ErrorString[1], "domestic");

    if (treeDataDebugFile.empty())
        return;

    FILE *stream = fopen(treeDataDebugFile.c_str(), "w");
    if (stream == NULL)
    {
        throw ModelException(routine, "Unable to open file \""+treeDataDebugFile+"\" for writing");
    }
    int count=0;
    fprintf(stream, "Following zero banks registered with the engine:\n\n");
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

    fprintf(stream, "\nDomestic Currency %s\n\n", mCurrency.c_str());

    Fix3_Print_FIX3_TREE_DATA(stream, &mTreeData);
    Print_MKTVOL_DATA(stream, &mMktVolData);

    fclose(stream);
    
}

//------------------------------------------------------------------------------
// Delete parts of the model that have been created by build.
//------------------------------------------------------------------------------
void Fix3::clear()
{
    if (!mTreeBuilt)
        return;
    mTreeBuilt = false;

    ESL_LOG(Log::debug) << __FUNCTION__ << endl;

    Fix3_Dev_Free (&mDevData, &mTreeData);

    ASSERT(zeroBankMode == ZeroBankMode::ZEROBANK);
    { // free zerobank
        map<string, NamedZeroBank>::iterator zb;
        for (zb = mZeroBank.begin(); zb != mZeroBank.end(); zb++)
        {
            vector<double*> &bank = zb->second.zeroBank;

            for (size_t i = 0; i < bank.size(); i++) {
                Fix3_Free_Slice(bank[i], &mTreeData);
                bank[i] = NULL;
            }
            bank.clear();
            zb->second.zeroCritDates.clear();
            zb->second.zeroBank.clear();
        }
        mZeroBank.clear();

        map<string, ZeroBond>::iterator zero;
        for (zero = mZeroBond.begin(); zero != mZeroBond.end(); zero++)
        {
            Fix3_Free_Slice (zero->second.zeroSlice, &mTreeData);
            zero->second.zeroSlice = NULL;
        }
        mZeroBond.clear();
    }

    Fix3_Tree_Free(&mTreeData);
    mTpIdx = -1;
}


//------------------------------------------------------------------------------
// Update engine at each time step
//------------------------------------------------------------------------------
void Fix3::update(int t, FDProduct::UpdateType type)
{
    try {
        if (!mTreeBuilt)
            throw ModelException("Not able to update tree as the tree has "
                                        "not been built");

        // set current range
        range->limits.bot1 = mTreeData.Bottom1[t];
        range->limits.top1 = mTreeData.Top1[t];
        range->limits.bot2 = mTreeData.Bottom2[t];
        range->limits.top2 = mTreeData.Top2[t];
        range->limits.bot3 = mTreeData.Bottom3[t];
        range->limits.top3 = mTreeData.Top3[t];
        range->treeStep = t;

        int T = mTreeData.NbTP;

        /* 
        *  'Update' tree.
        */
        if( Fix3_Lattice(&mDevData,
                        t,
                        T,
                        &mMktVolData,
                        &mTreeData ) != SUCCESS )
        {
            throw ModelException(IrConverter::slog.pop());
        }
        

        mTpIdx = t;

        // update the zerobanks
        ASSERT(zeroBankMode == ZeroBankMode::ZEROBANK);

        map<string, NamedZeroBank>::iterator zb;
        // loop through each stored named zero bank
        for (zb = mZeroBank.begin(); zb != mZeroBank.end(); zb++)
        {
            int devCurve = getCrvIdx(zb->second.curveName);
            if (Fix3_Zero_Bank(&zb->second.zeroBank[0],
                                &zb->second.zeroCritDates[0],
                                &zb->second.nbCurrentZeros,
                                zb->second.nbZeros,
                                0,       // not reset
                                getDate(t).toIrDate(),
                                t,
                                T,
                                devCurve,
                                &mDevData,
                                &mTreeData) != SUCCESS)
            {
                throw ModelException(IrConverter::slog.pop());
            }
        }

        map<string, ZeroBond>::iterator zero;
        // update any of the zero bonds maintained by the tree
        for (zero = mZeroBond.begin(); zero != mZeroBond.end(); zero++)
        {
            int devCurve = getCrvIdx(zero->second.curveName);
            if (Fix3_Zero_t(zero->second.zeroSlice,
                            0,          /* no reset now   */
                            t,
                            T,
                            devCurve,    // ??? these will have to be zerobond slice dependent
                            &mDevData,
                            &mTreeData) != SUCCESS)
            {
                throw ModelException(IrConverter::slog.pop());
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

//--------------------------------------
//
//--------------------------------------

void Fix3::updateFinalize(int t, FDProduct::UpdateType type)
{
    int status = FAILURE;

    DateTime currentDate = getDate(t);

    if (!mTreeBuilt)
        throw ModelException(__FUNCTION__, "Not able to update tree as the tree has "
                                      "not been built");

    // if zero bond reset date, reset zero bond (ie. set slice to 1)
    if (zeroBankMode == ZeroBankMode::ZEROBANK)
    {
        int i;
        map<string, NamedZeroBank>::iterator zb;
        map<string, ZeroBond>::iterator zero;
        // loop through each stored named zero bank
        for (zb = mZeroBank.begin(); zb != mZeroBank.end(); zb++)
        {
            bool zeroMatDate = false;
            int T = mTreeData.NbTP;
            long currentDateL = currentDate.toIrDate();

            // check if today is a zero bank maturity date  
            // ?? could use find but would be superceeded by event driver to be done
            for (i = 0; i < zb->second.zeroMatDateList.size(); i++)
            {
                if (currentDate == zb->second.zeroMatDateList[i])
                {
                    zeroMatDate = true;
                    break;
                }
            }

            if (zeroMatDate)
            {
                if (zb->second.zeroBank.empty())
                    throw ModelException(__FUNCTION__, "Named zero bank " + zb->second.name +
                                                  " zero bank slice collection is not defined");

                int devCurve = getCrvIdx(zb->second.curveName);
                if (Fix3_Zero_Bank(&zb->second.zeroBank[0],
                                   &zb->second.zeroCritDates[0],
                                   &zb->second.nbCurrentZeros,
                                   zb->second.nbZeros,
                                   1,       // reset zerobank slice
                                   currentDateL,
                                   t,
                                   T,
                                   devCurve,
                                   &mDevData,
                                   &mTreeData) != SUCCESS)
                {
                    goto RETURN;
                }
            }
        }

        for (zero = mZeroBond.begin(); zero != mZeroBond.end(); zero++)
        {
            bool zeroMatDate = false;
            int T = mTreeData.NbTP;

            // check if today is a zero bond maturity date  
            for (i = 0; i < zero->second.matDates.size(); i++)
            {
                if (currentDate == zero->second.matDates[i])
                {
                    zeroMatDate = true;
                    break;
                }
            }

            if (zeroMatDate)
            {
                if (Fix3_Zero_t(zero->second.zeroSlice,
                                1,  // force zero bond reset
                                t,
                                T,
                                mTreeData.CvDisc,
                                &mDevData,
                                &mTreeData) != SUCCESS)
                {
                    goto RETURN;
                }
            }
        }
    }

    status = SUCCESS;

    RETURN:

    if (status != SUCCESS)
        throw ModelException(__FUNCTION__, IrConverter::slog.pop());
}

//------------------------------------------------------------------------------
void Fix3::sliceExpand(TreeSlice& treeSlice, const string& curveName) const {
    try {
        TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);
        slice.expand(nbFactors);
    }
    catch (exception& e) {
        throw ModelException(e, "Fix3::sliceExpand");
    }
}

//------------------------------------------------------------------------------
// Perform one backward induction step by moving the slice from t+1 to t
//------------------------------------------------------------------------------
void Fix3::sliceDev(TreeSlice& treeSlice, int curveIdx) const
{
    try {
        TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);

        // nothing to do at the back of the tree 
        if (mTpIdx == mTreeData.NbTP)
            return;

        if (Fix3_Dev (slice.getValuePtr(), mTpIdx, mTreeData.NbTP, curveIdx, &mDevData, &mTreeData)
             != SUCCESS)
        {
            throw ModelException(IrConverter::slog.pop());
        }
    }
    catch (exception& e) {
        throw ModelException(e, "Fix3::sliceDev");
    }

}


//------------------------------------------------------------------------------
// Perform one backward induction step by moving the slice from t+1 to t - no discounting
//------------------------------------------------------------------------------
void Fix3::sliceEv(TreeSlice& treeSlice) const
{
    TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);
    // nothing to do at the back of the tree 
    if (mTpIdx == mTreeData.NbTP)
    {
        return;
    }

    if (slice.getDim() < nbFactors)
        throw ModelException("Fix3::sliceEv()","Wrong dimension for slice "+slice.name);

    Fix3_Ev (slice.getValuePtr(), mTpIdx, mTreeData.NbTP, &mDevData, &mTreeData);
}

/************************************************************/
// for class loading
bool Fix3Load(void) {
    return (Fix3::TYPE != 0);
}

//----------------------------------------------------
//
//----------------------------------------------------
void Fix3::getIRIndex(const IndexSpecIR& rateSpec, DateTime currentDate, TreeSlice& treeSlice)
{
    int status = FAILURE;
    TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);
    double* sliceL=0;

    if (slice.getDim() != nbFactors)
        slice.allocDim(nbFactors);

    string name = rateSpec.getName();
    char dcc = IrConverter::dccTo035A(*rateSpec.dcc);
    char freq = IrConverter::toEslFreq(rateSpec.frequency.get());

    ASSERT(zeroBankMode == ZeroBankMode::ZEROBANK);

    map<string, NamedZeroBank>::iterator zb = mZeroBank.find(name);
    if (zb == mZeroBank.end())
        throw ModelException(__FUNCTION__, "Unable to find named zero bank " + name);

    if (zb->second.zeroBank.empty())
    {
        DR_Error("Internal error - No zero bank slices have been allocated in the internal "
                    "tree named zero bank %s", name.c_str());
        goto RETURN;
    }
    if (zb->second.zeroCritDates.empty())
    {
        DR_Error("Internal error - No zero bank critical dates have been allocated in the "
                    "internal named zero bank %s", name.c_str());
        goto RETURN;
    }

    sliceL = slice.getValuePtr();

    if (Fix3_Par_Yield_t(sliceL,
                            zb->second.nbCurrentZeros,
                            &zb->second.zeroBank[0],
                            &zb->second.zeroCritDates[0],
                            1,  // definite reset
                            currentDate.toIrDate(),
                            currentDate.toIrDate(),  // spot starting
                            rateSpec.tenor->toMonths(),
                            dcc, 
                            freq,
                            0,  // float coupon spread = 0
                            mTpIdx,
                            mTreeData.NbTP,
                            mTreeData.CvDisc,
                            &mDevData,
                            &mTreeData) != SUCCESS)
    {
        goto RETURN;
    }

    slice.treeStep = range->treeStep;

    status = SUCCESS;

    RETURN:

    if (status != SUCCESS)
        throw ModelException(__FUNCTION__, IrConverter::slog.pop());

}


// Here we always return a slice of value 1.0 in order to allow product
// logic to be generalized
void Fix3::getFXIndex(TreeSlice& treeSlice)
{
    const char* routine = "Fix3::getFXIndex";
    try {
        TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);

        if (slice.getDim() != nbFactors)
            slice.allocDim(nbFactors);

        slice = 1.0;
        slice.treeStep = range->treeStep;
    }
    catch (exception& e) {
        throw ModelException(e, routine);
    }
}


/** Create a MarketDataFetcher which will be used for retrieving market data etc */
MarketDataFetcherSP Fix3::createMDF() const {
    MarketDataFetcherSP mdf(new RateTreeMDF(nbFactors));
    MDFUtil::setUseCurrencyBasis(*mdf, true);
    return mdf;
}

DRLIB_END_NAMESPACE
