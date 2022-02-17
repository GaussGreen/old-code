//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : Fix3CB.cpp
//
//   Description : temporary version of fix3 to work on the claim/bank and
//                 tree starting today without cluttering up the existing Fix3
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Fix3CB.hpp"
#include "edginc/IRVolPair.hpp"
#include "edginc/IrConverter.hpp"
#include "edginc/RatesUtils.hpp"
#include "edginc/IndexSpecIR.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/Format.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/Results.hpp"
#include "edginc/BenchmarkDate.hpp"
#include "edginc/MDFUtil.hpp"
#include "esl_log.h"

#include <stdio.h>

DRLIB_BEGIN_NAMESPACE


/************************************** Fix3CB *************************************/

void Fix3CB::load(CClassSP& clazz) {

    clazz->setPublic();
    REGISTER(Fix3CB, clazz);
    SUPERCLASS(RateTree);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(nbFactors,"Number of interest rate factors for tree");
    FIELD(IRParams, "Collection of interest rate model parameters")
    FIELD(treeDataDebugFile, "Print tree debug data to named output file");
    FIELD_MAKE_OPTIONAL(treeDataDebugFile);
    FIELD(zeroInterpStyle, "zero curve interpolation style.  Defaults to "
                                  "FLAT_FWD if not supplied");
    FIELD_MAKE_OPTIONAL(zeroInterpStyle);

    FIELD(cetSmoothing,"");
    FIELD_MAKE_OPTIONAL(cetSmoothing);

    // transient fields - copied but not exported
    FIELD(today, "today/reference date");
    FIELD_MAKE_TRANSIENT(today);
    FIELD(valueDate, "actual value date");
    FIELD_MAKE_TRANSIENT(valueDate);
    FIELD(curveValueDate, "spot/value date of the curve");
    FIELD_MAKE_TRANSIENT(curveValueDate);
    IrConverter::checkSlog();
}

// constructor
Fix3CB::Fix3CB(const CClassConstSP &type) :
RateTree(type), treeBuilt(false), nbFactors(0), 
zeroInterpStyle(ZeroInterpStyle::FLAT_FWD), 
cetSmoothing(false) 
{

    memset(&treeData,0,sizeof(treeData));
    memset(&devData,0,sizeof(devData));
    memset(treeCurves,0,sizeof(treeCurves));
    memset(&mktVolData,0,sizeof(mktVolData));
    memset(&today2ValDateZeros,0,sizeof(today2ValDateZeros));
    Fix3_Tree_Init(&treeData);
    Fix3_Dev_Init(&devData);
}

// destructor
Fix3CB::~Fix3CB() {
    clear();
}

// interface reflection/class loading mechanism
CClassConstSP const Fix3CB::TYPE = CClass::registerClassLoadMethod(
    "Fix3CB", typeid(Fix3CB), Fix3CB::load);

void Fix3CB::getMarket(const MarketData*  market, IInstrumentCollectionSP instruments) {
    try {
        today = market->GetReferenceDate();
        // ??? need final decision by business what to do here
        // set valueDate = today for now, but today may be SOD or EOD but we
        // want valueDate to be EOD no matter what time today is - this means
        // payments are always dropped anytime today
        valueDate = DateTime(today.getDate(), DateTime::END_OF_DAY_TIME);
        
        int time = today.getTime();

        // tree only supports today as SOD or EOD - no other times
        if (time != DateTime::START_OF_DAY_TIME &&
            time != DateTime::END_OF_DAY_TIME) 
            throw ModelException("Model only supports today date times of START_OF_DAY_TIME "
                                 "and END_OF_DAY_TIME.  DateTime time received = " +
                                 Format::toString(time));

        IRParams->getData(this, market);  // required for market wrapper

        // Get Engine
        if (!engineTable.isEmpty())
            engineTable.getData(this, market);

        // recurse the supplied instruments to retrieve the first YieldCurve object
        // type encountered.  This should be the top level instrument discount
        // curve object.
        class RetrieveYC : public ObjectIteration::IAction {
            const string & name;
            YieldCurveConstSP & yc;
        public:
            RetrieveYC(const string & name, YieldCurveConstSP& yc) :
                name(name),
                yc(yc)
            {}

            virtual bool invoke(ObjectIteration::State& state, IObjectConstSP obj) {
                const YieldCurveConstSP & ryc = YieldCurveConstSP::dynamicCast(obj);
                if (! yc && name == ryc->getName())
                    yc = ryc;

                // don't recurse inside the YieldCurve
                return false;
            }
        };
        RetrieveYC retrieveYC(getDomesticYCName(), discYC);
        ObjectIteration retrieveYCIter(YieldCurve::TYPE);
        retrieveYCIter.recurse(retrieveYC, instruments);

        if (!discYC)
            throw ModelException("Instrument does not contain discount YieldCurve type - "
                                 "unable to populate FDModel::discYC field");

        // The top level instrument must supply discount factor curve in order to define the
        // pricing/PnL currency of the trade.  This is retrieved from the IModel perspective
        // through a market data fetcher, and KComponent implementing the IInstrumentCollection
        // interface function...
        // virtual string discountYieldCurveName() const = 0
        // To double check internal consistency, ensure that the name of the discount curve 
        // retrieved in the above instrument object recursion is the same as that returned by 
        // the IModel::getDomesticYCName
        string discYCName = getDomesticYCName();
        if (discYCName.empty())
            throw ModelException("instrument must supply a discount Curve");
        if (discYCName != discYC->getName())
            throw ModelException("Internal consistency error - name of model discount curve "
                                 "returned from instrument iteration " + discYC->getName() + 
                                 " does not equal model discount curve derived from the internal "
                                 "IModel::getDomesticYCName() function " + discYCName);

        // if model discount curve supplied, check curve same as derived from instrument
        if (IRParams->curveToDiscount.get()) {  
            discYCName = discYC->getName();
            string modelDiscYCName = IRParams->curveToDiscount->getName();
            if (discYCName != modelDiscYCName)
                throw ModelException("Name of curveToDiscount field set in Fix3CB model " + 
                                     modelDiscYCName +
                                    " must be same as defined on the top level instrument " 
                                    + discYCName);
        }

        // recurse the instrument components to find and store all the IMarketFactor types found
        // in the exported fields - these are yield curve types for fix3
        class RetrieveFactors : public ObjectIteration::IAction {
            FDModel * model;
            const MarketData * market;
            IMarketFactorArray & factors;

        public:
            RetrieveFactors(FDModel* model, const MarketData* market, IMarketFactorArray& factors) :
                model(model), market(market), factors(factors) {factors.clear();}

            virtual bool invoke(ObjectIteration::State& state, IObjectConstSP obj) {

                IMarketFactor * factor = dynamic_cast< IMarketFactor * >(state.getObject().get());

                string name = factor->getName();
                string type = factor->getClass()->getName();
                int i;
                for (i = 0; i < factors.size(); ++i) {
                    if (name == factors[ i ]->getName() && 
                        type == factors[ i ]->getClass()->getName())
                        break;
                }
                if (i >= factors.size() && model->acceptFactor(factor)) {
                    factors.push_back(IMarketFactorSP::attachToRef(factor));
                    return false;
                }
                else
                    return true;
            }
        };
        RetrieveFactors retrieveFactors(this, market, factors);
        ObjectIteration iteration(IMarketFactor::TYPE);
        iteration.recurse(retrieveFactors, instruments);
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void Fix3CB::initTreeData(void)
{
    try {
        // ??? for now reuse the fix3 OverWriteString mechanism until model/market rewrite
        char OverWriteString[6][MAXBUFF];
        RateTree::IRModelParams mp;

        treeData.NbFactor = nbFactors;

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

        treeData.NbSigmaMax = mp.nbStdDevs;

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
        treeData.CvDiff = 0;
        treeData.CvIdx1 = 1;
        treeData.CvIdx2 = 2;
        treeData.CvDisc = 1;

        if (Fix3_Param_Input (&mktVolData, &treeData, treeData.NbFactor, 
                              OverWriteString,"") != SUCCESS)
            throw ModelException("Fix3_Param_Input falied: "+IrConverter::slog.pop());
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

int Fix3CB::getCurveIdx( const string & curveName ) const {
    return curveName.empty() ? -1 : getCrvIdx( curveName );
}

// map between the product defined curve name and the 
// internal Fix3 T_CURVE array index
int Fix3CB::getCrvIdx(const string& curveName) const {

    if (curveName.empty())
        throw ModelException(__FUNCTION__, "Curve id not supplied");
    
    int i, index=-1;
    for (i = 0; i < 3; i++) {
        if (curveName == this->curveName[i]) {
            index=i;
            break;
        }
    }
    if (index == -1)
        throw ModelException(__FUNCTION__, "Unable to find curveName " +
                             curveName +
                             " in list of 3 fix3 model registered curve names (" +
                             this->curveName[0] + ", " + this->curveName[1] + ", " +
                             this->curveName[2] + ")");

    return index;
}


// register product discounting slice requirements
void Fix3CB::registerZero(DateTime useDate, DateTime matDate, 
                          string curveName) {
    try {
        if (matDate < useDate)
            throw ModelException("matDate < useDate ! ("
            +matDate.toString()+" < "+useDate.toString()+")");

        map<string, NamedClaimBank>::iterator cb;

        // any payments are dropped on or before the value date
        if (matDate <= valueDate)
            return;

        if (matDate < useDate) {
            throw ModelException("Cannot register zero with useDate "
            +useDate.toString()+" > matDate "+matDate.toString());
        }

        cb = claimBank.find(curveName);

        // if claimBank doesn't exist, create a new one and insert into map
        if (cb == claimBank.end()) {
            NamedClaimBank newClaimBank;
            newClaimBank.curveName = curveName;
            claimBank[curveName] = newClaimBank;
            cb = claimBank.find(curveName);
        }

        DateTime useDateL = max(today, useDate);
        cb->second.critZeroUseDates.push_back(useDateL.toIrDate());
        cb->second.critZeroMatDates.push_back(matDate.toIrDate());
        cb->second.critZeroLabel.push_back(string("ZERO/DISCOUNT"));
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

// retrieve zero/discounting slice from the tree
void Fix3CB::getZero(DateTime useDate, DateTime matDate, 
                     string curveName, TreeSliceSP &slice) {
    try {
        map<string, NamedClaimBank>::iterator cb;
        DateTime currentDate = DateTime::fromIrDate(treeData.TPDate[mTpIdx]);

        if (!dynamic_cast<TreeSliceRates*>(slice.get())) {
            slice = RateTree::createSimpleSlice(curveName);
            slice->name = "Zero "+useDate.toString()+" to "+matDate.toString();
        }
        TreeSliceRates& sliceL = dynamic_cast<TreeSliceRates&>(*slice);

        if (matDate <= valueDate)
            throw ModelException("Requested zero maturity date " + matDate.toString() +
                                 " must be greater than the model valueDate " +
                                 valueDate.toString());

        cb = claimBank.find(curveName);
        if (cb == claimBank.end())
            throw ModelException("Unable to find curveName " + curveName +
                                 "in the list of model registered claimBanks");

        if (useDate != currentDate)
            throw ModelException("requested zero observation date " + useDate.toString() +
                                 " must be the same as the model's current date " +
                                 currentDate.toString());

        CLAIM_BANK const* zeroBank = NULL;
        double* slicePtr = NULL;

        zeroBank = &cb->second.zeroBank;
        sliceL.allocDim(nbFactors);
        slicePtr = sliceL.getValuePtr();

        // ??? to keep consistent with hyb3
        if (currentDate == matDate) {

            // set slice value as 1.0
            if (Fix3_Set_Slice(slicePtr,
                               1.0,
                               mTpIdx,
                               &treeData) != SUCCESS)
                throw ModelException("Fix3_SetSlice failed: " + IrConverter::slog.pop());
        }
        else {
            double* zeroBkPt = NULL;

            // pointer to internal zero slice maturity - zeroBkPt is read only and is not freed
            zeroBkPt = Fix3_ZbkReadZero((CLAIM_BANK*)zeroBank,
                                        matDate.toIrDate(),
                                        FALSE,  // ??? no interp, do we want to relax this
                                        currentDate.toIrDate(),
                                        mTpIdx,
                                        &treeData);

            if (zeroBkPt == NULL)
                throw ModelException("Fix3_ZbkReadZero failed: " + IrConverter::slog.pop());

            // copy results of ReadZero function call to output slice
            if (Fix3_Copy_Slice(slicePtr,
                                zeroBkPt,
                                mTpIdx,
                                &treeData) != SUCCESS)
                throw ModelException("Fix3_CopySlice failed: " + IrConverter::slog.pop());
        }

        slice->treeStep = range->treeStep;  // set internal slice timestep counter
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

// adjust parYield rate for any offsets where we are requesting/observing the
// rate no a date other than the reset date - uses deterministic ratio method
void Fix3CB::dateAdjustIRIndex(const IndexSpecIR& rateSpec, DateTime currentDate, 
                               DateTime resetDate, TreeSlice& treeSlice) {
    try {
        TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);

        string name = rateSpec.getName();
        char dcc = IrConverter::dccTo035A(*(rateSpec.dcc));
        char freq = IrConverter::toEslFreq(rateSpec.frequency.get());
        string curveName = rateSpec.getFactor()->getName();
        DateTime indexStartDate;
        double* sliceL = slice.getValuePtr();
        int curveIdx = getCrvIdx(curveName);
        double ratio;

        if (currentDate == resetDate)
            return;  // no adjustment required if they're the same dates
        else if (currentDate < resetDate)
            throw ModelException("Reset dates " + resetDate.toString() +
                                 "greater than the currentDate " + currentDate.toString() +
                                 "are not adjusted - product logic should cope with this case");

        if (rateSpec.fwdRateOffset.get())
            // ??? indexStartDate should be holiday adjusted
            indexStartDate = rateSpec.fwdRateOffset->toDate(resetDate);
        else
            indexStartDate = resetDate;  //  no offset

        // ??? think more about this later
        int nbFwdDays = indexStartDate.daysDiff(resetDate);
        if (nbFwdDays > maxFwdPeriodIRAdj)
            throw ModelException("Reset date adjustment algorithm only supports small forward "
                "starting rate offsets (<= " + Format::toString(maxFwdPeriodIRAdj) + 
                " days) as the deterministic parYieldRatio function does not currently account "
                "for forward starting rates.  Offset requested in days = " + 
                Format::toString(nbFwdDays));

        // if the reset date is before the currentDate when the rate is being asked for, 
        // approximate the rate by taking the ratio of the deterministic parYield rates 
        // between the actual reset date and the current date, then scale the stochastic 
        // payYield value by this ratio.  Essentially this is matching the first
        // moment.  We do not at this stage make volatility adjustments to the rate
        int nbOffsetDays = currentDate.daysDiff(resetDate);
        if (nbOffsetDays > maxIRDateOffsetAdj)
            throw ModelException("Current date " + currentDate.toString() + " to resetDate " +
                resetDate.toString() + " parYield approximate rate adjustment is only supported "
                "if the date offset is < " + Format::toString(maxIRDateOffsetAdj) + 
                " days due to issues around scaling the rate over too larger a time offset. "
                "Requested offset = " + Format::toString(nbOffsetDays));

        // returns ratio of (rate at resetDate) / (rate and current date)
        // ??? need a new function to take forward starting payYields rather than assuming 
        // the rate starts today - either spot offset or a significant forward starting offset
        if (ParYieldRatio(&ratio,
                          resetDate.toIrDate(), 
                          currentDate.toIrDate(),
                          0.0,  // no spread
                          &treeCurves[curveIdx],
                          rateSpec.tenor->toMonths(),
                          dcc,
                          freq) != SUCCESS)
            throw ModelException("ParYieldRatio failed: " + IrConverter::slog.pop());

        // scale stochastic rate by deterministic ratio of parYields
        if (Fix3_MultiplyScalar(sliceL, ratio, mTpIdx, &treeData) != SUCCESS)
            throw ModelException("Fix3_MultiplyScalar failed: " + IrConverter::slog.pop());
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__, "indexSpecName = " + rateSpec.getName());
    }
}



/** retrieving market data */
void Fix3CB::retrieveFactor() {
    try{
        if (!discYC)
            throw ModelException("FDModel::discYC field is not populated - must "
                                 "contain the top insrument's discount curve. Internal "
                                 "error in FDModel/product code");

        string currency = discYC->getCcy();

        if (!IRParams->curveToDiffuse)
            throw ModelException("Must supply curveToDiffuse to fix3 model");

        // check curves are same currency
        if (IRParams->curveToDiffuse->getCcy() != currency)
            throw ModelException("Currency of model curve to diffuse " + 
                                 IRParams->curveToDiffuse->getCcy() +
                                 " must be the same as currency of instrument " +
                                 currency);

        // for consistency, check model curve to discount if supplied is that same as
        // the instrument curve to discount
        if (!IRParams->curveToDiscount.isEmpty()) {
            if (IRParams->curveToDiscount->getName() != discYC->getName())
                throw ModelException("model curveToDiscount supplied (name = " +
                                     IRParams->curveToDiscount->getName() +
                                     ") is not same name as discount curve defined by the product " +
                                     discYC->getName());
        }
        // by convention set 0=diff curve, 1=discount curve, 2=other
        // assume vols are associated with diffusion curve
        curveName[0] = IRParams->curveToDiffuse->getName();
        IrConverter::to_T_CURVE(treeCurves[0], IRParams->curveToDiffuse.get(),
                                false);


        curveName[1] = discYC->getName();
        IrConverter::to_T_CURVE(treeCurves[1], discYC.get(), false);

        // all factors supplied must be yield curves and in same currency for fix3.
        // fix3 only supports max 3 curves, so check if this limit is exceeded
        bool thirdCurve=false;
        int i;
        for (i = 0; i < factors.size(); i++) {
            IMarketFactorConstSP &mf = factors[i];
            const YieldCurve* yc = dynamic_cast<const YieldCurve*>(mf.get());\

            if (!yc)
                throw ModelException("Fix3 only supports yield curve type market objects - "
                                     "type supplied = " +
                                     mf->getClass()->getName());
            // now check currency
            if (yc->getCcy() != currency)
                throw ModelException("Yield curve defined in instrument of currency " +
                                     yc->getCcy() +
                                     ". Fix3 is single currency model in currency - " +
                                     currency);
            // check if third curve (ie. not discount or diffusion)
            if (yc->getName() == curveName[0] ||
                yc->getName() == curveName[1])
                continue;
            else {
                // error if thirdCurve is already defined
                if (thirdCurve)
                    throw ModelException("Fix3 engine supports maximum of three yield curves "
                                         "and supplied yield curve/factor " + yc->getName() +
                                         "is a fourth.  Current 3 curves registered are " +
                                         curveName[0] + ", " + curveName[1] + ", " +
                                         curveName[2]);
                curveName[2] = yc->getName();
                IrConverter::to_T_CURVE(treeCurves[2], yc, false);
                thirdCurve = true;
            }
        }
        if (thirdCurve == false) {
            // assign 3rd curve same as second curve - rates default behaviour
            curveName[2] = curveName[1];
            treeCurves[2] = treeCurves[1];
        }

        // check that the curves have the same value date, and set this to the
        // spot date of the IR curves
        curveValueDate = DateTime::fromIrDate(treeCurves[0].ValueDate);
        for (i = 1; i < 3; i++) {
            if (treeCurves[i].ValueDate != curveValueDate.toIrDate())
                throw ModelException("ValueDate " + 
                                     DateTime::fromIrDate(treeCurves[i].ValueDate).toString() +
                                     " defined in zero curve index [" + Format::toString(i) + "] is not the same "
                                     "as the valueDate defined in the other zero curve(s) " +
                                     curveValueDate.toString());
        }

        IRVolSelector volSelector(getIRVolRaw(IRParams), treeCurves[0],
                                  IRParams->volCalibIndex,
                                  cetSmoothing, IRParams);
        volSelector.getMktVolData(mktVolData);
    }
    catch (exception& e){
        throw ModelException(e, __FUNCTION__);
    }
}

// Implementation of IRVegaPointwise::ISensitivePoints for smart vega tweaking
// and volatility exposure reporting.
IRGridPointAbsArraySP Fix3CB::getSensitiveIRVolPoints(OutputNameConstSP  outputName, 
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


void Fix3CB::configureClaimBank(const CRIT_DATE* critDate, int nbCritDate) {
    try {
        int i;

        // store zero/discount factor between today and the value date to
        // in order to forward the price back to the value date
        for (i=0; i<NBCRV; ++i) {
            today2ValDateZeros[i] = ::pow(1.0 + treeCurves[i].Zero[0], 
                                    -Daysact(today.toIrDate(), curveValueDate.toIrDate())/365.0);

            // adjust the zero dates from the value to today by scaling all of 
            // the zeroRates in the TCurves       
            RateTree::ExtendTreeCurve(&treeCurves[i], today);
        }

        // add critical claim bank zero dates for each named bank
        map<string, NamedClaimBank>::iterator cb;

        // produce zero bank dates and optimize
        for (cb = claimBank.begin(); cb != claimBank.end(); cb++) {

            int nbCrit = 0;
            long *critZeroMatDL = NULL;
            long *critZeroUseDL = NULL;
            IrConverter::AutoDeleteDRArray<long>
                critZeroMatDelete(&critZeroMatDL, LONG, &nbCrit);
            IrConverter::AutoDeleteDRArray<long>
                critZeroUseDelete(&critZeroUseDL, LONG, &nbCrit);

            // optimize non-critical dates
            int     nbMatDates = 0;
            long*   matDates   = NULL;
            int     nbUseDates = 0;
            long*   useDates   = NULL;
            int     j;

            // take copy of registered critical dates to use when printing out
            // the dates in the Fix3CB::print function, along with the description labels
            cb->second.critZeroMatDatesUnsorted = cb->second.critZeroMatDates;
            cb->second.critZeroUseDatesUnsorted = cb->second.critZeroUseDates;

            // get pointer to required array to make debugging easier
            IRDate* optMatDates = &(cb->second.optZeroMatDates[0]);
            IRDate* optUseDates = &(cb->second.optZeroUseDates[0]);

            if (ZbkOptDates(cb->second.optZeroMatDates.size(),
                            optMatDates,
                            cb->second.optZeroUseDates.size(),
                            optUseDates,
                            nbCritDate,
                            const_cast<CRIT_DATE*>(critDate),  // remove const
                            valueDate.toIrDate(),
                            &nbMatDates,
                            &matDates,
                            &nbUseDates,
                            &useDates) != SUCCESS)
                throw ModelException(IrConverter::slog.pop());
                        
            // combine optimized dates with critical dates
            for (j=0; j<nbUseDates; ++j) {
                cb->second.critZeroUseDates.push_back(useDates[j]);
                cb->second.critZeroMatDates.push_back(matDates[j]);
            }

            // ??? for now, take a copy of the dateList to reuse cbkProcessDL function
            // which reallocs the arrays, so I have to use the C array rather than 
            // C++ vector
            nbCrit = cb->second.critZeroUseDates.size();
            critZeroMatDL = (long *) DR_Array (LONG, 0, nbCrit-1);
            critZeroUseDL = (long *) DR_Array (LONG, 0, nbCrit-1);

            if (critZeroMatDL == NULL || critZeroUseDL == NULL) {
                throw ModelException(
                    "Unable to allocate memory for claim bank critical dates. "
                    +IrConverter::slog.pop());
            }
            for (i = 0; i < nbCrit; i++) {
                critZeroUseDL[i] = cb->second.critZeroUseDates[i];
                critZeroMatDL[i] = cb->second.critZeroMatDates[i];
            }

            // sort and remove duplicates from the critical date list
            if (CbkProcessDL(&nbCrit,
                            &critZeroMatDL,
                            &nbCrit,
                            &critZeroUseDL) != SUCCESS)
                throw ModelException(IrConverter::slog.pop());

            // copy result back to the claimBank vectors and resize
            cb->second.critZeroMatDates.resize(nbCrit);
            cb->second.critZeroUseDates.resize(nbCrit);

            for (i = 0; i < nbCrit; i++) {
                cb->second.critZeroMatDates[i] = critZeroMatDL[i];
                cb->second.critZeroUseDates[i] = critZeroUseDL[i];
            }

            // set the critical dates index for the last date
            cb->second.critDatesIdx = MAX(0, nbCrit - 1);
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


/**  collect model initialisation data, set up timeline  */
void Fix3CB::initModel(void) {
    try {
        treeBuilt = true;

        EslSetZeroInterpolation(zeroInterpStyle == ZeroInterpStyle::LINEAR ?
            ESL_INTERP_LINEAR : ESL_INTERP_FLATFWD);

        int i, j;
        CRIT_DATE* critDate=NULL;
        int nbCritDate = 0;
        IrConverter::AutoDeleteDRArray<CRIT_DATE>
            critDateDelete(&critDate, CRITDATE, &nbCritDate);

        map<string, NamedClaimBank>::iterator cb;

        initTreeData();

        // needed to setup fix3 timeline structures - have to register 
        // critical dates that are >= today
        critDate = (CRIT_DATE *) DR_Array (CRITDATE, 0, 0);
        if (critDate == NULL)
            throw ModelException(IrConverter::slog.pop());
        
        // add tree startDate as critical date
        if (Add_To_DateList(&nbCritDate, 
                            &critDate,
                            today.toIrDate(),
                            0, 0, 0, 0, 0, 0, 0, 0, 0) != SUCCESS)
            throw ModelException(IrConverter::slog.pop());

        // add valueDate as critical date
        if (Add_To_DateList(&nbCritDate, 
                            &critDate,
                            valueDate.toIrDate(),
                            0, 0, 0, 0, 0, 0, 0, 0, 0) != SUCCESS)
            throw ModelException(IrConverter::slog.pop());

        // add product specific critical dates
        for (j = 0; j < critDates.size(); j++) {
            DateTime critDateL = critDates[j];
            if (critDateL > today) {
                if (Add_To_DateList(&nbCritDate, 
                                    &critDate,
                                    critDateL.toIrDate(),
                                    0, 0, 0, 0, 0, 0, 0, 0, 0) != SUCCESS)
                    throw ModelException(IrConverter::slog.pop());
            }
        }

        // setup claim banks and add these critical dates to list
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
                    throw ModelException(IrConverter::slog.pop());
                if (Add_To_DateList(&nbCritDate, 
                                    &critDate,
                                    useDate,
                                    0,0, 0, 0, 
                                    0, 0, 0, 0, 0) != SUCCESS)
                    throw ModelException(IrConverter::slog.pop());
            }
        }

        // remove duplicates, sort critical dates and store them to
        // help debugging by printing them out in the debug data file
        if (Sort_CritDate (nbCritDate, critDate) != SUCCESS)
            throw ModelException(IrConverter::slog.pop());

        sortedCritDates.resize(1);
        sortedCritDates[0] = critDate[0].CritDate;
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
        // construct the tree timeline - tree starts from today, but
        // assumes product values itself to valueDate
        if (Fix3_Time_Line(today.toIrDate(),
                        nbCritDate,
                        critDate,
                        'I',                            
                        &treeData) != SUCCESS) {
            throw ModelException(IrConverter::slog.pop());
        }

        // numerical volatility calibration routine
        if (Fix3_Cet_Main(FALSE, treeCurves, &mktVolData, &treeData) != SUCCESS) {
            throw ModelException(IrConverter::slog.pop());
        }

        if (Fix3_Build_Tree(treeCurves, &mktVolData, &treeData) != SUCCESS) {
            throw ModelException(IrConverter::slog.pop());
        }

        if (treeData.NbFactor == 1) {
            range.reset(new TreeSliceRates::Range(
                -treeData.HalfWidth[0], 
                treeData.HalfWidth[0]));
        }
        else if (treeData.NbFactor == 2) {
            range.reset(new TreeSliceRates::Range(
                -treeData.HalfWidth[0],
                treeData.HalfWidth[0],
                -treeData.HalfWidth[1], 
                treeData.HalfWidth[1]));
        }
        else if (treeData.NbFactor == 3) {
            range.reset(new TreeSliceRates::Range(
                -treeData.HalfWidth[0],
                treeData.HalfWidth[0],
                -treeData.HalfWidth[1], 
                treeData.HalfWidth[1],
                -treeData.HalfWidth[2], 
                treeData.HalfWidth[2]));
        }

        // allocate claim bank internal memory
        for (cb = claimBank.begin(); cb != claimBank.end(); cb++) {
            Fix3_CbkInit(&cb->second.zeroBank);
            int nbCrit = cb->second.critZeroMatDates.size();

            if (nbCrit > 0) {
                if (Fix3_CbkAlloc(&cb->second.zeroBank, nbCrit, &treeData) != SUCCESS) {
                    throw ModelException(IrConverter::slog.pop()
                        +" Failed to allocate claim bank "+cb->second.curveName);
                }
            }
        }

        if (Fix3_Dev_Alloc(&devData, &treeData) != SUCCESS) {
            throw ModelException(IrConverter::slog.pop());
        }

        mTpIdx = treeData.NbTP;  // set internal step counter to final node number

        // print debug tree file if requested
        print();
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


DateTime Fix3CB::getToday() const {

    return today;
}

//------------------------------------------------------------------------------
// Get date on the timeline pointed by the index
//------------------------------------------------------------------------------
DateTime Fix3CB::getDate(int dateIdx) const
{
    if (treeData.NbTP == -1)  // degenerate tree with no critical date (=0 price)
        return DateTime();

    if (dateIdx < 0 || dateIdx > treeData.NbTP) {    
        throw ModelException(__FUNCTION__, 
        "Invalid time point index (" + Format::toString(dateIdx) 
        + "). Must be in the range 0..." + Format::toString(treeData.NbTP));
    }
    return DateTime::fromIrDate(treeData.TPDate[dateIdx]);
}

DateTimeArray Fix3CB::getDates() const
{
    DateTimeArray dates( treeData.NbTP + 1 );
    for( int i = 0; i <= treeData.NbTP; ++i )
        dates[ i ] = DateTime::fromIrDate( treeData.TPDate[ i ] );
    return dates;
}

/** get last step (total num of steps) on time line */
int Fix3CB::getLastStep() const
{
    return treeData.NbTP;
}

DateTime Fix3CB::getCurveValueDate(string curveName) const {
    try {
        getCrvIdx(curveName);  // just to check if valid curveName
        return valueDate;
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


void Fix3CB::print() {

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
    fprintf(stream, "today:     %s\n", today.toString().c_str());
    fprintf(stream, "valueDate: %s\n", valueDate.toString().c_str());
    fprintf(stream, "curveValueDate: %s\n", curveValueDate.toString().c_str());

    fprintf(stream, "\nFollowing critical dates registered with the engine:\n\n");
    for (i = 0; i < (int)sortedCritDates.size(); i++) {
        fprintf(stream, "%3d   %ld\n", i, sortedCritDates[i]);
    }

    fprintf(stream, "\nFollowing claim banks registered with the engine:\n");

    for (cb = claimBank.begin(); cb != claimBank.end(); cb++) {
        fprintf(stream, "\n%d, name: %s, curveName: %s, critDatesIdx: %d\n",
                count,
                cb->first.c_str(),
                cb->second.curveName.c_str(),
                cb->second.critDatesIdx);

        fprintf(stream, "\n\tCritical Zero Dates used by tree\n");
        fprintf(stream, "\t[idx]  (zero Use/Obs Date)   (zero Maturity Date)\n");

        for (i = 0; i < (int)cb->second.critZeroMatDates.size(); i++) {
            fprintf(stream, "\t%3d   %ld     %ld\n",
                    i,
                    cb->second.critZeroUseDates[i],
                    cb->second.critZeroMatDates[i]);
        }

        fprintf(stream, "\n\tOptional Zero Dates\n");
        fprintf(stream, "\t[idx]  (zero Use/Obs Date)   (zero Maturity Date)\n");

        for (i = 0; i < (int)cb->second.optZeroMatDates.size(); i++) {
            fprintf(stream, "\t%3d   %ld     %ld\n",
                    i,
                    cb->second.optZeroUseDates[i],
                    cb->second.optZeroMatDates[i]);
        }
        fprintf(stream, "\n\tUnsorted Critical Zero Dates (as registered by products "
                        "before internal sorting and duplicate removal)\n");
        fprintf(stream, "\t[idx]  (zero Use/Obs Date)   (zero Maturity Date)    (label)\n");

        for (i = 0; i < (int)cb->second.critZeroMatDatesUnsorted.size(); i++) {
            fprintf(stream, "\t%3d   %ld     %ld    %s\n",
                    i,
                    cb->second.critZeroUseDatesUnsorted[i],
                    cb->second.critZeroMatDatesUnsorted[i],
                    cb->second.critZeroLabel[i].c_str());
        }
    }

    fprintf (stream, "\n");

    fprintf(stream, "Domestic Currency: %s\n\n", currency.c_str());

    Fix3_Print_FIX3_TREE_DATA(stream, &treeData);
    Print_MKTVOL_DATA(stream, &mktVolData);

    fclose(stream);    
}


void Fix3CB::clear() {

    if (!treeBuilt)
        return;
    treeBuilt = false;

    Fix3_Dev_Free(&devData, &treeData);

    map<string, NamedClaimBank>::iterator cb;
    for (cb = claimBank.begin(); cb != claimBank.end(); cb++) {
        Fix3_CbkFree(&cb->second.zeroBank, &treeData);
    }

    Fix3_Tree_Free(&treeData);

    for (size_t i=0; i<NBCRV; ++i)
        today2ValDateZeros[i] = -1.0;

    mTpIdx = -1;
}


void Fix3CB::update(int t, FDProduct::UpdateType type) {
    try {
        map<string, NamedClaimBank>::iterator cb;

        if (!treeBuilt)
            throw ModelException("Not able to update tree as the "
                                "tree has not been built");

        // set current range
        range->limits.bot1 = treeData.Bottom1[t];
        range->limits.top1 = treeData.Top1[t];
        range->limits.bot2 = treeData.Bottom2[t];
        range->limits.top2 = treeData.Top2[t];
        range->limits.bot3 = treeData.Bottom3[t];
        range->limits.top3 = treeData.Top3[t];
        range->treeStep = t;

        int T = treeData.NbTP;

        // update tree and probabilities
        if (Fix3_Lattice(&devData,
                        t,
                        T,
                        &mktVolData,
                        &treeData) != SUCCESS)
            throw ModelException(IrConverter::slog.pop());
        
        mTpIdx = t;

        // update the claim banks
        for (cb = claimBank.begin(); cb != claimBank.end(); cb++) {

            int  addFlag = 0;
            long useDate = 0;
            long matDate = 0;
            int  idxActive = cb->second.critDatesIdx;
            string curveName = cb->second.curveName;
            long currentDate = getDate(mTpIdx).toIrDate();

            if (cb->second.critZeroMatDates.size() == 0)
                continue;

            if (idxActive >= 0 &&
                cb->second.critZeroMatDates[idxActive] >= currentDate) {

                addFlag = 1;
                useDate = cb->second.critZeroUseDates[idxActive];
                matDate = cb->second.critZeroMatDates[idxActive];

                cb->second.critDatesIdx--;
            }

            if (Fix3_ZbkUpdate(&cb->second.zeroBank,
                            addFlag,
                            currentDate,
                            useDate,
                            mTpIdx,
                            treeData.NbTP,
                            getCrvIdx(curveName),
                            &devData,
                            &treeData) != SUCCESS) {
                throw ModelException(
                    "zero bank update failed for curve "+curveName
                    +". "+IrConverter::slog.pop());
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


// ?? nothing to do in claim bank mode - remove when zerobankMode depricated
void Fix3CB::updateFinalize(int t, FDProduct::UpdateType type) {

    if (!treeBuilt)
        throw ModelException(__FUNCTION__, "Not able to update tree as "
                             "the tree has not been built");
}

//------------------------------------------------------------------------------
void Fix3CB::sliceExpand(TreeSlice& treeSlice, const string& curveName) const {
    try {
        TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);
        slice.expand(nbFactors);
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}
//------------------------------------------------------------------------------

void Fix3CB::sliceDev(TreeSlice& treeSlice, int curveIdx) const {
    try {
        TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);

        // nothing to do at the back of the tree 
        if (mTpIdx == treeData.NbTP)
            return;

        if (slice.treeStep != mTpIdx+1)
            throw ModelException("Unable to DEV Slice (name = " + slice.name +
                                ") from tree step " + Format::toString(mTpIdx+1) + 
                                ", date " + getDate(mTpIdx+1).toString() + ", to tree step " +
                                Format::toString(mTpIdx) + ", date " + getDate(mTpIdx).toString()  +
                                " as internal slice time step counter = " +
                                Format::toString(slice.treeStep));

        if (Fix3_Dev (slice.getValuePtr(), mTpIdx, treeData.NbTP, curveIdx, &devData, &treeData)
            != SUCCESS) {
                throw ModelException(IrConverter::slog.pop());
        }

        slice.treeStep = mTpIdx;
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__ + string(", tree slice name = ") + treeSlice.name);
    }
}

// record a lot of additional debug information - for now always store
// ??? need to turn this on/off from a control setting so have to specify what
// that control setting is 
void Fix3CB::recordOutput(Control* ctrl, Results* results) const {
    try {
        int i, j;

        if (!ctrl->requestsOutput(OutputRequest::DBG))
            return; // no debug info requested

        // record today and valueDates
        DateTimeSP todaySP(new DateTime(today));
        results->storeGreek(todaySP,
                            Results::DEBUG_PACKET,
                            OutputNameConstSP(new OutputName("TODAY")));

        DateTimeSP valueDateSP(new DateTime(valueDate));
        results->storeGreek(valueDateSP,
                            Results::DEBUG_PACKET,
                            OutputNameConstSP(new OutputName("VALUE_DATE")));

        DateTimeSP curveValueDateSP(new DateTime(curveValueDate));
        results->storeGreek(curveValueDateSP,
                            Results::DEBUG_PACKET,
                            OutputNameConstSP(new OutputName("IR_SPOT_DATE")));

        // store tree critical dates
        DateTimeArraySP treeDates(new DateTimeArray(treeData.NbTP+1));
        for (i = 0; i <= treeData.NbTP; i++)
            (*treeDates)[i] = DateTime::fromIrDate(treeData.TPDate[i]);

        results->storeGreek(treeDates,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("TREE_DATES")));

        // store number of days offsets from start of tree for each critical date
        IntArraySP daysOffsets(new IntArray(treeData.NbTP+1));
        for (i = 0; i <= treeData.NbTP; i++) {

            double days = Daysact(treeData.TPDate[0], treeData.TPDate[i]);
            (*daysOffsets)[i] = (int)days;
        }

        results->storeGreek(daysOffsets,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("TREE_DATE_OFFSETS")));

        /*
        // store tree drift - drift of centre node
        DoubleArraySP treeDrift(new DoubleArray(treeData.NbTP+1));
        for (i = 0; i <= treeData.NbTP; i++)
            (*treeDrift)[i] = treeData.ZCenter[i];

        results->storeGreek(treeDrift,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("TREE_DRIFT")));  */

        // store zero rates
        CDoubleMatrixSP domZeros(new CDoubleMatrix(3, treeData.NbTP+1));
        for (i = 0; i < 3; i++)
            for (j = 0; j <= treeData.NbTP; j++)
                (*domZeros)[i][j] = treeData.ZeroRate[i][j];

        results->storeGreek(domZeros,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("DOMESTIC_IR_ZEROS")));

        /*
        DoubleArraySP oneFacIRSpotVol(new DoubleArray(treeData.NbTP+1));
        for (i = 0; i <= treeData.NbTP; i++)
            (*oneFacIRSpotVol)[i] = treeData.Aweight[0][i];

        results->storeGreek(oneFacIRSpotVol,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("1FACTOR_IR_SPOT_VOLS"))); */

        if (nbFactors > 1) {
            // ??? record spot vols/weights for 2nd factor
        }
        if (nbFactors > 2) {
            // ??? record spot vols/weighs for 3rd factor
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


// Perform one backward induction step by moving slice from t+1 to t - no discounting
void Fix3CB::sliceEv(TreeSlice& treeSlice) const {

    TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);
    // nothing to do at the back of the tree 
    if (mTpIdx == treeData.NbTP)
        return;   

    if (slice.getDim() != nbFactors)
        slice.allocDim(nbFactors);

    if (Fix3_Ev(slice.getValuePtr(), mTpIdx, treeData.NbTP, &devData, &treeData)
        != SUCCESS)
        throw ModelException(__FUNCTION__, IrConverter::slog.pop());
}


void Fix3CB::getIRIndex(const IndexSpecIR& rateSpec, DateTime currentDate, DateTime resetDate, 
                        TreeSlice& treeSlice) {
    try {
        TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);

        if (slice.getDim() != nbFactors)
            slice.allocDim(nbFactors);

        string name = rateSpec.getName();
        char dcc = IrConverter::dccTo035A(*rateSpec.dcc);
        char freq = IrConverter::toEslFreq(rateSpec.frequency.get());

        string curveName = rateSpec.getFactor()->getName();
        DateTime indexStartDate;

        if (rateSpec.fwdRateOffset.get())
            // ??? indexStartDate should be holiday adjusted
            indexStartDate = rateSpec.fwdRateOffset->toDate(currentDate);
        else
            indexStartDate = currentDate;  //  no offset supplied

        if (indexStartDate < currentDate)
            throw ModelException("IR Index swap start date " + 
                                 indexStartDate.toString() +
                                 " can not be before the reset date " + currentDate.toString() +
                                 ".  Index name = " + name);

        if (resetDate < today)
            throw ModelException("Tree can not calculate past rate index values - resetDate = "
                                 + resetDate.toString() + " is before today (= " + today.toString() +
                                 ").  (NOTE: time of day may be the issue).  An internal  logic error " +
                                 "exists as any past fixing requests should have been intercepted in the "
                                 "index FDProduct code and this function should not have been called");

        if (currentDate < resetDate)
            throw ModelException("Currently future reset dates (" + resetDate.toString() + 
                                 ") greater than the currentDate (" + currentDate.toString() + 
                                 ") are not adjusted - the product logic should be able to observe "
                                 "and use the actual reset value observed on the reset date ");

        map<string, NamedClaimBank>::iterator cb = claimBank.find(curveName);
        if (cb == claimBank.end())
            throw ModelException("Unable to find named claim bank " + curveName +
                                 " to calculate IR index payYield " + name);

        CLAIM_BANK const* zeroBank = &cb->second.zeroBank;

        // annuity slice is not used, but required and populated for function call
        // so construct local slice
        TreeSliceRates annuitySlice(*range, "", -1);  //
        annuitySlice.allocDim(nbFactors);
        double* sliceL = slice.getValuePtr();
        double* annuityL = annuitySlice.getValuePtr();

        if (Fix3_ZbkParYield_t(sliceL,
                               annuityL,
                               &cb->second.zeroBank,
                               currentDate.toIrDate(),
                               indexStartDate.toIrDate(),
                               rateSpec.tenor->toMonths(),
                               dcc,
                               freq,
                               0.0,  // spread
                               mTpIdx,
                               &treeData) != SUCCESS)
            throw ModelException("Fix3_ZbkParYield_t failed: " + 
                                 IrConverter::slog.pop());

        slice.treeStep = range->treeStep;

        // adjust the stochastic rate if resetDate != currentDate - ie. estimating
        // the rate on a date that is not the actual reset date
        dateAdjustIRIndex(rateSpec, currentDate, resetDate, treeSlice);

    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


// Here we always return a slice of value 1.0 in order to allow product
// logic to be generalized
void Fix3CB::getFXIndex(TreeSlice& treeSlice) {
    try {
        TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);

        if (slice.getDim() != nbFactors)
            slice.allocDim(nbFactors);

        slice = 1.0;
        slice.treeStep = range->treeStep;
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


/** Create a MarketDataFetcher which will be used for retrieving market data etc */
MarketDataFetcherSP Fix3CB::createMDF() const {
    MarketDataFetcherSP mdf(new RateTreeMDF(nbFactors));
    MDFUtil::setUseCurrencyBasis(*mdf, true);
    return mdf;
}

/** forward the price back to the valueDate as the convention is to report values to
    the valueDate, NOT today */
double Fix3CB::getPrice0(const TreeSlice& price) const {
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
        // vdPrice /= today2ValDateZeros[1];

        // ??? new method:
        // as today=valueDate in the current implementation, no need to forward value
        // from today to valueDate

        return vdPrice;
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

bool Fix3CBLoad(void) {
    return (Fix3CB::TYPE != 0);
}

DRLIB_END_NAMESPACE
