//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : SingleCcyTree.cpp
//
//   Description : temporary version of fix3 to work on the claim/bank and
//                 tree starting today without cluttering up the existing Fix3
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SingleCcyTree.hpp"
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


/************************************** SingleCcyTree *************************************/

/** Create a MarketDataFetcher which will be used for retrieving market data etc */
MarketDataFetcherSP SingleCcyTree::createMDF() const {
    MarketDataFetcherSP mdf(new RateTreeMDF(nbFactors));
    MDFUtil::setUseCurrencyBasis(*mdf, true);
    return mdf;
}

void SingleCcyTree::getMarket(const MarketData*  market, IInstrumentCollectionSP instruments) {
    try {
        // ??? these may be modified later on to curveValueDate if model in legacyMode
        today = market->GetReferenceDate();
        valueDate = today;

        {// check time validity
            int time = today.getTime();

            // tree only supports today as SOD or EOD - no other times
            if (time != DateTime::START_OF_DAY_TIME &&
                time != DateTime::END_OF_DAY_TIME) 
                throw ModelException("Model only supports today date times of START_OF_DAY_TIME "
                                    "and END_OF_DAY_TIME.  DateTime time received = " +
                                    Format::toString(time));
        }
        IRParams->getData(this, market);  // required for market wrapper

        // Get Engine
        if (!engineTable.isEmpty())
            engineTable.getData(this, market);

        discYC = collectYC(getDomesticYCName(), instruments);

        {// check discount curve consistency

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
                string modelDiscYCName = IRParams->curveToDiscount->getName();
                if (discYCName != modelDiscYCName)
                    throw ModelException("Name of curveToDiscount field set in SingleCcyTree model " + 
                                        modelDiscYCName +
                                        " must be same as defined on the top level instrument " 
                                        + discYCName);
            }
        }

        // recurse the instrument components to find and store all the IMarketFactor types found
        // in the exported fields - these are yield curve types for fix3
        collectFactors(this, instruments, factors);
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void SingleCcyTree::populateTreeCurves() {
    try {
        int i;
        string currency = discYC->getCcy();
        IRParams->validate(currency, discYC->getName());

        // by convention set 0=diff curve, 1=discount curve, 2=other
        // assume vols are associated with diffusion curve
        curveNames[0] = IRParams->curveToDiffuse->getName();
        IrConverter::to_T_CURVE(treeCurves[0], IRParams->curveToDiffuse.get(), false);

        curveNames[1] = discYC->getName();
        IrConverter::to_T_CURVE(treeCurves[1], discYC.get(), false);

        // all factors supplied must be yield curves and in same currency for fix3.
        // fix3 only supports max 3 curves, so check if this limit is exceeded
        bool thirdCurve=false;
        for (i = 0; i < factors.size(); i++) {
            const YieldCurve* yc = dynamic_cast<const YieldCurve*>(factors[i].get());

            if (!yc) {
                throw ModelException(
                    getClass()->getName()+" only supports factors of type YieldCurve, not " 
                    + factors[i]->getClass()->getName());
            }
            // now check currency
            if (yc->getCcy() != currency) {
                throw ModelException(
                    "Instrument references a YieldCurve with currency " + yc->getCcy() +
                    ". "+getClass()->getName()+" is single currency model and the currency "
                    "of IRParams->curveToDiffuse is different (" + currency + ")");
            }
            // check if third curve (ie. not discount or diffusion)
            if (yc->getName() == curveNames[0] ||
                yc->getName() == curveNames[1])
            {
                continue;
            }
            else {
                // error if thirdCurve is already defined
                if (thirdCurve) {
                    throw ModelException(getClass()->getName()+" engine supports maximum "
                                         "of three yield curves "
                                         "and supplied yield curve/factor " + yc->getName() +
                                         "is a fourth.  Current 3 curves registered are " +
                                         curveNames[0] + ", " + curveNames[1] + ", " +
                                         curveNames[2]);
                }
                thirdCurve = true;
                curveNames[2] = yc->getName();
                IrConverter::to_T_CURVE(treeCurves[2], yc, false);
            }
        }
        if (!thirdCurve) {
            // assign 3rd curve same as second curve - rates default behaviour
            curveNames[2] = curveNames[1];
            treeCurves[2] = treeCurves[1];
        }

        // set the interpolation type on each curve, as esl uses this now instead of
        // the global interp type flag;
        // ??? add this to irconverter::to_T_CURVE
        for (i = 0; i < NBCRV; i++)
            treeCurves[i].InterpType = zeroInterpStyle;

        // check that the curves have the same value date, and set this to the
        // spot date of the IR curves
        curveValueDate = DateTime::fromIrDate(treeCurves[0].ValueDate);
        for (i = 1; i < NBCRV; i++) {
            if (treeCurves[i].ValueDate != curveValueDate.toIrDate()) {
                throw ModelException("Curve valueDate " + 
                    DateTime::fromIrDate(treeCurves[i].ValueDate).toString() +
                    " defined in zero curve index [" + Format::toString(i) + "] is not the same "
                    "as the valueDate defined in the other zero curve(s) " +
                    curveValueDate.toString());
            }
        }

        // ??? temporary switch to enable legacy mode where valueDate = curveValueDate and
        // today = curveValueDate
        if (legacyMode) {
            // maintain time of today when changing it to valueDate
            DateTime tmpDate(curveValueDate.getDate(), today.getTime());
            today = tmpDate;
            valueDate = tmpDate;
        }
    }
    catch (exception& e){
        throw ModelException(e, __FUNCTION__);
    }
}

void SingleCcyTree::initModel(void) {
    try {
        treeBuilt = true;

        if (Fix3_Model_Interface_Init(mktVolData.ModelChoice) != SUCCESS)
            throw ModelException(IrConverter::slog.pop());

        EslSetZeroInterpolation(zeroInterpStyle == ZeroInterpStyle::LINEAR ?
            ESL_INTERP_LINEAR : ESL_INTERP_FLATFWD);

        initTreeData();

        { // [begin] critical dates
            RatesUtils::LocalArray<CRIT_DATE, CRITDATE> critDatesL;

            RatesUtils::AddDateToList(critDatesL, today);
            RatesUtils::AddDateToList(critDatesL, valueDate);

            // add product specific critical dates
            for (int j = 0; j < critDates.size(); j++) {
                DateTime d = critDates[j];
                if (today < d) {
                    RatesUtils::AddDateToList(critDatesL, d);
                }
            }

            // setup claim banks and add these critical dates to list
            configureClaimBank(critDatesL.array, critDatesL.nb);

            // add claim bank critical dates to list
            map<string, NamedClaimBank>::iterator cbIter;
            for (cbIter = claimBank.begin(); cbIter != claimBank.end(); cbIter++) {
                NamedClaimBank &cb = cbIter->second;

                // add combined critical zero bank dates
                for (size_t j=0; j < cb.critZeroUseDates.size(); ++j) {
                    RatesUtils::AddDateToList(critDatesL, DateTime::fromIrDate(cb.critZeroMatDates[j]));
                    RatesUtils::AddDateToList(critDatesL, DateTime::fromIrDate(cb.critZeroUseDates[j]));
                }
            }

            // remove duplicates, sort critical dates and store them to
            // help debugging by printing them out in the debug data file
            if (Sort_CritDate (critDatesL.nb, critDatesL.array) != SUCCESS)
                throw ModelException(IrConverter::slog.pop());

            sortedCritDates.resize(1);
            sortedCritDates[0] = critDatesL.array[0].CritDate;
            for (int i = 1; i < critDatesL.nb; i++) {
                if (sortedCritDates.back() != critDatesL.array[i].CritDate) {
                    sortedCritDates.push_back(critDatesL.array[i].CritDate);
                }
            }
            if (sortedCritDates.size()<=1) {
                range.reset(new TreeSliceRates::Range(0,0,0,0,0,0));
                treeData.NbTP = -1;
                return;
            }
            // construct the tree timeline - tree starts from today, but
            // assumes product values itself to valueDate
            if (Fix3_Time_Line_Nmr(today.toIrDate(),
                            &critDatesL.nb,
                            &critDatesL.array,
                            'I',
                            &mktVolData,                            
                            &treeData) != SUCCESS)
            {
                throw ModelException(IrConverter::slog.pop());
            }
        } // [end] critical dates

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
        map<string, NamedClaimBank>::iterator cbIter;
        for (cbIter = claimBank.begin(); cbIter != claimBank.end(); cbIter++) {

            NamedClaimBank &cb = cbIter->second;
            int nbCrit = cb.critZeroMatDates.size();
            if (nbCrit > 0) {
                if (Fix3_CbkAlloc(&cb.zeroBank, nbCrit, &treeData) != SUCCESS) {
                    throw ModelException(IrConverter::slog.pop()
                        +" Failed to allocate claim bank "+cb.curveName);
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

// used for single currency trees that support standard 2Q and non-term structure
// VNFM paramters
void SingleCcyTree::initTreeData(void) {
    try {

        // ??? for now reuse the fix3 OverWriteString mechanism until model/market rewrite
        char OverWriteString[6][MAXBUFF];
        RateTree::IRModelParams mp;

        treeData.NbFactor = nbFactors;

        // retrieve model parameters from market cache - must be marketTable style
        if (!IRParams->irCalib.isEmpty())
            throw ModelException("Model does not support legacy IRCalib structure "
                "(field name irCalib) to define model and smile parameters - must "
                "use new market structures");

        IRExoticParamTable* pSmileTable = IRParams->smileTable.get();
        IRExoticParamTable* pModelTable = IRParams->modelTable.get();

        if (engineTable.isEmpty())
            throw ModelException("must supply engineTable when supplying other tables");
        if (!engineTable.get())
            throw ModelException("engineTable market wrapper pointer is NULL");

        // Use new method using MarketTable objects
        IrConverter::to_MODELPARAMETERS_DATA(mp, 
                                            IRParams->smileSet, 
                                            IRParams->modelSet,
                                            engineSet,
                                            *(engineTable.get()),
                                            *pSmileTable,
                                            *pModelTable);
        if (nbFactors != mp.nbFactors) {
            throw ModelException("nbFactors defined in tree model object (" +
                Format::toString(nbFactors) + ") does not equal number of factors "
                "defined in named modelParameters set " + IRParams->modelSet + 
                ", which defines " + Format::toString(mp.nbFactors) + " factors");
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

void SingleCcyTree::configureClaimBank(const CRIT_DATE* critDate, int nbCritDate) {
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
        // produce zero bank dates and optimize
        map<string, NamedClaimBank>::iterator cbIter;
        for (cbIter = claimBank.begin(); cbIter != claimBank.end(); cbIter++) {

            NamedClaimBank &cb = cbIter->second;

            cb.saveUnsorted();

            RatesUtils::zbkOptDates(
                cb.optZeroMatDates,
                cb.optZeroUseDates,
                nbCritDate,
                critDate,
                valueDate,
                cb.critZeroMatDates,
                cb.critZeroUseDates);

            RatesUtils::cbkProcessDL(cb.critZeroMatDates,
                                     cb.critZeroUseDates);

            // set the critical dates index for the last date
            cb.critDatesIdx = MAX(0, cb.critZeroMatDates.size() - 1);
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void SingleCcyTree::clear() {

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

void SingleCcyTree::registerZero(DateTime useDate, DateTime matDate, 
                          string curveName) 
{
    try {
        if (matDate < useDate) {
            throw ModelException("Cannot register zero with matDate "
                + matDate.toString() + " < useDate " +useDate.toString());
        }

        // any payments are dropped on or before the value date
        if (matDate <= valueDate)
            return;

        // find the corresponding claimbank
        map<string, NamedClaimBank>::iterator cb;
        cb = claimBank.find(curveName);

        // if claimBank doesn't exist, create a new one and insert into map
        if (cb == claimBank.end()) {
            NamedClaimBank newClaimBank;
            newClaimBank.curveName = curveName;
            cb = claimBank.insert(
                pair<string, NamedClaimBank>(curveName, newClaimBank)).first;
        }
  
        // add the dates to the claim bank
        DateTime useDateL = max(today, useDate);
        cb->second.critZeroUseDates.push_back(useDateL.toIrDate());
        cb->second.critZeroMatDates.push_back(matDate.toIrDate());
        cb->second.critZeroLabel.push_back("ZERO/DISCOUNT");
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

/** retrieving market data */
// Implementation of IRVegaPointwise::ISensitivePoints for smart vega tweaking
// and volatility exposure reporting.
IRGridPointAbsArraySP SingleCcyTree::getSensitiveIRVolPoints(OutputNameConstSP  outputName, 
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


//------------------------------------------------------------------------------
DateTime SingleCcyTree::getToday() const {
    
    return today;
}

DateTime SingleCcyTree::getDate(int dateIdx) const {

    if (dateIdx < 0 || dateIdx > treeData.NbTP) {
        if (treeData.NbTP <= 0)  // degenerate tree with no critical date
            return DateTime();
        throw ModelException(__FUNCTION__, 
            "Invalid time point index (" + Format::toString(dateIdx) 
            + "). Must be in the range 0..." + Format::toString(treeData.NbTP));
    }
    return DateTime::fromIrDate(treeData.TPDate[dateIdx]);
}


int SingleCcyTree::getLastStep() const {

    return treeData.NbTP;
}


DateTime SingleCcyTree::getCurveValueDate(string curveName) const {
    try {
        getCurveIdxEx(curveName);  // check if valid curveName
        return valueDate;
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

int SingleCcyTree::getCurveIdx( const string & curveName ) const {
    if (curveName.empty()) {
        return -1;
    }
    return getCurveIdxEx( curveName );
}

// map between the product defined curve name and the 
// internal Fix3 T_CURVE array index
int SingleCcyTree::getCurveIdxEx(const string& curveName) const {
    try {
        if (curveName.empty())
            throw ModelException("Curve id not supplied");
        
        for (int i = 0; i < NBCRV; ++i) {
            if (curveName == curveNames[i]) {
                return i;
            }
        }
        // perpare error message
        string curvesAvailable;
        for (int i = 0; i < NBCRV; ++i) {
            if (curveNames[i].size()) {
                if (curvesAvailable.size()) {
                    curvesAvailable += ", ";
                }
                curvesAvailable += curveNames[i];
            }
        }
        throw ModelException(
            "Unable to find curveName " 
            + curveName + ". Registered curve names are: " 
            + curvesAvailable + ".");
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


// register product discounting slice requirements

void SingleCcyTree::print() {

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
    // ??? fprintf(stream, "Model type = %s\n", *this->TYPE->getName().c_str());
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

    fprintf(stream, "Domestic Currency: %s\n\n", discYC->getCcy().c_str());

    Fix3_Print_FIX3_TREE_DATA(stream, &treeData);
    Print_MKTVOL_DATA(stream, &mktVolData);

    fclose(stream);    
}


//------------------------------------------------------------------------------
void SingleCcyTree::update(int t, FDProduct::UpdateType type) {
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

        // numeraire calculation if required
        if (Fix3_NmrToCcy(&devData,
                          t,
                          T,
                          &mktVolData,
                          &treeData) != SUCCESS)
            throw ModelException(IrConverter::slog.pop());

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
                            getCurveIdxEx(curveName),
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

void SingleCcyTree::updateFinalize(int t, FDProduct::UpdateType type) {
    try {
        if (!treeBuilt)
            throw ModelException(__FUNCTION__, "Not able to update tree as "
                                "the tree has not been built");

        int T = treeData.NbTP;

        if (Fix3_CcyToNmr(&devData,
                          t,
                          T,
                          &mktVolData,
                          &treeData) != SUCCESS)
            throw ModelException(IrConverter::slog.pop());
    } 
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void SingleCcyTree::sliceExpand(TreeSlice& treeSlice, const string& curveName) const {
    try {
        TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);
        slice.expand(nbFactors);
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}
void SingleCcyTree::sliceDev(TreeSlice& treeSlice, int curveIdx) const {
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
void SingleCcyTree::sliceEv(TreeSlice& treeSlice) const {

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


// Perform one backward induction step by moving slice from t+1 to t - no discounting
void SingleCcyTree::getIRIndex(const IndexSpecIR& rateSpec, DateTime currentDate, DateTime resetDate, 
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

// adjust parYield rate for any offsets where we are requesting/observing the
// rate no a date other than the reset date - uses deterministic ratio method
void SingleCcyTree::dateAdjustIRIndex(const IndexSpecIR& rateSpec, DateTime currentDate, 
                               DateTime resetDate, TreeSlice& treeSlice) {
    try {
        TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);
        string   name = rateSpec.getName();
        char     dcc = IrConverter::dccTo035A(*(rateSpec.dcc));
        char     freq = IrConverter::toEslFreq(rateSpec.frequency.get());
        string   curveName = rateSpec.getFactor()->getName();
        DateTime indexStartDate;
        double*  sliceL = slice.getValuePtr();
        int      curveIdx = getCurveIdxEx(curveName);
        double   ratio;

        if (currentDate == resetDate) {
            return;  // no adjustment required if they're the same dates
        }
        if (currentDate < resetDate) {
            throw ModelException("Reset dates " + resetDate.toString() +
                                 "greater than the currentDate " + currentDate.toString() +
                                 "are not adjusted - product logic should cope with this case");
        }
        if (rateSpec.fwdRateOffset.get())
            // ??? indexStartDate should be holiday adjusted
            indexStartDate = rateSpec.fwdRateOffset->toDate(resetDate);
        else
            indexStartDate = resetDate;  //  no offset

        // ??? think more about this later
        int nbFwdDays = indexStartDate.daysDiff(resetDate);
        if (maxFwdPeriodIRAdj < nbFwdDays)
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




// Here we always return a slice of value 1.0 in order to allow product
// logic to be generalized
// retrieve zero/discounting slice from the tree
void SingleCcyTree::getZero(DateTime useDate, DateTime matDate, 
                     string curveName, TreeSliceSP &slice) {
    try {
        DateTime currentDate = DateTime::fromIrDate(treeData.TPDate[mTpIdx]);

        if (matDate <= valueDate) {
            throw ModelException("Requested zero maturity date " + matDate.toString() +
                                 " must be greater than the model valueDate " +
                                 valueDate.toString());
        }

        if (useDate != currentDate) {
            throw ModelException("requested zero observation date " + useDate.toString() +
                                 " must be the same as the model's current date " +
                                 currentDate.toString());
        }

        // find the claim bank
        map<string, NamedClaimBank>::iterator cb;
        cb = claimBank.find(curveName);
        if (cb == claimBank.end()) {
            throw ModelException("Unable to find curveName " + curveName +
                                 "in the list of model registered claimBanks");
        }

        // allocate destination slice if needed
        if (!dynamic_cast<TreeSliceRates*>(slice.get())) {
            slice = RateTree::createSimpleSlice(curveName);
            slice->name = "Zero "+useDate.toString()+" to "+matDate.toString();
        }
        TreeSliceRates& sliceL = dynamic_cast<TreeSliceRates&>(*slice);
        sliceL.allocDim(nbFactors);


        CLAIM_BANK const &zeroBank = cb->second.zeroBank;
        double* slicePtr = sliceL.getValuePtr();

        // ??? to keep consistent with hyb3
        if (currentDate == matDate) {
            sliceL = 1.;
        }
        else {
            double* zeroBkPt = NULL;

            // pointer to internal zero slice maturity - zeroBkPt is read only and is not freed
            zeroBkPt = Fix3_ZbkReadZero((CLAIM_BANK*)&zeroBank,
                                        matDate.toIrDate(),
                                        FALSE,  // ??? no interp, do we want to relax this
                                        currentDate.toIrDate(),
                                        mTpIdx,
                                        &treeData);

            if (zeroBkPt == NULL) {
                throw ModelException("Fix3_ZbkReadZero failed: " + IrConverter::slog.pop());
            }
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

void SingleCcyTree::getFXIndex(TreeSlice& treeSlice) {
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


//------------------------------------------------------------------------------
void SingleCcyTree::recordOutput(Control* ctrl, Results* results) const {
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

        DateTimeSP IRSpotDateSP(new DateTime(curveValueDate));
        results->storeGreek(IRSpotDateSP,
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
                            OutputNameConstSP(new OutputName("IR_ZEROS")));

        // store ir discount factors
        CDoubleMatrixSP discFact(new CDoubleMatrix(3, treeData.NbTP+1));
        for (i = 0; i < 3; i++)
            for (j = 0; j <= treeData.NbTP; j++)
                (*discFact)[i][j] = treeData.ZeroCoupon[i][j];

        results->storeGreek(discFact,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("IR_DISCOUNT_FACTORS")));

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


/** forward the price back to the valueDate as the convention is to report values to
    the valueDate, NOT today */
double SingleCcyTree::getPrice0(const TreeSlice& price) const {
    try {
        double vdPrice;
        const TreeSliceLayer * layer = dynamic_cast< const TreeSliceLayer * >( &price );
        if (layer)
            vdPrice = layer->getPrice0();
        else
            vdPrice = price.getCentre();

        // now forward to valueDate using precalculate discount factors between today/valueDate
        // tree is setup so curveToDiffuse = [1], which is what we are reporting the price in
        // so we use that curve's discountFactor to forward back to valueDate
        if (today.getDate() != valueDate.getDate()) {
            vdPrice /= today2ValDateZeros[1];
        }

        return vdPrice;
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


//------------------------------------------------------------------------------
SingleCcyTree::SingleCcyTree(const CClassConstSP &type) :
    RateTree(type), 
    nbFactors(0), 
    zeroInterpStyle(ZeroInterpStyle::FLAT_FWD), 
    cetSmoothing(false),
    legacyMode(false),
    treeBuilt(false)
{
    memset(&treeData,0,sizeof(treeData));
    memset(&devData,0,sizeof(devData));
    memset(&mktVolData,0,sizeof(mktVolData));
    memset(treeCurves,0,sizeof(treeCurves));
    memset(&today2ValDateZeros,0,sizeof(today2ValDateZeros));
    Fix3_Tree_Init(&treeData);
    Fix3_Dev_Init(&devData);
}

SingleCcyTree::~SingleCcyTree() {
    clear();
}

void SingleCcyTree::load(CClassSP& clazz) {

    clazz->setPublic();
    REGISTER(SingleCcyTree, clazz);
    SUPERCLASS(RateTree);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(nbFactors,"Number of interest rate factors for tree");
    FIELD(IRParams, "Collection of interest rate model parameters")
    FIELD(treeDataDebugFile, "Print tree debug data to named output file");
    FIELD_MAKE_OPTIONAL(treeDataDebugFile);
    FIELD(zeroInterpStyle, "Zero curve interpolation style.  Defaults to "
                                  "FLAT_FWD if not supplied");
    FIELD_MAKE_OPTIONAL(zeroInterpStyle);
    FIELD(legacyMode, "starts tree from curve spot date, and thus doesn't extend curves "
          "between valueDate and today.  FALSE by default");
    FIELD_MAKE_OPTIONAL(legacyMode);

    FIELD(cetSmoothing,"");
    FIELD_MAKE_OPTIONAL(cetSmoothing);

    // transient fields - copied but not exported
    FIELD(today, "today/reference date");
    FIELD_MAKE_TRANSIENT(today);
    FIELD(valueDate, "actual value date");
    FIELD_MAKE_TRANSIENT(valueDate);
}

// interface reflection/class loading mechanism
CClassConstSP const SingleCcyTree::TYPE = CClass::registerClassLoadMethod(
    "SingleCcyTree", typeid(SingleCcyTree), SingleCcyTree::load);

bool SingleCcyTreeLoad(void) {
    return (SingleCcyTree::TYPE != 0);
}

DRLIB_END_NAMESPACE
