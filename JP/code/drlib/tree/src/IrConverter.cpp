#include "edginc/config.hpp"
#include "edginc/IrConverter.hpp"
#include "assert.h"
#include "edginc/Actual360.hpp"
#include "edginc/Actual365.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/ActualActual.hpp"
#include "edginc/B30360.hpp"
#include "edginc/B30360F.hpp"
#include "edginc/BootstrappedYieldCurve.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/IRCalib.hpp"
#include "edginc/IRVol.hpp"
#include "edginc/IRVolBase.hpp"
#include "edginc/IRVolPair.hpp"
#include "edginc/IRModelConfig.hpp"
#include "edginc/IREngineTree.hpp"
#include "edginc/IRModelVNFM.hpp"
#include "edginc/IRExoticParam.hpp"
#include "edginc/IRSmile2Q.hpp"
#include "edginc/IRSmileMQ.hpp"
#include "edginc/MaturityTimePeriod.hpp"
#include "edginc/SRMFXVol.hpp"
#include "edginc/Surface.hpp"
#include "edginc/StubBond.hpp"
#include "edginc/StubFactory.hpp"
#include "edginc/StubNone.hpp"
#include "edginc/StubSimple.hpp"
#include "edginc/UntweakableYC.hpp"
#include "edginc/VolCalibInterp.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/YieldCurve.hpp"
#include "irx/zc3.h"
#include "irx/calendar.h"
#include "irx/dateutils.h"
#include "irx/zerocurve.h"
#include "irx/rate.h"
#include "edginc/SRMEQVol.hpp"
#include "edginc/SimpleEquity.hpp"
//#include "edginc/Equity.hpp"
#include "esl_log.h"

DRLIB_BEGIN_NAMESPACE


/***************************** IrConverter ***********************/

StringLogger IrConverter::slog(eslLog());

void IrConverter::DateTimeArrayToIr(const DateTimeArray &dates, vector<long> &irDates) {
    irDates.resize(dates.size());
    for (int i=0; i<dates.size(); ++i) {
        irDates[i] = dates[i].toIrDate();
    }
}

void IrConverter::IrToDateTimeArray(long *irDates, int nbIrDates, DateTimeArray &dates) {
    dates.resize(nbIrDates);
    for (int i=0; i<nbIrDates; ++i) {
        dates[i] = DateTime::fromIrDate(irDates[i]);
    }
}


void IrConverter::checkSlog() {
    // We should not have to do that but for whatever 
    // reason IrConverter::slog gets deleted fron eslLog. Why?
    if (!eslLog().active())
        eslLog().insert(&slog);
}

const char* IrConverter::toIrIndexName(char freq)
{
    switch( freq )
    {
        case 'A': return "12m";
        case 'S': return "6m";
        case 'Q': return "3m";
        case 'M': return "1m";
        default : throw ModelException(__FUNCTION__, "Unsupported index frequency: "
                                        + Format::toString(freq));
    }
}

int IrConverter::eslToQlibFreq(char freq)
{
    switch( freq )
    {
        case 'A': return 1;
        case 'S': return 2;
        case 'Q': return 4;
        case 'M': return 12;
        default : throw ModelException(__FUNCTION__, "Unsupported frequency: "
                                        + Format::toString(freq));
    }
}

void IrConverter::to_NEW_T_CURVE(T_CURVE &tc, const YieldCurve *yc, bool growthCurve)
{
    try {
#if 0

        const UntweakableYC* untweakableYC;
        const BootstrappedYieldCurve *csc;

        untweakableYC = dynamic_cast<const UntweakableYC *>(yc);
        if (untweakableYC)
        {
            bool useGrowthCurve = growthCurve && untweakableYC->growthZeroCurve.get();
            CashFlowArraySP ratesAndDates = (useGrowthCurve ? 
                untweakableYC->growthZeroCurve : untweakableYC->zeroCurve)->getRatesAndDates();

            
            IrxTDate baseDate = yc->getSpotDate().getDate();
            int numDates = ratesAndDates->size();
            
            
            int i,j;
            if (yc->getToday()==(*ratesAndDates)[0].date) {j=1; --numDates;} else j=0;

            IrxTDate* zeroDates = new IrxTDate[numDates];
            double*   zeroRates = new double[numDates];

            for( i = 0; i < numDates; ++i, ++j ) {
                zeroDates[i] = (*ratesAndDates)[j].date.getDate();
                zeroRates[i] = (*ratesAndDates)[j].amount;
            }
           
            // Should be taken from ZeroCurve, like
            // IrxTRateType rateType = zeroCurve->getBasis();
            // IrxTDayCountConv dcc = zeroCurve->getDCC();
            IrxTRateType rateType = IRX_ANNUAL_RATE;
            IrxTDayCountConv dcc = IRX_ACT_365F;

            irxZeroCurveMakeFromRates(baseDate, numDates,
                                      zeroDates, zeroRates, rateType, dcc);
            delete zeroDates;
            delete zeroRates;

            return;
            
        } else {
            csc = dynamic_cast<const BootstrappedYieldCurve *>(yc);
            if (!csc) 
                throw ModelException("YieldCurve should be of BootstrappedYieldCurve or UntweakableYC type");
            
            int i;

            IrxTBootstrapMethod method = IRX_BOOTSTRAP_SMOOTH;
            IrxTBool            reuseDiscountCurve = TRUE;
            IrxTMarketConv      marketConv;
            IrxTDateInterval *  interval;
            /** Date interval for the fixed leg */ // Should have a faster conversion function
            interval = irxDateIntervalMakeFromString(csc->fixedIvl->toString().c_str());
            marketConv.fixedIvl = *interval;
            free(interval);

            /** Date interval for the floating leg */
            interval = irxDateIntervalMakeFromString(csc->floatIvl->toString().c_str());
            marketConv.floatIvl = *interval;
            free(interval);

            /** Date interval for currency basis swaps (into USD) */
            interval = irxDateIntervalMakeFromString(csc->basisIvl->toString().c_str());
            marketConv.cbsIvl = *interval;
            free(interval);

            /** Day count convention for the fixed leg */
            irxDayCountConvFromString(csc->fixedDcc->toString().c_str(), &(marketConv.fixedDcc));
            /** Day count convention for the floating leg */
            irxDayCountConvFromString(csc->floatDcc->toString().c_str(), &(marketConv.floatDcc));
            /** Day count convention for currency basis swaps */
            irxDayCountConvFromString(csc->basisDcc->toString().c_str(), &(marketConv.cbsDcc));
            /** Day count convention for the money market payment */
            irxDayCountConvFromString(csc->moneyMarketDayCount->toString().c_str(), &(marketConv.mmDcc));
            /** Bad day convention for payments */
            irxBadDayConvFromString(csc->badDayConvention->toString().c_str(),&(marketConv.paymentBdc));
            /** Bad day convention for accruals */
            marketConv.accrualBdc = marketConv.paymentBdc;
            /** Bad day convention for reset of floating leg */
            marketConv.resetBdc = marketConv.paymentBdc;
            /** Bad day convention for money market payments.

            Note that if this results in a payment adjustment to on or before the
            start date, then it will automatically switch to use following.
            (Prime example - the 1D rate when the start date is the last business day
            of the month, but not the last calendar day of the month). */
            marketConv.mmBdc = marketConv.paymentBdc;
            /** Number of days to spot */
            marketConv.daysToSpot = csc->spotOffset;

            IrxTCalendar*       calendar = IrConverter::makeCalendar(csc->hols.get());
            IrxTDate            today = csc->today.getDate();

            IrxTDate            baseDate = csc->valueDate.getDate();

            const ZeroCurveBenchmarkArray& benchmarks = csc->benchmarks;
            int                 numInstruments = benchmarks.size();

            const char * *      instTypes = new (const char *[numInstruments]);
            IrxTDate*           maturityDates = new IrxTDate[numInstruments];
            double*             rates = new double[numInstruments];
            double*             cbsSpreads = new double[numInstruments];
            long*                includeFlags = new long[numInstruments];

            for (i=0; i<numInstruments; ++i) {
                // we should have a function ZeroCurveBenchmark::getType()
                ZeroCurveBenchmark *benchmark = benchmarks[i].get();
                if (dynamic_cast<MoneyMarketBenchmark *>(benchmark))
                    instTypes[i] = "M";
                else if (dynamic_cast<SwapBenchmark *>(benchmark))
                    instTypes[i] = "S";
                else 
                    instTypes[i] = "UNSUPPORTED";
                maturityDates[i] = benchmark->getBenchmarkDate(csc->valueDate).getDate();
                rates[i] = benchmark->getRate();
                cbsSpreads[i] = 0.0;
                includeFlags[i] = benchmark->getIncludeFlag();
            }

            IrxTCalendar*       cbsCalendar = NULL;

            const CurrencyBasis * cb = csc->ccyBasis.get();
            if (cb) {
                for (i=0; i<numInstruments; ++i) {
                    // needs review 
                    cbsSpreads[i] = cb->cashBasis(*(benchmarks[i]->getEnd()));
                }

                cbsCalendar = IrConverter::makeCalendar(cb->getHolidays().get());
                //Holiday * holidays = cb->hols.get();

            }


            IrxTSwapZeroCurve * swapCurve = irxZeroCurve3CB(
                method,
                reuseDiscountCurve,
                &marketConv,
                today,
                baseDate,
                numInstruments,
                const_cast<char **>(instTypes),
                maturityDates,
                rates,
                cbsSpreads,
                NULL, // prices,
                includeFlags,
                NULL, // adjustments,
                calendar,
                cbsCalendar,
                0,    // numShortRates,
                NULL, // shortMaturityDates,
                NULL, // shortRates,
                NULL, // shortCbsSpreads,
                NULL, // shortIncludeFlags,
                FALSE,// firstFloatFixed,
                0.0,
                IRX_LONG_FRONT_STUB,    // stubLocation, not used ?
                DateTime::fromIrDate(Nxtmth( DateTime(baseDate,0).toIrDate(), 1200L, 1L )).getDate() // extrapDate
                // We require at least 100 years of zero curve for swaption vol bootstrapping
                );

            delete instTypes;
            delete maturityDates;
            delete rates;
            delete cbsSpreads;
            delete includeFlags;

            irxCalendarFree(calendar);
            irxCalendarFree(cbsCalendar);



        }

        

        // if (growthCurve)
        //     tc = swapCurve->indexCurve;
        // else 
        //     tc = swapCurve->discountCurve;

#endif
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__, "(name="+yc->getName()+")");
    }
}

// Used in IrConverter::to_T_CURVE(...)
int getDayCountDenom(const DayCountConvention* dcc)
{
    int denom;

    if (Actual360::TYPE->isInstance(dcc))
    {
        denom = 360;
    }
    else if (B30360::TYPE->isInstance(dcc))
    {
        denom = 360;
    }
    else if (Actual365::TYPE->isInstance(dcc))
    {
        denom = 365;
    }
    else if (Actual365F::TYPE->isInstance(dcc))
    {
        denom = 365;
    }
    else if (ActualActual::TYPE->isInstance(dcc))
    {
        denom = 365;
    }
    else
    {
        throw ModelException("Unknown MMBasis");
    }

    return denom;
}

// Used in IrConverter::to_T_CURVE(...)
void getDayCountDenomStr(const DayCountConvention* dcc, char* const& denom)
{
    if (Actual360::TYPE->isInstance(dcc))
    {
        ::strcpy(denom, "360");
    }
    else if (B30360::TYPE->isInstance(dcc))
    {
        ::strcpy(denom, "360");
    }
    else if (Actual365::TYPE->isInstance(dcc))
    {
        ::strcpy(denom, "365");
    }
    else if (Actual365F::TYPE->isInstance(dcc))
    {
        ::strcpy(denom, "365");
    }
    else if (ActualActual::TYPE->isInstance(dcc))
    {
        ::strcpy(denom, "ACT");
    }
    else
    {
        throw ModelException("Unknown swapDayCount");
    }
}

void IrConverter::to_T_CURVE(T_CURVE &tc, const YieldCurve *yc, bool growthCurve)
{
    try {
        CashFlowArraySP ratesAndDates;

        // Start date
        tc.Today = tc.ValueDate = yc->getSpotDate().toIrDate();
        tc.SpotDays = 0; // since tc.Today = tc.Valuedate

        const UntweakableYC* untweakableYC;
        const BootstrappedYieldCurve *bsYieldCurve;

        untweakableYC = dynamic_cast<const UntweakableYC *>(yc);
        if (untweakableYC)
        {
            // Money Market basis (360 or 365)
            tc.MMB = untweakableYC->moneyMarketDayCount->daysPerYear();

            // Annual or semi-annual curve ("A" or "S")
            tc.SwapFreq = toEslFreq(untweakableYC->swapFrequency);

            // Year basis for benchmark swaps ("ACT", "365" or "360")
            ::strcpy( tc.SwapDCC,
                (B30360::TYPE->isInstance(untweakableYC->swapDayCount.get()) ? "ACT" /*[sic]*/ :
                Actual360::TYPE->isInstance(untweakableYC->swapDayCount.get()) ? "360" : "365"));
            //getDayCountDenomStr(untweakableYC->swapDayCount.get(), tc.SwapDCC);

            bool useGrowthCurve = growthCurve && untweakableYC->growthZeroCurve.get();
            ratesAndDates = (useGrowthCurve ? untweakableYC->growthZeroCurve :
                                              untweakableYC->zeroCurve)->getRatesAndDates();
        } else {
            bsYieldCurve = dynamic_cast<const BootstrappedYieldCurve *>(yc);
            if (!bsYieldCurve) throw ModelException("YieldCurve should be of BootstrappedYieldCurve or UntweakableYC type");

            // Money Market basis (360 or 365)
            tc.MMB = bsYieldCurve->moneyMarketDayCount->daysPerYear();

            // Annual or semi-annual curve ("A" or "S")

            tc.SwapFreq = toEslFreq(bsYieldCurve->fixedIvl->annualFrequency());

            // Year basis for benchmark swaps ("ACT", "365" or "360")
            ::strcpy( tc.SwapDCC,
                (B30360::TYPE->isInstance(bsYieldCurve->fixedDcc.get()) ? "ACT" /*[sic]*/ :
                Actual360::TYPE->isInstance(bsYieldCurve->fixedDcc.get()) ? "360" : "365"));
            //getDayCountDenomStr(bsYieldCurve->fixedDcc.get(), tc.SwapDCC);

            ratesAndDates = bsYieldCurve->get( bsYieldCurve->isIndexCurve() || growthCurve )->getRatesAndDates();
        }

        tc.NbZero = ratesAndDates->size();
        {
            int i,j;
            if (yc->getToday()==(*ratesAndDates)[0].date) {j=1; --tc.NbZero;} else j=0;

            if( tc.NbZero > MAXNBDATE )
                throw ModelException("Number of rates in zero curve exceeds maximum");

            for( i = 0; i < tc.NbZero; ++i, ++j ) {
                tc.ZeroDate[i] = (*ratesAndDates)[j].date.toIrDate();
                tc.Zero[i] = (*ratesAndDates)[j].amount;
            }
        }

        /* We require at least 100 years of zero curve for swaption vol bootstrapping */
        long lastDate = Nxtmth( tc.ValueDate, 1200L, 1L );

        /* If the zero curve does not extend up to last date we add an extra point */
        if( tc.ZeroDate[ tc.NbZero - 1 ] < lastDate ) {
            if( tc.NbZero < MAXNBDATE ) {
                tc.ZeroDate[ tc.NbZero ] = lastDate;
                /* Flat zero curve */
                tc.Zero[ tc.NbZero ] = tc.Zero[ tc.NbZero - 1 ];
                ++tc.NbZero;
            } else tc.ZeroDate[ tc.NbZero - 1 ] = lastDate;
        }
        if (Term_Check_W(&tc) != SUCCESS ) {
            throw ModelException("Term_Check_W falied: "+IrConverter::slog.pop());
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__, "(name="+yc->getName()+")");
    }
}

double IrConverter::convExpiryToYears(const Expiry *expiry)
{
    try {
        if (!dynamic_cast<const MaturityPeriod*>(expiry)) {
            throw ModelException("Expecting labels (MaturityPeriod), not actual dates (BenchmarkDate)");
        }
        DateTime dateStart(0,0);
        DateTime dateEnd(dateStart);
        return dateStart.yearFrac(expiry->toDate(dateEnd));
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }

}

unsigned IrConverter::convExpiryToMonths(const Expiry *expiry) {
    double x = convExpiryToYears(expiry);
    return (unsigned)(0.5+12*x);
}

void IrConverter::to_SWAPVOL_DATA(SWAPVOL_DATA &sv, ExpiryArray &selectedTenors,
                                  ExpiryArray &selectedExpiries,
                                  const IRVolBase *volBase)
{
    const IRVol *volP = dynamic_cast<const IRVol*>(volBase);
    try {
        if (!volP) {
            if (!volBase) throw ModelException("volBase is NULL");
            throw ModelException("Expected IRVol type, not \""
                +volBase->getClass()->getName()+"\" for IRVolBase ");
        }
        const IRVol &vol = *volP;

        int nbExpir = vol.expiry->size();
        int nbTenor = vol.tenor->size();

        if (nbExpir > NBSWAPTION) 
            throw ModelException("Nb of expiries ("+Format::toString(nbExpir)
            +" rows) exceeds maximum: "+ Format::toString(NBSWAPTION));
        if (nbTenor > NBSWAPTION) 
            throw ModelException("Nb of tenors ("+Format::toString(nbTenor)
            +" columns) exceeds maximum: "+ Format::toString(NBSWAPTION));

        {
            int tenorShifted=0;
            int expiryOffset=0;
            int i;

            // skip expiries smaller than 1 month
            for (expiryOffset = 0; expiryOffset<nbExpir; ++expiryOffset) {
                double months = 12*convExpiryToYears((*vol.expiry)[expiryOffset].get());
                if (months>(29.0/30.0)) break;
            }
            nbExpir -= expiryOffset;
            sv.NbSwaptionExpiries = nbExpir;

            for (i = 0 ; i<nbTenor; ++i) {
                // copy tenor labels
                sv.SwapTenors[tenorShifted] = convExpiryToMonths( (*vol.tenor)[i].get() );

                // skip tenors < 1 year (SWAP vol)
                if (sv.SwapTenors[tenorShifted]<12) continue;

                selectedTenors.push_back((*vol.tenor)[i]);

                // copy the expiries for the current tenor
                for (int j = 0; j<nbExpir; ++j) {
                    sv.VolMatrix[j][tenorShifted] = vol.matrix[i][j+expiryOffset];
                }
                ++tenorShifted;
            }
            sv.NbSwapTenors = tenorShifted;

            // copy expiries labels
            for (i = 0; i<nbExpir; ++i) {
                sv.SwaptionExpiries[i] = convExpiryToMonths((*vol.expiry)[i+expiryOffset].get());
                selectedExpiries.push_back((*vol.expiry)[i+expiryOffset]);
            }

        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__, "IRVolBase="
            +(volP ? volP->getName() : string("?")));
    }
}

void IrConverter::to_BASEVOL_DATA(BASEVOL_DATA &bv,
                                  ExpiryArray &selectedTenors,
                                  ExpiryArray &selectedExpiries,
                                  DateTime& volBaseDate,
                                  const IRVolBase *volBase,
                                  const string &calibIndex)
{
    const IRVol *volP = dynamic_cast<const IRVol*>(volBase);
    try {
        // convert to IRVol
        if (!volP) {
            if (!volBase) throw ModelException("volBase is NULL");
            throw ModelException("Expected IRVol type, not \""
                +volBase->getClass()->getName()+"\" for IRVolBase ");
        }
        const IRVol &vol = *volP;

        volBaseDate = vol.baseDate;

        // check array sizes
        int nbExpir = vol.expiry->size();
        if (nbExpir > MAXNBDATE) throw ModelException("Nb of expiries (rows) exceeds maximum");
        if (vol.tenor->size()<1) throw ModelException("No tenor");

        // get expected tenor
        int calibTenor = atoi(calibIndex.c_str());
        {
            int n= (calibTenor>=10 ? 2:1);
            if (calibTenor>=12
                || (int)calibIndex.size()!=n+1
                || calibIndex[n]!='m') 
            {
                calibTenor = -1;
            }
        }

        // find the right column
        int col=0;
        int nbTenor = vol.tenor->size();
        while (col < nbTenor && (int)convExpiryToMonths((*vol.tenor)[col].get())!=calibTenor) ++col;

        // fill BASEVOL_DATA fields

        if (calibTenor > 0 && col < nbTenor) { // we use base vol
            selectedTenors.push_back((*vol.tenor)[col]);
            bv.Frequency = toEslFreq(12/calibTenor);
            bv.NbVols = nbExpir;
            for (int i = 0; i<nbExpir; ++i) {
                bv.VolDates[i] = (*vol.expiry)[i]->toDate(vol.baseDate).toIrDate();
                bv.Vols[i] = vol.matrix[col][i];
                selectedExpiries.push_back((*vol.expiry)[i]);
            }
        }
        else
        {  
            if (calibTenor != -1) {
                throw ModelException("The basevol provided is not appropriate to "+calibIndex+" calibration");
            }
            // base vol is not used
            bv.Frequency = 0;
            bv.NbVols = 0;
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__, "IRVolBase="
            +(volP ? volP->getName() : string("?")));
    }
}

void IrConverter::to_FXVOLATILITY_DATA(FXVOLATILITY_DATA &fv, const FXAsset *asset)
{
    try {
        if (!asset->fxVol.get())
            throw ModelException("fxVol field (name = " + asset->fxVol->getName() + 
                                 ") of FXAsset name (name = " + asset->getName() + 
                                 ") is NULL pointer - market wrapper may not have populated fxVol correctly");

        const SRMFXVol *vol = dynamic_cast<const SRMFXVol*>(asset->fxVol.get());
        if (!vol) 
            throw ModelException("FX vol market structure must be of type SRMFXVol, type supplied = " +
                                 asset->fxVol.get()->getClass()->getName() + " and name " +
                                 asset->fxVol.get()->getName());

        DateTime today = vol->today;
        if (today != asset->today) 
            throw ModelException("\"today\" not matching in FX vol structure " +
                                 today.toString() + " and asset "+ asset->today.toString());
        fv.ValueDate = today.toIrDate();
        fv.FXSpotRate = asset->spotFX;
        fv.BaseVolFreq=0; // unused;
        int i,nb;

        nb = fv.NbBaseVols = vol->compVolExpiry->size();
        for (i=0; i<nb; ++i) {
            fv.BaseVolDates[i] = (*vol->compVolExpiry)[i]->toDate(today).toIrDate();
            fv.BaseVols[i] = 100*vol->compVol[i];
        }

        if (vol->spotVolExpiry.get()) {
            nb = fv.NbSpotVols = vol->spotVolExpiry->size();
            for (i=0; i<nb; ++i) {
                fv.SpotVolDates[i] = (*vol->spotVolExpiry)[i]->toDate(today).toIrDate();
                fv.SpotVols[i] = 100*vol->spotVol[i];
            }
        } 
        else 
            fv.NbSpotVols=0;
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void IrConverter::to_FXSMILE_DATA(FXSMILE_DATA &fs, const FXAsset *asset)
{
    try {
        const SRMFXVol *vol = dynamic_cast<const SRMFXVol*>(asset->fxVol.get());
        if (!vol) 
            throw ModelException("volBase must be of type SRMFXVol");

        int i, nb;
        nb = fs.NbSmileDates = vol->smileDate.size();
        for (i=0; i<nb; ++i) {
            fs.SmileDates[i] = vol->smileDate[i].toIrDate();
            fs.A1[i] = vol->smileA1[i];
            fs.A2[i] = vol->smileA2[i];
            fs.A3[i] = vol->smileA3[i];
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void IrConverter::fill_FIX3_TREE_DATA(FIX3_TREE_DATA &tree_data)
{
    Fix3_Tree_Init( &tree_data );

    // Model: Nb of standard deviations used to cut the tree
    tree_data.NbSigmaMax = 6;
    // Model: Number of factors
    tree_data.NbFactor = 1;

    /* Hardcoded assigment of zero curves for the engine */
    tree_data.CvDiff = 0;  /* zero.dat is index curve */
    tree_data.CvIdx1 = 1;
    tree_data.CvIdx2 = 2;
    tree_data.CvDisc = 1;  /* disczero.dat is disc curve */
}

void IrConverter::to_MKTVOL_DATA(MKTVOL_DATA &mktvol_data,
                                 BASEVOL_EXPOSURE_DATA *selectedBV,
                                 SWAPVOL_EXPOSURE_DATA *selectedSV,
                                 const string &calibIndex, const T_CURVE &t_curve,
                                 const BASEVOL_DATA *bv, const SWAPVOL_DATA *sv,
                                 bool smoothing)
{
    try {
        const int maxLen=100;
        char index[maxLen];
        strncpy(index, calibIndex.c_str(), maxLen);

        if (EslMktVolCalibrationAndRecording(&mktvol_data, selectedBV,
                selectedSV, index, &t_curve, bv, sv) != SUCCESS) {
            throw ModelException("EslMktVolCalibrationAndRecording falied: "+IrConverter::slog.pop());
        }
        mktvol_data.SmoothingFlag = (smoothing ? 'Y' : 'N');
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void IrConverter::to_MKTVOL_DATA(
    MKTVOL_DATA&           mktvol_data,
    BASEVOL_EXPOSURE_DATA* selectedBV,
    SWAPVOL_EXPOSURE_DATA* selectedSV,
    const T_CURVE&         t_curve,
    bool                   smoothing,
    const IRVolBase*       pVolBase,
    IRVol::VolType         calibVolType,
    IRVol::CalibType       calibType,
    ExpirySP               calibTenor,
    CDoubleSP              calibVolOverride)
{
    static const string method = "IrConverter::to_MKTVOL_DATA";
    try
    {
        // Get Vol
        const IRVol* pVol = dynamic_cast<const IRVol*>(pVolBase);
        if (!pVol)
        {
            if (!pVolBase)
                throw ModelException(method, " - IRVol object is missing!");
            throw ModelException(method, " - Expected IRVol type, not \""
                + pVolBase->getClass()->getName() + "\" for IRVolBase.");
        }

        // Reset number of selected points
        if (selectedBV)
            selectedBV->NbPoints = 0;
        if (selectedSV)
            selectedSV->NbPoints = 0;

        // Initialise market vol data
        MktVol_Init(&mktvol_data);
        mktvol_data.CalibFlag = TRUE;       // Set calib flag
        mktvol_data.SkipFlag = FALSE;       // Vol points skipping not allowed
        mktvol_data.SmoothingFlag = (smoothing ? 'Y' : 'N');    // Set smoothing flag

        // Read conventions from t_curve
        mktvol_data.BaseDate = t_curve.Today;
        mktvol_data.Freq     = t_curve.SwapFreq;

        // Set Daycount number
        if (!strcmp(t_curve.SwapDCC, "360"))
            mktvol_data.DCC = '0';
        else if (!strcmp(t_curve.SwapDCC, "365"))
            mktvol_data.DCC = '5';
        else
            mktvol_data.DCC = '3';

        // Perform interpolation
        VolProcessedCalibSP ptrVolProcessedCalib = VolCalibInterp::interpVolsForCalibration(
                    pVol->getExpiries(), pVol->getTenors(), pVol->getVolMatrix(),
                    pVol->getBaseDate(), calibType, calibTenor, calibVolOverride);

        // Get results
        DoubleArraySP       ptrVols(ptrVolProcessedCalib->getVols());
        DateTimeArraySP     ptrSwapStartDates(ptrVolProcessedCalib->getSwapStartDates());
        DateTimeArraySP     ptrSwapMaturityDates(ptrVolProcessedCalib->getSwapMaturityDates());
        ExpiryArraySP       ptrTenors(ptrVolProcessedCalib->getSelectedTenors());
        ExpiryArraySP       ptrExpiries(ptrVolProcessedCalib->getSelectedExpiries());
        IntArraySP          ptrTenorIndices(ptrVolProcessedCalib->getSelectedTenorIndices());
        IntArraySP          ptrExpiryIndices(ptrVolProcessedCalib->getSelectedExpiryIndices());
        DateTimeArraySP     ptrSelectedSwapMaturityDates(ptrVolProcessedCalib->getSelectedSwapMaturityDates());

        // Get number of vol points
        int i, j, ix;
        mktvol_data.NbVol = ptrVols->size();

        // Setup vol matrix in market vol data
        for (i = 0; i < mktvol_data.NbVol; ++i)
        {
            //volDate = (*ptrExpiries)[i]->toDate(pVol->getBaseDate());
            mktvol_data.VolDate[i] = (*ptrSwapStartDates)[i].toIrDate(); //volDate.toIrDate();
            mktvol_data.SwapSt[i]  = (*ptrSwapStartDates)[i].toIrDate();
            mktvol_data.SwapMat[i] = (*ptrSwapMaturityDates)[i].toIrDate();
            mktvol_data.Vol[i]     = (*ptrVols)[i];
            mktvol_data.VolUsed[i] = 1; // TRUE
        }

        // Setup selected swap vols (for sensitivities)
        for (ix = 0; selectedSV && ix < ptrExpiryIndices->size(); ++ix)
        {
            i = (*ptrExpiryIndices)[ix];
            j = (*ptrTenorIndices)[ix];
            selectedSV->SwaptionExpiryIndices[ix] = i;
            selectedSV->SwaptionExpiryDates[ix] = mktvol_data.VolDate[i];
            selectedSV->SwapTenorIndices[ix] = j;
            selectedSV->SwapMatDates[ix] = (*ptrSelectedSwapMaturityDates)[i].toIrDate();
            selectedSV->NbPoints++;
        }

        // Setup selected base vols (for sensitivities)
        for (ix = 0; selectedBV && ix < ptrExpiryIndices->size(); ++ix)
        {
            i = (*ptrExpiryIndices)[ix];
            selectedBV->VolIndices[ix] = i;
            selectedBV->VolDates[ix] = mktvol_data.VolDate[i];
            selectedBV->NbPoints++;
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, __FUNCTION__);
    }
}


void IrConverter::to_MKTVOL_DATA_TMX_SMILE(
    MKTVOL_DATA&            mktvol_data,
    RateTree::CcyIRParamsSP ccyIRParams)
{
    static const string method = "IrConverter::to_MKTVOL_DATA_TMX_SMILE";
    try
    {
        // Handle TMX Smile Data
        // Lookup Smile
        MarketObjectSP ptrSmile = ccyIRParams->smileTable->getMarketObject(ccyIRParams->smileSet);
        MarketObject*  pSmile = ptrSmile.get();
        if (!pSmile)
            throw ModelException(method, " - Smile object " + ccyIRParams->smileSet + " not found in Smile Table " + ccyIRParams->smileTable.getName() + "!");

        // Handle case for MultiQ
        const IRSmileMQ* pIRSmileMQ = dynamic_cast<const IRSmileMQ*>(pSmile);
        if (!pIRSmileMQ)
            throw ModelException(method, " - Expected IRSmileMQ type, not \""
                + pSmile->getClass()->getName() + "\" for IRExoticParam.");
        if (pIRSmileMQ)
        {                
            // Get Vol
            IRVolRawSP ptrVol = RateTree::getIRVolRaw(ccyIRParams);
            if (!ptrVol.get())
                throw ModelException(method, " - IRVol object is missing!");
            const IRVol* pVol = (ccyIRParams->calibVolType == IRVol::BASE_VOL)
                                ? ptrVol->getBaseVol().get() : ptrVol->getSwapVol().get();
            if (!pVol)
                throw ModelException(method, " - Expected IRVol type, not \""
                    + ptrVol->getClass()->getName() + "\" for IRVolBase.");

            // Get MQ Smile vols: Skew, VolOfVol, BBRP, BBVP, DeltaLeft, DeltaRight, TauLeft, TauRight, LiqudityFlag
            SurfaceConstSP surfaceArray[9];
            surfaceArray[0] = pIRSmileMQ->getSurfaceSkew();
            surfaceArray[1] = pIRSmileMQ->getSurfaceVolOfVol();
            surfaceArray[2] = pIRSmileMQ->getSurfaceBBRP();
            surfaceArray[3] = pIRSmileMQ->getSurfaceBBVP();
            surfaceArray[4] = pIRSmileMQ->getSurfaceDeltaLeft();
            surfaceArray[5] = pIRSmileMQ->getSurfaceTauLeft();
            surfaceArray[6] = pIRSmileMQ->getSurfaceDeltaRight();
            surfaceArray[7] = pIRSmileMQ->getSurfaceTauRight();
            surfaceArray[8] = pIRSmileMQ->getSurfaceLiquidityFlag();

            // Perform interpolation
            VolProcessedCalibSP ptrVolProcessedCalib = VolCalibInterp::interpVolsForCalibration(
                        pVol->getExpiries(), pVol->getTenors(), pVol->getVolMatrix(),
                        pVol->getBaseDate(), ccyIRParams->calibType, ccyIRParams->calibTenor, ccyIRParams->calibVolOverride);

            // Get results
            DoubleArraySP ptrVols = ptrVolProcessedCalib->getVols();
            int expectedSize = ptrVols->size();
            int i, j;

            // Set results for ir vol only (first column)
            for (i = 0; i < expectedSize; ++i)
                mktvol_data.Smile[0][i] = (*ptrVols)[i];

            // Loop through each smile distribution
            for (j = 0; j < 9; ++j)
            {
                if (surfaceArray[j].get())
                {
                    // Perform interpolation
                    VolProcessedCalibSP ptrVolProcessedCalib = VolCalibInterp::interpVolsForCalibration(
                                pIRSmileMQ->getExpiries(), pIRSmileMQ->getTenors(), *(surfaceArray[j]->Zs().get()),
                                pVol->getBaseDate(), ccyIRParams->calibType, ccyIRParams->calibTenor);

                    // Get results
                    ptrVols = ptrVolProcessedCalib->getVols();

                    // Check size
                    if (ptrVols->size() != expectedSize)
                        throw ModelException(method, "Incorrect size received from interpolating the "
                               + surfaceArray[j]->toString() + " distribution in the smile object " + pIRSmileMQ->getName());
                }

                // Setup vols in mktVolData
                if (j == 8)
                {   // Setup liquidity flags
                    for (i = 0; i < expectedSize; ++i)
                        mktvol_data.SmlLiqDate[i] = (ptrVols.get() && (*ptrVols)[i] > 0.0) ? 1 : 0;
                }
                else
                {
                    for (i = 0; i < expectedSize; ++i)
                        mktvol_data.Smile[j+1][i] = (*ptrVols)[i];
                }
            }

            // Normal Cutoff and NCK (MQ numerical configuration)
            mktvol_data.NckMQ = pIRSmileMQ->getNCK();
            mktvol_data.NbSigmaMQ = pIRSmileMQ->getNormalCutoff();
        }                
    }
    catch (exception& e)
    {
        throw ModelException(e, __FUNCTION__);
    }
}

// returns -1 if element (qlib arrays appear to have no find functionality) 
int IrConverter::findNamedElementIndex(const StringArray& elementArray, const string& elementName, const string &errMsg) {

    for (int i = 0; i < elementArray.size(); i++) {
        if (elementArray[i] == elementName)
            return i;
    }
    if (errMsg.size()) 
        throw ModelException(__FUNCTION__, errMsg + elementName);

    return -1;
}


// simple wrapper until removing HYTreeModel as it's easier this way than rewriting it
void IrConverter::to_MODELPARAMETERS_DATA(MODELPARAMETERS_DATA &mp, const string &smileStyle, const string &modelSet, IRCalib &irCalib)
{

    RateTree::IRModelParams mpTmp;
    to_MODELPARAMETERS_DATA(mpTmp, smileStyle, modelSet, irCalib);

    mp.nbFactors = mpTmp.nbFactors;
    mp.QLeft = mpTmp.QLeft;
    mp.QRight = mpTmp.QRight;
    mp.FwdShift = mpTmp.FwdShift;
    mp.CetNbIter = mpTmp.CetNbIter;

    // One factor
    mp.OneFactorPPY = mpTmp.PPY;
    mp.OneFactorMR = mpTmp.FactorMR[0];
    mp.OneFactorVol = mpTmp.FactorVol[0];
   
    // Two factor
    mp.TwoFactorPPY = mpTmp.PPY;
    mp.TwoFactorMR1 = mpTmp.FactorMR[0];
    mp.TwoFactorMR2 = mpTmp.FactorMR[1];
    mp.TwoFactorVol1 = mpTmp.FactorVol[0];
    mp.TwoFactorVol2 = mpTmp.FactorVol[1];
    mp.TwoFactorCorr = mpTmp.FactorCorr[0];

    // Three factor
    mp.ThreeFactorPPY = mpTmp.PPY;
    mp.ThreeFactorMR1 = mpTmp.FactorMR[0];
    mp.ThreeFactorMR2 = mpTmp.FactorMR[1];
    mp.ThreeFactorMR2 = mpTmp.FactorMR[2];
    mp.ThreeFactorVol1 = mpTmp.FactorVol[0];
    mp.ThreeFactorVol2 = mpTmp.FactorVol[1];
    mp.ThreeFactorVol2 = mpTmp.FactorVol[2];
    mp.ThreeFactorCorr12 = mpTmp.FactorCorr[0];
    mp.ThreeFactorCorr13 = mpTmp.FactorCorr[1];
    mp.ThreeFactorCorr23 = mpTmp.FactorCorr[2];
}

void IrConverter::to_MODELPARAMETERS_DATA(RateTree::IRModelParams &mp, const string &smileStyle, const string &modelSet, IRCalib &irCalib)
{
    try
    {
        IRCalib::SmileRequest smileRequest(smileStyle);
        CVolProcessedSP volProcessed(irCalib.getProcessedVol(&smileRequest,0));
        IRCalib::VolProcessed *volData = dynamic_cast<IRCalib::VolProcessed*>(volProcessed.get());
        if (!volData) 
            throw ModelException("volProcessed should be of type IRCalib::VolProcessed");
        const DoubleArray& smileParams = volData->getParams();
        if (smileParams.size()!=4) 
            throw ModelException("Number of ir vol smile params supplied = " + 
                                 Format::toString(smileParams.size()) +
                                 " must supply 4 (QLeft, QRight, FwdShift, nbCETIters)");
        
        mp.QLeft = smileParams[0];
        mp.QRight = smileParams[1];
        mp.FwdShift = smileParams[2];
        mp.CetNbIter = (int)smileParams[3];
        
        IRCalib::ModelRequest modelRequest(modelSet);
        CVolProcessedSP volProcessed2(irCalib.getProcessedVol(&modelRequest,0));
        volData = dynamic_cast<IRCalib::VolProcessed*>(volProcessed2.get());
        if (!volData) 
            throw ModelException("volProcessed should be of type IRCalib::VolProcessed");
        
        const DoubleArray& modelParams = volData->getParams();
        const StringArray& paramLabel = volData->getParamLabel();
        
        if (paramLabel.empty())
            throw ModelException("Must supply optional paramLabel array field for IRCalib::Model[1-3]FL "
                                 "when using rates trees IR model parameters");
        
        mp.nbFactors = 1;
        if (mp.nbFactors == 1)
        {
            int index;
            const string errMsg = "Unable to find IR model parameter: ";

            index = findNamedElementIndex(paramLabel, "meanReversion", errMsg);
            mp.FactorMR[0] = modelParams[index];
            
            index = findNamedElementIndex(paramLabel, "factorWeight", errMsg);
            mp.FactorVol[0] = modelParams[index];

            index = findNamedElementIndex(paramLabel, "ppy", errMsg);
            mp.PPY = (int)modelParams[index];

            index = findNamedElementIndex(paramLabel, "nbStdDevs", errMsg);
            mp.nbStdDevs = (int)modelParams[index];

            index = findNamedElementIndex(paramLabel, "backBone", errMsg);
            mp.backBone = modelParams[index];
        }
        else 
            throw ModelException("Currently only one factor IR model parameters supported");
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


void IrConverter::to_MODELPARAMETERS_DATA(
    RateTree::IRModelParams& mp,
    const string&            smileKey,
    const string&            modelKey,
    const string&            engineKey,
    IRModelConfigTable&       engineTable,
    IRExoticParamTable&        smileTable,
    IRExoticParamTable&        modelTable)
{
    static const string method = "IrConverter::to_MODELPARAMETERS_DATA";
    try
    {
        // Lookup Engine
        MarketObjectSP ptrEngine = engineTable.getMarketObject(engineKey);
        MarketObject*  pEngine = ptrEngine.get();
        if (!pEngine)
            throw ModelException(method, " - Engine object " + engineKey + " not found in Engine Table " + engineTable.getName() + "!");
        const IREngineTree* pIREngineTree = dynamic_cast<const IREngineTree*>(pEngine);
        if (!pIREngineTree)
            throw ModelException(method, " - Object " + engineKey + " is not of type IREngineTree!");

        // Populate engine parameters
        mp.nbStdDevs = pIREngineTree->nbStdDevs;
        mp.nbStateVariables = pIREngineTree->nbStateVars;
        mp.stateVarStdDevs = pIREngineTree->nbStateVarStdDevs;
        mp.CetNbIter = pIREngineTree->nbCET;

        // Lookup Smile
        MarketObjectSP ptrSmile = smileTable.getMarketObject(smileKey);
        MarketObject*  pSmile = ptrSmile.get();
        if (!pSmile)
            throw ModelException(method, " - Smile object " + smileKey + " not found in Smile Table " + smileTable.getName() + "!");

        // Smile2Q
        const IRSmile2Q* pIRSmile2Q = dynamic_cast<const IRSmile2Q*>(pSmile);
        if (pIRSmile2Q)
        {
            // Populate smile parameters
            mp.QLeft = pIRSmile2Q->qLeft();
            mp.QRight = pIRSmile2Q->qRight();
            mp.FwdShift = pIRSmile2Q->fwdShift();
        }
        else if (!dynamic_cast<const IRSmileMQ*>(pSmile))
            throw ModelException(method, " - Object " + smileKey + " is not of type IRSmile2Q or IRSmileMQ!");

        // Lookup Model
        MarketObjectSP ptrModel = modelTable.getMarketObject(modelKey);
        MarketObject*  pModel = ptrModel.get();
        if (!pModel)
            throw ModelException(method, " - Model object " + modelKey + " not found in Model Table " + modelTable.getName() + "!");
        const IRModelVNFM* pIRModelVNFM = dynamic_cast<const IRModelVNFM*>(pModel);
        if (!pIRModelVNFM)
            throw ModelException(method, " - Object " + modelKey + " is not of type IRModelVNFM!");
        
        // Populate model parameters
        const DoubleMatrix& mMR = pIRModelVNFM->meanReversion();
        const DoubleMatrix& mWeight = pIRModelVNFM->weight();
        const DoubleMatrix& mCorr = pIRModelVNFM->correlation();
        mp.nbFactors = pIRModelVNFM->numFactors();
        mp.backBone = pIRModelVNFM->backbone();
        mp.PPY = pIREngineTree->nbPPY;
        
        int i;
        for (i = 0; i < mp.nbFactors; ++i)
        {
            mp.FactorMR[i] = mMR[i][0];
            mp.FactorVol[i] = mWeight[i][0];
        }
        for (i = mp.nbFactors + 1; i < 3; ++i)      // to -999 for unused bits
        {
            mp.FactorMR[i] = -999.;
            mp.FactorVol[i] = -999.;
        }
        for (i = 0; i < mCorr.numCols(); ++i)
            mp.FactorCorr[i] = mCorr[i][0];
        for (i = mCorr.numCols() + 1; i < 3; ++i)
            mp.FactorCorr[i] = -999.;
    }
    catch (exception& e) {
        throw ModelException(e, "IrConverter::to_MODELPARAMETERS_DATA");
    }
}

void IrConverter::checkEnum(const string &oneCharEnum, const char *allowedValues, const char *varName) {
    if (oneCharEnum.size()!=1) {
        throw ModelException(varName+string(" must be 1 character long"));
    }
    if (oneCharEnum.find_first_of(allowedValues)==string::npos) {
        throw ModelException(varName + string("=") + oneCharEnum
        + " is not one of the following allowed values: " + allowedValues);
    }
}

const IRVolBase* IrConverter::getVol(const YieldCurve *yc, bool baseVol) {
    try {
        VolRequestRaw vrr;
        IRVolRawSP vol0(dynamic_cast<IRVolRaw*>(yc->getProcessedVol(&vrr)));
        if (!vol0) throw ModelException("IRVolPair required in market for yc "+yc->getName());
        return (baseVol ? vol0->getBaseVol().get() : vol0->getSwapVol().get());

        /*const IRVolBase *ivb;
        const UntweakableYC *uyc = dynamic_cast<const UntweakableYC *>(yc);
        if (!uyc) {
            const BootstrappedYieldCurve *csc = dynamic_cast<const BootstrappedYieldCurve *>(yc);
            if (!csc) throw ModelException("Expected UntweakableYC or BootstrappedYieldCurve type");
            ivb = csc->irVol.get();
        } else ivb = uyc->irVol.get();
        const IRVolPair *ivp = dynamic_cast<const IRVolPair *>(ivb);
        if (!ivp) throw ModelException("IRVol in "+yc->getName()+" must be of type IRVolPair");
        return (baseVol ? ivp->baseVol.get() : ivp->swapVol.get());
        */
    }
    catch (exception &e) {
        throw ModelException(e, __FUNCTION__);
    }
}

char IrConverter::dccTo035A(const string &dcc) {
    auto_ptr<DayCountConvention> dccObj(DayCountConventionFactory::make(dcc));
    return dccTo035A(*dccObj);
}

char IrConverter::dccTo035A(const DayCountConvention &dcc) {
    if (Actual360::TYPE->isInstance(dcc)) return '0';
    if (B30360::TYPE->isInstance(dcc)) return '3';
    if (Actual365::TYPE->isInstance(dcc)) return '5';
    if (Actual365F::TYPE->isInstance(dcc)) return '5';   // ??? temporary, need to support properly in tree
    if (ActualActual::TYPE->isInstance(dcc)) return 'A';

    throw ModelException(__FUNCTION__,"Unknown DCC \""+dcc.toString()+"\". "
        "Expecting Act/360, 30/360, Act/365, Act/365F or Act/Act");
}

string IrConverter::dccFrom035A(char dcc) {
    switch (dcc) {
        case '0': return "Actual/360";
        case '3': return "B30/360";
        case '5': return "Actual/365";
        case 'A': return "Actual/Actual";
        default: 
            checkEnum(string(dcc,1),"035A", "dccFrom035A(char dcc)");
            return ""; // to remove compiler warning (never reached)
    }
}

char IrConverter::compoundBasisToCoS(int cb) {
    if (cb==CompoundBasis::CONTINUOUS) return 'C';
    if (cb==CompoundBasis::SIMPLE) return 'S';
    throw ModelException(__FUNCTION__,
        "Compount basis \""+Format::toString(cb)+"\"must be CompoundBasis::CONTINUOUS (=5000) "
        "or CompoundBasis::SIMPLE (=0)");
}

char IrConverter::stubConvToBNS(const string &stub) {
    auto_ptr<Stub> stubObj(StubFactory::make(stub));
    if (StubBond::TYPE->isInstance(*stubObj)) return 'B';
    if (StubNone::TYPE->isInstance(*stubObj)) return 'N';
    if (StubSimple::TYPE->isInstance(*stubObj)) return 'S';

    throw ModelException(__FUNCTION__,"Unknown stub \""+stub+"\". "
        "Expecting Bond, None or Simple");
}

char IrConverter::toEslFreq(int freq)
{
    switch(freq) {
        case   1: return ESL_FREQ_ANNUAL; // 'A';
        case   2: return ESL_FREQ_SEMI_ANNUAL; // 'S';
        case   4: return ESL_FREQ_QUARTERLY; // 'Q';
        case  12: return ESL_FREQ_MONTHLY; // 'M';
        case  52: return ESL_FREQ_WEEKLY; // 'W';
        case 365: return ESL_FREQ_DAILY; // 'D';
        default : throw ModelException(__FUNCTION__, 
                      "Unsupported swap frequency: " + Format::toString(freq));
    }
}

char IrConverter::toEslFreq(Expiry *expiry) {
    MaturityPeriod *mp = dynamic_cast<MaturityPeriod *>(expiry);
    if (!mp) throw ModelException(__FUNCTION__, "expiry is not of type MaturityPeriod");
    return toEslFreq(mp->annualFrequency());
}

bool IrConverter::optimizeCriticalDates(
            IRDate vd,
            vector<IRDate> const& cd,
            vector<IRDate> const& cm, vector<IRDate> const& cu,
            vector<IRDate> const& om, vector<IRDate> const& ou,
            vector<IRDate>&       zm, vector<IRDate>&       zu,
            bool optimize)
{
    static char const*  routine = "RateTree::optimizeCriticalDates";

    ESL_LOG(Log::debug) << routine << endl;

    // consistency checks
    assert(cm.size() == cu.size());
    assert(om.size() == ou.size());

    // critical and non-critical schedule sizes
    size_t cSize = cm.size();
    size_t oSize = om.size();

    // require that target vectors are empty
    assert(zm.size() == 0);
    assert(zu.size() == 0);

    // allocate empty critical date array
    int nbCD = 0;
    CRIT_DATE* CD = (CRIT_DATE *)malloc(sizeof(CRITDATE));
    if (!CD)
    {
        ESL_LOG(Log::error) << routine << " - failed to allocate critical dates object";
        return false;
    }

    // add critical dates
    size_t i;
    for (i=0; i<cd.size(); ++i)
    {
        if (Add_To_DateList(&nbCD, &CD, cd[i], 0,0,0,0,0,0,0,0,0) != SUCCESS)
        {
            free(CD);
            ESL_LOG(Log::error) << routine << " - failed to add date to critical dates";
            return false;
        }
    }

    // add critical zero bank dates
    for (i=0; i<cSize; ++i)
    {
        if (Add_To_DateList(&nbCD, &CD, cm[i], 0,0,0,0,0,0,0,0,0) != SUCCESS ||
            Add_To_DateList(&nbCD, &CD, cu[i], 0,0,0,0,0,0,0,0,0) != SUCCESS)
        {
            free(CD);
            ESL_LOG(Log::error) << routine << " - failed to add date to critical dates";
            return false;
        }
    }

    // optimize non-critical dates
    int         nbMat    = 0;
    IRDate*     matDates = NULL;
    int         nbUse    = 0;
    IRDate*     useDates = NULL;

    if (om.size() && ZbkOptDates(om.size(), (long*)&om[0],
                                 ou.size(), (long*)&ou[0],
                                 nbCD, CD,
                                 vd,
                                 &nbMat, &matDates,
                                 &nbUse, &useDates) != SUCCESS)
    {
        free(CD);
        ESL_LOG(Log::error) << routine << " - failed to optimize non-critical dates";
        return false;
    }
                    
    // consistency checks
    assert(nbMat == nbUse);

    int             optSize   = nbMat;
    IRDate const*   optMDates = matDates;
    IRDate const*   optUDates = useDates;

    if (!optimize)
    {
        // with no optimization copy non-critical dates as-is
        optSize   = oSize;
        optMDates = &om[0];
        optUDates = &ou[0];
    }

    // allocate date lists for sorting and use date optimization
    IRDate* matTmp = (IRDate*)malloc(sizeof(IRDate) * (cSize + optSize));
    if (!matTmp)
    {
        free(CD);
        ESL_LOG(Log::error) << routine << " - failed to allocate dates array";
        return false;
    }
    IRDate* useTmp = (IRDate*)malloc(sizeof(IRDate) * (cSize + optSize));
    if (!matTmp)
    {
        free(CD);
        free(matTmp);
        ESL_LOG(Log::error) << routine << " - failed to allocate dates array";
        return false;
    }

    // combine critical zero dates and optimized non-critical dates
    std::copy(cm.begin(), cm.end(), matTmp);
    std::copy(cu.begin(), cu.end(), useTmp);

    std::copy(optMDates, optMDates + optSize, matTmp + cSize);
    std::copy(optUDates, optUDates + optSize, useTmp + cSize);

    // set counters to the total number of dates - CbkProcessDL needs two
    nbMat += cSize;
    nbUse += cSize;

    // sort and remove duplicates from the critical date list
    if (CbkProcessDL(&nbMat, &matTmp, &nbUse, &useTmp) != SUCCESS)
    {
        free(CD);
        free(matTmp);
        free(useTmp);
        ESL_LOG(Log::error) << routine << " - failed to optimize critical dates";
        return false;
    }

    // put combined dates to one vector
    for (int j=0; j<nbMat; ++j)
    {
        zm.push_back(matTmp[j]);
        zu.push_back(useTmp[j]);
    }

    free(CD);
    free(matTmp);
    free(useTmp);

    if (matDates)
        free(matDates);
    if (useDates)
        free(useDates);

    return true;
}

void  IrConverter::to_EQ_DATA(EQ_DATA &eq_data, const SimpleEquity &se) {
    try {
        if (!se.vol)
            throw ModelException("SimpleEquity::vol not populated in market for equity " +
                                 se.getName());
        const SRMEQVol *srmEqVol = dynamic_cast<const SRMEQVol *>(se.vol.get());
        if (!srmEqVol) 
            throw ModelException("SimpleEquity::vol expected to be of type SRMEQVol");

        const DateTimeArray& smileBMDates = srmEqVol->getSmileDate();
        eq_data.NbSmilePt = smileBMDates.size();

        if ((eq_data.NbSmilePt < 0) ||
            (eq_data.NbSmilePt > MAXNBDATE))
            throw ModelException("Nb EQ smile param lines OWS supplied ( " + 
                                Format::toString(eq_data.NbSmilePt) + ")is out of range!");

        for (int i = 0; i < eq_data.NbSmilePt; i++)
        {
            eq_data.SmileDate[i] = smileBMDates[i].toIrDate();
            eq_data.a1[i] = srmEqVol->smileA1[i];
            eq_data.a2[i] = srmEqVol->smileA2[i];
            eq_data.a3[i] = srmEqVol->smileA3[i];
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

IrxTCalendar* makeCalendar(const Holiday *holidays) 
{ 
    return NULL; 
}

void IrConverter::q3GetZerosData(DateTimeArray &dates, DoubleArray &rates, const YieldCurve &yc)
{
    try {
	    const UntweakableYC* untweakableYC;
	    untweakableYC = dynamic_cast<const UntweakableYC *>(&yc);
        if (!untweakableYC)
            throw ModelException("YieldCurve "+yc.getName()+" is not of type UntweakableYC");

	    CashFlowArraySP ratesAndDates = untweakableYC->zeroCurve->getRatesAndDates();
	    int numDates = ratesAndDates->size();
	    int i,j;
        /*if (untweakableYC->getToday()==(*ratesAndDates)[0].date) {j=1; --numDates;} else*/ j=0;

	    dates.resize(numDates);
	    rates.resize(numDates);

	    for( i = 0; i < numDates; ++i, ++j ) 
	    {
		    dates[i] = (*ratesAndDates)[j].date;
		    rates[i] = (*ratesAndDates)[j].amount;
	    }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

// retrieve the zeroDates and rates, and extend curve to today if required
// ??? upgrade to allow cash/swap
void IrConverter::q3GetandExtendZerosData(DateTimeArray &dates, DoubleArray &rates, DateTime &spotDate,
                                          const YieldCurve &yc, const DateTime& today)
{
    try {
        spotDate = yc.getSpotDate();

	    const UntweakableYC* untweakableYC;
	    untweakableYC = dynamic_cast<const UntweakableYC *>(&yc);
        if (!untweakableYC)
            throw ModelException("YieldCurve "+yc.getName()+" is not of type UntweakableYC, "
                "as q3 currently only supports zero Curves - soon to be upgraded");

	    CashFlowArraySP ratesAndDates = untweakableYC->zeroCurve->getRatesAndDates();
	    int numDates = ratesAndDates->size();
	    int i,j;
        if (untweakableYC->getToday()==(*ratesAndDates)[0].date) {
            j=1; 
            --numDates;
        } 
        else 
            j=0;

        if (numDates <= 0)
            throw ModelException("Yield curve " + yc.getName() + " must have at least one "
                "active date and rate entry");
	    dates.resize(numDates);
	    rates.resize(numDates);

        if (today > spotDate)
            throw ModelException("supplied today date " + today.toString() + " must be <= "
                "supplied spot date " + spotDate.toString());
        double stmZero;
        double T2Vrate = (*ratesAndDates)[j].amount;  // use first rate in curve

        // calculate discount factor between today and spot date
        if (today < spotDate)
            stmZero = ::pow((1.0 + T2Vrate), -Daysact(today.toIrDate(), spotDate.toIrDate())/365.0);
        else 
            stmZero = 1.0;

	    for( i = 0; i < numDates; ++i, ++j ) 
	    {
            DateTime tmpDate = (*ratesAndDates)[j].date;
            double T = Daysact(spotDate.toIrDate(), tmpDate.toIrDate())/365.;
            double AZero = ::pow( 1 + (*ratesAndDates)[j].amount , -T);
            
            if (AZero <= TINY)
                throw ModelException("discount factor calculated between spotDate " + spotDate.toString() +
                    " and curve date " + tmpDate.toString() + " is invalid (ie. < TINY)");

            AZero *= stmZero;   // convert to discount factor from today

            // convert back to rate
            double ARate = ::pow(AZero, -365.0/Daysact(today.toIrDate(), tmpDate.toIrDate())) - 1.0;  

		    dates[i] = (*ratesAndDates)[j].date;
		    rates[i] = ARate;
	    }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

string IrConverter::StringFromIRDate(IRDate date)
{
    char buffer[12];
    ::StringFromIRDate(date, buffer);
    return string(buffer);
}

DRLIB_END_NAMESPACE
