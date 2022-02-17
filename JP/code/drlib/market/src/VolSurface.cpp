//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolSurface.cpp
//
//   Description : Surface Based Implementation of Vol Interface
//
//   Author      : Mark A Robson
//
//   Date        : 21 Jan 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_VOLSURFACE_CPP
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/LinearStrikeSpreadVolRequest.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/CliquetVolRequest.hpp"
#include "edginc/BenchmarkDate.hpp"
#include "edginc/CompositeVol.hpp"
#include "edginc/Addin.hpp"
#include "edginc/PDFDefaultLNStrike.hpp"
#include "edginc/VolSpline.hpp"
#include "edginc/VolRequestTime.hpp"
#include "edginc/VolProcessedDispatch.hpp"
#include "edginc/IStrikeMap.hpp"

DRLIB_BEGIN_NAMESPACE

const string VolSurface::ABSOLUTE = "Absolute";
const string VolSurface::SPOTMONEYNESS = "SpotMoneyness";
const string VolSurface::FWDMONEYNESS = "FwdMoneyness";

    // InterpList contains the information required to interpolate on 
    // each benchmark at some constant function of (strike, T = option maturity) 
    // as defined by the passed strikeMap. The information is stored as 
    // (lower idx, frac, upper idx) for each benchmark.
    // For example, if strikeMap is a FwdMoneyness, InterpList will define 
    // the point on each strike axis in order to interpolate at
    // FwdMoneyness = K/F(T) across benchmarks
    // Associate with VolSurface?
    class MARKET_DLL InterpList : public CObject {
    public:
        static  CClassConstSP const TYPE;
        IntArraySP lower;
        IntArraySP upper;
        DoubleArraySP frac;
        //double strike;

        static void getLinearInterpFrac(double level,               // I
                                        const DoubleArray& levels,  // I
                                        int& lower,                 // O
                                        int& upper,                 // O
                                        double& frac)               // O
        {
            if (level <= levels[0]) {
                lower = 0;
                upper = 0;
                frac  = 0.0;
            } else if (level >= levels.back()) {
                lower = levels.size()-1;
                upper = levels.size()-1;
                frac  = 0.0;
            } else {
                /* Note we have already checked for case when strike is above
                the upper bound of the surface => loop terminates okay */
                int j;
                for (j = 1; level > levels[j]; j++)
                {
                    ; /* empty loop */
                }
        
                /* interpolate between i and i - 1. Value is y1 +
                (x - x1)*(y2 - y1)/(x2 - x1). Can cache part of this */
                frac  = (level - levels[j-1])/(levels[j] - levels[j-1]);
                lower = j - 1;
                upper = j;
            }
        }

        // Builds up a list of interpolation indicies based on the mapped strike
        InterpList(const VolSurface* surface, 
                   IStrikeMapSP strikeMap, 
                   const DateTime& maturity,
                   double strike) : CObject(TYPE) {
            
            double level = strikeMap->fromStrike(strike, maturity);
            DateTimeArray dates = surface->getDates();
            DoubleArray strikes = surface->getStrikes();
            DoubleArray levels(strikes.size());
            upper = IntArraySP(new IntArray(dates.size()));
            lower = IntArraySP(new IntArray(dates.size()));
            frac = DoubleArraySP(new DoubleArray(dates.size()));
            int i, j;

            for (i=0; i<dates.size(); i++) {
                // Map strikes at benchmark i
                for (j=0; j<strikes.size(); j++) {
                    levels[j] = strikeMap->fromStrike(strikes[j], dates[i]);
                }

                InterpList::getLinearInterpFrac(level,
                                                levels,
                                                (*lower)[i],
                                                (*upper)[i],
                                                (*frac)[i]);
            }
        }
    private:
        InterpList::InterpList() : CObject(TYPE) {}
        static void load(CClassSP& clazz) {
            REGISTER(InterpList, clazz);
            SUPERCLASS(CObject);
            EMPTY_SHELL_METHOD(defaultInterpList);
            FIELD_NO_DESC(lower);
            FIELD_MAKE_TRANSIENT(lower);
            FIELD_NO_DESC(upper);
            FIELD_MAKE_TRANSIENT(upper);
            FIELD_NO_DESC(frac);
            FIELD_MAKE_TRANSIENT(frac);
          //  FIELD_NO_DESC(strike);
          //  FIELD_MAKE_TRANSIENT(strike);
        }
        static IObject* defaultInterpList() {
            return new InterpList();
        }
    };

    typedef smartConstPtr<InterpList> InterpListConstSP;
    typedef smartPtr<InterpList> InterpListSP;

CClassConstSP const InterpList::TYPE = 
    CClass::registerClassLoadMethod(
    "InterpList", 
    typeid(InterpList), load);

class VolSurface::Interp: public CVolProcessedBS{
public:
    static CClassConstSP const TYPE;
    VolSurfaceConstSP   volSurf;    /* vol (sometimes a ref to original). */
    bool           fwdStartInterp;  /* derived value */
    double         interpLevel;  /* either absolute strike or
                                    moneyness if fwd starting */
    double         spot;         /* if fwdStart */
    DateTime       startDate;    /* 'option' start date if fwd start */
    DateTime       endDate;      /* 'option' end date if fwd start */
    bool           allowNegativeFwdVar; /* true: don't fail if fwd variance
                                           negative (unless doing vols) */
    bool           isManualSpread; /* if true
                                    then spread below is added regardless
                                    of whether fwd starting or not */
    double         spread;       /* shift for when manualSpread = true */
    // strikeMap is used for mapping a function of strike (e.g. FwdMoneyness)
    // into strike and vice versa. Needed to interpolate
    // along some constant quantity = F(strike, option maturity) in time 
    IStrikeMapSP   strikeMap;   

    void  validatePop2Object(){
        // Default strikeMap
        if (!strikeMap) {
            strikeMap = IStrikeMapSP(new StrikeMapNone());
        }
    }

//    Interp(const VolSurface*  vol):
//        CVolProcessedBS(TYPE), volSurf(VolSurfaceConstSP::attachToRef(vol)),
//        spot(0.0), 
//        allowNegativeFwdVar(false), isManualSpread(false),
//        spread(0.0), strikeMap(new StrikeMapNone()) {}

    Interp(const VolSurface*  vol, IStrikeMapSP strikeMap):
        CVolProcessedBS(TYPE), volSurf(VolSurfaceConstSP::attachToRef(vol)),
        spot(0.0), 
        allowNegativeFwdVar(false), isManualSpread(false),
        spread(0.0), strikeMap(strikeMap) {}
    
    /** identifies the market data name of the volatility */
    virtual string getName() const {
        return volSurf->getName();
    }

    /** retrieve time measure for the vol */
    virtual TimeMetricConstSP GetTimeMetric()const{
        return volSurf->metric;
    }

    /** calculates the trading time between two dates */
    virtual double calcTradingTime(const DateTime &date1, 
                                   const DateTime &date2) const {
        return volSurf->metric->yearFrac(date1, date2);
    }

    /** Populate a CompositeVol with data it requires (essentially
        benchmark dates and vols) */
    virtual void populateCompositeVol(CompositeVol* compositeVol) const{
        if (!compositeVol){
            throw ModelException("Interp::populateCompositeVol",
                                 "NULL compositeVol");
        }
        if (fwdStartInterp){
            /* forward starting - don't bother calculating the vols now */
            compositeVol->addDates(volSurf->baseDate, volSurf->expiries.get());
        } else {
            const CDoubleMatrix& vol = *volSurf->vol;
            CDoubleArray  volsAtExpiries(volSurf->expiries->size());
            /* now just do linear interpolation at benchmark dates */
            TInterp        interpCache;
            strikeInterpLinearDouble(volSurf->strikes,
                                     interpLevel,
                                     allowNegativeFwdVar,
                                     &interpCache);
            for (int i = 0; i < volSurf->expiries->size(); i++) {
                volsAtExpiries[i] = volAtStrikeExpiry(interpLevel,
                                                      i,
                                                      endDate);
                //volsAtExpiries[i] = vol[interpCache.lower][i] + 
                //    interpCache.frac *
                //    (vol[interpCache.upper][i] - vol[interpCache.lower][i]);
            }
            // then supply expiries and vols to composite vol
            compositeVol->addVols(volSurf->baseDate,
                                  volSurf->expiries.get(), volsAtExpiries);
        }
    }

    /* interpolate vol */
    void CalcVol(const DateTimeArray& dateList, 
                 TCalcType            calcType, 
                 CDoubleArray&        vols) const{
        const static string routine = "Interp::CalcVol";
        try{
            if (calcType == forward) {
                if (dateList.size() < 2){
                    throw ModelException(routine, 
                                         "Must supply at least two dates");
                }
                if (vols.size() < dateList.size()-1){
                    throw ModelException(routine, 
                                         "Double array too short");
                }
                if (fwdStartInterp) {
                    /* do the forward starting version */
                    fwdStartInterpVol(dateList,
                                      true,    /* do vols */
                                      allowNegativeFwdVar,
                                      vols);
                } else {
                    surfInterpVol(interpLevel,
                                  allowNegativeFwdVar,
                                  dateList,
                                  true,    /* do vols */
                                  vols);
                }
            }
            else {
                CVolProcessedBS::CalcVol(dateList, calcType, vols);
            }

        } catch (exception& e){
            throw ModelException(e, routine);
        }
        return;
    }


    void CalcVar(const DateTimeArray& dateList,
                 TCalcType            calcType, 
                 CDoubleArray&        vars) const{
        const static string routine = "Interp::CalcVar";
        try{
            if (dateList.size() < 2){
                throw ModelException(routine, 
                                     "Must supply at least two dates");
            }
            if (vars.size() < dateList.size()-1){
                throw ModelException(routine, 
                                     "Double array too short");
            }
            if (fwdStartInterp) {
                /* do the forward starting version */
                fwdStartInterpVol(dateList,
                                  false,    /* do variance */
                                  allowNegativeFwdVar,
                                  vars);
            } else {
                surfInterpVol(interpLevel,
                              allowNegativeFwdVar,
                              dateList,
                              false,    /* do variance */
                              vars);
            }
            
            switch (calcType){
                int i; /* MSVC is broken - doesn't support separate
                          variables in loop */
            case forward:
                // do nothing
                break;
            case fromFirst:
                for (i = 1; i < dateList.size()-1; i++){
                    vars[i] += vars[i-1];
                }
                break;
            case toLast:
                for (i = dateList.size()-2; i >= 0; i--){
                    vars[i] += vars[i+1];
                }
                break;
            default:
                throw ModelException(routine, "Unknown calculate type");
            }
        } catch (exception& e){
            throw ModelException(&e, "Interp::CalcVar");
        }
    }

    /** Calculates volatility between 2 dates */
    double CalcVol(const DateTime& date1, const DateTime& date2) const{
        // might be a bit slower than possible
        DateTimeArray dateList(2);
        CDoubleArray  vols(1);
        dateList[0] = date1;
        dateList[1] = date2;
        CalcVol(dateList, forward, vols);
        return vols[0];
    }
                
    /** Calculates variance between 2 dates */
    double CalcVar(const DateTime& date1,
                   const DateTime& date2) const{
        // might be a bit slower than possible
        DateTimeArray dateList(2);
        CDoubleArray  vars(1);
        dateList[0] = date1;
        dateList[1] = date2;
        CalcVar(dateList, forward, vars);
        return vars[0];
    }

    /* interpolate a vol surface - approach 1 */
    static Interp* interpolateLinearStrike(
        const VolSurface*             surface,
        const LinearStrikeVolRequest* interp,// (I) 
        const CAsset*                 asset){ // (I) access to asset market info

        IStrikeMapSP noMap = IStrikeMapSP(new StrikeMapNone());
        smartPtr<Interp> interpVol(new Interp(surface, noMap));
        interpVol->allowNegativeFwdVar = interp->negativeFwdVarAllowed();
        interpVol->fwdStartInterp = interp->getStartDate().
            isGreater(surface->baseDate);
        interpVol->isManualSpread = false;
        if (interpVol->fwdStartInterp) {
            interpVol->startDate = interp->getStartDate();
            interpVol->endDate   = interp->getEndDate();
            /* get spot price while we have access to the asset */
            interpVol->spot = asset->getSpot();
        }
        interpVol->interpLevel = interp->getInterpLevel(interpVol->
                                                        fwdStartInterp);
        return interpVol.release();
    }

public:
    /* interpolate a vol surface - approach 1 simplified */
    static Interp* interpolateAbsLinearStrike(const VolSurface* surface,
                                              double            absoluteStrike){
        IStrikeMapSP noMap = IStrikeMapSP(new StrikeMapNone());
        Interp* interpVol = new Interp(surface, noMap);
        interpVol->interpLevel = absoluteStrike;
        interpVol->fwdStartInterp = false;
        interpVol->isManualSpread = false;
        return interpVol;
    }

    /* interpolate a vol surface - approach 2 */
    static Interp* interpolateATM(
        const VolSurface*    surface,
        const ATMVolRequest* request,
        const CAsset*        asset){// (I) access to asset market info

        IStrikeMapSP noMap = IStrikeMapSP(new StrikeMapNone());
        smartPtr<Interp> interpVol(new Interp(surface, noMap));
        // interpolate at the asset's spot
        interpVol->interpLevel = asset->getSpot();
        interpVol->fwdStartInterp = false;
        interpVol->isManualSpread = false;
        return interpVol.release();
    }
    
    /* interpolate a vol surface - approach 3 */
    static Interp* interpolateWithTermStructure(
        const VolSurface*               surface,
        const LinearStrikeTSVolRequest* interp, 
        const CAsset*                   asset){  // access to asset market info
        IStrikeMapSP noMap = IStrikeMapSP(new StrikeMapNone());
        smartPtr<Interp> interpVol(new Interp(surface, noMap));        
        interpVol->isManualSpread = false;
        interpVol->allowNegativeFwdVar = interp->negativeFwdVarAllowed();

        interpVol->fwdStartInterp = interp->getStartDate().
            isGreater(surface->baseDate);
        interpVol->interpLevel = interp->getInterpLevel(interpVol->
                                                        fwdStartInterp);
        // deal with the non-forward starting case first
        if (!interpVol->fwdStartInterp)
        {
            // just return a reference to the original surface
            interpVol->volSurf = VolSurfaceConstSP::attachToRef(surface);
        }
        else   // forward starting
        {
            /* currently (for ease) build a new vol if fwd starting -
               this is slow.  New vol is a curve and just use normal
               interpolation on it */
            /* get spot price while we have access to the asset */
            interpVol->spot = asset->getSpot(); 
            interpVol->startDate = interp->getStartDate();
            interpVol->endDate = interp->getEndDate();

            interpVol->volSurf = interpVol->interpolateFSWithTermStructure(
                interpVol->interpLevel,
                interpVol->spot,
                interpVol->allowNegativeFwdVar,
                interpVol->startDate);
        }
        
        return interpVol.release();
    }

    static Interp* interpolateCliquet(
        const VolSurface*          surface,
        const CliquetVolRequest*   interp,   // (I) 
        const CAsset*              asset) {  // (I) access to asset market info

        static const string method = "Interp::interpolateCliquet";

        int          numNewBMDates, numNewBMDatesCliq;
        DateTime     intervalEndDate;
        int          idx, cliqIdx;
        double       yearsToStart;
        int          firstCliqIdx;
        double       varToStart;
        IStrikeMapSP noMap = IStrikeMapSP(new StrikeMapNone());
        smartPtr<Interp> interpVol(new Interp(surface, noMap));
        interpVol->isManualSpread = false;
        interpVol->allowNegativeFwdVar = interp->negativeFwdVarAllowed();
        /* interp level is actually irrelevant */
        interpVol->interpLevel = 0.0;

        /* currently (for ease) build a new vol 
           New vol is a curve and just use normal interpolation on it */
        interpVol->fwdStartInterp = false;

        double spot = asset->getSpot();

        const DateTimeArray& liveCliqStartDates = interp->getCliqStartDates();
        bool firstCliqStarted = surface->getBaseDate().
            isGreaterOrEqual(liveCliqStartDates[0]);
        const DoubleArray&   interpLevels = 
            interp->getInterpLevels(!firstCliqStarted);

        /* allocate more than enough space */
        numNewBMDates = liveCliqStartDates.size() * 
            (surface->getDates().size()+1);
        DoubleArray vols(numNewBMDates);
        DateTimeArray bmDates(numNewBMDates);

        const DateTimeArray& surfaceDates = surface->getDates();

        if (firstCliqStarted) {
            /* firstCliqHasStarted - treat separately since vols are
             computed differently, and there's special account for
             baseDate */
            intervalEndDate = (liveCliqStartDates.size() ==1)?
                surfaceDates.back(): liveCliqStartDates[1];

            /* build bmDates for this interval stopping after
               intervalEndDate */
            ExpiryArrayConstSP  bmExpiriesCliq(
                interpVol->bmExpiriesForFSWithTermStructure(surface->getBaseDate(),
                                                 intervalEndDate,
                                                 numNewBMDatesCliq));
            DateTimeArray bmDatesCliq(0);

            /* intervalStartDate is baseDate (required to take advantage of 
             * BMDatesForFSWithTermStructure), but bmdates should
             * not include the base date, so filter out here.
             * Started interp method - relies on interpLevel being absolute */
            for (idx = 0; idx < bmExpiriesCliq->size()-1; ++idx) {
                bmDatesCliq.push_back((*bmExpiriesCliq)[idx]->
                                      toDate(surface->getBaseDate()));
            }

            /* cut at end of first cliquet and add a BM date there */
            if ( (*bmExpiriesCliq)[numNewBMDatesCliq-1]->
                 toDate(surface->getBaseDate()).isGreater(intervalEndDate) ) {
                bmDatesCliq.push_back(intervalEndDate);
            } else {
                bmDatesCliq.push_back((*bmExpiriesCliq)[numNewBMDatesCliq-1]->toDate(surface->getBaseDate()));
            }

            interpVol->surfInterpVol(interpLevels[0], 
                          interpVol->allowNegativeFwdVar,
                          bmDatesCliq,
                          true, /* want vols */
                          vols);

            for (idx = 0; idx < numNewBMDatesCliq-1; ++idx) {
                bmDates[idx] = bmDatesCliq[idx+1];
            }


            numNewBMDates = numNewBMDatesCliq-1; /* start the true count */
            firstCliqIdx = 1;
        } else {
            numNewBMDates = 0; /* start the true count */
            firstCliqIdx = 0;
        }

        for(cliqIdx=firstCliqIdx;cliqIdx<liveCliqStartDates.size();cliqIdx++)
        {
            intervalEndDate = (cliqIdx==liveCliqStartDates.size()-1)?
                surfaceDates.back():liveCliqStartDates[cliqIdx+1];

            /* build bmDates for this interval stopping after intervalEndDate */
            ExpiryArrayConstSP  bmExpiriesCliq(interpVol->bmExpiriesForFSWithTermStructure(liveCliqStartDates[cliqIdx],
                                                      intervalEndDate,
                                                      numNewBMDatesCliq));

            DateTimeArray bmDatesCliq(0);
            for (int i=0 ; i < bmExpiriesCliq->size() ; ++i) {
                bmDatesCliq.push_back((*bmExpiriesCliq)[i]->toDate(surface->baseDate));
            }

            /* Extend the full curve across this interval. Translate into
             * standard spot vols as we go. Note vols[] at each cliquet date
             * is defined from the previous interval. */
            DoubleArraySP volsCliq(
                interpVol->fwdVolsForFSWithTermStructure(interpLevels[cliqIdx], spot, 
                                              interpVol->allowNegativeFwdVar,
                                              liveCliqStartDates[cliqIdx],
                                              numNewBMDatesCliq, bmDatesCliq));

            if (cliqIdx==0) { 
                /* first benchmark date at first start date */
                bmDates[0]    = bmDatesCliq[0];
                vols[0]       = (*volsCliq)[0];
                numNewBMDates = 1;
            } else if (numNewBMDates<1) {
                /* continue from end of previous cliquet - which should exist! */
                throw ModelException(method, "INTERNAL ERROR : numNewBMDates is " + 
                    Format::toString(numNewBMDates) + " < 1!\n");
            }

            /* calculate variance between base and start date - need calendar time
               from base date to start date */
            yearsToStart = surface->metric->yearFrac(surface->baseDate,
                                                     liveCliqStartDates[cliqIdx]);
            varToStart = vols[numNewBMDates-1] * vols[numNewBMDates-1] * yearsToStart;

            for (idx = 1; idx < numNewBMDatesCliq; ++idx) {
                double  yearFrac = surface->metric->yearFrac(liveCliqStartDates[cliqIdx], 
                                                             bmDatesCliq[idx]);
                vols[numNewBMDates+idx-1] = sqrt((varToStart + (*volsCliq)[idx] * (*volsCliq)[idx] * yearFrac)/
                                                 (yearFrac + yearsToStart));
                bmDates[numNewBMDates+idx-1] = bmDatesCliq[idx];
            }
            numNewBMDates += numNewBMDatesCliq-1;

            /* if cliq ends not on benchmark truncate what we've got and correct the vol */
            idx = numNewBMDates-1; /* reuse for brevity&clarity */
            if (bmDates[idx].isGreater(intervalEndDate)) {
                double  t1 = surface->metric->yearFrac(surface->baseDate, bmDates[idx-1]);
                double  var1 = vols[idx-1] * vols[idx-1] * t1;
                double  t2 = surface->metric->yearFrac(surface->baseDate, bmDates[idx]);
                double  var2 = vols[idx] * vols[idx] * t2;
                double  t = surface->metric->yearFrac(surface->baseDate, intervalEndDate);
                /* linearly interpolate variance then extract vol */
                vols[idx] = sqrt((var1*(t2-t)+var2*(t-t1))/(t*(t2-t1)));
                bmDates[idx] = intervalEndDate;
            }
        }

        // we have possibly overallocated memory - reduce size to number of benchmark dates
        vols.resize(numNewBMDates);

        DoubleArray strike(1);
        strike[0] = interpLevels[0];
        DoubleMatrix volsAsMatrix(vols);

        ExpiryArray expiryBMDates(0);

        for (int i=0 ; i<numNewBMDates ; ++i) {
            expiryBMDates.push_back(ExpirySP(new BenchmarkDate(bmDates[i])));
        }


        /* "strike" is irrelevant since this is a curve but possibly useful for error messages */
        interpVol->volSurf = VolSurfaceConstSP(new VolSurface(surface->getName(),
                                                              surface->metric.get(),
                                                              strike,
                                                              volsAsMatrix,
                                                              &expiryBMDates,
                                                              surface->baseDate));

        return interpVol.release();
    }

    /* interpolate a vol surface - approach 5 */
    static Interp* interpolateLinearStrikeSpread(
        const VolSurface*                   surface,
        const LinearStrikeSpreadVolRequest* interp,// (I) 
        const CAsset*                       asset){// (I) access to asset market info

        IStrikeMapSP noMap = IStrikeMapSP(new StrikeMapNone());
        smartPtr<Interp> interpVol(new Interp(surface, noMap));
        interpVol->interpLevel = interp->getInterpLevel(true);
        interpVol->fwdStartInterp = true; // so we use our spread
        interpVol->isManualSpread = true;
        interpVol->allowNegativeFwdVar = interp->negativeFwdVarAllowed();
        interpVol->spread = interp->getSpread();
        interpVol->startDate = interp->getStartDate();
        interpVol->endDate   = interp->getEndDate();
        /* get spot price while we have access to the asset */
        interpVol->spot = asset->getSpot();
        return interpVol.release();
    }

    /* spline surface then call getProcessedVol on that */
    static IVolProcessed* interpolateSpline(
        const VolSurface*      surface,
        const CVolRequestDVF*  request,// (I) 
        const CAsset*          asset){// (I) access to asset market info
        CVolBaseSP spline(new VolSpline(*surface, asset->getSpot()));
        return spline->getProcessedVol(request, asset);
    }
    /** just want access to trading time */
    static IVolProcessed* interpolateTime(
        const VolSurface*      surface,
        const VolRequestTime*  request,// (I) 
        const CAsset*          asset){// (I) access to asset market info
        return VolRequestTime::createVolProcessed(surface->getName(),
                                                  surface->metric);
    }

public:
    /* structure to hold genuine cached strike interpolation data - only to be
       used in this file with structure placed on the stack */
    typedef struct _TInterp
    {
        int     lower;      /* lower idx of axis */
        int     upper;      /* upper idx of axis */
        double  frac;       /* fraction between above two */
        double  strike;     /* only for error reporting */
        bool    allowNegativeFwdVar; // true: don't fail on negative fwd var
    } TInterp;

    /* get different numbers on NT optimised to NT debug and solaris 
       debug/optimised. Not worth the performance gain so turn off */
#if defined(_MSC_VER)
#pragma optimize("", off)
#endif

    /* returns bounding points of a value on a double axis (not
       checked for) together with the fraction the point is between these
       bounding points. If interpLevel is off bottom of axis lowerIdx =
       upperIdx = 0. If off top then lowerIdx = upperIdx = last value of
       array index on axis. In these two cases fraction = 0. If axis only
       has one point then lowerIdx = upperIdx = 0 */
    static void strikeInterpLinearDouble(
        const CDoubleArray& strikes, /* (I) */
        double              level,   /* (I) where to interpolate */
        bool                allowNegativeFwdVar, // (I) don't fail if var<0
        TInterp            *interp){ /* (O) interp data for strike axis */

        interp->strike = level;  /* record for error reporting */
        interp->allowNegativeFwdVar = allowNegativeFwdVar;

        InterpList::getLinearInterpFrac(level,
                                        strikes,
                                        interp->lower,
                                        interp->upper,
                                        interp->frac);
        return;
    }

private:

    /* determine the spread to apply to forward vols when
     * interpolating a forward starting vol
     */
    double fwdStartVolSpread() const{   /* (I)  */
        double        spreadToUse;
        
        if (isManualSpread){
            /* pre-supplied spread */
            spreadToUse = spread;
        } else if (Maths::equals(interpLevel, 1.0)){
            /* struck at-the-money, so zero spread */
            spreadToUse = 0.0;
        } else {
            const VolSurface *surface = volSurf.get();
            double        spotLevel;
            double        fwdStrike;
            double        spotLevelVol;
            double        strikeLevelVol;

            /* need to interpolate at the 'effective maturity'
             * i.e. if it's a 6mo option, get the 6mo vol */
            DateTime interpDate = surface->baseDate.add(
                endDate.subtract(startDate));

            /* first interpolate at the "fwd strike", i.e. at the
             * moneyness * spot price */
            spotLevel = spot;
            fwdStrike = interpLevel * spotLevel;
            
            try{
                DateTimeArray dates(2);
                DoubleArray   vols(1);
                dates[0] = surface->baseDate;
                dates[1] = interpDate;
                surfInterpVol(fwdStrike,
                              allowNegativeFwdVar,
                              dates,
                              true,    /* do vols */
                              vols);
                strikeLevelVol = vols[0];

                /* now interpolate at spot */
                surfInterpVol(spotLevel,
                              allowNegativeFwdVar,
                              dates,
                              true,    /* do vols */
                              vols);
                spotLevelVol = vols[0];
            } catch (exception& e){
                throw ModelException(&e, "Interp::fwdStartVolSpread");
            }
            spreadToUse = strikeLevelVol - spotLevelVol;
        }
        return spreadToUse;
    }

    /* Computes new benchmark dates with correct intervals from start
     * date enforcing a prescribed last benchmark date */
    ExpiryArray* bmExpiriesForFSWithTermStructure(
        const DateTime&            startDate,
        const DateTime&            lastBMDate,
        int&                       numNewBMDates){
     
        static const string method = 
            "Interp::bmExpiriesForFSWithTermStructure";

        ExpiryArraySP bmExpiries(new ExpiryArray());
        ExpirySP firstExpiry(new BenchmarkDate(startDate));
        bmExpiries->push_back(firstExpiry);
        
        if (startDate.isGreaterOrEqual(lastBMDate))
        {
            throw ModelException(method,
                                 "Option start date " +
                                 startDate.toString() +
                                 " must be before last bench mark date " +
                                 lastBMDate.toString());
                                 
        }

        DateTime::Interval interval = startDate.subtract(volSurf->baseDate);
        numNewBMDates = 1;
        do
        {
            ExpirySP expiry(new BenchmarkDate(
                (*volSurf->dates)[numNewBMDates-1].add(interval)));
            bmExpiries->push_back(expiry);
            numNewBMDates++;
        }
        while (lastBMDate.isGreater((*bmExpiries)[numNewBMDates-1]->toDate(startDate)));
        
        // numNewBMDates now contains correct count

        return bmExpiries.release();
    }

    /* Computes the "fwd vols" which measure from start date as an internal part
       of the "FSWithTermStructure" method. Also needed for Cliquet. 
       "vols" is an array of length numNewBMDates */
    DoubleArray* fwdVolsForFSWithTermStructure(
        double                   interpLevel,
        double                   spot,
        bool                     allowNegativeFwdVar,
        const                    DateTime& startDate,
        int                      numNewBMDates,
        const DateTimeArray&     bmDates){

        static const string method = 
            "MatrixVolInterp::fwdVolsForFSWithTermStructure";

        DoubleArraySP vols(new DoubleArray(numNewBMDates-1));
        
        // calculate ATM vols on new dates 
        surfInterpVol(spot,
                      allowNegativeFwdVar,
                      bmDates,
                      true,
                      *vols);

        // insert the first vol value i.e baseDate to startDate
        DateTimeArray baseToStart(2);
        DoubleArray baseToStartVol(1);
        baseToStart[0] = volSurf->baseDate;
        baseToStart[1] = startDate;
        surfInterpVol(spot,
                      allowNegativeFwdVar,
                      baseToStart,
                      true,
                      baseToStartVol);

        /* insert the base to start date vol at the front
           NOTE: if any iterator code is added after this point
           you must reserve extra memory on the vols doubleArray
           to ensure the iterators remain valid */
        vector<double>::iterator idx = vols->begin();
        vols->insert(idx, baseToStartVol[0]);

        // interpolate at strike
        TInterp interpStrike,interpATM; 
        strikeInterpLinearDouble(volSurf->strikes,
                                 spot * interpLevel,
                                 allowNegativeFwdVar,
                                 &interpStrike);

        // interpolate atm
        strikeInterpLinearDouble(volSurf->strikes,
                                 spot,
                                 allowNegativeFwdVar,
                                 &interpATM);
        
        /* Note the assumption that surface->vol[] provides vols at the same 
           **intervals** as bmDates[] */
        double volAtStrike, volATM; 
        const CDoubleMatrix& vol = *volSurf->vol;
        int i = 0;
        for (i = 0; i < numNewBMDates-1; i++)
        {
            //calculate spot smile
            volAtStrike = volAtStrikeExpiry(spot*interpLevel,
                                            i,
                                            endDate);
            //volAtStrike = vol[interpStrike.lower][i] +
            //    interpStrike.frac * (vol[interpStrike.upper][i] - 
            //                         vol[interpStrike.lower][i]);
            volATM = volAtStrikeExpiry(spot,
                                       i,
                                       endDate);
            //volATM = vol[interpATM.lower][i] + interpATM.frac *
            //    (vol[interpATM.upper][i] - vol[interpATM.lower][i]);
            
            // add smile to fwd vol
            (*vols)[i+1] += volAtStrike - volATM;
            if (i == 0)
            {
                (*vols)[0] += volAtStrike - volATM; 
            }
        }

        DoubleArray years(numNewBMDates);
        // compute years from start to each benchmark
        for (i = 0; i < numNewBMDates; i++) {
            years[i] = volSurf->metric->yearFrac(startDate, bmDates[i]);
        }

        // have to ensure neither -ve vol nor -ve variance
        for (i = 0; i < numNewBMDates; i++)
        {
            // floor it
            (*vols)[i] = Maths::max((*vols)[i], 
                                    VolSurface::FWD_START_MIN_FWD_VOL);

            // now deal with -ve variance
            if (i > 1) {
				// set the minimum fwd variance to 2 bps
				double volFloor = ((*vols)[i-1]*VolSurface::FWD_START_MIN_FWD_VARIANCE)*sqrt(years[i-1]/years[i]);

                (*vols)[i] = Maths::max((*vols)[i], volFloor);
            }
        } 

        return vols.release();
    }

    VolSurfaceSP interpolateFSWithTermStructure(
        double            interpLevel,
        double            spot,
        bool              allowNegativeFwdVar,
        const DateTime&   startDate){
        
        int numNewBMDates = 0;

        ExpiryArraySP bmExpiries(bmExpiriesForFSWithTermStructure(
            startDate,
            volSurf->dates->back(),
            numNewBMDates));
        
        // convert the expiries to dates for the vol calculations
        DateTimeArraySP bmDatesSP(new DateTimeArray(numNewBMDates));
        DateTimeArray&  bmDates = *bmDatesSP;
        for (int i = 0; i < bmExpiries->size(); i++)
        {
            // the DateTime parameter passed to 'toDate' has no affect
            bmDates[i] = (*bmExpiries)[i]->toDate(startDate);
        }

        DoubleArraySP vols(fwdVolsForFSWithTermStructure(interpLevel,
                                                         spot,
                                                         allowNegativeFwdVar,
                                                         startDate,
                                                         numNewBMDates,
                                                         bmDates));
        
        /* calculate variance between base and start date - need calendar time
           from base date to start date */
        double yearsToStart = volSurf->metric->yearFrac(volSurf->baseDate,
                                                        startDate);
        // save yearFrac's for use in VolSurface constructor
        DoubleArraySP myTradYears(new DoubleArray(bmDates.size()));
        (*myTradYears)[0] = yearsToStart; // first expiry = startDate
        // build standard spot vols so can call constructor 
        double varToStart = (*vols)[0] * (*vols)[0] * yearsToStart;
        int idx;
        for (idx = 1; idx < numNewBMDates; idx++)
        {
            double yearFrac = volSurf->metric->
                yearFrac(bmDates[idx-1], bmDates[idx]);
            (*myTradYears)[idx] = (*myTradYears)[idx-1] + yearFrac;
            double yearFracFromStart = (*myTradYears)[idx] - yearsToStart;
            (*vols)[idx] =  sqrt((varToStart + (*vols)[idx] * 
                                  (*vols)[idx] * yearFracFromStart)/
                                 (yearFracFromStart + yearsToStart));
        }
            
        /* spec says last date of new benchmark dates is the last original 
           benchmark date, so cut if necessary. Note that numNewBMDates>=2 by
           construction */   
        idx = numNewBMDates-1; /* reuse for brevity&clarity */
        DateTime lastSurfaceDate = volSurf->dates->back();
        if (bmDates[idx].isGreater(lastSurfaceDate)) {
            double  t1 = volSurf->metric->yearFrac(volSurf->baseDate,
                                                   bmDates[idx-1]);
            double  var1 = (*vols)[idx-1] * (*vols)[idx-1] * t1;
            double  t2 = volSurf->metric->yearFrac(volSurf->baseDate, 
                                                   bmDates[idx]);
            double  var2 = (*vols)[idx] * (*vols)[idx] * t2;
            double  t = volSurf->metric->
                yearFrac(volSurf->baseDate,
                         lastSurfaceDate);
            
            // linearly interpolate variance then extract vol 
            (*vols)[idx] =  sqrt((var1*(t2-t)+var2*(t-t1))/(t*(t2-t1))); 
            bmDates[idx] = lastSurfaceDate;
        }
        
        // use special constructor for performance
        VolSurfaceSP fwdVolSurface(new VolSurface(volSurf.get(),
                                                  interpLevel,
                                                  *vols,
                                                  bmExpiries,
                                                  bmDatesSP,
                                                  myTradYears));
        
        return fwdVolSurface;
    }

    /* Calculates the 'forward' vol between two dates.
       calculates (sigma(t2)^2 * t2 - sigma(t1)^2 * t1)/(t2-t1) where
       sigma is the interpolated value on the surface and t1 and t2 are
       adjacent points on the time on the axis. t1 and t2 represent
       trading time (calculated using EdgTimeMetricYearFrac). The forward
       vol is given by the sqrt of the above formula.
       If timeIdx is the last point on the axis then the routine works as 
       if timeIdx-1 had been used - but yearFracToPt remains the same and 
       yearFracDiff is undefined */
    double surfCalcSquareDiff(
        TInterp*          interp,       /* (I) interp data for strike axis */
        int               timeIdx) const {      /* (I) index for t1 */

        const DateTimeArray&   dates     = *volSurf->dates;
        const DoubleArray&     tradYears = *volSurf->tradYears;
        const DoubleMatrix&    vol       = *volSurf->vol;
        double           loVol;
        double           hiVol;
        /* this means we extrapolate using flat fwd vols */
        int              i = timeIdx == dates.size()-1? timeIdx-1: timeIdx;
        double           ttYearFracToPt; /* trading time year frac from 
                                            origin to timeIdx */
        double           ttYearFracDiff;
    
        /* get time between benchmarks (from precomputed cache) - guaranteed to
           be > 0 */
        ttYearFracDiff = tradYears[i+1] - tradYears[i];

        ttYearFracToPt = tradYears[timeIdx];
        if (timeIdx == dates.size() - 1)
        {
            /* need to adjust yearFracToPt in this situation */
            ttYearFracToPt -= ttYearFracDiff;
        }

        loVol = vol[interp->lower][i] + interp->frac *
            (vol[interp->upper][i] - vol[interp->lower][i]);

        hiVol = vol[interp->lower][i+1] + interp->frac *
            (vol[interp->upper][i+1] - vol[interp->lower][i+1]);
        
        double output = (hiVol * hiVol * (ttYearFracToPt + ttYearFracDiff) -
                         loVol * loVol * ttYearFracToPt)/(ttYearFracDiff);

        if (Maths::isZero(output)){
            output = 0.0;
        }

        if (!interp->allowNegativeFwdVar && Maths::isNegative(output)){
                string m = "Negative variance ("+Format::toString(output)+")"+
                    " found for strike "+Format::toString(interp->strike)+
                    "\nbetween "+dates[i].toString()+" and "+
                    dates[i+1].toString()+" on " + volSurf->getName() + 
                    " volatility surface";
                throw ModelException("MatrixInterpVol::surfCalcSquareDiff", m);

        }
        return output;
    }

    /* calculate the variance between origin and time point index (ie
       integral of square of surface between origin and T). Interpolation
       between time points is via constant fwd vols. Routine works if
       nextTimeIdx = num Dates on axis. If nextTimeIdx = 0 then varToIdx = 0.
       Does not support time axis with single point */
    void surfIntegrateSquare(
        TInterp*          interp,     /* interp data for strike axis */
        double            level,
        int               nextTimeIdx,/* (I) idx of next time point 
                                         on/after timept */
        double           *fwdVolSq,   /* (O) square of fwd vol in timePt
                                         interval */
        double           *varToIdx) const {   /* (O) variance from origin to 
                                          nextTimeIdx-1*/

        const DoubleArray&     tradYears = *volSurf->tradYears;
        const DoubleMatrix&    vol       = *volSurf->vol;
        double interpVol;

        if (nextTimeIdx == 0) {
            /* timePt is before first point on axis */

            interpVol = volAtStrikeExpiry(level, 0, endDate);
            *fwdVolSq = interpVol * interpVol;
            *varToIdx = 0.0;
        } else {
            int    timeIdx = nextTimeIdx - 1;
            
            /* need fwd vol from previous time point */
            *fwdVolSq = surfCalcSquareDiff(interp, timeIdx);
            /* can calculate variance directly */
            interpVol = volAtStrikeExpiry(level, timeIdx, endDate);
            *varToIdx = interpVol * interpVol * tradYears[timeIdx];
        }
        return;
    }

    void surfMoveToNextInterval(
        TInterp*          interp,    /* (I) */
        int               newIdx,    /* (I) where we are now */
        double           *varSum,    /* (O) variance to idx */
        double           *fwdVolSq) const { /* (O) fwd vol squared for new interval */

        const DoubleArray&  ttYears = *volSurf->tradYears;
        const DoubleMatrix& vol     = *volSurf->vol;
        double      interpVol;
    
        /* refresh variance up to newIdx */
        interpVol = vol[interp->lower][newIdx] + interp->frac *
           (vol[interp->upper][newIdx] - vol[interp->lower][newIdx]);
        *varSum = interpVol * interpVol * ttYears[newIdx];
    
        *fwdVolSq = surfCalcSquareDiff(interp, newIdx);
    }

    // Obtains vol(strike K, expiryIdx t) by interpolating at a function of strike K
    // at optionMat T specified by stri
    // This is achieved by linearly interpolating in strike for inv_f(t, f(T, K))
    // An alternative method would be to map each strike at expiryIdx and linearly interpolate
    // in mapped strike. The two approaches are equivalent if f() is a linear function of
    // strike.
    double volAtStrikeExpiry(double strike, int expiryIdx, const DateTime& optionMat) const {
        const DoubleMatrix& vol = *volSurf->vol;
        const DoubleArray& strikes = volSurf->getStrikes();
        double mappedStrike = strikeMap->fromStrike(strike, optionMat);     // f(T, K)
        double targetStrike = strikeMap->toStrike(mappedStrike, (volSurf->getDates())[expiryIdx]); // finv(t, f(T, K))
        int lowerIdx, upperIdx;
        double interpFrac;
        InterpList::getLinearInterpFrac(targetStrike,
                                        strikes,
                                        lowerIdx,
                                        upperIdx,
                                        interpFrac);
        return vol[lowerIdx][expiryIdx] + interpFrac *
            (vol[upperIdx][expiryIdx] - vol[lowerIdx][expiryIdx]);
    }

    double volAtExpiry(int idx, TInterp* interp) const {
        const DoubleMatrix& vol = *volSurf->vol;
        double volAtExpiry = vol[interp->lower][idx] + interp->frac *
                (vol[interp->upper][idx] - vol[interp->lower][idx]);
        return volAtExpiry;
    }

    /* only one date on time axis - so simple interpolation */
    void simpleInterpVol(
        TInterp*             interp,     /* (I) */
        double               level,      /* (I) */
        bool                 calcVols,   /* (I) */
        const DateTimeArray& dates,      /* (I) */
        double               spread,     /* (I) */
        DoubleArray&         output) const {     /* (O) */
        static string     routine = "Interp::simpleInterpVol";
        try{
            const DoubleMatrix& vol = *volSurf->vol;
            int              i;
            DateTime         loDate;
            const DateTime&  dateFrom = dates[0];
            int              numDatesTo = dates.size()-1;
            /* only one date on time axis - so simple interpolation */
            double interpVol = volAtStrikeExpiry(level, 0, endDate);

            if (spread < -interpVol)
            {
                string m = "Vol spread ("+Format::toString(spread)+") will "
                    "make vol ("+Format::toString(interpVol)+") negative";
                throw ModelException(routine, m);
            }

            interpVol += spread;

            if (calcVols) {
                /* Don't need to back out the vols from the var since the vol
                   is just the interpVol in this case */
                for (i = 0; i < numDatesTo; i++)
                {
                    output[i] = interpVol;
                }
            } else {
                interpVol *= interpVol;

                for (i = 0, loDate = dateFrom; 
                     i < numDatesTo; loDate = dates[i+1], i++)
                {
                    double years = volSurf->metric->yearFrac(loDate,
                                                             dates[i+1]);
                    /* Calculate variance */
                    output[i] = interpVol * years;
                }
            }
    
        } catch (exception& e){
            throw ModelException(&e, routine);
        }
    }


#if defined(_MSC_VER)
#pragma optimize("", on)
#endif

    /** Does surfIntegrateSquare and then adds supplied spread to
        fwdVolSq - checks for negative spread etc */
    void surfIntegrateSquareWithSpread(
        double            spread,
        TInterp*          interp,   
        double            level,
        int               timePtIdx,
        double&           fwdVolSq,
        double&           var) const { 

        /* add variance after last benchmark up to dateTo */
        surfIntegrateSquare(interp, 
                            level,
                            timePtIdx, 
                            &fwdVolSq, 
                            &var);
        if (Maths::isZero(fwdVolSq)){
            fwdVolSq = 0.0;
        }
        if (Maths::isZero(spread)){
            // do nothing
        } else if (Maths::isNegative(fwdVolSq)){
            // cunning: switch flag to disallow negative variance and redo
            interp->allowNegativeFwdVar = false;
            surfIntegrateSquare(interp, 
                                level,
                                timePtIdx, 
                                &fwdVolSq, 
                                &var);    
            // should never get here
            throw ModelException("addSpread", "internal error");
        } else {
            /* add spread to fwd vol */
            fwdVolSq = sqrt(fwdVolSq);
            fwdVolSq += spread;
            fwdVolSq *= fwdVolSq;
        }
    }
        
    /* Uses the level to determine position on strike axis.  On date
       axis, from datesTo[0] to each of datesTo[i>0], either the
       'average' value is returned (calcVols = TRUE) or the 'integrated'
       value is returned (calcVols = FALSE).  The convention is that the
       points in the surface (for a fixed strike) are constant going
       forward in time. The 'integrated' value is the value squared * time
       interval (in years). The 'average' value is calculated by backing
       out from the integrated value.  This routine has no concept of the
       idea of forward starting interpolation All vols in surface are spot
       vols - this routine does not work for a fwd vol surface. */
public:
    // as we do this all the time
    static void throwNegVariance(const string&   method,
                                 double          var,
                                 const DateTime& t1,
                                 const DateTime& t2,
                                 double          k) {
        throw ModelException(method,
                             "Variance ("+ Format::toString(var)+ ") " 
                             "is <= zero between " + t1.toString()+
                             " and " + t2.toString() + " at strike " +
                             Format::toString(k));
    }
                                 
private:

    void surfInterpVol(
        double               level,      /* (I) where to interpolate on 
                                            strike axis */
        bool                 allowNegativeFwdVar, // (I) don't fail if var<0
        const DateTimeArray& datesTo,    /* (I) measure to these points */
        bool                 calcVols,   /* (I) TRUE = calc vols,
                                            else calc variance */
        DoubleArray&         output) const {     /* (O) vols or variance */

        static string   routine = "surfInterpVol";

        try{
            volSurf->checkCache();
            const DateTimeArray&  dates   = *volSurf->dates;
            const DoubleArray&    strikes = volSurf->strikes;
            const DoubleArray&    years   = *volSurf->tradYears;
            TInterp         interp  = {0};
            const DateTime& dateFrom = datesTo[0];
            int             numDatesTo = datesTo.size()-1;

            if (volSurf->baseDate.isGreater(dateFrom)) {
                string m = "From date ("+dateFrom.toString()+") is before "
                    "origin of vol surface ("+volSurf->baseDate.toString()+")";
                throw ModelException(routine, m);
            }

            if (numDatesTo > 0 && dateFrom.isGreater(datesTo[1])) {
                string m = "First date ("+datesTo[1].toString()+
                    ") is before date from ("+dateFrom.toString()+")";
                throw ModelException(routine, m);
            }

            /* interpolate on strike axis using level */
            strikeInterpLinearDouble(strikes, level, allowNegativeFwdVar, 
                                     &interp);

            if (dates.size() == 1) {
                simpleInterpVol(&interp, 
                                level,
                                calcVols,
                                datesTo, 
                                0.0,       /* no spread */
                                output);
            } else {
                double varDateFrom;
                int    i;
                double fwdVolSq;
                double varSum;
                double yearFracTo;
                double yearFracFrom;
                double yearDateFrom;
                int    timePtIdx;
                int    lastIdx;
                int    dateFromIdx;
                double varSumSoFar;
                bool   withinSameBMint;

                /* get variance to dateFrom - first need to find place on 
                   date axis */
                for (timePtIdx = 0; timePtIdx < dates.size() &&
                         !dates[timePtIdx].isGreater(dateFrom);
                     timePtIdx++)
                {
                    ; /* empty loop */
                }
                dateFromIdx = timePtIdx - 1;  

            /* then get data for corresponding time index */
                surfIntegrateSquare(&interp, 
                                    level,
                                    timePtIdx, 
                                    &fwdVolSq, 
                                    &varSum);
        
                /* calc time from timePtIdx-1 to dateFrom */
                if (dateFrom.equals(volSurf->baseDate)) {
                    yearFracFrom = 0.0;
                    varDateFrom  = 0.0;
                    yearDateFrom = 0.0;
                } else {
                    yearFracFrom = volSurf->metric->yearFrac(
                        timePtIdx==0? volSurf->baseDate: dates[timePtIdx-1],
                        dateFrom);

                    varDateFrom  = varSum + fwdVolSq * yearFracFrom;
                    yearDateFrom = (timePtIdx == 0 ? 0.0: 
                                    years[timePtIdx-1]) + yearFracFrom;       
                }

                /* now loop through dateTimes */
                for (i = 0, varSumSoFar = 0.0; i < numDatesTo; i++) {
                    withinSameBMint = false;
                    if (i > 0 && 
                        datesTo[i].isGreater(datesTo[i+1])) {
                        throw ModelException(routine, 
                                             "DateTimes are not increasing\n");
                    }
            
                    /* starting at dateFromIdx, search for index for date 
                       in datesTo */
                    for (lastIdx = timePtIdx; /* timePtIdx has value already */
                         timePtIdx < dates.size() &&
                             !dates[timePtIdx].isGreater(datesTo[i+1]);
                         timePtIdx++) {
                        ; /* empty loop */
                    }
            
                    /* Check if 'from' and 'to' dates are within same
                       BM interval */
                    if (timePtIdx - 1  == dateFromIdx) {
                        withinSameBMint = true;
                    }

                    /* refresh cached values for this period */
                    if (timePtIdx != lastIdx) {
                        surfMoveToNextInterval(&interp,
                                               timePtIdx-1, 
                                               &varSum, 
                                               &fwdVolSq);
                        lastIdx = timePtIdx;
                    }
                    
                    /* Calculate yearFrac1 from nextTimePtIdx-1 to date */
                    yearFracTo = volSurf->metric->yearFrac(
                        timePtIdx==0? volSurf->baseDate: dates[timePtIdx-1],
                        datesTo[i+1]);

                    /* Before calculating vols, first Calculate variance using
                       varSum and interp vol at that pt plus yearFrac1 and
                       fwdVolSq */
                    // volatile to ensure numerical consistency
                    double volatile tmp = varSum + fwdVolSq * yearFracTo;
                    output[i] = tmp - varDateFrom;
            
                    /* do we need to check for negative variance here?  In
                       theory have only added positive numbers together so
                       should be okay.  However, on intel with fp
                       registers of 10 bytes can get a tiny negative
                       number here */
                    if (Maths::isZero(output[i])) {
                        output[i] = 0.0;
                    } else if (Maths::isNegative(output[i])) {
                            throwNegVariance(routine,
                                             output[i],
                                             dateFrom,
                                             datesTo[i+1],
                                             level);
                    }

                    if (calcVols) {
                        /* vols are for trading time so divide variance by
                           trading time i.e timeMetricYearFrac. Trad time
                           between dateFrom and dateTo can be calculated
                           using the tradYears cache. */
                        double tradYearFrac;  
                                     
                        if (withinSameBMint) {
                            tradYearFrac = yearFracTo - yearFracFrom;
                        } else {
                            tradYearFrac = (years[timePtIdx-1] - 
                                            (dateFromIdx<0? 0: 
                                             years[dateFromIdx]))-
                                yearFracFrom + yearFracTo;
                        }

                        /* here we need to check for zero time interval */
                        if (Maths::isZero(tradYearFrac)){
                            // be generous and just use the current vol
                            output[i] = sqrt(fwdVolSq);
                        } else if (tradYearFrac < DBL_EPSILON) {
                            string m = "Year fraction ("+
                                Format::toString(tradYearFrac)+ ") " 
                                "is < zero\nbetween "+dateFrom.toString()+
                                " and "+datesTo[i+1].toString();
                            throw ModelException(routine, m);
                        } else {
                            output[i] = sqrt(output[i]/tradYearFrac);
                        }
                    } else {
                        // for variances, we're returning variance per interval
                        double temp = varSumSoFar;
                        varSumSoFar = output[i];
                        output[i] -= temp;

                        // floor it to zero as we seem to get numerical error. 
                        // Motivated from rollTheta12.xml !!!
                        if(!Maths::isPositive(output[i]) && 
                           Maths::isPositive(output[i] + 10.0 * DBL_EPSILON)){
                            output[i] = 0.0;
                        }

                        // output[i] is now a Fwd Variance
                        else if (Maths::isNegative(output[i]) && 
                                 !allowNegativeFwdVar) {
                           throwNegVariance(routine,
                                             output[i],
                                             datesTo[i],
                                             datesTo[i+1],
                                             level);  
                        }
                    }
                }
            }
        } catch (exception& e){
            throw ModelException(&e, routine);
        }
    }

    void fwdStartInterpVol(
        const DateTimeArray&   datesTo,    /* (I) measure to these points */
        bool                   calcVols,   /* (I) TRUE = calc vols,
                                              else calc variance */
        bool                   allowNegativeFwdVar, /* (I) */
        DoubleArray&           output) const {    /* (O) vols or variance */

        static string    routine = "fwdStartInterpVol";

        volSurf->checkCache(); // must be done as very first thing
        const DateTimeArray&    dates   = *volSurf->dates;
        const DoubleArray&      strikes = volSurf->strikes;
        const DoubleArray&      years   = *volSurf->tradYears;
        TInterp           interp  = {0};
        double            spread;
        const DateTime&   dateFrom = datesTo[0];
        int               numDatesTo = datesTo.size()-1;


        if (volSurf->baseDate.isGreater(dateFrom)) {
            string m = "From date ("+dateFrom.toString()+") is before "
                "origin of vol surface ("+volSurf->baseDate.toString()+")";
            throw ModelException(routine, m);
        }

        if (numDatesTo > 0 && dateFrom.isGreater(datesTo[1])) {
            string m = "First date ("+datesTo[1].toString()+
                ") is before date from ("+dateFrom.toString()+")";
            throw ModelException(routine, m);
        }

        /* get spread to apply to spot level forward vols */
        spread = fwdStartVolSpread();

        /* interpolate on strike axis using spotlevel */
        strikeInterpLinearDouble(strikes, spot, allowNegativeFwdVar,
                                 &interp);

        if (dates.size() == 1) {
            simpleInterpVol(&interp, 
                            spot,
                            calcVols,
                            datesTo, 
                            spread,
                            output);
        }
        else
        {
            double varDateFrom;
            int    i;
            double fwdVolSq;
            double var;        /* value not used */
            double varSum;
            double yearFracTo;
            double yearFracFrom;
            double yearDateFrom;
            double yearsToPrevious;
            int    timePtIdx;
            int    dateFromIdx;
            double varSumSoFar;
            bool   withinSameBMint;

            /* get variance to dateFrom - first need to find place on 
               date axis */
            timePtIdx = 0;
            varSum    = 0.0;
            yearsToPrevious = 0.0;

            /* spread is applied to forward vol, so have to build up variance
             * to dateFrom benchmark by benchmark - can't apply spread to
             * 'average' variance
             */

            while (timePtIdx < dates.size() && 
                   !dates[timePtIdx].isGreater(dateFrom))
            {
                /* get fwd vol for corresponding time index */
                surfIntegrateSquareWithSpread(spread,
                                              &interp, 
                                              spot,
                                              timePtIdx, 
                                              fwdVolSq, 
                                              var);
                
                varSum += fwdVolSq *(years[timePtIdx]- yearsToPrevious);
            
                /* remember current benchmark */
                yearsToPrevious = years[timePtIdx];
                timePtIdx++;
            }
            dateFromIdx = timePtIdx - 1;

            /* get variance up to dateFrom from the last benchmark prior to 
             * dateFrom (or from base date if dateFrom lies inside first
             * benchmark)
             */
            surfIntegrateSquareWithSpread(spread,
                                          &interp, 
                                          spot,
                                          timePtIdx, 
                                          fwdVolSq, 
                                          var);

            /* calc time from timePtIdx-1 to dateFrom */
            if (dateFrom.equals(volSurf->baseDate)) {
                yearFracFrom = 0.0;
                varDateFrom  = 0.0;
                yearDateFrom = 0.0;
            } else {
                yearFracFrom = volSurf->metric->yearFrac(
                    timePtIdx==0? volSurf->baseDate: dates[timePtIdx-1],
                    dateFrom);
                varDateFrom  = varSum + fwdVolSq * yearFracFrom;
                yearDateFrom = (timePtIdx == 0 ? 0.0: years[timePtIdx-1]) + 
                    yearFracFrom;           
            }

            /* now loop through dateTimes */
            for (i = 0, varSumSoFar = 0.0; i < numDatesTo; i++)
            {
                withinSameBMint = false;
                if (i > 0 && 
                    datesTo[i].isGreater(datesTo[i+1]))
                {
                    throw ModelException(routine,
                                         "DateTimes are not increasing");
                }
            
                /* this is poor performance wise - as we recalculating the
                   variance from timePtIdx = 0 to the date in question each
                   time through the loop cf surfInterpVol */
                timePtIdx = 0;
                varSum = 0.0;
                yearsToPrevious = 0.0;

                /* in a similar way, accumulate variance up to each date in
                 * datesTo by getting fwd vol between benchmarks and adding
                 * the spread
                 */

                while (timePtIdx < dates.size() &&
                       !dates[timePtIdx].isGreater(datesTo[i+1]))
                {
                    /* get fwd vol for corresponding time index */
                    surfIntegrateSquareWithSpread(spread,
                                                  &interp, 
                                                  spot,
                                                  timePtIdx, 
                                                  fwdVolSq, 
                                                  var);
                    varSum += fwdVolSq * (years[timePtIdx] - yearsToPrevious);

                    yearsToPrevious = years[timePtIdx];
                    timePtIdx++;
                }

                /* Check if 'from' and 'to' dates are within same
                   BM interval */
                if (timePtIdx - 1  == dateFromIdx) {
                    withinSameBMint = true;
                }

                /* add variance after last benchmark up to dateTo */
                surfIntegrateSquareWithSpread(spread,
                                              &interp, 
                                              spot,
                                              timePtIdx, 
                                              fwdVolSq, 
                                              var);

                /* Calculate yearFrac1 from nextTimePtIdx-1 to date */
                yearFracTo = volSurf->metric->yearFrac(
                    timePtIdx==0? volSurf->baseDate: dates[timePtIdx-1],
                    datesTo[i+1]);

                /* Before calculating vols, first calculate variance using
                   varSum and interp vol at that pt plus yearFrac1 and
                   fwdVolSq */
                // volatile to ensure numerical consistency
                double volatile tmp = varSum + fwdVolSq * yearFracTo;
                output[i] = tmp - varDateFrom;

                /* do we need to check for negative variance here?  In
                   theory have only added positive numbers together so
                   should be okay.  However, on intel with fp
                   registers of 10 bytes can get a tiny negative
                   number here */
                if (Maths::isZero(output[i])) {
                    output[i] = 0.0;
                } else if (Maths::isNegative(output[i])) {
                    string m = "Variance ("+
                        Format::toString(output[i])+ ") " 
                        "is <= zero\nbetween "+dateFrom.toString()+
                        " and "+datesTo[i+1].toString();
                    throw ModelException(routine, m);
                }
                if (calcVols) {
                    /* vols are for trading time so divide variance by
                       trading time i.e timeMetricYearFrac. Trad time
                       between dateFrom and dateTo can be calculated
                       using the tradYears cache. */
                    double tradYearFrac;
                
                    if (withinSameBMint) {
                        tradYearFrac = yearFracTo - yearFracFrom;
                    } else  {
                        tradYearFrac = (years[timePtIdx-1] - 
                                        (dateFromIdx<0? 0: 
                                         years[dateFromIdx])) -
                            yearFracFrom + yearFracTo;
                    }

                    /* here we need to check for zero time interval */
                    if (Maths::isZero(tradYearFrac)){
                        // be generous and just use the current vol
                        output[i] = sqrt(fwdVolSq);
                    } else if (tradYearFrac < DBL_EPSILON) {
                        throw ModelException(routine, "Year fraction ("+
                                             Format::toString(tradYearFrac)+
                                             ") is < zero");
                    } else {
                        // we checked that output[i] is not -ve above
                        output[i] = sqrt(output[i]/tradYearFrac);
                    }
                } else {
                    /* for variances, we're returning variance per interval */
                    // output[i] is now Total Variance
                    double temp = varSumSoFar;
                    varSumSoFar = output[i];
                    output[i] -= temp;
                    
                    // output[i] is now Total Variance
                    /* 1st check/adjust for numerical noise. Then for -ve var*/
                    if (Maths::isZero(output[i])) {
                        output[i] = 0.0;
                    } else if (Maths::isNegative(output[i]) && 
                               !allowNegativeFwdVar) {
                        string m = "Variance ("+
                            Format::toString(output[i])+ ") " 
                            "is <= zero\nbetween "+datesTo[i].toString()+
                            " and "+datesTo[i+1].toString();
                        throw ModelException(routine, m);
                    }
                }

            }
        }
    }
private:
    /* for reflection */
    Interp(): CVolProcessedBS(TYPE), spot(0.0), 
        allowNegativeFwdVar(false), spread(0.0) {};
    
private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(Interp, clazz);
        SUPERCLASS(CVolProcessedBS);
        EMPTY_SHELL_METHOD(defaultInterp);
        FIELD(volSurf, "Vol Surface");
        FIELD(fwdStartInterp, "is it forward starting interpolation");
        FIELD(interpLevel, "interpolation level");
        FIELD(spot, "spot price");
        FIELD(startDate, "interpolation start date");
        FIELD(endDate, "interpolation end date");
        FIELD(isManualSpread, "use manual spread");
        FIELD(spread, "manual spread");
        FIELD(allowNegativeFwdVar, "Allow negative forward variance");
        FIELD_MAKE_OPTIONAL(allowNegativeFwdVar); // default false
        FIELD(strikeMap, "Strike Map");
        FIELD_MAKE_OPTIONAL(strikeMap);
    }

    static IObject* defaultInterp(){
        return new Interp();
    }

};

CClassConstSP const VolSurface::Interp::TYPE = 
CClass::registerClassLoadMethod(
    "VolSurface::Interp", typeid(VolSurface::Interp), load);

/** minimum fwd vol for fwd starting interpolation */
const double VolSurface::FWD_START_MIN_FWD_VOL = 0.02;


/** minimum fwd variance for fwd starting interpolation, set to 2 bps */
const double VolSurface::FWD_START_MIN_FWD_VARIANCE = 1.0002;

/** Returns name of vol */
string VolSurface::getName() const{
    return name;
}


/** this gets called after an object is constructed from a data dictionary.
    Not after an object has been copied (see override of clone method below) */
void VolSurface::validatePop2Object(){
    static const char routine[] = "VolSurface::validatePop2Object";
    try {
        /* validate that the number of strikes and number of dates matches up
           with double matrix */
        int numStrikes = strikes.size();
        int numDates = expiries->size();

        if (numStrikes < 1 || numDates < 1){
            throw ModelException(routine, getName()+" - at least one strike and "
                                 "benchmark date must be supplied");
        }
        if (strikes.size() != vol->numCols()) {
            throw ModelException(routine, getName()+
                                 " - mismatch between number of strikes ("+
                                 Format::toString(strikes.size())+") and\n"
                                 "number of columns in matrix ("+
                                 Format::toString(vol->numCols())+")");
        }
        if (expiries->size() != vol->numRows()) {
            throw ModelException(routine, getName()+
                                 " - mismatch between number of benchmark dates ("+
                                 Format::toString(expiries->size())+") and\n"
                                 "number of rows in matrix (" + 
                                 Format::toString(vol->numRows()) + ")");
        }
        /* check that strikes are > 0 && strictly increasing - although allow
           through curves with zero strike (repercussions for vega skew though) */
        if (numStrikes > 1 && !Maths::isPositive(strikes[0])){
            throw ModelException(routine, 
                                 getName()+ " - first strike is not positive");
        }
        for (int idx = 1; idx < numStrikes; idx++){
            if (!Maths::isPositive(strikes[idx] - strikes[idx-1])) {
                throw ModelException(routine, getName()+ " - "
                                     "Strikes "+Format::toString(strikes[idx-1])+
                                     " and "+Format::toString(strikes[idx])+
                                     " need to be strictly increasing");
            }
        }

        // flag up if we have any absolute benchmarks (for time rolling)
        for (int i = 0; i < expiries->size(); i++) {
            bool isAbsolute = BenchmarkDate::TYPE->isInstance((*expiries)[i].get()); 

            if (i > 0 && isAbsolute) {
                // must have fixed BM before rolling ones (to prevent overlap)
                if (!BenchmarkDate::TYPE->isInstance((*expiries)[i-1].get())) {
                    throw ModelException(routine, "must have rolling  benchmarks "
                                         "after fixed ones - got " +
                                         (*expiries)[i]->toString() + " after " + 
                                         (*expiries)[i-1]->toString());
                }
            }                                
        }
        zapPastBM(false);  // remove historic fixed benchmarks
        numDates = expiries->size(); // reset incase we removed any
        for (int iStrike = 0; iStrike < numStrikes; iStrike++){
            for (int iExpiry = 0; iExpiry < numDates; iExpiry++){
                if ((*vol)[iStrike][iExpiry] < -DBL_EPSILON){
                    throw ModelException(routine, 
                                         "the volatility with strike "
                                         + Format::toString(strikes[iStrike])
                                         + " and expiry "
                                         + (*expiries)[iExpiry]->toString() 
                                         + " ("
                                         + Format::toString((*vol)[iStrike][iExpiry])
                                         + ") is negative");
                                     
                }
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, routine, "validation failed for " + getName());
    }
}

/* build up a cache of year fractions from surface 'base date' to
 * each benchmark date  */
void VolSurface::buildYearFracCache() const{   /* (M) */
    DoubleArraySP   myTradYears(new DoubleArray(expiries->size()));
    DateTimeArraySP myDates(new DateTimeArray(expiries->size()));
    for (int i = 0; i < myDates->size(); i++) {
        //calc benchmark date from expiry
        (*myDates)[i] = (*expiries)[i]->toDate(baseDate);
    }
    // calculate yearFractions
    metric->yearFrac(baseDate, *myDates, *myTradYears);
    // validate them
    for (int j = 0; j < myTradYears->size(); j++){
        double yearFracSoFar = j == 0? 0.0: (*myTradYears)[j-1];
        double yearFracDiff = (*myTradYears)[j] - yearFracSoFar;
        if (yearFracDiff < DBL_EPSILON) {
            const DateTime& startDate = j == 0? baseDate: (*myDates)[j-1];
            string m = "Trading time ("+Format::toString(yearFracDiff)+
                ") between benchmark dates ("+startDate.toString()+") and ("+
                (*myDates)[j].toString()+") is <= 0";
            throw ModelException("buildYearFracCache", m);
        }
    }            
    tradYears = myTradYears;
    dates = myDates;
    gotCache = true;
}

//// wrapper to checkNonNegative - also avoids strange crash on solaris
static void matrixHelper(const CDoubleMatrix& matrix, double shiftSize,
                         const string&        name,   const string& routine){
    try{
        matrix.checkNonNegative();
    } catch (exception& e){
        string m = "Vol shift ("+ Format::toString(shiftSize)+
            ") results in negative vol for " + name;
        throw ModelException(e, routine, m); 
    }
}

/** Returns name identifying vol for vega parallel */
string VolSurface::sensName(const VolParallel*) const {
    return getName();
}

/** Shifts the object using given shift */
TweakOutcome VolSurface::sensShift(const PropertyTweak<VolParallel>& tweak){
    static const string routine = "VolSurface::sensShift<VolParallel>";
    if (!Maths::isZero(tweak.coefficient)){
        vol->scalarAdd(tweak.coefficient);
        if (tweak.coefficient < 0){
            try{
                matrixHelper(*vol, tweak.coefficient, getName(), routine);
            } catch (exception& e){
                // if we fail the restore method is not called
                sensRestore(tweak);
                throw ModelException(e, routine);
            }
        }
    }
    return TweakOutcome(tweak.coefficient,
                        false); // none of our components has a vega type sensitivity
}

/** Restores the object to its original form */
void VolSurface::sensRestore(const PropertyTweak<VolParallel>& tweak){
    if (!Maths::isZero(tweak.coefficient)){
        vol->scalarAdd(-tweak.coefficient);
    }
}

/** Returns name identifying vol for vega pointwise */
string VolSurface::sensName(const VolPointwise*) const{
    return getName();
}

/** Returns the array of expiries (ie maturities/benchmark dates) that
    need to be tweaked for this vol */
ExpiryWindowArrayConstSP VolSurface::sensQualifiers(const VolPointwise*) const{
    return ExpiryWindow::series(expiries);
}

/** Shifts the object using given shift */
TweakOutcome VolSurface::sensShift(const PropertyTweak<VolPointwise>& shift){
    static const string routine = "VolSurface::sensShift";
    try{
        if (!Maths::isZero(shift.coefficient)){
            int expiryIdx = shift.qualifier->expiry->search(expiries.get());
            vol->rowAdd(expiryIdx, shift.coefficient);
            if (shift.coefficient < 0){
                matrixHelper(*vol, shift.coefficient, getName(), routine);
            }
        }
        return TweakOutcome(shift.coefficient,
                            false); // none of our components has a vega type sensitivity
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

/** Returns name identifying vol for vega pointwise */
string VolSurface::sensName(VegaMatrix* shift) const{
    return getName();
}

/** Returns the array of expiries (ie maturities/benchmark dates) that
    need to be tweaked for this vol */
ExpiryArrayConstSP VolSurface::sensExpiries(VegaMatrix* shift) const{
    return expiries;
}

/** Shifts the object using given shift */
bool VolSurface::sensShift(VegaMatrix* shift) {
    static const string routine = "VolSurface::sensShift";


    /* retrieve shift information from MatrixShift */
    const ExpiryArray&  shiftExpiries = *(shift->getExpiries());
    const DoubleArray&  shiftAxis     = *(shift->getXAxisValues());
    int                 expiryIdx     = shift->getExpiryIdx();

    /* build a vol surface with just the sensitive strikes in it */
    CDoubleMatrixSP newVol(new DoubleMatrix(shiftAxis.size(), 
                                            shiftExpiries.size()));
    DoubleMatrix& volMatrix = *newVol;

    /* check that the expiry to bump matches ours (it always should since
       we supplied the expiries */
    if (expiryIdx >= expiries->size() ||
        !shiftExpiries[expiryIdx]->equals((*expiries)[expiryIdx].get())){
        throw ModelException(routine, "Maturity "+
                             shiftExpiries[expiryIdx]->toString()+
                             " unexpected");
    }

    /* now just do linear interpolation at benchmark dates for each of the
       strikes */
    int j;
    for (j = 0; j < shiftAxis.size(); j++){
        VolSurface::Interp::TInterp   interpCache;
        VolSurface::Interp::strikeInterpLinearDouble(strikes,
                                                     shiftAxis[j],
                                                     false,
                                                     &interpCache);
        double* volsAtExpiries = volMatrix[j];
        for (int i = 0; i < expiries->size(); i++) {
            volsAtExpiries[i] = (*vol)[interpCache.lower][i] + 
                interpCache.frac *
                ((*vol)[interpCache.upper][i] - (*vol)[interpCache.lower][i]);
        }
    }

     // copy the new strikes over
    strikes = shiftAxis; // structure copy
    // set the new vols  
    vol = newVol;  

    /* tweak points given by strike & date.
       If we have adjacent strikes that are extremely close together
       then tweak all of these points at the same time defined by the lower
       and upper indices contained in the shift. Only one sensitivity will 
       be reported for the lowest strike in the series */
    double shiftSize = shift->getShiftSize();
    int    upperIdx  = shift->getUpperIdx();
    for (j = shift->getLowerIdx(); j <= upperIdx; j++) {
        volMatrix[j][expiryIdx] += shiftSize;
        if (shiftSize < 0 && Maths::isNegative(volMatrix[j][expiryIdx])) {
             throw ModelException(routine, getName() + " - shift " + 
                                  Format::toString(shiftSize) + " makes vol at "
                                  + Format::toString(shiftAxis[j]) + 
                                  "negative.");
         }            
     }

     return false; // nothing else to tweak
 }



 /** Returns name identifying vol for root time vega */
 string VolSurface::sensName(RootTimeVega* shift) const{
     return getName();
 }

 /** apply a shift to each point on the surface equal to 
     shift / SQRT(TTi)  where TTi is the trading time year fraction
     up to the benchmark date Ti */
 bool VolSurface::sensShift(RootTimeVega* shift) {
     double shiftSize = shift->getShiftSize();
     if (!Maths::isZero(shiftSize)){
         /* only bother if non zero */
         for (int j = 0; j < expiries->size(); j++) {
             /* Get trading time to the ith benchmark from cache*/
             double ttYearFrac = (*tradYears)[j]; 
             double rtShift = shift->rtVegaShift(ttYearFrac);
             for (int i = 0; i < strikes.size(); i++)
             {
                 (*vol)[i][j] += rtShift;
             }
         }

         if (shiftSize < 0){
             try{
                 vol->checkNonNegative();
             } catch (exception& e){
                 string m = "Vol shift ("+ Format::toString(shiftSize)+
                     ") results in negative vol for " + getName();
                 throw ModelException(&e, "VolSurface::sensShift", m);
             }
         }
     }
     return false; // none of our components has a vega type sensitivity
 }

 /* Apply a skew shift for a benchmark date. Vols for all strikes are
  * shifted by -shift * log(strike/spot)/log(1.1) */
 void VolSurface::skewShift(int expiryIndex, VegaSkewParallel* shift) {
     // curves with zero strike are currently allowed
     if (strikes.size() == 1 && Maths::isZero(strikes[0])){
         throw ModelException("VolSurface::skewShift", getName()+" - strike "
             "for vol curve is zero");
     }
     for (int i = 0; i < strikes.size(); i++){
         // strikes should have already been checked already to ensure > 0
         /* skew moved by shift amount at 110% strikes */
         double volShift = shift->skewShift(strikes[i]);
         (*vol)[i][expiryIndex] += volShift;
         if ((*vol)[i][expiryIndex] < -DBL_EPSILON) {
             // floor at zero rather than fail
             (*vol)[i][expiryIndex] = 0.0;
         }
     }
 }

 /* Apply a skew shift for a benchmark date. Vols for all strikes are
  * shifted by -shift * log(strike/spot)/log(1.1) */
 void VolSurface::skewShift(int expiryIndex, VegaSkewPointwise* shift) {
     // curves with zero strike are currently allowed
     if (strikes.size() == 1 && Maths::isZero(strikes[0])){
         throw ModelException("VolSurface::skewShift", getName()+" - strike "
             "for vol curve is zero");
     }
     for (int i = 0; i < strikes.size(); i++){
         // strikes should have already been checked already to ensure > 0
         /* skew moved by shift amount at 110% strikes */
         double volShift = shift->skewShift(strikes[i]);
         (*vol)[i][expiryIndex] += volShift;
         if ((*vol)[i][expiryIndex] < -DBL_EPSILON) {
             // floor at zero rather than fail
             (*vol)[i][expiryIndex] = 0.0;
         }
     }
 }


string VolSurface::sensName(VegaSkewParallel* shift) const{
    return getName();
}

bool VolSurface::sensShift(VegaSkewParallel* shift){
    double shiftSize = shift->getShiftSize();
    if (!Maths::isZero(shiftSize)){
        /* only bother if non zero */
        for (int j = 0; j < expiries->size(); j++){
            skewShift(j, shift);
        }
    }
     return false; // no more tweaking required here
}

string VolSurface::sensName(VegaSkewPointwise* shift) const{
    return getName();
}

ExpiryArrayConstSP VolSurface::sensExpiries(VegaSkewPointwise* shift) const{
    return expiries;
}


/* Apply a skew shift to a specific benchmarh. Vols for all strikes are
 * shifted by -shift * log(strike/spot)/log(1.1) */
bool VolSurface::sensShift(VegaSkewPointwise* shift){
    double shiftSize = shift->getShiftSize();
    if (!Maths::isZero(shiftSize)){
        int expiryIdx = shift->getExpiry()->search(expiries.get());
        skewShift(expiryIdx, shift);
    }
    return false; // no more tweaking required here
}

/** Returns name identifying this object for VolLevel */
string VolSurface::sensName(VolLevel* shift) const {
    return getName();
}

/** Shifts the object using given shift (see VolLevel::Shift)*/
bool VolSurface::sensShift(VolLevel* shift) {
    static const string method = "VolSurface::sensShift";
    try {
        double shiftSize = shift->getShiftSize();
        if (Maths::isNegative(shiftSize)) {
            throw ModelException(method,
                                 "vol level (" + 
                                 Format::toString(shiftSize) + 
                                 ") is -ve");
        }

        for (int j = 0; j < expiries->size(); j++) {
            for (int i = 0; i < strikes.size(); i++)
            {
                (*vol)[i][j] = shiftSize;
            }
        }
    }
    catch (exception&e ) {
        throw ModelException(e, method, 
                             "VolLevel scenario failed for " + getName());
    }
    return false;  // all done
}


/** Returns name identifying this object for VolParallelShift */
string VolSurface::sensName(VolParallelShift* shift) const {
    return getName();
}

/** Shifts the object using given shift (see VolParallelShift::Shift)*/
bool VolSurface::sensShift(VolParallelShift* shift)
{
    static const string method = "VolSurface::sensShift";
    try {
        double shiftSize = shift->getShiftSize();
        if (!Maths::isZero(shiftSize)) {

            for (int j = 0; j < expiries->size(); j++) {
                for (int i = 0; i < strikes.size(); i++)
                {
                    double shiftVol = (*vol)[i][j] + shiftSize;

                    if( Maths::isNegative(shiftSize) ) {
                        // Floor a downshift, but don't increase a low vol
                        if (shiftVol < VolParallelShift::MIN_SPOT_VOL) {
                            (*vol)[i][j] =
                                Maths::min((*vol)[i][j],
                                           VolParallelShift::MIN_SPOT_VOL);
                        }
                        else {
                            (*vol)[i][j] = shiftVol;
                        } 
                    } else {
                        // Don't floor an upshift
                        (*vol)[i][j] = shiftVol;
                    }
                }
            }

            // now remove -ve variance
            safeFwdVol();
        }
    }
    catch (exception&e ) {
        throw ModelException(e, method, 
                             "VolParallelShift scenario failed for "+getName());
    }
    return false; // all done
}

/** Returns name identifying this object for VolBenchmarkShift */
string VolSurface::sensName(VolBenchmarkShift* shift) const {
    return getName();
}
/** Shifts the object using given shift (see VolBenchmarkShift::Shift)*/
bool VolSurface::sensShift(VolBenchmarkShift* shift) {
    static const string method = "VolSurface::sensShift";
    try {
        double shiftSize = shift->getShiftSize();
        if (!Maths::isZero(shiftSize)) {
            int  i = 0;
            bool found = false;

            while (i < expiries->size() && !found) {
                found = shift->expiryEquals((*expiries)[i].get());
                i++;
            }

            if (!found) {
                throw ModelException(method, 
                                     "benchmark not found on vol surface " + 
                                     getName());
            }

            i--;  // we've stepped over it
                    
            // shift all strikes for this expiry
            for (int j = 0; j < strikes.size(); j++) {
                double shiftVol = (*vol)[j][i] + shiftSize;

                if( Maths::isNegative(shiftSize) ) {
                    // Floor a downshift, but don't increase a low vol
                    if (shiftVol < VolParallelShift::MIN_SPOT_VOL) {
                        (*vol)[j][i] = Maths::min((*vol)[j][i], 
                                                  VolParallelShift::MIN_SPOT_VOL);
                    }
                    else {
                        (*vol)[j][i] = shiftVol;
                    } 
                } else {
                    // Don't floor an upshift
                    (*vol)[j][i] = shiftVol;
                }
            }

        }

        // now remove -ve variance
        if (shift->lastShift()) {
            safeFwdVol();
        }
    }
    catch (exception&e ) {
        throw ModelException(e, method, "VolBenchmarkShift scenario "
                             "failed for " + getName());
    }
    return false; // all done
}

/** Returns name identifying this object for PowerVega */
string VolSurface::sensName(PowerVega* shift) const {
    return getName();
}
/** Shifts the object using given shift (see PowerVega::Shift)*/
bool VolSurface::sensShift(PowerVega* shift) {
    static const string method = "VolSurface::sensShift";
    try {
        double shiftSize = shift->getShiftSize();
        if (!Maths::isZero(shiftSize)){
            /* only bother if non zero */
            
            for (int j = 0; j < expiries->size(); j++) {
                /* Get trading time to the ith benchmark from cache*/
                double ttYearFrac = (*tradYears)[j]; 
                double rtShift = shift->powerShift(ttYearFrac);
                for (int i = 0; i < strikes.size(); i++)
                {
                    double shiftVol = (*vol)[i][j] + rtShift;

                    if( Maths::isNegative(rtShift) ) {
                        // Floor a downshift, but don't increase a low vol
                        if (shiftVol < VolParallelShift::MIN_SPOT_VOL) {
                            (*vol)[i][j] = 
                                Maths::min((*vol)[i][j], 
                                           VolParallelShift::MIN_SPOT_VOL);
                        }
                        else {
                            (*vol)[i][j] = shiftVol;
                        } 
                    } else {
                        // Don't floor an upshift
                        (*vol)[i][j] = shiftVol;
                    }
                }
            }

            // now remove -ve variance
            safeFwdVol();
        }
    }
    catch (exception&e ) {
        throw ModelException(e, method, 
                             "PowerVega scenario failed for " + getName());
    }
    return false; // all done
}

// remove fixed benchmarks in the past 
// useTradTime reflects that fact that we don't always have all the market
// data to hand (e.g. hols in time metric) - via Pyramid we can rely on
// getMarket being called, on s/sheet you only get calendar checking
// - compromise on real world vs making everything mutable
void VolSurface::zapPastBM(bool useTradTime) {
    static const string method("VolSurface::zapPastBM");
    try {
        // fixed always come before relative, so if first one is fixed ..
        if (BenchmarkDate::TYPE->isInstance((*expiries)[0].get())) {
            int i;
            int j;
            int skip = 0;
            // if a fixed benchmark is in the past, need to drop them out 
            // and resize the matrix
            for (i = 0; i < expiries->size(); i++) {
                Expiry*  expiry = (*expiries)[i].get();  // shorthand
                DateTime bmDate = expiry->toDate(baseDate);
                // if it's historic in calendar terms OR in trading time terms
                if (bmDate <= baseDate || 
                    (useTradTime && 
                     !Maths::isPositive(metric->yearFrac(baseDate, bmDate)))) {
                    if (!BenchmarkDate::TYPE->isInstance(expiry)) {
                        throw ModelException(method,
                                             "expiry " + 
                                             expiry->toString() + 
                                             " is historic (today is " + 
                                             baseDate.toString() + ")");
                    }
                    skip++;
                }
            }
            if (skip) {
                if (skip == expiries->size()) {
                    throw ModelException(method,
                                         "rolling to " + baseDate.toString() +
                                         " makes all vol benchmarks historic");
                }

                //uh-oh work to do
                CDoubleMatrixSP zapVol(new DoubleMatrix(vol->numCols(),
                                                        vol->numRows()-skip));

                // slice off skip rows of matrix
                for (i = 0; i < zapVol->numCols(); i++) {
                    for (j = 0; j < zapVol->numRows(); j++) {
                        (*zapVol)[i][j] = (*vol)[i][j+skip];
                    }
                }
                   
                // and ripple the expiries up
                for (i = 0; i < expiries->size() - skip; i++) {
                    (*expiries)[i] = (*expiries)[i+skip];
                }
                expiries->resize(expiries->size() - skip);

                vol = zapVol;
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method); 
    }
}

/** Shifts the object using given shift. */
bool VolSurface::sensShift(Theta* shift)
{
    static const string method("VolSurface::sensShift (theta)");
    try
    {
        baseDate = shift->rollDate(baseDate);       
        zapPastBM(true);  // remove historic fixed benchmarks
        buildYearFracCache();
    }
    catch (exception& e)
    {
        throw ModelException(e, method); 
    }
    return false;  // non of the components have Theta sensitivity
}

string VolSurface::sensName(DeltaSurface* shift) const{
    return getName();
}

bool VolSurface::sensShift(DeltaSurface* shift){
    double shiftSize = shift->getShiftSize();
    if (!Maths::isZero(shiftSize)){
        // only bother if non zero 

        checkCache();

        // what to do - spot moves to spot + shift
        // interpolate ATM curve using new spot
        // create new surface witk k-> k (S+dS)/S
        // preserve smile around new ATM vols
        double spot    = shift->getSpot();
        double newSpot = spot * (1.0 + shiftSize);

        CVolProcessedBSSP oldATM(getProcessedVol(spot));
        CVolProcessedBSSP newATM(getProcessedVol(newSpot));

        CDoubleMatrixSP newVols(new DoubleMatrix(strikes.size(),
                                                 expiries->size()));

        int i;

        for (int j = 0; j < expiries->size(); j++){
            double atmVol = oldATM->CalcVol(baseDate, (*dates)[j]);
            double twkVol = newATM->CalcVol(baseDate, (*dates)[j]);
            for (i = 0; i < strikes.size(); i++) {
                // measure smile at each strike
                double smile  = (*vol)[i][j] - atmVol;

                (*newVols)[i][j] = twkVol + smile;               
            }                        
        }

        // scale up all the strikes
        for (i = 0; i < strikes.size(); i++) {
            strikes[i] *= (newSpot/spot);
        }

        // and use our new vols
        vol = newVols;
    }
     return false; // no more tweaking required here
}

/** Returns name identifying this object for VolRelativeShift */
string VolSurface::sensName(VolRelativeShift* shift) const {
    return getName();
}   

/** Shifts the object using given shift (see VolRelativeShift::IShift)*/
bool VolSurface::sensShift(VolRelativeShift* shift) {
    static const string method("VolSurface::sensShift");
    try {
        for (int i = 0; i < expiries->size(); i++) {
            double shiftSize = shift->shiftSize(baseDate,(*expiries)[i]->toDate(baseDate));

            // shift all strikes for this expiry
            for (int j = 0; j < strikes.size(); j++) {
                double shiftVol = (*vol)[j][i] + shiftSize;

                if( Maths::isNegative(shiftSize) ) {
                    // Floor a downshift, but don't increase a low vol
                    if (shiftVol < VolParallelShift::MIN_SPOT_VOL) {
                        (*vol)[j][i] = Maths::min((*vol)[j][i], 
                                                  VolParallelShift::MIN_SPOT_VOL);
                    }
                    else {
                        (*vol)[j][i] = shiftVol;
                    } 
                } else {
                    // Don't floor an upshift
                    (*vol)[j][i] = shiftVol;
                }
            }
        }
         
        // now remove -ve variance
        safeFwdVol();
    }
    catch (exception& e) {
        throw ModelException(e, method, 
                             "VolRelativeShift scenario failed for "+getName());
    }
    return false; // all done
}

/** Returns name identifying this object for VolAbsoluteShift */
string VolSurface::sensName(VolAbsoluteShift* shift) const {
    return getName();
}   

/** Shifts the object using given shift (see VolAbsoluteShift::IShift)*/
bool VolSurface::sensShift(VolAbsoluteShift* shift) {
    static const string method("VolSurface::sensShift");
    try {
        for (int i = 0; i < expiries->size(); i++) {
            double shiftSize = shift->shiftSize(baseDate,(*expiries)[i]->toDate(baseDate));
            
            // shift all strikes for this expiry
            for (int j = 0; j < strikes.size(); j++) {
                double shiftVol = (*vol)[j][i] + shiftSize;

                if( Maths::isNegative(shiftSize) ) {
                    // Floor a downshift, but don't increase a low vol
                    if (shiftVol < VolParallelShift::MIN_SPOT_VOL) {
                        (*vol)[j][i] = Maths::min((*vol)[j][i], 
                                                  VolParallelShift::MIN_SPOT_VOL);
                    }
                    else {
                        (*vol)[j][i] = shiftVol;
                    } 
                } else {
                    // Don't floor an upshift
                    (*vol)[j][i] = shiftVol;
                }
            }
        }
         
        // now remove -ve variance
        safeFwdVol();
    }
    catch (exception& e) {
        throw ModelException(e, method, 
                             "VolAbsoluteShift scenario failed for "+getName());
    }
    return false; // all done
}

// adjust surface such that the fwd vol between benchmarks is at least 5%
void  VolSurface::safeFwdVol() {
    double v2 = VolParallelShift::MIN_FWD_VOL * VolParallelShift::MIN_FWD_VOL;
    for (int i = 0; i < strikes.size(); i++) {
        for (int j = 0; j < expiries->size()-1; j++) {

            double vol1 = (*vol)[i][j];
            double vol2 = (*vol)[i][j+1];

            double var1 = vol1*vol1*(*tradYears)[j];
            double var2 = vol2*vol2*(*tradYears)[j+1];

            double minVar = v2 * ((*tradYears)[j+1] - (*tradYears)[j]);

            if (var2 - var1 < minVar) {
                (*vol)[i][j+1] = sqrt((var1+minVar)/(*tradYears)[j+1]);
            }
        }
    }
}


/** Combines market and instrument data together to give a
    Processed Vol */
CVolProcessed* VolSurface::getProcessedVol(
    const CVolRequest* volRequest,
    const CAsset*      asset) const{
    checkCache();
    // we delegate the work of doing a 'double dispatch' to VolProcessedDispatch
    // Here 'double dispatch' means choosing the method based upon the type
    // of the vol and the request
    return VolProcessedDispatch::dispatch(this, volRequest, asset);
}

/** Same as getProcessedVol(CVolRequest, CAsset) but takes smart pointer
    to surface and solves memory ownership problems. Note static nethod.
    Method is now redundant because refCount is inside object. */
CVolProcessed* VolSurface::getProcessedVol(
    const VolSurfaceConstSP&  volSurface,
    const CVolRequest*        volRequest,
    const CAsset*             asset)
{
    return volSurface->getProcessedVol(volRequest, asset);
}    

/** Combines market and instrument data together to give a
    Processed Vol. Here the processed volatility is a processed
    struck volatility ie it reflects the combination of this
    CVolBase together with the supplied FX asset and the
    correlation between this CVolBase and the vol of the
    FX. */
CVolProcessed* VolSurface::getProcessedVol(
    const CVolRequest* volRequest,
    const CAsset*      eqAsset,
    const FXAsset*     fxAsset,
    const Correlation* eqFXCorr) const
{
    static const string routine("CVolProcessed::getProcessedVol");
    checkCache();
    if (!fxAsset){
        throw ModelException(routine, "NULL fx asset");
    }
    if (!eqFXCorr){
        throw ModelException(routine, "NULL correlation");
    }
    /* vol(struck)= sqrt(vol(equity) * vol(equity) + vol(FX) * vol(FX)
       + 2 * correlation * vol(equity) * vol(FX))
    */
    try{
        // create an atm vol request for the fx side of things
        CVolRequestSP atmInterp(new ATMVolRequest());
        CVolProcessedSP interpFXVol(fxAsset->getProcessedVol(atmInterp.get()));
        // cast to the type of vol we need
        CVolProcessedBSSP volFXBS = 
            CVolProcessedBSSP::dynamicCast(interpFXVol);

        // then calculate vol at each expiry - unfortunately vol method
        // can't handle everything in one go
        double firstFxVol = volFXBS->CalcVol(baseDate, (*dates)[0]);
        DoubleArray fxVols(dates->size()-1);
        volFXBS->CalcVol(*dates, CVolProcessedBS::fromFirst, fxVols);

        double fxSpot = fxAsset->getSpot();
        double corr = eqFXCorr->getCorrelation();
        // create new vol to store it in
        VolSurfaceSP struckVol(copy(this));
        // loop over strikes
        for (int i = 0; i < strikes.size(); i++){
            // scale strikes by spot fx
            struckVol->strikes[i] *= fxSpot;
            double eqVol = (*vol)[i][0];
            (*struckVol->vol)[i][0] = sqrt(eqVol * eqVol +
                                           firstFxVol * firstFxVol + 
                                           2.0 * corr * eqVol * firstFxVol);
            
            for (int j = 1; j < dates->size(); j++){
                double eqVol = (*vol)[i][j];
                double fxVol = fxVols[j-1];
                (*struckVol->vol)[i][j] = sqrt(eqVol * eqVol + fxVol * fxVol + 
                                               2.0 * corr * eqVol * fxVol);
            }
        }
        // this takes care of the struck vol - ensures that it is not freed
        // until the interpVol is freed
        return struckVol->getProcessedVol(volRequest, eqAsset);
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

/** Builds vol surface using existing values from surface but with
    new single strike, vols and bm dates. Note no clones of incoming
    data made */
VolSurface::VolSurface(const VolSurface*         surface,
                       double                    strike,
                       const DoubleArray&        vols,
                       const ExpiryArraySP&      expiries,
                       const DateTimeArraySP&    dates,  // must match expiries
                       const DoubleArraySP&      tradTime): // must match dates
    CVolBase(TYPE), name(surface->getName()), metric(surface->metric), 
    strikes(1, strike), vol(new DoubleMatrix(vols)), expiries(expiries), 
    baseDate(surface->baseDate), gotCache(true), tradYears(tradTime),
    dates(dates), strikeMode(VolSurface::ABSOLUTE) {
    validatePop2Object();
}

/** Creates a vol surface from an existing surface but uses the
    supplid strikes and vols - this is optimised for performance */
VolSurface::VolSurface(const VolSurface*    surface,
                       const DoubleArray&   strikes,
                       const CDoubleMatrix& matrix):
    CVolBase(TYPE), name(surface->getName()), metric(surface->metric), 
    strikes(strikes), vol(new DoubleMatrix(matrix)),
    expiries(surface->expiries), 
    baseDate(surface->baseDate), gotCache(true), tradYears(surface->tradYears),
    dates(surface->dates), strikeMode(VolSurface::ABSOLUTE) {
    validatePop2Object();
}


/** Constructor - in general though object is created via data dict.
    Note that the metric must already be populated with its market
    data (ie its holidays) */
VolSurface::VolSurface(const string&        volName,
                       const TimeMetric*    metric,
                       const DoubleArray&   strikes,
                       const DoubleMatrix&  vol,
                       const ExpiryArray*   expiries,
                       const DateTime&      baseDate):
    CVolBase(TYPE), name(volName), strikes(strikes), vol(new DoubleMatrix(vol)),
    baseDate(baseDate),
    gotCache(false), strikeMode(VolSurface::ABSOLUTE) {
    try{
        this->metric = TimeMetricSP(copy(metric));
        this->expiries = ExpiryArraySP(copy(expiries));
        validatePop2Object();
        buildYearFracCache();
    } catch (exception& e){
        throw ModelException(&e, "VolSurface::VolSurface");
    }
}

//// Like above constructor but takes variances rather than vols
VolSurface* VolSurface::createFromVars(const string&        volName,
                                       const TimeMetric*    metric,
                                       const DoubleArray&   strikes,
                                       const DoubleMatrix&  vars,
                                       const ExpiryArray*   expiries,
                                       const DateTime&      baseDate){
    // create surface
    VolSurfaceSP surf(new VolSurface(volName, metric, strikes,
                                     vars, expiries, baseDate));
    // then turn variances into vols
    for (int iStrike = 0; iStrike < surf->vol->numCols(); iStrike++){
        double* vols = (*surf->vol)[iStrike];
        for (int iDate = 0; iDate < surf->vol->numRows(); iDate++){
            vols[iDate] = sqrt(vols[iDate]/(*surf->tradYears)[iDate]);
        }
    }
    return surf.release();
}
            
/** Updates this vol surface with the supplied strikes and vols. Note
    that no copy is made of the vols. (Also note that the vol surface
    may at any point overwrite or alter them.) The return value is the
    new VolSurface to use - it may or may not be the same as 'this'.
    Use smart pointers to handle memory management.
    If validate is true then the strikes are validated to be increasing
    and distinct and the vols are validated to be positive. Only set
    to this false if you're confident that the data will always be
    good */
VolSurface* VolSurface::update(const CDoubleArray&    strikes,
                               const CDoubleMatrixSP& vols,
                               bool                   validate){
    if (getRefCount() == 1){
        // can reuse this vol
        vol = vols;
        this->strikes = strikes; // structure copy
        if (validate){
            validatePop2Object();
        }
        return this;
    }
    // must create a new copy
    return new VolSurface(this, strikes, *vols);
}

TimeMetricSP VolSurface::getTimeMetric()const{
    return metric;
}

void VolSurface::acceptValueDateCollector(const VolSurface*    volSurface,
                                          CValueDateCollector* collector)
{
    collector->valueDateValidate(volSurface->baseDate,
                                 volSurface->getName());
}

/** given a current spot level, get the next strike on the vol surface */
double VolSurface::getNextStrike(const double& strike,
                                 bool          isUp,
                                 bool&         offSurface) const
{
    int  idx          = 0;
    double nextStrike = 0.0;
	offSurface = false;

    // find first strike which is >= than strike 
    while ( idx < strikes.size()     && 
            strikes[idx]   < strike )
    {
        ++idx;
    }

    if ( isUp )
    {
        if ( idx == strikes.size() ) {
            offSurface = true;
        } else {
            nextStrike = strikes[idx];
        }
    }
    else
    {
        if ( idx == 0 ) {
            offSurface = true;
        } else {
            nextStrike = strikes[idx-1];
        }
    }
    return nextStrike;
}

/** returns the base date within the vol - needed as there is a tendency
    to use the VolSurface as a data holder into the library. Probably
    better to be able to get a new VolCurve class from a vol surface */
const DateTime& VolSurface::getBaseDate() const{
    return baseDate;
}


/** returns a double list of all strikes on the vol surface */
DoubleArraySP VolSurface::getAllStrikes() const
{
    DoubleArraySP doubleArray(new DoubleArray(strikes.size()));

    for(int i=0 ; i<strikes.size() ; ++i )
    {
        (*doubleArray)[i] = strikes[i];
    }
    return doubleArray;
}

/** handles the delta shift size adjustment */
void VolSurface::adjustDeltaShiftSize(ShiftSizeCollector* collector,
                                      const string assetName,
                                      double spot) const {
    collector->adjustShiftSize(getAllStrikes(), assetName, spot);
}

/* simple linear (non forward starting) interpolation at absolute strike */
CVolProcessedBS* VolSurface::getProcessedVol(double absoluteStrike) const{
    checkCache();
    return VolSurface::Interp::interpolateAbsLinearStrike(this, 
                                                          absoluteStrike);
}

/** Returns the array of expiries in the vol suurface. See comments
    under getBaseDate. To review in conjunction with parameterised vols */
ExpiryArrayConstSP VolSurface::getExpiries() const{
    return expiries;
}

/** Returns the array of expiries in the vol suurface converted to
    actual dates */
const DateTimeArray& VolSurface::getDates() const{
    return *dates;
}

/** Returns the array of strikes in the vol surface */
const DoubleArray& VolSurface::getStrikes() const{
    return strikes;
}

/** Returns the array of strikes in the vol surface */
const DoubleMatrix& VolSurface::getVolMatrix() const{
    return *vol;
}

/** This routine highlights the problem with using a market data cache - we
    have to check to see if we have our market data before we do any work */
void VolSurface::checkCache() const{
    if (!gotCache){
        buildYearFracCache();
    }
}

/** populate from market cache */
void VolSurface::getMarket(const IModel* model, const MarketData* market) {
    try{
        market->GetReferenceDate(baseDate);
        metric->getMarket(model, market);
        zapPastBM(true);  // remove historic fixed benchmarks
        buildYearFracCache();
    } catch (exception& e){
        throw ModelException(e, "VolSurface::getMarket" + name);
    }
}

PDFCalculator* VolSurface::getPDFCalculator(
    const PDFRequest* request,
    const CAsset*     asset) const {
    static const string method("VolSurface::getPDFCalculator");
    try {
        return new PDFDefaultLNStrike(baseDate,asset,request);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/* for reflection */
VolSurface::VolSurface(): CVolBase(TYPE), gotCache(false), strikeMode(VolSurface::ABSOLUTE) {}


class VolSurfaceHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VolSurface, clazz);
        SUPERCLASS(CVolBase);
        IMPLEMENTS(IVolatilityBS);
        IMPLEMENTS(IRestorableWithRespectTo<VolParallel>);
        IMPLEMENTS(VegaMatrix::IShift);
        IMPLEMENTS(ITweakableWithRespectTo<VolPointwise>);
        IMPLEMENTS(VegaSkewPointwise::IShift);
        IMPLEMENTS(VegaSkewParallel::IShift);
        IMPLEMENTS(RootTimeVega::IShift);
        IMPLEMENTS(VolLevel::Shift);
        IMPLEMENTS(VolParallelShift::Shift);
        IMPLEMENTS(VolBenchmarkShift::Shift);
        IMPLEMENTS(PowerVega::Shift);
        IMPLEMENTS(Theta::IShift);
        IMPLEMENTS(IPDFCalculator);
        IMPLEMENTS(DeltaSurface::IShift);
        IMPLEMENTS(VolRelativeShift::IShift);
        IMPLEMENTS(VolAbsoluteShift::IShift);
        EMPTY_SHELL_METHOD(defaultVolSurface);
        FIELD(name, "Vol identifier");
        FIELD(metric, "metric for trading time");
        FIELD(strikes, "strikes");
        FIELD(vol, "volatility matrix");
        FIELD(expiries, "benchmark dates");
        FIELD(baseDate, "base date");
        FIELD_MAKE_OPTIONAL(baseDate);
        FIELD(tradYears, "cached trading years");
        FIELD_MAKE_TRANSIENT(tradYears); // hide from dd interface
        FIELD(dates, "cached trading years");
        FIELD_MAKE_TRANSIENT(dates); // hide from dd interface
        FIELD(gotCache, "built our cache yet");
        FIELD_MAKE_TRANSIENT(gotCache); // hide from dd interface
        FIELD(strikeMode, "Strike Representation");
        FIELD_MAKE_OPTIONAL(strikeMode);
        // register the flavours of vol request that we support
        VolProcessedDispatch::registerVolMethod(
            VolSurface::Interp::interpolateWithTermStructure);
        VolProcessedDispatch::registerVolMethod(
            VolSurface::Interp::interpolateLinearStrike);
        VolProcessedDispatch::registerVolMethod(
            VolSurface::Interp::interpolateLinearStrikeSpread);
        VolProcessedDispatch::registerVolMethod(
            VolSurface::Interp::interpolateATM);
        VolProcessedDispatch::registerVolMethod(
            VolSurface::Interp::interpolateCliquet);
        VolProcessedDispatch::registerVolMethod(
            VolSurface::Interp::interpolateSpline);
        VolProcessedDispatch::registerVolMethod(
            VolSurface::Interp::interpolateTime);
        

        ClassSetAcceptMethod(VolSurface::acceptValueDateCollector);
        Addin::registerConstructor("VOL_SURFACE",
                                   Addin::MARKET,
                                   "Creates a handle to a vol surface",
                                   VolSurface::TYPE);

    }

    static IObject* defaultVolSurface(){
        return new VolSurface();
    }
};

VolSurface::~VolSurface(){}


CClassConstSP const VolSurface::TYPE = CClass::registerClassLoadMethod(
    "VolSurface", typeid(VolSurface), VolSurfaceHelper::load);


DRLIB_END_NAMESPACE
