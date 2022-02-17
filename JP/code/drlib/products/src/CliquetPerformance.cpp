//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : TemporalPerformance.cpp
//
//   Description : Measures performance of an asset through time
//
//   Date        : 7 Nov 2001
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/Performance.hpp"
#include "edginc/CliquetVolRequest.hpp"

DRLIB_BEGIN_NAMESPACE

/** Cliquet performance where perfType, strike, participation, and 
    sampleDates are the same across time */
class CliquetPerformance: public CObject,
                          virtual public IPerformance{
private:
    // all fields apply to all time periods
    string             perfType;         
    double             strike;           
    /* ideally we'd have a DateTimeCluster but IMS/EAS can't support this
       yet */
    DateTimeArray      sampleDates;    
    IntArray           dateTypes; // 1 = final cliquet date, 0 = ordinary date
    IntArray           numDatesInCliquet;  // transient field
public:
    static CClassConstSP const TYPE;

    static const int   CLIQUET_END_DATE; // = 1
    static const int   NORMAL_DATE;      // = 0

    class MCPerf: virtual public IPerformance::IMCPerf{
    private:
        /// fields ///
        const CliquetPerformance* performance;
        mutable DoubleArray       sumOut;
        mutable DoubleArray       cliqRefLevels;
        IntArray                  dateMap;
        bool                      isTrivialMap;
        mutable int               dateIndex; // where we're starting from
        mutable int               cliquetIndex; // where we're starting from
    public:
        /** Returns the number of performances that this object will be
            calculating ie the size of the perf parameter in calcPerf() */
        virtual int numPerformances() const{
            return performance->numDatesInCliquet.size();
        }

        /** Calculates the performance of each asset returning the results
            in perf (which must of the required size) */
        virtual void calcPerf(const IMCPathGenerator*  pathGen,
                              int                      iPath,
                              DoubleArray&             perf) const{
            perf = sumOut; // structure copy
            // record ref level for first cliquet
            cliqRefLevels[0] = pathGen->refLevel(0 /* iAsset */, iPath);
            const IntArray& dateTypes = performance->dateTypes;
            int iCliquet = cliquetIndex;
            int iDate = dateIndex;
            int endIdx = pathGen->end(0 /* iAsset */);
            int iStep = pathGen->begin(0 /* iAsset */);
            const double* path = pathGen->Path(0 /* iAsset */, iPath);
            if (isTrivialMap){ // for performance
                for (; iStep < endIdx; iStep++, iDate++){
                    perf[iCliquet] += path[iStep];
                    if (dateTypes[iStep] == CLIQUET_END_DATE){
                        iCliquet++; // move to next cliquet
                        cliqRefLevels[iCliquet] = path[iStep]; // save level
                    }
                }
            } else {
                // loop over our dates (a subset of path)
                for (iStep += dateMap[iStep]; iStep < endIdx; 
                     iStep++, iStep += dateMap[iStep], iDate++){
                    perf[iCliquet] += path[iStep];
                    if (dateTypes[iDate] == CLIQUET_END_DATE){
                        iCliquet++;
                        cliqRefLevels[iCliquet] = path[iStep];
                    }
                }
            }
            if (pathGen->doingPast()){ // preserve values
                sumOut = perf;
                dateIndex = iDate;
                cliquetIndex = iCliquet;
            }
            const IntArray& numDatesInCliquet = performance->numDatesInCliquet;
            // when doing past avoid cliquets in the future
            int   numCliquetsToDo = iCliquet;
            for (iCliquet = 0; iCliquet < numCliquetsToDo; iCliquet++){
                double cliqRefLevel = cliqRefLevels[iCliquet];
                double cliqPerf = (perf[iCliquet]/numDatesInCliquet[iCliquet])/
                    cliqRefLevel - performance->strike;
                perf[iCliquet] = 
                    IPerformance::Util::calcPerf(performance->perfType,
                                                 cliqPerf);
            }
        }

        // for the LogNormal path generator
        CVolRequestLN* getVolInterp(const IMCProduct*        mcProduct,
                                    const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
            int numCliquets = performance->numDatesInCliquet.size();
            const DateTime&  today = mcProduct->getToday();
            const IRefLevel* refLevel = mcProduct->getRefLevel();
            const DateTime&  startDate = refLevel->getAllDates().front();

            // get hold of the future strike dates
            int numFutureStrikeDates = numCliquets - cliquetIndex;
            DateTimeArray futStrikeDates(numFutureStrikeDates);
            int lastDatePos = 0;
            // loop over cliquets identifying start dates
            for (int iCliquet = 0; iCliquet < numCliquets; iCliquet++){
                if (iCliquet >= cliquetIndex){
                    // a live cliquet (ie not in the past)
                    futStrikeDates[iCliquet-cliquetIndex] = iCliquet == 0?
                        startDate: performance->sampleDates[lastDatePos-1];
                }
                lastDatePos += performance->numDatesInCliquet[iCliquet];
            }

            // same strike per cliquet (but may need to adjust first one)
            DoubleArray  strikes(futStrikeDates.size(), performance->strike);
            bool fwdStarting = startDate.isGreater(today);
            if (!fwdStarting && strikes.size() > 0){
                // need to set first level to 'current' strike - adjusted
                // additionally for any average out samples for this cliquet
                double cliqRefLevel = cliquetIndex == 0?
                    pathGen->refLevel(iAsset, 0): // in case there is no past
                    cliqRefLevels[cliquetIndex];
                int numInCliquet = 
                    performance->numDatesInCliquet[cliquetIndex];
                int numRemaining = 1; // for when dateTypes[iDate] == 1
                for (int iDate = dateIndex; 
                     performance->dateTypes[iDate] == NORMAL_DATE;
                     iDate++){
                    numRemaining++;
                }
                strikes[0] = (numInCliquet * cliqRefLevel * performance->strike
                              - sumOut[cliquetIndex])/ numRemaining;
            }
            return new CliquetVolRequest(fwdStarting, 
                                         futStrikeDates, 
                                         performance->sampleDates.back(),
                                         strikes);
        }

        /** Returns the maximum factor by which this Performance
            scales the delta for asset i wrt all paths. This is
            typically participation/ref level. Note that the path
            generator is the past one. To do: review this (fnc needed
            for QuickGreeks) */
        virtual double maxDeltaScalingFactor(
            const IMCPathGenerator* pastPathGen, 
            int                     iAsset) const{
            double refLevel = pastPathGen->refLevel(iAsset,
                                                    0 /* path irrelevant*/ );
            return (1.0/refLevel); // participation is 1
        }

        MCPerf(const CliquetPerformance* performance, 
               const DateTime&           today,
               const SimSeries*          allDates): 
            performance(performance), 
            sumOut(performance->numDatesInCliquet.size()),
            cliqRefLevels(performance->numDatesInCliquet.size()+1),
            dateIndex(0), cliquetIndex(0) {
            static const string routine("CliquetPerformance::MCPerf");
            try{
                /* get mapping for dates (tells us our how dates map onto all
                   simulation dates) */
                dateMap = allDates->createMap(0,// iAsset(same dates per asset)
                                              performance->sampleDates, 
                                              isTrivialMap);
            } catch (exception& e){
                throw ModelException(e, routine);
            }
        }
    };

    virtual PerformanceSVSP getStateVarGen(
        IRefLevel::IStateVarGenSP refLevelGen, int nbAssets) {
        static const string routine = "CliquetPerformance::getStateVarGen";

        throw ModelException(routine,"State variable representation not supported");
    }

    
    friend class MCPerf;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(CliquetPerformance, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IPerformance);
        EMPTY_SHELL_METHOD(defaultCliquetPerformance);
        FIELD(perfType, "Performance Type (for all cliquets)");
        FIELD(strike,   "Strike (for all cliquets)");
        FIELD(sampleDates,  "Cliquet averaging out dates");
        FIELD(dateTypes, "Indicates final averaging date in "
                     "each cliquet");
        FIELD(numDatesInCliquet, "");
        FIELD_MAKE_TRANSIENT(numDatesInCliquet);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
     
    // for reflection
    CliquetPerformance(): CObject(TYPE), strike(0.0) {}

    static IObject* defaultCliquetPerformance(){
        return new CliquetPerformance();
    }

public:
    IMCPerf* getMCPerf(const DateTime&    today,
                       const SimSeries*   allDates) const{
        return new MCPerf(this, today, allDates);
    }

    /** Records the dates for which this object needs points for 
        within the MC simulation */
    virtual void recordDates(SimSeries* simSeries) const{
        simSeries->addDates(sampleDates);
    }
        
    // validation
    void validatePop2Object(){
        static const string routine("CliquetPerformance::validatePop2Object");
        IPerformance::Util::validatePerfFlags(perfType, strike, sampleDates);
        if (sampleDates.size() != dateTypes.size()){
            throw ModelException(routine, Format::toString(sampleDates.size())+
                                 " sample dates but "+
                                 Format::toString(dateTypes.size())+
                                 " date types. They must be equal in length");
        }
        // populate numDatesInCliquet array
        numDatesInCliquet = IntArray(1);
        for (int i = 0; i < sampleDates.size(); i++){
            numDatesInCliquet.back()++;
            if (dateTypes[i] == CLIQUET_END_DATE){
                numDatesInCliquet.push_back(0);
            } else if (dateTypes[i] != NORMAL_DATE){
                throw ModelException(routine, "Cliquet Flags must be 0 or 1");
            }
        }
        if (numDatesInCliquet.back() != 0){
            throw ModelException(routine,
                                 "Last average cliquet date must be of"
                                 " type 1");
        }
        numDatesInCliquet.pop_back(); // remove bogus final element
    }
};
CClassConstSP const CliquetPerformance::TYPE =
CClass::registerClassLoadMethod("CliquetPerformance", 
                                typeid(CliquetPerformance), load);

const int CliquetPerformance::CLIQUET_END_DATE = 1;
const int CliquetPerformance::NORMAL_DATE      = 0;

// force linker to include this file (avoid having header file) */
bool CliquetPerformanceLoad() {
    return (CliquetPerformance::TYPE != 0);
}

DRLIB_END_NAMESPACE

    
