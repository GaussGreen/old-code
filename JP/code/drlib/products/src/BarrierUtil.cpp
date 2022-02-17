//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BarrierUtil.cpp
//
//   Description : 
//
//   Date        : Q Hou. Feb 2005
//
//
//   $Log: BarrierUtil.cpp,v $
//   Revision 1.10  2005/04/27 08:55:49  qhou
//   BarrierMgr knownCF bug fix w.r.t. sticky barrier. change TrivialPay to link obsDates with last prev mon date
//
//   Revision 1.9  2005/04/25 04:01:10  qhou
//   relevantMonDt fnction for barUtil/Pay
//
//   Revision 1.8  2005/04/18 09:52:11  qhou
//   redefine class names IBarrierPay, IBarrierMgr and BarrierUtil's preprocess()
//
//   Revision 1.7  2005/04/15 04:22:14  qhou
//   add interface to collect knownCf/phyDelivery. add initialHitVal function to sticky barrier
//
//   Revision 1.6  2005/04/01 08:25:21  qhou
//   missed 1 spot in prev bug fix
//
//   Revision 1.5  2005/03/31 10:34:30  qhou
//   fix Unix compile error need IBarrierUtil::ICanHaveSticky instead of ICanHaveSticky
//
//   Revision 1.4  2005/03/30 06:43:14  qhou
//   replace IBarrier by IBarrierUtil interface and move sticky/coupler/reverser etc into cpp
//
//   Revision 1.3  2005/03/02 16:50:00  qhou
//   fix too strict validation in BarrierCoupler
//
//   Revision 1.2  2005/03/02 12:02:28  qhou
//   fix UMR msg from purify
//
//   Revision 1.1  2005/02/25 09:39:56  qhou
//   IBarrier and IBarrierPayObject classes
//
//
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/Algorithm.hpp"
#include "edginc/Barrier.hpp"
#include "edginc/PhysicalDelivery.hpp"
#include "edginc/Maths.hpp"


DRLIB_BEGIN_NAMESPACE


//////////////////////////////////////////////////////////////////////////////
//
// Several concrete instances of IBarrierUtil used as wrapper on basic barriers
//
//////////////////////////////////////////////////////////////////////////////

       
class BarrierReverser : public IBarrierUtil, virtual public IBarrierUtil::ICanHaveSticky {
public:
    virtual bool getHitDateAndStep(int& step, DateTime& hitDate) const{
        return barrier->getHitDateAndStep(step, hitDate);
    }
    virtual bool getEvents(EventResults* events, bool useIsHitFlags, const string& barrierName, StringArray* assNames, const DateTime& breachDate, int maxHits) const{
        return barrier->getEvents(events, useIsHitFlags, barrierName, assNames, breachDate, maxHits);
    }
    virtual bool getHitDataForBreachEvent(double& refLevel, double& barrierLevel, double& assetlevel,string& monType) const{
        return barrier->getHitDataForBreachEvent(refLevel,barrierLevel,assetlevel, monType);
    }

    friend class BarrierFactory; // so can access barrier component

    BarrierReverser(IBarrierUtilSP barrier) : barrier(barrier)
    {
        if( !barrier )
            throw ModelException("BarrierReverser", "Underlying barrier is empty");
    }
    ~BarrierReverser(){};

    // relevant dates from this barrier itself
    virtual void relevantDates(DateTimeArray &relevantDts, bool &hasContinuous) const
    { barrier->relevantDates(relevantDts, hasContinuous); }

    // linkage, memory management w.r.t. simulation timeline
    virtual void preprocess(DateTime valueDate, DateTimeArrayConstSP simDates)
    { barrier->preprocess(valueDate, simDates); }

    virtual const DateTimeArray& relevantMonDates() const
    { return barrier->relevantMonDates(); }

    virtual void pathUpdated(const IMCPathGenerator* pathGenIn)
    { barrier->pathUpdated(pathGenIn); }

    virtual double hitValue(int obsStep) const
    { return 1.0 - barrier->hitValue(obsStep); }

    virtual bool isSticky() const
    { 
        IBarrierUtil::ICanHaveSticky *s = dynamic_cast<IBarrierUtil::ICanHaveSticky*>(barrier.get());
        return s&&s->isSticky();
    }

    virtual double initialHitVal() const
    {
        IBarrierUtil::ICanHaveSticky *s = dynamic_cast<IBarrierUtil::ICanHaveSticky*>(barrier.get());
        if( !s || !s->isSticky() )
            throw ModelException("BarrierReverser::initialState", "Internal error. Barrier not sticky");
        
        return 1.0 - s->initialHitVal();
    }

    virtual int periodEndStep(int obsStep) const
    { 
        if( !isSticky() )
            throw ModelException("BarrierReverser::periodEndStep", "Internal error. Not sticky");
        return dynamic_cast<IBarrierUtil::ICanHaveSticky*>(barrier.get())->periodEndStep(obsStep);
    }

protected:
    IBarrierUtilSP  barrier;
};
   
class BarrierReverserBend : public BarrierReverser, public IBarrierUtil::ICanHaveBend {
public:
    BarrierReverserBend(IBarrierUtilSP barrier) : BarrierReverser(barrier){};
    ~BarrierReverserBend(){};
    
    virtual void bendDates(DateTimeArray &startDates, DateTimeArray &endDates) const
    { dynamic_cast<IBarrierUtil::ICanHaveBend*>(barrier.get())->bendDates(startDates, endDates); }

    virtual double hitValueBend(int obsStep, int payObsStep) const
    { return 1.0 - dynamic_cast<IBarrierUtil::ICanHaveBend*>(barrier.get())->hitValueBend(obsStep, payObsStep); }
};

 
// hit value is sticky. once hit, the hit value carries to end of sim dates
// for this type, no need to enforce that hitValue() is called only for relevant obs step
class BarrierSticky :   public IBarrierUtil, 
                        virtual public IBarrierUtil::ICanHaveSticky 
{
public:

    virtual bool getHitDateAndStep(int& step, DateTime& hitDate) const{
        return barrier->getHitDateAndStep(step, hitDate);
    }
    virtual bool getEvents(EventResults* events, bool useIsHitFlags, const string& barrierName, StringArray* assNames, const DateTime& breachDate, int maxHits) const{
        return barrier->getEvents(events, useIsHitFlags, barrierName, assNames, breachDate, maxHits);
    }

    virtual bool getHitDataForBreachEvent(double& refLevel, double& barrierLevel, double& assetLevel,string& monType) const{
        return barrier->getHitDataForBreachEvent(refLevel, barrierLevel, assetLevel, monType);
    }

    BarrierSticky(IBarrierUtilSP barrier) 
        : barrier(barrier), lastCalcStep(-1)
    {
        if( !barrier )
            throw ModelException("BarrierSticky", "Underlying barrier is empty");
        if( dynamic_cast<BarrierReverser *>(barrier.get()) )
            throw ModelException("BarrierSticky", "Underlying barrier can not be BarrierReverser");
    }
    ~BarrierSticky(){};

    virtual void relevantDates(DateTimeArray &relevantDts, bool &hasContinuous) const
    { 
        barrier->relevantDates(relevantDts, hasContinuous); 
    }

    virtual void preprocess(DateTime valueDate, DateTimeArrayConstSP simDates)
    {
        static string method = "BarrierSticky::preprocess";

        barrier->preprocess(valueDate, simDates);
        
        // create mapping to relevant dates
        // make sure it's part of sim date and create mapping
        const DateTimeArray& relevantMonDts = barrier->relevantMonDates();
        if( !DateTime::isSubset(*simDates, relevantMonDts) )
            throw ModelException(method, ". relevant monitor dates must be subset of sim dates");

        bool dummy;
        relevantMap = DateTime::createMapping(*simDates, relevantMonDts, dummy);

        hitValues.resize(simDates->size());
    }

    // return all sim dates after the first mon date of the underlyer
    virtual const DateTimeArray& relevantMonDates() const
    { return barrier->relevantMonDates(); }

    // reset 'lastCalcStep' member
    virtual void pathUpdated(const IMCPathGenerator* pathGenIn)
    {
        barrier->pathUpdated(pathGenIn);
        lastCalcStep = Maths::max(-1, Maths::min(lastCalcStep, pathGenIn->begin(0)-1));
    }
    
    virtual double hitValue(int obsStep) const
    {
        // if sticky, but already calculated, return cached value
        if( obsStep <= lastCalcStep )
            return hitValues[obsStep];

        // if not already calculated, calculate from end of last calc step
        int i = lastCalcStep;
        double val = (i<0?0.0:hitValues[i]);
        while(i<obsStep)
        {
            i++;
            if ( relevantMap[i] == 0 )
                val += (1.0-val) * barrier->hitValue(i);

            hitValues[i] = val;
            if( Maths::isZero(val - 1.0) )
            {
                // propagate once hit
                while((++i)<hitValues.size()) 
                    hitValues[i] = val;
                i--;
            }
        }
        lastCalcStep = i;

        return val;
    }

    virtual bool isSticky() const
    { return true; }

    virtual double initialHitVal() const
    { return 0.0; }

    virtual int periodEndStep(int obsStep) const
    { return hitValues.size()-1; }

protected:
    IBarrierUtilSP              barrier;

    // *** transient fields below *** 
    // array with same size as sim dates to map sim step to relevant dates
    IntArray                    relevantMap;

    // array of cached hitValues
    mutable DoubleArray         hitValues;

    // need to record end of calculation already finished
    mutable int                 lastCalcStep;
};


class BarrierStickyBend : public BarrierSticky, virtual public IBarrierUtil::ICanHaveBend  {
public:
    BarrierStickyBend(IBarrierUtilSP barrier) : BarrierSticky(barrier)
    {
        if( !dynamic_cast<IBarrierUtil::ICanHaveBend *>(barrier.get()) )
            throw ModelException("BarrierStickyBend", "barrier must be of type IBarrierUtil::ICanHaveBend");
    }
    ~BarrierStickyBend(){};

    virtual void preprocess(DateTime valueDate, DateTimeArrayConstSP simDates)
    {
        static string method = "BarrierStickyBend::preprocess";
        BarrierSticky::preprocess(valueDate, simDates);

        // create mapping to bend period start/end dates
        DateTimeArray bendStartDts, bendEndDts;
        bendDates(bendStartDts, bendEndDts);
        if( bendEndDts.size() ==0 )
        {
            bendEndMap.resize(simDates->size());
            for(int i=0; i<simDates->size(); i++) bendEndMap[i] = -1;
            return;
        }

            DateTime::ensureIncreasing(bendEndDts, method + ". barrrier bend end dates must in strict ascending order", true);
            if( !DateTime::isSubset(*simDates, bendEndDts) )
                throw ModelException(method, ". barraier bend period end dates must be subset of sim dates");

            bool dummy;
            bendEndMap = DateTime::createMapping(*simDates, bendEndDts, dummy);

            bendEndStep.resize(bendEndDts.size());
            bendStartStep.resize(bendEndDts.size());
            int i, j;
            for(i=0, j=0; i<simDates->size(); i++)
            {
                if( bendEndMap[i] == 0 )
                {
                    // look for start of bend period
                    bendStartStep[j] = i;
                    while( bendStartStep[j] >= 0 && 
                        (*simDates)[bendStartStep[j]] >= bendStartDts[j] )
                        bendStartStep[j]--;
                    bendStartStep[j]++;

                    bendEndStep[j] = i;
                    bendEndMap[i] = j;
                    j++;
                
                    if( j== bendEndDts.size() )
                    {
                        i++;
                        while(i<simDates->size())
                        {
                            bendEndMap[i++] = -1;
                        }
                    }
                }
                else
                    bendEndMap[i] = j;
            }

            lastCalcStepBend.resize(bendEndStep.size());
            hitValuesBend.resize(bendEndStep.size());
            for(i=0; i<bendEndStep.size(); i++)
                hitValuesBend[i].resize(bendEndStep[i]+1);


    }
    virtual void pathUpdated(const IMCPathGenerator* pathGenIn)
    {
        BarrierSticky::pathUpdated(pathGenIn);
        for(int i=0; i<lastCalcStepBend.size(); i++)
            lastCalcStepBend[i] = Maths::max(-1, Maths::min(lastCalcStepBend[i], pathGenIn->begin(0)-1));
    }

    virtual void bendDates(DateTimeArray &startDates, DateTimeArray &endDates) const
    { dynamic_cast<IBarrierUtil::ICanHaveBend*>(barrier.get())->bendDates(startDates, endDates); }

    virtual double hitValueBend(int obsStep, int payObsStep) const
    {
        static string method = "BarrierStickyBend::hitValueBend";

        if( obsStep > payObsStep )
            throw ModelException(method, "obsStep > payObsStep");

        // find bend period, if is before start of bend period, call base schedule hitValue
        int bendPeriodIdx = bendEndMap[payObsStep];
        if( bendPeriodIdx < 0 || obsStep < bendStartStep[bendPeriodIdx] )
            return hitValue(obsStep);

        // if sticky, but already calculated, return cached value
        if( obsStep <= lastCalcStepBend[bendPeriodIdx] )
            return hitValuesBend[bendPeriodIdx][obsStep];

        // if not already calculated, start from base schedule hit value at start of step (which may invoke calculation)
        // or prev finished bend step, whichever is later
        int i = lastCalcStepBend[bendPeriodIdx];
        double val;
        if(i<bendStartStep[bendPeriodIdx])
        {
            i = bendStartStep[bendPeriodIdx]-1;
            val = i<0?0.0:hitValue(i);
        }
        else
            val = hitValuesBend[bendPeriodIdx][i];

        while(i<obsStep)
        {
            i++;
            if ( relevantMap[i] == 0 )
                val += (1.0-val) *  dynamic_cast<IBarrierUtil::ICanHaveBend*>(barrier.get())->hitValueBend(i, payObsStep);
            
            hitValuesBend[bendPeriodIdx][i] = val;
            if( Maths::isZero(val - 1.0) )
            {
                // propagate once hit
                while((++i)<=bendEndStep[bendPeriodIdx])
                    hitValuesBend[bendPeriodIdx][i] = val;
                i--;
            }
        }            
        lastCalcStepBend[bendPeriodIdx] = i;

        return val;
    }

    virtual int periodEndStep(int obsStep) const
    {
        int bendPeriodIdx = bendEndMap[obsStep];
        return (bendPeriodIdx>=0)?bendEndStep[bendPeriodIdx]:BarrierSticky::periodEndStep(obsStep);
    }

protected:
    // *** transient fields below ***
    // mapping between bend period end date's step and cardinal of period
    IntArray            bendEndMap;

    // array of bend period start/end index, both inclusive.
    IntArray            bendStartStep;
    IntArray            bendEndStep;

    // array of hitValue arrays, one for each bend period
    mutable DoubleArrayArray    hitValuesBend;

    // keep track of steps already calculated. 
    mutable IntArray            lastCalcStepBend; 
};

class BarrierCoupler : virtual public IBarrierUtil {
public:
    virtual bool getHitDateAndStep(int& step, DateTime& hitDate) const{
        //we are cheating here, we know that when a KI and KO barriers are coupled together,
        //the KI barrier is the last one and we can use that for the event
         return barriers[barriers.size()-1]->getHitDateAndStep(step, hitDate);
    }
    virtual bool getEvents(EventResults* events, bool useIsHitFlags, const string& barrierName, StringArray* assNames, const DateTime& breachDate, int hits) const{
        static const string method = "BarrierCoupler::getEvents";
        try {
            if (useIsHitFlags) {

                StringArraySP assetNames(new StringArray(0));
                DoubleArraySP assetLevels(new DoubleArray(0));
                DoubleArraySP barrLevels(new DoubleArray(0));
                
                assetNames->push_back("from breach flag");
                assetLevels->push_back(0.0);
                barrLevels->push_back(0.0);
                events->addEvent(new BarrierBreach(breachDate,
                                barrierName,
                                "from breach flag",
                                isOut?BarrierBreach::KNOCK_OUT : BarrierBreach::KNOCK_IN,
                                isUp,
                                0,
                                assetNames,
                                assetLevels,
                                barrLevels));
                return false;
            }
            bool processEvent = false;
            for (unsigned int iAsset=0; iAsset<barriers.size();++iAsset) {
                processEvent = processEvent || barriers[iAsset]->getEvents(events, useIsHitFlags, barrierName, assNames, breachDate, hitsRecorded);
            }

            if (processEvent) {
                StringArraySP assetNames(new StringArray(0));
                DoubleArraySP assetLevels(new DoubleArray(0));
                DoubleArraySP barrLevels(new DoubleArray(0));
                vector<int> AssetIndices(0);

                string monType;
                for(unsigned int iAsset=0; iAsset<barriers.size(); iAsset++)
                {
                    double refLevel = 0;
                    double assetLevel = 0.0;
                    double barrierLevel = 0.0;
                    if (barriers[iAsset]->getHitDataForBreachEvent(refLevel, barrierLevel,assetLevel,monType)) {
                        AssetIndices.push_back(iAsset);
                        assetNames->push_back((*assNames)[iAsset]);
                        barrLevels->push_back(barrierLevel);
                        assetLevels->push_back(assetLevel);
                    }
                }
                if (AssetIndices.size()>0) {
                    events->addEvent(new BarrierBreach(breachDate,
                                                       barrierName,
                                                       monType,
                                                       isOut?BarrierBreach::KNOCK_OUT : BarrierBreach::KNOCK_IN,
                                                       isUp,
                                                       hits--,
                                                       assetNames,
                                                       assetLevels,
                                                       barrLevels));
                }
            }
            return false;
        }catch (exception& e) {
            throw ModelException(e, method);
        }
    } 
    virtual bool getHitDataForBreachEvent(double& refLevel, double& barrierLevel, double& assetlevel,string& monType) const{
        return false; 
    }

    BarrierCoupler(vector<IBarrierUtilSP> barriers, int numHit, bool isOut, bool isUp)
        : barriers(barriers), numHit(numHit), isOut(isOut),hitsRecorded(0), isUp(isUp)
    {
        if( barriers.size() < 1 )
            throw ModelException("BarrierCoupler", "Must contain more than 1 barriers as components");
        if( (int)barriers.size() < numHit || numHit < 1)
            throw ModelException("BarrierCoupler", "numHit must between 1 and " + Format::toString((int)barriers.size()));
        for(unsigned int i=0; i<barriers.size(); i++)
        {
            if( !barriers[i] )
                throw ModelException("BarrierCoupler", "barrier " + Format::toString((int)i) + " is empty");
        }

        hitValTemp.resize(barriers.size());
    }
    ~BarrierCoupler(){};

    // relevant dates from this barrier itself
    virtual void relevantDates(DateTimeArray &relevantDts, bool &hasContinuous) const
    {
        barriers[0]->relevantDates(relevantDts, hasContinuous);
        for(unsigned int i=1; i<barriers.size(); i++)
        {
            bool elemIsCont;
            DateTimeArray elemDates;
            barriers[i]->relevantDates(elemDates, elemIsCont);
            relevantDts = DateTime::merge(relevantDts, elemDates);
            hasContinuous = (hasContinuous || elemIsCont);
        }
    }

    // linkage, memory management w.r.t. simulation timeline
    virtual void preprocess(DateTime valueDt, DateTimeArrayConstSP simDts)
    {
        relevantMonDts.clear();
        for(unsigned int i=0; i<barriers.size(); i++)
        {
            barriers[i]->preprocess(valueDt, simDts);
            relevantMonDts = DateTime::merge(relevantMonDts, barriers[i]->relevantMonDates());
        }
    }

    virtual const DateTimeArray& relevantMonDates() const
    {
        return relevantMonDts;
    }

    virtual double hitValue(int obsStep) const
    {
        for(unsigned int i=0; i<barriers.size(); i++)
        {
            hitValTemp[i] = barriers[i]->hitValue(obsStep);
        }
    
        // sorted in ascending order (in the sense of rainbow, ie.
        // index 0 is highest value!) and use the numHits'th highest value;
        Algorithm::shellSort(hitValTemp);
        return hitValTemp[numHit-1];    
    }
    
    virtual void pathUpdated(const IMCPathGenerator* pathGenIn)
    {
        for(unsigned int i=0; i<barriers.size(); i++)
            barriers[i]->pathUpdated(pathGenIn);
    }
        
protected:
   // component barriers
    vector<IBarrierUtilSP>  barriers;
    
    // threshold used to combine barrier's hitValues into 1 hitValue
    int                     numHit;

    // transient field for temp memory
    DateTimeArray           relevantMonDts;
    mutable DoubleArray     hitValTemp;
    bool                    isOut;
    bool                    isUp;
    mutable int             hitsRecorded;
};

class BarrierCouplerBend : public BarrierCoupler, virtual public IBarrierUtil::ICanHaveBend   {
public:
    BarrierCouplerBend(vector<IBarrierUtilSP> barriers, int numHit, bool isOut, bool isUp) 
        : BarrierCoupler(barriers, numHit, isOut, isUp){};
    ~BarrierCouplerBend(){};

    virtual void bendDates(DateTimeArray &startDates, DateTimeArray &endDates) const
    {
        DateTimeArray sDts, eDts;
        startDates.clear();
        endDates.clear();
        for(unsigned int i=0; i<barriers.size(); i++)
        {
            IBarrierUtil::ICanHaveBend* bendView = dynamic_cast<IBarrierUtil::ICanHaveBend*>(barriers[i].get());
            if( bendView )
            {
                bendView->bendDates(sDts, eDts);
                if( sDts.size() != 0 && eDts.size() != 0 )
                {
                    if( startDates.size() == 0 )
                    {
                        startDates = sDts;
                        endDates = eDts;
                    }
                    if( !DateTime::equals(sDts, startDates) || !DateTime::equals(eDts, endDates) )
                    throw ModelException("BarrierCouplerBend::bendDates", 
                            "Currently does not support different bend period among barrier components");
                }
            }
        }
    }
    
    virtual double hitValueBend(int obsStep, int payObsStep) const 
    {
        for(unsigned int i=0; i<barriers.size(); i++)
        {
            IBarrierUtil::ICanHaveBend* bendView = dynamic_cast<IBarrierUtil::ICanHaveBend*>(barriers[i].get());
            if( bendView )
                hitValTemp[i] = bendView->hitValueBend(obsStep, payObsStep);
            else
                hitValTemp[i] = barriers[i]->hitValue(obsStep);
        }
    
        // sorted in ascending order & use the numHits'th highest value;
        Algorithm::shellSort(hitValTemp);
        return hitValTemp[numHit-1];
    }
};


/****************************************************************************/

// provide access to barrier hit value
class TrivialBarrierPay : public IBarrierPay
{
public:
    TrivialBarrierPay(const DateTimeArray& obsDts) : obsDates(obsDts){};
    ~TrivialBarrierPay(){};

    void preprocess(DateTime valueDt, DateTimeArrayConstSP simDts, const DateTimeArray& monitorDates, const YieldCurve* discount)
    {
        if( monitorDates.front() > obsDates.front() )
            throw ModelException("TrivialBarrierPay:preprocess", "Obs date " + obsDates.front().toString() + 
                                 " has no monitor date associated with it");

        // match obs dates with last prev (inclusive) mon dates
        relevantMonDts.clear();
        IntArray counts;
        int i, j;
        for(i=0, j=0; i<obsDates.size(); i++)
        {
            while( j<monitorDates.size() && monitorDates[j] <= obsDates[i] )
                j++;
            j--;
            
            if( relevantMonDts.size()==0 || monitorDates[j] != relevantMonDts.back() )
            {
                relevantMonDts.push_back(monitorDates[j]);
                counts.push_back(0);
            }
            counts.back()++;
        }

        payPerStep.resize(simDts->size());
        for(i=0, j=0; i<simDts->size(); i++)
        {
            payPerStep[i] = (j<relevantMonDts.size() && (*simDts)[i]==relevantMonDts[j])?counts[j++]:0.0;
        }
    };

    bool payOnce() const 
    { return false; }

    bool payIfHitOne() const 
    { return true; }

    void payObsDates(DateTimeArray &obsDts, bool &hasContinuous) const
    { 
        obsDts = obsDates; 
        hasContinuous = false;
    }

    virtual const DateTimeArray& relevantMonDates() const
    { return relevantMonDts; }

    double value(int step, DateTime valueDate, const IMCPathGenerator* pathGen)
    { return payPerStep[step]; }
    
    void collect(int obsStep, const IMCPathGenerator* pathGenIn, double scalingFactor,
                 CashFlowArray& knownCFs, PhysicalDeliveryArray& phyDs)
    {};
 
private:
    DateTimeArray     obsDates;
    DateTimeArray     relevantMonDts;
    DoubleArray       payPerStep;
};

              
/****************************************************************************/
// base class that takes care of generic payoff calculation over a path
// taking into account of barrier bending periods
// allow customization w.r.t. payoff calculation at a local step
// to be derived by instances under StateVar context and pre-StateVar implementation
/****************************************************************************/

class BarrierMgrBase : virtual public IBarrierMgr {
public:
    // checks if there have been barrier breaches
    // and generates an event for each one
    // used by RDK
    void getEvents(const  IMCPathGenerator*  pathGenIn,
        DateTime valueDate,
        EventResults* events,
        const string& barrierName,
        StringArray* assNames) {
        double hitVal;

        static const string method = "BarrierMgrBase::getEvents";
        try
        {
            int beginStep = pathGenIn->begin(0), endStep = pathGenIn->end(0) - 1; 
            IBarrierUtil::ICanHaveBend* bendView = dynamic_cast<IBarrierUtil::ICanHaveBend*>(barrier.get());
            hitVal = bendView?bendView->hitValueBend(endStep, endStep):barrier->hitValue(endStep);
            hitVal = pay->payIfHitOne()?hitVal:(1.0-hitVal);
            if (Maths::isPositive(hitVal)) {
                DateTime hitDate;
                int step;
                bool usingHitFlag = barrier->getHitDateAndStep(step, hitDate);
                
                barrier->getEvents(events,usingHitFlag,barrierName,assNames,hitDate, 0);
            }
        }
        catch (exception& e)
        {
            throw ModelException(e,method);        	
        }

    }

    BarrierMgrBase(IBarrierUtilSP barrier, IBarrierPaySP pay) 
        : barrier(barrier), pay(pay)
    {
        if( !barrier || !pay )
            throw ModelException("BarrierMgrBase", "Barrier and/or pay is empty");

        if( pay->payOnce() )
        {
            IBarrierUtil::ICanHaveSticky* s = dynamic_cast<IBarrierUtil::ICanHaveSticky*>(barrier.get());
            if( !s || !s->isSticky() )
                throw ModelException("BarrierMgrBase", "If payOnce, must have sticky barrier");
        }
    }
    ~BarrierMgrBase(){};
    
    // create mapping between pay observe and sim dates
    // and for obs steps, assign corresponding payment date
    void preprocess(DateTime valueDt, DateTimeArrayConstSP simDts, const YieldCurve* discount)
    {
        static string method = "BarrierMgrBase::preprocess";
    
        valueDate = valueDt;
        simDates = simDts;

        // just in case barrier was not called
        barrier->preprocess(valueDt, simDts);
        const DateTimeArray& monDates = barrier->relevantMonDates();
        if( !DateTime::isSubset(*simDts, monDates) )         // sanity check
            throw ModelException(method, "barrier relevant monitor dates must be subset of sim dates");
        pay->preprocess(valueDt, simDts, monDates, discount);

        // create mapping to relevant monitor dates
        // make sure it's part of sim date in case of sticky (disallow mon dates between sim dates for now)
        // and part of mon date in case of non-sticky
        const DateTimeArray& payMonDts = pay->relevantMonDates();
        if( payMonDts.empty() )
            throw ModelException(method, "No pay relevant mon dates returned");
        IBarrierUtil::ICanHaveSticky* s = dynamic_cast<IBarrierUtil::ICanHaveSticky*>(barrier.get());
        if( s && s->isSticky() )
        {
            if( !DateTime::isSubset(*simDts, payMonDts) || payMonDts[0] < monDates[0] )
                throw ModelException(method, "pay observe dates must be subset of sim dates");
        } 
        else 
        {
            if( !DateTime::isSubset(monDates, payMonDts) )
                throw ModelException(method, "pay observe dates must be subset of barrier mon dates");
        }
        bool dummy;
        relMonDtMap = DateTime::createMapping(*simDts, payMonDts, dummy);
    }

    // loop through relevant dates to check barrier status at each date,
    // if hitValue=payIfHit, then
    //      if payOnce=true, only pay if hitValue changed from prev relevant date
    //      otherwise, pay given by value()
    // payoffs are cached in case it's used by value(). eg. step up coupon
    //
    // For bending, assume barrier bend period end dates are ascend order Ti, i=1 to n.
	// payoff within (Ti-1, Ti] uses hitValue for the i'th bend period
	// if payOnce=true, only pay if hit step is within (Ti-1, Ti] 
    //
    // To be consistent, if payOnce=true, must have payIfHit=true and sticky barrier
    virtual double calculate(DateTime valueDt, const IMCPathGenerator* pathGenIn)
    {
        double hitVal, hitValOld, payVal=0;

        int beginStep = pathGenIn->begin(0), endStep = pathGenIn->end(0); 
        IBarrierUtil::ICanHaveBend* bendView = dynamic_cast<IBarrierUtil::ICanHaveBend*>(barrier.get());
        IBarrierUtil::ICanHaveSticky* s = dynamic_cast<IBarrierUtil::ICanHaveSticky*>(barrier.get());

        // if payOnce, need to go through each bending period, 
        // check that not already hit at end of last bend period
        //
        int i = beginStep, iOld = i-1;
        i += relMonDtMap[i];
        if( pay->payOnce() )
        {
            // we know that barrier is sticky
            for(; i<endStep; i++, i+=relMonDtMap[i])
            {
                hitValOld = iOld<0?s->initialHitVal():bendView?bendView->hitValueBend(iOld, i):barrier->hitValue(iOld);
                hitVal = bendView?bendView->hitValueBend(i, i):barrier->hitValue(i);
                payVal +=  (hitVal - hitValOld) * pay->value(i, valueDt, pathGenIn);

                if( Maths::isZero(pay->payIfHitOne()?(hitVal-1.0):hitVal))
                {
                    i = s->periodEndStep(i);
                }
                iOld = i;
            }

            if( !pay->payIfHitOne() ) payVal = -payVal;

        }
        else
        {
            for(; i<endStep; i++, i+=relMonDtMap[i])
            {
                hitVal = bendView?bendView->hitValueBend(i, i):barrier->hitValue(i);
                hitVal = pay->payIfHitOne()?hitVal:(1.0-hitVal);
                payVal +=  hitVal * pay->value(i, valueDt, pathGenIn);
            }
        }

        return payVal;
    }
    
    // obtain knowCF and phyD from historical hits. ignore bending since it's event handling
    void getKnowCFPhyD(DateTime valueDate, const IMCPathGenerator* pathGenIn, double scalingFactor,
                       CashFlowArray& knownCFs, PhysicalDeliveryArray& phyDs)
    {
        static string method = "BarrierMgrBase::getKnowCFPhyD";

        knownCFs.clear();
        phyDs.clear();
        if( !pathGenIn->doingPast() )
            throw ModelException(method, "Internal error. Can only call if doingPast");

        double hitVal;
        bool   noRelevantDate = true;       // to avoid UMR of hitVal.
        int i, beginStep = pathGenIn->begin(0), endStep = pathGenIn->end(0);
        for(i=beginStep, i+=relMonDtMap[i]; i<endStep; i++, i+=relMonDtMap[i])
        {
            noRelevantDate = false;
            hitVal = barrier->hitValue(i);
            
            // historical hit value must be binary
            if( !Maths::isZero(hitVal) && !Maths::isZero(hitVal-1.0) )
                throw ModelException(method, "Internal error. Historical hitValue must be 0 or 1");
            
            if( Maths::isZero(pay->payIfHitOne()?(hitVal-1.0):hitVal) )
            {
                pay->collect(i, pathGenIn, scalingFactor, knownCFs, phyDs);

                if( pay->payOnce() )
                    return; // terminate
            }
        }

        // if sticky and not payOnce, we roll forward to get all FUTURE knownCF/PhyD as well, 
        // even beyond endStep
        IBarrierUtil::ICanHaveSticky* s = dynamic_cast<IBarrierUtil::ICanHaveSticky*>(barrier.get());
        if (!noRelevantDate){
            if (Maths::isZero(pay->payIfHitOne()?(hitVal-1.0):hitVal) &&
                s && s->isSticky() && !Maths::isZero(hitVal - s->initialHitVal()) )
            {
                for(; i<simDates->size(); i++, i+=relMonDtMap[i])
                    pay->collect(i, pathGenIn, scalingFactor, knownCFs, phyDs);
            }
        }
    }


protected:
    // barrier. also get simDate and path from.
    IBarrierUtilSP          barrier;
    IBarrierPaySP           pay;
        
    // *** transient field ***
    DateTime                valueDate;
    DateTimeArrayConstSP    simDates;

    // array with same size as sim dates to map sim step to dates when pay off is determined
    IntArray                relMonDtMap;

};

typedef refCountPtr<BarrierMgrBase> BarrierMgrBaseSP;


/****************************************************************************/
//
// BarrierFactory functions to convert barriers to sticky/reverser/pay
//
/****************************************************************************/

IBarrierUtilSP BarrierFactory::makeSticky(IBarrierUtilSP barrier)
{
    IBarrierUtil::ICanHaveSticky* s = dynamic_cast<IBarrierUtil::ICanHaveSticky*>(barrier.get());
    if( s && s->isSticky() )
        return barrier;

    // if is reverser, need to make the component sticky, then get reverser
    BarrierReverser *r = dynamic_cast<BarrierReverser*>(barrier.get());
    if( r )
        return makeReverser(makeSticky(r->barrier));

    if( dynamic_cast<IBarrierUtil::ICanHaveBend*>(barrier.get()) )
        return IBarrierUtilSP(new BarrierStickyBend(barrier));
    else
        return IBarrierUtilSP(new BarrierSticky(barrier));
}

IBarrierUtilSP BarrierFactory::makeCoupler(vector<IBarrierUtilSP> barriers, int numHit, bool isOut, bool isUp)
{
    for(unsigned int i=0; i<barriers.size(); i++)
        if( dynamic_cast<IBarrierUtil::ICanHaveBend*>(barriers[i].get()) )
            return IBarrierUtilSP(new BarrierCouplerBend(barriers, numHit,isOut, isUp));

    return IBarrierUtilSP(new BarrierCoupler(barriers, numHit, isOut, isUp));
}

IBarrierUtilSP BarrierFactory::makeReverser(IBarrierUtilSP barrier)
{
    // if is reverser, return the component
    BarrierReverser *r = dynamic_cast<BarrierReverser*>(barrier.get());
    if( r )
        return r->barrier;

    if( dynamic_cast<IBarrierUtil::ICanHaveBend*>(barrier.get()) )
        return IBarrierUtilSP(new BarrierReverserBend(barrier));
    else
        return IBarrierUtilSP(new BarrierReverser(barrier));
}

IBarrierPaySP BarrierFactory::makeTrivialPay(const DateTimeArray& obsDts)
{
    return IBarrierPaySP(new TrivialBarrierPay(obsDts));
}

IBarrierMgrSP BarrierFactory::makeManager(IBarrierUtilSP barrier, IBarrierPaySP pay)
{
    return IBarrierMgrSP(new BarrierMgrBase(barrier, pay));
}

DRLIB_END_NAMESPACE

