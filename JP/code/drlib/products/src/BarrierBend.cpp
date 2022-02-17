//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BarrierBend.cpp
//
//   Description :
//
//   Date        : Q Hou. Nov 2004
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/BarrierBend.hpp"
#include "edginc/Barrier.hpp"

DRLIB_BEGIN_NAMESPACE


void IBarrierBendMaker::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER_INTERFACE(IBarrierBendMaker, clazz);
    EXTENDS(IObject);
}

CClassConstSP const IBarrierBendMaker::TYPE = CClass::registerInterfaceLoadMethod(
    "IBarrierBendMaker", typeid(IBarrierBendMaker), load);



////////////////////////////////////////////////////
//
// ScheduleBarrierBend and its maker
//
////////////////////////////////////////////////////

class ScheduleBarrierBendMaker;

class ScheduleBarrierBend : virtual public IBarrierBend
{
public:
    // basically data holder for barrier bending info for MC
    class BendInfo : public CObject
    {
    public:
        static CClassConstSP const TYPE;
        BendInfo();
        ~BendInfo(){};

        int getNbBendAmt() const; // return number of bending amounts

        bool isTrivial() const;

        static void load(CClassSP& clazz);
        static IObject *defaultBendInfo();

        BendInfo(CClassConstSP& clazz);
        BendInfo(const BendInfo & clazz);

    public:
        DateTime    startDate; // $unregistered
        DateTime    endDate; // $unregistered

        int         startIdx;       // index of start date in the simDates $unregistered
        int         endIdx;         // index of end date in the simDates $unregistered
        DoubleArray bendAmts;       // array of size (endIdx - startIdx + 1) $unregistered
    };

    typedef smartPtr<BendInfo> BendInfoSP;
    typedef array<BendInfoSP, BendInfo> BendInfoArray;

public:
    ~ScheduleBarrierBend(){};

    // methods to create barrier bending info
    ScheduleBarrierBend(
        const ScheduleBarrierBendMaker *bend,
        const DateTime&         valueDate,
        const DateTimeArray&    simDates);

    BendInfoSP getBendInfo(DateTime myEndDate);

    int getMaxOverlap() const;

    bool isTrivial() const;

    virtual double getBending(int bendDateIdx,
                              int simDateIdx) const;

    virtual bool getBending(int bendDateIdx,
                           DoubleArray &bendAmt) const;

    virtual IntArray getBendingEndByStart(int startDateIdx) const;

    virtual IntArray getBendingEndByCover(int startDateIdx) const;

    virtual int getBendingStartByEnd(int endDateIdx) const;

    virtual IntArray getBendingStart() const;

    virtual IntArray getBendingEnd() const;

public:
    BendInfoArray   bendInfos;
};


class ScheduleBarrierBendMaker : public CObject,
                            virtual public IBarrierBendMaker {
public:
    static CClassConstSP const TYPE;
    friend class ScheduleBarrierBend;
    friend class ScheduleBarrierBendMakerHelper;

    virtual IBarrierBendSP getBarrierBend(const DateTime& valueDate, const DateTimeArray& simDates) const;

    virtual bool trivial() const;

    virtual void validatePop2Object();

protected:
    // schedule to specify bending amounts and length of bending periods
    DateTimeArray       dates;
    DoubleArray         amounts;
    DoubleArray         lengths;
    bool                startAscending;
    bool                useHist;
private:
    ScheduleBarrierBendMaker(); // for reflection
    ScheduleBarrierBendMaker(const ScheduleBarrierBendMaker& rhs); // not implemented
    ScheduleBarrierBendMaker& operator=(const ScheduleBarrierBendMaker& rhs); // not implemented

};


ScheduleBarrierBend::ScheduleBarrierBend(
        const ScheduleBarrierBendMaker *bend,
        const DateTime&         valueDate,
        const DateTimeArray&    simDates)
{
    static string method = "ScheduleBarrierBend::ScheduleBarrierBend";

    int nbObs = simDates.size();
    if( nbObs == 0 )
        throw ModelException(method, "No sim dates passed in");

    // no need to do any thing if no bend date is supplied
    const DateTimeArray &bendDates = bend->dates;
    int nbBendInfo = bendDates.size();
    if( nbBendInfo == 0 ) return;

    DateTime::ensureIncreasing(simDates, method + "sim dates must in strict ascending order", true);
    DateTime::ensureIncreasing(bendDates, method + "bend dates must in strict ascending order", true);
    if( !DateTime::isSubset(simDates, bendDates) )
        throw ModelException(method, "bend dates must be subset of sim dates");

    // find out location of bendDates inside sim dates for later lookup
    bool dummy;
    IntArray endMap = DateTime::createMapping(simDates, bendDates, dummy);

    // allocate memory for the holding area for calculations
    // calculate bending amounts for the immediate horizon before each required dates
    DoubleArray yrFrac(nbObs);
    yrFrac[0] = 0;
    int i, j;
    for(i=1; i<nbObs; i++)
        yrFrac[i] = simDates[0].yearFrac( simDates[i] );

    bendInfos.clear();
    for(i=0, j=endMap[0]; j<nbObs; i++, j++, j+=endMap[j])
    {
        // if bend does not apply to historical dates, ignore fully historical bending period
        if( !bend->useHist && bendDates[i] <= valueDate )
            continue;

        // allocate memory
        BendInfoSP bendInfo(new BendInfo);

        // interpolate to find the bending info for the bend date required
        // may not fail if can't interpolate bending info
        double bendAmt=0, bendLen=0;
        bendAmt = bend->amounts[i];
        bendLen = bend->lengths[i]/Barrier::BUS_DAYS_IN_YEAR_BA; //convert into year fraction

        // ignore trivial bendings
        if( Maths::isZero(bendAmt) )
            continue;

        // find position of the bend date, ie. end inddex, within the observation dates
        bendInfo->endDate = simDates[j];
        bendInfo->endIdx = j;

        // find position of start index
        // if bend does not apply to historical dates, ignore historical portion of partial hist bend period
        double startT = yrFrac[j] - bendLen;
        while( j >=0 && yrFrac[j] >= startT  && (bend->useHist || simDates[j] >= valueDate) )
        {
            j--;
        }
        j++;
        bendInfo->startDate = simDates[j];
        bendInfo->startIdx = j;

        // calc linear interpolated bending amounts to be applied as multiplicative factor to barrier
        while( j <= bendInfo->endIdx )
        {
            bendInfo->bendAmts.push_back( 1.0 + bendAmt * (Maths::isZero(bendLen)?1.0:((yrFrac[j]-startT)/bendLen)) );
            j++;
        }
        j--;

        // may require the bended schedules must start in increasing order
        if( bend->startAscending && bendInfos.size() > 0 && bendInfo->startIdx < bendInfos.back()->startIdx )
            throw ModelException(method, "Barrier bend start dates must be of increasing order");

        bendInfos.push_back(bendInfo);
    }
}

// max number of overlapped bending period.
// used to determine nb of pricing slice (tree) or nb of monitor flags (MC)
int ScheduleBarrierBend::getMaxOverlap() const
{
    int i, maxDIdx = -1;
    for(i=0; i<bendInfos.size(); i++)
        if( maxDIdx < bendInfos[i]->getNbBendAmt() ) maxDIdx = bendInfos[i]->getNbBendAmt();
        return maxDIdx;
}

bool ScheduleBarrierBend::isTrivial() const
{
    return (bendInfos.size()==0);
}

double ScheduleBarrierBend::getBending(int endDateIdx,
                              int simDateIdx) const
{
    throw ModelException("Not implemented");
}

bool ScheduleBarrierBend::getBending(int endDateIdx,
                           DoubleArray &bendAmt) const
{
    for(int i=0; i<bendInfos.size(); i++)
        if( bendInfos[i]->endIdx == endDateIdx )
        {
            bendAmt = bendInfos[i]->bendAmts;
            return true;
        }
    return false;
}

IntArray ScheduleBarrierBend::getBendingEndByStart(int startDateIdx) const
{
    IntArray endIdxs;
    int i;
    for(i=0; i<bendInfos.size(); i++)
    {
        if( bendInfos[i]->startIdx == startDateIdx )
            endIdxs.push_back(bendInfos[i]->endIdx);
    }
    sort(endIdxs.begin(), endIdxs.end());
    return endIdxs;
}


IntArray ScheduleBarrierBend::getBendingEndByCover(int startDateIdx) const
{
    throw ModelException("Not implemented");
}

IntArray ScheduleBarrierBend::getBendingStart() const
{
    IntArray startIdxs(bendInfos.size());
    int i;
    for(i=0; i<bendInfos.size(); i++)
    {
        startIdxs[i] = bendInfos[i]->startIdx;
    }
    sort(startIdxs.begin(), startIdxs.end());
    // remove duplicate
    vector<int>::iterator iter = startIdxs.begin() + 1;
    while( iter < startIdxs.end() )
    {
        if ( *iter == *(iter-1) )
            iter = startIdxs.erase(iter);
        else
            ++iter;
    }
    return startIdxs;
}

int ScheduleBarrierBend::getBendingStartByEnd(int endDateIdx) const
{
    for(int i=0; i<bendInfos.size(); i++)
    {
        if( bendInfos[i]->endIdx == endDateIdx ) return bendInfos[i]->startIdx;
    }
    return -1;
}

IntArray ScheduleBarrierBend::getBendingEnd() const
{
    IntArray endIdxs(bendInfos.size());
    for(int i=0; i<bendInfos.size(); i++)
    {
        endIdxs[i] = bendInfos[i]->endIdx;
    }
    return endIdxs;
}


ScheduleBarrierBend::BendInfo::BendInfo()
: CObject(TYPE), startIdx(-1), endIdx(-1)
{}

ScheduleBarrierBend::BendInfo::BendInfo(CClassConstSP& clazz)
: CObject(clazz), startIdx(-1), endIdx(-1)
{}

int ScheduleBarrierBend::BendInfo::getNbBendAmt() const
{
    return bendAmts.size();
}

bool ScheduleBarrierBend::BendInfo::isTrivial() const
{
    return bendAmts.empty();
}

void ScheduleBarrierBend::BendInfo::load(CClassSP& clazz){
        clazz->setPrivate();
        REGISTER(ScheduleBarrierBend::BendInfo, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultBendInfo);
}

IObject *ScheduleBarrierBend::BendInfo::defaultBendInfo()
{
    return new BendInfo();
}


CClassConstSP const ScheduleBarrierBend::BendInfo::TYPE = CClass::registerClassLoadMethod(
    "ScheduleBarrierBend::BendInfo", typeid(ScheduleBarrierBend::BendInfo), load);

// initialise type for array of assets
// work around for msvc 7 bug
typedef ScheduleBarrierBend::BendInfoArray ScheduleBarrierBendBendInfoArray;
DEFINE_TEMPLATE_TYPE_WITH_NAME("ScheduleBarrierBend::BendInfoArray", ScheduleBarrierBendBendInfoArray);



ScheduleBarrierBendMaker::ScheduleBarrierBendMaker()
: CObject(TYPE), startAscending(true), useHist(false)
{}

IBarrierBendSP ScheduleBarrierBendMaker::getBarrierBend(const DateTime& valueDate, const DateTimeArray &simDates) const
{
    IBarrierBendSP bend;
    if( dates.size() && (dates.back() > valueDate || useHist) )
        bend = IBarrierBendSP(new ScheduleBarrierBend(this, valueDate, simDates));
    return bend;
}

bool ScheduleBarrierBendMaker::trivial() const
{
    return dates.size()==0;
}

void ScheduleBarrierBendMaker::validatePop2Object()
{
    static string method = "ScheduleBarrierBendMaker::validatePop2Object";
    try {
        if( dates.size() != amounts.size() ||
            dates.size() != lengths.size() )
            throw ModelException(method, "Dates, amount, lengths array mismatch");

        int i;
        for(i=0; i<dates.size(); i++)
        {
            if( i>0 && dates[i-1].isGreaterOrEqual(dates[i]) )
                throw ModelException(method, "Dates must be in strict ascend order");
            if( Maths::isNegative(lengths[i]) )
                throw ModelException(method, "Length must be non negative");
        }
    }
    catch (exception &e)
    {
        throw ModelException(e, method);
    }
}

class ScheduleBarrierBendMakerHelper{
public:
    static IObject* defaultScheduleBarrierBendMaker(){
        return new ScheduleBarrierBendMaker();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ScheduleBarrierBendMaker, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IBarrierBendMaker);
        FIELD(dates,     "bending schedule: dates for term structure");
        FIELD(amounts,   "bending schedule: amount of bending (multiplicative)");
        FIELD(lengths,   "bending schedule: length for bending (in days)");
        FIELD(startAscending,    "if enforce bent period starts are in ascending order. default is true");
        FIELD_MAKE_OPTIONAL(startAscending);
        FIELD(useHist,    "if bend historical portion of barrier. default is false");
        FIELD_MAKE_OPTIONAL(useHist);
        EMPTY_SHELL_METHOD(defaultScheduleBarrierBendMaker);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

CClassConstSP const ScheduleBarrierBendMaker::TYPE =
CClass::registerClassLoadMethod("ScheduleBarrierBendMaker",
                                typeid(ScheduleBarrierBendMaker), ScheduleBarrierBendMakerHelper::load);



DRLIB_END_NAMESPACE

