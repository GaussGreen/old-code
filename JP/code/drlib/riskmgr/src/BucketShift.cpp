//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BucketShift.cpp
//
//   Description : Bucket shift interface 
//
//   Author      : Stephen Hope
//
//   Date        : 14 Mar 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ExpiryResult.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/BenchmarkDate.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/MuSpecial.hpp"

DRLIB_BEGIN_NAMESPACE

// a struct to help with sorting expiry objects
struct dateSortObj
{
    DateTime date;
    int      arrayIdx;
};

static int dateCompare(
    const void* ptr1,   /*(I) Pointer to first date to compare */
    const void* ptr2)   /*(I) Pointer to second date to compare */
                    
{
    
    DateTime date1 = ((dateSortObj*)ptr1)->date;      
    DateTime date2 = ((dateSortObj*)ptr2)->date;

    if (date1.isLess(date2))
    {
        return(-1);
    }
    else if (date1.isGreater(date2))
    {
        return(1);
    }
    else
    {
        return(0);
    }
}

/** Returns the expiries which are to be tweaked */
ExpiryArrayConstSP BucketShift::calculateExpiries()
{
    ExpiryArraySP filteredExpiries(new ExpiryArray(0));
    for (int i = 0; i < expiries->size(); ++i) {
        if ((*expiries)[i].get()) {
            filteredExpiries->push_back((*expiries)[i]);
        }
    }
    return filteredExpiries;
}

/** Returns the expiry which is currently being tweaked */
ExpiryConstSP BucketShift::getExpiry()const 
{
    if (!expiry) {
       throw ModelException("BucketShift::getExpiry", "Expiry is Null");
    }
    return expiry;
}

/** Returns the expiries which are to be tweaked */
ExpiryArrayConstSP BucketShift::getExpiries()const{
    return expiries;
}

/** returns a sorted expiryArray which excludes any past date */
void BucketShift::getRevisedExpiries(const DateTime& valueDate,
                                     ExpiryArraySP& revisedList,
                                     ExpiryArraySP& pastDateList)
{
    ExpiryArrayConstSP expiries = calculateExpiries();
    int i;

    dateSortObj* sortArray = new dateSortObj[expiries->size()];
    for (i=0 ; i<expiries->size() ; ++i)
    {
        sortArray[i].date     = (*expiries)[i]->toDate(valueDate);
        sortArray[i].arrayIdx = i;
    }

    qsort((void *)(sortArray),
                   expiries->size(),
                   sizeof(dateSortObj),
                   dateCompare);

    revisedList  = ExpiryArraySP(new ExpiryArray(0));
    pastDateList = ExpiryArraySP(new ExpiryArray(0));

    DateTime currentDate(0,0);
    DateTime previousDate(0,0);
    for (i=0 ; i<expiries->size() ; ++i)
    {
        currentDate = (*expiries)[sortArray[i].arrayIdx]->toDate(valueDate);
        if (  currentDate.isGreater(valueDate) && 
             !currentDate.equals(previousDate) )
        {
            revisedList->push_back((*expiries)[sortArray[i].arrayIdx]);
        }
        else if ( !currentDate.isGreater(valueDate)  &&
                  !currentDate.equals(previousDate)  )
        {
            pastDateList->push_back((*expiries)[sortArray[i].arrayIdx]);
        }
                   
        previousDate = currentDate;
    }

    delete [] sortArray;
    // return revisedList;
}

BucketShift::BucketShift(const CClassConstSP& clazz,
                         const string&        outputName,
                         const double&        shiftSize,
                         const ExpiryArray*   expiries): 
    ScalarShift(clazz, outputName, shiftSize), expiries(copy(expiries)){}

/** implements a one sided bucket derivative for each instance of the
    market data which is sensitive to this SensControl - same as vector 
    shift with some additions to set up the front and back stubs */
void BucketShift::calculate(TweakGroup*  tweakGroup,
                            CResults*    results){
    try{
        // get list of names to calculate result for
        OutputNameArrayConstSP names(this->names(tweakGroup, results));

        if (names->empty()) {
            results->storeNotApplicable(this);            
        }
    
        // see if the instrument has a last sens date method
        LastSensDate* lsd = dynamic_cast<LastSensDate*>(tweakGroup->
                                                        getInstrument());
        DateTime      endDate;
        DateTime      valueDate;
        int           idx, jdx;

        valueDate = tweakGroup->getInstrument()->getValueDate();

        ExpiryArraySP expiries;
        ExpiryArraySP pastExpiries;


        getRevisedExpiries(valueDate, expiries, pastExpiries);

        // add lastSensDate if it's after the last bucket date
        if ( lsd ) {
            endDate   = lsd->endDate(this);
            if ( !MuSpecial::TYPE->isInstance(this) ) {
                if ( expiries->size() == 0 || 
                     expiries->back()->toDate(valueDate).isLess(endDate))
                {
                    ExpirySP expyDate(new BenchmarkDate(endDate));
                    expiries->push_back(expyDate);
                }
            }
        }

        if ( expiries->size() > 0 || pastExpiries->size() > 0) {
            for (idx = 0; idx < names->size(); idx++){
                // store the name of what we want to shift
                setMarketDataName((*names)[idx]);
                /* skip over where result has been calculated already */
                if (!results->exists(this)){

                    try {
                        ExpiryResultArraySP tweaks(new ExpiryResultArray(0));
                        for (int pdx = 0; pdx < pastExpiries->size(); pdx++){
                            // store 0.0 result in array
                            tweaks->push_back(
                                ExpiryResult((*pastExpiries)[pdx], 0.0));
                        }

                        if ( expiries->size() > 0 ) {
                            // convert the expiries to dates
                            DateTimeArray bucketDates(0);
                            if (valueDate.isLess((*expiries)[0]->
                                                 toDate(valueDate)))
                            {
                                bucketDates.push_back(valueDate);
                            }

                            int expyIdx;
                            for ( expyIdx=0 ; expyIdx < expiries->size() ;
                                  ++expyIdx )
                            {
                                DateTime nextExpiryDate =
                                    (*expiries)[expyIdx]->toDate(valueDate);
                                bucketDates.push_back(nextExpiryDate);
                            }

                            if ( lsd )
                            {
                                DateTime lastExpiryDate(
                                    (*expiries)[expiries->size()-1]->
                                    toDate(valueDate));
                                if ( endDate.isGreater(lastExpiryDate))
                                {
                                    bucketDates.push_back(endDate);
                                }
                            }

                            // then loop over the expiries
                            bool expired = false;
                            for (jdx = 0; jdx < expiries->size(); jdx++){
                                // store the expiry which we want to tweak
                                expiry          = (*expiries)[jdx];
                                bucketStartDate = bucketDates[jdx];
                                bucketEndDate   = bucketDates[jdx+1];

                                // calculate sens (if not expired)
                                double firstDeriv = 0.0;
                                if (!expired) {
                                    firstDeriv = 
                                        calcOneSidedFirstDeriv(tweakGroup, 
                                                               results);
                                }                           
                                // store result in array
                                tweaks->push_back(ExpiryResult(expiry, 
                                                               firstDeriv));

                                // do we need to tweak anymore ?
                                if (lsd) {
                                    expired = expiry->
                                        toDate(valueDate).isGreater(endDate);
                                }
                            }
                        }
                        // and store it
                        results->storeGreek(tweaks, this);
                    }
                    catch (exception& e) {
                        results->storeGreek(IObjectSP(new Untweakable(e)),
                                            this);
                    }
                }
            }
        }

    } catch (exception& e){
        throw ModelException(&e,  "BucketShift::calculate");
    }
}

/** returns the start- and end-dates for the current bucket */
void BucketShift::getBucketDates(DateTime& startDate, 
                                 DateTime& endDate) const
{
    startDate = bucketStartDate;
    endDate   = bucketEndDate;
}

/** for reflection */
BucketShift::BucketShift(const CClassConstSP& clazz,
                         const string&        sensName):
    ScalarShift(clazz, sensName){}

class BucketShiftHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(BucketShift, clazz);
        SUPERCLASS(ScalarShift);
        FIELD_NO_DESC(expiry);
        FIELD_MAKE_TRANSIENT(expiry);
        FIELD(expiries, "Expiry dates");
        FIELD_NO_DESC(bucketStartDate);
        FIELD_MAKE_TRANSIENT(bucketStartDate);
        FIELD_NO_DESC(bucketEndDate);
        FIELD_MAKE_TRANSIENT(bucketEndDate);
    }
};

CClassConstSP const BucketShift::TYPE = CClass::registerClassLoadMethod(
    "BucketShift", typeid(BucketShift), BucketShiftHelper::load);

DRLIB_END_NAMESPACE
