//   Filename    : SampleList.cpp 
//
//   Description : SampleList representation
//
//   Author      : Stephen Hope
//
//   Date        : 27 Feb 2001
//
//
//--------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_SAMPLELIST_CPP
#include "edginc/SampleList.hpp"
#include "edginc/VolProcessedBS.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SimSeries.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/Addin.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/MaturityTimePeriod.hpp"
#include "edginc/MultiFactors.hpp"

DRLIB_BEGIN_NAMESPACE

// this is about as good as we get from Pyramid 
#define EDG_CLOSE_ENOUGH_TO_ONE 1e-9

// PUBLIC methods

/** Constructor */
SampleList::SampleList(const DateTimeArray&  dates,
                       const DoubleArray&    values,
                       const DoubleArray&    weights):
    CObject(TYPE), dates(dates), values(values), weights(weights),
    startSimAtFirstAvgDate(false){
  //    try { // fix for gcc 2.95 on solaris opt
        validatePop2Object();
	//    } catch (exception&){
	//        throw;
	//    }
}
/* weight scaling shouldn't be greater than order of magnitude 10 (1000%) */
#define MAX_WEIGHT_SCALING 10.0

/** Create a sample using a generated array of dates formed from first and
    last date and rules for generating the intervening dates. The
    resultant dates (exluding the specified first and last dates)
    are then business date adjusted. */
SampleListSP SampleList::create(
    const DateTime&                valueDate, // to define historical samples
    const DateTime&                firstDate,
    const DateTime&                lastDate,
    const ExpiryConstSP&           matPeriod,
    bool                           stubAtEnd,
    const BadDayConventionConstSP& badDayConv,
    const HolidayConstSP&          hols,
    double                         histAvgToDate, // historic average to date
    double                         weightScaling) /* weights = weightScaling/
                                                     num samples*/
{
    static const string routine("SampleList::create");
    if (!Maths::isPositive(weightScaling)){
        throw ModelException(routine, "Weight scaling ("+
                             Format::toString(weightScaling)+
                             ") must be > 0.0");
    }
    if (weightScaling > MAX_WEIGHT_SCALING)
    {
        throw ModelException(routine, "Weight scaling factor ("+
                             Format::toString(weightScaling)+") must be a "
                             "decimal not a percentage (1.0 not 100.0%)");
    }
    
    
	DateTimeArraySP dates(generateDates(firstDate, lastDate, matPeriod,
                                        stubAtEnd, badDayConv, hols));
    SampleListSP sample(new SampleList());
    sample->dates = *dates;
    sample->values = DoubleArray(dates->size());
    sample->weights = DoubleArray(dates->size());
    for (int i = 0; i < dates->size(); i++){
        sample->weights[i] = weightScaling/dates->size();
        if (valueDate.isGreaterOrEqual((*dates)[i])){
            sample->values[i] = histAvgToDate;
        }
    }
    return sample;
}

void SampleList::GetMarket(const IModel* model, const CMarketDataSP market) {
}

/** Destructor */
SampleList::~SampleList()
{
    // empty
}

/** SampleList validation
    Conditions:
    1. dates, values or weights cannot be empty
    2. dates, values and weights arrays must be equal length
    3. dates must be in strictly increasing order
    4. weights and values must be positive */
void SampleList::validatePop2Object()
{
    static const string method = "SampleList::validatePop2Object";
    try {
        int numDates = dates.size();
        int i;
        // Condition 1
        if (numDates == 0 || values.size() == 0 || weights.size() == 0) {
            throw ModelException(method,
                                 "Sample list dates, values or weights arrays"
                                 " cannot be empty");
        }
        
        
        // Condition 2
        if (values.size() != numDates || weights.size() != numDates) {
            throw ModelException(method,
                                 "Sample list dates, values and weights "
                                 " arrays must be of equal length");
            
        }
        
        // Condition 3
        for (i=1; i < numDates; i++) {
            if (dates[i-1].isGreaterOrEqual(dates[i])) {
                throw ModelException(method,
                                     "sample list dates not in "
                                     "ascending order.\n"
                                     + dates[i-1].toString() +
                                     " occurs before " +
                                     dates[i].toString() +
                                     " in the dates list.");
            }
        }
        
        // Condition 4
        for (i = 0; i < numDates; i++) {
            if (Maths::isNegative(weights[i])) {
                throw ModelException(method,
                                     "Sample list weights cannot be "
                                     "negative (" +
                                     Format::toString(weights[i]) + ")");
            }
            if (Maths::isNegative(values[i])) {
                throw ModelException(method,
                                     "Sample list values cannot be "
                                     "negative (" +
                                     Format::toString(values[i]) + ")");
            }
        }
    }  catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** How to get a IMCStateVarGen (because we implement IRefLevel) */
IRefLevel::IStateVarGen* SampleList::createStateVarGen(
    const IMultiFactors* multiFactors,
    const DateTime&      valueDate) const{
    throw ModelException("SampleList::createStateVarGen", "Not yet done");
}


/** returns the weighted sum of samples so far */
double SampleList::sumToDate(const DateTime& today)const
{
    double sumSoFar = 0.0;
    int i = 0;
    
    while (i < dates.size() &&
           today.isGreaterOrEqual(dates[i])) {
        if (!Maths::isPositive(values[i])) {
            throw ModelException("SampleList::sumToDate",
                                 "sample at (" + 
                                 dates[i].toString() + 
                                 ") is historic [today is ("+
                                 today.toString() + 
                                 ")] and has not been set");
        }
        sumSoFar += values[i] * weights[i];
        i++;
    }

    return sumSoFar;
}

/** returns the average of the samples so far */
double SampleList::averageToDate(const DateTime& today)const
{
    double avgSoFar = 0.0;
    double wgtSoFar = 0.0;
    int i = 0;

    while (i < dates.size() &&
           today.isGreaterOrEqual(dates[i]))
    {
        if (!Maths::isPositive(values[i])) {
            throw ModelException("SampleList::averageToDate",
                                 "sample at (" + 
                                 dates[i].toString() + 
                                 ") is historic [today is ("+
                                 today.toString() + 
                                 ")] and has not been set");
        }

        avgSoFar += values[i] * weights[i];
        wgtSoFar += weights[i];
        i++;
    }

    if (!(Maths::isZero(wgtSoFar)))
    {
        avgSoFar = avgSoFar/wgtSoFar;
    }

    return avgSoFar;
}

/** returns the expected life of a sample:
    expected life = 
    [ sum for i = 1 to N (weight(i) * time(i))]
    -----------------------------------------
    [sum for i = 1 to N (weight(i)]
    
    where N is number of future weights. */
double SampleList::expectedLife(const DateTime& today)const
{
    static const string method = "SampleList::expectedLife";
    double weighting = 0.0;
    double expectedLife = 0.0;
    double years;

    int numFutureSamples = countFutureSamples(today);

    if (numFutureSamples > 0)
    {
        int i = dates.size() - numFutureSamples;
        
        while (i < dates.size())
        {
            years = today.yearFrac(dates[i]);
            expectedLife += weights[i] * years;
            weighting += weights[i];
            i++;
        }
        
        if (Maths::isZero(weighting))
        {
            throw ModelException(method,
                                 "sum of future weights must be > 0");
        }
        
        expectedLife /= weighting;
    }
    return expectedLife;
}

/** returns the expected date of a sample (cf expectedLife) */
DateTime SampleList::expectedEndDate(const DateTime& today)const
{
    DateTime expectedEndDate;
    int numFutureSamples = countFutureSamples(today);
    int i = 0;

    if (numFutureSamples > 0)
    {
        int numPastSamples = dates.size() - numFutureSamples;
        DateTimeArray futDates(numFutureSamples);
        DoubleArray futWeights(numFutureSamples);

        for (i = 0; i < numFutureSamples ; i++)
        {
            futDates[i]   = dates[i + numPastSamples];
            futWeights[i] = weights[i + numPastSamples];
        }
        expectedEndDate = DateTime::expectedFutureDate(today,
                                                       futDates,
                                                       futWeights);
    }
    else
    {
        // all samples are set 
        expectedEndDate = dates[dates.size() -1];
    }
    
    return expectedEndDate;
}

/** returns the number of samples remaining in the future */
int SampleList::countFutureSamples(const DateTime& today)const
{
    int count = 0;
    int i = 0;
    int historic = 0;

    while (i < dates.size() &&
           today.isGreaterOrEqual(dates[i]))
    {
        historic++;
        i++;
    }
    count = dates.size() - historic;

    return count;
}

/** returns the weighting of samples that remain in the future */
double SampleList::futureWeight(const DateTime& today)const
{
    double weight = 0.0;

    int numFutureSamples = countFutureSamples(today);

    int i = dates.size() - numFutureSamples;

    while (i < dates.size())
    {
        weight += weights[i];
        i++;
    }

    return weight;
}

/** Returns the first date in the sample list */
const DateTime& SampleList::getFirstDate() const{
    return dates[0];
}

/** Returns the first date in the sample list */
const DateTime& SampleList::getLastDate() const{
    return dates[dates.size()-1];
}

/** return the first and last dates in a sample list */
void SampleList::getBoundingDates(DateTime& loDate,
                                  DateTime& hiDate)const
{
    // validation ensure at least one date in SampleList
    loDate = dates[0];
    hiDate = dates[dates.size()-1];
}

/** returns TRUE if all the weights sum to 1.0 
    (within a sensible tolerance) ? */
bool SampleList::weightsSumToOne()const
{
    double weight = 0.0;
    int i = 0;

    for (i = 0 ; i < dates.size(); i++)
    {
        weight += weights[i];
    }

    return (Maths::areEqualWithinTol(weight, 1.0, EDG_CLOSE_ENOUGH_TO_ONE));
}

/** Calculates fwdPrices on dates[i] to dates[dates.size()-1] */
static void fwdPrice(const CAsset*        asset,
                     const DateTimeArray& dates,
                     int                  startingIndex,
                     CDoubleArray&        fwdPrices){
    int numFutureSamples = dates.size() - startingIndex;
    if (numFutureSamples > 0){
        /* estimate any future sample levels - note need to
           optimise call to forward price */
        if (startingIndex == 0){
            asset->fwdValue(dates, fwdPrices);
        } else {
            CDateTimeArray fwdDates(numFutureSamples);
            for (int j = 0; j < numFutureSamples; j++){
                fwdDates[j] = dates[startingIndex+j];
            }
            asset->fwdValue(fwdDates, fwdPrices);
        }
    }
}

/** returns an estimate of the weighted sum of future samples */
double SampleList::futureSampleSum(const CAsset* asset,
                                   const DateTime& today)const
{
    try{
        double sumSoFar = 0.0;
        int    numFutureSamples = countFutureSamples(today);
        if (numFutureSamples > 0){
            int futureSamplesOffset = dates.size() - numFutureSamples;
            /* estimate any future sample levels - note need to
               optimise call to forward price */
            CDoubleArray fwdPrices(numFutureSamples);
            fwdPrice(asset, dates, futureSamplesOffset, fwdPrices);

            for (int j = 0; j < numFutureSamples; j++){
                sumSoFar += fwdPrices[j] * weights[futureSamplesOffset+j];
            }
        }
        return sumSoFar;
    } catch (exception& e){
        throw ModelException(e, "SampleList::futureSampleSum");
    }
}

/** returns the expected average */
double SampleList::expectedAverage(const CAsset* asset,
                                   const DateTime& today)const
{
    static const string routine("SampleList::expectedAverage");
    double average = 0.0;
    double weight = 0.0;
    int    i;

    for (i = 0; i < dates.size() && today.isGreaterOrEqual(dates[i]); i++){
        if (!Maths::isPositive(values[i])) {
            string m("sample at (" +  dates[i].toString() + ") is historic "
                     "[today is ("+ today.toString() + ")] "
                     "and has not been set");
            throw ModelException(routine, m);
        }

        average += values[i] * weights[i];
        weight += weights[i];
    }

    // estimate any future sample levels
    int numFutureSamples = dates.size() - i;
    if (numFutureSamples > 0){
        try{
            /* estimate any future sample levels - note need to
               optimise call to forward price */
            CDoubleArray fwdPrices(numFutureSamples);
            fwdPrice(asset, dates, i, fwdPrices);
            
            for (int j = 0; j < numFutureSamples; j++){
                double thisWeight = weights[i+j];
                average += fwdPrices[j] * thisWeight;
                weight += thisWeight;
            }
        } catch(exception& e){
            throw ModelException(e,routine);
        }
    }
    // work out overall average 
    average = Maths::isZero(weight)? 0.0: (average/weight);
    return average;
}

/** returns the average mode level ( S * exp (mu*t - 1/2 * sigma ^2 t) 
// Not yet Tested.....
double SampleList::averageModeLevel(const CAsset* asset,
                                    const CVolProcessedBS* procVol,
                                    const DateTime& today)const
{
    static const string method("SampleList::averageModeLevel");
    try {
        int numFutureSamples = countFutureSamples(today);
        int numWeights = 0;
        double weight = 0.0;
        double variance = 0.0;
    
        double avgMode = 0.0;
        int    i, j;

        for (i = 0; i < dates.size() && today.isGreaterOrEqual(dates[i]); i++){
            if (!Maths::isPositive(values[i])) {
                string m("sample at (" +  dates[i].toString() + ") is historic "
                         "[today is ("+ today.toString() + ")] "
                         "and has not been set");
                throw ModelException(method, m);
            }

            avgMode += values[i] * weights[i];
            weight += weights[i];
        }

        // estimate any future average mode levels
        if (numFutureSamples > 0){
            CDoubleArray fwdPrices(numFutureSamples);
            fwdPrice(asset, dates, i, fwdPrices);
        
            int numPastSamples = dates.size() - numFutureSamples;
            DateTimeArray futDates(numFutureSamples+1);
            DoubleArray var(numFutureSamples);
            
            // variances will be calculated from today
            futDates[0] = today;
            for (i = 0; i < numFutureSamples; i++)
            {
                futDates[i+1] = dates[i + numPastSamples];
            }

            procVol->CalcVar(futDates,
                             procVol->fromFirst,
                             var);

            for (j = 0; j < numFutureSamples; j++){
                double thisWeight = weights[numPastSamples+j];
                avgMode += thisWeight * fwdPrices[j] * exp( - 0.5*var[j]);
                weight += thisWeight;
            }
        }
        return avgMode;
    } 
    catch(exception& e) {
        throw ModelException(e, method);
    }
}
*/

/** returns the average variance */
double SampleList::averageVariance(const CVolProcessedBS* procVol,
                                   const DateTime& today,
                                   bool  usePastWeights)const
{
    static const string method("SampleList::averageVariance");
    try {
        int numFutureSamples = countFutureSamples(today);
        int numWeights = 0;
        double weight = 0.0;
        double variance = 0.0;
    
        // interpolate variance at future samples 
        if (numFutureSamples > 0)
        {
            int numPastSamples = dates.size() - numFutureSamples;
            DateTimeArray futDates(numFutureSamples+1);
            DoubleArray var(numFutureSamples);
            int i = 0, j = 0;

            // variances will be calculated from today
            futDates[0] = today;
            for (i = 0; i < numFutureSamples; i++)
            {
                futDates[i+1] = dates[i + numPastSamples];
            }

            procVol->CalcVar(futDates,
                             procVol->fromFirst,
                             var);

            // we now have total variance from today to futDates[1], futDates[2] etc

            if (!usePastWeights)
            {
                weight = futureWeight(today);
                numWeights = numFutureSamples;
            }
            else
            {
                numWeights = dates.size();
                for (i = 0; i < dates.size(); i++)
                {
                    weight += weights[i];
                }
            }
        
            // see if all future weights are the same 
            bool equalWeights = true;
            i = 0;

            while (i < numFutureSamples && equalWeights)
            {
                if(!(Maths::equals(weights[i + numPastSamples],
                                   weights[numPastSamples])))
                {
                    equalWeights = false;
                }
                i++;
            }

            if (equalWeights)
            {
                for (i = 0; i < numFutureSamples; i++)
                {
                    variance += (2.0 * numFutureSamples - 2 * (i+1) + 1) *var[i];
                }
            
                variance /= numWeights * numWeights;
            }
            else
            {
                for (i = 0; i < numFutureSamples; i++)
                {
                    for (j = 0; j < numFutureSamples; j++)
                    {
                        variance += weights[i + numPastSamples] *
                            weights[j + numPastSamples] *
                            Maths::min(var[i], var[j]);
                    }
                }

                variance /= weight * weight;
            }
        }
        else
        {
            // all samples are set
            variance = 0.0;
        }
    
        return variance;
    } 
    catch(exception& e) {
        throw ModelException(e, method);
    }
}

/** returns the co-variance between two sample dates */
double SampleList::averageCovariance(const SampleList* other,
                                     const CVolProcessedBS* procVol,
                                     const DateTime& today)const
{
    static const string method("SampleList::averageCovariance");
    try {
        int numPastThis = 0, numPastOther = 0;
        int i = 0, j = 0;

    
        // get totals of weights in both samples 
        double weightThis = 0.0;
        for (i = 0; i < this->dates.size(); i++)
        {
            weightThis += this->weights[i];
        }

        double weightOther = 0.0;
        for (i = 0; i < other->dates.size(); i++)
        {
            weightOther += other->weights[i];
        }
    
        int numFutThis = this->countFutureSamples(today);
        DoubleArray varThis(numFutThis);

        // interpolate variance at future samples
        if (numFutThis > 0)
        {
            numPastThis = this->dates.size() - numFutThis;
            DateTimeArray futDatesThis(numFutThis+1);
            futDatesThis[0] = today;
        
            for (i = 0; i < numFutThis; i++)
            {
                futDatesThis[i+1] = this->dates[i + numPastThis];
            }

            procVol->CalcVar(futDatesThis,
                             procVol->fromFirst,
                             varThis);
        }

        int numFutOther = other->countFutureSamples(today);
        DoubleArray varOther(numFutOther);

        // interpolate variance at future samples
        if (numFutOther > 0)
        {
            numPastOther = other->dates.size() - numFutOther;
            DateTimeArray futDatesOther(numFutOther+1);
            futDatesOther[0] = today;
        
            for (i = 0; i < numFutOther; i++)
            {
                futDatesOther[i+1] = other->dates[i + numPastOther];
            }

            procVol->CalcVar(futDatesOther,
                             procVol->fromFirst,
                             varOther);
        }

        double covariance = 0.0;

        for (i = 0; i < numFutThis; i++)
        {
            for (j = 0; j < numFutOther; j++)
            {
                covariance += this->weights[i+numPastThis] *
                    (other->weights[j+numPastOther]) *
                    Maths::min(varThis[i], varOther[j])/(weightThis * weightOther);
            }
        }
    
        return covariance;
    } 
    catch(exception& e) {
        throw ModelException(e, method);
    }
}

bool SampleList::getPreviousSample(
    const DateTime &today, 
    DateTime&       previousDate, 
    double &        value, 
    double &        weight) const
{
    int idx=dates.size()-countFutureSamples(today)-1;
    if((idx<0)||(idx>dates.size()))
    {
        return false;
    }
    else
    {
        previousDate=dates[idx];
        value = values[idx];
        weight= weights[idx];
        return true;
    }
}

/** populate future samples up to 'rollDate' with either the spot or the
    forward price of an asset
*/
void SampleList::roll(const Asset*    asset,
                      const DateTime& today,
                      const DateTime& rollDate,
                      bool            useSpot) {
    static const string method = "SampleList::roll";
    try {
        double spot=0.0;

        if (useSpot) {
            spot = asset->getSpot();
        }

        bool pastRollDate = false;
        int  i = 0;
        while ( !pastRollDate && i<dates.size() ) {
            if ( ( dates[i].isGreater(today) && !dates[i].isGreater(rollDate)) ||
                 ( dates[i].equals(today)    && Maths::isZero(values[i]))    ) {
                if (useSpot) {
                    values[i] = spot;
                }
                else {
                    values[i] = asset->fwdValue(dates[i]);
                }
            } else if ( dates[i].isGreater(rollDate) ) {
                pastRollDate = true;
            }
            ++i;
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Rolls the PastValues through time populating any samples as they
    become historic */
void SampleList::roll(const Theta::Util&         thetaUtil,
                      int                        iAsset,
                      const Asset*               asset){
    const DateTime& origDate = thetaUtil.getOriginalValueDate();
    const DateTime& newDate = thetaUtil.getNewValueDate();
    roll(asset, origDate, newDate, !thetaUtil.useAssetFwds());
}

/** rolls past values for all assets inside a MultiFactor */
void SampleList::roll(const Theta::Util&           thetaUtil,
                      const IMultiMarketFactors*   multiFactor){
    static const string method = "SampleList::roll";
    try {
        const DateTime& today = thetaUtil.getOriginalValueDate();
        const DateTime& rollDate = thetaUtil.getNewValueDate();
        IMarketFactorConstSP factor(multiFactor->getFactor(0));
        if (!IGeneralAsset::TYPE->isInstance(factor)){
            throw ModelException(method, "Past values only support for "
                                 "IGeneralAssets");
        }
        const IGeneralAsset* asset = DYNAMIC_CAST(IGeneralAsset, factor.get());
        bool pastRollDate = false;
        int  i = 0;
        while ( !pastRollDate && i < dates.size()) {
            if ((dates[i].isGreater(today) && !dates[i].isGreater(rollDate)) ||
                (dates[i].equals(today)    && Maths::isZero(values[i]))) {
                double newSpot = thetaUtil.useAssetFwds()?
                    asset->fwdValue(dates[i]): asset->getSpot();
                values[i] = newSpot;
            } else if ( dates[i].isGreater(rollDate) ) {
                pastRollDate = true;
            }
            ++i;
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

void SampleList::validate(const DateTime& today,
                          const DateTimeArrayArray& requiredDates,
                          IMultiFactors* assets,
                          SamplingConvention* dateAdjust) {
    // nothing for now
    return;
}

// Retrieves all past samples used by an instrument
// For past dates the provided samples are used and if missing
// samples are retrieved from the centralised asset history
void SampleList::pastSamplesEvents(const DateTime& today,
                                   const DateTimeArrayArray& requiredDates,
                                   IMultiFactors* assets,
                                   SamplingConvention* dateAdjust,
                                   EventResults* events) {
}


/** populate samples between date1 < date2, assuming that the spot level is
	value1 and value2 on these dates, and using linear interpolation. */
void SampleList::rollLinear(const DateTime& date1, double value1,
				const DateTime& date2, double value2)
{
	if (date1 >= date2)
	{
		throw ModelException("SampleList::rollLinear",
							 "date1 must be < date2");
	}

	bool pastDate2 = false;
	int  i = 0;
    while ( !pastDate2 && i < dates.size()) {
		// note we only compare the date for date2, so we will populate
		// an EOD sample on date2.
        if ( dates[i] > date1 && dates[i].getDate() <= date2.getDate() )
		{
			values[i] = value1 + 
				(double)(dates[i].getDate() - date1.getDate()) / (double)(date2.getDate() - date1.getDate()) *
				(value2 - value1);

		} else if ( dates[i].getDate() > date2.getDate() ) {
            pastDate2 = true;
        }
        i++;
    }
}


/** Returns the number of assets for which this object has data for */
int SampleList::getNumAssets() const{
    return 1;
}

/** Returns the dates in the sample - needed until the MC supports
    averaging */
const DateTimeArray& SampleList::getDates() const{
    return dates;
}

/** Returns the weights in the sample. Need to review interaction
    of MC with this class */
const DoubleArray& SampleList::getWeights() const{
    return weights;
}

/** returns simulation dates (of all assets) which are
    strictly in the future */
DateTimeArray SampleList::getFutureDates(const DateTime& valueDate) const{
    return  valueDate.getFutureDates(dates);
}

/** returns all the simulation dates (of all assets) - here there is only
    one asset */
const DateTimeArray& SampleList::getAllDates() const{
    return dates;
}
    
/** Returns all the simulation dates of the given asset */
const DateTimeArray& SampleList::getDates(int iAsset) const{
    return dates;
}

const DoubleArray& SampleList::getValues() const{
    return values;
}

/** returns the number of future dates for given asset */
int SampleList::numFutureDates(const DateTime& valueDate,
                               int             iAsset) const{
    return valueDate.numFutureDates(dates);
}
    
/** returns the number of dates (past and future) for given asset */
int SampleList::numDates(int iAsset) const{
    return dates.size();
}

/** returns the date when a model should start simulating - will not
    return a date before today */
const DateTime& SampleList::getSimStartDate(const DateTime& today) const{
    if (startSimAtFirstAvgDate){
        // start simulation on start date
        return (dates.front().isGreater(today)? dates.front(): today);
    }
    // otherwise start simulation immediately - not at first average date
    return today;
}

/** Returns true if each asset has the same set of simulation dates */
bool  SampleList::sameDatesPerAsset() const{
    return true;
}

/** Generate an array of dates given a first and last date and rules
    for generating the intervening dates. The resultant dates
    (exluding the specified first and last dates) are then business
    date adjusted. Not sure where this method belongs but can go here
    for now */
DateTimeArraySP SampleList::generateDates(
    const DateTime&                firstDate,
    const DateTime&                lastDate,
    const ExpiryConstSP&           expiry,
    bool                           stubAtEnd,
    const BadDayConventionConstSP& badDayConv,
    const HolidayConstSP&          hols){
    static const string routine("SampleList::generateDates");
    // special handling for trivial case
    if (firstDate.equals(lastDate)){
        return DateTimeArraySP(new DateTimeArray(1, firstDate));
    }
    if (firstDate.isGreaterOrEqual(lastDate)){
        throw ModelException(routine, "First date ("+
                             firstDate.toString()+") must be before last "
                             "date ("+lastDate.toString()+")");
    }

    // swap tool could do with being upgraded to take an Expiry
    int count;
    string interval;
    const MaturityPeriod*     matPeriod = 0;
    const MaturityTimePeriod* matTPeriod = 0;
    if (MaturityPeriod::TYPE->isInstance(expiry)){
        matPeriod = dynamic_cast<const MaturityPeriod*>(expiry.get());
    } else if (MaturityTimePeriod::TYPE->isInstance(expiry)){
        matTPeriod = dynamic_cast<const MaturityTimePeriod*>(expiry.get());
        matPeriod = matTPeriod->getMaturityPeriod().get();
    } else {
        throw ModelException(routine, "Expiry of type "+
                             expiry->getClass()->getName()+" is not supported."
                             "Use MaturityPeriod or MaturityTimePeriod");
    }
        
    matPeriod->decompose(count, interval);
    // first generate unmodified dates
    DateTimeArraySP dates(SwapTool::dateArray(firstDate, lastDate,
                                              count, interval,
                                              stubAtEnd));
    // then bad date adjust - excluding first and last dates
    badDayConv->adjust(dates->begin()+1, dates->end()-1, hols.get());
    // ensure all dates are within bounds
    DateTime::removeOutliers(*dates, firstDate, lastDate);
    // then remove duplicates
    DateTime::removeDuplicates(*dates, true /* ignore time of day */);
    // finally force time of day if specified
    if (matTPeriod){
        int time = matTPeriod->getTime();
        for (int i = 1; i < dates->size()-1; i++){
            if ((*dates)[i].getTime() != time){
                (*dates)[i] = DateTime((*dates)[i].getDate(), time);
            }
        }
    }
    return dates;
}

class SampleList::MCPath: virtual public IRefLevel::IMCPath{
private:
    const SampleList*  sample;
    int                numFutureValues; 
    double             sumSoFar; 
    IntArray           trivialIndex;
public:
    // constructor for when we start simulation at first avg date and
    // value date is before first avg date
    MCPath(const SampleList* sample, const DoubleArray& spotsAtSimStart):
        sample(sample), 
        numFutureValues(sample->dates.size()-1), // value at first date is known
        trivialIndex(1, 0){
        sumSoFar = spotsAtSimStart[0] * sample->weights[0];
    }

    // constructor for when we start simulation today
    MCPath(const SampleList*  sample, const DateTime& valueDate):
        sample(sample), numFutureValues(sample->numFutureDates(valueDate, 0)),
        sumSoFar(0), trivialIndex(1, 0){
        int numPastValues =  sample->dates.size() - numFutureValues;
        // calculate sumSoFar
        for (int i = 0; i < numPastValues; i++){
            sumSoFar += sample->values[i] * sample->weights[i];
        }
    }
       
    /** Returns an array of integers denoting which assets have a simulation
        date given the index into getAllDates(). The number of the assets
        (0, 1, 2, ..., numAssets-1) is the same as the order in which the
        series is supplied to the constructor */
    const IntArray& assetsOnDate(int index) const{
        return trivialIndex;
    }

    /** Main method. Returns refLevel given values for future path of
        specified asset. The array must run from the next [future]
        date in the RefLevel to the last date */
    virtual double refLevel(int            iAsset,
                            const double*  futurePath) const{
        double sum = sumSoFar;
        int numPastValues =  sample->dates.size() - numFutureValues;
        for (int i = 0; i < numFutureValues; i++){
            sum += futurePath[i] * sample->weights[numPastValues+i];
        }
        return sum;
    }

    /** Used for estimating refLevel when future path depends upon
        refLevel (ie vol interp level for LN) */
    virtual double refLevel(int                  iAsset,
                            const IMultiFactors* multiFactor) const{
        double sum = sumSoFar;
        if (numFutureValues > 0){
            // extrapolate just using spot
            double spot = multiFactor->assetGetSpot(0);
            int numPastValues =  sample->dates.size() - numFutureValues;
            double remainingWeight = 0.0;
            for (int i = 0; i < numFutureValues; i++){
                remainingWeight += sample->weights[numPastValues+i];
            }
            sum += spot * remainingWeight;
        }
        return sum;
    }            
        
    //// utility method - returns the number of future points
    virtual int numFutureDates(int iAsset) const{
        return numFutureValues;
    }
};
    

/** Creates an object implementing the IMCPath interface */
IRefLevel::IMCPath* SampleList::createMCPath(
    const DateTime&    valueDate,
    const DoubleArray& spotsAtSimStart, // not needed since sim start = today
    const IPastValues* pastValues) const{ // ignored here as we have them
    if (startSimAtFirstAvgDate && valueDate.isLess(dates.front())){
        return new MCPath(this, spotsAtSimStart);
    }
    // otherwise since our simStartDate starts immediately we don't need
    // to use spotsAtSimStart
    return new MCPath(this, valueDate);
}

/** This method is to be retired once we move to 'state variables' for the
    Monte Carlo. This is a hack to get around something that was designed in
    namely the fact that the ref level is responsible for choosing the
    sim start date */
void SampleList::setSimStartDateToFirstRefDate(){
    startSimAtFirstAvgDate = true;
}

/** Returns the total number of past values */
int SampleList::getNbPastValues(int iAsset) const{
    return numDates(iAsset);
}

/** Returns the total number of past values 
    up to and including 'date' */
int SampleList::getNbPastValues(const DateTime& valueDate,
                                int             iAsset) const{
    return valueDate.numPastDates(dates);
}

/** Returns DoubleArray containing historic values for the given asset.
    Array returned is of length equal to the number of past values */
DoubleArray SampleList::getPastValues(const DateTimeArray& userDates,
                                      int                  assetIndex,
                                      const DateTime&      today) const{
    // need to put this in a utility class - similar code in RefLevel.cpp
    const static string routine("SampleList::getPastValues");
    if (assetIndex != 0){ // 'C' style index
        throw ModelException(routine, "No data for asset number "+
                             Format::toString(assetIndex+1));
    }
        
    DoubleArray past;
    int j = 0; // index to dates
    for (int i = 0; 
         i < userDates.size() && today.isGreaterOrEqual(userDates[i]);
         i++){
        // find dates[i] in allDates
        while (j < dates.size() && !dates[j].equals(userDates[i])){
            j++;
        }
        if (j == dates.size()){
            throw ModelException(routine, "No past value for "+
                                 userDates[i].toString());
        }
        double value = values[j];
        if (Maths::isZero(value)){
            throw ModelException(routine, "Past value for asset number "+
                                 Format::toString(assetIndex+1)+ " on "+
                                 userDates[i].toString()+" is 0");
        }
        past.push_back(value);
    }
    return past;
}

/** Records the dates for which this object needs points for 
    within a simulation */
void SampleList::recordDates(SimSeries* simSeries) const{
    simSeries->addDates(dates);
}



// PRIVATE METHODS

SampleList::SampleList(): CObject(TYPE), startSimAtFirstAvgDate(false){}


class SampleListHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(SampleList, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IRefLevel);
        IMPLEMENTS(IPastValues);
        clazz->enableCloneOptimisations();
        EMPTY_SHELL_METHOD(defaultSampleList);
        FIELD(dates, "Sample dates");
        FIELD(values, "Sample values");
        FIELD(weights, "Sample weights");
        FIELD_NO_DESC(startSimAtFirstAvgDate);
        FIELD_MAKE_TRANSIENT(startSimAtFirstAvgDate);
        Addin::registerConstructor("SAMPLE",
                                   Addin::UTILITIES,
                                   "Creates a Sample",
                                   SampleList::TYPE);
        
    }

    static IObject* defaultSampleList(){
        return new SampleList();
    }
};

class SampleListCreateAddin: public CObject{
    static CClassConstSP const TYPE;

    DateTime            valueDate; // to define historical samples
    DateTime            firstDate;
    DateTime            lastDate;
    ExpirySP            matPeriod;
    bool                stubAtEnd;
    string              badDayConv;
    HolidayConstSP      hols;
    double              histAvgToDate; // historic average to date
    double              weightScaling; // weights = weightScaling/ num samples
    
    static IObjectSP create(SampleListCreateAddin* params){
        BadDayConventionSP badDay(
            BadDayConventionFactory::make(params->badDayConv));
        return SampleList::create(params->valueDate, params->firstDate,
                                  params->lastDate,
                                  params->matPeriod, params->stubAtEnd,
                                  badDay,
                                  params->hols, params->histAvgToDate,
                                  params->weightScaling);
    }

    SampleListCreateAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(SampleListCreateAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultSampleListCreateAddin);
        FIELD(valueDate, "value Date");
        FIELD(firstDate, "first date for date generation");
        FIELD(lastDate, "last date for date generation");
        FIELD(matPeriod, "MaturityPeriod eg 1W");
        FIELD(stubAtEnd, "stubAtEnd");
        FIELD(badDayConv, "bad day convention");
        FIELD(hols, "hols");
        FIELD(histAvgToDate, "Historical average so far");
        FIELD(weightScaling, "weightScaling");

        Addin::registerClassObjectMethod("SAMPLE_BETWEEN_DATES",
                                         Addin::UTILITIES,
                                         "Creates a sample with "
                                         "generated dates",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)create);
    }
    static IObject* defaultSampleListCreateAddin(){
        return new SampleListCreateAddin();
    }
    
};

CClassConstSP const SampleListCreateAddin::TYPE = 
CClass::registerClassLoadMethod(
    "SampleListCreateAddin", typeid(SampleListCreateAddin), load);

CClassConstSP const SampleList::TYPE = CClass::registerClassLoadMethod(
    "SampleList", typeid(SampleList), SampleListHelper::load);

DRLIB_END_NAMESPACE
