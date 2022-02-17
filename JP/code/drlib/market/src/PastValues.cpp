//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PastValues.cpp
//
//   Description : Handles idea of a reference level for a series of assets
//
//   Author      : Mark A Robson
//
//   Date        : 19 Oct 2001
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_PASTVALUES_CPP
#include "edginc/PastValues.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Asset.hpp"
#include "edginc/SimSeries.hpp"
#include "edginc/PastObservations.hpp"
#include "edginc/RevertTypeConvert.hpp"

DRLIB_BEGIN_NAMESPACE

/** Returns observation types if they exist 
default implementation returns dummies*/
ObservationTypeArray* IPastValues::getObsTypes(IMultiFactors* assets) {
    int nbAssets = getNumAssets();
    ObservationTypeArray* obsTypes = new ObservationTypeArray(nbAssets);
    for (int i = 0; i < nbAssets; i++) {
        (*obsTypes)[i] = ObservationType::make("NotUsed");
    }
    return obsTypes;
}

/** Returns observation types if they exist 
default implementation returns dummies*/
ObservationSourceArray* IPastValues::getSources(IMultiFactors* assets) {
    int nbAssets = getNumAssets();
    ObservationSourceArray* sources = new ObservationSourceArray(nbAssets);
    for (int i = 0; i < nbAssets; i++) {
        (*sources)[i] = ObservationSourceSP(IMarketObservable::getDefaultObsSource());
    }
    return sources;
}
/** Indicates whether this PastValues object uses sampling rules to construct a history.
 *  By default, this indicates that they are not constructed this way. */
bool IPastValues::hasSamplingConvention() {
    return false;
}

IPastValues::MissingSampleException::~MissingSampleException() throw (){}

/** Instances of this class are thrown when a zero historic value 
    is tried to be used */
IPastValues::MissingSampleException::MissingSampleException(
    const string&    routine,
    int              assetIndex,
    const DateTime&  sampleDate):
    ModelException(routine, "Past value for asset number "+
                   Format::toString(assetIndex+1)+ " on "+
                   sampleDate.toString()+" is 0.0"),
    assetIndex(assetIndex), sampleDate(sampleDate){}


/** returns the index of the asset which had the missing sample */
int IPastValues::MissingSampleException::getAssetIndex() const{
    return assetIndex;
}

/** records the asset name whose sample was missing */
void IPastValues::MissingSampleException::setAssetName(const string& name){
    setDescription("Past value for asset "+name+
                   " on "+ sampleDate.toString()+" is missing");
}

IPastValues::MissingSampleException::MissingSampleException(
    const MissingSampleException& e): 
    ModelException(create(e)), assetIndex(e.assetIndex), sampleDate(e.sampleDate) {}

IPastValues::MissingSampleException::MissingSampleException(
    const ModelException& e,
    int                   assetIndex,
    const DateTime&       sampleDate):
    ModelException(create(e)), assetIndex(assetIndex), sampleDate(sampleDate){}

/** creates a [deep] copy of the exception */
ModelException* IPastValues::MissingSampleException::clone() const{
    MissingSampleException* e = new MissingSampleException(*this,
                                                           assetIndex,
                                                           sampleDate);
    return e;
}

/** indicates whether this exception is derived from ModelException -
    used to drive whether this is stored as a 'cause' when a new
    exception is created */
bool IPastValues::MissingSampleException::isDerived() const{
    return true;
}

/** Returns the 'cause' exception if it's a MissingSampleException
    otherwise returns null. Equivalent to
    dynamic_cast<MissingSampleException*>e.getCause() */
IPastValues::MissingSampleException* 
IPastValues::MissingSampleException::getInstance(
    ModelException& e){
    return dynamic_cast<MissingSampleException*>(e.getCause());
}

bool IPastValues::PositivePastValueCheck::isValid(double value){
    return Maths::isPositive(value);
}

bool IPastValues::NonNegativePastValueCheck::isValid(double value){
    return !Maths::isNegative(value);
}

/** Class implementing IPastValues where each asset has the same set of
    "averaging in" dates */
template <class ValidPastValueCheck>
class IPastValues_Simple: public CObject,
                          virtual public IPastValues, 
						  virtual public IRevertTypeConvert{
private:
    DateTimeArray allDates;
    DoubleMatrix  pastValues; // dereferened as [assetIndex][dateIndex]
    int           numAssets;  /* optional parameter - useful at the spreadsheet
                                 when there is no historic data */
    bool          throwMissingFutureDatesException;
public:
    static CClassConstSP const TYPE;

    /** Creates a IPastValues where each asset has the same set of simulation
        dates. */
    IPastValues_Simple(const DateTimeArray& allDates,
                       const CDoubleMatrix& pastValues): 
    CObject(TYPE), 
    allDates(allDates), 
    pastValues(pastValues), 
    numAssets(0),
    throwMissingFutureDatesException(false) {
        validatePop2Object();
    }

    /** Returns the total number of past values */
    virtual int getNbPastValues(int iAsset) const{
        return allDates.size();
    }

    /** Returns the total number of past values 
        up to and including 'date' */
    virtual int getNbPastValues(const DateTime& valueDate,
                                int             iAsset) const{
        return valueDate.numPastDates(allDates);
    }

    /** Returns the number of assets for which this object has data for */
    virtual int getNumAssets() const{
        return numAssets; // pastValues.numCols();
    }

    /** Returns an array holding the values on the requested dates. The
        length of the array is equal to the number of historic dates.
        An exception is thrown if there are no values for any of the dates
        or any historic values are 0 or if iAsset is out of bounds */
    virtual DoubleArray getPastValues(const DateTimeArray& dates,
                                      int                  iAsset,
                                      const DateTime&      today) const{
        const static string routine("IPastValues::Simple::getPastValues");
        if (iAsset < 0 || iAsset >= numAssets){
            throw ModelException(routine, "No data for asset number "+
                                 Format::toString(iAsset+1)+
                                 " (index out of bounds)");
        }

        DoubleArray past;
        int j = 0; // index to allDates
        for (int i = 0; 
             i < dates.size() && today.isGreaterOrEqual(dates[i]);
             i++){
            // find dates[i] in allDates
            while (j < allDates.size() && !allDates[j].equals(dates[i])){
                j++;
            }
            if (j == allDates.size()){
                throw MissingSampleException(routine, iAsset, dates[i]);
            }
            double value = pastValues[iAsset][j];
            if (!ValidPastValueCheck::isValid(value)){
                throw MissingSampleException(routine, iAsset, dates[i]);
            }
            past.push_back(value);
        }
        return past;
    }

    /** rolls past values for all assets inside a MultiFactor */
    virtual void roll(const Theta::Util&         thetaUtil,
                      const IMultiMarketFactors* multiFactor){
        static const string method("PastValues::Util::roll");
        try{
            const DateTime& origDate = thetaUtil.getOriginalValueDate();
            const DateTime& newDate = thetaUtil.getNewValueDate();
            bool pastRollDate = false;
            for (int i = 0; !pastRollDate && i < allDates.size(); i++){
                for (int iAsset = 0; iAsset < getNumAssets(); iAsset++){
                    IMarketFactorConstSP factor(multiFactor->getFactor(iAsset));
                    if (!IGeneralAsset::TYPE->isInstance(factor)){
                        throw ModelException(method, "Past values only support "
                                             "for IGeneralAssets");
                    }
                    const IGeneralAsset* asset = DYNAMIC_CAST(IGeneralAsset,
                                                              factor.get());
                    if ((allDates[i].isGreater(origDate) && 
                         !allDates[i].isGreater(newDate)) ||
                        (allDates[i].equals(newDate)    
                         && !ValidPastValueCheck::isValid(pastValues[iAsset][i]))) {
                        double newSpot = thetaUtil.useAssetFwds()?
                            asset->fwdValue(allDates[i]): asset->getSpot();
                        pastValues[iAsset][i] = newSpot;
                    } else if (allDates[i].isGreater(newDate) ) {
                        pastRollDate = true;
                    }
                }
            }
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

    /** Rolls the PastValues through time populating any samples as they
        become historic */
    virtual void roll(const Theta::Util&  thetaUtil,
                      int                 iAsset,
                      const Asset*        asset){
        try{
            const DateTime& origDate = thetaUtil.getOriginalValueDate();
            const DateTime& newDate = thetaUtil.getNewValueDate();
            bool pastRollDate = false;
            for (int i = 0; !pastRollDate && i < allDates.size(); i++){
                if ((allDates[i].isGreater(origDate) && 
                     !allDates[i].isGreater(newDate)) ||
                    (allDates[i].equals(newDate)    
                     && !ValidPastValueCheck::isValid(pastValues[iAsset][i]))){                    
                    pastValues[iAsset][i] = 
                        asset->getThetaSpotOnDate(thetaUtil.getTheta(),
                                                  allDates[i]);
                } else if (allDates[i].isGreater(newDate) ) {
                    pastRollDate = true;
                }
            }
        } catch (exception& e){
            throw ModelException(e, "PastValues::Util::roll");
        }
    }                            

    // validates and fills in past samples if possible
    // ensures all dates are there and all past dates have valid samples
    // For past dates the provided samples are used and if missing
    // samples are retrieved from the centralised asset history
    virtual void validate(const DateTime& today,
                          const DateTimeArrayArray& requiredDates,
                          IMultiFactors* assets,
                          SamplingConvention* dateAdjust) {
/*        static const string routine("IPastValues::Simple::validate");

        int j = 0;
        vector<DateTime>::iterator iter(allDates.begin());
        for(int i=0; i<requiredDates.size(); i++) {
           bool found = false;
            for(; !found && iter != allDates.end()
                            && (*iter)<=requiredDates[i]; ++iter, ++j) {
                if ((*iter)==requiredDates[i]) {
                    found = true;
                }
            }
            if (!found) {
                // need to add the missing date in the list
                // but note cannot get sample as Simple PastValues has
                // no ability to do centralised sampling
                // if it's in the past it will get picked up later as 
                // a dodgy past sample
                iter = allDates.insert(iter, requiredDates[i]);
                pastValues.rowAdd(j-1, 0.0);
            }
        }*/
    }

    virtual void validatePop2Object(){
        static const string routine("IPastValues::Simple:validatePop2Object");
        DateTime::ensureIncreasing(allDates, "Sample dates", false);
        // sort out the numAssets field
        int numCols = pastValues.numCols();
        if (numCols > 0){
            // have a proper matrix
            if (numAssets == 0){
                // numAssets not specified
                numAssets = numCols;
            } else if (numAssets != numCols){
                throw ModelException(routine, "Inconsistent numAssets with "
                                     "num of columns in matrix");
            }
        } else if (numAssets > 0){
            // resize matrix to have correct number of cols
            pastValues = DoubleMatrix(numAssets, 0);
        }
        if (pastValues.numRows() < allDates.size()){
            // be generous and fill in with zeros - this happens on the
            // s/sheet as the interface ignores blank cells
            DoubleMatrix origVals(pastValues);
            pastValues = DoubleMatrix(numAssets, allDates.size());
            for (int i = 0; i < numAssets; i++){
                for (int j = 0; j < origVals.numRows(); j++){
                    pastValues[i][j] = origVals[i][j];
                }
            }
        }   
        if (allDates.size() != pastValues.numRows()){
            throw ModelException(routine,
                                 "Mismatch between number of dates (" + 
                                 Format::toString(allDates.size()) +
                                 ") and number of past values (" + 
                                 Format::toString(pastValues.numRows()) +
                                 ")");
        } 
    }

	// for the IMS interface (see IRevertTypeConvert)
	virtual IObjectSP revert(const string& interfaceType) const {
		static const string method = "IPastValues::Simple::revert";
		try {
			if (interfaceType != IRevertTypeConvert::PYRAMID) {
				throw ModelException( method, 
					"Cannot convert a IPastValues::Simple for"
					" the interface " + interfaceType);               
			}

			CashFlowClusterSP cashFlows(new CashFlowCluster(numAssets));
			int numRows = pastValues.numRows();
			DoubleArraySP dbl(new DoubleArray(numRows));//createCashFlowArray copies the values anyways
			for (int i=0;i<numAssets;++i) {
				for (int j=0;j<numRows;++j) {
					(*dbl)[j]=pastValues[i][j];
				}
				cashFlows->set(i,CashFlow::createCashFlowArray(allDates,*dbl));
			}
			return cashFlows;
		}
		catch (exception& e) {
			throw ModelException(e, method);
		}
	}

    // Retrieves all past samples used by an instrument
    // For past dates the provided samples are used and if missing
    // samples are retrieved from the centralised asset history
    virtual void pastSamplesEvents(const DateTime& today,
                          const DateTimeArrayArray& requiredDates,
                          IMultiFactors* assets,
                          SamplingConvention* dateAdjust,
                          EventResults* events) {
    }

private:
    /* for reflection */
    IPastValues_Simple(): 
    CObject(TYPE),
    numAssets(0),
    throwMissingFutureDatesException(false){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(IPastValues_Simple<ValidPastValueCheck>, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IPastValues);
		IMPLEMENTS(IRevertTypeConvert);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(allDates, "Sample dates for all assets");
        FIELD(pastValues, "Historic values for each asset");
        FIELD(numAssets, "Number of assets");
        FIELD_MAKE_OPTIONAL(numAssets);
        FIELD_NO_DESC(throwMissingFutureDatesException);  // DO NOT USE
        FIELD_MAKE_OPTIONAL(throwMissingFutureDatesException);
    }

    static IObject* defaultCtor(){
        return new IPastValues_Simple<ValidPastValueCheck>();
    }
};

template<> CClassConstSP const IPastValues_Simple<IPastValues::PositivePastValueCheck>::TYPE =
CClass::registerClassLoadMethod(
    "IPastValues::Simple",  // for backward compatibility
    typeid(IPastValues_Simple<PositivePastValueCheck>), load);

template<> CClassConstSP const IPastValues_Simple<IPastValues::NonNegativePastValueCheck>::TYPE =
CClass::registerClassLoadMethod(
    "IPastValues::Simple<NonNegativePastValueCheck>", 
    typeid(IPastValues_Simple<NonNegativePastValueCheck>), load);

/** Creates a PastValues where each asset has the same set of simulation
    dates. */
IPastValues* IPastValues::Util::makeSimple(
    const DateTimeArray&  commonSimDates,
    const CDoubleMatrix&  pastValues){ // numCols = numAssets
    return new IPastValues_Simple<PositivePastValueCheck>(commonSimDates, pastValues);
}

/** Creates a PastValues where each asset has the same set of sample
    dates using a NonNegativePastValue policy */
IPastValues* IPastValues::Util::makeSimple(
    const DateTimeArray&       commonSimDates,
    const CDoubleMatrix&       pastValues,
    NonNegativePastValueCheck  notUsed){
    return new IPastValues_Simple<NonNegativePastValueCheck>(commonSimDates, pastValues);
}

/** Creates a PastValues where each asset 'averages' in on a single
    day */
IPastValues* IPastValues::Util::makeTrivial(
    const DateTime&      singleDate,
    const DoubleArray&   pastValues){ // one entry per asset
    DateTimeArray commonSimDates(1, singleDate);
    DoubleMatrix  matrixPastValues(pastValues);
    matrixPastValues.transpose();
    return makeSimple(commonSimDates, matrixPastValues);
}

        
/** Creates a SimSeries for a single asset which 'averages' in
    on a single day */
IPastValues* IPastValues::Util::makeTrivial(
    const DateTime&      singleDate,
    const double         pastValue){ // one entry per asset
    DoubleArray pastValues(1, pastValue);
    return makeTrivial(singleDate, pastValues);
}

/** Class implementing IPastValues where each asset has the same set of
    "averaging in" dates */
class IPastValues::General: public CObject,
                            virtual public IPastValues{
private:
    CashFlowCluster valuesByAsset; // array of cash flow arrays
    // source and obs type for looking up each asset's past samples
    // form the centralised asset history
    StringArray             sources; 
    ObservationTypeArray    obsTypes;
    StringArray             fxSources; 
    ObservationTypeArray    fxObsTypes;
    bool            hasSamplingInfo;
    bool            throwMissingFutureDatesException; // NOT USED
    // very lame and temporary handling of overrides for event purposes
    IntArrayArray           overrideFlags;
public:
    static CClassConstSP const TYPE;

    /** Creates a IPastValues where each asset has a different set of
        simulation dates. */
    General(const CashFlowCluster& samples):
        CObject(TYPE), valuesByAsset(samples),
        hasSamplingInfo(false),
        throwMissingFutureDatesException(false) { 
        validatePop2Object();
    }

    /** Creates a IPastValues where each asset has a different set of
        simulation dates and we pass in source/obs type for each asset*/
    General(const PastObservationsArray& samples): CObject(TYPE),
                            hasSamplingInfo(true),
                            throwMissingFutureDatesException(false) {
        int numAssets = samples.size();
        valuesByAsset = CashFlowCluster(numAssets);
        sources = StringArray(numAssets);
        obsTypes = ObservationTypeArray(numAssets);
        fxSources = StringArray(numAssets);
        fxObsTypes = ObservationTypeArray(numAssets);
       
        for (int i = 0; i < numAssets; i++) {
            valuesByAsset[i] = samples[i]->samples;
            sources[i] = samples[i]->source;
            obsTypes[i] = samples[i]->obsType;
            fxSources[i] = samples[i]->fxSource;
            fxObsTypes[i] = samples[i]->fxObsType;
        }
        validatePop2Object();
    }

    /** Returns the total number of past values */
    virtual int getNbPastValues(int iAsset) const{
        return valuesByAsset[iAsset].size();
    }

    /** Returns the total number of past values 
        up to and including 'date' */
    virtual int getNbPastValues(const DateTime& valueDate,
                                int             iAsset) const{
        // copied from DateTime::numPastDates
        int numDates = valuesByAsset[iAsset].size();
        int i;
        for (i = 0; i < numDates && valueDate.isGreaterOrEqual(valuesByAsset[iAsset][i].date); i++); // empty
        return i;
    }

    /** Returns the number of assets for which this object has data for */
    virtual int getNumAssets() const{
        return valuesByAsset.size();
    }

    /** Returns an array holding the values on the requested dates. The
        length of the array is equal to the number of historic dates.
        An exception is thrown if there are no values for any of the dates
        or any historic values are 0 or if iAsset is out of bounds */
    virtual DoubleArray getPastValues(const DateTimeArray& dates,
                                      int                  iAsset,
                                      const DateTime&      today) const{
        const static string routine("IPastValues::General::getPastValues");
        if (iAsset < 0 || iAsset >= valuesByAsset.size()){
            throw ModelException(routine, "No data for asset number "+
                                 Format::toString(iAsset+1)+
                                 " (index out of bounds) ");
        }

        DoubleArray past;
        const CashFlowArray& cfArray = valuesByAsset[iAsset];
        int j = 0; // index to cfArray
        for (int i = 0; 
             i < dates.size() && today.isGreaterOrEqual(dates[i]);
             i++){
            // find dates[i] in cfArray
            while (j < cfArray.size() && !cfArray[j].date.equals(dates[i])){
                j++;
            }
            if (j == cfArray.size()){
                throw MissingSampleException(routine, iAsset, dates[i]);
            }
            double value = cfArray[j].amount;
            if (Maths::isZero(value)){
                throw MissingSampleException(routine, iAsset, dates[i]);
            }
            past.push_back(value);
        }
        return past;
    }

    /** rolls past values for all assets inside a MultiFactor */
    virtual void roll(const Theta::Util&         thetaUtil,
                      const IMultiMarketFactors* multiFactor){
        static const string method("PastValues::Util::roll");
        try{
            const DateTime& origDate = thetaUtil.getOriginalValueDate();
            const DateTime& newDate = thetaUtil.getNewValueDate();
            for (int iAsset = 0; iAsset < getNumAssets(); iAsset++){
                IMarketFactorConstSP factor(multiFactor->getFactor(iAsset));
                if (!IGeneralAsset::TYPE->isInstance(factor)){
                    throw ModelException(method, "Past values only support "
                                         "for IGeneralAssets");
                }
                const IGeneralAsset* asset = DYNAMIC_CAST(IGeneralAsset,
                                                          factor.get());
                bool pastRollDate = false;
                CashFlowArray& cfArray = valuesByAsset[iAsset];
                for (int i = 0; !pastRollDate && i < cfArray.size(); i++){
                    const DateTime& currDate = cfArray[i].date;
                    if ((currDate.isGreater(origDate) && 
                         !currDate.isGreater(newDate)) ||
                        (currDate.equals(newDate)    
                         && Maths::isZero(cfArray[i].amount)))
                    {
                        double newSpot = thetaUtil.useAssetFwds()?
                            asset->fwdValue(currDate): asset->getSpot();
                        cfArray[i].amount = newSpot;
                    } else if (currDate.isGreater(newDate) ) {
                        pastRollDate = true;
                    }
                }
            }
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

    /** Rolls the PastValues through time populating any samples as they
        become historic */
    virtual void roll(const Theta::Util&  thetaUtil,
                      int                 iAsset,
                      const Asset*        asset){
        try{
            const DateTime& origDate = thetaUtil.getOriginalValueDate();
            const DateTime& newDate = thetaUtil.getNewValueDate();
            bool pastRollDate = false;
            CashFlowArray& cfArray = valuesByAsset[iAsset];
            for (int i = 0; !pastRollDate && i < cfArray.size(); i++){
                const DateTime& currDate = cfArray[i].date;
                if ((currDate.isGreater(origDate) && 
                     !currDate.isGreater(newDate)) ||
                    (currDate.equals(newDate) && 
                     Maths::isZero(cfArray[i].amount))){
                    cfArray[i].amount = 
                        asset->getThetaSpotOnDate(thetaUtil.getTheta(), 
                                                  currDate);
                } else if (currDate.isGreater(newDate) ) {
                    pastRollDate = true;
                }
            }
        } catch (exception& e){
            throw ModelException(e, "PastValues::Util::roll");
        }
    }                            
      
    /** Returns observation types if they exist*/
    ObservationTypeArray* getObsTypes(IMultiFactors* assets) {
        int nbAssets = getNumAssets();
        ObservationTypeArray* obs = new ObservationTypeArray(nbAssets);
        if (hasSamplingInfo) {
            for (int i = 0; i < nbAssets; ++i) {
                if (assets->assetGetCcyTreatment(i)== CAsset::CCY_TREATMENT_STRUCK){ 
                    (*obs)[i] =
                        ObservationTypeSP(new StruckObservationType(obsTypes[i].get(),
                                                                    fxObsTypes[i].get()));
                } else {
                    (*obs)[i] = obsTypes[i];
                }
            }
        } else {
            for (int i = 0; i < nbAssets; i++) {
                (*obs)[i] = ObservationType::make("NotUsed");
            }
        }
        return obs;
    }

    /** Returns observation types if they exist*/
    ObservationSourceArray* getSources(IMultiFactors* assets) {
        int nbAssets = getNumAssets();
        ObservationSourceArray* obsSources = new ObservationSourceArray(nbAssets);
        if (hasSamplingInfo) {
            for (int i = 0; i < nbAssets; ++i) {
                if (assets->assetGetCcyTreatment(i)== CAsset::CCY_TREATMENT_STRUCK){ 
                    (*obsSources)[i] = 
                        ObservationSourceSP(new StruckObservationSource(sources[i],
                                                                        fxSources[i]));
                } else {
                    (*obsSources)[i] = 
                        ObservationSourceSP(new ObservationSource(sources[i]));
                }
            }
        } else {
            for (int i = 0; i < nbAssets; i++) {
                (*obsSources)[i] = ObservationSourceSP(IMarketObservable::getDefaultObsSource());
            }
        }
        return obsSources;
    }

    /** Indicates whether this PastValues object uses sampling rules to construct a history. */
    bool hasSamplingConvention() {
        // If we have samplingInfo, then the history is constructed using its rules.
        return hasSamplingInfo;
    }

    // validates and fills in past samples if possible
    // ensures all dates are there and all past dates have valid samples
    // For past dates the provided samples are used and if missing/zero
    // samples are retrieved from the centralised asset history
    // note samples are for SAMPLING dates not OBSERVATION dates so if we have 
    // a sample for a sampling date it will be used even if the one for
    // the obs date is different
    virtual void validate(const DateTime& today,
                          const DateTimeArrayArray& requiredDates,
                          IMultiFactors* assets,
                          SamplingConvention* dateAdjust) {
        static const string routine("IPastValues::General::validate");
        try {
            int numAssets = valuesByAsset.size();
            int iAsset;

            // if we're going to do any centralised sampling we first need
            // to set up some objects for struck assets
            ObservationSourceArraySP srcs(getSources(assets));
            ObservationTypeArraySP obs(getObsTypes(assets));

            // now fill in all the missing dates
            for(iAsset=0; iAsset<numAssets; iAsset++) {
                CashFlowArray& cfArray = valuesByAsset[iAsset];
                IntArraySP overrideArray = overrideFlags[iAsset]; 

                vector<CashFlow>::iterator iter(cfArray.begin());
                vector<int>::iterator overrideIter(overrideArray->begin());
                for(int i=0; i<requiredDates[iAsset].size(); i++) {
                    bool found = false;
                    for( ; !found && iter != cfArray.end()
                                    && (*iter).date <= requiredDates[iAsset][i]; ) {
                        if ((*iter).date==requiredDates[iAsset][i]) {
                            found = true;
                        } else {
                            ++iter;
                            ++overrideIter;
                        }
                    }
                    if (!found || (requiredDates[iAsset][i] < today  &&
                                    !Maths::isPositive((*iter).amount))) {
                        double sample = 0.0;
                        if (requiredDates[iAsset][i] < today && hasSamplingInfo) {
                            try {
                                sample = assets->getAsset(iAsset).pastValue(requiredDates[iAsset][i],
                                                                    (*obs)[iAsset].get(),
                                                                    (*srcs)[iAsset].get(),
                                                                    0, // fixing type
                                                                    0, // overrides
                                                                    dateAdjust);
                            } catch (exception& e) {
                                throw (ModelException(e, "Past value for asset " +
                                            assets->getAsset(iAsset).getTrueName() +
                                            " for date " + 
                                            requiredDates[iAsset][i].toString() +
                                            " is missing or is <= 0.0 and "
                                            "cannot be retrived from the "
                                            "centralised asset history"));
                            }
                        }
                        if (requiredDates[iAsset][i] < today  && !Maths::isPositive(sample)) {
    #if 1
    // recent change to MissingSampleException may fix a vanishing call stack so leave
    // this in while we see. Otherwise the 'else' clause is the fix. 
                            MissingSampleException mse(routine, iAsset, requiredDates[iAsset][i]);
                            mse.setAssetName(assets->getName(iAsset));
                            throw mse;
    #else
                            throw ModelException(routine, "Past value for asset "+
                                                assets->assetGetTrueName(iAsset)+
                                                " on "+ requiredDates[iAsset][i].toString()+" is missing");
    #endif
                        }
                        if (!found) {
                            // need to add the missing date in the list
                            iter = cfArray.insert(iter,  CashFlow(requiredDates[iAsset][i], sample));
                            overrideIter = overrideArray->insert(overrideIter, 0);
                        } else {
                            // just need to overwrite the zero sample
                            (*iter).amount = sample;
                            (*overrideIter) = 0; // it's ot an override now
                        }
                    }
                }
            }
        } catch (exception& e) {
            throw ModelException(e, routine);
        }   
    }

    virtual void validatePop2Object(){
        static const string routine("IPastValues::General:validatePop2Object");
        int i;
        int numAssets = valuesByAsset.size();
        for (i = 0; i < numAssets; i++){
            CashFlow::ensureDatesIncreasing(valuesByAsset[i], 
                                            "Sample dates for asset "+
                                            Format::toString(i+1), 
                                            false);
        }
        // check to see if we have sampling info
        // and check the 4 arrays are the same length 
        // no real way of missing this from currrent interfaces but let's do it
        // as this is how xml files sent by Analytics will appear to us
        // i.e. a IPastValues::General objects rather than arrays of
        // CashFlowArrays or PastObservations
        if (sources.size() > 0) {
            hasSamplingInfo = true;
        }
        if (sources.size() != obsTypes.size()) {
            throw ModelException("IPastValues::General::validatePop2Object",
                "Must be as many sources as there are observation types");
        }
        if (sources.size() != fxSources.size()) {
            throw ModelException("IPastValues::General::validatePop2Object",
                "Must be as many sources as there are FX sources");
        }
        if (sources.size() != fxObsTypes.size()) {
            throw ModelException("IPastValues::General::validatePop2Object",
                "Must be as many sources as there are FX observation types");
        }
        // now set up the rather lame arrays of flags telling you which
        // levels are overrides (because we fill in the gaps with non-overrides)
        // at this stage everything present is an override!!
        overrideFlags = IntArrayArray(valuesByAsset.size());
        for (int i = 0; i < numAssets; ++i) {
            overrideFlags[i] = IntArraySP(new IntArray(valuesByAsset[i].size(), 1));
        }
    }

    // Retrieves all past samples used by an instrument
    // For past dates the provided samples are used and if missing
    // samples are retrieved from the centralised asset history
    virtual void pastSamplesEvents(const DateTime& today,
                          const DateTimeArrayArray& requiredDates,
                          IMultiFactors* assets,
                          SamplingConvention* dateAdjust,
                          EventResults* events) {
        static const string routine("IPastValues::General:validatePop2Object");
        try {
            int numAssets = valuesByAsset.size();
            int iAsset;

            // if we're going to do any centralised sampling we first need
            // to set up some objects for struck assets
            ObservationSourceArraySP srcs(getSources(assets));
            ObservationTypeArraySP obs(getObsTypes(assets));
            PastSamplesCollectorSP collector =
                        PastSamplesCollectorSP(new PastSamplesCollector(today));

            for(iAsset=0; iAsset<numAssets; iAsset++) {
                CashFlowArray& cfArray = valuesByAsset[iAsset];
                const CAsset* thisAsset = &(assets->getAsset(iAsset));
                AssetFixTypeSP fixType = 
                        AssetFixTypeSP(new AssetFixType(thisAsset->getName()));
                DateTimeArray assetDates = requiredDates[iAsset];

                vector<CashFlow>::iterator iter(cfArray.begin());
                vector<int>::iterator overrideIter(overrideFlags[iAsset]->begin());
                for(int i=0; i<assetDates.size() &&
                            assetDates[i] <= today; i++) {
                    // find the sample - we know it's there
                    // ISS actually if PastValues ONLY had the ones we wanted
                    // we could just step through the list
                    bool found = false;
                    for( ; !found && iter != cfArray.end()
                                    && (*iter).date <= assetDates[i]; ) {
                        if ((*iter).date==assetDates[i]) {
                            found = true;
                        } else {
                            ++iter;
                            ++overrideIter;
                        }
                    }
                    if ((*overrideIter) == 1) {
                        // add overriden sample 
                        collector->addSample(assetDates[i],
                                             assetDates[i], // sample = obs
                                             (*srcs)[iAsset].get(),
                                             (*obs)[iAsset].get(),
                                             (*iter).amount,
                                             fixType.get(),
                                             true); // override
                    } else {
                        // we've done validation so we must be able to do this
                        thisAsset->addPastSampleEvent(assetDates[i],
                                                (*obs)[iAsset].get(),
                                                (*srcs)[iAsset].get(),
                                                fixType.get(),
                                                0, // overrides use eventually
                                                dateAdjust,
                                                collector.get());
                    }
                }
            }
            // now put all the events in the EventsResults object
            collector->finaliseResults(events);
        } catch (exception& e) {
            throw ModelException(e, routine);
        }   
    }

    static IObjectSP fromCashFlowArrayArray(const IObjectSP& object, 
                                            CClassConstSP    requiredType) {
        CashFlowArrayArray& samples =
            dynamic_cast<CashFlowArrayArray&>(*object);
        return IObjectSP(new General(samples));
    }

    // ability to construct a General PastValues from an array of
    // PastObservations objects (each contains cashflow array and source/obs type)
    static IObjectSP fromPastObsArray(const IObjectSP& object, 
                                      CClassConstSP    requiredType) {
        PastObservationsArray& samples =
            dynamic_cast<PastObservationsArray&>(*object);
        return IObjectSP(new General(samples));
    }

private:

    /* for reflection */
    General(): CObject(TYPE), hasSamplingInfo(false),
               throwMissingFutureDatesException(false){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(General, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IPastValues);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(valuesByAsset, "Averaging in dates and historic values"
                     " per asset");
        FIELD_NO_DESC(sources);
        FIELD_MAKE_OPTIONAL(sources);
        FIELD_NO_DESC(obsTypes);
        FIELD_MAKE_OPTIONAL(obsTypes);
        FIELD_NO_DESC(fxSources);
        FIELD_MAKE_OPTIONAL(fxSources);
        FIELD_NO_DESC(fxObsTypes);
        FIELD_MAKE_OPTIONAL(fxObsTypes);
        FIELD_NO_DESC(hasSamplingInfo);
        FIELD_MAKE_TRANSIENT(hasSamplingInfo);
        FIELD_NO_DESC(overrideFlags);
        FIELD_MAKE_TRANSIENT(overrideFlags);
        FIELD(throwMissingFutureDatesException, "DO NOT USE ME"); // DO NOT USE
        FIELD_MAKE_OPTIONAL(throwMissingFutureDatesException);
        registerObjectFromArrayMethod(CashFlowArrayArray::TYPE,
                                      TYPE,
                                      &fromCashFlowArrayArray);
        registerObjectFromArrayMethod(PastObservationsArray::TYPE,
                                      TYPE,
                                      &fromPastObsArray);
    }

    static IObject* defaultCtor(){
        return new General();
    }
};

/** Creates a PastValues where each asset can have a different
    set of sample dates. */
IPastValues* IPastValues::Util::makeGeneral(
    const CashFlowCluster&  samplesPerAsset){
    return new General(samplesPerAsset);
}

CClassConstSP const IPastValues::General::TYPE =
CClass::registerClassLoadMethod(
    "IPastValues::General", typeid(IPastValues::General), load);

    

void IPastValues::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER_INTERFACE(IPastValues, clazz);
    EXTENDS(IObject);
}    

CClassConstSP const IPastValues::TYPE = CClass::registerInterfaceLoadMethod(
    "IPastValues", typeid(IPastValues), load);


DRLIB_END_NAMESPACE
