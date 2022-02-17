//----------------------------------------------------------------------------
//
//   Group       : Credit Derivatives Research
//
//   Filename    : SCIDCalibparameters.cpp
//
//   Date        : August, 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SCIDCalibparameters.hpp"
#include "edginc/Format.hpp"

DRLIB_BEGIN_NAMESPACE

/* Make sure the class links */
bool  SCIDCalibparametersLoad() {
    return (SCIDCalibParameters::TYPE != 0);
   }

SCIDCalibParameters::SCIDCalibParameters(const CClassConstSP& clazz) 
    : MarketObject(clazz)
{
}

SCIDCalibParameters::SCIDCalibParameters()
	: MarketObject(SCIDCalibParameters::TYPE),
	name("sCID calibration parameter default Name")
{
}

SCIDCalibParameters::~SCIDCalibParameters()
{
}

string SCIDCalibParameters::getName() const
{
    return name;
}

/** populate from market cache - default implementation provided */
void SCIDCalibParameters::getMarket(const IModel* model, const MarketData* market) {
    static const string method = "SCIDCalibParameters::getMarket";
	try
	{
        market->GetReferenceDate(today);
        ELCalibrateMaturities.resize(ELCalibrateExpiries->size());
		for (int k=0; k<ELCalibrateMaturities.size(); k++)
			ELCalibrateMaturities[k] = (*ELCalibrateExpiries)[k]->toDate(today);
                
        // populates trancheQuotes
		for (int i = 0; i < trancheQuotes->size(); ++i) {
            (*(*trancheQuotes)[i]).getData(model, market);
		}
        
        for (int i=0; i< trancheQuotes->size(); ++i)
            QLIB_VERIFY((*(*trancheQuotes)[i])->getMaturityDates()->size() == (*trancheWeights)[i].size(), 
                method + " CDOTrancheQuotes and trancheWeights do not have the same number of maturities for " + (*(*trancheQuotes)[i])->getName());                   
        
        for (int i=0; i< trancheQuotes->size(); ++i)
            QLIB_VERIFY((*(*trancheQuotes)[i])->getMaturityDates()->size() == (*trancheStartDates)[i].size(), 
                method + " CDOTrancheQuotes and trancheStartDates do not have the same number of maturities for " + (*(*trancheQuotes)[i])->getName());                   
        
        // builds strikesToTrancheQuotesIdxMap
        initStrikesToTrancheQuotesIdxMap();
    }
    catch (exception& e)
	{
		throw ModelException(e, method);
	}
}

/** Init strikesToTrancheQuotesIdxMap */
void SCIDCalibParameters::initStrikesToTrancheQuotesIdxMap() {
    // clear map
    strikesToTrancheQuotesIdxMap.clear();

    // build strikesToTrancheQuotesIdxMap
    for (int i = 0; i < trancheQuotes->size(); ++i) {
        strikesToTrancheQuotesIdxMap[pair<double, double>(
            (*(*trancheQuotes)[i])->getLowStrike(),
            (*(*trancheQuotes)[i])->getHighStrike())] = i;
    }
}

void SCIDCalibParameters::getTrancheAttachementPoints(DoubleArray& _kmin, DoubleArray& _kmax) {
    _kmin.resize(strikesToTrancheQuotesIdxMap.size());
    _kmax.resize(strikesToTrancheQuotesIdxMap.size());
    int size = trancheQuotes->size();
    for (int i=0;  i < size; ++i) {
        _kmin[i] = (*(*trancheQuotes)[i])->getLowStrike();
        _kmax[i] = (*(*trancheQuotes)[i])->getHighStrike();
    }
}

DateTimeArraySP SCIDCalibParameters::getTrancheMaturities(double _kmin, double _kmax) {
    // Retrieve tranche quotes index
    map<pair<double, double>, int>::const_iterator iter =
        strikesToTrancheQuotesIdxMap.find(pair<double, double>(_kmin,_kmax));
    if (iter == strikesToTrancheQuotesIdxMap.end()) 
        throw ModelException(
            "SCIDCalibParameters::getTrancheMaturities", "No maturities available for tranche " +
            Format::toString(_kmin) + "-" + Format::toString(_kmax));
    int trancheQuoteIdx = iter->second;
    
    return (*(*trancheQuotes)[trancheQuoteIdx])->getMaturityDates();
}

DoubleArraySP SCIDCalibParameters::getTrancheWeights(double _kmin, double _kmax) {
    // Retrieve tranche quotes index
    map<pair<double, double>, int>::const_iterator iter =
        strikesToTrancheQuotesIdxMap.find(pair<double, double>(_kmin,_kmax));
    if (iter == strikesToTrancheQuotesIdxMap.end()) 
        throw ModelException(
            "SCIDCalibParameters::getTrancheMaturities", "No maturities available for tranche " +
            Format::toString(_kmin) + "-" + Format::toString(_kmax));
    int trancheQuoteIdx = iter->second;
    
    return DoubleArraySP(new DoubleArray((*trancheWeights)[trancheQuoteIdx]));
}

DateTimeArray SCIDCalibParameters::getTrancheStartDates(double _kmin, double _kmax) {
    // Retrieve tranche quotes index
    map<pair<double, double>, int>::const_iterator iter =
        strikesToTrancheQuotesIdxMap.find(pair<double, double>(_kmin,_kmax));
    if (iter == strikesToTrancheQuotesIdxMap.end()) 
        throw ModelException(
            "SCIDCalibParameters::getTrancheStartDates", "No maturities available for tranche " +
            Format::toString(_kmin) + "-" + Format::toString(_kmax));
    int trancheQuoteIdx = iter->second;
    
    return (*trancheStartDates)[trancheQuoteIdx];
}

DoubleArraySP SCIDCalibParameters::getTranchePrices(double _kmin, double _kmax) {
    // Retrieve tranche quotes index
    map<pair<double, double>, int>::const_iterator iter =
        strikesToTrancheQuotesIdxMap.find(pair<double, double>(_kmin,_kmax));
    if (iter == strikesToTrancheQuotesIdxMap.end()) 
        throw ModelException(
            "SCIDCalibParameters::getTranchePrices", "No weights available for tranche " +
            Format::toString(_kmin) + "-" + Format::toString(_kmax));
    int trancheQuoteIdx = iter->second;
    
    DoubleArraySP result(new DoubleArray());
    DateTimeArraySP maturities = (*(*trancheQuotes)[trancheQuoteIdx])->getMaturityDates();
    for(DateTimeArray::iterator pos = maturities->begin(); pos != maturities->end(); ++pos)
        result->push_back((*(*trancheQuotes)[trancheQuoteIdx])->getSpread(*pos));
    return result;        
}

DoubleArraySP SCIDCalibParameters::getTrancheUpFrontPremiums(double _kmin, double _kmax) {
    // Retrieve tranche quotes index
    map<pair<double, double>, int>::const_iterator iter =
        strikesToTrancheQuotesIdxMap.find(pair<double, double>(_kmin,_kmax));
    if (iter == strikesToTrancheQuotesIdxMap.end()) 
        throw ModelException(
            "SCIDCalibParameters::getTranchePrices", "No weights available for tranche " +
            Format::toString(_kmin) + "-" + Format::toString(_kmax));
    int trancheQuoteIdx = iter->second;
    
    DoubleArraySP result(new DoubleArray());
    DateTimeArraySP maturities = (*(*trancheQuotes)[trancheQuoteIdx])->getMaturityDates();
    for(DateTimeArray::iterator pos = maturities->begin(); pos != maturities->end(); ++pos)
        result->push_back((*(*trancheQuotes)[trancheQuoteIdx])->getUpfront(*pos));
    return result;        
}

/**Check stuff here */
void SCIDCalibParameters::validatePop2Object() 
{
    static const string routine("SCIDCalibParameters::validatePop2Object");

	try 
    {
        //verify that CDS quotes and weights have the same size
        QLIB_VERIFY(trancheQuotes->size() == trancheWeights->size(), routine+ " CDOTrancheQuotes and trancheWeights do not have the same size");
	}
	catch (exception& e)
	{
		throw ModelException(e, routine);
	}

}

/** Specific clone method to copy "strikesToTrancheQuotesIdxMap" field */
IObject* SCIDCalibParameters::clone() const {
    SCIDCalibParameters* copy = DYNAMIC_CAST(SCIDCalibParameters, MarketObject::clone());
    copy->strikesToTrancheQuotesIdxMap = strikesToTrancheQuotesIdxMap;
    return copy;
}

void SCIDCalibParameters::load(CClassSP& clazz)
{
    clazz->setPublic(); 
    REGISTER(SCIDCalibParameters, clazz);
    SUPERCLASS(MarketObject);
	EMPTY_SHELL_METHOD(defaultSCIDCalibParameters);

	FIELD(name, "Name of the scid calibration parameter class");
	FIELD(TrancheCoupon,      "Tranche Coupon");
	FIELD( dcc,      				 "Day Count Convention");
	FIELD(today,  "Today");
	FIELD_MAKE_TRANSIENT(today);
    FIELD(ELCalibrateExpiries, "expiries to calibrate EL");
    FIELD(ELCalibrateMaturities, "maturities to calibrate EL");
    FIELD_MAKE_TRANSIENT(ELCalibrateMaturities);
	FIELD(ELCalibrateWeights, "weights to calibrate EL");
    FIELD(trancheQuotes, "Tranche quotes");
    FIELD(trancheWeights, "Tranche weights");
    FIELD(trancheStartDates, "Tranche start dates");
}

IObject* SCIDCalibParameters::defaultSCIDCalibParameters(){
    return new SCIDCalibParameters();
}

CClassConstSP const SCIDCalibParameters::TYPE = CClass::registerClassLoadMethod(
   "SCIDCalibParameters", typeid(SCIDCalibParameters), load);

DEFINE_TEMPLATE_TYPE(SCIDCalibParametersWrapper);
DRLIB_END_NAMESPACE

