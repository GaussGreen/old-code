//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : SCIDCalibParameters.hpp
//
//   Description : Parameters needed for the calibration of sCID model 
//
//   Date        : July 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_sCID_CALIB_PARAMETERS_HPP
#define QLIB_sCID_CALIB_PARAMETERS_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/CDOTrancheQuotes.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE_WRAPPER(SCIDCalibParameters)
FORWARD_DECLARE(MarketData)
FORWARD_DECLARE(Expiry)

typedef smartPtr<SCIDCalibParameters> SCIDCalibParametersSP;
typedef smartConstPtr<SCIDCalibParameters> SCIDCalibParametersConstSP;

class MARKET_DLL SCIDCalibParameters: public MarketObject
{
public:

    static CClassConstSP const TYPE;
    virtual ~SCIDCalibParameters();

    virtual string getName() const;

	  /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);
    virtual void validatePop2Object();

	int getTrancheCoupon() { return TrancheCoupon; }
	DayCountConventionSP getTrancheDCC() {  return dcc;}
    DoubleArrayConstSP getPELCalibrationWeights(){return ELCalibrateWeights;}
    DateTimeArray getPELCalibrationMaturities(){return ELCalibrateMaturities;}
   
    /** methods which use CDOTrancheQuotes **/
    void getTrancheAttachementPoints(DoubleArray& _kmin, DoubleArray& _kmax);
	DateTimeArraySP getTrancheMaturities(double _kmin, double _kmax);
	DateTimeArray getTrancheStartDates(double _kmin, double _kmax);
	DoubleArraySP getTrancheWeights(double _kmin, double _kmax);
	DoubleArraySP getTranchePrices(double _kmin, double _kmax);
	DoubleArraySP getTrancheUpFrontPremiums(double _kmin, double _kmax);
    
    
    /** Specific clone method to copy "strikesToTrancheQuotesIdxMap" field */
    virtual IObject* clone() const;
protected:  
    SCIDCalibParameters(const CClassConstSP& clazz);
	SCIDCalibParameters(const SCIDCalibParameters& v);
	SCIDCalibParameters();

    static void load(CClassSP& clazz);
	static IObject* defaultSCIDCalibParameters();

    /** the name of SCID calibration parameters */
    string               name;
private: 

	/** frequency of the coupon, in Months */
	int TrancheCoupon;
	/** Day count convention used for the tranche legs computations*/
	DayCountConventionSP			dcc;

	DateTime today;
    /** Expiries to calibrate EL **/
    ExpiryArrayConstSP  ELCalibrateExpiries;
    
    /** Maturities to calibrate EL **/
    DateTimeArray  ELCalibrateMaturities;
    
    /** Weights to calibrate EL **/
	DoubleArrayConstSP ELCalibrateWeights;

	/** Tranche quotes and weights for calibration[Mandatory] */
	CDOTrancheQuotesWrapperArraySP trancheQuotes;
   
    /** Tranche weights **/
    DoubleArrayArraySP  trancheWeights;
    //** Tranche Start Dates **/
    DateTimeClusterSP   trancheStartDates;
    /**
     * Internal map between (low strike, high strike) pair and
     * tranche quotes and weights index [non exposed]
     * */
    map<pair<double, double>, int> strikesToTrancheQuotesIdxMap;     // $unregistered
    
    /** Init strikesToTrancheQuotesIdxMap */
    void initStrikesToTrancheQuotesIdxMap();
};


DRLIB_END_NAMESPACE

#endif
