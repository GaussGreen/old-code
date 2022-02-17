//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CreditSupport.cpp
//
//   Description : Credit support class
//
//   Author      : Ning Shen
//
//   Date        : 10 June 2002
//
//
//----------------------------------------------------------------------------

#ifndef CREDIT_SUPPORT_HPP
#define CREDIT_SUPPORT_HPP


#define MIN_EQUITY_PRICE 1E-6


#include "edginc/Instrument.hpp"
#include "edginc/Model.hpp"
#include "edginc/Asset.hpp"
#include "edginc/FlatVol.hpp"

DRLIB_BEGIN_NAMESPACE

/** credit underlier that spot scan be set */
class CREDIT_DLL CreditUnd
{
public:
    virtual ~CreditUnd() {};

    virtual string getName() const = 0;
    virtual string getTrueName() const = 0;
    virtual CVolProcessed * getProcessedVol(
        const CVolRequest * volRequest ) const = 0;
    virtual void fwdValue(
        const DateTimeArray & dateList,
        CDoubleArray        & result ) const = 0;
    virtual double getSpot() const = 0;

    virtual void setSpot( double spot ) = 0;
};
typedef refCountPtr< CreditUnd > CreditUndSP;

class CREDIT_DLL AssetUnd : public CreditUnd
{
public:
    AssetUnd( IGeneralAsset * asset );
    virtual string getName() const;
    virtual string getTrueName() const;
    virtual CVolProcessed * getProcessedVol(
        const CVolRequest * volRequest ) const;
    virtual void fwdValue(
        const DateTimeArray & dateList,
        CDoubleArray        & result ) const;
    virtual double getSpot() const;
    virtual void setSpot( double spot );
protected:
    IGeneralAsset * m_asset;
};

class CreditSupport;
typedef refCountPtr<CreditSupport> CreditSupportSP;

/** Base credit support class  */
class CreditSupport
{
public:
    class CREDIT_DLL Interface
    {
    public:
        virtual CreditSupportSP createCreditSupport(CMarketDataSP market) = 0;
    };

    /** preprocess instrument for a given set of path dates */
    virtual void preProcess(const DateTimeArray	&dates, 
							const DoubleArray	&atmFwd,
							const DoubleArray	&atmVar) = 0;

    /** calculate values for a given path */
    virtual void calcPathValues(DoubleArray& result, const DateTimeArray& dates, 
                                const double* spots, double spotRef) = 0;

    /** return model for this instrument */
    virtual IModelSP getModel() = 0;

    /** return instrument ccy ISO code */
    virtual string getInstCcyCode() const = 0;

    /** return instrument's last exposure date */
    virtual DateTime getInstLastExposureDate() const = 0;

    ////////// implemented methods //////////////

    /** get equity for spot shift */
    virtual CreditUndSP getUnderlier() const = 0;

    CREDIT_DLL CreditSupport(){};
    CREDIT_DLL virtual ~CreditSupport(){};

    /** shift value date to newDate. perform intermediate theta shift for fixingDates in between */
    /* static */ CREDIT_DLL void shiftValueDate(IObjectSP& inst,  const DateTime& valDate, const DateTime& newDate); 

    // static helper functions 
    /** put fixng dates into groups according to sample date intervals */
    CREDIT_DLL static void groupFixingDates(const DateTimeArray& fixingDatesInput, 
                                 const DateTimeArray& pathDates,
                                 vector<DateTimeArray>& fixDates);


	/** compute price cache */
    CREDIT_DLL /* static */ void computePriceCache
	(
		IObjectSP& inst,
		const DateTimeArray& dates,
		CreditSupport *creditSupport,
		const IModelSP& model,
		DoubleArray &FwdCache					// output
	);
	
	/** compute delta cache */
	/* static */ CREDIT_DLL void computeDeltaCache
	(
		IObjectSP& inst,	
		const DateTimeArray& dates,
		CreditSupport *creditSupport,
		const IModelSP &model,
		double tweakSize,
		const DoubleArray &priceCache,					
		DoubleArray &deltaCache					// output
	);

	/** compute scaling factor for forward starting */
    CREDIT_DLL double computeScaleFactor
	(
		DateTime			startDate,
		const DateTimeArray	&dates,
		const double		*spots,
		DateTime			valueDate,
		double				spotRef
	);

	/** set up spot grid spanning around the equity atm fwd price */
    CREDIT_DLL virtual void setupSpotGrid
	(
		const DateTimeArray	&dates, 
		const DoubleArray	&atmFwd,
		const DoubleArray	&atmVar,
		int					numFactor,
		const double		*factor,
		DoubleArrayArray	&spotGrids	// output
	);

	/** returns linearly interpolated value on otherDate given 
    function goes through (date1, value1) and (date2, value2) */
    CREDIT_DLL static double interpolateDates(
		const DateTime& date1, double value1, 
		const DateTime& date2, double value2,
		const DateTime& otherDate);

	
    /** populates sampleLevels using linear interpolation.  Similar to rollLinear, 
        only applied to flat data, for example as seen in the Equity Swap. */
    CREDIT_DLL static void populateSampleArray(
		const DateTimeArray& sampleDates,
			  DoubleArray&   sampleLevels,
      		  const DateTime& date1, double value1, 
			  const DateTime& date2, double value2);



    /** make a copy of an instrument */
    template<class T> void copyInstrument(smartPtr<T>& copy, T* original)
    {
        IObject* ptr = smartPtr<T>::attachToRef(original)->clone();
        copy =  smartPtr<T> (dynamic_cast<T*>(ptr));
    }

    CREDIT_DLL void setValueDateShift(bool useFwdRate);

private:
    bool useThetaFwdRate;
};

DRLIB_END_NAMESPACE

#endif
