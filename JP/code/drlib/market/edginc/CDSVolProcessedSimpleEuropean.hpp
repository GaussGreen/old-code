//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : CDSVolProcessedSimpleEuropean.hpp
//
//   Description : Processed vol for a simple european CDSOption of the
//                 same type as the underlying BM CDS of the vol cube.
//
//   Author      : Charles Morcom
//
//   Date        : 20 December 2005
//
//
//----------------------------------------------------------------------------

#ifndef CDSVOLPROCESSEDSIMPLEEUROPEAN_HPP
#define CDSVOLPROCESSEDSIMPLEEUROPEAN_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/TimeMetric.hpp"
#include "edginc/VolProcessed.hpp"
#include "edginc/MultiQDistribution.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(CDSVolRequestSimpleEuropean);

/** defines base object used for holding the result of combining instrument
    specific data with the volatility market data */

class MARKET_DLL CDSVolProcessedSimpleEuropean: public CObject,
                                    public virtual IVolProcessed {
public:
    static CClassConstSP const TYPE;
    /** identifies the market data name of the volatility */
    virtual string getName() const;
    /** calculates the trading time between two dates */
    virtual double calcTradingTime(const DateTime &date1, 
                                   const DateTime &date2) const;
    /** retieve time measure for the vol */
    virtual TimeMetricConstSP GetTimeMetric() const;

    virtual ~CDSVolProcessedSimpleEuropean();

//     CDSVolProcessedSimpleEuropean(
// 		const string& name, 
// 		TimeMetricSP metric, 
// 		double px, 
// 		double impVol, 
// 		double atmV, 
// 		double moneyness,
// 		double delta,
// 		double gamma,
// 		double vega,
// 		MultiQDistributionSP mq);
   CDSVolProcessedSimpleEuropean(
		const string& name, 
		TimeMetricSP metric, 
		double px, 
		double impVol, 
		double atmV, 
		double moneyness,
		double delta,
		double gamma,
		double vega,
		MultiQDistributionSP mq,
		double beta = 0.000001);
    /* Note that this is the forward option price at the exercise date!*/
    virtual double optionPrice() const;
    /* Black-Scholes implied vol of the option*/
    virtual double impliedVolatility() const;
    /* ln(F/K)/sigma sqrt(t) where sigma is the ATM vol.*/
    virtual double moneyness() const;
    /* Implied vol for an ATM option */
    virtual double atmVolatility() const;
	virtual double delta() const;
	virtual double gamma() const;
	virtual double vega() const;
  MultiQDistributionSP multiQ() const;
  double meanReversion() const;

protected:
    CDSVolProcessedSimpleEuropean();
private:
    string name;
    TimeMetricSP metric;
    /**Price of the option*/
    double price;
    /**Implied vol of the option*/
    double impliedVol;
    /**ATM implied vol */
    double atmVol;
    /**Moneyness of option*/
    double mny;
	/**Delta of option using BS implied vol*/
	double dlta;
	/**Gamma of option using BS implied vol*/
	double gmma;
	/**Vega of option using BS implied vol */
	double vga;
  /**MultiQDistribution used for pricing*/
  MultiQDistributionSP mq;
  double beta;

    friend class CDSVolProcessedSimpleEuropeanHelper;
};

typedef smartConstPtr<CDSVolProcessedSimpleEuropean> CDSVolProcessedSimpleEuropeanConstSP;
typedef smartPtr<CDSVolProcessedSimpleEuropean> CDSVolProcessedSimpleEuropeanSP;

DRLIB_END_NAMESPACE

#endif
