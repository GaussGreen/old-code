//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SRMEnergyUtil.hpp
//
//   Description : Helper for Energy SRM - used for holding intermediate data
//
//   Author      :
//
//   Date        :
//
//
//----------------------------------------------------------------------------

#ifndef _SRM_ENERGY_UTIL_HPP_
#define _SRM_ENERGY_UTIL_HPP_

#include "edginc/SRMEnergyDiffuse.hpp"
#include "edginc/EnergyFuturesCurve.hpp"
#include "edginc/DECLARE.hpp"
#include <cassert>

DRLIB_BEGIN_NAMESPACE

class SRMEnergyUtilBase:    public virtual VirtualDestructorBase
{
public:
    ///// constructs and populates SRMEnergyUtilBase object
    SRMEnergyUtilBase(
		int	_nbFactors,
		const DateTime & _today,
		EnergyFuturesCurveConstSP _futureCurve,
		const string & _corrInstrStart,
		const string & _corrInstrMaturity
		);

	~SRMEnergyUtilBase() {};

	// Utility type functions:

    // populates sqrtYearFracs with sqrt of year fracs between sim dates.
	// todayIdx is the index of today in the diffusionDates array.
	static void computeSqrtYearFrac(
		int todayIdx,
		const DateTimeArray & diffusionDates,
		vector<double>& sqrtYearFracs
		);

	// Same as above except we assume the first diffusion date is the
	// base date: diffusionDates[0] == baseDate
	static void computeSqrtYearFrac(const DateTimeArray & diffusionDates,
									vector<double>& sqrtYearFracs);

	// compute step-wise k factors (from diffusion time ti to ti+1):
	// used in states diffusion
	static void computeKFactors(vector<double> & _vK,
		double _alpha,
		const vector<double> & _dtSqrts);

	// compute cumulative (from 0 to Ti) k factors:
	// used in computing lnFwd after diffusion is performed
	static void computeCumKFactors(vector<double> & _vKT,
		double _alpha,
		const DateTime & _today,
		const DateTimeArray & _fwdMaturities);

	// compute the gamma's as in Dickinson's paper for the corr instrument
	// (e.g. an energy future with a particular contract maturity)
	virtual vector<double> getVectorGamma() const = 0;

	// compute the matrix A as in Dickinson's paper for the corr instrement
	virtual DoubleMatrix getMatrixAByPCA() const = 0;


    // Accessors:

	/** returns the CDS curve that is diffused */
    EnergyFuturesCurveConstSP getFutureCurve() const { return futureCurve; }

    /** returns today */
    const DateTime& getToday() const { return today; }

	/** returns number of factors */
	int numFactors() const {return nbFactors; }

	/** returns calibrated model vols */
	const DoubleMatrix & getSigmas() const {return sigmas; }

	/** extend calibrated model vols from benchmark dates to timeline.
		note: today is assumed to be in the timeline, and we only
		consider extending the model vols to diffusion dates that
		are strictly > today */
	virtual void setTimeLine(DateTimeArrayConstSP timeline) = 0;

	/** returns calibrated model vols that are extended to the timeline */
	const DoubleMatrix & getExtendedSigmas() const {return extSigmas; }

	/** handy... */
	double getAlpha() {return alpha; }
    double getBeta() {return beta; }

protected:

	/** computes total volatility rates at simulation dates from calibrated
	model vols by assuming option expires = underlying future maturities.

	The volatility rate at simulation dates s_i will be like the option
	implied vol rate, with both the option maturity and underlying future
	maturity equal s_i (equivalent to market quotes).
	**/
	virtual DoubleArray calcTotalVolRatesAtSimDates(
		const DateTimeArray & bmDates, // vol benchmark dates
		const DateTimeArray & timeline, // simulation timeline
        const DateTimeArray & maturities // future maturity associated with each timeline date
		) const = 0;

    const DateTime & today;
	int nbFactors;

	EnergyFuturesCurveConstSP futureCurve;
	EnergyUnderlyerConstSP underlyer;
	EnergyInstVolBaseConstSP instVolBase;

	const string & corrInstrStart; // correlation instrument start label
	const string & corrInstrMaturity; // correlation instrument maturity label
	DoubleMatrix sigmas; // calibrated model vols
	DoubleMatrix extSigmas; // extended sigmas to the full timeline
	double alpha; // all energy model has alpha
    double beta;

}; // end SRMEnergyUtilBase

DECLARE(SRMEnergyUtilBase);


class SRMEnergyUtilOil : public SRMEnergyUtilBase
{
public:
	///// constructs and populates SRMEnergyUtilBase object
	SRMEnergyUtilOil(
		int	_nbFactors,
		const DateTime & _today,
		EnergyFuturesCurveConstSP _futureCurve,
		const string & _corrInstrStart,
		const string & _corrInstrMaturity
		);

	~SRMEnergyUtilOil() {};

	// compute the gamma's as in Dickinson's paper for the corr instrument
	// (e.g. an energy future with a particular contract maturity)
	vector<double> getVectorGamma() const;

	// compute the matrix A (by PCA) as in Dickinson's paper
	DoubleMatrix getMatrixAByPCA() const;

	/** extend calibrated model vols from benchmark dates to timeline.
		note: today is assumed to be in the timeline, and we only
		consider extending the model vols to diffusion dates that
		are strictly > today */
	virtual void setTimeLine(DateTimeArrayConstSP timeline);

	/** handy */
	double getX() const {return x; }
	double getY() const {return y; }
	double getZ() const {return z; }

protected:

	/** computes volatility rates at simulation dates from calibrated
	model vols (piecewise constant with jump times equal benchmark
	dates T_i).

	The volatility rate at simulation dates s_i will be like the
	option implied vol with both the option maturity and underlying
	future maturity equal s_i (equivalent to market quotes).
	**/
	virtual DoubleArray calcTotalVolRatesAtSimDates(
		const DateTimeArray & bmDates, // vol benchmark dates
		const DateTimeArray & timeline, // simulation time line
        const DateTimeArray & maturities // future maturity associated with each timeline date
		) const;

private:

	double x, y, z; // model const calibrated from vol and corr ratios

    // calibration functions:
    // equivalent to EnergyInstVolCalibrated::getNormalizedSigma1
    double calibX(
        const DoubleArray & alphas,
        int firstContract,
        int secondContract,
        double instVolRatio,
        double instCorr) const;

    // equivalent to EnergyInstVolCalibrated::getNormalizedSigma2bar
    double calibZ(
        double x,
        const DoubleArray & alphas,
        int firstContract,
        int secondContract,
        double instVolRatio) const;

    double calibY() const {return 1.0;}

    // equivalent to EnergyInstVolRegular::deriveSigmas
    DoubleMatrix calibSigmas(
        const DateTime & baseDate,
        const DateTimeArray & fFutureMaturityDates,
        const DateTimeArray & fOptionExpiryDates,
        const DoubleArray & vols, // atm vols,
        const DoubleArray & alphas,
        double x,
        double z) const;

    // equivalent to EnergyInstVolBase::getInstVolRatio
    double getInstVolRatio(
        const double alpha,
        const double beta,
        const double v1Bar,
        const double v1,
        const double v2Bar,
        const double v2,
        const double alphaFactor1,
        const double betaFactor1,
        const double alphaFactor2,
        const double betaFactor2) const;

}; // end SRMEnergyUtilOil


class SRMEnergyUtilGas : public SRMEnergyUtilBase
{
public:
    ///// constructs and populates SRMEnergyUtilBase object
    SRMEnergyUtilGas(
        int	_nbFactors,
        const DateTime & _today,
        EnergyFuturesCurveConstSP _futureCurve,
        const string & _corrInstrStart,
        const string & _corrInstrMaturity
        );

    ~SRMEnergyUtilGas() {};

    // compute the gamma's as in Dickinson's paper for the corr instrument
    // (e.g. an energy future with a particular contract maturity)
    vector<double> getVectorGamma() const;

    // compute the matrix A (by PCA) as in Dickinson's paper
    DoubleMatrix getMatrixAByPCA() const;

    /** extend calibrated model vols from benchmark dates to timeline.
    note: today is assumed to be in the timeline, and we only
    consider extending the model vols to diffusion dates that
    are strictly > today */
    virtual void setTimeLine(DateTimeArrayConstSP timeline) {};

    /** handy */
    const DoubleArray & getX() const {return x; }
    const DoubleArray & getZ() const {return z; }

protected:

    /** computes volatility rates at simulation dates from calibrated
    model vols (piecewise constant with jump times equal benchmark
    dates T_i).

    The volatility rate at simulation dates s_i will be like the
    option implied vol with both the option maturity and underlying
    future maturity equal s_i (equivalent to market quotes).
    **/
    virtual DoubleArray calcTotalVolRatesAtSimDates(
        const DateTimeArray & bmDates, // vol benchmark dates
        const DateTimeArray & timeline, // simulation time line
        const DateTimeArray & maturities // future maturity associated with each timeline date
        ) const {};

private:

    DoubleArray x;  // calibrated seasonality consts
    DoubleArray z;  // calibrated seasonality consts

    // calibration functions:
    // equivalent to EnergyInstVolCalibrated::getNormalizedSigma1
    double calibX() const;

    // equivalent to EnergyInstVolCalibrated::getNormalizedSigma2bar
    double calibZ() const;

    // equivalent to EnergyInstVolRegular::deriveSigmas
    DoubleMatrix calibSigmas() const;

    // equivalent to EnergyInstVolSeasonal::getInstVolRatio
    double getInstVolRatio(const double alpha,
        const double beta,
        const double normV1_1,
        const double normV1_2,
        const double normV2Bar,
        const double normV2,
        const double alphaFactor1,
        const double betaFactor1,
        const double alphaFactor2,
        const double betaFactor2) const;

}; // end SRMEnergyUtilGas


/** Create the appropriate tier 1 energy util class */
SRMEnergyUtilBaseSP SRMEnergyUtilCreate(
	string				        _modelType,
	int							_nbFactors,
	const DateTime &			_today,
	EnergyFuturesCurveConstSP	_futureCurve,
	const string &				_corrInstrStart,
	const string &				_corrInstrMaturity
	);


/***************************************************************************/
/******************** Tier 2 Energy Util Stuffs ****************************/
/***************************************************************************/
class SRMEnergyUtilTier2 : public SRMEnergyUtilBase
{
public:
	///// constructs and populates SRMEnergyUtilTier2 object
	SRMEnergyUtilTier2(
		int	_nbFactors,
		const DateTime & _today,
		EnergyFuturesCurveConstSP _futureCurve,
		const string & _corrInstrStart,
		const string & _corrInstrMaturity
		);

	virtual ~SRMEnergyUtilTier2() {};

	// compute the gamma's as in Dickinson's paper for the corr instrument
	// (e.g. an energy future with a particular contract maturity)
    vector<double> getVectorGamma() const;

	// compute the matrix A (by PCA) as in Dickinson's paper
    DoubleMatrix getMatrixAByPCA() const;

	/** not implemented for tier 2 energy util */
    virtual void setTimeLine(DateTimeArrayConstSP timeline);

protected:

	/** not implemented for tier 2 energy util **/
	virtual DoubleArray calcTotalVolRatesAtSimDates(
		const DateTimeArray & bmDates, // vol benchmark dates
		const DateTimeArray & timeline, // simulation time line
        const DateTimeArray & maturities // future maturity associated with each timeline date
		) const { return DoubleArray(0); }

private:

}; // end SRMEnergyUtilTier2


DRLIB_END_NAMESPACE
#endif
