//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SRMSwaption.hpp
//
//   Description : Helper for SRM - used for calibration against a swaption
//
//   Author      : Mark A Robson
//
//   Date        : 11 June 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_SRMSWAPTION_HPP
#define EDR_SRMSWAPTION_HPP
#include "edginc/SRMSwaptionPricer.hpp"
#include "edginc/SRMRatesHJMUtil.hpp"
#include "edginc/RootFinder.hpp"
//#include "edginc/SVGenIRSwap.hpp"


DRLIB_BEGIN_NAMESPACE

/** Helper for SRM - for calibrating against swaption vols */
class MCARLO_DLL SRMSwaption :  public virtual VirtualDestructorBase
{

#define SIGMADELTA    0.0001 


    //Used in SRMSwaption::calibVol2Q
    class SRMSwaptionFunctor : public Func1D::NoDeriv
    {
        const SRMSwaption &     m_swaption;
        const SRMRatesHJMUtil & m_ratesUtil;
        const DateTimeArray &   m_mergedList;
        vector<double>      &   m_extendedSpotVol;
        double                  m_targetSwaptionVol;
    public:
        SRMSwaptionFunctor( const SRMSwaption &swaption, 
            const SRMRatesHJMUtil &ratesUtil,
            const DateTimeArray &mergedList,
            vector<double>      &extendedSpotVol,
            double targetVol) : m_swaption(swaption),
            m_ratesUtil(ratesUtil),
            m_mergedList(mergedList),
            m_extendedSpotVol(extendedSpotVol),
            m_targetSwaptionVol(targetVol)
        {

        }

        virtual double operator()(double x) const
        {
            double impVol = m_swaption.impVol2Q(m_ratesUtil, m_mergedList,
                x, m_extendedSpotVol);
            double diff = impVol - m_targetSwaptionVol;
            return diff;        
        }

        double approxDerivative(double x, double f, double delta)
        {
            double sigma_2 = x + delta;
            double impVol = f + m_targetSwaptionVol;

            double impVol_2 = m_swaption.impVol2Q(m_ratesUtil, m_mergedList,
                sigma_2, m_extendedSpotVol);

            /* approximate derivative of the implied vol with respect to sigma */
            double derivative = (impVol_2 - impVol) / delta;
            return derivative;
        }
    };

    
public:

    /** Constructor using explicit data. Use today for prevSwaptionExpiry
        for first swaption */
    SRMSwaption(const SRMRatesHJMUtil&  ratesHJMUtil,
                const DateTime&         swaptionExpiry,
                const DateTime&         swapStart,
                const DateTime&         swapMat,
                const DateTime&         prevSwaptionExpiry,
                const DateTimeArray&    dates);

    /** Implied vol of a call(?) swaption. Returns < 0 if implied vol cannot
        be found */
    double impliedVol(double swaptionPrice) const;

    /* Calibrates the Swaption to the target implied vol target
     * From swapvol::CalibSwaption */
    void calibrate(const SRMRatesHJMUtil& ratesHJMUtil,
                   double target,
                   vector<double>& SpotVol /* (M) */ ) const;

    /* Main calibration routine. From SRM3::CalibVol2Q */
    static DateTimeArray calibVol2Q(
						 const SRMRatesHJMUtil& ratesHJMUtil,
 						 bool skipFlag,
						 DoubleArray& SpotVol); // O

    /* Inner calibration step: access is needed by functor used by Newton-Raphson **/
    double impVol2Q( // from SRM3::SwpnVol2Q
        const SRMRatesHJMUtil&     ratesHJMUtil,
        const DateTimeArray&  mergedList, // zero dates + vol dates
        double                sigma,           /* spotvol prevailing between
                                               expiry(i-1) and expiry(i) */
                                               vector<double>&       extendedSpotVol) const; /* : dicretized spotVol */

    /* MAR: This is the old calibration method - used currently for quanto
    adjustment to CDS Par Spreads */
    static DoubleArray spotVolVF(
					  const SRMRatesHJMUtil&  ratesHJMUtil,
					  bool  skipFlag);

    /** Construct a Pricer which can be used to price paths for this swaption */
    SRMSwaptionPricer* createPricer(const SRMRatesHJMUtil& ratesHJMUtil) const;

private:
    SRMSwaption();
    SRMSwaption(const SRMSwaption& rhs);
    SRMSwaption& operator=(const SRMSwaption& rhs);
    void muFactor(const SRMRatesHJMUtil& ratesHJMUtil,
                  const DateTime& currDate,
                  const DateTime& swapStart,
                  const DateTime& swapMat,
                  double&         muFacSquare1, /* (O) */
                  double&         muFacSquare2 /* (O) */) const;

    DateTime findTau(const SRMRatesHJMUtil&  ratesHJMUtil,
                     const DateTime& swapStart,
                     const DateTime& swapMat,
                     const DateTime& today) const;

    void swapDetails(const SRMRatesHJMUtil& ratesHJMUtil,
                     const DateTime& swapStart,
                     const DateTime& swapMat,
                     bool            stubAtEnd,
                     bool            useDiffusedCurve);

    static double meanFSquare(
        double t,    /* (I) time t : cumulative variance v(t) (see report)*/
        double qR,   /* (I) QRight                                        */
        double qL);  /* (I) QLeft                                         */

    void populate(const SRMRatesHJMUtil& ratesHJMUtil,
                  const DateTime& swaptionExpiry,
                  const DateTime& swapStart,
                  const DateTime& swapMat,
                  int             prevIdx,
                  const DateTimeArray& dates);

    static DateTimeArray extendTau(
        const SRMRatesHJMUtil& ratesHJMUtil,
        const DateTimeArray& tau); /* (I) : tau value at each vol dates  */

    void correlMuF(
        const SRMRatesHJMUtil&     ratesHJMUtil,
        const DateTimeArray&  mergedList, // zero dates + vol dates
        const DateTime&       swapStart,
        const DateTime&       swapMat,
        const vector<double>& extendedSpotVol, /* : descretized spotVol */
              vector<double>& correl)const;//(O)correlation between mu^2 and F^2

    DateTime       swapStart;      // swap start date
    DateTime       swapMat;        // swap maturity date
    DateTime       swapExp;        // swaption expiry date
    int            idx;            /*index of the swaption we are calibrating */
    double         expiry;         /*Swaption expiry */
    /* Discretization variables */
    vector<double> variance;      /*Cumulative variance */
    vector<double> rbart;         /* array of rbar values at dicrete times */
    vector<double> del_t;         /*time step array for the discretization */
    vector<double> muFacSquare1;  /* mu factor squared for lognormal assumption
                                     as specified in the doc*/
    vector<double> muFacSquare2;  /* mu factor squared for normal assumption
                                     as specified in the doc*/
    DateTime       tau;           /* tau as specified in the doc for a given
                                     swaption                  */
    vector<vector<double> > mr;   /* integration of the mean-reversion term
                                     mr[date][factor] */
    int            NumCpns;       /* number of coupon payment dates */
    DateTimeArray  CpnPayDates;   /* coupon payment dates */
    vector<double> weights;       /* deterministic frozen weights */
    int            prevIdx;       /* called current_index in SRM3 code
                                     'idx' of the last swaption calibrated */
    /* Swap details */
    double         ParYield;      /* Swap par Yield */
    double         Annuity;       /* Swap annuity */
};

typedef smartPtr<SRMSwaption> SRMSwaptionSP;

DRLIB_END_NAMESPACE
#endif
