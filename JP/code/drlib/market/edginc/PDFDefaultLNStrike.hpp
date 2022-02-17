//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PDFDefaultLNStrike.hpp
//
//   Description : Implementation of PDFCalculator for parameterised vols
//
//   Author      : Mark A Robson
//
//   Date        : 24 March 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDR_PDF_DEFAULT_LNSTRIKE_HPP
#define EDR_PDF_DEFAULT_LNSTRIKE_HPP

#include "edginc/PDFCalculator.hpp"
#include "edginc/VolBase.hpp"
#include "edginc/Asset.hpp"
#include "edginc/PDFRequestLNStrike.hpp"
#include "edginc/VolParam.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL PDFDefaultLNStrike : public PDFCalculator {
public:
    static CClassConstSP const TYPE;

    virtual ~PDFDefaultLNStrike();

    /** calculates transient fields */
    virtual void validatePop2Object();

    /** awaiting spec */
    virtual void probabilities(const DoubleArray& strikes,
                               const DateTime&    maturity,
                               DoubleArray&       probs) const;

    /** awaiting spec */
    virtual void probabilities(const CLatticeDouble&    strikes,
                               const DateTimeArray&     maturities,
                               CLatticeDouble&          probs) const;
    
    /** awaiting spec */
    virtual void localDensity(const DoubleArray& strikes,
                              const DateTime&    maturity,
                              DoubleArray&       density) const;

    /** awaiting spec */
    virtual void integratedDensity(const DoubleArray& strikes,
                                   const DateTime&    maturity,
                                   DoubleArray&       density) const;

    /** constructor - to do alter type of request to PDFRequestLNStrike */
    PDFDefaultLNStrike(const DateTime&   valueDate,
                       const CAsset*     asset,
                       const PDFRequest* request);  

protected:
    /** constructor for derived classes */
    PDFDefaultLNStrike(const CClassConstSP&             clazz,
                       const DateTime&                  valueDate,
                       const CAssetConstSP&             asset,
                       const PDFRequestLNStrikeConstSP& request);  

    /* same as next spreads method but takes in strikes as a DoubleArray.
       Just calls next spreads() method having mapped data types */
    virtual void spreads(const DoubleArray& strikes,
                         const DateTime&    maturity,
                         double&            adjEpsilon,
                         DoubleArray&       probs) const;
    
    /** Caclulates spreads using specified strike(1+epsilon) and
        strike(1-epsilon) via call to spreads routine below */
    virtual void spreads(const CLatticeDouble&    strikes,
                         const DateTimeArray&     maturities,
                         double&                  adjEpsilon,
                         CLatticeDouble&          spread) const;

    /** core routine that actually does the work. Derived classes will
        probably want to override this method */
    virtual void spreads(const CLatticeDouble&    loStrikes,
                         const CLatticeDouble&    hiStrikes,
                         const DateTimeArray&     maturities,
                         const DoubleArray&       fwds, // at maturities
                         CLatticeDouble&          spread) const;

    /** Populates supplied lostrikes and histrikes with strikes * (1.0
        - epsilon) and strikes * (1.0 + epsilon) where epsilon is
        modified so that none of the strikes overlap */
    static void constructStrikes(
        const CSliceDouble&     strikes,
        double&                 epsilon, // note modified
        CSliceDouble&           lostrikes,
        CSliceDouble&           histrikes);
    
    /** validates that the strikes for specified maturity are ok */
    static void validateStrikes(
        const CSliceDouble& strikes,
        const DateTime&     maturity);

    PDFDefaultLNStrike(const CClassConstSP& clazz);

    // fields
    DateTime                    valueDate;
    CAssetConstSP               asset;
    PDFRequestLNStrikeConstSP   lnRequest;
    double                      epsilon;
    double                      accuracy;
    // following field are transient
    // DateTime                    startDate;
    mutable bool                fwdStarting;
    bool                        volFwdStarting;
    mutable double              spotAtStart;
private:
    PDFDefaultLNStrike();
    class Helper;
    friend class Helper;

};

typedef smartConstPtr<PDFDefaultLNStrike> PDFDefaultLNStrikeConstSP;
typedef smartPtr<PDFDefaultLNStrike> PDFDefaultLNStrikeSP;

DRLIB_END_NAMESPACE
#endif
