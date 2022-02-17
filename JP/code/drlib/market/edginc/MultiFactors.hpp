//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MultiFactors.hpp
//
//   Description : Interface that exposes some of the internals of a collection of assets
//
//   Date        : Oct 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_OPEN_ASSETS_HPP
#define EDG_OPEN_ASSETS_HPP
#include "edginc/MultiMarketFactors.hpp"
#include "edginc/VolBase.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/OutputName.hpp"
#include "edginc/Correlation.hpp"
#include "edginc/CorrelationTerm.hpp"
#include "edginc/LocalCorrSqueeze.hpp"
#include "edginc/TimeMetric.hpp"
#include "edginc/ObservationType.hpp"
#include "edginc/ObservationSource.hpp"
#include "edginc/ObservationOverride.hpp"
#include "edginc/SamplingConvention.hpp"

DRLIB_BEGIN_NAMESPACE
class SensitiveStrikeDescriptor;
class PDFCalculator;
class PDFRequest;
class Results;

/* Provides a view of a collection of assets as a cooperating set :
   while there may be M "natural" assets we allow that they may be
   viewed as a collection of N "component" assets, or factors.  We are
   opening some of the internals to public view.  Example being n
   assets viewed as M=1 (an XCB) or N=n correlated assets (the XCB
   components).  So for pricing we'd express any payoff as a function
   of the M assets, but for modelling we would simulate N.  */
class MARKET_DLL IMultiFactors : virtual public IMultiMarketFactors {
public:
    static CClassConstSP const TYPE;

    // "Asset" view 

    //// Returns the number of assets
    virtual int NbAssets() const = 0;

    //// Returns the name of the selected asset
    virtual string assetGetName(int iAsset) const = 0;

    //// Returns the pure equity name of the selected asset (e.g. not "_P_" synthetic)
    virtual string assetGetTrueName(int iAsset) const = 0;

    //// Returns the ccy treatment for selected asset
    virtual string assetGetCcyTreatment(int iAsset) const = 0;

    //// Returns the spot value of the selected asset
    virtual double assetGetSpot(int iAsset) const = 0;

    //// calculates the forward value of the selected asset on the supplied
    //// date
    virtual double assetFwdValue(int             iAsset,
                                 const DateTime& date) const = 0;

    //// calculates the forward value of the selected asset on the supplied
    //// dates
    virtual void assetFwdValue(int                  iAsset,
                               const DateTimeArray& dates,
                               CDoubleArray&        result) const = 0;

    
    // "Factor" view
    //// Returns the number of factors
    virtual int NbFactors() const = 0;

    //// Returns the name of the asset selected by its factor index
    virtual string factorGetName(int iFactor) const = 0;

    /** Returns an processed vol - which combines the vol market data with the
        instrument data in the volRequest for given factor */
    virtual IVolProcessed * factorGetProcessedVol(
        int                iFactor,
        const CVolRequest* volRequest) const = 0;

    //// Returns the spot value of the asset selected by its factor index
    virtual double factorGetSpot(int iFactor) const = 0;

    /** Calculates the expected spot price of the given factor at the
        given date */
    virtual double factorFwdValue(int             iFactor,
                                  const DateTime& date) const = 0;

    /** Calculates the expected spot price of the given factor at each of the
        given dates */
    virtual void factorFwdValues(int                  iFactor,
                                 const DateTimeArray& dates,
                                 CDoubleArray&        result) const = 0;

    // Is there a correlation class to generalise+hide the matrix?
    // This is NbFactors x NbFactors 
    virtual CDoubleMatrixConstSP factorsCorrelationMatrix() const = 0;

    /** Pull out the data which is necessary for CorrTS */
    virtual BoolArray               getSkipFwdCorrArray() const = 0;
    virtual CorrelationCommonArray        getCorrObjArray() const = 0;
    virtual CorrelationTermArray    getCorrTermObjArray() const = 0;
    virtual LocalCorrSqueezeArray     getLocalCorrSqueezeArray() const = 0;
    virtual IntArray                  getLocalCorrSqueezeIndexArray() const = 0;
    virtual TimeMetricArray         getTimeMetricArray() const = 0;
        
    /** A utility to record forwards at maturity (for each factor) */
    virtual void recordFwdAtMat(OutputRequest*  request,
                                Results*        results,
                                const DateTime& maturityDate) const = 0;

    /* "Low level" method to translate from "factors" to "assets".
       e.g. N paths of factors need to be translated into M for use in payoff
       Uses preallocated memory.
     */
    virtual void CollapseFactors(
        int             len,                     // length of each array
        const double**  factorLevels,            // [NbFactors] arrays
        double**        assetLevels) const = 0;  // [NbAssets] arrays

    /** Queries selected asset for its sensitive strikes given the
        supplied volRequest. The sensitive strikes are added to the
        sensitiveStrikes parameter */
    virtual void assetGetSensitiveStrikes(
        int                               iAsset,
        const CVolRequest*                volRequest,
        const OutputNameConstSP&          outputName,
        const SensitiveStrikeDescriptor&  sensStrikeDesc,
        const DoubleArraySP&              sensitiveStrikes) const = 0;

    /** Returns a list of asset indexes indicating which asset is
        sensitive to supplied SensControl. If crossAssetSensitivities
        is true then this only covers sensitivities such as phi which
        can live at this level otherwise these are excluded. The
        relevant name of what is being tweaked must be stored in the
        SensControl */
    virtual IntArray getSensitiveAssets(
        const Sensitivity* sens,
        bool               crossAssetSensitivities) const = 0;

    /** Returns a reference to the i-th Asset (which of course potentially
        breaks the whole asset/factor view if you use this in conjunction
        with the factorXXXX methods eg in a path generator) */
    virtual const CAsset& getAsset(int iAsset) const = 0;
    
    /** returns PDF calculator for iAsset */
    virtual PDFCalculator* assetPdfCalculator(
        const PDFRequest* request,
        int               iAsset) const = 0;

    /** Helper for centralised sampling/isda adjustment in the n-factor case
    gets a list of dates for which the given sample will finally be known for 
    each asset. Results are put into obsDates array and the BoolArray tracks
    those assets for which date is omitted*/
    virtual void observationDates(const DateTime&                sampleDate,
                                  const ObservationSourceArray&  sources,
                                  const SamplingConvention*      sampleRule,
                                  DateTimeArray&                 obsDates,
                                  BoolArray&                     sampled) const;

    /** Helper for parametrised (SPI) dates in the n-factor case
    return a bool indicating whether any of the assets or their components
    have a sampling holiday on the given date */
    virtual bool isHoliday(const DateTime&                sampleDate,
                           const ObservationSourceArray&  sources) const;

protected:
    IMultiFactors() {};
};

typedef smartConstPtr<IMultiFactors> IMultiFactorsConstSP;
typedef smartPtr<IMultiFactors> IMultiFactorsSP;

DRLIB_END_NAMESPACE
#endif
