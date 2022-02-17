//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CDFMapping.hpp
//
//   Description : CDFMapping
//
//   Date        : 03 Jun 04
//
//
//   $Log: CDFMapping.hpp,v $
//   Revision 1.8  2004/11/30 18:18:14  nshen
//   use assetMap
//
//   Revision 1.7  2004/11/24 17:01:28  nshen
//   added assets DD for storing market (LN) asset for pdf.
//
//   Revision 1.6  2004/11/23 22:07:25  nshen
//   registered ImpliedSample parameters here. Added model for pdf.
//
//   Revision 1.5  2004/11/11 22:05:53  rgaragnon
//   change of interface due to total re-writing of CDFMapping class
//
//   Revision 1.4  2004/06/24 19:58:45  flu
//   change the differents methods.
//
//   Revision 1.3  2004/06/21 23:04:26  nshen
//   updates.
//
//   Revision 1.2  2004/06/17 18:37:19  flu
//   this class do the mapping spot model to spot market
//
//   Revision 1.1  2004/06/11 22:17:10  nshen
//   modified design for general model.
//
//   Revision 1.1  2004/06/10 20:13:58  flu
//   First version of Mapping to market implied vol surface
//
//----------------------------------------------------------------------------

#ifndef EDR_CDFMapping_HPP
#define EDR_CDFMapping_HPP
#include "edginc/Class.hpp"
#include "edginc/Asset.hpp"
#include "edginc/ImpliedSample.hpp"
#include "edginc/PDFCalculator.hpp"
#include "edginc/PDFDefaultLNStrike.hpp"
#include "edginc/VolProcessedBS.hpp"
#include "edginc/Spline.hpp"


DRLIB_BEGIN_NAMESPACE

/** this class contains MarketCDF and modelCDF objects and perform mapping of
model spot to market spot using respective cdf's */
class MARKET_DLL CDFMapping: public CObject
{
private:
    CDFMapping():CObject(TYPE){}
        
    /** this stores the dates at which the mapping is made */       
    DateTimeArray maturities; // $unregistered
    /** this stores the model spots / probabilities at maturity dates
        and the corresponding market spots at maturity dates */
    DoubleArrayArray spotModel; // $unregistered
    DoubleArrayArray probModel; // $unregistered
    DoubleArrayArray spotMarket; // $unregistered
    IntArray         gridSize; // $unregistered
        
    /** this stores the interpolants to go from model spot to market spot */
    InterpolantArray interpModelMarket; // $unregistered
        
    /** parameters used when creating the market spots grid */
    // registered i/o fields
    IModelConstSP model; // for model CDF computation

    int     numMidStrikesPerLogStrike;
    int     minMidStrikes;
    int     numTailStrikesPerLogStrike;
    int     minTailStrikes;
    double  numMidStdDevs;
    double  numTailStdDevs;
    double  DEBUG_propMidStrikes;
    bool    DEBUG_failProbs;

    // transient but tweakable
    MarketObjectArray assetObjects;
        
public:

    CDFMapping(CClassConstSP clazz): CObject(clazz){}

    // for testing only
    CDFMapping(DateTimeArray maturities,
               const CDFMapping* param): CObject(TYPE),
                                         maturities(maturities)
    {
        // calculate the nber of maturities
        int nberMat = maturities.size();

        // resize the different arrays used in CDFMapping
        spotModel.resize(nberMat);
        probModel.resize(nberMat);
        spotMarket.resize(nberMat);
        gridSize.resize(nberMat);

        interpModelMarket.resize(nberMat);

        numMidStrikesPerLogStrike = param->numMidStrikesPerLogStrike;
        minMidStrikes  = param->minMidStrikes; 
        numTailStrikesPerLogStrike  = param->numTailStrikesPerLogStrike;
        minTailStrikes  = param->minTailStrikes;
        numMidStdDevs  = param->numMidStdDevs; 
        numTailStdDevs  = param->numTailStdDevs;
        DEBUG_propMidStrikes  = param->DEBUG_propMidStrikes;
        DEBUG_failProbs  = param->DEBUG_failProbs;
    };
        
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);
    static IObject* defaultCDFMapping(){
        return new CDFMapping();
    }

    void setAsset(MarketObjectSP asset);
    CAssetConstSP getAsset(const string& name);

    /** set the grid of model spots and the corresponding grid of market spots */
    void setSpotsGrids(const DateTimeArray& maturitiesInp,
                       const CAssetConstSP& asset,
                       const VolRequestLNStrikeSP& volRequest,
                       const DateTime& today,
                       const DateTime& startDate,
                       const PDFDefaultLNStrikeSP& PDFMarket,
                       const PDFCalculatorSP& PDFModel);

    /** build the spline interpolant from model spots to market spots */
    void buildSplineInterpolant();

    /** build the linear interpolant from model spots to market spots */
    void buildLinearInterpolant();

    /** map inSpotPath to outSpotPath for Monte-Carlo engine */
    void impliedMCSpot(double* outSpotPath, const double* inSpotPath, int startStep, int endStep);

    /** map inSpotGrid to outSpotGrid for Finite Differences / Tree engine */
    void impliedSliceSpot(double* outSpotSlice, const double* inSpotSlice, int lo, int hi, int step);

    /** return the market spots grid */
    const DoubleArrayArray getSpotMarket() const;

    /** return the model spots grid */
    const DoubleArrayArray getSpotModel() const;

    /** return model for CDF */
    IModelConstSP getModel(){return model;}
    
    /** override clone() to copy assets which is unregistered*/
    // IObject* CDFMapping:: clone() const;
};

typedef smartPtr<CDFMapping> CDFMappingSP;
typedef smartConstPtr<CDFMapping> CDFMappingConstSP;
typedef array<CDFMappingSP, CDFMapping> CDFMappingArray;

DRLIB_END_NAMESPACE
#endif
