//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FD2DEQCR.hpp
//
//   Description : two factor finite difference engine for credit
//                    asset factor    dX = drift*dt + sigma*dW
//                    credit factor        dp = k(p_0 - p)dt + vol*sqrt(p)*dZ
//                    Z and W have a correlation rho.
//
//   Author      : Bruno Melka
//
//   Date        : 10 Jun 2005
//
//----------------------------------------------------------------------------

#ifndef FD2DEQCR_HPP
#define FD2DEQCR_HPP

#include "edginc/FD2D.hpp"
#include "edginc/VolProcessedBS.hpp"
#include "edginc/MarketDataFetcherLN.hpp"
#include "edginc/CDSParSpreads.hpp"
#include "edginc/ICDSVol.hpp"

DRLIB_BEGIN_NAMESPACE

class TREE_DLL FD2DEQCR : public FD2D{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& );
    static IObject* defaultFD2DEQCR(){
        return new FD2DEQCR();
    }
    
    friend class FD2DSolver;

    FD2DEQCR();
    virtual ~FD2DEQCR();

    virtual MarketObjectSP GetMarket(const MarketData*    market,
                                    const string&        name,
                                    const CClassConstSP& type) const;

    /** interpolate/integrate for prices at t=0 with asset=asset(t=0) etc. */
    virtual double getPrice0( const TreeSlice & price ) const;

    /** Create a MarketDataFetcher which will be used for retrieving market data etc */
    MarketDataFetcherSP createMDF() const;

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * Returns riskMappingIrrelevant --- though you can maybe imagine
     * otherwise?  See IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

    /** get the initial conditions and corresponding indexes in the slices
        for forward induction */
    virtual void getInitialConditions(IntArray& initialIndex, DoubleArray& initialValue) const{};

protected:

    /** retrieving market data from MDF */
    virtual void retrieveFactor();


    /**-----------------------------------------------------*/
    /** Compute Long Run term structure for CIR model on
    /// the curve and then assumes it piecewise constant
    /// on the simulation dates
    ///-----------------------------------------------------*/

    void createTermStructureCIR(CashFlowArraySP defRatesDates);

    /**-----------------------------------------------------*/
    /** model initialisation (first part)
    /// get vol 
    ///-----------------------------------------------------*/

    virtual void initModel();

    /**-----------------------------------------------------*/
    /** model initialisation  (second part)
    /// most memory init should be here
    /// finalise model initialisation once after initModel() and product init() 
    /// called by validate 
    ///-----------------------------------------------------*/

    virtual void finaliseModel(CControl*    control);

    /**-----------------------------------------------------*/
    /** each derived model can set diff. boundaries 
    /// based on the dynamics of underlying,
    /// alpha is the input truncation1D (nb of std) 
    /// outLowB, outUpB are low and up boundaries for FD 
    ///-----------------------------------------------------*/

    virtual  void setFdBounds(double& volForBound1, double alpha1, double& outLowB1, double& outUpB1,
                            double& volForBound2, double alpha2, double& outLowB2, double& outUpB2);
    
    /**----------------------------------------------------------------*/
    /**  setup the coefficients of the PIDE needed to be solved
    ///  Assume the most general PDE is 
    ///  U_t + a*U_x + b*U_y + c*U_xx + d*U_yy + e*U_xy + f*U = 0;
    ///  this ft is called at each time step, so they're time depd.
    ///----------------------------------------------------------------*/

    /** assuming that DimX (vertical) is stock and DimY is Variance*/

    void pdeCoeff(int step, double*** coeff,  
                    int bot1, int top1, int bot2, int top2 );

    // ------------------ data ------------------------------------------

    ICDSParSpreadsConstSP            credit; // $unregistered
    CVolProcessedBSSP                volLN; // $unregistered
    DoubleArray                        vols; // $unregistered

    double correl; // $unregistered
    double vol_CIR; // $unregistered
    double mr_CIR; // $unregistered
    double init_CIR; // $unregistered
    DoubleArray lr_CIR; // $unregistered
};

typedef smartPtr<FD2DEQCR> FD2DEQCRSP;

DRLIB_END_NAMESPACE
#endif
