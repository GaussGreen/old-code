//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FD2DSV.hpp
//
//   Description : two factor finite difference engine for SVCJ vol
//                    asset factor   dX = drift*dt + sqrt(V) dW + Merton jump in X
//                    vol factor  dV = k(V_0 - V)dt + vol*sqrt(V) dZ + Merton jump in V
//                    V_0 is allowed to be a gamma distribution of mean v0.
//
//   Author      : Ning Shen
//                   Xiaolan Zhang
//
//   Date        : November 29, 2004
//
//----------------------------------------------------------------------------

#ifndef FD2DSV_HPP
#define FD2DSV_HPP
#include "edginc/VolSV.hpp"
#include "edginc/FD2D.hpp"
#include "edginc/MarketDataFetcherLN.hpp"

DRLIB_BEGIN_NAMESPACE

class TREE_DLL FD2DSV : public FD2D{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& );
    static IObject* defaultFD2DSV(){
        return new FD2DSV();
    }
    
    friend class FD2DSolver;
    friend class FD2DSolverJumps;

    //FD2DSV();
    FD2DSV( const CClassConstSP & type = TYPE );

    virtual ~FD2DSV();

    /** Create a MarketDataFetcher which will be used for retrieving market data etc */
    MarketDataFetcherSP createMDF() const;
    
    
    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * Returns allowRiskMapping ? riskMappingAllowed : riskMappingDisallowed.
     * See IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;


    /** interpolate/integrate for prices at t=0 with asset=asset(t=0) etc. */
    virtual double getPrice0( const TreeSlice & price ) const;


    /** get the initial conditions and corresponding indexes in the slices
        for forward induction */
    virtual void getInitialConditions(IntArray& initialIndex, DoubleArray& initialValue) const;

protected:

    /** retrieving market data from MDF */
    virtual void retrieveFactor();

    
    /**-----------------------------------------------------
           model initialisation  (second part)
           most memory init should be here
           finalise model initialisation once after initModel() and product init() 
           called by validate 
           -----------------------------------------------------*/

    virtual void finaliseModel(CControl*    control);

    /**-----------------------------------------------------
           each derived model can set diff. boundaries 
           based on the dynamics of underlying,
           alpha is the input truncation1D (nb of std) 
           outLowB, outUpB are low and up boundaries for FD 
           -----------------------------------------------------*/

    virtual  void setFdBounds(double& volForBound1, double alpha1, double& outLowB1, double& outUpB1,
                              double& volForBound2, double alpha2, double& outLowB2, double& outUpB2);
    
    //virtual   void SetSpaceGridData();
    
    /**----------------------------------------------------------------
       setup the coefficients of the PIDE needed to be solved
       Assume the most general PDE is 
       U_t + a*U_x + b*U_y + c*U_xx + d*U_yy + e*U_xy + f*U + Jump(x,y) = 0;
       this ft is called at each time step, so they're time depd.
       ----------------------------------------------------------------*/
    
    /** assuming that DimX (vertical) is stock and DimY is Variance*/
    
    void pdeCoeff(int step, double*** coeff,  
                  int bot1, int top1, int bot2, int top2 );
    
    /**-------------------------------------
       for SVCJ, not in the base FD2D model 
       -------------------------------------*/
    //jumps parts
    double jumpProDen(double x, double y);
    
    VolSVConstSP          vol; // $unregistered
    bool                    allowRiskMapping;
};

typedef smartPtr<FD2DSV> FD2DSVSP;

DRLIB_END_NAMESPACE
#endif
