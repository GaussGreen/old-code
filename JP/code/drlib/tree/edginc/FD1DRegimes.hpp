#ifndef FD1DRegimes_HPP
#define FD1DRegimes_HPP

#include "edginc/FD1DMulti.hpp"
#include "edginc/MarketDataFetcherLN.hpp"
#include "edginc/MultiFactors.hpp"
#include "edginc/VolProcessedBS.hpp"
#include "edginc/MultiRegimeFactor.hpp"

DRLIB_BEGIN_NAMESPACE

class TREE_DLL FD1DRegimes : public FD1DMulti
{
public:
    static const CClassConstSP TYPE;
    static void load( CClassSP & type );

    FD1DRegimes( const CClassConstSP & type = TYPE );
    ~FD1DRegimes();

    static IObject * defaultFD1DRegimes()
    {
        return new FD1DRegimes();
    }

    /** validate some inputs */
    virtual void validatePop2Object();

    /** Override default createMDF in order to set the right MDF */
    virtual MarketDataFetcherSP createMDF() const;

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;
    
    /** get the initial conditions and corresponding indexes in the slices
        for forward induction */
    virtual void getInitialConditions(IntArray& initialIndex, DoubleArray& initialValue) const{};

	//temporary, to make barrier option work, to review
    virtual int GetStepVol(int step, vector<double>& vol, const double* s, int start, int end);
    
protected:
    /** retrieving market data from MDF */    
    virtual void retrieveFactor();
    
    virtual void initModel();
    
    virtual void finaliseModel(CControl*    control);

    /**-----------------------------------------------------
        each derived model can set diff. boundaries 
        based on the dynamics of underlying,
        alpha is the input truncation1D (nb of std) 
        outLowB, outUpB are low and up boundaries for FD 
    -----------------------------------------------------*/
    virtual void setFdBounds(
        double alpha,
        double & outLowB,
        double & outUpB ) const;

	virtual void pdeCoeff(
        int step,
        DoubleArray & a,
        DoubleArray & c,
        DoubleArray & f,
        DoubleArray & g,
		DoubleArray & q,
		DoubleArray & jump) const;

    // log-normal vol
    string volType; // $unregistered

private:
	/** Invoked after instrument has got its market data. */
    virtual void getMarket(const MarketData* market, IInstrumentCollectionSP instruments);  

	//transient fields ...
	MultiRegimeFactorSP		stochFactor;
};

DRLIB_END_NAMESPACE

#endif
