//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolProcessedDVF.hpp
//
//   Description : interface for processed local vol object
//
//   Date        : 30 Oct 2001
//
//
//----------------------------------------------------------------------------

#ifndef VOL_PROCESSED_DVF_HPP
#define VOL_PROCESSED_DVF_HPP

#include "edginc/VolProcessed.hpp"
#include "edginc/Lattice.hpp"

DRLIB_BEGIN_NAMESPACE

/** Type of CVolProcessed. Returned following a LocalVolRequest */
class MARKET_DLL CVolProcessedDVF: public CObject,
                        public virtual IVolProcessed{
public:
    static CClassConstSP const TYPE;
    
    virtual double computeImpVol(const DateTime& maturity,
                                 double          strike) const = 0;
    
    /** volatility calculator object using a cache of forward values */
    class MARKET_DLL IVolCalculator {
    public:
        /** Calculate local volatility at a particular time point */
        virtual void CalcLocVol(
            const CSliceDouble& strikes,        // Slice of strikes
            int                 iDate,          // index in the array of maturities
            CSliceDouble&       locVol) = 0;    // output
        
        virtual ~IVolCalculator(){}
    };  
    
    /** create a volatility calculator with a cache of forward values, storing the array of maturity dates */
    virtual IVolCalculator* CreateVolCalculator(
        const DateTimeArray&   maturity, // the time axis of the lattice
        bool                   isIntraDayInterp = true) const = 0;      
    
    virtual double CalcLocVol(const DateTime& maturity,
                              double          strike,
                              bool			  isIntraDayInterp=true) const = 0;
    
    virtual void CalcLocVol(CLatticeDouble*     strikes,
                            DateTimeArray&      maturities, // the time axis of the lattice above
                            CLatticeDouble*	locVol,
                            bool	        isIntraDayInterp=true) const = 0;
    
    /* old versions that do not use the cache (quite slow) */
    virtual void CalcLocVar(CLatticeDouble* 	strikes,
                            DateTimeArray& 	maturities, // the time axis of the lattice above
                            CLatticeDouble*	locVar,
                            bool                isIntraDayInterp=true) const = 0;
    
protected:
    CVolProcessedDVF(const CClassConstSP& clazz);
private:
    CVolProcessedDVF(const CVolProcessedDVF &rhs);
    CVolProcessedDVF& operator=(const CVolProcessedDVF& rhs);
};

typedef smartConstPtr<CVolProcessedDVF> CVolProcessedDVFConstSP;
typedef smartPtr<CVolProcessedDVF> CVolProcessedDVFSP;

DRLIB_END_NAMESPACE
#endif

