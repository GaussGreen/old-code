//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : Tranche Option.hpp
//
//   Description : A pricer for SCID model 
//
//   Author      : Adrian Bozdog
//
//----------------------------------------------------------------------------

#ifndef QR_TRANCHEOPTION_HPP
#define QR_TRANCHEOPTION_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/CDOIndexOption.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(TrancheOption)

class PRODUCTS_DLL TrancheOption :    public CDOIndexOption, 
                     public virtual SCID::IIntoProduct {
public:
    static CClassConstSP const TYPE;

    virtual ~TrancheOption();
    /*=========================================================================
     * I/CInstrument Interface
     *=======================================================================*/
    void Validate();

    /*=========================================================================
     * Instrument specific functions
     *=======================================================================*/
    void setFastMCParameters(SCIDparametersSP &sCIDparamSP);
    double getLossAndLegs(long expiryIndexInSim, DateTimeArray & simDates, DateTimeArray& maturities, DateTime& expiry, DoubleArray& RA, DoubleArray& DL, SCID * model); 


private:
    TrancheOption(const TrancheOption& rhs);
    TrancheOption& operator=(const TrancheOption& rhs);
	TrancheOption(CClassConstSP clazz = TYPE);

    /*=========================================================================
     * DATA FIELDS
     *=======================================================================*/
    double    lowStrike;        /** lowStrike and highStrike defines portion of
                                    portfolio that payoff depends upon */
    double    highStrike;

    /*=========================================================================
     * FRIENDS
     *=======================================================================*/
    /** for reflection **/
    friend class TrancheOptionHelper;
};

DRLIB_END_NAMESPACE
#endif




