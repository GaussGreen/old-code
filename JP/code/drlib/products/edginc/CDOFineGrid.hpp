//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : CDOFineGrid.hpp
//
//   Description :  concrete CDO fine grid class responsible for
//                 "extending" market quotes
//
//   Author      : Matthias Arnsdorf
//
//   Date        : January 2006
//
//----------------------------------------------------------------------------


#ifndef CDO_FINE_GRID_HPP
#define CDO_FINE_GRID_HPP

#include "edginc/ICDOFineGrid.hpp"
#include "edginc/ImpliedLossModel.hpp"

DRLIB_BEGIN_NAMESPACE

/** Generic definition of a grid */
class PRODUCTS_DLL CDOFineGrid:	public CObject,
					public virtual ICDOFineGrid
 {
public:

	
	/** Destructor */
    virtual ~CDOFineGrid();

	/** Called immediately after object constructed */
    virtual void validatePop2Object();    

	/** get market */
	virtual void getMarket(MarketData *market);

    /**
     * Extend the input market quotes into a new set of quotes
     * corresponding to this CDO fine grid
     * */
    virtual CDOQuotesSP extendQuotes(CDOQuotesConstSP marketQuotes);

    /** TYPE for reflection */
    static CClassConstSP const TYPE;

protected:
	
	/** private constructor */
	CDOFineGrid();

private:

	// For reflection
    static void load (CClassSP& clazz);

	static IObject* defaultCDOFineGrid();


	/** default value dor running coupon */
	static const double DEFAULT_COUPON;

	// FIELDS ---------------------------------------

	/** expiries of fine grid quotes */
	ExpiryArraySP expiries;

	/** low strikes of fine qrid quotes */
	DoubleArraySP lowStrikes;

	/** high strikes of fine grid quotes */
	DoubleArraySP highStrikes;

	/** output types for quote (SPREAD or UPFRONT) */
	StringArraySP outputTypes;

	/** running copuon if outputType is UPFRONT */
	DoubleArraySP coupons;

	/** index swap spreads for CDO to interpolate */
	ICDSParSpreadsWrapper indexSpreads;

	/** interpolation model */
	ImpliedLossModelSP model;
};



DRLIB_END_NAMESPACE

#endif /*CDO_FINE_GRID_HPP*/

