//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : ICDOFineGrid.hpp
//
//   Description : Interface for "CDO fine grids" responsible for
//                 "extending" market quotes
//
//   Author      : Antoine Gregoire
//
//   Date        : January 2006
//
//----------------------------------------------------------------------------


#ifndef I_CDO_FINE_GRID_HPP
#define I_CDO_FINE_GRID_HPP

#include "edginc/CDOQuotes.hpp"

DRLIB_BEGIN_NAMESPACE

/** Generic definition of a grid */
class PRODUCTS_DLL ICDOFineGrid: public virtual IObject {
public:
    /**
     * Extend the input market quotes into a new set of quotes
     * corresponding to this CDO fine grid
     * */
    virtual CDOQuotesSP extendQuotes(CDOQuotesConstSP marketQuotes) = 0;

	/** get market */
	virtual void getMarket(MarketData *market) = 0;


    /** TYPE for IBootstrappable */
    static CClassConstSP const TYPE;
};

typedef smartPtr<ICDOFineGrid> ICDOFineGridSP;

DRLIB_END_NAMESPACE

#endif /*I_CDO_FINE_GRID_HPP*/

