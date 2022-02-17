//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : ImpliedLossAddIn.hpp
//
//   Description :	Add in that returns par spreads, annuities or losses generated 
//					by ImpliedLossModel
//                 
//   Author      : Matthias Arnsdorf
//
//   Date        : January 2006
//
//
//----------------------------------------------------------------------------

#ifndef IMPLIED_LOSS_ADD_IN_HPP
#define IMPLIED_LOSS_ADD_IN_HPP

#include "edginc/ImpliedLossModel.hpp"

DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL ImpliedLossAddIn : public CObject
{
public:
    static CClassConstSP const TYPE;

    
    virtual IObjectSP run();

private:

	
    static void load(CClassSP& clazz);

    static IObject* defaultImpliedLossAddIn();

    //for reflection
    ImpliedLossAddIn();


    // FIELDS -------------------------------

	/** market data */
	MarketDataSP market;

	/** CDO quotes */
	CDOQuotesWrapper cdoQuotes;

	/** Index swap spreads */
	ICDSParSpreadsWrapper indexSpreads;

	/** interpolation model */
	ImpliedLossModelSP model;

	/**Output expiries */
	ExpiryArraySP expiries;

	/** Output low strikes */
	DoubleArraySP lowStrikes;

	/** Output high strikes */
	DoubleArraySP highStrikes;

	/** outputTypes = "SPREAD", "UPFRONT", "ANNUITY" , "RISKY_ZERO" or "CONT_LEG" */
	StringArraySP outputTypes;

	/** coupons for upfront */
	DoubleArraySP coupons;
};

DRLIB_END_NAMESPACE

#endif



