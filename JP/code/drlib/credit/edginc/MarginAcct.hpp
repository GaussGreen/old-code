//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MarginAcct.hpp
//
//   Description : Parameters for a Margin Collateral Account
//
//   Author      : Jay Blumenstein
//
//   Date        : 13 Sep 2002
//
//
//----------------------------------------------------------------------------

#ifndef CREDIT_MARGIN_ACCT__H
#define CREDIT_MARGIN_ACCT__H

#include "edginc/config.hpp"
#include "edginc/Object.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/Model.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/Holiday.hpp"

DRLIB_BEGIN_NAMESPACE

//////////////////////////////////////////////////////////////////////
//
//  class MarginAcct
//
//  contains: parameters for a Margin Collateral Account.  Cf. officialbuild\credit\src\ecvgrid.c
//  methods: 
//////////////////////////////////////////////////////////////////////
class CREDIT_DLL MarginAcct: public CObject
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);
	
    static IObject* defaultMarginAcct(){
        return new MarginAcct();
    }

    virtual void validatePop2Object();

    double        upfrontMargin;       /**  var. margin in account already  */
    double        initialMargin;       /**  init margin acct. Never repaid  */
    double        thresholdMarginCpty; /**  unsecured margin counter party  */
    double        thresholdMarginJPM;  /**  unsecured margin JPM is never asked for */
    double        minMarginCall;       /**  (applied after threshold)       */
    int           marginCallDelay;     /**  days btn margin call and receipt*/
    bool          rehypothecate;       /**  if counter party can use margin */
	HolidaySP	  marginCallHol;	   /** optional field: used to find the last MTM date for an exposure date */

private:

	MarginAcct();
    MarginAcct(const MarginAcct& rhs);
    MarginAcct& operator=(const MarginAcct& rhs);

};

typedef smartPtr<MarginAcct> MarginAcctSP;

DRLIB_END_NAMESPACE

#endif
