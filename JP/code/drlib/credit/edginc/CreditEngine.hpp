//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CreditEngine.cpp
//
//   Description : Credit engine
//
//   Author      : Ning Shen
//
//   Date        : 10 June 2002
//
//
//----------------------------------------------------------------------------

#ifndef CREDIT_ENGINE_HPP
#define CREDIT_ENGINE_HPP

#include "edginc/CreditPathValuesIn.hpp"
#include "edginc/CreditPathValuesOut.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/RegressionTest.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/CreditSupport.hpp"
#include "edginc/MarginAcct.hpp"
#include "edginc/CreditSpreadCurve.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

/** credit engine */
class CREDIT_DLL CreditEngine : public CObject,
                 public ClientRunnable,
                 public IRegressionTest
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);

    static IObject* defaultCreditEngine(){
        return new CreditEngine();
    }

    // EdrAction version of addin
    virtual IObjectSP run();
    /** addin for calculate options */
    static IObjectSP computePathValuesAddin(const CreditEngine*);
    /** for regression run */
    virtual IObjectSP runTest() const;

    /** takes credit data, group instruments according to underliers,
        call calcOneUndPath() with an array of indicators that signal the instruments using the same underlier. */
    CreditPathValuesOutSP computePathValues() const;

protected:
    /** generate paths for one assset and call price() on all instruments with the same underlying underlier
        results are added to target */
    void calcOneUndPath(CreditPathValuesOutSP target,
                        map<string, int> & usedUnd,
                        const CIntArray& instIndex,
                        int seed, int market, double corr,
		                const DoubleArray& positions) const;

    /** create paths dates that require sample generation */
    void createPathDates(const DateTime & valueDate, DateTimeArraySP& pathDates, CIntArraySP& pathDateTypes) const;

    /** calculate fx fwds  */
    void calcFXFwd(const string& creditCcy, const string& instCcy,
                   const DateTimeArray& dates, DoubleArray& fx) const;

    /** get a credit support object  */
    CreditSupportSP getCreditSupport(CInstrument* inst) const;

    /** generate spot paths */
    void genPaths(CreditUndSP und,
                  const DateTimeArray& pathDates,
                  int numOfPaths,
                  const DateTime& valueDate,
                  int seed,
                  int marketIndex,
                  double corr,
                  DoubleArray& spot1,
                  DoubleArray& spot2,
                  DoubleArray& atmFwd,
                  DoubleArray& var_dt) const;

    /** returns the maximum of each instrument's last exposure date. */
    DateTime getLastExposureDateFromImnts() const;

private:
    CreditEngine();

    // input data
    CreditPathValuesInSP creditData;
    CMarketDataSP		 market;            /** market cache */
    MarginAcctSP         marginAcct;        /** margin collateral account */

    double               percentile;        /** percentile used for peak computation */
    string               creditCcyName;      /** name of currency we wish to compute exposure in */
	CreditSpreadCurveSP	 creditSpreads;		/** optional field: required if cvr is to be computed. */

	bool				 outputGrids;		/** if true, returns all values on all paths */

    bool                 useThetaFwdRate;   /** true=fwd rate used when shifting value date, false=whole yiled curve shift */

    IModelSP               model;             /** used for getting yield curve from market */
    // unregistered
    mutable vector< vector< float > > randMarket; // $unregistered

    // to do: do we need this?
    // for convenience
    mutable FXAssetSP     fxAsset; // $unregistered

};

DECLARE_REF_COUNT(CreditEngine);

DRLIB_END_NAMESPACE

#endif
