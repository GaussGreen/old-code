//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolMertonProcessed.hpp
//
//   Description : interface for VolMertonProcessed object
//
//   Date        : 20 Fev 2004
//
//
//----------------------------------------------------------------------------

#ifndef VolMertonLV_PROCESSED_HPP
#define VolMertonLV_PROCESSED_HPP

#include "edginc/VolProcessed.hpp"
#include "edginc/VolProcessedDVF.hpp"
#include "edginc/VolSurface.hpp"
#include "edginc/VolParam.hpp"
#include "edginc/VolProcessedDVFParam.hpp"
#include "edginc/VolMertonLV.hpp"

DRLIB_BEGIN_NAMESPACE

/********* class VolMertonLVProcessed ****/
class MARKET_DLL VolMertonLVProcessed: public CVolProcessedDVFParam{
public:
    static CClassConstSP const TYPE;
	static void load(CClassSP& clazz);
    static IObject* defaultVolMertonLVProcessed(){
        return new VolMertonLVProcessed();
    }

    //// construct a vol processed DVF for parameterised vol
    VolMertonLVProcessed(const CVolBase*         vol,
                          const CVolParamConstSP&  volParam,
                          const CVolRequestDVF*    volRequest, 
                          const CAsset*            asset,
                          const VolSurface*        VolSurf,
						  VolMertonLVSP			   mertonLV);

    double getMertonParam(const string& param_name) const;
	VolSurfaceSP getMertonImpliedSurface() const;

	//////////form Olivier Brochaus Merton Model////
	void Quantile(
    const DateTime &date1,
    const DateTime &date2,
    double         epsilon,
    int            maxJumps,
    double*        quantile,
    int*           nbJumps ) const;

    /** Calculates jump factor between dates. */
    double CalcJump(double noise, int numJumps) const;

protected:
    VolMertonLVProcessed(const CClassConstSP& clazz);
    // all fields are transient
    VolMertonLVSP			myVolMertonLV;
    AssetConstSP            asset;
	CVolParamConstSP		myParamVol;
	YieldCurveConstSP		equityYC;

    DateTime startDate;

private:

    virtual bool computeLocV2(double  yrsToMat,
		                      //const DateTime&  maturity,
                              double  strike,
                              double  forward,
                              double  growthRate,
                              SImpV&  impV,
                              double  &v2) const;

    friend class VolMertonLV;
    friend class VolMertonLVCalib;
    VolMertonLVProcessed();
    VolMertonLVProcessed(const VolMertonLVProcessed &rhs);
    VolMertonLVProcessed& operator=(const VolMertonLVProcessed& rhs);
};

typedef smartConstPtr<VolMertonLVProcessed> VolMertonLVProcessedConstSP;
typedef smartPtr<VolMertonLVProcessed> VolMertonLVProcessedSP;

DRLIB_END_NAMESPACE
#endif
