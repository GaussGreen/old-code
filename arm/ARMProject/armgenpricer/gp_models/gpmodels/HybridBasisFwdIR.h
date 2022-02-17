/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file HybridBasisFwdIR.h
 *
 *  \brief
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date April 2004
 */


#ifndef _INGPMODELS_HYBRIDBASISFWDIR_H
#define _INGPMODELS_HYBRIDBASISFWDIR_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
#include "gpbase/env.h"
#include "gpbase/port.h"

#include "MultiAssets.h"
#include "typedef.h"

/// forward declaration in the global namespace
class ARM_Forex;



CC_BEGIN_NAMESPACE( ARM )


//-----------------------------------------------------------------------------
// \class HybridBasisFwdIR
// \brief
//  Hybrid model for basis swap margin and two interest rate curves.
//  It is based on a single reference IR model (stochastic).
//  The other IR model uses a forward margin model (deterministic) and the
//  reference IR model
//  The basis margin also uses a forward margin and one of both IR models as
//  its reference model
//  The forex used a forward forex model. It gives deterministic forex value
//  equal to the forward forex at asOfDate.
//-----------------------------------------------------------------------------

class ARM_HybridBasisFwdIR :    public ARM_MultiAssetsModel
{
public:
    enum modelsAlias
    {
        RefModel=0,         // the stochastic IR model
        IrMarginModel,      // the forward margin model to generated the 2nd IR model
        BasisMarginModel,   // the forward margin model for basis swap
        ForexModel,         // the forex model

        NbModels
    };

    ARM_StringVectorVector* ValidatePricingModels(const ARM_StringVector& names,
        const vector< ARM_PricingModelPtr >& models,modelsAlias refBasisModelIdx,
        modelsAlias domesticModelIdx, modelsAlias foreignModelIdx);

	void UpdateLinks();

public:
	/// constructor, copy constructor, destructor, assignment operator
    ARM_HybridBasisFwdIR(const ARM_StringVector& names, const vector< ARM_PricingModelPtr >& models,
        modelsAlias refBasisModelIdx,modelsAlias domesticModelIdx, modelsAlias foreignModelIdx);
	ARM_HybridBasisFwdIR(const ARM_HybridBasisFwdIR& rhs);
	virtual ~ARM_HybridBasisFwdIR() {}
	ARM_HybridBasisFwdIR& operator = (const ARM_HybridBasisFwdIR& rhs);

    /// Utility to change the reference model
    void ReplaceRefModel(ARM_PricingModel* newRefModel);
    virtual void SetModelMap(ARM_ModelNameMap* RefModelMap);
	virtual void SetModelMapNoClone(ARM_ModelNameMap* RefModelMap);

    /// Update curves (clones the curve and uses smart pointors!)
	void UpdateCurves( const ARM_ZeroCurvePtr& refCurve, const ARM_ZeroCurvePtr irMarginCurve, 
		const ARM_ZeroCurvePtr& basisMarginCurve, const ARM_Forex& forex);

	virtual ARM_PricingStatesPtr Init( const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos );
	virtual ARM_PricingStatesPtr Induct(ARM_PricingStatesPtr& states,double toTime);
	virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;


	/// --------- access to the stochastic Model
    const ARM_PricingModel* GetStochasticModel() const;
    ARM_PricingModel* GetStochasticModel();


    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_HybridBasisFwdIR(*this); }
	virtual string toString(const string& indent="",const string& nextIndent="") const;
	virtual string ExportShortName() const { return "LBFIR";}
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
