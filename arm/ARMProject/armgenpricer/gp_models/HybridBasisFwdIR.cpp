/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file HybridBasisFwdIR.cpp
 *
 *  \brief
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date April 2004
 */


#include "gpbase/removeidentifiedwarning.h"
#include "gpmodels/hybridbasisfwdir.h"

/// gpbase
#include "gpbase/autocleaner.h"

/// gpinfra
#include "gpinfra/modelnamemap.h"
#include "gpinfra/pricingstates.h"


/// gpmodels
#include "gpmodels/forwardmargin.h"
#include "gpmodels/forwardforex.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : HybridBasisFwdIR
///	Routine: Constructor with a pricing models vector
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_HybridBasisFwdIR::ARM_HybridBasisFwdIR( 
	const ARM_StringVector& names, 
	const vector< ARM_PricingModelPtr>& models,
    modelsAlias refBasisModelIdx, 
	modelsAlias domesticModelIdx, 
	modelsAlias foreignModelIdx )
:
	ARM_MultiAssetsModel()
{
    /// Check linked models consistency
    if(models.size() != NbModels)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : not enough models to built a Hybrid Basis/Fwd IR model");

    double asOfDate=models[RefModel]->GetZeroCurve()->GetAsOfDate().GetJulian();
    if( models[IrMarginModel]->GetZeroCurve()->GetAsOfDate().GetJulian() != asOfDate ||
        models[BasisMarginModel]->GetZeroCurve()->GetAsOfDate().GetJulian() != asOfDate)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Only one AsOfDate is allow for all linked models");

    SetName(ARM_BASISFWDIRMODEL);

    /// Check models consistency and build links
    ARM_StringVectorVector* modelLinks = ValidatePricingModels(names,models,
        refBasisModelIdx,domesticModelIdx,foreignModelIdx);

    ARM_AutoCleaner< ARM_StringVectorVector > HoldSVV(modelLinks);

    /// Built the model name map (no clone of the models)
    ARM_ModelNameMap* modelMap = new ARM_ModelNameMap(names,models,*modelLinks);
	(*modelMap)[IrMarginModel]->SetUsedInPricing( false );
	(*modelMap)[BasisMarginModel]->SetUsedInPricing( false );
	(*modelMap)[ForexModel]->SetUsedInPricing( false );
	SetModelMapNoClone(modelMap);

    /// The numerical method is set to the IR reference model one
    SetNumMethod(models[RefModel]->GetNumMethod());

    /// The Zc curve is set to the IR reference model one
    SetZeroCurve(models[RefModel]->GetZeroCurve());

	/// update links
	UpdateLinks();
	ARM_PricingModel* refModel = &*(*modelMap)[RefModel]->Model();
	SetRefModel( refModel );

	/// Get the numeraire from the refModel
	ARM_PricingModel::SetNumeraire( refModel->GetNumeraire() );
}



////////////////////////////////////////////////////
///	Class  : ARM_HybridBasisFwdIR
///	Routine: Induct
///	Returns: 
///	Action : Call the reference model
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_HybridBasisFwdIR::Induct(ARM_PricingStatesPtr& states,double toTime)
{
	ARM_ModelNameMap* modelMap = GetModelMap();
    return ((*modelMap)[RefModel])->Model()->Induct(states,toTime);
}


////////////////////////////////////////////////////
///	Class  : ARM_HybridBasisFwdIR
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_HybridBasisFwdIR::ARM_HybridBasisFwdIR(const ARM_HybridBasisFwdIR& rhs)
:	ARM_MultiAssetsModel(rhs)
{ 
	ARM_ModelNameMap* modelMap = GetModelMap();

    /// In forward margin models, update links their reference model
    ARM_ForwardMargin* irMargin = static_cast< ARM_ForwardMargin* >(&*((*modelMap)[IrMarginModel]->Model()));
    irMargin->SetRefModel(&*(*modelMap)[RefModel]->Model()); /// linked model is always the reference model

    ARM_ForwardMargin* basisMargin = static_cast< ARM_ForwardMargin* >(&*((*modelMap)[BasisMarginModel]->Model()));
    basisMargin->SetRefModel(&*(*modelMap)[(*modelMap)[BasisMarginModel]->OtherModelRefNb()[0]]->Model());
}



////////////////////////////////////////////////////
///	Class  : HybridBasisFwdIR
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_HybridBasisFwdIR& ARM_HybridBasisFwdIR::operator=(const ARM_HybridBasisFwdIR& rhs)
{
	if(this != &rhs)
	{
		ARM_MultiAssetsModel::operator=(rhs);
		UpdateLinks(); 
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_HybridBasisFwdIR
///	Routine: ValidatePricingModels
///	Returns: 
///	Action : Check the models consistency
////////////////////////////////////////////////////
ARM_StringVectorVector* ARM_HybridBasisFwdIR::ValidatePricingModels(
	const ARM_StringVector& names,
    const vector<ARM_PricingModelPtr>& models, 
	modelsAlias refBasisModelIdx,
    modelsAlias domesticModelIdx, 
	modelsAlias foreignModelIdx)
{
    /// Check the model consistency
    if(models.size() != NbModels || names.size() != NbModels )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : the Hybrid Basis Forward & IR model needs 3 models");

    if(refBasisModelIdx >= NbModels)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : the reference model of the basis forward model is not available");

    if( !dynamic_cast< ARM_PricingModelIR* >(&*(models[RefModel])) ||
        !dynamic_cast< ARM_ForwardMargin* >(&*(models[IrMarginModel])) ||
        !dynamic_cast< ARM_ForwardMargin* >(&*(models[BasisMarginModel])) ||
        !dynamic_cast< ARM_ForwardForex* >(&*(models[ForexModel])) )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : the model needs an IR model, two forward margin models (2nd IR curve & basis swap) and a forward forex model");

    /// Place the nested model
    static_cast< ARM_ForwardMargin& >(*(models[IrMarginModel])).SetRefModel( &*models[RefModel] );
    static_cast< ARM_ForwardMargin& >(*(models[BasisMarginModel])).SetRefModel( &*models[refBasisModelIdx] );

    /// Build links to auxillary models (but auxillary models are not used
    /// because forward margin or forward forex models already get a pointor to linked model)
    ARM_StringVectorVector* modelLinks  = new ARM_StringVectorVector(NbModels);
    (*modelLinks)[IrMarginModel]        = ARM_StringVector(1,names[RefModel]);
    (*modelLinks)[BasisMarginModel]     = ARM_StringVector(1,names[refBasisModelIdx]);
    (*modelLinks)[ForexModel]           = ARM_StringVector(2);
    (*modelLinks)[ForexModel][0]        = names[domesticModelIdx];
    (*modelLinks)[ForexModel][1]        = names[foreignModelIdx];

    return modelLinks;
}




////////////////////////////////////////////////////
///	Class  : ARM_HybridBasisFwdIR
///	Routine: UpdateLinks
///	Returns: void
///	Action : udpate the various links for a given model map
////////////////////////////////////////////////////

void ARM_HybridBasisFwdIR::UpdateLinks()
{
	ARM_ModelNameMap* modelMap = GetModelMap();

	/// update links
	ARM_PricingModelPtr& newRefModel = (*modelMap)[RefModel]->Model();

	ARM_ForwardMargin* irMargin = static_cast< ARM_ForwardMargin* >(&*((*modelMap)[IrMarginModel]->Model()));
	irMargin->SetRefModel(&*newRefModel); /// linked model is always the reference model
	
	ARM_ForwardMargin* basisMargin = static_cast< ARM_ForwardMargin* >(&*((*modelMap)[BasisMarginModel]->Model()));
	if( strcmp( basisMargin->GetZeroCurve()->GetCurrencyUnit()->GetCcyName(), irMargin->GetZeroCurve()->GetCurrencyUnit()->GetCcyName() ) == 0 )
		basisMargin->SetRefModel(irMargin);
	else
		basisMargin->SetRefModel(&*newRefModel);

	/// update the numMethod
	ARM_NumMethodPtr numMethod = newRefModel->GetNumMethod();
	ARM_PricingModel::SetNumMethod( numMethod );
}


////////////////////////////////////////////////////
///	Class  : ARM_HybridBasisFwdIR
///	Routine: SetModelMap,SetModelMapNoClone
///	Returns: nothing
///	Action : ReplaceModelMap and updates links
////////////////////////////////////////////////////
void ARM_HybridBasisFwdIR::SetModelMap(ARM_ModelNameMap* RefModelMap)
{
	ARM_MultiAssetsModel::SetModelMap(RefModelMap);
	UpdateLinks();
}


void ARM_HybridBasisFwdIR::SetModelMapNoClone(ARM_ModelNameMap* RefModelMap)
{
	ARM_MultiAssetsModel::SetModelMapNoClone(RefModelMap);
	UpdateLinks();
}


////////////////////////////////////////////////////
///	Class  : ARM_HybridBasisFwdIR
///	Routine: ReplaceRefModel
///	Returns: nothing
///	Action : Replace the reference model by the input one
////////////////////////////////////////////////////
void ARM_HybridBasisFwdIR::ReplaceRefModel(ARM_PricingModel* newRefModel)
{
	ARM_ModelNameMap* modelMap = GetModelMap();
    ARM_PricingModelPtr oldRefModel = (*modelMap)[RefModel]->Model();


	newRefModel->SetModelName(oldRefModel->GetModelName());
    (*modelMap)[RefModel]->Model() = ARM_PricingModelPtr(newRefModel);

    ARM_ForwardMargin* irMargin = static_cast< ARM_ForwardMargin* >(&*((*modelMap)[IrMarginModel]->Model()));
    irMargin->SetRefModel(newRefModel);		/// linked model is always the reference model

	/// compared to update links, we need to know if we need to chang the basis Margin model reference model
    ARM_ForwardMargin* basisMargin = static_cast< ARM_ForwardMargin* >(&*((*modelMap)[BasisMarginModel]->Model()));
    if(basisMargin->GetRefModel() == &*oldRefModel)
        basisMargin->SetRefModel(newRefModel);
}


////////////////////////////////////////////////////
///	Class  : ARM_HybridBasisFwdIR
///	Routine: UpdateCurves
///	Returns: nothing
///	Action : Update zero curves of the model and forex
////////////////////////////////////////////////////
void ARM_HybridBasisFwdIR::UpdateCurves(
	const ARM_ZeroCurvePtr& refCurve,
	const ARM_ZeroCurvePtr irMarginCurve,
	const ARM_ZeroCurvePtr& basisMarginCurve,
    const ARM_Forex& forex)
{
	ARM_ModelNameMap* modelMap = GetModelMap();

    SetZeroCurve(refCurve);
    (*modelMap)[RefModel]->Model()->SetZeroCurve(refCurve);
    (*modelMap)[IrMarginModel]->Model()->SetZeroCurve(irMarginCurve);
    (*modelMap)[BasisMarginModel]->Model()->SetZeroCurve(basisMarginCurve);

    static_cast< ARM_ForwardForex&>( *((*modelMap)[ForexModel]->Model()) ).SetForex(forex);

    /// Set forex domestic curve
    static_cast< ARM_ForwardForex&>( *((*modelMap)[ForexModel]->Model()) ).SetCurves(
        (*modelMap)[(*modelMap)[ForexModel]->OtherModelRefNb()[0]]->Model()->GetZeroCurve(),
        (*modelMap)[(*modelMap)[ForexModel]->OtherModelRefNb()[1]]->Model()->GetZeroCurve() );
}


////////////////////////////////////////////////////
///	Class  : ARM_HybridBasisFwdIR
///	Routine: UpdateCurves
///	Returns: nothing
///	Action : Update zero curves of the model and forex
////////////////////////////////////////////////////
string ARM_HybridBasisFwdIR::toString(const string& indent, const string& nextIndent) const
{
	const ARM_ModelNameMap* const modelMap = GetModelMap();

    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "Forward Basis Margin, Forward Rate Margin & IR Hybrid Model\n";
    os << indent << "-----------------------------------------------------------\n";

    ARM_PricingModelPtr refModel = (*modelMap)[RefModel]->Model();

    /// Set model flags to avoid multiple full dumps of the reference model
    ARM_PricingModelPtr irMarginModel = (*modelMap)[IrMarginModel]->Model();
    bool prevIrMarginDump(static_cast< ARM_ForwardMargin& >(*irMarginModel).GetRefModelDump());
    static_cast< ARM_ForwardMargin& >(*irMarginModel).SetRefModelDump(false);

    ARM_PricingModelPtr basisMarginModel = (*modelMap)[BasisMarginModel]->Model();
    bool prevBasisMarginDump(static_cast< ARM_ForwardMargin& >(*basisMarginModel).GetRefModelDump());
    if( static_cast< ARM_ForwardMargin& >(*basisMarginModel).GetRefModel() == &*refModel )
        static_cast< ARM_ForwardMargin& >(*basisMarginModel).SetRefModelDump(false) ;

    /// Dump the model map (and its linked models)
    os << indent << modelMap->toString(indent,nextIndent);

    /// Restore previous flags
    static_cast< ARM_ForwardMargin& >(*irMarginModel).SetRefModelDump(prevIrMarginDump);
    static_cast< ARM_ForwardMargin& >(*basisMarginModel).SetRefModelDump(prevBasisMarginDump);

    return os.str();
}

////////////////////////////////////////////////////
///	Class  : ARM_HybridBasisFwdIR
///	Routine: Init
///	Returns: ARM_PricingStatesPtr
///	Action : Init the model
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_HybridBasisFwdIR::Init( const string& payModelName, 
	const ARM_TimeInfoPtrVector& timeInfos )
{
	const ARM_ModelNameMap* const modelMap = GetModelMap();
    ARM_PricingModel* refModel = &*(*modelMap)[payModelName]->Model();
	SetRefModel( refModel );
	UpdateLinks();
	return refModel->Init( payModelName, timeInfos );
}


////////////////////////////////////////////////////
///	Class  : ARM_HybridBasisFwdIR
///	Routine: GetStochasticModel
///	Returns: const ARM_PricingModel*
///	Action : Returns the stochastic model
////////////////////////////////////////////////////
const ARM_PricingModel* ARM_HybridBasisFwdIR::GetStochasticModel() const
{
	const ARM_ModelNameMap* const modelMap = GetModelMap();
	return &*(*modelMap)[RefModel]->Model();
}

////////////////////////////////////////////////////
///	Class  : ARM_HybridBasisFwdIR
///	Routine: GetStochasticModel
///	Returns: const ARM_PricingModel*
///	Action : Returns the stochastic model
////////////////////////////////////////////////////
ARM_PricingModel* ARM_HybridBasisFwdIR::GetStochasticModel()
{
	ARM_ModelNameMap* modelMap = GetModelMap();
	return &*(*modelMap)[RefModel]->Model();
}

////////////////////////////////////////////////////
///	Class  : ARM_HybridBasisFwdIR
///	Routine: GetStochasticModel
///	Returns: const ARM_PricingModel*
///	Action : Returns the stochastic model
////////////////////////////////////////////////////
void ARM_HybridBasisFwdIR::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
	const ARM_ModelNameMap* const modelMap = GetModelMap();
	(*modelMap)[RefModel]->Model()->MCModelStatesFromToNextTime(states,timeIndex);
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

