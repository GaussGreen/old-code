
/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	 /
 *																							 *
 *			Class Hybrid Inf Ir Leg															 *
 *																							 *
 *			This class builds a hybrid inf ir leg from swap legs and  Inf Swap Leg			 *									 *
 *																							 *
 *			Author	: Mathieu Bernardo														 *
 *			Date	: April, 2nd 2007														 *																											 *
 *			Copyright (c) NatIxis															 *
 *			Version	: 1.0																	 *
 *																							 *
 *	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

#ifndef _INGPINFLATION_HYBRIDINFIRPRICER_H
#define _INGPINFLATION_HYBRIDINFIRPRICER_H

#pragma warning(disable : 4786) 
#include <gpinflation/hybridinfleg.h>
#include <gpinflation/infvisitor.h>
#include <gpinflation/infpayoffvisitor.h>
#include <gpinflation/infmodelvisitor.h>

CC_BEGIN_NAMESPACE( ARM )

 
class ARM_HybridInfLegPricer: public ARM_InfPricer{

public:
	ARM_HybridInfLegPricer( ARM_HybridInfIrLeg* ins, ARM_InfHybridPayOff* payOff, ARM_InfModel* model=NULL );
	ARM_HybridInfLegPricer(	const ARM_HybridInfLegPricer & );
	virtual ~ARM_HybridInfLegPricer(){};
	
	ASSIGN_OPERATOR		( ARM_HybridInfLegPricer			);

	virtual ARM_Object* Clone() const { return new ARM_HybridInfLegPricer(*this); }
 	virtual string toString(	const string& indent="", const string& nextIndent="") const;

	virtual void Compute();

private:
	ARM_MAP_Double	CptFwd( const int & i );
	ARM_MAP_VolPar	CptVol( const int & i );
	ARM_MAP_PairDb	CptCor( const int & i );

	void			CptFlows();

	void			Map	( ARM_InfPayOffValuePtr& , ARM_InfModelValuePtr& );

private:
	bool									isComputed;
	ARM_CountedPtr<ARM_InfHybridPayOff>		itsPayOff;
	ARM_CountedPtr<ARM_InfModel>			itsModel;
	ARM_CountedPtr<ARM_HybridInfIrLeg>		itsSecurity;
};




CC_END_NAMESPACE()

#endif
/*--------------------------------------------------------------------------*/
/*---- End of file ----*/


