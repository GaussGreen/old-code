/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file Q1F_Fx.cpp
 *
 *  \brief Q model 1 factor FX version
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date December 2004
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpmodels/Q1F_Fx.h"

/// gpinfra
#include "gpinfra/zccrvfunctor.h"
#include "gpinfra/pricingmodeltype.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/modelnamemap.h"
#include "gpinfra/pricingcontext.h"
#include "gpinfra/nummethod.h"
#include "gpinfra/discretisationscheme.h"

/// gpbase
#include "gpbase/utilityport.h"
#include "gpbase/interpolatorvector.h"
#include "gpbase/vectormanip.h"
#include "gpbase/comparisonfunctor.h"

/// gpmodel
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/ModelParamsQ1F.h"
#include "gpmodels/Q1F_PricingContext.h"
#include "gpmodels/Q1F.h"
#include "gpmodels/QModelAnalytics.h"
#include "gpmodels/ForwardMargin.h"
#include "gpmodels/multiassets.h"
#include "gpmodels/bs_modelparams.h"

/// gpnummethods
#include "gpnummethods/pde3Fnumericalschemes.h"
#include "gpnummethods/pdemethod.h"
#include "gpnummethods/tree1D.h"
#include "gpnummethods/markoviandriftsampler.h"
#include "gpnummethods/meanrevertingsampler.h"

/// gpclosedforms
#include "gpclosedforms/vanille_bs_formula.h"
#include "gpclosedforms/gaussian_integrals.h"

/// gpnumlib
#include "gpnumlib/gaussiananalytics.h"

/*** for file dump ***
#include "gpnummethods/treebase.h"
#include "gpnummethods/slice.h"
#include "gpnummethods/treeindex.h"
*** for file dump ***/

/// ARM Kernel
#include <inst/portfolio.h>
#include <inst/optionportfolio.h>

CC_BEGIN_NAMESPACE( ARM )

const size_t GL_PY_NBPOINTS         = 2;
const size_t GL_MIN_NBPOINTS        = 2;
const size_t DECAY_NBSTEP_PER_YEAR  = 4;

const double QFX_LIMIT      = 1.0e-3;
const double LNVOLFX_LIMIT  = 2.0; // 200%
const double QFXDIFF_LIMIT	= 1.0e-1;



////////////////////////////////////////////////////
///	Class  : ARM_QModel1F_Fx
///	Routine: Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_QModel1F_Fx::ARM_QModel1F_Fx(const ARM_ZeroCurvePtr& zc, 
								 ARM_ModelParamsQ1F_Fx* modelParam,
								 const ARM_CurveMatrix& correlMatrix,
								 bool integratedVersion,
								 bool isForexDiffusion)
:	ARM_EqFxBase(zc,modelParam,correlMatrix),
	itsPrecomputedFwds(NULL),
	itsIntegratedVersion(integratedVersion),
	itsIsForexDiffusion(isForexDiffusion),
	itsIRDomModel(NULL),
	itsIRForModel(NULL),
    itsIsImpliedVolCalc(false),
    itsFwdFxVol(0),
	itsLocalFunctional(NULL)
{
	if (!zc.IsNull() && modelParam)
		ARM_EqFxBase::Init();
}


////////////////////////////////////////////////////
///	Class  : ARM_QModel1F_Fx
///	Routine: Copy Constructor
///	Returns: 
///	Action : Copy Constructor
////////////////////////////////////////////////////

ARM_QModel1F_Fx::ARM_QModel1F_Fx( const ARM_QModel1F_Fx& rhs )
:	
	ARM_EqFxBase(rhs),
	itsPrecomputedFwds(rhs.itsPrecomputedFwds),
	itsIntegratedVersion(rhs.itsIntegratedVersion),
	itsIsForexDiffusion(rhs.itsIsForexDiffusion),
	itsIRDomModel( rhs.itsIRDomModel ),
	itsIRForModel( rhs.itsIRForModel ),
    itsIsImpliedVolCalc( rhs.itsIsImpliedVolCalc ),
	itsFwdFxVol( rhs.itsFwdFxVol ),
	itsLocalFunctional(rhs.itsLocalFunctional)
{}

////////////////////////////////////////////////////
///	Class  : ARM_QModel1F_Fx
///	Routine: destructor
///	Returns: 
///	Action : destructor
////////////////////////////////////////////////////
ARM_QModel1F_Fx::~ARM_QModel1F_Fx()
{
}

////////////////////////////////////////////////////
///	Class   : ARM_QModel1F_Fx
///	Routines: toString
///	Returns :
///	Action  : stringify the object
////////////////////////////////////////////////////
string ARM_QModel1F_Fx::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "1F Q FX Model\n";
    os << indent << "-------------\n";
    os << ARM_PricingModel::toString(indent);

	if (!itsLocalFunctional.IsNull())
		os << itsLocalFunctional->toString();

    return os.str();
}

////////////////////////////////////////////////////
///	Class  : ARM_QModel1F_Fx
///	Routine: GetQ1FFXModelParams
///	Returns: ARM_ModelParamsQ1F_Fx*
///	Action : get the fxm model param
////////////////////////////////////////////////////

ARM_ModelParamsQ1F_Fx* ARM_QModel1F_Fx::GetQ1FFXModelParams()
{
	ARM_ModelParams* modelParam= GetModelParams();
	ARM_ModelParamsQ1F_Fx* fxModelParams = (ARM_ModelParamsQ1F_Fx*) modelParam;

#if defined(__GP_STRICT_VALIDATION)
	if( !dynamic_cast<const ARM_ModelParamsQ1F_Fx*>(modelParam) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": model params is not ARM_ModelParamsQ1F_Fx*" );
#endif
	return fxModelParams;
}

////////////////////////////////////////////////////
///	Class  : ARM_QModel1F_Fx
///	Routine: SetIRForeignModel
///	Returns: void
///	Action : sets the interest rate foreign model
////////////////////////////////////////////////////

void ARM_QModel1F_Fx::SetIRForeignModel( const ARM_PricingModelPtr& irModel )
{
	ValidateIRModel( irModel );
	itsIRForModel = irModel;
	GetQ1FFXModelParams()->SetForCurve( irModel->GetZeroCurve()  );
}


////////////////////////////////////////////////////
///	Class  : ARM_QModel1F_Fx
///	Routine: SetIRForeignModel
///	Returns: void
///	Action : sets the interest rate domestic model
////////////////////////////////////////////////////
void ARM_QModel1F_Fx::SetIRDomesticModel( const ARM_PricingModelPtr& irModel )
{
	ValidateIRModel( irModel );
	itsIRDomModel = irModel;
	GetQ1FFXModelParams()->SetDomCurve( irModel->GetZeroCurve() );
}


////////////////////////////////////////////////////
///	Class  : ARM_QModel1F_Fx
///	Routine: MarkovianDriftData
///	Returns: ARM_GP_VectorPtr
///	Action : Compute data for markovian drift calculation
////////////////////////////////////////////////////
void ARM_QModel1F_Fx::MarkovianDriftData( size_t timeIdx,
    double& qFx, double& fx0, double* qLinearCoef, double* fx0Deriv, double* dt ) const
{
    double evalTime = GetNumMethod()->GetTimeStep(timeIdx);
	const ARM_ModelParam& fxQModelParam	= GetModelParams()->GetModelParam(ARM_ModelParamType::QParameter);

    /// Get q(t) and S(0,t)
	qFx     = fxQModelParam.GetValue(evalTime);
	fx0     = (*GetPrecomputedFwds())[timeIdx];

    if(qLinearCoef)
    {
        /// Get q'(t)
        double qShift   = 0.1;
	    double qFxRd1	= (fxQModelParam.GetValue(evalTime+qShift) - qFx)/qShift;
	    double qFxLd1	= (qFx - fxQModelParam.GetValue(evalTime-qShift))/qShift;
        double qFxd1;
        if(fabs(qFxRd1-qFxLd1) > K_NEW_DOUBLE_TOL)
        {
            /// Not a simple linear interpolation or just around a breakPointTime
            const ARM_CurveModelParam& fxCurveQModelParam=fxQModelParam.ToCurveModelParam();
            ARM_Interpolator<double>* interpolator = fxCurveQModelParam.GetCurve()->GetInterpolator();
            if( dynamic_cast< ARM_StepUpLeftOpenCstExtrapol<double>* >(interpolator) ||
                dynamic_cast< ARM_StepUpRightOpenCstExtrapol<double>* >(interpolator) )
                qFxd1 = 0.0;
            else
                /// Linear interpolator => average left & right slopes
                qFxd1 = 0.5*(qFxLd1+qFxRd1);
        }
        else qFxd1 = qFxLd1;
        *qLinearCoef = (fabs(qFxd1) > K_NEW_DOUBLE_TOL && fabs(qFx) > K_NEW_DOUBLE_TOL ? qFxd1/qFx : 0.0);
    }

    if(fx0Deriv)
    {
        int nextTimeIdx;
        if(timeIdx+1 == GetNumMethod()->GetTimeSteps()->size())
            nextTimeIdx = timeIdx-1;
        else
            nextTimeIdx = timeIdx+1;

	    *fx0Deriv    = (*GetPrecomputedFwds())[nextTimeIdx];
	    *dt   	    = (GetNumMethod()->GetTimeStep(nextTimeIdx)-evalTime)/K_YEAR_LEN;
	    *fx0Deriv	= ((*fx0Deriv)-fx0)/(*dt);
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_QModel1F_Fx
///	Routine: ComputeForwardFxQ
///	Returns: double
///	Action : Compute between [t,T], the equivalent
///          fwd Fx Q using Piterbarg's method
///          assuming deterministic interest rates
///          then qFwdFx = Integ{s=t->T,w(s).q(s).ds})
////////////////////////////////////////////////////
double ARM_QModel1F_Fx::ComputeForwardFxQ(double t, double T, const ARM_GP_Vector& times, const ARM_GP_Vector& sigmas, const ARM_GP_Vector& qs, const ARM_GP_Vector& dqs, double& qDrift)
{
    /// Find Tn(t)-1 < t <= Tn(t) and Tn(T)-1< T <= Tn(T) for vol & Q
	int nt = lower_boundPosWithPrecision(times,t);
	int nT = lower_boundPosWithPrecision(times,T);

    size_t i,ii,nbTimes=times.size();
    double lastU=t/K_YEAR_LEN,U,dt;
    double vol,qU,q,dq,A,B,dA,var;
    double varTot=0.0,sumWQ=0.0;
    qDrift=0.0;
    for(i=nt;i<=nT;++i)
    {
		U=(i<nT ? times[i] : T)/K_YEAR_LEN;

        ii  = i<nbTimes ? i :nbTimes-1;

        vol = sigmas[ii];
		if(fabs(vol)<ARM_ModelParamsHW::VOL_LIMIT)
			vol = ARM_ModelParamsHW::VOL_LIMIT;
        qU  = times[ii]/K_YEAR_LEN;
        q   = qs[ii];
        dq  = dqs[i<nbTimes ? i :nbTimes];

        dq = 0.0;

        A   = dq;
        B   = dq*(lastU-qU)+q;

        dt  = U-lastU;
        dA  = A*dt;
        var = vol*vol*dt;

        sumWQ   +=2*var*(dA*var/3.0 + 0.5*(dA*varTot+B*var) + B*varTot);
        varTot  += var;
        qDrift  += var*(0.5*dA+B);

        lastU = U;
    }

    if(varTot < 1.0e-15)
    {
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": total variance null in equivalent fwd Q computation" );
    }
    else
        return sumWQ/(varTot*varTot);
}


////////////////////////////////////////////////////
///	Class  : ARM_QModel1F_Fx
///	Routine: IntegratedMarkovianDrift
///	Returns: ARM_GP_MatrixPtr
///	Action : add the MarkovianDrift to an existing drift
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_QModel1F_Fx::IntegratedMarkovianDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates, const ARM_GP_VectorPtr& driftCorrection ) const
{
	const ARM_ModelParamsQ1F* fxParams = static_cast< const ARM_ModelParamsQ1F*>(GetModelParams());

	double meanRevFx = fxParams->GetModelParam(ARM_ModelParamType::MeanReversion).GetValue(0);
    if(meanRevFx != 0.0)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Forex model MRS need to be set to 0" );

#if defined(__GP_STRICT_VALIDATION)
	if( itsIRDomModel == ARM_PricingModelPtr(NULL) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": itsIRDomModel == ARM_PricingModelPtr(NULL)" );
	if( itsIRForModel == ARM_PricingModelPtr(NULL) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": itsIRForModel == ARM_PricingModelPtr(NULL)" );
#endif	

	size_t i,nbStates = numMethodStates->cols();
	ARM_GP_MatrixPtr result	= ARM_GP_MatrixPtr( new ARM_GP_Matrix( numMethodStates->rows() , nbStates, 0.0 ) );
	size_t fxModelIdx = GetModelNb();

    const ARM_Curve& qFxCurve = * (fxParams->GetModelParam(ARM_ModelParamType::QParameter).ToCurveModelParam().GetCurve());
    if(qFxCurve == 1.0)
    {
        /// Compute integrated part of zero-coupon risk neutral drifts
        /// Curve dependent part disappears
        /// Other determinsitic parts will be computed at 3D spot FX calibration time
	    ARM_GP_VectorPtr integDriftDom = itsIRDomModel->IntegratedRiskNeutralDrift( timeIdx, numMethodStates );
	    ARM_GP_VectorPtr integDriftFor = itsIRForModel->IntegratedRiskNeutralDrift( timeIdx, numMethodStates );

#if defined(__GP_STRICT_VALIDATION)
	if( integDriftDom->size() != nbStates )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": integDriftDom->size() != nbStates" );
	if( integDriftFor->size() != nbStates )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": integDriftFor->size() != nbStates" );
#endif

	    for( i=0; i<nbStates; ++i )
	    {
            /// Only markovian terms appear here
		    (*result)(fxModelIdx,i) = (*integDriftDom)[i] - (*integDriftFor)[i];
	    }
    }
    else
    {
        double evalTime     = GetNumMethod()->GetTimeStep(timeIdx);
        double expiryTime   = GetNumMethod()->GetTimeStep(timeIdx+1);

/****
        // Initial method : not so accurate because it implies that fwd Fx S(t,t+dt)
        // is sampled by the lattice and the effective fwd Fx volatility
        // between t and t+dt  must depends of S(t,t+dt) => an ACP by state should be used !!

        /// First we use Fwd(t,t+dt) = Spot(t).ZcFor(t,t+dt)/ZcDom(t,t+dt)
        /// then Spot(t+dt)+m = Fwd(t+dt,t+dt)+m = (Fwd(t,t+dt)+m)*exp(drift(t,t+dt))*MG(volFwdFx(t,t+dt))
        /// with drift(t=t+dt) = probability change : fwd neutral dom t+dt -> dom risk neutral
        /// MG(volFwdFx(t,t+dt) = exponential martingale from t -> t+dt with fwd fx volatility
        /// Using Piterbarg's shifted LN forward FX approximation the forward forex vol and then
        /// its drift are both deterministic (domestic rates are assumed to be gaussian).
        /// So E[S(t+dt)] = (Fwd(t,t+dt)+m)*exp(drift(t,t+dt))-m

        /// Get S(0,t), qFx(t), S(0,t+dt), qFx(t+dt)
        double qFx,fx0;
        MarkovianDriftData(timeIdx,qFx,fx0);

        double qFxNext,fx0Next;
        MarkovianDriftData(timeIdx+1,qFxNext,fx0Next);

        /// Compute S(t,t+dt). Note that the calibrated deterministic drift
        /// at current slice is already added to states to correctly compute the fwd Fx
        ARM_PricingStatesPtr states(new ARM_PricingStates);
        states->SetModelStates(numMethodStates);

        ARM_GP_VectorPtr fwdFx = Forward(GetModelName(),evalTime,expiryTime,expiryTime,expiryTime,states);

       /// Compute drift correction (only spot vol part is used in fwd vol)
	    double domFxCorr	= (*GetCorrelMatrix())(DomModel,FxModel);
	    const ARM_ModelParamsHW1FStd* const domZcParams   = static_cast<const ARM_ModelParamsHW1FStd* const>( itsIRDomModel->GetModelParams() );
        const ARM_CurveModelParam& fxVolParam = static_cast< const ARM_CurveModelParam& >(fxParams->GetModelParam(ARM_ModelParamType::QVol));
        double absDrift = - domFxCorr * ARM_ModelParamsHW1F::HW1FEqFxZcCovariance(fxVolParam,domZcParams,evalTime,expiryTime,expiryTime);

        /// Compute the equivalent fwd Fx Q between both slices
        ARM_GP_Vector times,sigmas,qs,dqs;
        size_t nbTimes1 = fxParams->GetModelCurves(times,sigmas,qs,dqs);

        double qDrift,xDrift;
        double qFwdFx = ComputeForwardFxQ(evalTime,expiryTime,times,sigmas,qs,dqs,qDrift);

        double relDrift = exp(qFwdFx*absDrift);
        double xFx,fwdFxMean,xMean,fwdfFxShift;

        /// Link the conditional expectation of Spot(t+dt) using a Taylor expansion
        /// of the mapping function around the conditional expectation of X(t+dt)
        /// then compute the markovian drift E[X(t+dt)]-X(t)


//FILE* f=fopen("c:\\temp\\dumpQ1F_FX.txt","a");
//fprintf(f,"Idx=%3d (%8.1lf)\n",timeIdx,evalTime);
//fclose(f);
//ARM_TreeBase* tree = dynamic_cast< ARM_TreeBase* >(&*(GetNumMethod()));
//ARM_SliceNDBase* slice = dynamic_cast< ARM_SliceNDBase* >((*(tree->GetSlices()))[timeIdx]);
//ARM_IntVector treeMin(slice->GetMin()),treeMax(slice->GetMax());
//ARM_TreeIndex treeIdx(treeMin,treeMax);
//int nbPb=0;

        bool status;
	    for( i=0; i<nbStates; ++i )
	    {
            xFx = (*numMethodStates)(fxModelIdx,i);

            if(fabs(qFwdFx)>K_NEW_DOUBLE_TOL)
            {
                fwdfFxShift = (*fwdFx)[i]*(1/qFwdFx - 1);
                fwdFxMean   = ( (*fwdFx)[i] + fwdfFxShift )*relDrift - fwdfFxShift;
            }
            else
                fwdFxMean   = (*fwdFx)[i]*(1.0 + absDrift);

            xMean       = InverseMappingFunction( fwdFxMean, fx0Next, qFxNext, status );
            if(status)
            {
                /// Inversion succeded
		        xDrift = xMean - xFx;

                /// Deterministic drift part due to lognormality of
                /// the Q model between the adjacent slices must be removed.
                /// Last determinsitic drift part (already calibrated for timeIndex slice)
                /// is assumed to be identical at this stage for timeIndex+1 slice
                /// (then disappears automatically by difference)
                if(fabs(qFwdFx)>K_NEW_DOUBLE_TOL)
                    xDrift += 0.5*qDrift;
            }
            else
            {
                /// May occur only with very low probability

                xDrift = 0.0;

//++nbPb;
//f=fopen("c:\\temp\\dumpQ1F_FX.txt","a");
//while (treeIdx.GetPosition() != i) ++treeIdx;
//if(! treeIdx.IsOnHedge() && slice->GetArrowDebreuPrices(i) > 1.0e-5)
//    fprintf(f,"Node=%6d :(%2d/[%2d,%2d],%2d/[%2d,%2d],%2d/[%2d,%2d])\tfwdFx=%15.10lf\tAD=%15.10lf\n",i,treeIdx[0],treeMin[0],treeMax[0],
//            treeIdx[1],treeMin[1],treeMax[1],treeIdx[2],treeMax[2],treeMin[1],(*fwdFx)[i],slice->GetArrowDebreuPrices(i));
//fclose(f);

            }

            (*result)(fxModelIdx,i) = xDrift;

	    }


//f=fopen("c:\\temp\\dumpQ1F_FX.txt","a");
//fprintf(f,"TotPb=%6.2lf%%\n\n",nbPb*100.0/nbStates);
//fclose(f);

         /// Substract drift correction
        size_t j,nbFactors = driftCorrection->size();
	    for(j=0;j<nbFactors;++j)
            for(i=0;i<nbStates;++i)
                (*numMethodStates)(j,i) -= (*driftCorrection)[j];

****/

/****/
        /// The underlying process X(t) is diffused and markovian drift is integrated
        /// by freezing X(t)=X(ti) in integrals

        /// Compute S(ti,ti+1) for further use if necessary
        /// Note that the calibrated deterministic drift at current slice
        /// is already added to states to correctly compute the fwd Fx
        ARM_PricingStatesPtr states(new ARM_PricingStates);
        states->SetModelStates(numMethodStates);
        ARM_GP_VectorPtr fwdFx = Forward(GetModelName(),evalTime,expiryTime,expiryTime,expiryTime,states);


        /// Substract drift correction to compute correcly integrated RN
        size_t j,nbFactors = driftCorrection->size();
	    for(j=0;j<nbFactors;++j)
            for(i=0;i<nbStates;++i)
                (*numMethodStates)(j,i) -= (*driftCorrection)[j];

        /// Get integrated RN drifts  (IR curve independent part) =  
        /// ThetaDom(ti,ti+1) + BetaDom(ti,ti+1)Xdom(ti) - ThetaFor(ti,ti+1) - BetaFor(ti,ti+1)Xfor(ti)
        bool isDriftAdded = true;
	    ARM_GP_VectorPtr integDriftDom = itsIRDomModel->IntegratedRiskNeutralDrift( timeIdx, numMethodStates, isDriftAdded );
	    ARM_GP_VectorPtr integDriftFor = itsIRForModel->IntegratedRiskNeutralDrift( timeIdx, numMethodStates, isDriftAdded );

#if defined(__GP_STRICT_VALIDATION)
	    if( integDriftDom->size() != nbStates )
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": integDriftDom->size() != nbStates" );
	    if( integDriftFor->size() != nbStates )
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": integDriftFor->size() != nbStates" );
#endif


        /// Quanto correction induced by foreign short rate integral
		ARM_GP_Matrix correlMatrix = GetCorrelMatrix().Interpolate(expiryTime);
	    double forFxCorr = correlMatrix(ForModel,FxModel);
	    const ARM_ModelParamsHW1FStd* const forZcParams   = static_cast<const ARM_ModelParamsHW1FStd* const>( itsIRForModel->GetModelParams() );
	    const ARM_ModelParamsQ1F* const fxParams   = static_cast<const ARM_ModelParamsQ1F* const>( GetModelParams() );
        const ARM_CurveModelParam& fxVolParam = static_cast< const ARM_CurveModelParam& >(fxParams->GetModelParam(ARM_ModelParamType::QVol));
        double quantoDrift = 0.0;
        if(forFxCorr!=0.0)
            quantoDrift = forFxCorr * ARM_ModelParamsQ1F::Q1FStateZcCovariance(fxParams,forZcParams,evalTime,expiryTime,expiryTime,expiryTime);
        double residualXDrift = -0.5*ARM_ModelParamsQ1F::Q1FStateCovariance(fxParams,fxParams,evalTime,expiryTime,expiryTime,false);

       /// Compute drift correction and effective fwd Fx Q (only spot vol part is used in fwd vol)
	    double domFxCorr = correlMatrix(DomModel,FxModel);
	    const ARM_ModelParamsHW1FStd* const domZcParams   = static_cast<const ARM_ModelParamsHW1FStd* const>( itsIRDomModel->GetModelParams() );
        double fwdToCashDrift;
        double qDrift,qFwdFx,mFwdFx;
        bool isFwdToCashDrift = (domFxCorr != 0.0);
        if(isFwdToCashDrift)
        {
            fwdToCashDrift = - domFxCorr * ARM_ModelParamsQ1F::Q1FStateZcCovariance(fxParams,domZcParams,evalTime,expiryTime,expiryTime,expiryTime);

            /// Compute the equivalent fwd Fx Q between both slices
            ARM_GP_Vector times,sigmas,qs,dqs;
            size_t nbTimes1 = fxParams->GetModelCurves(times,sigmas,qs,dqs);

            qFwdFx = ComputeForwardFxQ(evalTime,expiryTime,times,sigmas,qs,dqs,qDrift);
        }

        /// Compute Integ{ti->ti+1, q^2(t).volXS^2(t).dt}
        double qVar = ARM_ModelParamsQ1F::Q1FStateCovariance(fxParams,fxParams,evalTime,expiryTime,expiryTime);


        /// Get S(0,ti), q(ti) & q(ti+1)
        double qFx,fx0,mFx;
        MarkovianDriftData(timeIdx,qFx,fx0);
        bool isQFxTiny;
        if((isQFxTiny=(fabs(qFx)<QFX_LIMIT)))
            qFx = (qFx > 0 ? QFX_LIMIT : -QFX_LIMIT);
        mFx = fx0*(1/qFx-1);
        double qFxNext,fx0Next,mFxNext;
        MarkovianDriftData(timeIdx+1,qFxNext,fx0Next);
        if(fabs(qFxNext)<QFX_LIMIT)
            qFxNext = (qFxNext > 0 ? QFX_LIMIT : -QFX_LIMIT);
        mFxNext = fx0Next*(1/qFxNext-1);

        double qFxRatio = log(fabs(qFxNext/qFx));

        /// Compute a maximum value of (1+mFx/fx) due to limitation
        /// on LN spot FX vol ( = q.volXS.(1+mFx/fx) )
        if(qVar < K_DOUBLE_TOL)
            qVar = K_DOUBLE_TOL;
        double maxDriftTerm = LNVOLFX_LIMIT * sqrt((expiryTime-evalTime) / (qVar*K_YEAR_LEN));

        double fx,xFx,yFx,zFx,stochCoef1,stochCoef2,driftTerm;
        double xDrift,fxNext,xFxNext;
        bool inversionSuccess;
	    for( i=0; i<nbStates; ++i )
	    {
            xFx = (*numMethodStates)(fxModelIdx,i) + (*driftCorrection)[fxModelIdx];
		    fx  = MappingFunction(xFx,fx0,qFx);

            if(isQFxTiny)
            {
                stochCoef1  = fx/fx0;
                stochCoef2  = 0.0;
            }
            else
            {
		        yFx         = 1.0/( ((fx/fx0)-1.0)*qFx+1.0 ); /// or exp(-qFx*xFx)
                zFx         = (1.0-yFx)/qFx;
                stochCoef1  = yFx + zFx;
                stochCoef2  = zFx - xFx;
            }

            driftTerm = 1+mFx/fx;
            if(driftTerm < -maxDriftTerm)
                driftTerm = - maxDriftTerm;
            else
            {
                if(driftTerm > maxDriftTerm)
                    driftTerm = maxDriftTerm;
            }

            /// -0.5*Integ(ti->ti+1, q(t).volXS(t)^2) not added because
            /// will be computed at 3D spot FX calibration time
		    xDrift = stochCoef2 * qFxRatio
                + stochCoef1 * ( (*integDriftDom)[i] - (*integDriftFor)[i] - driftTerm*quantoDrift );

            if(evalTime/K_YEAR_LEN >= 10.0)
            {
                /// Drift may diverge, switch to do better way using S(ti,ti+1);
                if(isFwdToCashDrift)
                {
                    if(fabs(qFwdFx)>K_NEW_DOUBLE_TOL)
                    {
                        mFwdFx      = (*fwdFx)[i]*(1/qFwdFx - 1);
                        fxNext   = ( (*fwdFx)[i] + mFwdFx )*exp(qFwdFx*driftTerm*fwdToCashDrift) - mFwdFx;
                    }
                    else
                        fxNext   = (*fwdFx)[i]*(1 + qFwdFx*driftTerm*fwdToCashDrift);
                }
                else
                    fxNext = (*fwdFx)[i];

                xFxNext = InverseMappingFunction(fxNext,fx0Next,qFxNext,inversionSuccess);
                if(inversionSuccess)
                    /// Inversion succeded
		            xDrift = xFxNext - xFx - residualXDrift;
                else
                {
                    /// May occur only with very low probability
                    xDrift = 0.0;

                }
            }

            (*result)(fxModelIdx,i) = xDrift;
        }

/****/

    }

    return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_QModel1F_Fx
///	Routine: MarkovianDrift
///	Returns: ARM_GP_MatrixPtr
///	Action : add the MarkovianDrift to an existing drift
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_QModel1F_Fx::MarkovianDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates ) const
{
	double meanRevFx = GetModelParams()->GetModelParam(ARM_ModelParamType::MeanReversion).GetValue(0);

#if defined(__GP_STRICT_VALIDATION)
	if( itsIRDomModel == ARM_PricingModelPtr(NULL) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": itsIRDomModel == ARM_PricingModelPtr(NULL)" );
	if( itsIRForModel == ARM_PricingModelPtr(NULL) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": itsIRForModel == ARM_PricingModelPtr(NULL)" );
    if(meanRevFx != 0.0)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Forex model MRS need to be set to 0" );
#endif	

	size_t nbStates		    = numMethodStates->cols();
	ARM_GP_MatrixPtr result	= ARM_GP_MatrixPtr( new ARM_GP_Matrix( numMethodStates->rows() , nbStates, 0.0 ) );
	size_t i;

	size_t fxModelIdx = GetModelNb();

    /// Collect markovian drift datas
    double qFx,qLinearCoef,fx0,fx0Deriv,dt;
    MarkovianDriftData(timeIdx,qFx,fx0,&qLinearCoef,&fx0Deriv,&dt);

    /// Get instantaneous RN drifts (drift corrections are aleardy added in numMethodStates)
	ARM_GP_VectorPtr instDriftDom = itsIRDomModel->RiskNeutralDrift( timeIdx, numMethodStates );
	ARM_GP_VectorPtr instDriftFor = itsIRForModel->RiskNeutralDrift( timeIdx, numMethodStates );

#if defined(__GP_STRICT_VALIDATION)
	if( instDriftDom->size() != nbStates )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": instDriftDom->size() != nbStates" );
	if( instDriftFor->size() != nbStates )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": instDriftFor->size() != nbStates" );
#endif

	double fx,stochCoef;
    double eqMinusXt=1.0;
	for( i=0; i<nbStates; ++i )
	{
		fx = MappingFunction( (*numMethodStates)(fxModelIdx,i), fx0, qFx );
        if(fabs(qFx) <= K_NEW_DOUBLE_TOL)
            stochCoef   = fx/fx0;
        else
        {
		    eqMinusXt   = 1.0/( ((fx/fx0)-1.0)*qFx+1.0 );
            stochCoef   = eqMinusXt+(1.0-eqMinusXt)/qFx;
        }

        /// Only markovian terms appear here
		(*result)(fxModelIdx,i) =
            ( ( (*instDriftDom)[i]-(*instDriftFor)[i] -fx0Deriv/fx0 ) * stochCoef
              + qLinearCoef/qFx * ( (1-eqMinusXt)/qFx - (*numMethodStates)(fxModelIdx,i) )
			)*dt;
	}

	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_QModel1F_Fx
///	Routine: MarkovianDriftPDE
///	Returns: ARM_GP_MatrixPtr
///	Action : Compute the markbkovian drift
////////////////////////////////////////////////////
void ARM_QModel1F_Fx::MarkovianDriftPDE( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates, ARM_GP_VectorPtr result ) const
{
	double meanRevFx = GetModelParams()->GetModelParam(ARM_ModelParamType::MeanReversion).GetValue(0);

#if defined(__GP_STRICT_VALIDATION)
	if( itsIRDomModel == ARM_PricingModelPtr(NULL) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": itsIRDomModel == ARM_PricingModelPtr(NULL)" );
	if( itsIRForModel == ARM_PricingModelPtr(NULL) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": itsIRForModel == ARM_PricingModelPtr(NULL)" );
    if(meanRevFx != 0.0)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Forex model MRS need to be set to 0" );
#endif	

	size_t nbStates		    = numMethodStates->cols();
	size_t i;

	size_t fxModelIdx = GetModelNb();

    /// Collect markovian drift datas
    double qFx,qLinearCoef,fx0,fx0Deriv,dt;
    MarkovianDriftData(timeIdx,qFx,fx0,&qLinearCoef,&fx0Deriv,&dt);

	ARM_GP_VectorPtr instDriftDom;
	ARM_GP_VectorPtr instDriftFor;

	double t,T, DomForCorr, domTheta, forTheta;
	
	if(GetNumeraire()->GetType() == ARM_Numeraire::Cash)
	{
		/// Get instantaneous RN drifts (drift corrections are aleardy added in numMethodStates)
		instDriftDom = itsIRDomModel->RiskNeutralDrift( timeIdx, numMethodStates, true );
		instDriftFor = itsIRForModel->RiskNeutralDrift( timeIdx, numMethodStates, true );
	}
	else if(GetNumeraire()->GetType() == ARM_Numeraire::TerminalZc)
	{
		 /// Get instantaneous FN drifts (drift corrections are aleardy added in numMethodStates)
		instDriftDom = itsIRDomModel->RiskNeutralDrift( timeIdx, numMethodStates, true );//pour le moment
		instDriftFor = itsIRForModel->RiskNeutralDrift( timeIdx, numMethodStates, true );
	    const ARM_ModelParamsHW1FStd* const modelParamsDom   = static_cast<const ARM_ModelParamsHW1FStd* const>( itsIRDomModel->GetModelParams() );
	    const ARM_ModelParamsHW1FStd* const modelParamsFor   = static_cast<const ARM_ModelParamsHW1FStd* const>( itsIRForModel->GetModelParams() );
 
        t = GetNumMethod()->GetTimeStep(timeIdx);
		ARM_GP_Matrix correlMatrix = GetCorrelMatrix().Interpolate(t);
		T=GetNumeraire()->GetMaturity();
		DomForCorr = correlMatrix(DomModel,ForModel);
		domTheta=ARM_ModelParamsHW1F::HW1FStateZcCovariance( modelParamsDom, modelParamsDom, 0.0, t, t, T );
		forTheta=DomForCorr*ARM_ModelParamsHW1F::HW1FStateZcCovariance( modelParamsFor, modelParamsDom, 0.0, t, t, T );

		*instDriftDom= domTheta + *instDriftDom;
		*instDriftFor= forTheta + *instDriftFor;

	}
	else
	{
	    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": not Cash or TerminalZC numeraire!");
	}

#if defined(__GP_STRICT_VALIDATION)
	if( instDriftDom->size() != nbStates )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": instDriftDom->size() != nbStates" );
	if( instDriftFor->size() != nbStates )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": instDriftFor->size() != nbStates" );
#endif

	double fx,stochCoef;
    double eqMinusXt=1.0;
	for( i=0; i<nbStates; ++i )
	{
		fx = MappingFunction( (*numMethodStates)(fxModelIdx,i), fx0, qFx );
        if(fabs(qFx) <= K_NEW_DOUBLE_TOL)
            stochCoef   = fx/fx0;
        else
        {
		    eqMinusXt   = 1.0/( ((fx/fx0)-1.0)*qFx+1.0 );
            stochCoef   = eqMinusXt+(1.0-eqMinusXt)/qFx;
        }

        /// Only markovian terms appear here
		(*result)(i) = 
			  ( (*instDriftDom)[i]-(*instDriftFor)[i]  ) * stochCoef
              + qLinearCoef/qFx * ( (1-eqMinusXt)/qFx - (*numMethodStates)(fxModelIdx,i) );
	}

}

////////////////////////////////////////////////////
///	Class   : ARM_QModel1F_EqFxBase
///	Routines: LocalDiscounts
///	Returns : void
///	Action  : Computes the LocalDiscounts
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QModel1F_Fx::LocalDiscounts(
	size_t timeIdx, 
	double dt, 
	const ARM_PricingStatesPtr& states) const
{
	size_t stateSize= states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
	const ARM_GP_Vector* const timeSteps = GetNumMethod()->GetTimeSteps();
	double startTime = (*timeSteps)[timeIdx];

	ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();
	double zct	= ZcCurve->DiscountPrice(startTime/K_YEAR_LEN);
	double zcT	= ZcCurve->DiscountPrice((startTime+dt)/K_YEAR_LEN);

	return ARM_VectorPtr( new ARM_GP_Vector( stateSize, zcT/zct ) );
}


////////////////////////////////////////////////////
///	Class  : ARM_QModel1F_EqFxBase
///	Routine: LocalDrifts
///	Returns: void
///	Action : local drifts
////////////////////////////////////////////////////
void ARM_QModel1F_Fx::IntegratedLocalDrifts(
	const ARM_GP_Vector& timeSteps,
	ARM_GP_MatrixPtr& relativeDrifts,
	ARM_GP_MatrixPtr& absoluteDrifts) const
{
	size_t nbSteps	= timeSteps.size();
    double step		= timeSteps[0],nextStep;
	relativeDrifts	= ARM_GP_MatrixPtr( new ARM_GP_Matrix( nbSteps-1, 1, 0.0 ) );
	absoluteDrifts	= ARM_GP_MatrixPtr(NULL);

	for(size_t i=0;i<nbSteps-1;++i)
	{
		nextStep=timeSteps[i+1];
		/// [i] => local variance from ti->ti+1
		(*relativeDrifts)(i,0) = ((const ARM_ModelParamsQ1F* const) GetModelParams())->StateLocalDrift(step,nextStep);
		step=nextStep;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_QModel1F_Fx
///	Routine: ModelStateLocalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute local and global variances
///          (from 0) of the state variable between
///          each time step
////////////////////////////////////////////////////
void ARM_QModel1F_Fx::ModelStateLocalVariances(
	const ARM_GP_Vector& timeSteps,
	ARM_MatrixVector& localVariances ) const
{
	size_t nbSteps	= timeSteps.size();
	size_t modelNb	= GetModelNb();
    double step		= timeSteps[0],
		   nextStep;
	size_t offsetIndex	= (nbSteps-1)*modelNb;
	localVariances.resize((nbSteps-1)*(modelNb+1));
	size_t i;

	for(i=0;i<nbSteps-1;++i)
	{
		nextStep=timeSteps[i+1];
		
		/// [i] => local variance from ti->ti+1
		localVariances[offsetIndex+i] = new ARM_GP_TriangularMatrix(1,((const ARM_ModelParamsQ1F* const) GetModelParams())->StateLocalVariance(step,nextStep,nextStep));
	
		step=nextStep;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_QModel1F_Fx
///	Routine: FirstPricingStates
///	Returns: ARM_PricingStatesPtr
///	Action  : create the first pricing state
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_QModel1F_Fx::FirstPricingStates( size_t bucketSize ) const
{
	/// ARM_PricingStates(nbStates = bucketSize, nbModelStates = 1F , nbPayoffs = 0)
	ARM_PricingStatesPtr initStates( new ARM_PricingStates(bucketSize,1,0,1) );
    double initValue = 0.0;
	if( itsIRForModel != ARM_PricingModelPtr(NULL) && itsIRDomModel != ARM_PricingModelPtr(NULL) &&
        itsIsForexDiffusion )
        initValue = ComputeFwdAtTime(0.0);

    for(size_t i=0;i<bucketSize;++i)
        initStates->SetModelState(i,0,initValue);

    return initStates;
}


////////////////////////////////////////////////////
///	Class  : ARM_QModel1F_Fx
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : move model states for timeIndex to timeIndex+1
////////////////////////////////////////////////////
void ARM_QModel1F_Fx::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{

#if defined(__GP_STRICT_VALIDATION)
	if(timeIndex<0)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index is negative!");
	if( timeIndex >= GetNumMethod()->GetTimeSteps()->size()-1 )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index bigger than max size!");
#endif

	if( itsIRForModel != ARM_PricingModelPtr(NULL) && itsIRDomModel != ARM_PricingModelPtr(NULL) )
	{
		const ARM_NumMethodPtr numMethod= GetNumMethod();
        size_t nbSteps = numMethod->GetTimeSteps()->size();
		double evalTime	= numMethod->GetTimeStep(timeIndex);
		double nextTime	= numMethod->GetTimeStep(timeIndex+1);
		size_t i,nbStates = states->size();
	    size_t modelNb	= GetModelNb();

        bool isLnFwdFx;
        bool isDriftAdded;

		ARM_GP_MatrixPtr numMethodStates    = states->GetNumMethodStates();
		ARM_GP_MatrixPtr modelStates        = states->GetModelStates();

		if( GetIntegratedVersion() )
		{
	        const ARM_ModelParamsQ1F* modelParams = static_cast< const ARM_ModelParamsQ1F*>(GetModelParams());
            const ARM_Curve& qFxCurve = * (modelParams->GetModelParam(ARM_ModelParamType::QParameter).ToCurveModelParam().GetCurve());
            isLnFwdFx = qFxCurve == 1.0;

            ARM_GP_VectorPtr domDf(NULL),forDf(NULL);
            ARM_GP_VectorPtr integDriftDom(NULL),integDriftFor(NULL);
            double qti,qtk,ti,S0ti,mti,Bdom0ti,Bfor0ti,tk,S0tk,mtk,Bdom0tk,Bfor0tk,dt;
            double domFwdRate,forFwdRate;
            double qDrift;

            if(isLnFwdFx)
            {
                /// Pure lognormal forward FX S(t,ti+1) is diffused from ti to ti+1 and
                /// then S(ti+1)=S(ti+1,ti+1)

			    domDf = itsIRDomModel->DiscountFactor( itsIRDomModel->GetModelName(), evalTime, nextTime, states );
			    forDf = itsIRForModel->DiscountFactor( itsIRForModel->GetModelName(), evalTime, nextTime, states );
            }
            else
            {
                ti = GetNumMethod()->GetTimeStep(timeIndex);

                qti=qFxCurve.Interpolate(ti);

                S0ti = ComputeFwdAtTime(ti);
                mti=S0ti*(1.0/qti-1.0);

                domDf = itsIRDomModel->DiscountFactor( itsIRDomModel->GetModelName(), 0.0, ti, states );
                Bdom0ti = (*domDf)[0];
                forDf = itsIRForModel->DiscountFactor( itsIRForModel->GetModelName(), 0.0, ti, states );
                Bfor0ti = (*forDf)[0];

                tk = GetNumMethod()->GetTimeStep(timeIndex+1);

                qtk=qFxCurve.Interpolate(tk);

                S0tk = ComputeFwdAtTime(tk);
                mtk=S0tk*(1.0/qtk-1.0);

                dt = (tk-ti)/K_YEAR_LEN;

                domDf = itsIRDomModel->DiscountFactor( itsIRDomModel->GetModelName(), 0.0, tk, states );
                Bdom0tk = (*domDf)[0];
                forDf = itsIRForModel->DiscountFactor( itsIRForModel->GetModelName(), 0.0, tk, states );
                Bfor0tk = (*forDf)[0];

                domFwdRate = log(Bdom0ti/Bdom0tk);
                forFwdRate = log(Bfor0ti/Bfor0tk);

                isDriftAdded = true;

	            integDriftDom = itsIRDomModel->IntegratedRiskNeutralDrift( timeIndex, modelStates, isDriftAdded );
	            integDriftFor = itsIRForModel->IntegratedRiskNeutralDrift( timeIndex, modelStates, isDriftAdded );
            }

            /// Compute lognormal drift of diffused variable (S(t,ti+1) or S(t))
            const ARM_GP_Matrix& localVar = *(GetModelStateLocalVars()[timeIndex]);
            double lnDrift = -0.5*localVar(modelNb,modelNb);

            double spotFx,spotFx1,spotFx2,markovDrift;

			for(i=0;i<nbStates;++i )
            {
                if(isLnFwdFx)
                {
                    /// Compute S(ti,ti+1)
			        double fwdFx = (*modelStates)(modelNb,i)*(*forDf)[i]/(*domDf)[i];

                    /// Diffuse fwd Fx to ti+1
                    (*modelStates)(modelNb,i) = MappingFunction((*numMethodStates)(modelNb,i) + lnDrift,fwdFx,1.0);

                }
                else
                {
                    /// Get S(ti)
                    spotFx = (*modelStates)(modelNb,i);

                    /// Becareful : if S(t)+m(t) is close to 0,
                    /// dS = (rd(t)-rf(t))S(t)dt + q(t).vol(t)(S(t)+m(t))dW
                    /// has no more volatility !!
                    bool isZeroSpot = fabs(spotFx)<0.01;
                    bool isNoVol    = !isZeroSpot && fabs(1+mti/spotFx)<0.005;

                    /// Predictor : compute S1(ti+1) assuming S(t)=S(ti)
                    if(isNoVol)
                    {
                        /// No more volatility on shifted log spot Fx => dS/S = (rdom-rfor).dt
                        markovDrift = ( domFwdRate - forFwdRate + (*integDriftDom)[i] - (*integDriftFor)[i] );
                        spotFx1 = spotFx * exp(markovDrift);
                    }
                    else
                    {
						qDrift = (mtk-mti)/(spotFx+mti);
						markovDrift = (isZeroSpot ? 0.0 : ( domFwdRate - forFwdRate + (*integDriftDom)[i] - (*integDriftFor)[i] )/(1.0 + mti/spotFx)) + qDrift;
						spotFx1  = (spotFx+mti) * exp( markovDrift + (*numMethodStates)(modelNb,i) + lnDrift ) - mtk;
                    }

                    isZeroSpot = fabs(spotFx1)<0.01;
                    isNoVol    = !isZeroSpot && fabs(1+mtk/spotFx1)<0.005;

                    /// Corrector : compute S2(ti+1) assuming S(t)=S1(ti+1)
                    if(isNoVol)
                    {
                        markovDrift = (domFwdRate - forFwdRate + (*integDriftDom)[i]  - (*integDriftFor)[i] );
                        spotFx2 = spotFx * exp(markovDrift);
                    }
                    else
                    {
						qDrift = (mtk-mti)/(spotFx1+mtk);
						markovDrift = (isZeroSpot ? 0.0 : ( domFwdRate - forFwdRate + (*integDriftDom)[i] - (*integDriftFor)[i] )/(1.0 + mtk/spotFx1)) + qDrift;
						spotFx2  = (spotFx+mti) * exp( markovDrift + (*numMethodStates)(modelNb,i) + lnDrift ) - mtk;
                    }

                    (*modelStates)(modelNb,i) = 0.5*(spotFx1+spotFx2);

                }

            }
		}
		else if(itsIsForexDiffusion)
		{
			/// Euler discretisation scheme on spot Fx
			/// spot(t+dt) = spot(t)*(1 + (rd(t)-rf(t))*dt) + volNor(t)*dWt
            double t		= GetNumMethod()->GetTimeStep(timeIndex);
            double nextt	= GetNumMethod()->GetTimeStep(timeIndex+1);
            double dt		= (nextt-t)/K_YEAR_LEN;

			double qVol	= GetModelParams()->GetModelParam(ARM_ModelParamType::QVol).GetValue(t);
			double q	= GetModelParams()->GetModelParam(ARM_ModelParamType::QParameter).GetValue(t);

            double fx0t = ComputeFwdAtTime(t);

            isDriftAdded = true;

			ARM_GP_VectorPtr rd = itsIRDomModel->RiskNeutralDrift(timeIndex,modelStates,isDriftAdded);
            ARM_GP_VectorPtr df = itsIRDomModel->DiscountFactor(itsIRDomModel->GetModelName(),0.0,t,states);
            double zc0t = (*df)[0];
            df = itsIRDomModel->DiscountFactor(itsIRDomModel->GetModelName(),0.0,nextt,states);
            double fd0t = log(zc0t / (*df)[0]);

			ARM_GP_VectorPtr rf = itsIRForModel->RiskNeutralDrift(timeIndex,modelStates,isDriftAdded);
            df = itsIRForModel->DiscountFactor(itsIRForModel->GetModelName(),0.0,t,states);
            zc0t = (*df)[0];
            df = itsIRForModel->DiscountFactor(itsIRForModel->GetModelName(),0.0,nextt,states);
            double ff0t = log(zc0t / (*df)[0]);

            double fx,volNorFx;
			double volRelShift = q * qVol;
			double volAbsShift = (qVol - volRelShift) * fx0t;
			for(i=0;i<nbStates;++i)
			{
				fx			= (*modelStates)(modelNb,i);
				volNorFx	= volRelShift * fx + volAbsShift;
				(*modelStates)(modelNb,i) = fx * ( 1.0 + fd0t - ff0t + ((*rd)[i] - (*rf)[i])*dt )
											+ volNorFx * (*numMethodStates)(modelNb,i);
			}
		}
		else
		{
			/// add the drift correction due to stochastic model
			ARM_GP_MatrixPtr markovianDrift = MarkovianDrift(timeIndex,states->GetNumMethodStates());

			const ARM_MatrixVector& localVar	= GetModelStateLocalVars();
			const ARM_MatrixVector& localStdDev = GetModelStateLocalStdDevs();

			double currentState,stdDev;
		
			for( size_t i=0;i<nbStates; ++i )
			{
				currentState = (*markovianDrift)(modelNb,i) + states->GetModelState(i,modelNb);
				stdDev			= localStdDev[ OffsetTimeIndexWithModelNb(timeIndex,modelNb) ]->Elt(0,0);
				currentState  += stdDev*states->GetNumMethodState(i,modelNb);
				states->SetModelState(i,modelNb,currentState);
			}
		}
	}
	else 
	{
		double nextTime = GetNumMethod()->GetTimeStep(timeIndex+1);
		size_t factorsNb= FactorCount();
		size_t statesNb = states->size();
		double currentState;
		size_t modelNb	= GetModelNb();
		
		for( size_t i=0;i<statesNb; ++i )
		{
			for( size_t j=0;  j<factorsNb; ++j )
			{
				currentState   = 0.0;
				for( size_t k =0; k<=j; ++k )
				{
					double gaussian = states->GetNumMethodState(i,modelNb+k);
					currentState  += gaussian;
				}
				states->SetModelState(i,j+modelNb,currentState);
			}
		}
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_QModel1F_Fx
///	Routine: SetNumMethod
///	Returns: void 
///	Action : Set the num method
////////////////////////////////////////////////////
void ARM_QModel1F_Fx::SetNumMethod(const ARM_NumMethodPtr& numMethodPtr)
{
	if( !GetFromMultiFactor() )
	{
		const ARM_TreeBase* treeBase = dynamic_cast<const ARM_TreeBase*>(&*numMethodPtr);
		if(treeBase)
		{
			if( dynamic_cast<const ARM_Tree1D*>(treeBase) )
			{
				if( !dynamic_cast<const ARM_MarkovianDriftSampler1D*>(treeBase->GetSampler()) &&
					!dynamic_cast<const ARM_MeanRevertingSampler1D*>(treeBase->GetSampler()))
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
					" Q model requires Markovian drift or Mean revering samplers !" );
			}
			else
			{
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " tree is not 1D compatible" );
			}
		}
	}
	ARM_PricingModel::SetNumMethod( numMethodPtr );
}

////////////////////////////////////////////////////
///	Class   : ARM_QModel1F_Fx
///	Routine : SupportAnalyticMarginal
///	Returns : bool
///	Action  : cannot support analytical marginal as the Markovian drift
///				imposes to compute using Euler scheme the Markovian drift correction
////////////////////////////////////////////////////
bool ARM_QModel1F_Fx::SupportAnalyticMarginal() const
{
#if defined(__GP_STRICT_VALIDATION)
	if( itsIRDomModel == ARM_PricingModelPtr(NULL ) &&  itsIRForModel != ARM_PricingModelPtr(NULL ) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": the two ir models are either deterministic or stochastic" );
	if( itsIRDomModel != ARM_PricingModelPtr(NULL ) &&  itsIRForModel == ARM_PricingModelPtr(NULL ) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": the two ir models are either deterministic or stochastic" );
#endif

	if( itsIRDomModel == ARM_PricingModelPtr(NULL ) && itsIRForModel == ARM_PricingModelPtr(NULL ))
		return true;
	else
		return GetIntegratedVersion();
}


////////////////////////////////////////////////////
///	Class   : ARM_QModel1F_Fx
///	Routine : ARM_PricingStatesPtr
///	Returns : Init
///	Action  : initialise the model
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_QModel1F_Fx::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{
    int nbEvents=timeInfos.size();
    bool isSpotUse = nbEvents == 0 || (nbEvents==1 && timeInfos[0]->GetEventTime() <= K_NEW_DOUBLE_TOL);

    if(!isSpotUse)
    {
        // Numerical method and numeraire are needed to go on
        ARM_NumMethodPtr numMethod=GetNumMethod();
		if( numMethod == ARM_NumMethodPtr(NULL))
		{
			CC_Ostringstream os;
			os << ARM_USERNAME << ": numerical method not set in the Q mode!";
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
		}

        ARM_NumerairePtr numeraire=GetNumeraire();
        if( numeraire == ARM_NumerairePtr(NULL) )
		{
			CC_Ostringstream os;
			os << ARM_USERNAME << ": numeraire not set in the Qmode!";
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
            
		}

        if( numeraire->GetType() != ARM_Numeraire::Cash )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
            ": only Cash numeraire supported by Q model at the moment!");
		
		/// creates the model schedule (smart pointor for exception safety!)
		ARM_DiscretisationScheme& discretisationScheme = ARM_EventTime();
		CC_NS(std,auto_ptr)<ARM_GP_Vector> ptimeSteps( discretisationScheme.ModelTimesFromTimesInfo(timeInfos, *this) );

        //// Initialise the numeraire
		numeraire->Init(*(GetDiscountFunctor()),payModelName,timeInfos);

		/// Set the basic schedule in the numerical method and...
		numMethod->SetTimeSteps(*ptimeSteps);

		double firstInductTime = timeInfos[0]->GetEventTime();

		/// ...initialise it
		return numMethod->Init(*this,firstInductTime);
	}
	else
	{
        // Compute a single model states set to (0.0,...,0.0)
        ARM_PricingStatesPtr initStates(new ARM_PricingStates(1,1,0));
		initStates->EnqueuePricingStatesContext( ARM_PricingStatesContextPtr( new ARM_Q1FPricingStatesContext() ) );
        for(size_t i=0;i<1;++i)
            initStates->SetModelState(0,i,0.0);
        return initStates;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_QModel1F_Fx
///	Routine: PostInit
///	Returns: nothing
///	Action : Function to init after the numerical method
///			 has been initialised.. non const
////////////////////////////////////////////////////
void ARM_QModel1F_Fx::PostInit()
{
    ARM_GP_Vector* timeSteps = GetNumMethod()->GetTimeSteps();
	size_t timeStepSize = timeSteps->size();
	itsPrecomputedFwds = ARM_GP_VectorPtr(new ARM_GP_Vector(timeStepSize,0.0));
	
	for( size_t i=0; i<timeStepSize; ++i )
		(*itsPrecomputedFwds)[i] = ComputeFwdAtTime( (*timeSteps)[i] );
}


////////////////////////////////////////////////////
///	Class   : ARM_QModel1F_Fx
///	Routine	: GetStochasticModel (static function)
///	Returns : ARM_QModel1F*
///	Action  : the stochastic model for the interest rate model
////////////////////////////////////////////////////

ARM_QModel1F* ARM_QModel1F_Fx::GetStochasticModel( const ARM_PricingModelPtr& model )
{
	if( model != ARM_PricingModelPtr(NULL) )
	{
		if( dynamic_cast<ARM_ForwardMargin*>( &*model ) )
			return dynamic_cast<ARM_QModel1F*>( model->GetRefModel() );
		else
			return dynamic_cast<ARM_QModel1F*>( &*model );
	}
	else
		return NULL;
}

////////////////////////////////////////////////////
///	Class  : ARM_QModel1F_Fx
///	Routine: Forward
///	Returns: a vector of forward (t,T)
///	Action : computes the forward of fx model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QModel1F_Fx::Forward(
		const string& modelName, 
        double evalTime,
	    double expiryTime,
	    double settlementTime,
	    double payTime,
        const ARM_PricingStatesPtr& states) const
{
	//Get itsLambda which enables to choose the PDE and the reconstruction formula(centred, non centred, Qmapping, SpotFX)
	double lambda=0.0;
	const ARM_NumMethodPtr& numMethod = GetNumMethod();
	ARM_PDEMethod* pdemethod = dynamic_cast<ARM_PDEMethod*>( &*numMethod );
	if( pdemethod )
	{
		ARM_PDENumericalScheme* scheme=pdemethod->GetPDENumericalScheme();
		int NumSchemeId=scheme->getNumericalSchemeType();
		scheme->getNumericalSchemeInstanceById(NumSchemeId,-1,-1,-1,-1.0,-1.0,-1.0,-1,0.0,-1,ARM_GP_Matrix(0));
		ARM_PDE3FCraigSneydNumericalScheme* CS3Fscheme = dynamic_cast<ARM_PDE3FCraigSneydNumericalScheme*>( &*scheme );
		lambda=CS3Fscheme->GetLambda();
	}
    
	if( !IsFwdModels() || GetNumMethod() == ARM_NumMethodPtr(NULL)
        || states == ARM_PricingStatesPtr(NULL) || (evalTime<=K_NEW_DOUBLE_TOL))
	{
		bool isFwdAndInteg = GetNumMethod() != ARM_NumMethodPtr(NULL) && 
		GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_FWDLOOKING && 
		itsIntegratedVersion && IsFwdModels();

		double forwardValueMaturityTime	= ComputeFwdAtTime( settlementTime );
		double qParameter = GetModelParams()->GetModelParam( ARM_ModelParamType::QParameter).ToCurveModelParam().GetCurve()->Interpolate(expiryTime);
		
		if( evalTime<=K_NEW_DOUBLE_TOL
			|| states == ARM_PricingStatesPtr(NULL) )
		{
			if(isFwdAndInteg && !states->GetPricingStatesContextVector().IsNull() && states->GetPricingStatesContextVector()->size())
				return states->GetPricingStatesContext( GetModelRank() )->ToQ1FPricingStatesContext()->GetFwdValue();
			else
				return ARM_VectorPtr( new ARM_GP_Vector(1,forwardValueMaturityTime) );
		}
		else
		{
			size_t stateSize = states->size();
			ARM_GP_VectorPtr fwdValue;
			if(isFwdAndInteg)
				fwdValue = states->GetPricingStatesContext( GetModelRank() )->ToQ1FPricingStatesContext()->GetFwdValue();
			else
				fwdValue = ARM_GP_VectorPtr( new ARM_GP_Vector(stateSize,forwardValueMaturityTime) );

			ARM_GP_VectorPtr fwdVector = ARM_GP_VectorPtr( new ARM_GP_Vector(stateSize, 0.0 ));
			size_t modelNb = GetModelNb();

			if( fabs(qParameter)>K_NEW_DOUBLE_TOL)
			{
				double stateLocalVar = static_cast< const ARM_ModelParamsQ1F* const>(GetModelParams())->StateLocalVariance(0.0,evalTime,evalTime);
				double drift = 0.5*qParameter*stateLocalVar;
				for( size_t i=0; i<stateSize; ++i )
					(*fwdVector )[i] = MappingFunction( states->GetModelState(i,modelNb)-drift,(*fwdValue)[i],qParameter);
			}
			else
			{
				for( size_t i=0; i<stateSize; ++i )
					
					(*fwdVector )[i] = MappingFunction( states->GetModelState(i,modelNb), (*fwdValue)[i],qParameter);
			}

			return fwdVector;
		}
	}

    size_t i,nbStates = states->size();
    size_t modelNb = GetModelNb();
    ARM_GP_VectorPtr domDf,forDf;
    ARM_GP_VectorPtr fwdFxValues;

    if( (GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_FWDLOOKING)||(GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_FWDBCKWDLOOKING)||(abs(lambda)>1) )//abs(lambda)>1 implies no Qmaping
    {
        /// Spot Fx is path-dependently already diffused at current time
        if(settlementTime <= evalTime + K_NEW_DOUBLE_TOL)
        {
            domDf = ARM_GP_VectorPtr( new ARM_GP_Vector(nbStates,1.0) );
            forDf = domDf;
        }
        else
        {
            domDf = itsIRDomModel->DiscountFactor(itsIRDomModel->GetModelName(),evalTime,settlementTime,states);
            forDf = itsIRForModel->DiscountFactor(itsIRForModel->GetModelName(),evalTime,settlementTime,states);
        }
		fwdFxValues = ARM_GP_VectorPtr(new ARM_GP_Vector(nbStates, 0.0));

		double drift = 0.0;
		if (GetNumeraire()->GetType() == ARM_Numeraire::TerminalZc)
		{
			double numTime = GetNumeraire()->GetMaturity();

			ARM_GP_Matrix correlMatrix = GetCorrelMatrix().Interpolate(evalTime);
			double corrDomFor    = correlMatrix(DomModel,ForModel);
			double corrDomFx    = correlMatrix(DomModel,FxModel);
			double corrForFx    = correlMatrix(ForModel,FxModel);

			const ARM_ModelParamsHW1FStd* const domParam   = static_cast<const ARM_ModelParamsHW1FStd*>( itsIRDomModel->GetModelParams() );
			const ARM_ModelParamsHW1FStd* const forParam   = static_cast<const ARM_ModelParamsHW1FStd*>( itsIRForModel->GetModelParams() );
			const ARM_CurveModelParam& fxVolParam = static_cast< const ARM_CurveModelParam& >(GetModelParams()->GetModelParam(ARM_ModelParamType::QVol));

			double betatT = domParam->BetatT(expiryTime,numTime);

			double covarFXDom = ARM_ModelParamsHW1FStd::HW1FEqFxStateCovariance(fxVolParam,domParam,0.0,expiryTime,expiryTime)*corrDomFor;
			double covarZcForDom = ARM_ModelParamsHW1FStd::HW1FStateZcCovariance(domParam,forParam,0.0,expiryTime,expiryTime,expiryTime)*corrDomFx;
			double covarZcDomDom = ARM_ModelParamsHW1FStd::HW1FStateZcCovariance(domParam,domParam,0.0,expiryTime,expiryTime,expiryTime);

			drift = exp(betatT*(covarZcDomDom-covarZcForDom+covarFXDom));
		}

		if (!itsLocalFunctional.IsNull())
		{
			if (settlementTime != evalTime)
				ARM_THROW(ERR_INVALID_ARGUMENT, "ARM_QModel1F_Fx::Forward : We cannot use the forward keyword in the case settlement time != eval time with LocalFunctional.");

			for( i=0; i<nbStates; ++i )
				(*fwdFxValues )[i] = states->GetModelState(i,modelNb);
			
			size_t index = IdxFromValue(itsLocalFunctional->ResetTimes(),evalTime);

			double resetTime = evalTime/K_YEAR_LEN;
			double vol = itsLocalFunctional->GetVol(index);
			double fwd = itsLocalFunctional->GetFwd(index);
			
			for( i=0; i<nbStates; ++i )
				(*fwdFxValues)[i] = (log((*fwdFxValues )[i]/fwd)+0.5*vol*vol*resetTime)/vol/sqrt(resetTime);

			fwdFxValues = itsLocalFunctional->Func(evalTime,fwdFxValues);

		}
		else
		{
			for( i=0; i<nbStates; ++i )
				(*fwdFxValues )[i] = states->GetModelState(i,modelNb) * (*forDf)[i] / (*domDf)[i];
		}
    }
    else
    {
        /// Yield curves are stochastic : spot is computed using the mapping function then
        /// calibrated domestic/foreign zero-coupons give the forward Fx

#if defined( __GP_STRICT_VALIDATION )
        if(expiryTime < evalTime - K_NEW_DOUBLE_TOL)
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": try to compute a forward Fx expiring strictly before evalTime" );
#endif

        if(settlementTime <= evalTime + K_NEW_DOUBLE_TOL)
        {
            domDf = ARM_GP_VectorPtr( new ARM_GP_Vector(nbStates,1.0) );
            forDf = domDf;
        }
        else
        {
            domDf = itsIRDomModel->DiscountFactor(itsIRDomModel->GetModelName(),evalTime,settlementTime,states);
            forDf = itsIRForModel->DiscountFactor(itsIRForModel->GetModelName(),evalTime,settlementTime,states);
        }


	    double fwdFx0 = ComputeFwdAtTime(evalTime);
        double qFx = GetModelParams()->GetModelParam( ARM_ModelParamType::QParameter).ToCurveModelParam().GetCurve()->Interpolate(evalTime);

		fwdFxValues = ARM_GP_VectorPtr(new ARM_GP_Vector(nbStates, 0.0));

        /// No drift is computed here. The integrated correction should use the
        /// time dependency of the q parameter so it is simpler to add it through
        /// the markovian drift correction
        /// Convexity correction not computed yet (settlementTime may differ from payTime)
        double spotFx;
		double drift = 0.0;
		double lagDrift = 1.0;

		if (GetNumeraire()->GetType() == ARM_Numeraire::TerminalZc)
		{
			double numTime = GetNumeraire()->GetMaturity();

			ARM_GP_Matrix correlMatrix = GetCorrelMatrix().Interpolate(evalTime);
			double corrDomFor    = correlMatrix(DomModel,ForModel);
			double corrDomFx    = correlMatrix(DomModel,FxModel);
			double corrForFx    = correlMatrix(ForModel,FxModel);

			const ARM_ModelParamsHW1FStd* const domParam   = static_cast<const ARM_ModelParamsHW1FStd*>( itsIRDomModel->GetModelParams() );
			const ARM_ModelParamsHW1FStd* const forParam   = static_cast<const ARM_ModelParamsHW1FStd*>( itsIRForModel->GetModelParams() );
			const ARM_ModelParamsHW1FStd* const fxParam		= static_cast<const ARM_ModelParamsHW1FStd*>( GetModelParams() );
			const ARM_CurveModelParam& fxVolParam = static_cast< const ARM_CurveModelParam& >(fxParam->GetModelParam(fxParam->GetVolatilityType()));

			double varFx       = fxParam->StateLocalVariance(0.0,evalTime,evalTime);
			double varZcDom    = ARM_ModelParamsHW1F::HW1FZcCovariance( domParam, domParam, 0.0,evalTime,evalTime );
			double varZcFor    = ARM_ModelParamsHW1F::HW1FZcCovariance( forParam, forParam, 0.0,evalTime,evalTime );

			double covarZcDomZcFor = ARM_ModelParamsHW1F::HW1FZcCovariance( domParam, forParam, 0.0,evalTime,evalTime );
			double covarZcDomFx    = ARM_ModelParamsHW1F::HW1FEqFxZcCovariance( fxVolParam, domParam, 0.0,evalTime,evalTime );
			double covarZcForFx    = ARM_ModelParamsHW1F::HW1FEqFxZcCovariance( fxVolParam, forParam, 0.0,evalTime,evalTime );

			drift = varFx + varZcDom + varZcFor
			- 2*corrDomFor*covarZcDomZcFor
			- 2*corrDomFx*covarZcDomFx
			+ 2*corrForFx*covarZcForFx;

			drift = 0.5*drift;

			double betatT = domParam->BetatT(expiryTime,numTime);

			double covarFXDom = ARM_ModelParamsHW1FStd::HW1FEqFxStateCovariance(fxVolParam,domParam,0.0,expiryTime,expiryTime)*corrDomFx;
			double covarZcForDom = ARM_ModelParamsHW1FStd::HW1FStateZcCovariance(domParam,forParam,0.0,expiryTime,expiryTime,expiryTime)*corrDomFor;
			double covarZcDomDom = ARM_ModelParamsHW1FStd::HW1FStateZcCovariance(domParam,domParam,0.0,expiryTime,expiryTime,expiryTime);

			lagDrift = exp(betatT*(covarZcDomDom-covarZcForDom-covarFXDom));
		}

		for( i=0; i<nbStates; ++i )
        {
            spotFx = MappingFunction(states->GetModelState(i,modelNb)-drift,fwdFx0,qFx);
			(*fwdFxValues )[i] = spotFx * lagDrift * (*forDf)[i] / (*domDf)[i];
        }
    }

    return fwdFxValues;
}


////////////////////////////////////////////////////
///	Class   : ARM_QModel1F_Fx
///	Routine	: VarianceFwdFx
///	Returns : a double
///	Action  : Computes the variance in [a,b] of the forward FX
///           S(.T)=S(.)Bf(.,T)/Bd(.,T) (T = fwd Fx settlement time)
///           S follows a Q model dynamics and the naïve approximation
///           at time 0, S(0,T).Bf(0,t,T)/Bd(0,t,T), is used in the
///           exact volatility formula for forward FX
///	          If needed, it also computes the covariance in [a,b]
///           of the forward Zc Bd(.,T)/Bf(.,T) and the forward FX
////////////////////////////////////////////////////
double ARM_QModel1F_Fx::VarianceFwdFx(double a, double b, double settlementTime,
                                      const ARM_QModel1F* domRefModel,
                                      const ARM_QModel1F* forRefModel,
                                      const ARM_ModelParamsQ1F* fxParams,
                                      bool isFwdZcFwdFxCov,
                                      double& fwdZcFwdFxCov) const
{
    const ARM_ModelParamsQ1F* domZcParams   = static_cast< const ARM_ModelParamsQ1F*>(domRefModel->GetModelParams());
    const ARM_ModelParamsQ1F* forZcParams   = static_cast< const ARM_ModelParamsQ1F* >(forRefModel->GetModelParams());

	int forModelIdx = forRefModel->GetModelNb();
	int fxModelIdx = GetModelNb();

	double spotFxVar = fxParams->StateLocalVariance(a,b,b);

	ARM_CurveMatrix curveCorrel = GetCorrelMatrix();
	ARM_GP_Matrix correlMatrix;

	ARM_GP_Vector correlTimes = curveCorrel.GetTimes();
	int startCorrelTimeIdx = 0;
	int endCorrelTimeIdx = 0;
	int i;
	ARM_GP_Vector timeSteps;
	
	while ((startCorrelTimeIdx < correlTimes.size())&&(correlTimes[startCorrelTimeIdx] <= a))
				startCorrelTimeIdx++;

	endCorrelTimeIdx = startCorrelTimeIdx;

	while ((endCorrelTimeIdx < correlTimes.size())&&(correlTimes[endCorrelTimeIdx] < b))
		endCorrelTimeIdx++;

	timeSteps.push_back(a);

	for (i = startCorrelTimeIdx; i < endCorrelTimeIdx; ++i)
		timeSteps.push_back(correlTimes[i]);

	timeSteps.push_back(b);

	correlMatrix = curveCorrel.Interpolate(a);

    double domForCorr = correlMatrix(DomModel,forModelIdx);
	double domFxCorr = correlMatrix(DomModel,fxModelIdx);
	double forFxCorr = correlMatrix(forModelIdx,fxModelIdx);

    double fwdFxVar=0.0,domZcVar=0.0,forZcVar=0.0,domForZcCovar=0.0,domZcFxCovar=0.0,forZcFxCovar=0.0;
    double domFor,domFx,forFx;
    fwdZcFwdFxCov = 0.0;

    if(a + K_NEW_DOUBLE_TOL < b)
    {
        /// Restore fxVolCurve (underlying gaussian process volatility in fact)
        const ARM_CurveModelParam& fxVolParam = static_cast< const ARM_CurveModelParam& >(fxParams->GetModelParam(ARM_ModelParamType::QVol));

        if(domRefModel->IsDegenerateInHW() && forRefModel->IsDegenerateInHW())
        {	

			fwdFxVar = spotFxVar;
			for (i = 0; i < timeSteps.size()-1; ++i)
			{
				correlMatrix = curveCorrel.Interpolate(timeSteps[i]);

				domForCorr = correlMatrix(DomModel,forModelIdx);
				domFxCorr = correlMatrix(DomModel,fxModelIdx);
				forFxCorr = correlMatrix(forModelIdx,fxModelIdx);

				/// Both IR models are degenerated to H&W 1F
				domZcVar        = ARM_ModelParamsHW1F::HW1FZcCovariance(domZcParams,domZcParams,timeSteps[i],timeSteps[i+1],settlementTime);
				forZcVar        = ARM_ModelParamsHW1F::HW1FZcCovariance(forZcParams,forZcParams,timeSteps[i],timeSteps[i+1],settlementTime);
				domForZcCovar   = ARM_ModelParamsHW1F::HW1FZcCovariance(domZcParams,forZcParams,timeSteps[i],timeSteps[i+1],settlementTime);
				domZcFxCovar    = ARM_ModelParamsHW1F::HW1FEqFxZcCovariance(fxVolParam,domZcParams,timeSteps[i],timeSteps[i+1],settlementTime);
				forZcFxCovar    = ARM_ModelParamsHW1F::HW1FEqFxZcCovariance(fxVolParam,forZcParams,timeSteps[i],timeSteps[i+1],settlementTime);

				domFor      = domForCorr*domForZcCovar;
				domFx       = domFxCorr*domZcFxCovar;
				forFx       = forFxCorr*forZcFxCovar;

				///  Var(fwdFx(.,T))
				fwdFxVar    += forZcVar + domZcVar + 2*(forFx - domFx - domFor);

				/// Covar(Bd(.T)/Bf(.T),fwdFx(.,T))
				if(isFwdZcFwdFxCov)
					fwdZcFwdFxCov -= (forZcVar + domZcVar - 2*domFor + forFx - domFx);
			}
        }
        else
        {
            /// The interval [evalTime,maturityTime] is cut in a predetermined schedule to optimise
            /// the double integration needed to compute all covariances
            /// At the moment this schedule is rebuilt each time and decay factors recomputed !
            double yfAToB = (b-a)/K_YEAR_LEN;
            size_t i,nbSteps = static_cast<size_t>(floor(yfAToB*DECAY_NBSTEP_PER_YEAR));
            if(nbSteps < 1)
                nbSteps = 1;
            double step = (b-a)/nbSteps, yfStep = yfAToB/nbSteps;
            double yfFwdTime = settlementTime/K_YEAR_LEN;

            double domDecay,forDecay;
            double fwdFxVar = spotFxVar;
            double fwdZcFwdFxCov = 0.0;
            double lastTime = a, time = a+step, yfTime = time/K_YEAR_LEN;
            for(i=0;i<nbSteps;++i)
            {
                domDecay        = domRefModel->DecayFactor(yfTime, yfFwdTime, yfTime, domZcParams );
                forDecay        = forRefModel->DecayFactor(yfTime, yfFwdTime, yfTime, forZcParams );
                domZcVar        = ARM_ModelParamsHW1F::HW1FStateCovariance(domZcParams, domZcParams, lastTime, time, time );
                forZcVar        = ARM_ModelParamsHW1F::HW1FStateCovariance(forZcParams, forZcParams, lastTime, time, time );
                domForZcCovar   = ARM_ModelParamsHW1F::HW1FStateCovariance(domZcParams, forZcParams, lastTime, time, time );
                domZcFxCovar    = ARM_ModelParamsHW1F::HW1FStateCovariance(fxParams, domZcParams, lastTime, time, time );
                forZcFxCovar    = ARM_ModelParamsHW1F::HW1FStateCovariance(fxParams, forZcParams, lastTime, time, time );

                domFor      = domDecay*(domDecay*domZcVar - 2*domForCorr*forDecay*domForZcCovar) + forDecay*forDecay*forZcVar;
                domFx       = domDecay*domFxCorr*domZcFxCovar;
                forFx       = forDecay*forFxCorr*forZcFxCovar;

                ///  Var(fwdFx(.,T))
                fwdFxVar += domFor + 2*domFx - 2*forFx; /// sign are inverted because zc vol is - sig(t).Integ{u=t->T,f(0,u)exp(-mrs(u-t))du}

                /// Covar(Bd(.T)/Bf(.T),fwdFx(.,T))
                if(isFwdZcFwdFxCov)
                    fwdZcFwdFxCov -= domFor + domFx - forFx;

                lastTime    = time;
                time        += step;
                yfTime      += yfStep;
            }
        }
    }
    else
        fwdFxVar = 0.0;

    return fwdFxVar;
}


////////////////////////////////////////////////////
///	Class   : ARM_QModel1F_Fx
///	Routine	: CallVectorial
///	Returns : a vector of call vectorial
///	Action  : computes fx call
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QModel1F_Fx::CallVectorial(
		const string& modelName,
        double evalTime,
	    double expiryTime,
	    double settlementTime,
		const ARM_GP_Vector& strikePerState,
		int callPut,
	    double payTime,
        const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context) const
{
    /// Restore reference IR models
	ARM_QModel1F* domRefModel = GetStochasticModel( itsIRDomModel );
	ARM_QModel1F* forRefModel = GetStochasticModel( itsIRForModel );

    if(expiryTime < evalTime - K_NEW_DOUBLE_TOL)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": the fx option resets before event time" );

    bool isIVCalc;
    if(isIVCalc = (expiryTime <= evalTime + K_NEW_DOUBLE_TOL))
        expiryTime = evalTime;  // force evaluation at reset time

	if( expiryTime<=K_NEW_DOUBLE_TOL || !forRefModel || !domRefModel || GetCorrelMatrix().empty() )
    {
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": CallVectorial is not implemented when the forward rate are deterministic." );
    }


	size_t i,nbStates = (evalTime > K_NEW_DOUBLE_TOL && states != ARM_PricingStatesPtr(NULL))
                        ? states->size(): 1;

#if defined( __GP_STRICT_VALIDATION )
	if( strikePerState.size() != nbStates )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " strike vector size != states size" );
#endif

    double yfEval   = evalTime/K_YEAR_LEN;
    double yfExpiry = expiryTime/K_YEAR_LEN;

    const ARM_ModelParamsQ1F* fxParams = static_cast< const ARM_ModelParamsQ1F*>(GetModelParams());

    double fwdFxQ,zcFxCov;
    double fwdFxVarTtoM = VarianceFwdFx(evalTime,expiryTime,settlementTime,domRefModel,forRefModel,fxParams,false,zcFxCov);

	const ARM_Curve& spotFxQCurve = *(fxParams->GetModelParam( ARM_ModelParamType::QParameter).ToCurveModelParam().GetCurve());
    if(spotFxQCurve == 1.0)
    {
        /// Standard case : pure lognormal forward FX vol
        fwdFxQ = 1.0;
    }
    else if(!isIVCalc)
    {
        /// Compute a correction to get an approximated displaced forward FX volatility
        /// and shift (Piterbarg's "FX volatility skew" 2005 paper)

        /// Merge IR & FX vol schedules
        const ARM_ModelParamsQ1F* domParams   = static_cast< const ARM_ModelParamsQ1F*>(domRefModel->GetModelParams());
        const ARM_GP_Vector& domTimes     = ((ARM_CurveModelParam&) domParams->GetModelParam(domParams->GetVolatilityType())).GetCurve()->GetAbscisses();

        const ARM_ModelParamsQ1F* forParams   = static_cast< const ARM_ModelParamsQ1F* >(forRefModel->GetModelParams());
        const ARM_GP_Vector& forTimes     = ((ARM_CurveModelParam&) forParams->GetModelParam(forParams->GetVolatilityType())).GetCurve()->GetAbscisses();

        const ARM_GP_Vector& fxTimes     = ((ARM_CurveModelParam&) fxParams->GetModelParam(fxParams->GetVolatilityType())).GetCurve()->GetAbscisses();

		const ARM_GP_Vector correlTimes = GetCorrelMatrix().GetTimes();

        ARM_GP_VectorPtr times1(MergeSortedVectorNoDuplicates(domTimes,forTimes));
		ARM_GP_VectorPtr times2(MergeSortedVectorNoDuplicates(*times1,fxTimes));
        ARM_GP_VectorPtr times(MergeSortedVectorNoDuplicates(*times2,correlTimes));

        /// Find n(a) & n(b) such that Tn(a)-1 < evalTime <= Tn(a) and Tn(b)-1< expiryTime <= Tn(b)
	    int na = lower_boundPosWithPrecision(*times,evalTime);
	    int nb = lower_boundPosWithPrecision(*times,expiryTime);

		ARM_GP_Matrix correlMatrix;

	    double domFxCorr	= 0.0;
	    double forFxCorr	= 0.0;

        double fwdFxVarTtot,fwdZcFwdFxCovTtot;
        double lastt=evalTime,yfLastt=yfEval,t,yft,tt,lasttt,zcFxCov,w;
        double volXSt,khiRatiot,spotFxQt,etat,volDomZct,volForZct;
        fwdFxQ = 0.0;
        fwdFxVarTtot = 0.0;
        fwdZcFwdFxCovTtot = 0.0;

        size_t n,nbTimes=times->size();
        size_t i,nbPoints;
        for(n=na;n<=nb;++n)
        {
            t   = (n<nbTimes ? ((*times)[n]>expiryTime ? expiryTime : (*times)[n]): expiryTime);
            yft = t/K_YEAR_LEN;

			correlMatrix = GetCorrelMatrix().Interpolate(t);

			domFxCorr = correlMatrix(DomModel,FxModel);
			forFxCorr = correlMatrix(ForModel,FxModel);

            /// A numerical integration is needed to get the equivalent shift on forward FX
            nbPoints=static_cast<int>(ceil(yft-yfLastt)*GL_PY_NBPOINTS);
            if(nbPoints<GL_MIN_NBPOINTS)
                nbPoints = GL_MIN_NBPOINTS;

            GaussLegendre_Coefficients GLPoints(nbPoints,yfLastt,yft);

            /// Gauss-Legendre summation between each [Ti,Ti+1]
            lasttt = lastt;
            for(i=0;i<nbPoints;++i)
            {
                tt= GLPoints.get_point(i)*K_YEAR_LEN;
                w = GLPoints.get_weight(i);

                fwdFxVarTtot        += VarianceFwdFx(lasttt,tt,settlementTime,domRefModel,forRefModel,fxParams,true,zcFxCov);
                fwdZcFwdFxCovTtot   += zcFxCov;
                khiRatiot           = fwdZcFwdFxCovTtot/fwdFxVarTtot;

	            spotFxQt        = fxParams->GetModelParam( ARM_ModelParamType::QParameter).ToCurveModelParam().GetCurve()->Interpolate(tt);
	            volXSt          = fxParams->GetModelParam( ARM_ModelParamType::QVol).ToCurveModelParam().GetCurve()->Interpolate(tt);
                etat            = -volXSt*(1+khiRatiot)*(1-spotFxQt);

                volDomZct       = domRefModel->VolZc(tt,settlementTime);
                volForZct       = forRefModel->VolZc(tt,settlementTime);

                fwdFxQ          += w * fwdFxVarTtot*etat*(volXSt + volForZct*forFxCorr - volDomZct*domFxCorr);

                lasttt = tt;
            }

            if(n<nb)
            {
                fwdFxVarTtot        += VarianceFwdFx(lasttt,t,settlementTime,domRefModel,forRefModel,fxParams,true,zcFxCov);
                fwdZcFwdFxCovTtot   += zcFxCov;
            }
            else
                fwdFxVarTtot        += VarianceFwdFx(lasttt,t,settlementTime,domRefModel,forRefModel,fxParams,false,zcFxCov);

            lastt   = t;
            yfLastt = yft;
        }
		double error = fabs(fwdFxVarTtoM-fwdFxVarTtot);
		if(fwdFxVarTtoM > 1.0e-1)
			error/=fwdFxVarTtoM;

        if(fabs(error)>K_NEW_DOUBLE_TOL)
    		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
            " : inconsistency in variances computation to get an equivalent forward FX Q parameter" );

        fwdFxQ = 1.0 + 2.0*fwdFxQ/(fwdFxVarTtoM*fwdFxVarTtoM);
    }

	/// Compute option values without memory allocation by replacing fwd by option values
    ARM_VectorPtr zcPay = itsIRDomModel->DiscountFactor(itsIRDomModel->GetModelName(),evalTime,payTime,states);
	ARM_VectorPtr fwdFx = Forward(modelName, evalTime, expiryTime, settlementTime, payTime, states );

    ARM_GP_Vector* fxOption = new ARM_GP_Vector(nbStates);

    double yfTtoM = yfExpiry-yfEval;
	double fwdFxVol=0.0;
    if(isIVCalc)
    {
        /// Just compute the option intrinsic value
	    for(i=0; i<nbStates; ++i )
		    (*fxOption)[i] = CC_Max<double>(callPut*((*fwdFx)[i]-strikePerState[i]),0.0)*(*zcPay)[i];

        if(itsIsImpliedVolCalc)
        {
            itsFwdFxVol.resize(nbStates);
	        for(i=0; i<nbStates; ++i )
                itsFwdFxVol[i]=0.0;
        }
    }
    else
    {
        /// Compute a Q model price formula
        fwdFxVol = sqrt(fwdFxVarTtoM/yfTtoM);
        if(itsIsImpliedVolCalc)
        {
            itsFwdFxVol.resize(nbStates);
            double accuracy,fwdFxValue,sTtoM=sqrt(yfTtoM);
	        for(i=0; i<nbStates; ++i )
            {
                fwdFxValue = (*fwdFx)[i];
		        (*fxOption)[i] = QModelAnalytics::BSQFunction( fwdFxValue, strikePerState[i], fwdFxVol, yfTtoM, fwdFxQ, (*zcPay)[i], callPut, fwdFxValue );

                if(fwdFxQ != 1.0)
                {
                    accuracy = 1.0e-4*(*fxOption)[i];
                    if(accuracy < 1.0e-8)
                        accuracy = 1.0e-8;
                    itsFwdFxVol[i] = ARM_CF_BS_Formula::callimplicit_totalvolatility(fwdFxValue, (*zcPay)[i],
                        strikePerState[i],callPut,(*fxOption)[i],accuracy) / sTtoM;
                }
                else
                    itsFwdFxVol[i]=fwdFxVol;
            }
        }
        else
        {
	        for(i=0; i<nbStates; ++i )
		        (*fxOption)[i] = QModelAnalytics::BSQFunction( (*fwdFx)[i], strikePerState[i], fwdFxVol, yfTtoM, fwdFxQ, (*zcPay)[i], callPut, (*fwdFx)[i] );
        }
    }

//FILE* f=fopen("c:\\temp\\dumpQ1F_FX.txt","a");
//fprintf(f,"t=%5.0lf\tfwdFx=%15.10lf\tfwdQ=%15.10lf\tfwdQvol=%15.10lf\n",expiryTime,(*fwdFx)[0],fwdFxQ,fwdFxVol);
//fclose(f);

    if(context)
    {
        /// Save pricing context
        ARM_PricingCallContext* callContext = context->ToCallContext();
        if(fwdFxQ == 1.0)
        {
            callContext->SetVol(fwdFxVol);
            callContext->SetShift(1.0);
        }
        else
        {
            callContext->SetVol(fwdFxVol);
            callContext->SetShift(fwdFxQ);
        }
        callContext->SetForward(fwdFx);
    }

	return ARM_VectorPtr(fxOption);
}

////////////////////////////////////////////////////
///	Class   : ARM_QModel1F_Fx
///	Routine	: CallVectorial
///	Returns : a vector of call vectorial
///	Action  : computes fx call
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QModel1F_Fx::ComputeFwdFXModelParam(
		double evalTime,
	    double expiryTime,
	    double settlementTime,
		ARM_PricingContext* context) const
{
    /// Restore reference IR models
	ARM_QModel1F* domRefModel = GetStochasticModel( itsIRDomModel );
	ARM_QModel1F* forRefModel = GetStochasticModel( itsIRForModel );

    if(expiryTime < evalTime - K_NEW_DOUBLE_TOL)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": the fx option resets before event time" );

    bool isIVCalc;
    if(isIVCalc = (expiryTime <= evalTime + K_NEW_DOUBLE_TOL))
        expiryTime = evalTime;  // force evaluation at reset time

	if( expiryTime<=K_NEW_DOUBLE_TOL || !forRefModel || !domRefModel || GetCorrelMatrix().empty() )
    {
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": CallVectorial is not implemented when the forward rate are deterministic." );
    }

	double yfEval   = evalTime/K_YEAR_LEN;
    double yfExpiry = expiryTime/K_YEAR_LEN;

    const ARM_ModelParamsQ1F* fxParams = static_cast< const ARM_ModelParamsQ1F*>(GetModelParams());

    double fwdFxQ,zcFxCov;
    double fwdFxVarTtoM = VarianceFwdFx(evalTime,expiryTime,settlementTime,domRefModel,forRefModel,fxParams,false,zcFxCov);

	const ARM_Curve& spotFxQCurve = *(fxParams->GetModelParam( ARM_ModelParamType::QParameter).ToCurveModelParam().GetCurve());
    if(spotFxQCurve == 1.0)
    {
        /// Standard case : pure lognormal forward FX vol
        fwdFxQ = 1.0;
    }
    else if(!isIVCalc)
    {
        /// Compute a correction to get an approximated displaced forward FX volatility
        /// and shift (Piterbarg's "FX volatility skew" 2005 paper)

        /// Merge IR & FX vol schedules
        const ARM_ModelParamsQ1F* domParams   = static_cast< const ARM_ModelParamsQ1F*>(domRefModel->GetModelParams());
        const ARM_GP_Vector& domTimes     = ((ARM_CurveModelParam&) domParams->GetModelParam(domParams->GetVolatilityType())).GetCurve()->GetAbscisses();

        const ARM_ModelParamsQ1F* forParams   = static_cast< const ARM_ModelParamsQ1F* >(forRefModel->GetModelParams());
        const ARM_GP_Vector& forTimes     = ((ARM_CurveModelParam&) forParams->GetModelParam(forParams->GetVolatilityType())).GetCurve()->GetAbscisses();

        const ARM_GP_Vector& fxTimes     = ((ARM_CurveModelParam&) fxParams->GetModelParam(fxParams->GetVolatilityType())).GetCurve()->GetAbscisses();

		const ARM_GP_Vector correlTimes = GetCorrelMatrix().GetTimes();

        ARM_GP_VectorPtr times1(MergeSortedVectorNoDuplicates(domTimes,forTimes));
		ARM_GP_VectorPtr times2(MergeSortedVectorNoDuplicates(*times1,fxTimes));
        ARM_GP_VectorPtr times(MergeSortedVectorNoDuplicates(*times2,correlTimes));

        /// Find n(a) & n(b) such that Tn(a)-1 < evalTime <= Tn(a) and Tn(b)-1< expiryTime <= Tn(b)
	    int na = lower_boundPosWithPrecision(*times,evalTime);
	    int nb = lower_boundPosWithPrecision(*times,expiryTime);

		ARM_GP_Matrix correlMatrix;

	    double domFxCorr	= 0.0;
	    double forFxCorr	= 0.0;

        double fwdFxVarTtot,fwdZcFwdFxCovTtot;
        double lastt=evalTime,yfLastt=yfEval,t,yft,tt,lasttt,zcFxCov,w;
        double volXSt,khiRatiot,spotFxQt,etat,volDomZct,volForZct;
        fwdFxQ = 0.0;
        fwdFxVarTtot = 0.0;
        fwdZcFwdFxCovTtot = 0.0;

        size_t n,nbTimes=times->size();
        size_t i,nbPoints;
        for(n=na;n<=nb;++n)
        {
            t   = (n<nbTimes ? ((*times)[n]>expiryTime ? expiryTime : (*times)[n]): expiryTime);
            yft = t/K_YEAR_LEN;

			correlMatrix = GetCorrelMatrix().Interpolate(t);

			domFxCorr = correlMatrix(DomModel,FxModel);
			forFxCorr = correlMatrix(ForModel,FxModel);

            /// A numerical integration is needed to get the equivalent shift on forward FX
            nbPoints=static_cast<int>(ceil(yft-yfLastt)*GL_PY_NBPOINTS);
            if(nbPoints<GL_MIN_NBPOINTS)
                nbPoints = GL_MIN_NBPOINTS;

            GaussLegendre_Coefficients GLPoints(nbPoints,yfLastt,yft);

            /// Gauss-Legendre summation between each [Ti,Ti+1]
            lasttt = lastt;
            for(i=0;i<nbPoints;++i)
            {
                tt= GLPoints.get_point(i)*K_YEAR_LEN;
                w = GLPoints.get_weight(i);

                fwdFxVarTtot        += VarianceFwdFx(lasttt,tt,settlementTime,domRefModel,forRefModel,fxParams,true,zcFxCov);
                fwdZcFwdFxCovTtot   += zcFxCov;
                khiRatiot           = fwdZcFwdFxCovTtot/fwdFxVarTtot;

	            spotFxQt        = fxParams->GetModelParam( ARM_ModelParamType::QParameter).ToCurveModelParam().GetCurve()->Interpolate(tt);
	            volXSt          = fxParams->GetModelParam( ARM_ModelParamType::QVol).ToCurveModelParam().GetCurve()->Interpolate(tt);
                etat            = -volXSt*(1+khiRatiot)*(1-spotFxQt);

                volDomZct       = domRefModel->VolZc(tt,settlementTime);
                volForZct       = forRefModel->VolZc(tt,settlementTime);

                fwdFxQ          += w * fwdFxVarTtot*etat*(volXSt + volForZct*forFxCorr - volDomZct*domFxCorr);

                lasttt = tt;
            }

            if(n<nb)
            {
                fwdFxVarTtot        += VarianceFwdFx(lasttt,t,settlementTime,domRefModel,forRefModel,fxParams,true,zcFxCov);
                fwdZcFwdFxCovTtot   += zcFxCov;
            }
            else
                fwdFxVarTtot        += VarianceFwdFx(lasttt,t,settlementTime,domRefModel,forRefModel,fxParams,false,zcFxCov);

            lastt   = t;
            yfLastt = yft;
        }
		double error = fabs(fwdFxVarTtoM-fwdFxVarTtot);
		if(fwdFxVarTtoM > 1.0e-1)
			error/=fwdFxVarTtoM;

        if(fabs(error)>K_NEW_DOUBLE_TOL)
    		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
            " : inconsistency in variances computation to get an equivalent forward FX Q parameter" );

        fwdFxQ = 1.0 + 2.0*fwdFxQ/(fwdFxVarTtoM*fwdFxVarTtoM);
    }

	double fwdFxVol = sqrt(fwdFxVarTtoM / (yfExpiry-yfEval));

	ARM_GP_Vector* fwdFXModelParam = new ARM_GP_Vector(2);
	(*fwdFXModelParam)[0] = fwdFxQ; 
	(*fwdFXModelParam)[1] = fwdFxVol;
    
    return ARM_VectorPtr(fwdFXModelParam);
}

////////////////////////////////////////////////////
///	Class  : ARM_QModel1F_Fx
///	Routine: DigitalVectorial
///	Returns: a vector of digital vectorial
///	Action : computes equity or fx digital
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QModel1F_Fx::DigitalVectorial(
	const string& modelName,
    double evalTime,
	double expiryTime,
	double settlementTime,
	const ARM_GP_Vector& strikePerState,
	double notional,
	int callPut,
	double payTime,
	ARM_DigitType digitType,
	double epsilon,
    const ARM_PricingStatesPtr& states,
    ARM_PricingContext* context) const
{
    /// FIX FIX : review formula with payTime !!

	ARM_GP_VectorPtr result = Forward(modelName, evalTime, expiryTime, settlementTime, payTime, states );

#if defined(__GP_STRICT_VALIDATION)
	if( strikePerState.size() != result->size() )
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
                ": strikePerState.size() != forward.size()" );
	}
#endif

	if(	expiryTime<=K_NEW_DOUBLE_TOL )
	{
		/// get the intrinsic value
		size_t payoffSize = result->size();
		for( size_t i=0; i<payoffSize; ++i )
			(*result)[i] = callPut*(*result)[i]>callPut*strikePerState[i];
	}
	else
	{
		ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();
		size_t nbStates = states!=ARM_PricingStatesPtr(NULL)? states->size(): 1;

	    double qParameter	= GetModelParams()->GetModelParam( ARM_ModelParamType::QParameter).ToCurveModelParam().GetCurve()->Interpolate(expiryTime);
		double qSigma		= static_cast< const ARM_ModelParamsQ1F* const>(GetModelParams())->StateLocalVariance(evalTime,expiryTime,expiryTime);
		double time			= (expiryTime-evalTime)/K_YEAR_LEN;
		qSigma				= sqrt(qSigma/time);

		double zcT	= ZcCurve->DiscountPrice(settlementTime/K_YEAR_LEN);
		double zct	= ZcCurve->DiscountPrice(evalTime/K_YEAR_LEN);
		double df	= zcT/zct;

		/// get the intrinsic value
		size_t payoffSize = result->size();
		for( size_t i=0; i<payoffSize; ++i )
			(*result)[i] = QModelAnalytics::DigitalBSQFunction((*result)[i], strikePerState[i], qSigma, time, qParameter, df, callPut, (*result)[i] );

        if(context)
        {
            /// Save pricing context
            ARM_PricingCallContext* callContext = context->ToCallContext();
            callContext->SetVol(qSigma*qParameter);
	        ARM_GP_VectorPtr fwd0 = Forward(modelName, 0.0, expiryTime, settlementTime, payTime, ARM_PricingStatesPtr(NULL) );
            callContext->SetShift( (*fwd0)[0]*(1.0/qParameter-1.0) );
        }
    }

	return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_MarketHybridModel
///	Routine: HybridCallVectorial
///	Returns: price an hybrid call that is a spread
///			 between a forward equity or FX strip
///			 and an IR swap :
///			 Call = O1(t)*Et[max(EqFxStrip(T)/O1(T) - (SwapRate(T)-Strike),0)] 
///			 Put = O1(t)*Et[max(SwapRate(T)-Strike - EqFxStrip(T)/O1(T),0)] 
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QModel1F_Fx::HybridCallVectorial(
	const string& modelName,
    double evalTime,
	double expiryTime,
	int callPut,
	const ARM_GP_Vector& strikesPerState,

	/// Strip of forwards FX (or equity)
	const ARM_GP_Vector& fxExpiryTimes,
	const ARM_GP_Vector& fxSettlementTimes,
	const ARM_GP_Vector& fxPayTimes,
	const ARM_GP_Vector& fxNotionals,

	/// IR Swap
	double swapResetTime,
	const ARM_GP_Vector& fixNotionals,
	const ARM_GP_Vector& floatNotionals,
	double floatStartTime,
	double floatEndTime,
	const ARM_GP_Vector& floatResetTimes,
	const ARM_GP_Vector& floatStartTimes,
	const ARM_GP_Vector& floatEndTimes,
	const ARM_GP_Vector& floatIntTerms,
	const ARM_GP_Vector& fixPayTimes,
	const ARM_GP_Vector& fixPayPeriods,

	const ARM_PricingStatesPtr& states,
    ARM_PricingContext* context) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_QModel1F_Fx::HybridCallVectorial not maintained!" );

	return ARM_VectorPtr(new ARM_GP_Vector(0));
}


////////////////////////////////////////////////////
///	Class   : ARM_QModel1F_Fx
///	Routine : DigitalVectorial
///	Returns : a vector of digital vectorial
///	Action  : computes equity or fx digital
////////////////////////////////////////////////////
void ARM_QModel1F_Fx::UpdateLinks( const ARM_MultiAssetsModel& multiAssetsModel )
{
	const ARM_ModelNameMap* modelMap = multiAssetsModel.GetModelMap();
	/// find linked models
	ARM_IntVector itsOtherModels = (*modelMap)[GetModelName()]->OtherModelRefNb();

	if( !itsOtherModels.empty() ){
		// check that there is two and only two...
		if( itsOtherModels.size() != 2 )
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_QModel1F_Fx: wrong number of otherModels" );

		itsIRDomModel = (*modelMap)[ itsOtherModels[DomModel] ]->Model();
		itsIRForModel = (*modelMap)[ itsOtherModels[ForModel] ]->Model();
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_QModel1F_Fx
///	Routines: VolatilitiesAndCorrelations
///	Returns :
///	Action  : computes the volatilities its derivatives and the correlation
////////////////////////////////////////////////////
void ARM_QModel1F_Fx::VolatilitiesAndCorrelations( const ARM_GP_Vector& timeSteps, 
	ARM_GP_MatrixPtr& vols,
	ARM_GP_MatrixPtr& d1Vols,
	ARM_GP_MatrixPtr& correls,
	bool linearVol) const
{
	if (!linearVol)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "ARM_QModel1F_Fx::VolatilitiesAndCorrelations : only linearVol case is implemented");

	/// to speed up, compute in one go the vol and vol derivatives
	ARM_GP_Vector times = GetModelParams()->GetModelParam( ARM_ModelParamType::QVol).ToCurveModelParam().GetCurve()->GetAbscisses();
	ARM_GP_Vector values= GetModelParams()->GetModelParam( ARM_ModelParamType::QVol).ToCurveModelParam().GetCurve()->GetOrdinates();

    /// First step : shift values to be equivalent as a constant left interpolation
    size_t nbT=times.size();
    ARM_GP_Vector newTimes;
    ARM_GP_Vector newValues;
    if(times[0]>0.0)
    {
        newTimes.push_back(0.0);
        newValues.push_back(values[0]);
    }
    size_t i;
    for(i=0;i<nbT;++i)
    {
        newTimes.push_back(times[i]);
        newValues.push_back(values[i+1 < nbT ? i+1 : nbT-1]);
    }

    /// Second step : use special interpolation to compute vol values and 1st derivatives
    nbT = timeSteps.size();
	ARM_GP_Vector volsVec(nbT),d1VolsVec(nbT);
    for(i=0;i<nbT;++i)
        volsVec[i] = FunctionSpecialInterpolation(timeSteps[i],newTimes,newValues);
    for(i=0;i<nbT;++i)
        d1VolsVec[i] = DerivativeSpecialInterpolation(timeSteps[i],timeSteps,volsVec) * K_YEAR_LEN;

    /// factor by line and zero correl because in one factor!
	vols	= ARM_GP_MatrixPtr( new ARM_GP_Matrix(1,volsVec.size(),&volsVec[0]) );
	d1Vols	= ARM_GP_MatrixPtr( new ARM_GP_Matrix(1,d1VolsVec.size(),&d1VolsVec[0]));
	correls	= ARM_GP_MatrixPtr( NULL );
}


////////////////////////////////////////////////////
///	Class   : ARM_QModel1F_Fx
///	Routines: EulerLocalDrifts
///	Returns :
///	Action  : computes the relative and absolute drift
////////////////////////////////////////////////////
void ARM_QModel1F_Fx::EulerLocalDrifts( const ARM_GP_Vector& timeSteps,
		ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const
{
	relativeDrifts = ARM_GP_MatrixPtr( new ARM_GP_Matrix(timeSteps.size(),1,-GetModelParams()->GetModelParam( ARM_ModelParamType::MeanReversion).GetValueAtPoint(0) ) );
    for(size_t i=0;i<timeSteps.size()-1;++i)
        (*relativeDrifts)(i,0) *= (timeSteps[i+1]-timeSteps[i])/K_YEAR_LEN;
	absoluteDrifts = ARM_GP_MatrixPtr( new ARM_GP_Matrix(timeSteps.size(),1, 0.0 ) );
}


////////////////////////////////////////////////////
///	Class  : ARM_QModel1F_Fx
///	Routine: AdviseBreakPointTimes
///	Returns: void 
///	Action : sets the corresponding suggested break point times to the model param
////////////////////////////////////////////////////
void ARM_QModel1F_Fx::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* inputModelParam, size_t factorNb )
{
	ARM_CurveModelParam* modelParam = dynamic_cast<ARM_CurveModelParam*>(inputModelParam);
    double asOfDate			= GetAsOfDate().GetJulian();
    size_t portfolioSize	= portfolio->size();  
    ARM_GP_Vector  tmpdates;
    size_t i;
	
	/// FX option or IR/FX hybrid swaption allowed
	bool isFxOption;
	bool isIrFxSwaption;
	bool isIrSwaption;
	size_t nbOpts;
	ARM_OptionPortfolio* optPort;
	for( i=0; i<portfolioSize; ++i )
	{
		isFxOption = (portfolio->GetAsset(i)->GetName() == ARM_OPTION);
		if(!isFxOption)
		{
			isIrFxSwaption = (portfolio->GetAsset(i)->GetName() == ARM_OPTIONPORTFOLIO);
			if(isIrFxSwaption)
			{
				optPort = static_cast<ARM_OptionPortfolio*>(portfolio->GetAsset(i));
				nbOpts			= optPort->GetPtf()->GetSize();
				isFxOption		= (optPort->GetPtf()->GetAsset(0)->GetName() == ARM_OPTION);
				isIrSwaption	= (optPort->GetPtf()->GetAsset(nbOpts-1)->GetName() == ARM_SWAPTION);
				if(!isFxOption || !isIrSwaption)
				{
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
						": portfolio can only contain FX options or IR/FX hybrid swaptions" );
				}		
			}
		}
		else
		{
			double expiry = portfolio->GetAsset(i)->GetExpiryDate().GetJulian() - asOfDate;
			tmpdates.push_back(expiry);

//			if (expiry > FRMVOL_LAG_THRESHOLD)
		}
	}

	if (tmpdates.size())
	modelParam->UpdateValues(&tmpdates);


    /// For Fx option keep all reset dates even if they are very close to each other
    ///... so nothing to filter
}

////////////////////////////////////////////////////
///	Class   : ARM_QModel1F_Fx
///	Routines: UnderlyingCovariance
///	Returns : double
///	Action  : compute the local covariance between
///           fromTime -> toTime of forward FXs
///           with settlements at startTimes
////////////////////////////////////////////////////
double ARM_QModel1F_Fx::UnderlyingCovariance(
    string underlyingType,
    double	fromTime,
    double	toTime,
    double	startTime1,
    double  endTime1,
    double	startTime2,
    double  endTime2,
    double	startTime3,
    double  endTime3,
    double	startTime4,
    double  endTime4) const
{
    /// FIX FIX : only fwd Fx variance is allowed !
    if(startTime1 != startTime2 || endTime1 != endTime2 || startTime1 != endTime1)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " Q1F_FX only computes local variance of a given forward Fx");

	ARM_QModel1F* domRefModel           = GetStochasticModel( itsIRDomModel );
	ARM_QModel1F* forRefModel           = GetStochasticModel( itsIRForModel );
    const ARM_ModelParamsQ1F* fxParams  = static_cast< const ARM_ModelParamsQ1F*>(GetModelParams());

    double unsedZcFxCov;
    return VarianceFwdFx(fromTime,toTime,startTime1,domRefModel,forRefModel,fxParams,false,unsedZcFxCov);
}


////////////////////////////////////////////////////
///	Class   : ARM_QModel1F_Fx
///	Routine : ValidateModelParams
///	Returns : bool
///	Action  : validate model params
////////////////////////////////////////////////////
bool ARM_QModel1F_Fx::ValidateModelParams(const ARM_ModelParams& params) const
{
	const ARM_ModelParamsQ1F_Fx* ModelParamsQ1F = dynamic_cast<const ARM_ModelParamsQ1F_Fx*>(&params);
	if( !ModelParamsQ1F )
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, "modelparams is not of type ARM_ModelParamsQ1F_Fx" );
	}
	return true;
}



////////////////////////////////////////////////////
///	Class   : ARM_QModel1F_Fx
///	Routines: MappingFunction
///	Returns :
///	Action  : mapping function of the Q model
////////////////////////////////////////////////////

double ARM_QModel1F_Fx::MappingFunction( double x, double x0, double q0 ) const
{
	double value;
	if( fabs(q0)<K_NEW_DOUBLE_TOL )
		value = x0*(1+x);
	else
		value = x0*(1.0+(exp(q0*x)-1.0)/q0);
	return value;
}

////////////////////////////////////////////////////
///	Class   : ARM_QModel1F_Fx
///	Routines: InverseMappingFunction
///	Returns :
///	Action  : inverse of the mapping function of the Q model
////////////////////////////////////////////////////

double ARM_QModel1F_Fx::InverseMappingFunction( double fx, double fx0, double q0, bool& status ) const
{
	double value = fx/fx0-1.0;
	if( fabs(q0)>K_NEW_DOUBLE_TOL )
    {
        if( status = (q0*value+1.0 > 0.0) )
            value = log(q0*value+1.0)/q0;
    }
    else
        status = true;

	return value;
}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

