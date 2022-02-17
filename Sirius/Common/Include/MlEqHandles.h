//	MlEqHandles.h :			 Contains all the MlEq[...]Handle typedefs
//
//	Author :				 David Cuin
/////////////////////////////////////////////////////////////////////////////

#ifndef _MLEQHANDLES_H_
#define _MLEQHANDLES_H_

template<class T> class RCPtr;
template<class T> class RCConstPtr;


/////////////////////////////////////////////////////////////////////////////
//	Strike stuff
//
class MlEqStrike;
typedef RCPtr<MlEqStrike>				MlEqStrikeHandle;
typedef RCConstPtr<MlEqStrike>		    MlEqConstStrikeHandle;

class MlEqStrikes;
typedef RCPtr<MlEqStrikes>				MlEqStrikesHandle;
typedef RCConstPtr<MlEqStrikes>			MlEqConstStrikesHandle;

/////////////////////////////////////////////////////////////////////////////
//	VolatilityStructure stuff
//
class MlEqVolData;
typedef RCPtr<MlEqVolData>				MlEqVolDataHandle;

class MlEqJumpWingVolData;
typedef RCPtr<MlEqJumpWingVolData>		MlEqJumpWingVolDataHandle;

class MlEqJumpWingFitVolData;
typedef RCPtr<MlEqJumpWingFitVolData>	MlEqJumpWingFitVolDataHandle;

class SvlFit;
typedef RCPtr<SvlFit>					SvlFitHandle;

class MlSvlVolData;
typedef RCPtr<MlSvlVolData>				MlSvlVolDataHandle;

class MLHullVolData;
typedef RCPtr<MLHullVolData>			MLHullVolDataHandle;

class MLRamVolData;
typedef RCPtr<MLRamVolData>				MLRamVolDataHandle;

class MLHermiteVolData;
typedef RCPtr<MLHermiteVolData>			MLHermiteVolDataHandle;

class RamFit;
typedef RCPtr<RamFit>					RamFitHandle;

class MlARPropFit;
typedef RCPtr<MlARPropFit>				MlARPropFitHandle;

class MlARPropVolData;
typedef RCPtr<MlARPropVolData>			MlARPropVolDataHandle;

class MlEqVolMultiplierData;
typedef RCPtr<MlEqVolMultiplierData>	MlEqVolMultiplierDataHandle;

class MlEqVolHermiteData;
typedef RCPtr<MlEqVolHermiteData>		MlEqVolHermiteDataHandle;


/////////////////////////////////////////////////////////////////////////////
// Other object typedefs
//
class RCObject;
typedef RCObject				MlEqObject;
typedef RCPtr<RCObject>			MlEqObjectHandle;

class CForwardSkewMC;
typedef RCPtr<CForwardSkewMC>			ForwardSkewMCHandle;

class MlEqAsset;
typedef RCPtr<MlEqAsset>				MlEqAssetHandle;
typedef RCConstPtr<MlEqAsset>		    MlEqConstAssetHandle;

class MlEqDate;
typedef RCPtr<MlEqDate>					MlEqDateHandle;
typedef RCConstPtr<MlEqDate>			MlEqConstDateHandle;

class MlEqCorrelationMatrix;
typedef RCPtr<MlEqCorrelationMatrix>	MlEqCorrelationMatrixHandle;

class MlEqDividendSchedule;
typedef RCPtr<MlEqDividendSchedule>		MlEqDividendScheduleHandle;

class MlEqSpotSchedule;
typedef RCPtr<MlEqSpotSchedule>			MlEqSpotScheduleHandle;

class MlEqZeroCurve;
typedef RCPtr<MlEqZeroCurve>			MlEqZeroCurveHandle;	

class MlEqSwap;
typedef RCPtr<MlEqSwap>					MlEqSwapHandle;

class MlEqVolatilityStructure;
typedef RCPtr<MlEqVolatilityStructure>	MlEqVolatilityStructureHandle;

class MlEqStochBetaVolatilityStructure;
typedef RCPtr<MlEqStochBetaVolatilityStructure> MlEqStochBetaVolatilityStructureHandle;

class MlEqInterpolator;
typedef RCPtr<MlEqInterpolator>			MlEqInterpolatorHandle;
typedef RCConstPtr<MlEqInterpolator>	MlEqConstInterpolatorHandle;

class MlEq2DInterpolator;
typedef RCPtr<MlEq2DInterpolator>		MlEq2DInterpolatorHandle;

class randomGenerator;
typedef RCPtr<randomGenerator>			randomGeneratorHandle;

class MlEqHullAndWhite;
typedef RCPtr<MlEqHullAndWhite>			MlEqHullAndWhiteHandle;

class MlEqDateSchedule;
typedef RCPtr<MlEqDateSchedule>			MlEqDateScheduleHandle;

class MlEqMonteCarlo;
typedef RCPtr<MlEqMonteCarlo>			MlEqMonteCarloHandle;

class MlEqParameterList;
typedef RCPtr<MlEqParameterList>		MlEqParameterListHandle;

class MlEqArray;
typedef RCPtr<MlEqArray>				MlEqArrayHandle;

class MlEqMatrix;
typedef RCPtr<MlEqMatrix>				MlEqMatrixHandle;

class MlEqScenario;
typedef RCPtr<MlEqScenario>				MlEqScenarioHandle;

class MlEqPdeHelper;
typedef RCPtr<MlEqPdeHelper>			MlEqPdeHelperHandle;

class MlEqPdeDriver;
typedef RCPtr<MlEqPdeDriver>			MlEqPdeDriverHandle;

class knockOutBarrierHelper;
typedef RCPtr<knockOutBarrierHelper>	knockOutBarrierHelperHandle;

class DupireLocalVol;
typedef RCPtr<DupireLocalVol>	DupireLocalVolHandle;

class MlEqAnalyticCurveWithTanhWingEdge;
typedef RCPtr<MlEqAnalyticCurveWithTanhWingEdge>	MlEqAnalyticCurveWithTanhWingEdgeHandle;

class pdeLocalVolEffective;
typedef RCPtr<pdeLocalVolEffective>	pdeLocalVolEffectiveHandle;


#endif