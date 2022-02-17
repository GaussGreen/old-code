//---------------------------------------------------------------------------
//
//   Group       : Convertibles Research
//
//   Filename    : FDSolver1FGeneric.hpp
//
//   Description : A finite difference solver for generic 1-factor parabolic PDEs
//
//   Author      : André Segger
//
//   Date        : 30 September 2003
//
//---------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/FDSolver1FGeneric.hpp"
#include "edginc/FDUtils.hpp"

DRLIB_BEGIN_NAMESPACE

#define	TOL         1e-2


FDSolver1FGeneric::FDSolver1FGeneric() : FDSolver1F() {}

FDSolver1FGeneric::~FDSolver1FGeneric()
{
	// Clear memory.
	deleteMem();
}


int
FDSolver1FGeneric::allocateMem(int iMax)
{

	return FDSolver1F::allocateMem(iMax);
}

void
FDSolver1FGeneric::deleteMem()
{
	FDSolver1F::deleteMem();
    return;
}

int
FDSolver1FGeneric::getGridSpots(
	double		S1Min,
	double		S1Max,
	int			iMax,
	double*		spots,	// (O) spot vector to be populated
	double*		deltax,	// (O) physical step size of the grid
	int			coordinateType,	
	double		K)
{
	return FDSolver1F::getGridSpots(S1Min,S1Max,iMax,spots,deltax,coordinateType,K);
}

int
FDSolver1FGeneric::createGrid(
	double		S1Min,
	double		S1Max,
	int			iMax,
	int			coordinateType,
	double		K)
{
	return FDSolver1F::createGrid(S1Min,S1Max,iMax,coordinateType,K);
}


int
FDSolver1FGeneric::solve(
	int				    iMax,		        // Index of the last grid point
    FDTermStructureSP   driftTerm,          // (I) Vector of drift terms in the differential equation
    FDTermStructureSP   diffusionTerm,      // (I) Vector of diffusion terms in the differential equation
    FDTermStructureSP   discountTerm,       // (I) Vector of discounting terms in the differential equation 
    FDTermStructureSP   couponTerm,         // (I) Vector of coupon terms in the differential equation 
	double			    dt,
	double*			    Vnew,
	const double*	    V,
	double			    upBarrier1,
	double			    upBC1,
	double			    downBarrier1,
	double			    downBC1,
	const double*	    payoffMin,
	const double*	    Vold,		        // Old option prices for the three time level method
	int				    iMin
    ) const
{
	return FDSolver1F::solve(iMax,driftTerm,diffusionTerm,discountTerm,couponTerm,dt,Vnew,V,upBarrier1,upBC1,downBarrier1,
		downBC1,payoffMin,Vold,iMin);
}

int
FDSolver1FGeneric::solveQuick(
	int			iMax,
	double		r,
	double		divYld1,
	double		sigma1,
	double		coupon,		// Continuous coupon
	double		dt,
	double*		V,
	double		upBarrier1,
	double		upBC1,
	double		downBarrier1,
	double		downBC1,
	int			n)	const		// The number of iterations
{
	return FDSolver1F::solveQuick(iMax,r,divYld1,sigma1,coupon,dt,V,upBarrier1,upBC1,downBarrier1,
		downBC1,n);
}

int
FDSolver1FGeneric::solveBarrier(
	int			        m,
    FDTermStructureSP   driftTerm,          // (I) Vector of drift terms in the differential equation
    FDTermStructureSP   diffusionTerm,      // (I) Vector of diffusion terms in the differential equation
    FDTermStructureSP   discountTerm,       // (I) Vector of discounting terms in the differential equation 
    FDTermStructureSP   couponTerm,         // (I) Vector of coupon terms in the differential equation 
	double		        dt,
	double*		        Vnew,
	const double*		V,
	double		        upBarrierNew,
	double		        upPayoutNew,
	double		        upPayoutDeltaNew,
	double		        upBarrier,
	double		        upPayout,
	double		        upPayoutDelta,
	double		        downBarrierNew,
	double		        downPayoutNew,
	double		        downPayoutDeltaNew,
	double		        downBarrier,
	double		        downPayout,
	double		        downPayoutDelta,
	const double*		payoffMin,
	const double*		Vold,		// Old option prices for the three time level method
	double		        upBarrierOld,
	double		        upPayoutOld,
	double		        upPayoutDeltaOld,
	double		        downBarrierOld,
	double		        downPayoutOld,
	double		        downPayoutDeltaOld
   ) const
{
	return FDSolver1F::solveBarrier(m,driftTerm,diffusionTerm,discountTerm,couponTerm,dt,Vnew,V,upBarrierNew,
        upPayoutNew,upPayoutDeltaNew,
		upBarrier,upPayout,upPayoutDelta,downBarrierNew,downPayoutNew,downPayoutDeltaNew,downBarrier,downPayout,
		downPayoutDelta,payoffMin,Vold,upBarrierOld,upPayoutOld,upPayoutDeltaOld,downBarrierOld,downPayoutOld,
		downPayoutDeltaOld);
}


int 
FDSolver1FGeneric::interpolate(
	int				nOld,		// (I) size of the table to interpolate from
	const double*	xOld,		// (I) x values of the table to interpolate from
	const double*	yOld,		// (I) y values of the table to interpolate from
	int				n,			// (I) size of the output array
	const double*	x,			// (I) array of x values at which to interpolate
	double*			y,			// (O) array of interpolated values for x
	double*			y1,			// (O) array of first derivatives (don't calculate if NULL)
	double*			y2			// (O) array of second derivatives (don't calculate if NULL)
)
{
	double*	y1Temp=NULL;
	double* y2Temp=NULL;
	int rc=SUCCESS;

	try {

	if (y1==NULL && y2 == NULL)
		rc = FDInterpolation1F(nOld,(double*)xOld,(double*)yOld,n,(double*)x,y);
	else
	{
		double* y1Ptr;
		if (y1==NULL)
		{
			y1Temp = (double*) calloc(n, sizeof(double));
			y1Ptr = y1Temp;
		}
		else
			y1Ptr = y1;
		double* y2Ptr;
		if (y2==NULL)
		{
			y2Temp = (double*) calloc(n, sizeof(double));
			y2Ptr = y2Temp;
		}
		else
			y2Ptr = y2;
		
		rc = FDInterpolationD1F(nOld,(double*)xOld,(double*)yOld,n,(double*)x,y,y1Ptr,y2Ptr);
	}

	} catch(...) { rc = FAILURE; }

	free(y1Temp);
	free(y2Temp);
	return rc;
}

int	
FDSolver1FGeneric::gridInterpolate(
	int		        mOld,				// (I) Index of the last grid point in the old grid
	const double*	Sold,		        // (I) Spot vector in the old grid 
	const double*	Vold,		        // (I) Vector of option values in the old grid 
	int		        m,					// (I) Index of the last grid point in the new grid
	const double*	S,			        // (I) Spot vector in the new grid 
	double*	        V,					// (O) vector of option values
	double	        upBarrier,			// (I) Up barrier (<=0 - no barrier)
	double	        upPayout,			// (I) Payout coefficient at and above the up barrier
	double	        upPayoutDelta,		// (I) Payout delta coeffient at and above the up barrier
	double	        downBarrier,		// (I) Down barrier (<=0 - no barrier)
	double	        downPayout,			// (I) Payout coefficient at and below the down barrier
	double	        downPayoutDelta)	// (I) Payout delta coeffient at and below the down barrier
{
	return GridInterpolate1F(mOld,(double*)Sold,(double*)Vold,m,(double*)S,V,upBarrier,upPayout,upPayoutDelta,
	                         downBarrier,downPayout,downPayoutDelta);
}


DRLIB_END_NAMESPACE
