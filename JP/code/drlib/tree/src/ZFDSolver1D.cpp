// ===================================================================
// PROGRAM:		ZFDSolver1D.cpp
//
// AUTHOR:		Milan Kovacevic April 04, 2001.
//
// CONTAINS:	FD 2D Solver class
//
//=================================================================

//#include <zmarket.h>
//#include "drerrors.h"

//#include "mymalloc.h"
//#include "macros.h"
//#include <zerror.h>
#include "edginc/config.hpp" // needed for precompiled headers
#include "edginc/ZFDSolver1D.hpp"
using namespace std;

DRLIB_BEGIN_NAMESPACE

int ZFDSolver1D::triDiag1D(
	double*		a,
	double*		b,
	double*		c,
	double*		r,
	double*		u,
	int			n)
{
    int	status = FAILURE;
    double bet;
//	int j;
//	double*	gam=NULL;

//    gam = (double *)CALLOC(n + 1, sizeof(double));

//	if (b[0] == 0.0) goto general_failure;
/*
  u[0]=r[0]/(bet=b[0]);
  for (j=1;j<=n;j++) {
  gam[j]=c[j-1]/bet;
  bet=b[j]-a[j]*gam[j];
//		if (bet == 0.0)	goto general_failure;
u[j]=(r[j]-a[j]*u[j-1])/bet;
}
for (j=(n-1);j>=0;j--)
u[j] -= gam[j+1]*u[j+1];
*/
    double* priceBottomPtr;
    double* pricePtr = priceBottomPtr = &(u[0]);
    double* priceTopPtr = &(u[n]);
    double* rPtr = &(r[0]);
    double* aPtr = &(a[1]);
    double* bPtr = &(b[0]);
    double* cPtr = &(c[0]);
    double* gamPtr = &gam[1];
    *pricePtr = (*rPtr++)/(bet=(*bPtr++));
    while (pricePtr < priceTopPtr)
    {
        ++pricePtr;
        *gamPtr = (*cPtr++)/bet;
        bet = *bPtr++ - (*aPtr)*(*gamPtr++);
        *pricePtr = ((*rPtr++) - (*aPtr++)*(*(pricePtr-1)))/bet;
    }
    --gamPtr;

    while (pricePtr > priceBottomPtr)
    {
        *(pricePtr-1) -= (*gamPtr--)*(*pricePtr);
        --pricePtr;
    }

    status  = SUCCESS;
//general_failure:
//	::free(gam);
    return status;
}

/*
  int
  ZFDSolver1D::sor(
  double*		a,
  double*		b,
  double*		c,
  double*		d,
  double*		payoff,
  double*		u,
  int			n,
  double		tolerance)
  {
  int i,iter,maxIter;
  double newvar,normf,normg,normr,omega2,sizer;
  double		normSA;

  maxIter = 100;
  omega2 = 1.;

  iter = 1;
  if (tolerance<1e-10)
  {
//        DrError("Linear tolerance is too small.  tolerance < 1e-10.");
return FAILURE;
}

normSA = 0.;
for (i=0; i<=n; i++) {
normSA = normSA + a[i]*a[i] + b[i]*b[i] + c[i]*c[i];
}
normSA = sqrt(normSA / (double)(n + 1)) / sqrt(3.);

normg = 0.;
for (i=0; i<=n; i++) {
normg = normg + d[i]*d[i];
}
normg = sqrt(normg / (double)(n + 1));
while (1) {
gam[0] = d[0] - u[0]*b[0] - u[1]*c[0];
newvar = u[0] + (omega2*gam[0]) / b[0];

newvar = max(newvar,payoff[0]);
gam[0] = (newvar - u[0])*b[0];
u[0] = newvar;
for (i=1; i<=n - 1; i++)
{
gam[i] = d[i] - u[i - 1]*a[i] - u[i]*b[i] - u[i + 1]*c[i];
newvar = u[i] + (omega2*gam[i]) / b[i];

newvar = max(newvar,payoff[i]);
gam[i] = (newvar - u[i])*b[i];
u[i] = newvar;
}
gam[n] = d[n] - u[n - 1]*a[n] - u[n]*b[n];
newvar = u[n] + (omega2*gam[n]) / b[n];
newvar = max(newvar,payoff[n]);
gam[n] = (newvar - u[n])*b[n];
u[n] = newvar;
normr = 0.;
for (i=0; i<=n; i++) {
normr = normr + gam[i]*gam[i];
}
normr = sqrt(normr / (double)(n + 1));
normf = 0.;
for (i=0; i<=n; i++) {
normf = normf + u[i]*u[i];
}
normf = sqrt(normf / (double)(n + 1));
sizer = normg + normf*normSA;
if ((normr>sizer*tolerance&&iter<maxIter))
{
    iter++;
}
else
{
    if ((iter>=maxIter&&maxIter!=5))
    {
//                DrError"loop terminating because maximum number of iterations reached");
        return FAILURE;
    }
    return SUCCESS;
}
}
return SUCCESS;
}
*/

int
ZFDSolver1D::sor(
    double*		a,
    double*		b,
    double*		c,
    double*		d,
    double*		payoff,
    double*		u,
    int			n,
    double		tolerance)
{
//	double		tolerance = CHEM_EPSILON9;
    double		temp;
    double		norm;
    double		iter;
    double		maxiter=100;
    double		w=1.;
/*
  int			i;
  for (iter = 0; iter < maxiter; iter++)
  {
  norm = 0;
  for (i=1; i <=n-1; i++)
  {
  // -a[i]*u[i-1] + b[i]*u[i] - c[i]*u[i+1] = d[i]
//			temp = u[i] + w * (d[i] + a[i] * u[i - 1] - b[i] * u[i] + c[i] * u[i + 1]) / b[i];
temp = u[i] + w * (d[i] - a[i] * u[i - 1] - b[i] * u[i] - c[i] * u[i + 1]) / b[i];
//			temp = y - u[i];
//			norm = norm + temp * temp;
//			u[i] = u[i] + w * temp;
//			if (u[i] < payoff[i])
//				u[i] = payoff[i];
if (temp < payoff[i])
temp = payoff[i];
norm += (temp-u[i])*(temp-u[i]);
u[i] = temp;
}
iter++;
if (norm <= tolerance)
return SUCCESS;
}
return FAILURE;
*/

    double*	uPtr;
    double*	uBottomPtr = &(u[0]);
    double* uTopPtr = &(u[n]);
    double*	aPtr;
    double*	aBottomPtr = &(a[0]);
    double*	bPtr;
    double*	bBottomPtr = &(b[0]);
    double*	cPtr;
    double*	cBottomPtr = &(c[0]);
    double*	dPtr;
    double*	dBottomPtr = &(d[0]);
    double*	payoffPtr;
    double*	payoffBottomPtr=&(payoff[0]);
    double	temp1;

    for (iter = 0; iter < maxiter; iter++)
    {
        norm = 0;
        uPtr = uBottomPtr;
        aPtr = aBottomPtr;
        bPtr = bBottomPtr;
        cPtr = cBottomPtr;
        dPtr = dBottomPtr;
        payoffPtr = payoffBottomPtr;
        while (uPtr <= uTopPtr)
        {
            temp = *uPtr + w * ((*dPtr++) - (*aPtr++) * (*(uPtr-1)) - (*bPtr) * (*uPtr) - (*cPtr++) * (*(uPtr+1))) / (*bPtr);
            bPtr++;
            if (temp < *payoffPtr)
                temp = *payoffPtr;
            payoffPtr++;
            temp1 = temp - (*uPtr);
            norm += (temp1)*(temp1);
            *uPtr++ = temp;
        }
        iter++;
        if (norm <= tolerance)
            return SUCCESS;
    }
    return FAILURE;


}


int
ZFDSolver1D::sorMinMax(
	double*		a,
	double*		b,
	double*		c,
	double*		d,
	double*		payoffMin,
	double		payoffMax,
	double*		u,
	int			n,
	double		tolerance)
{
//	double		tolerance = CHEM_EPSILON9;
    double		temp;
    double		norm;
    double		iter;
    double		maxiter=100;
    double		w=1.;
/*
  int			i;
  for (iter = 0; iter < maxiter; iter++)
  {
  norm = 0;
  for (i=1; i <=n-1; i++)
  {
  // -a[i]*u[i-1] + b[i]*u[i] - c[i]*u[i+1] = d[i]
//			temp = u[i] + w * (d[i] + a[i] * u[i - 1] - b[i] * u[i] + c[i] * u[i + 1]) / b[i];
temp = u[i] + w * (d[i] - a[i] * u[i - 1] - b[i] * u[i] - c[i] * u[i + 1]) / b[i];
//			temp = y - u[i];
//			norm = norm + temp * temp;
//			u[i] = u[i] + w * temp;
//			if (u[i] < payoff[i])
//				u[i] = payoff[i];
if (temp < payoff[i])
temp = payoff[i];
norm += (temp-u[i])*(temp-u[i]);
u[i] = temp;
}
iter++;
if (norm <= tolerance)
return SUCCESS;
}
return FAILURE;
*/

    double*	uPtr;
    double*	uBottomPtr = &(u[0]);
    double* uTopPtr = &(u[n]);
    double*	aPtr;
    double*	aBottomPtr = &(a[0]);
    double*	bPtr;
    double*	bBottomPtr = &(b[0]);
    double*	cPtr;
    double*	cBottomPtr = &(c[0]);
    double*	dPtr;
    double*	dBottomPtr = &(d[0]);
    double*	payoffMinPtr;
    double*	payoffMinBottomPtr=&(payoffMin[0]);
    double	temp1;

    for (iter = 0; iter < maxiter; iter++)
    {
        norm = 0;
        uPtr = uBottomPtr;
        aPtr = aBottomPtr;
        bPtr = bBottomPtr;
        cPtr = cBottomPtr;
        dPtr = dBottomPtr;
        payoffMinPtr = payoffMinBottomPtr;
        while (uPtr <= uTopPtr)
        {
            temp = *uPtr + w * ((*dPtr++) - (*aPtr++) * (*(uPtr-1)) - (*bPtr) * (*uPtr) - (*cPtr++) * (*(uPtr+1))) / (*bPtr);
            bPtr++;
            if (payoffMax < *payoffMinPtr)
                return FAILURE;
            if (temp < *payoffMinPtr)
                temp = *payoffMinPtr;
            if (temp > payoffMax)
                temp = payoffMax;
            payoffMinPtr++;
            temp1 = temp - (*uPtr);
            norm += (temp1)*(temp1);
            *uPtr++ = temp;
        }
        iter++;
        if (norm <= tolerance)
            return SUCCESS;
    }
    return FAILURE;


}


ZFDSolver1D::ZFDSolver1D()
{
    // Clear memory.
    resetMem();
}

ZFDSolver1D::ZFDSolver1D(
	const ZFDSolver1D&	mh)
{
    // Initially we have no pointers
    resetMem();

    // Leverage the assignment operator!
    *this = mh;
}

ZFDSolver1D::~ZFDSolver1D()
{
    // Clear memory.
    deleteMem();
}

#if 0
// likely to core dump if ever used as the allocated memory is not copied
// the first one to be deleted corrupts any assignments
ZFDSolver1D&
ZFDSolver1D::operator=(const ZFDSolver1D& source)
{

    if(this == &source)
        return *this;

/*
  m_rhs = source.m_rhs;
  m_rhs1 = source.m_rhs1;
*/
    return *this;
}
#endif

void
ZFDSolver1D::resetMem()
{
    m_iMax=0;
    J1 = NULL;
    JdriftTerm = 0;
    JvolTermDown = 0;
    JvolTermUp = 0;
    m_rhs = 0;
    spot1 = 0;
    left1 = 0;
    diag1 = 0;
    right1 = 0;
    gam = 0;
    m_nonNegativeValues = TRUE;
}

int ZFDSolver1D::allocateMem(int iMax)
{
    static const string method = "ZFDSolver1D::allocateMem";
    int rc=SUCCESS;

    try
    {

        m_iMaxMem = iMax;

/*        m_rhs = (double *)CALLOC(iMax + 1, sizeof(double));
          J1 = (double *)CALLOC(iMax + 1, sizeof(double));
          JdriftTerm = (double *)CALLOC(iMax + 1, sizeof(double));
          JvolTermDown = (double *)CALLOC(iMax + 1, sizeof(double));
          JvolTermUp = (double *)CALLOC(iMax + 1, sizeof(double));

          spot1 = (double *)CALLOC(iMax + 1, sizeof(double));
          left1 = (double *)CALLOC(iMax + 1, sizeof(double));
          diag1 = (double *)CALLOC(iMax + 1, sizeof(double));
          right1 = (double *)CALLOC(iMax + 1, sizeof(double));
          gam = (double *)CALLOC(iMax + 1, sizeof(double));*/

        m_rhs = (double *)calloc(iMax + 1, sizeof(double));
        J1 = (double *)calloc(iMax + 1, sizeof(double));
        JdriftTerm = (double *)calloc(iMax + 1, sizeof(double));
        JvolTermDown = (double *)calloc(iMax + 1, sizeof(double));
        JvolTermUp = (double *)calloc(iMax + 1, sizeof(double));

        spot1 = (double *)calloc(iMax + 1, sizeof(double));
        left1 = (double *)calloc(iMax + 1, sizeof(double));
        diag1 = (double *)calloc(iMax + 1, sizeof(double));
        right1 = (double *)calloc(iMax + 1, sizeof(double));
        gam = (double *)calloc(iMax + 1, sizeof(double));


    } catch (exception& e){
        throw ModelException(&e, method);
    }


    if (rc != SUCCESS)
        deleteMem();

    return rc;
}

void
ZFDSolver1D::deleteMem()
{
    ::free(m_rhs);
    ::free(J1);
    ::free(JdriftTerm);
    ::free(JvolTermDown);
    ::free(JvolTermUp);
    ::free(spot1);
    ::free(left1);
    ::free(diag1);
    ::free(right1);
    ::free(gam);

    resetMem();
    return;

}

int
ZFDSolver1D::getGridSpots(
	double		S1Min,  // min stock price
	double		S1Max,  // max stock price
	int			iMax,   // no. stock intervals
	double*		spots,	// (O) spot vector to be populated
	double*		deltax,	// (O) physical step size of the grid
	int			coordinateType, // 0-sinh, 1-linear, 2-exponential
	double		K)      // strike - used for sinh grid type only
{
    double	temp;
    double	xMin,xMax;
    double	dx;
    int		i;

    switch (coordinateType)
    {

    case	0:	// sinh
        temp = S1Min/K;
        xMin = log(temp+sqrt(1+temp*temp));
        temp = S1Max/K;
        xMax = log(temp+sqrt(1+temp*temp));
        break;

    case	1:	// linear
        xMin = S1Min;
        xMax = S1Max;
        break;

    case	2:	// exponential
        xMax = log(S1Max);
        xMin = log(S1Min);
        break;

    default:
        return FAILURE;
    }

    dx = (xMax - xMin) / (double)(iMax);

    // Set upper and lower bound explicitly to avoid platform dependent rounding errors
    spots[0] = S1Min;
    spots[iMax] = S1Max;

    // Initialize coordinate arrays for S1
    switch (coordinateType)
    {

    case	0:	// sinh
        for (i=1; i<iMax; i++)
            spots[i] = K*sinh(dx*i + xMin);
        break;

    case	1:	// linear
        for (i=1; i<iMax; i++)
            spots[i] = dx*i + xMin;
        break;

    case	2:	// exponential
        for (i=1; i<iMax; i++)
            spots[i] = exp(dx*i + xMin);
        break;

    default:
        return FAILURE;
    }

    *deltax = dx;

    return SUCCESS;
}

int
ZFDSolver1D::createGrid(
	double		S1Min,
	double		S1Max,
	int			iMax,
	int			coordinateType,
	double		K)
{
    int		i;

    if (iMax > m_iMaxMem)
        return FAILURE;

    m_iMax = iMax;

    if (ZFDSolver1D::getGridSpots(S1Min,S1Max,iMax,spot1,&dx,coordinateType,K) != SUCCESS)
        return FAILURE;

    // Initialize step size
    //   ds1 = (s1Max - s1Min) / (double)(iMax);

    for (i=0; i<iMax; i++)
    {
        J1[i] = (spot1[i + 1] - spot1[i]) / dx;
    }
    // Precompute for performance
    JdriftTerm[0] = 1. / (J1[0]*dx);
    JdriftTerm[iMax] = 1. / (J1[iMax-1]*dx);
    for (i=1; i<iMax; i++)
    {
        JdriftTerm[i] = 1. / ((J1[i-1] + J1[i])*dx);
//		JvolTermDown[i] = 1./ ((J1[i-1]+J1[i])*J1[i-1]*dx*dx);
//		JvolTermUp[i] = 1./ ((J1[i-1]+J1[i])*J1[i]*dx*dx);
        JvolTermDown[i] = JdriftTerm[i]/(J1[i-1]*dx);
        JvolTermUp[i] = JdriftTerm[i]/(J1[i]*dx);
    }
    return SUCCESS;
}


int
ZFDSolver1D::solve(
	double		upBarrier1,
	double		upBC1,
	double		downBarrier1,
	double		downBC1,
	double		divYld1,
	int			iMax,
    int         iMin,
	double		r,
	double*		sigma1,
	double*		cs,
	double		coupon,		// Continuous coupon
	double		dt,
	double*		V,
	double*		Vnew,
	int			american,
	double*		payoffMin,
	double		toleranceSOR,
	double		payoffMax,
	double		theta,		// Discretization method: 0 - explicit, 0.5 - Crank - Nicholson, 1 - implicit, 2 - three time level
	double*		Vold,		// Old option prices for the three time level method
    int         varMethod
    )
{
    int i;

    double	discountTerm;
    double	driftTerm1,volTerm1;
    double	couponTerm;
    double	muTerm;
//	double	ds1;
    double	alphaUp=1.;
    double	alphaDown=1.;
//	double	theta=0.5;
    double	pup1,pdown1;
    int		hasUpBarrier1 = FALSE;
    int		hasDownBarrier1 = FALSE;
    int		threeTimeLevel=FALSE;

    double	temp;
    double	tempa;
    double  oTemp;
//	int		rc=0;
//	CT_TRY

    if (iMax > m_iMax)
        return FAILURE;

    if (theta > 1)
    {
        threeTimeLevel=TRUE;
        theta = 1;	// Do implicit
        if (Vold == NULL)
            return FAILURE;
    }

    if (upBarrier1 > 0)
    {
        if (DBL_EQUAL(upBarrier1,spot1[iMax]))
            alphaUp = 1;
        else if (upBarrier1 > spot1[iMax] || upBarrier1 <= spot1[iMax-1])
            return FAILURE;
        else
            // The boundary does not fall on a grid point (it is between nodes iMax-1 and iMax)
            alphaUp = (spot1[iMax]-spot1[iMax-1])/(upBarrier1-spot1[iMax-1]);
        hasUpBarrier1 = TRUE;
    }

    if (downBarrier1 > 0)
    {
        if (DBL_EQUAL(downBarrier1,spot1[0]))
            alphaDown = 1;
        else if (downBarrier1 < spot1[iMin] || downBarrier1 >= spot1[iMin+1])
            return FAILURE;
        else
            // The boundary does not fall on a grid point (it is between nodes 0 and 1)
            alphaDown = (spot1[iMin]-spot1[iMin+1])/(downBarrier1-spot1[iMin+1]);
        hasDownBarrier1 = TRUE;
    }

    // Set up the system

    if (cs)
        discountTerm = dt*(r+cs[0]);
    else
        discountTerm = dt*r;
    couponTerm = dt*coupon;	// continuous coupon (e.g., accrued coupon)
    muTerm = (r - divYld1)*dt;

    if (hasDownBarrier1)
    {
        // V[0] = downBC1*alphaDown + (1-alphaDown)*V[1]
        left1[iMin*1] = 0;
        diag1[iMin*1] = 1;
//		diag1[iMin] = alphaDown;
        right1[iMin*1] = -(1-alphaDown);
        m_rhs[iMin*1] = downBC1*alphaDown;
    }
    else
    {
        driftTerm1= muTerm*spot1[0]*JdriftTerm[iMin*1];
//		driftTerm1= dt*(r - divYld1)*spot1[0]/(J1[0]*dx);
        pup1 = driftTerm1;
        left1[iMin*1] = 0.;
        diag1[iMin*1] = 1.+ (discountTerm + pup1)*theta;
        right1[iMin*1] = - pup1*theta;
        if (threeTimeLevel)
        {
            diag1[iMin*1] += 0.5;
            m_rhs[iMin*1] = (4*V[iMin*1]-Vold[iMin*1])/2 + couponTerm;
        }
        else
        {
            m_rhs[iMin*1] = (1. - (discountTerm + pup1)*(1-theta))*V[iMin*1]
                + pup1*(1-theta)*V[iMin*1+1] + couponTerm;
        }
    }

    for (i=iMin*1+1; i<iMax; i++) {
        if (cs)
            discountTerm = dt*(r+cs[i]);
        if (varMethod <= 0) {
            temp =  sigma1[i]*spot1[i];
            tempa = dt*temp*temp;
            oTemp = 1.;
        } else if (varMethod >= 1) {
            temp = exp(sigma1[i]*sigma1[i]*dt);
            tempa = (temp-1.)*spot1[i]*spot1[i];
            oTemp = 1./temp;
        }
        driftTerm1 = muTerm*spot1[i]*JdriftTerm[i];
//		driftTerm1 = dt*spot1[i]*(r - divYld1)/((J1[i-1]+J1[i])*dx);
        volTerm1 = tempa*JvolTermDown[i];
//		volTerm1 = dt*sigma1[i]*sigma1[i]*spot1[i]*spot1[i]/((J1[i-1]+J1[i])*J1[i-1]*dx*dx);
        pdown1 = (volTerm1 - driftTerm1);
        volTerm1 = tempa*JvolTermUp[i];
//		volTerm1 = dt*sigma1[i]*sigma1[i]*spot1[i]*spot1[i]/((J1[i-1]+J1[i])*J1[i]*dx*dx);
        pup1 = (volTerm1 + driftTerm1);
        left1[i] = - oTemp*pdown1*theta;
        diag1[i] = 1.+ (discountTerm + oTemp*pdown1 + oTemp*pup1)*theta;
        right1[i] = - oTemp*pup1*theta;

        if (threeTimeLevel)
        {
            diag1[i] += 0.5;
            m_rhs[i] = (4*V[i]-Vold[i])/2 + couponTerm;
        }
        else
        {
            m_rhs[i] = (1. - (discountTerm + pdown1 + pup1)*(1-theta))*V[i]
                + pdown1*(1-theta)*V[i-1]
                + pup1*(1-theta)*V[i+1] + couponTerm;
        }
    }

    if (hasUpBarrier1)
    {
        // V[iMax] = upBC1*alphaUp + (1-alphaUp)*V[iMax-1]
        left1[iMax] = -(1-alphaUp);
        diag1[iMax] = 1;
        right1[iMax] = 0;
        m_rhs[iMax] = upBC1*alphaUp;
    }
    else
    {
        if (cs)
            discountTerm = dt*(r+cs[iMax]);
        driftTerm1= muTerm*spot1[iMax]/(J1[iMax-1]*dx);
        pdown1 = - driftTerm1;
        left1[iMax] = - pdown1*theta;
        diag1[iMax] = 1.+ (discountTerm + pdown1)*theta;
        right1[iMax] = 0;
        if (threeTimeLevel)
        {
            diag1[iMax] += 0.5;
            m_rhs[iMax] = (4*V[iMax]-Vold[iMax])/2 + couponTerm;
        }
        else
        {
            m_rhs[iMax] = (1. - (discountTerm + pdown1)*(1-theta))*V[iMax]
                + pdown1*(1-theta)*V[iMax-1] + couponTerm;
        }
    }


    // Solve using tridiagonal solvers

    triDiag1D(&left1[iMin],&diag1[iMin],&right1[iMin],&m_rhs[iMin],&Vnew[iMin],iMax-iMin);

/*
  // This is not needed since it is done implicitly inside the tridiagonal solver
  if (hasUpBarrier1)
  {
  if (alphaUp > 1)
  Vnew[iMax] = upBC1*alphaUp+Vnew[iMax-1]*(1-alphaUp);
  }

  if (hasDownBarrier1)
  {
  if (alphaDown > 1)
  Vnew[0] = downBC1*alphaDown+Vnew[1]*(1-alphaDown);
  }
*/
    if (m_nonNegativeValues)
    {
        // This can improve convergence - dampen oscillations around a barrier
        for (i=0; i<=iMax; i++)
            Vnew[i] = MAX(Vnew[i],0.);
    }
    if (american)
    {
//		for (int i=0; i<=iMax; i++)
//			Vnew[i] = MAX(Vnew[i],payoff[i]);

        if (hasUpBarrier1)
            payoffMin[iMax]=0;

        if (payoffMax>0)
            sorMinMax(left1,diag1,right1,m_rhs,payoffMin,payoffMax,Vnew,iMax,toleranceSOR);
        else
            sor(left1,diag1,right1,m_rhs,payoffMin,Vnew,iMax,toleranceSOR);

    }
    return SUCCESS;
}

int
ZFDSolver1D::solveBarrier(
	int			m,
	double		r,
	double		divYld1,
	double*		sigma1,
	double*		cs,
	double		coupon,		// Continuous coupon
	double		dt,
	double*		V,
	double*		Vnew,
	double		upBarrierNew,
	double		upPayoutNew,
	double		upPayoutDeltaNew,
	double		upBarrierOld,
	double		upPayoutOld,
	double		upPayoutDeltaOld,
	double		downBarrierNew,
	double		downPayoutNew,
	double		downPayoutDeltaNew,
	double		downBarrierOld,
	double		downPayoutOld,
	double		downPayoutDeltaOld,
	int			american,
	double*		payoffMin,
	double		toleranceSOR,
	double		payoffMax,
	double		theta,		// Discretization method: 0 - explicit, 0.5 - Crank - Nicholson, 1 - implicit, 2 - three time level
	double*		Vold,		// Old option prices for the three time level method
    int         varMethod
    )
{
	int		i;
	int		status;
	double	Voldtmp=0,Voldtmp1=0;
    double  VDoldtmp = 0, VDoldtmp1=0;
	int		iMax = m;
    int     iMin = 0;
//	int		hasUpBarrier = FALSE;
	double	upB = 0, downB = 0;
	double	alphaUp=1, alphaDown=1;
	int		numSaved=0, numDSaved=0;
	double	upBCNew=0;
	double	downBCNew=0;
//	int		rc=0;
//	CT_TRY


    if (r != 0. && divYld1 != 0. && dt == 0.) {
        double* stockFwd = new double [m+1];
        // r and divYld1 have different definitions when dt == 0
        // r is the ir forward factor
        // divYld1 is the ir forward factor divided by the equity forward factor
        for (i=m; i>= 0; i--) {
		    stockFwd[i] = spot1[i]*r/divYld1;
        }

        if (FDInterpolation(m+1, spot1, V, m+1, stockFwd, Vnew) == FAILURE) {
            throw ModelException("ZFDSolver1D::solveBarrier", "FDInterpolation failure");
        }

        for (i=m; i>= 0; i--) {
		    Vnew[i] /= r;
        }

        delete [] stockFwd;
        return SUCCESS;
    }

	iMax = m;
	if (upBarrierNew > 0)
	{
		// Find iMax for Vnew
		upB = upBarrierNew*0.9999999999;

        if (spot1[0] >= upB) { // all points are above the barrier
            for (i=m; i>= 0; i--) {
			    Vnew[i] = upPayoutNew+upPayoutDeltaNew*spot1[i];;
            }
            return SUCCESS;
        }

        while (spot1[iMax] >= upB)
			iMax--;

		if (iMax < m)
			iMax++;
		else
			upBarrierNew = -1; // spot1[m] < upB so the barrier will be ignored
	}

	if (upBarrierOld > 0)
	{
		// If the old  barrier is greater or equal than iMax, ignore it, otherwise
		// adjust option prices above the barrier up to spot1[iMax]
		if (spot1[iMax] >= upBarrierOld*0.9999999999)
		{
			double upBCOld = upPayoutOld + upPayoutDeltaOld*upBarrierOld;
			Voldtmp = V[iMax];	// Save the old value
			if (upBarrierOld >= spot1[iMax-1])
			{
				alphaUp = (spot1[iMax]-spot1[iMax-1])/(upBarrierOld-spot1[iMax-1]);
				V[iMax] = upBCOld*alphaUp + (1-alphaUp)*V[iMax-1];
				numSaved = 1;
			}
			else	// The barrier is below spot1[iMax-1]
			{
                if ( iMax <= 1) {
                    throw ModelException("ZFDSolver1D::solveBarrier",
                        "The barrier is very close to the lowest spot price, which can cause unstable numerical solutions. Please update the model paramters");
                }

				Voldtmp1 = V[iMax-1];	// Save the old value
				alphaUp = (spot1[iMax]-spot1[iMax-2])/(upBarrierOld-spot1[iMax-2]);
				V[iMax] = upBCOld*alphaUp + (1-alphaUp)*V[iMax-2];
				alphaUp = (spot1[iMax-1]-spot1[iMax-2])/(upBarrierOld-spot1[iMax-2]);
				V[iMax-1] = upBCOld*alphaUp + (1-alphaUp)*V[iMax-2];
				numSaved = 2;
			}
		}
	}

	if (downBarrierNew > 0)
	{
		// Find iMin for Vnew
		downB = downBarrierNew*1.00000001;

        if (spot1[iMax] <= downB) { // all points are below the barrier
            for (i=0; i<= iMax; i++) {
			    Vnew[i] = downPayoutNew + downPayoutDeltaNew*downBarrierNew;
            }
            return SUCCESS;
        }

        while (spot1[iMin] <= downB)
			iMin++;

		if (iMin > 0)
			iMin--;
		else
			downBarrierNew = -1; // spot1[0] > downB so the barrier will be ignored
	}

	if (downBarrierOld > 0)
	{
		// If the old  barrier is smaller or equal than iMin, ignore it, otherwise
		// adjust option prices above the barrier up to spot1[iMax]
		if (spot1[iMin] <= downBarrierOld*1.0000001)
		{
			double downBCOld = downPayoutOld + downPayoutDeltaOld*downBarrierOld;
			VDoldtmp = V[iMin];	// Save the old value
			if (downBarrierOld <= spot1[iMin+1])
			{
				alphaDown = (spot1[iMin]-spot1[iMin+1])/(downBarrierOld-spot1[iMin+1]);
				V[iMin] = downBCOld*alphaDown + (1-alphaDown)*V[iMin+1];
				numDSaved = 1;
			}
			else	// The barrier is below spot1[iMax-1]
			{
                if ( iMin >= m-1) {
                    throw ModelException("ZFDSolver1D::solveBarrier",
                        "The down barrier is very close to the highest spot price, which can cause unstable numerical solutions. Please update the model paramters");
                }

				VDoldtmp1 = V[iMin+1];	// Save the old value
				alphaDown = (spot1[iMin]-spot1[iMin+2])/(downBarrierOld-spot1[iMin+2]);
				V[iMin] = downBCOld*alphaDown + (1-alphaDown)*V[iMin+2];
				alphaDown = (spot1[iMin+1]-spot1[iMin+2])/(downBarrierOld-spot1[iMin+2]);
				V[iMin+1] = downBCOld*alphaDown + (1-alphaDown)*V[iMin+2];
				numDSaved = 2;
			}
		}
	}

	if (upBarrierNew > 0)
		upBCNew = upPayoutNew + upPayoutDeltaNew*upBarrierNew;
	if (downBarrierNew > 0)
		downBCNew = downPayoutNew + downPayoutDeltaNew*downBarrierNew;
	status = solve (upBarrierNew,upBCNew,downBarrierNew,downBCNew,
                    divYld1,iMax,iMin,r,sigma1,cs,coupon,dt,V,Vnew,american,payoffMin,toleranceSOR,payoffMax,theta,Vold,varMethod);

	if (status != SUCCESS)
		return status;

	// Restore old option values above the boundary
	if (numSaved > 0)
	{
		V[iMax] = Voldtmp;
		if (numSaved > 1)
			V[iMax-1] = Voldtmp1;
	}
	if (numDSaved > 0)
	{
		V[iMin] = VDoldtmp;
		if (numDSaved > 1)
			V[iMin+1] = VDoldtmp1;
	}

	if (upBarrierNew > 0)
	{
		// Set new options values which are above the boundary
		for (i=m; i>= iMax; i--)
			Vnew[i] = upPayoutNew+upPayoutDeltaNew*spot1[i];;
//			Vnew[i] = intrinsic[i];
	}

	if (downBarrierNew > 0)
	{
		// Set new options values which are above the boundary
		for (i=0; i<= iMin; i++)
			Vnew[i] = downPayoutNew+downPayoutDeltaNew*spot1[i];;
	}

/*
  CT_CATCH_RC
  if (rc == FAILURE)
  return rc;
*/
	return status;
}


DRLIB_END_NAMESPACE
