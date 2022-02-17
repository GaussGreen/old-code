//----------------------------------------------------------------------------
//
//   Group       : Convertibles Research
//
//   Filename    : FDSolver1F.cpp
//
//   Description : A base class for 1-factor FD solvers
//
//   Author      : André Segger
//
//   Date        : 13 October 2003
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/FDSolver1F.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

#define	TOL         1e-2
#define	FAILURE     -1
#define	SUCCESS     0
#define	EPSILON     1e-10


int FDSolver1F::triDiag1D(
	const double*		a, 
	const double*		b, 
	const double*		c, 
	const double*		r, 
	double*				u,
	int					n) const {
	int	status = FAILURE;
	double bet;
	double* priceBottomPtr;
	double* pricePtr = priceBottomPtr = &(u[0]);
	double* priceTopPtr = &(u[n]);
	double* rPtr = (double*) &(r[0]);
	double* aPtr = (double*) &(a[1]);
	double* bPtr = (double*) &(b[0]);
	double* cPtr = (double*) &(c[0]);
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
	return status;
}


int FDSolver1F::sor(
	const double*		a,
	const double*		b,
	const double*		c,
	const double*		d,
	const double*		payoff,
	double*				u,
	int					n,
	double				tolerance) {
	double		temp;
	double		norm;
	double		iter;
	double		maxiter=100;
	double		w=1.;

	double*	uPtr;
	double*	uBottomPtr = &(u[0]);
	double* uTopPtr = &(u[n]);
	double*	aPtr;
	double*	aBottomPtr = (double*) &(a[0]);
	double*	bPtr;
	double*	bBottomPtr = (double*) &(b[0]);
	double*	cPtr;
	double*	cBottomPtr = (double*) &(c[0]);
	double*	dPtr;
	double*	dBottomPtr = (double*) &(d[0]);
	double*	payoffPtr;
	double*	payoffBottomPtr= (double*) &(payoff[0]);
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


int FDSolver1F::sorMinMax(
	const double*		a,
	const double*		b,
	const double*		c,
	const double*		d,
	const double*		payoffMin,
	double		payoffMax,
	double*		u,
	int			n,
	double		tolerance) {
	double		temp;
	double		norm;
	double		iter;
	double		maxiter=100;
	double		w=1.;

	double*	uPtr;
	double*	uBottomPtr = &(u[0]);
	double* uTopPtr = &(u[n]);
	double*	aPtr;
	double*	aBottomPtr = (double* ) &(a[0]);
	double*	bPtr;
	double*	bBottomPtr = (double* ) &(b[0]);
	double*	cPtr;
	double*	cBottomPtr = (double* ) &(c[0]);
	double*	dPtr;
	double*	dBottomPtr = (double* ) &(d[0]);
	double*	payoffMinPtr;
	double*	payoffMinBottomPtr= (double*) &(payoffMin[0]);
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


FDSolver1F::FDSolver1F() {
	// Clear memory.
	resetMem();
}

FDSolver1F::~FDSolver1F() {
	// Clear memory.
	deleteMem();
}


void FDSolver1F::resetMem() {
	m_iMax=0;
	J1 = NULL;
	JdriftTerm = NULL;
	JvolTermDown = NULL;
	JvolTermUp = NULL;
	m_rhs = NULL;
	spot1 = NULL;
	left1 = NULL;
	diag1 = NULL;
	right1 = NULL;
	gam = NULL;
    Vtemp = NULL;
    Voldtemp = NULL;
	m_nonNegativeValues = true;
	m_coordinateType = 0;
	m_driftSpread = false;
	m_payoffMax = 0;
	m_theta = 0.5;
	m_cs = NULL;
}

int FDSolver1F::allocateMem(int iMax) {

	int rc=SUCCESS;

	m_iMaxMem = iMax;

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
    Vtemp = (double *)calloc(iMax + 1, sizeof(double));
    Voldtemp = (double *)calloc(iMax + 1, sizeof(double));

	return rc;
}

void FDSolver1F::deleteMem() {
	free(m_rhs);
	free(J1);
	free(JdriftTerm);
	free(JvolTermDown);
	free(JvolTermUp);
	free(spot1);
	free(left1);
	free(diag1);
	free(right1);
	free(gam);
	free(Vtemp);
	free(Voldtemp);

	resetMem();
    return;

}

int FDSolver1F::getGridSpots(
	double		S1Min,
	double		S1Max,
	int			iMax,
	double*		spots,	// (O) spot vector to be populated
	double*		deltax,	// (O) physical step size of the grid
	int			coordinateType,	
	double		K) {
	double	temp;
	double	xMin,xMax;
	double	dx;
//	double	dxExp;
	int		i;

	if (S1Min >= S1Max)
		return FAILURE;

	switch (coordinateType)
	{
	
		case	COORD_TYPE_SYNH:	// sinh
		temp = S1Min/K;
		xMin = log(temp+sqrt(1+temp*temp));
		temp = S1Max/K;
		xMax = log(temp+sqrt(1+temp*temp));
		break;

		case	COORD_TYPE_LINEAR:	// linear
		xMin = S1Min;
		xMax = S1Max;
		break;

		case	COORD_TYPE_EXP:	// exponential
		case	3:
		xMax = log(S1Max/K);
		xMin = log(S1Min/K);
		break;

		case	4:
		xMin = 0;
		xMax = 1;
		break;

		default:
		return FAILURE;
	}

    dx = (xMax - xMin) / (double)(iMax);

    // Initialize coordinate arrays for S1 
	switch (coordinateType)
	{
	
		case	COORD_TYPE_SYNH:	// sinh
		for (i=0; i<=iMax; i++)
			spots[i] = K*sinh(dx*i + xMin);

		break;

		case	COORD_TYPE_LINEAR:	// linear
		for (i=0; i<=iMax; i++) 
			spots[i] = dx*i + xMin;
		break;

		case	COORD_TYPE_EXP:	// exponential
		case	3:	
		for (i=0; i<=iMax; i++)
			spots[i] = K*exp(dx*i + xMin);
		break;

		case	4:	
		for (i=0; i<=iMax; i++) {
			spots[i] = ::pow(S1Min,1. - dx*i)*::pow(S1Max,dx*i);
		}
		break;

		default:
		return FAILURE;
	}

	*deltax = dx;

	return SUCCESS;
}

int FDSolver1F::createGrid(
	const double*		spots,
	int					iMax,
	int					noCheck,
    double              dxIn) {
	int		i;

	if (iMax > m_iMaxMem)
		return FAILURE;

	if (!noCheck)
	{

		for (i=0; i<iMax; i++)
			if (spots[i] >= spots[i+1])
				return FAILURE;
	}

	if (spots != spot1)
		memcpy(spot1,spots,(iMax+1)*sizeof(double));

    if (dxIn == 0)
        dx = 1. / iMax;     // Set default value
	else
		dx = dxIn;

    for (i=0; i<iMax; i++)
	{
        J1[i] = (spot1[i + 1] - spot1[i]) / dx;
	}
	// Precompute for performance
	JdriftTerm[0] = 1. / (J1[0]*dx);
	JdriftTerm[iMax] = 1. / (J1[iMax-1]*dx);
    double dx2 = 0.5*dx;
    for (i=1; i<iMax; i++)
	{
		JdriftTerm[i] = 1. / ((J1[i-1] + J1[i])*dx);
		JvolTermDown[i] = JdriftTerm[i]/(J1[i-1]*dx2);	
		JvolTermUp[i] = JdriftTerm[i]/(J1[i]*dx2);
	}
	m_iMax = iMax;

	return SUCCESS;

}


int FDSolver1F::createGrid(
	double		S1Min,
	double		S1Max,
	int			iMax,
	int			coordinateType,
	double		K) {

	if (iMax > m_iMaxMem)
		return FAILURE;

    // Generate vector of spots for the grid
	if (FDSolver1F::getGridSpots(S1Min,S1Max,iMax,spot1,&dx,coordinateType,K) != SUCCESS)
		return FAILURE;

    // Create the grid
	if (createGrid(spot1,iMax,1,dx) != SUCCESS)
		return FAILURE;

	return SUCCESS;

}

int FDSolver1F::insertGridPoint(
	double		S) {
	int idx,i;

	if (S > spot1[m_iMax])
		return FAILURE;
	for (idx=m_iMax-1; idx>= 0 ; idx--)
	{
		if (S > spot1[idx])
			break;
	}
	if (idx<0)
		return FAILURE;

	if (S==spot1[idx+1])
		return SUCCESS;

	// (spot1[idx] < S < spot1[idx+1]
	// Insert a point
	if (m_iMax == m_iMaxMem)
		return FAILURE;

	for (i=m_iMax; i> idx; i--)
	{
		spot1[i+1] = spot1[i];
	}
	spot1[idx+1] = S;

	m_iMax++;

    for (i=0; i<m_iMax; i++)
	{
        J1[i] = (spot1[i + 1] - spot1[i]) / dx;
	}
	// Precompute for performance
	JdriftTerm[0] = 1. / (J1[0]*dx);
	JdriftTerm[m_iMax] = 1. / (J1[m_iMax-1]*dx);
    double dx2 = 0.5*dx;
    for (i=1; i<m_iMax; i++)
	{
		JdriftTerm[i] = 1. / ((J1[i-1] + J1[i])*dx);
		JvolTermDown[i] = JdriftTerm[i]/(J1[i-1]*dx2);	
		JvolTermUp[i] = JdriftTerm[i]/(J1[i]*dx2);
	}
	return SUCCESS;
}


int FDSolver1F::solve(
	int			        iMax,		        // Index of the last grid point
    FDTermStructureSP   driftTerm,          // (I) Vector of drift terms in the differential equation
    FDTermStructureSP   diffusionTerm,      // (I) Vector of diffusion terms in the differential equation
    FDTermStructureSP   discountTerm,       // (I) Vector of discounting terms in the differential equation 
    FDTermStructureSP   couponTerm,         // (I) Vector of coupon terms in the differential equation 
	double		        dt,
	double*		        Vnew,
	const double*	    V,
	double		        upBarrier1,
	double		        upBC1,
	double		        downBarrier1,
	double		        downBC1,
	const double*	    payoffMin,
	const double*	    Vold,		        // Old option prices for the three time level method
	int			        iMin
    ) const {
    int i;

	double	discountTerm1;
	double	driftTerm1,volTerm1;
	double	alphaUp=1.;
	double	alphaDown=1.;
	double	pup1,pdown1;
	int		hasUpBarrier1   = false;
	int		hasDownBarrier1 = false;
	int		threeTimeLevel  = false;
	int		m;
	double	theta;
    double  coupTerm=0;

	m = iMax;
    if (iMax > m_iMax) {
        throw ModelException("FDSolver1F::solve",
                        "Number of stock steps exceeds maximum number of stock steps");
    }

    if (Maths::isZero(dt)) {
        throw ModelException("FDSolver1F::solve",
                        "Zero time between two time points is not supported at the moment");
    }

	if (m_theta > 1)
	{
		threeTimeLevel=true;
		theta = 1;	// Do implicit
		if (Vold == NULL)
			return FAILURE;
	}
	else
		theta=0.5;

	if (upBarrier1 > 0)
	{
		double upB = upBarrier1*0.9999999999;
//		double upB = upBarrier1*0.999;
		while (spot1[iMax] >= upB)
			iMax--;

		if (iMax < m)
			iMax++;
		else
			upBarrier1 = -1; // spot1[m] < upB so the barrier will be ignored

		alphaUp = (spot1[iMax]-spot1[iMax-1])/(upBarrier1-spot1[iMax-1]);
		if (alphaUp > 10)
			iMax--;

		hasUpBarrier1 = true;
	}

	if (downBarrier1 > 0)
	{
		double downB = downBarrier1*1.0000000001;
		while (spot1[iMin] <= downB)
			iMin++;
		if (iMin > 0)
			iMin--;
		else
			return FAILURE;

		hasDownBarrier1 = true;
	}

    if (couponTerm->isFlat())
        coupTerm = (*couponTerm)(0)*dt;

	// Set up the system  

    if ( discountTerm->isFlat() )
        discountTerm1 = (*discountTerm)(0)*dt;
    if ( driftTerm->isFlat() )
        driftTerm1 = (*driftTerm)(0)*dt;

	if (hasDownBarrier1)
	{
		for (i=0; i <= iMin; i++)
		{
			alphaDown = (spot1[iMin+1]-spot1[i])/(spot1[iMin+1]-downBarrier1);
			left1[i] = 0;
			diag1[i] = 1;
			right1[i] = -(1-alphaDown);
			m_rhs[i] = downBC1*alphaDown;
		}
	}
	else
	{
        if ( !discountTerm->isFlat() )
            discountTerm1 = (*discountTerm)(iMin)*dt;
        if ( !driftTerm->isFlat() )
		    driftTerm1 = (*driftTerm)(iMin)*dt/(J1[iMin]*dx);
        else
            driftTerm1 = (*driftTerm)(0)*dt/(J1[iMin]*dx);
		pup1 = driftTerm1;
		left1[iMin] = 0.;
		diag1[iMin] = 1.+ (discountTerm1 + pup1)*theta;
		right1[iMin] = - pup1*theta;
		if (threeTimeLevel)
		{
			diag1[iMin] += 0.5;
			m_rhs[iMin] = (4*V[iMin]-Vold[iMin])/2;
		}
		else
		{
			m_rhs[iMin] = (1. - (discountTerm1 + pup1)*(1-theta))*V[iMin] 
				+ pup1*(1-theta)*V[iMin+1];
		}

        if (!couponTerm->isFlat() )
            m_rhs[iMin] += (*couponTerm)(iMin)*dt;
        else if (coupTerm != 0)
            m_rhs[iMin] += coupTerm;

	}

	for (i=iMin+1; i<iMax; i++) 
    {
        if ( !discountTerm->isFlat() )
	        discountTerm1 = (*discountTerm)(i)*dt;
        if ( !driftTerm->isFlat() )
		    driftTerm1= (*driftTerm)(i)*dt*JdriftTerm[i];
        volTerm1 = (*diffusionTerm)(i)*dt*JvolTermDown[i];
//		volTerm1 = tempa*JvolTermDown[i];
		pdown1 = (volTerm1 - driftTerm1);
        volTerm1 = (*diffusionTerm)(i)*dt*JvolTermUp[i];
//		volTerm1 = tempa*JvolTermUp[i];
		pup1 = (volTerm1 + driftTerm1);
		left1[i] = - pdown1*theta;
		diag1[i] = 1.+ (discountTerm1 + pdown1 + pup1)*theta;
		right1[i] = - pup1*theta;

		if (threeTimeLevel)
		{
			diag1[i] += 0.5;
			m_rhs[i] = (4*V[i]-Vold[i])/2;
		}
		else
		{
			m_rhs[i] = (1. - (discountTerm1 + pdown1 + pup1)*(1-theta))*V[i] 
				+ pdown1*(1-theta)*V[i-1]
				+ pup1*(1-theta)*V[i+1];
		}
	}

    if ( !couponTerm->isFlat() ) {
	    for (i=iMin+1; i<iMax; i++)
           m_rhs[i] += (*couponTerm)(i)*dt;
    } else if (coupTerm != 0) {
	    for (i=iMin+1; i<iMax; i++)
           m_rhs[i] += coupTerm;
    }

	if (hasUpBarrier1)
	{
		for (i=iMax; i <= m; i++)
		{
			alphaUp = (spot1[i]-spot1[iMax-1])/(upBarrier1-spot1[iMax-1]);
			left1[i] = -(1-alphaUp);
			diag1[i] = 1;
			right1[i] = 0;
			m_rhs[i] = upBC1*alphaUp;
		}
	}
	else
	{
		// iMax == m in this case
        if ( !discountTerm->isFlat() )
            discountTerm1 = (*discountTerm)(iMax)*dt;
		driftTerm1= (*driftTerm)(iMax)*dt/(J1[iMax-1]*dx);
		pdown1 = - driftTerm1;
		left1[iMax] = - pdown1*theta;
		diag1[iMax] = 1.+ (discountTerm1 + pdown1)*theta;
		right1[iMax] = 0;
		if (threeTimeLevel)
		{
			diag1[iMax] += 0.5;
			m_rhs[iMax] = (4*V[iMax]-Vold[iMax])/2;
		}
		else
		{
			m_rhs[iMax] = (1. - (discountTerm1 + pdown1)*(1-theta))*V[iMax]
				+ pdown1*(1-theta)*V[iMax-1];
		}
        if (!couponTerm->isFlat() )
            m_rhs[iMax] += (*couponTerm)(iMax)*dt;
        else if (coupTerm != 0)
            m_rhs[iMax] += coupTerm;
	}


	// Solve using tridiagonal solvers

	triDiag1D(&(left1[iMin]),&(diag1[iMin]),&(right1[iMin]),&(m_rhs[iMin]),&(Vnew[iMin]),m-iMin);
	if (m_nonNegativeValues)
	{
		// This can improve convergence - dampen oscillations around a barrier
		for (i=iMin; i<=m; i++)
            Vnew[i] = Maths::max(Vnew[i],0.);
	}
	if (payoffMin)
	{
//		for (int i=0; i<=iMax; i++)
//			Vnew[i] = MAX(Vnew[i],payoff[i]);
		if (m_payoffMax>0)
			sorMinMax(&(left1[iMin]),&(diag1[iMin]),&(right1[iMin]),&(m_rhs[iMin]),&(payoffMin[iMin]),m_payoffMax,&(Vnew[iMin]),iMax-iMin,m_toleranceSOR);
		else
			sor(&(left1[iMin]),&(diag1[iMin]),&(right1[iMin]),&(m_rhs[iMin]),&(payoffMin[iMin]),&(Vnew[iMin]),iMax-iMin,m_toleranceSOR);

	}

	return SUCCESS;
}

int
FDSolver1F::solveQuick(
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
    int			n)	const {		// The number of iterations
    int i;

	double	discountTerm;
	double	driftTerm1,volTerm1;
	double	couponTerm;
	double	muTerm;
//	double	ds1;
	double	alphaUp=1.;
	double	alphaDown=1.;
//	double	theta=0.5;
	double	pup1,pdown1,pmid1;
	int		hasUpBarrier1   = false;
	int		hasDownBarrier1 = false;
//	int		threeTimeLevel  = false;
//	int		iMax;
	double	leftTerm=0,rightTerm=0,diagTerm=0;
	int		m;
	double	theta = m_theta;
//	int		rc=0;

	m = iMax;
	if (m > m_iMax)
		return FAILURE;

	if (upBarrier1 > 0)
	{
		double upB = upBarrier1*0.9999999999;
		while (spot1[iMax] >= upB)
			iMax--;
		if (iMax < m)
			iMax++;
		else
			return FAILURE;
		hasUpBarrier1 = true;
	}

	if (downBarrier1 > 0) {
        // got a few problems with DBL_EPSILON ( approx. 2e-16) here, thus using a bigger tolerance
        if ( fabs(downBarrier1-spot1[0]) < EPSILON) {
			alphaDown       = 1;
		    hasDownBarrier1 = true;
        } else if (downBarrier1 < spot1[0]) {
			alphaDown       = 0.0;
		    hasDownBarrier1 = false;
        } else if ( downBarrier1 >= spot1[1]) {
            throw ModelException("FDSolver1F::solveQuick", 
                                 "The down barrier must not be greater than the second spot price on the grid");
        } else {
			// The boundary does not fall on a grid point (it is between nodes 0 and 1)
			alphaDown       = (spot1[0]-spot1[1])/(downBarrier1-spot1[1]);
		    hasDownBarrier1 = true;
        }
	}


	// Set up the system  

	discountTerm = dt*r;
	couponTerm = dt*coupon;	// continuous coupon (e.g., accrued coupon)

	if (m_coordinateType != 3)
	{
		muTerm = (r - divYld1)*dt;

		if (hasDownBarrier1)
		{
			left1[0] = 0;
			diag1[0] = 1;
			right1[0] = -(1-alphaDown);
		}
		else
		{
			driftTerm1= muTerm*spot1[0]*JdriftTerm[0];
			pup1 = driftTerm1;
			left1[0] = 0.;
			diag1[0] = 1.+ (discountTerm + pup1)*theta;
			right1[0] = - pup1*theta;
		}

        double  dt2 = 0.5*dt;
		for (i=1; i<iMax; i++) {
			double	temp =  sigma1*spot1[i];
			double	tempa = dt2*temp*temp;
			driftTerm1 = muTerm*spot1[i]*JdriftTerm[i];
			volTerm1 = tempa*JvolTermDown[i];
			pdown1 = (volTerm1 - driftTerm1);
			volTerm1 = tempa*JvolTermUp[i];
			pup1 = (volTerm1 + driftTerm1);
			left1[i] = - pdown1*theta;
			diag1[i] = 1.+ (discountTerm + pdown1 + pup1)*theta;
			right1[i] = - pup1*theta;

		}

		if (hasUpBarrier1)
		{
			// V[iMax] = upBC1*alphaUp + (1-alphaUp)*V[iMax-1]
			for (i=iMax; i <= m; i++)
			{
				alphaUp = (spot1[i]-spot1[iMax-1])/(upBarrier1-spot1[iMax-1]);
				left1[i] = -(1-alphaUp);
				diag1[i] = 1;
				right1[i] = 0;
	//			m_rhs[i] = upBC1*alphaUp;
			}
		}
		else
		{
			// iMax == m in this case
			driftTerm1= muTerm*spot1[iMax]/(J1[iMax-1]*dx);
			pdown1 = - driftTerm1;
			left1[iMax] = - pdown1*theta;
			diag1[iMax] = 1.+ (discountTerm + pdown1)*theta;
			right1[iMax] = 0;
	//		m_rhs[iMax] = (1. - (discountTerm + pdown1)*(1-theta))*V[iMax]
	//			+ pdown1*(1-theta)*V[iMax-1] + couponTerm;
		}
	}
	else
	{
		driftTerm1 = 0.5*dt*(r - divYld1-sigma1*sigma1/2)/dx;
		volTerm1 = dt* ::pow(sigma1/dx,2)/2.;
		pdown1 = (volTerm1 - driftTerm1)*(1-theta);
		pup1 = (volTerm1 + driftTerm1)*(1-theta);
		pmid1 = 1.- (discountTerm*(1-theta) + pdown1 + pup1);
		leftTerm = -(volTerm1 - driftTerm1)*theta;
		rightTerm = - (volTerm1 + driftTerm1)*theta;
		diagTerm = 1 +  (discountTerm*theta - leftTerm - rightTerm);
	}
	// Solve using tridiagonal solvers
	for (int k=0; k< n; k++)
	{
		if (m_coordinateType != 3)
		{
			if (hasDownBarrier1)
				m_rhs[0] = downBC1*alphaDown;
			else
				m_rhs[0] = (2-diag1[0])*V[0] - right1[0]*V[1] + couponTerm;

			for (i=1; i<iMax; i++)
				m_rhs[i] = -left1[i]*V[i-1] + (2-diag1[i])*V[i] - right1[i]*V[i+1];

            if ( !Maths::equals(couponTerm,0))
			{
				for (i=1; i<iMax; i++)
					m_rhs[i] += couponTerm;
			}

			if (hasUpBarrier1)
				m_rhs[iMax] = upBC1*alphaUp;
			else
				m_rhs[iMax] = -left1[iMax]*V[iMax-1] + (2-diag1[iMax])*V[iMax] + couponTerm;

			triDiag1D(left1,diag1,right1,m_rhs,V,m);
		}
		else
		{
			for (i=1; i<iMax; i++)
				m_rhs[i] = leftTerm*V[i-1] + diagTerm*V[i] + rightTerm*V[i+1];

            if (!Maths::equals(coupon,0)) {
				for (i=1; i<iMax; i++)
					m_rhs[i] += couponTerm;
			}
			double frac = (spot1[1]-spot1[0])/(spot1[2]-spot1[1]);
			double frac1 = 1/frac;

			if (hasDownBarrier1)
			{
				// Adjust the co_efficients diag1[1] and rhs[1] due to 
				// elimination of V[0] since 
				// V[0] = downBC1*alphaDown + (1-alphaDown)*V[1]
				diag1[1] += left1[1]*(1-alphaDown);
  				m_rhs[1] -= left1[1]*downBC1*alphaDown;
			}
			else
			{
				// Adjust the co_efficients diag1[1] and right1[1] due to
				// elimination of V[0] using linearity condition equation
				// V[0]=V[1]-(V[2]-V[1])*frac
				diag1[1] += left1[1]*(1+frac);
				right1[1] -= left1[1]*frac;
			}
			// New left[1] is 0 since the V[0] is eliminated 
			left1[1] = 0.0;

			if (hasUpBarrier1)
			{
 				// Adjust the co_efficients diag1[iMax-1] and m_rhs[iMax-1] due to 
				// elimination of V[iMax] since 
				// V[iMax] = upBC1*alphaUp + (1-alphaUp)*V[iMax-1]) - old method
				diag1[iMax-1] += right1[iMax-1]*(1-alphaUp);
				m_rhs[iMax-1] -= right1[iMax-1]*upBC1*alphaUp;
			}
			else
			{
				// Adjust the co_efficients left1[iMax-1] and diag1[iMax-1] due to
				// elimination of V[iMax] using linearity condition equation
				// V[m] = V[m-1] + (V[m-1]- V[m-2])*frac 
				left1[iMax-1] -= right1[iMax-1]*frac1;
				diag1[iMax-1] += right1[iMax-1]*(1+frac1);
			}
			// New right[iMax-1] is 0 since the V[iMax] is eliminated 
			right1[iMax-1]= 0.0;
					
			triDiag1D(&(left1[1]),&(diag1[1]),&(right1[1]),&(m_rhs[1]),&(V[1]),m-2);

			if (hasDownBarrier1)
				V[0] = downBC1*alphaDown+V[1]*(1-alphaDown);
			else
				V[0] = V[1] - (V[2]- V[1])*frac;

			if (hasUpBarrier1)
				V[iMax] = upBC1*alphaUp+V[iMax-1]*(1-alphaUp); // Old method
			else
				V[iMax] = V[iMax-1] + (V[iMax-1]- V[iMax-2])*frac1;

		}

	}
    return SUCCESS;
}

int FDSolver1F::solveBarrier(
	int			      m,
    FDTermStructureSP driftTerm,                // (I) Vector of drift terms in the differential equation
    FDTermStructureSP diffusionTerm,            // (I) Vector of diffusion terms in the differential equation
    FDTermStructureSP discountTerm,             // (I) Vector of discounting terms in the differential equation 
    FDTermStructureSP couponTerm,               // (I) Vector of coupon terms in the differential equation 
	double		      dt,
	double*		      Vnew,
	const double*	  V,
	double		      upBarrierNew,
	double		      upPayoutNew,
	double		      upPayoutDeltaNew,
	double		      upBarrier,
	double		      upPayout,
	double		      upPayoutDelta,
	double		      downBarrierNew,
	double		      downPayoutNew,
	double		      downPayoutDeltaNew,
	double		      downBarrier,
	double		      downPayout,
	double		      downPayoutDelta,
	const double*	  payoffMin,
	const double*	  Vold,		// Old option prices for the three time level method
	double		      upBarrierOld,
	double		      upPayoutOld,
	double		      upPayoutDeltaOld,
	double		      downBarrierOld,
	double		      downPayoutOld,
	double		      downPayoutDeltaOld
   ) const {
	int		i;
	int		rc=0;
//	double	Vsaved[50];
//	double	VsavedOld[50];
//	double	Vsaved1[50];
//	double	VsavedOld1[50];
//	double	Vtemp=0,Vtemp1=0;
//	double	Voldtemp=0,Voldtemp1=0;
	int		iMax;
	int		iMaxOld;
	int		iMin;
	int		iMinOld;
//	int		hasUpBarrier = false;
//	double	upB = 0;
	double	alphaUp=1;
	double	alphaDown=1;
	int		numSaved=0;
	double	upBCNew=0;
	double	downBCNew=0;
	double	tol=TOL;
	double*	Vptr=NULL;
	double*	Voldptr=NULL;

	iMax = m;
	iMin = 0;
	if (upBarrierNew > 0)
	{
		// Find iMax for Vnew
		double upB = upBarrierNew*0.9999999999;

        if (spot1[0] >= upB) 
		{ // all points are above the barrier
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

//          return FAILURE;
//			upBarrierNew = -1; // spot1[m] < upB so the barrier will be ignored
	}
	if (downBarrierNew > 0)
	{
		// Find iMin for Vnew
		double downB = downBarrierNew*1.0000000001;
//		upB = upBarrierNew*0.999;
        if (spot1[m] <= downB) 
		{ // all points are below the barrier
            for (i=m; i>= 0; i--) {
			    Vnew[i] = downPayoutNew+downPayoutDeltaNew*spot1[i];;
            }
            return SUCCESS;
        }
		while (iMin < iMax && spot1[iMin] <= downB)
			iMin++;

		if (iMin > 0)
			iMin--;
		else
			downBarrierNew = -1; // spot1[m] < upB so the barrier will be ignored

//          return FAILURE;
//			upBarrierNew = -1; // spot1[m] < upB so the barrier will be ignored
	}

	if (upBarrier > 0 || downBarrier > 0)
	{
		memcpy(Vtemp,V,(m+1)*sizeof(double));
		Vptr = Vtemp;
	}
	else
		Vptr = (double*) V;

	iMaxOld = iMax;
    /*
	if (upBarrier > 0)
	{
		// Find iMax for V
		double upB = upBarrier*0.999;
//		upB = upBarrier*0.9999999999;
		while (spot1[iMaxOld] >= upB)
			iMaxOld--;

		if (iMaxOld < iMax)
		{
			iMaxOld++;
			alphaUp = (spot1[iMaxOld]-spot1[iMaxOld-1])/(upBarrier-spot1[iMaxOld-1]);
			if (alphaUp > 10)
				iMaxOld--;
			double upBC = upPayout + upPayoutDelta*upBarrier;
			if (fabs(upBarrier/spot1[iMaxOld-1]-1) < 1e-3)
				return FAILURE;
			for (i=iMaxOld; i <= iMax; i++)
			{
				alphaUp = (spot1[i]-spot1[iMaxOld-1])/(upBarrier-spot1[iMaxOld-1]);
				Vtemp[i] = upBC*alphaUp + (1-alphaUp)*Vtemp[iMaxOld-1];
			}
		}
	}
    */

    double Voldtmp  = 0.0;
    double Voldtmp1 = 0.0;
    if (upBarrier > 0) {
		// If the old  barrier is greater or equal than iMax, ignore it, otherwise
		// adjust option prices above the barrier up to spot1[iMax]
		if (spot1[iMax] >= upBarrier*0.9999999999)
		{
			double upBCOld = upPayout + upPayoutDelta*upBarrier;
			Voldtmp = Vptr[iMax];	// Save the old value
			if (upBarrier >= spot1[iMax-1]) {
				alphaUp = (spot1[iMax]-spot1[iMax-1])/(upBarrier-spot1[iMax-1]);
				Vptr[iMax] = upBCOld*alphaUp + (1-alphaUp)*Vptr[iMax-1];
				numSaved = 1;
            } else {
                // The barrier is below spot1[iMax-1]
                if ( iMax <= 1) {
                    throw ModelException("ZFDSolver1D::solveBarrier", 
                        "The barrier is very close to the lowest spot price, which can cause unstable numerical solutions. Please update the model paramters");
                }

				Voldtmp1 = Vptr[iMax-1];	// Save the old value
				alphaUp = (spot1[iMax]-spot1[iMax-2])/(upBarrier-spot1[iMax-2]);
				Vptr[iMax] = upBCOld*alphaUp + (1-alphaUp)*Vptr[iMax-2];
				alphaUp = (spot1[iMax-1]-spot1[iMax-2])/(upBarrier-spot1[iMax-2]);
				Vptr[iMax-1] = upBCOld*alphaUp + (1-alphaUp)*Vptr[iMax-2];
				numSaved = 2;
			}
		}
	}

	iMinOld = iMin;
	if (downBarrier > 0)
	{
		// Find iMin for V

		while (iMinOld < iMaxOld && spot1[iMinOld] <= downBarrier + (spot1[iMinOld+1]-spot1[iMinOld])*tol)
			iMinOld++;

		if (iMinOld > iMin)
		{
			iMinOld--;
			double downBC = downPayout + downPayoutDelta*downBarrier;
			if (fabs(downBarrier/spot1[iMinOld+1]) < tol)
				return FAILURE;
			for (i=iMinOld; i >= iMin; i--)
			{
				alphaDown = (spot1[iMinOld+1]-spot1[i])/(spot1[iMinOld+1]-downBarrier);
				Vtemp[i] = downBC*alphaDown + (1-alphaDown)*Vtemp[iMinOld+1];
			}
		}
	}

	if (m_theta > 1 && (upBarrier > 0 || downBarrier > 0))
	{
		memcpy(Voldtemp,Vold,(m+1)*sizeof(double));
		Voldptr = Voldtemp;
	}
	else
		Voldptr = (double*) Vold;

	iMaxOld = iMax;
	if (m_theta > 1 && upBarrierOld > 0)
	{
		if (Vold == NULL)
			return FAILURE;
		// Adjust Vold for the barrier
		// If the old  barrier is greater or equal than iMax, ignore it, otherwise
		// adjust option prices above the barrier up to spot1[iMax]
		double upB = upBarrierOld*0.999;
//		upB = upBarrierOld*0.9999999999;
		while (spot1[iMaxOld] >= upB)
			iMaxOld--;

		if (iMaxOld < iMax)
		{
			iMaxOld++;
			double upBC = upPayoutOld + upPayoutDeltaOld*upBarrierOld;
			if (fabs(upBarrierOld/spot1[iMaxOld-1]-1) < 1e-3)
				return FAILURE;
			for (i=iMaxOld; i <= iMax; i++)
			{
				alphaUp = (spot1[i]-spot1[iMaxOld-1])/(upBarrierOld-spot1[iMaxOld-1]);
				Voldtemp[i] = upBC*alphaUp + (1-alphaUp)*Voldtemp[iMaxOld-1];
			}
		}
	}

	iMinOld = iMin;
	if (m_theta > 1 && downBarrierOld > 0)
	{
		if (Vold == NULL)
			return FAILURE;
		// Adjust Vold for the barrier
		// If the old  barrier is greater or equal than iMax, ignore it, otherwise
		// adjust option prices above the barrier up to spot1[iMax]
		while (iMinOld < iMaxOld && spot1[iMinOld] <= downBarrierOld + (spot1[iMinOld+1]-spot1[iMinOld])*tol)
			iMinOld++;

		if (iMinOld > iMin)
		{
			iMinOld--;
			double downBC = downPayoutOld + downPayoutDeltaOld*downBarrierOld;
			double	alphaDownInv = (spot1[iMinOld+1]-downBarrierOld)/(spot1[iMinOld+1]-spot1[iMinOld]);
			if (fabs(alphaDownInv) < tol)
				return FAILURE;
			for (i=iMinOld; i >= iMin; i--)
			{
				alphaDown = (spot1[iMinOld+1]-spot1[i])/(spot1[iMinOld+1]-downBarrierOld);
				Voldtemp[i] = downBC*alphaDown + (1-alphaDown)*Voldtemp[iMinOld+1];
			}
		}
	}
	if (upBarrierNew > 0)
		upBCNew = upPayoutNew + upPayoutDeltaNew*upBarrierNew;
	if (downBarrierNew > 0)
		downBCNew = downPayoutNew + downPayoutDeltaNew*downBarrierNew;
	rc = solve (iMax,driftTerm,diffusionTerm,discountTerm,couponTerm,
            dt,Vnew,Vptr,upBarrierNew,upBCNew,downBarrierNew,downBCNew,
						payoffMin,Voldptr,iMin);
	if (rc != SUCCESS)
		return rc;

	// Restore old option values above the boundary
    if (numSaved > 0) {
		Vptr[iMax] = Voldtmp;
		if (numSaved > 1)
			Vptr[iMax-1] = Voldtmp1;
	}
    


	if (upBarrierNew > 0)
	{
		// Set new options values which are above the boundary
		for (i=m; i>= iMax; i--)
			Vnew[i] = upPayoutNew+upPayoutDeltaNew*spot1[i];
	}

	if (downBarrierNew > 0)
	{
		// Set new options values which are above the boundary
		for (i=0; i<= iMin; i++)
			Vnew[i] = downPayoutNew+downPayoutDeltaNew*spot1[i];
	}

	if (rc == FAILURE)
		return rc;


	return rc;
}


DRLIB_END_NAMESPACE
