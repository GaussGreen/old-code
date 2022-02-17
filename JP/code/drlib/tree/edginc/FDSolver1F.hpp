//---------------------------------------------------------------------------
//
//   Group       : Convertibles Research
//
//   Filename    : FDSolver.hpp
//
//   Description : A 1 factor FD solver base class
//
//   Author      : André Segger
//
//   Date        : 29 September 2003
//
//---------------------------------------------------------------------------

#ifndef FD_SOLVER_1FB_HPP
#define FD_SOLVER_1FB_HPP

#include <string>
#include "edginc/config.hpp"
#include "edginc/Object.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/FDTermStructure.hpp"

using namespace std;    // string

DRLIB_BEGIN_NAMESPACE

class TREE_DLL FDSolver1F
{
public:

	enum {COORD_TYPE_SYNH = 0, 
		COORD_TYPE_LINEAR = 1, 
		COORD_TYPE_EXP = 2
	};

	// Constructor
	FDSolver1F();

	// Allocate memory for the Finite Difference grid. The maximum amount of memory that will be used by the object
	// should be allocated. The input parameter m is the maximum index for a grid point in the spot dimension 
	// (e.g., S[0], S[1],...,S[m]. Therefore the grid size is m + 1). 
	int
	allocateMem(int iMax);

	// Free allocated memory
	void
	deleteMem();

	// Get the type of the coordinates used in the grid 
	int	getCoordinateType() const {return m_coordinateType; }

	// Methods for manipulating tolerance parameter for the Simulatanous Overrelaxation Method (SOR)
	double	getToleranceSOR() const {return m_toleranceSOR; }
	void setToleranceSOR (double toleranceSOR) {m_toleranceSOR = toleranceSOR; }

	// Discretization methods: 0 - Explicit, 0.5, - Crank Nicholson, 1 - Implicit, 2 - Three Time Level
	double getDiscretizationMethod() const {return m_theta; }
	void setDiscretizationMethod(double theta) {m_theta = theta; }

	// Enforces that values return by the solver are non-negative (
	void forceNonNegative(int	nonNegative=true) {m_nonNegativeValues = nonNegative ? true : false;}
	int	isNonNegative() const {return m_nonNegativeValues; } 


	// Creates an FD grid: creates a grid from a vector of spots prices created externally
    int
    createGrid(
	    const double*		spots,          // (I) vector of spot prices
	    int					m,              // (I) maximum index for a grid point 
	    int					noCheck=0,      // (I) don't check for validity of the grid if set
        double              dx=0);          // (I) physical step size of the grid (default = 1./iMax)

    // Creates an FD grid: creates spots vector internally
	int
	createGrid(
		double			spotMin,			// (I) mininum spot in the grid 
		double			spotMax,			// (I) maximum spot in the grid
		int				m,					// (I) maximum index for a grid point 
		int				coordinateType=0,	// (I) coordinate type to be used for the grid
		double			gridStrike=1.);

	// Get a pointer to the internally created vector of spot points in the grid (e.g., S[0], S[1],...,S[m])
	// Note: createGrid method needs to be called first to create a spot vector
	inline double*
		getSpotsPtr() const { return spot1;}

	// Finite difference solver method. 0 <= iMin < iMax <= m
	// Vnew represents option values for the current (k-th) time point and V values from the previous ((k+1)-th) 
	// time point. Vold are values from the (k+2)-th time point and need to be provided only when using
	// the three time level discretization method
	int
	solve(
		int				  iMax,			  // (I) Index of the last grid point
        FDTermStructureSP driftTerm,      // (I) Vector of drift terms in the differential equation
        FDTermStructureSP diffusionTerm,  // (I) Vector of diffusion terms in the differential equation
        FDTermStructureSP discountTerm,   // (I) Vector of discounting terms in the differential equation 
        FDTermStructureSP couponTerm,     // (I) Vector of coupon terms in the differential equation 
		double			  dt,			  // (I) Time step
		double*			  Vnew,			  // (O) New values
		const double*	  V,			  // (I) Previous values
		double			  upBarrier=0,	  // (I) Up barrier (<=0 - no barrier)
		double			  upPayout=0,	  // (I) Payout at the up barrier 
		double			  downBarrier=0,  // (I) Down barrier (<=0 - no barrier)
		double			  downPayout=0,	  // (I) Payout at the down barrier
		const double*	  payoffMin=NULL, // (I) Intrinsic value for the American options calculated using SOR method 
		const double*	  Vold=NULL,	  // (I) Old values for the three time level discretization method
		int				  iMin=0		  // (I) Index of the first grid point
    ) const;

	// An enhanced version of the solve method that can handle arbitrary barriers (e.g., constant or moving barriers).
	// The barriers must be within the range of grid points (S[0] <= downBarrier < upBarrier <= S[m]). The barrier 
	//payout is a linear function of spot (e.g., Vnew[i] = upPayoutDeltaNew*S[i] + upPayoutNew if S[i] >= upBarrierNew)
	int
	solveBarrier(
		int				  m,					 // (I) Index of the last grid point
        FDTermStructureSP driftTerm,             // (I) Vector of drift terms in the differential equation
        FDTermStructureSP diffusionTerm,         // (I) Vector of diffusion terms in the differential equation
        FDTermStructureSP discountTerm,          // (I) Vector of discounting terms in the differential equation 
        FDTermStructureSP couponTerm,            // (I) Vector of coupon terms in the differential equation 
		double			  dt,					 // (I) Time step
		double*			  Vnew,				     // (O) New values
		const double*	  V,					 // (I) Previous values
		double			  upBarrierNew=0,		 // (I) Up barrier at the current time point (<=0 - no barrier)
		double			  upPayoutNew=0,		 // (I) Payout coefficient at and above the up barrier at the current time point
		double			  upPayoutDeltaNew=0,	 // (I) Payout delta coeffient at and above the up barrier at the current time point
		double			  upBarrier=0,		     // (I) Up barrier at the previous time point (<=0 - no barrier)
		double			  upPayout=0,			 // (I) Payout coefficient at and above the up barrier at the previous time point
		double			  upPayoutDelta=0,	     // (I) Payout delta coeffient at and above the up barrier at the previous time point
		double			  downBarrierNew=0,	     // (I) Down barrier at the current time point (<=0 - no barrier)
		double			  downPayoutNew=0,	     // (I) Payout coefficient at and below the down barrier at the current time point
		double			  downPayoutDeltaNew=0,  // (I) Payout delta coeffient at and below the down barrier at the current time point
		double			  downBarrier=0,		 // (I) Down barrier at the previous time point (<=0 - no barrier)
		double			  downPayout=0,		     // (I) Payout coeffient at and below the down barrier at the previous time point
		double			  downPayoutDelta=0,	 // (I) Payout delta coeffient at and below the down barrier at the previous time point
		const double*	  payoffMin=NULL,		 // (I) Intrinsic value for the American options calculated using SOR method 
		const double*	  Vold=NULL,			 // (I) Old values for the three time level discretization method
		double			  upBarrierOld=0,		 // (I) Up barrier at the old time point (<=0 - no barrier)
		double			  upPayoutOld=0,		 // (I) Payout coefficient at and above the up barrier at the old time point
		double			  upPayoutDeltaOld=0,	 // (I) Payout delta coeffient at and above the up barrier at the old time point 
		double			  downBarrierOld=0,	     // (I) Down barrier at the old time point (<=0 - no barrier)
		double			  downPayoutOld=0,	     // (I) Payout coeffient at and below the down barrier at the old time point
		double			  downPayoutDeltaOld=0   // (I) Payout delta coeffient at and below the down barrier at the old time point
	) const;

	// A trimmed down and faster version of the solve method, which can be used when all the parameters are constant.
	// (e.g, intRate, divYld, vol, coupon, barriers and payouts).
	// It does not support local vols and three time level discretization method. 
	int
	solveQuick(
		int				m,				// (I) Index of the last grid point
		double			intRate,		// (I) Continuous interest rate
		double			divYld,			// (I) Continuous dividend yield
		double			vol,			// (I) Vector of volatilities
		double			coupon,			// (I) Continuous coupon
		double			dt,				// (I) Time step
		double*			V,				// (I/O) Initial (previous) and the final value
		double			upBarrier=0,	// (I) Up barrier (<=0 - no barrier)
		double			upPayout=0,		// (I) Payout at the up barrier 
		double			downBarrier=0,	// (I) Down barrier (<=0 - no barrier)
		double			downPayout=0,	// (I) Payout at the down barrier
		int				n=1				// (I) The number of time steps
	) const;	

	// virtual Destructor
	~FDSolver1F();

	// Populates a given spot vector with grid points. 
	static int
	getGridSpots(
		double			spotMin,			// (I) mininum spot in the grid 
		double			spotMax,			// (I) maximum spot in the grid
		int				m,					// (I) maximum index for a grid point 
		double*			spots,				// (O) spot vector to be populated
		double*			dx,					// (O) physical step size of the grid
		int				coordinateType=0,	// (I) coordinate type to be used for the grid
		double			gridStrike=1.		// (I) grid strike
	);

	// Tridiagonal solver
	int triDiag1D(
	    const double*		a, 
	    const double*		b, 
	    const double*		c, 
	    const double*		r, 
	    double*				u,
	    int					n) const;

protected:
	double		dx;
	int			m_iMaxMem;	// maximum memory allocated
	int			m_iMax;		
	double*		J1;
	double*		m_rhs;
	double*		spot1;
	double*		left1;
	double*		diag1;
	double*		right1;
	double*		gam;
	int			m_nonNegativeValues;
	int			m_coordinateType;
	int			m_driftSpread;
	double		m_toleranceSOR;
	double		m_payoffMax;
	double		m_theta;
	// Precompute for performance
	double*		JdriftTerm;
	double*		JvolTermDown;
	double*		JvolTermUp;
	double*		m_cs;
	double*		Vtemp;
	double*		Voldtemp;

	static int
	sor(
		const double*		a,
		const double*		b,
		const double*		c,
		const double*		d,
		const double*		payoff,
		double*				u,
		int					m,
		double				tolerance=1.0e-9
	);

	static int
	sorMinMax(
		const double*		a,
		const double*		b,
		const double*		c,
		const double*		d,
		const double*		payoffMin,
		double		payoffMax,
		double*		u,
		int			n,
		double		tolerance=1.0e-9
	);

	void
	resetMem();

	// This needs more testing
	int
	insertGridPoint(
		double		S);

private:

	// The following methods are currently not implemented and should not be used
	FDSolver1F(								// COPY constructor
		const FDSolver1F&);					// (I) Original object

	FDSolver1F& operator=(const FDSolver1F&); // Assignment operator

};


DRLIB_END_NAMESPACE
#endif
