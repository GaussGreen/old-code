  
 
#ifndef GRGPROBLEM_H
#define GRGPROBLEM_H
#include <vector>
#include <assert.h>
#include "CMatrix.h"

using namespace std;

//
// A.G. Forward declare the globals
//
struct LsgrgGlobals_t;

class LSGRGSolver
{	
public:
	//friend void p_gcomp( double* gVals,double *xval,void* vp);
	//friend void partialDerivs(double* x, int n, double* paij, int* iprow, int* ipcol, int *nnz,void* vp);

	void p_gcomp( double* gVals,double *xval);
	//void partialDerivs(double* x, int n, double* paij, int* iprow, int* ipcol, int *nnz);
	void partialDerivs(double* x, long n, double* paij, long* iprow, long* ipcol, long* nnz);	

protected:
	int		m_numberVariables;
	int		m_numberConstraints;
	iVector	m_linearVariableIndex;//[ilinearvar]
	double		m_InitialTolerance;
	double		m_finalTolerance;
	double		m_StoppingTolerance;
	
	//LsgrgInfo *m_info;
	//
	// A.G.
	//
	LsgrgGlobals_t *m_info;
	
	// columns [.][0]:lower bound. [.][1] higher bound
	CMatrix	m_BoundsOnConstraints;// [iconstraint][2]
	CMatrix	m_xValBounds;// [ivariable][2]
	CVector	m_xInitialGuess;
	iVector	m_nonZeroJacobianMapRows;//[nonzeroJacobians]->[iVariables]
	iVector	m_nonZeroJacobianMapCols;//[nonzeroJacobians]->[iConstraint]
	int		m_numberOfNonZeroJacobians;
	CVector	m_partialDerivatives;//[nonzeroJacobians]
	iMatrix	m_inverseJacobianMap;//[m_numberVariables][m_numberConstraints]->nonzeroJacobians
	CMatrix	m_NonZeroPartialDerivatives;//[m_numberVariables][m_numberConstraints]
	// 1 if partial is nonzero 0 else
	int		m_partialDerivsAreSet;
			
public:
	 class LSGRGSolverSolution
	 {
		  public:
		  CVector m_finalxValues;
		  CVector m_finalFunctionValues;
		  CVector m_finalLagrangeMultiplyers;
		  CVector m_finalReducedGradients;
	 }m_LSGRGSolverSolution;
	
	void	setPartialDerivative(double setVal,int ixVal,int iConstraint);
	double	getPartialDerivative(int ixVar,int iConstraint);
	
	void	initialize(CVector& initialGuess,iVector& linearVariableIndex,
					   CMatrix& xValBounds,CMatrix& ObjectiveBounds,
					   double InitialTolerance,double FinalTolerance,
					   double StoppingTolerance,
					   CMatrix& NonZeroPartialDerivativesSpecification,int outputFlag);
	
	
	virtual void  ConstraintFcn(double* contraintVals, double* xVals);
	virtual void  PartialDerivatives(CVector& xVal,CVector& partialDerivatives){};
	virtual void  ObjectiveFcn(double* gVals,double* xVals);
	
	LSGRGSolver();
	
	void solve( CVector& x, int maximizeFlag, int& returnCode);
};	

#ifdef USER_TERMINATION_ENABLED
	typedef long (LSGRGSolver::*P_GCOMP)  (double *, double *);
	typedef long (LSGRGSolver::*P_GCOMPX) (double *, double *, long*, long[]);
#else
	typedef void (LSGRGSolver::*P_GCOMP)  (double *, double *);
	typedef void (LSGRGSolver::*P_GCOMPX) (double *, double *, long*, long[]);
#endif
	typedef void (LSGRGSolver::*P_PARSH)  (double[],long, double[],long[],long[], long*);

//
// A.G. This requires the typedef's and the definition of LSGRGSolver
//
#include "lsinfo.h"

class RungeKuttaIntegrator
{
	
	friend void integratorFcn(double, double[], double[],void *vp);

	CVector m_res;

	double m_h_min;
	double m_eps;
	double m_bump;
	double m_h_init;
	double m_scale;

	public:

	void initialize(int nfunctions,CVector& integrationLimits,double eps=1e-3,double hInit=0.1,double hMin=1e-2,double bump=0.0008);
	virtual void integrands(CVector& res, double x){};// integrates vector of functions
	void integrate(CVector& integral);

	CVector m_IntegrationLimits;

};


enum enum_underdetermined_gradient{eUNDERDETERMINED_CALC_BEST=0,eUNDERDETERMINED_INDEX_MASK=16383,
eUNDERDETERMINED_SETUP_GRADIENT=16384,eUNDERDETERMINED_CALC_GRADIENT=32768};
	
class cSparseLinSolve
{	
	void linbcg(CVector& i_b);
	double snrm(int n,CVector & sx,int itol);
	CVector sx;
	void atimes(int n,CVector & x, CVector & y, int itrnsp);
	void asolve(int n,CVector& b,CVector& x,int itrnsp);
	void dsprstx(CVector& sa,iVector& ija,CVector & x,CVector& b,int n);
	void dsprsax(CVector& sa, iVector& t_ija,CVector& t_x,CVector& t_b,int n);
	void sprsin();
	iVector ija;
	CVector sa;
	int nmax;
	int ndim;
	CVector b;
	CVector x;
	const CMatrix* a;

			
	public:
			
// solves a x = b if a is a sparse matrix
				
	cSparseLinSolve(const CMatrix& i_a, CVector &i_b,CVector& x_guess,int i_tol,double tol,int i_itmax,double i_thresh);
	cSparseLinSolve(const CMatrix& i_a, int i_tol,double tol,int i_itmax,double i_thresh);

	int itol;	// convergence criteria: 1,2,3,4
	double tol; // convergence tolerance
	int itmax;	// max number of iterations
	int iter;	// iterations needed
	double error; // estimated error
	double thresh; // threshold to set matrixelements in i_a to zero
			
	void solve(CVector& solution,int& niter,double& error);
	void solve(CMatrix& solution,const CMatrix& rhs,CVector& x_guess,int& niter,double& error);			
};			
	
	
typedef int (*p_dev_underdet_func) (const CVector& x, void* vp, int calc_for_gradient, CVector& f);


void gs_rootfind_underdetermined_solve(const CVector& initial_x, const CVector& bump_x, const CVector& tolerances, 
								   p_dev_underdet_func fn, void* vp, int max_tries, 
								   int max_restarts, const CMatrix& weights, CVector& found_x);




#endif

