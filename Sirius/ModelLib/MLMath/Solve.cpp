
#include "StdAfx.h"
#include "Solve.h"
#include <iostream>
#include <fstream>
#include "MlEqMaths.h"
 
using namespace std;

#define SUCCESS  0
 

/***************************************************************
**	Class   : LSGRGSolver 
**	Function: ConstraintFcn
**	Returns : nothing
**	Comment : 
****************************************************************/


void  LSGRGSolver::ConstraintFcn(double* contraintVals, double* xVals)
{
	double z = 0.0;

/*	for ( int i = 0 ; i < m_numberVariables; i++ )
	{
		 z += xVals[i];
	}
*/
	contraintVals[0] = z;
//  default implementation: no constraints
}

/***************************************************************
**	Class   : LSGRGSolver 
**	Function: ObjectiveFcn
**	Returns : nothing
**	Comment : 
****************************************************************/

void  LSGRGSolver::ObjectiveFcn(double* gVals,double* xVals)
{
	double z = 0.0;
	for ( int i = 0 ; i < m_numberVariables; i++ )
	{	
		z += xVals[i]*xVals[i];
	}

	*gVals = z;
}

/***************************************************************
**	Class   :  
**	Function: p_gcomp
**	Returns : nothing
**	Comment : private function which is passed downstream into rootsolver code
****************************************************************/

void LSGRGSolver::p_gcomp( double* gVals,double *xval)
{
	// 
	// A.G. This is now a member function
	//

	//LSGRGSolver* lsrc = (LSGRGSolver*)vp;
	//lsrc->ConstraintFcn( gVals+1,xval+1);
	//lsrc->ObjectiveFcn(&(gVals[lsrc->m_numberConstraints+1]),xval+1);

	ConstraintFcn( gVals+1,xval+1);
	ObjectiveFcn(&(gVals[m_numberConstraints+1]),xval+1);
}

/***************************************************************
**	Class   :  
**	Function: partialDerivs
**	Returns : nothing
**	Comment : private function which is passed downstream into rootsolver code
****************************************************************/

void LSGRGSolver::partialDerivs(double* x, long n, double* paij, long* iprow, long* ipcol, long* nnz)
{
	// 
	// A.G. This is now a member function
	//
	//LSGRGSolver* lsrc = (LSGRGSolver*)vp;

	CVector xVal(n, x);
	PartialDerivatives(xVal, m_partialDerivatives);

	paij	= m_partialDerivatives.getPtr()-1;
	iprow	= (long*)m_nonZeroJacobianMapRows.getPtr()-1;    
	ipcol	= (long*)m_nonZeroJacobianMapCols.getPtr()-1;
	*nnz	= m_numberOfNonZeroJacobians;

	//PartialDerivatives(xVal,lsrc->m_partialDerivatives);
	/*
	paij	= lsrc->m_partialDerivatives.getPtr()-1;
	iprow	= lsrc->m_nonZeroJacobianMapRows.getPtr()-1;    
	ipcol	= lsrc->m_nonZeroJacobianMapCols.getPtr()-1;
	*nnz	= lsrc->m_numberOfNonZeroJacobians;
	*/
}

/***************************************************************
**	Class   : LSGRGSolver 
**	Function: constructor
**	Returns : nothing
**	Comment : 
****************************************************************/


LSGRGSolver::LSGRGSolver()
:m_info( LsgrgInitialize() )
{}




/***************************************************************
**	Class   : LSGRGSolver 
**	Function: initialization
**	Returns : nothing
**	Comment : 
****************************************************************/

void LSGRGSolver::initialize(CVector& initialGuess,iVector& linearVariableIndex,
							  CMatrix& xValBounds,CMatrix& BoundsOnConstraints,
							  double InitialTolerance,double FinalTolerance,double StoppingTolerance,
							  CMatrix& NonZeroPartialDerivatives,int outputFlag)
{

//	NonZeroPartialDerivatives[numberVariables][numberConstraints]

 if ( BoundsOnConstraints.rows() == 0 )
 {
	BoundsOnConstraints.resize(1,2);
	BoundsOnConstraints[0][0]	=	-1e99;
	BoundsOnConstraints[0][1]	=	+1e99;
 }

 m_numberVariables		= initialGuess.getsize();
 m_numberConstraints	= BoundsOnConstraints.rows(); 

 m_linearVariableIndex	= linearVariableIndex;
 m_xValBounds			= xValBounds;
 m_BoundsOnConstraints  = BoundsOnConstraints;
 m_xInitialGuess		= initialGuess;

 m_InitialTolerance		= InitialTolerance;
 m_finalTolerance		= FinalTolerance;
 m_StoppingTolerance	= StoppingTolerance;

 // m_NonZeroPartialDerivatives[m_numberVariables][m_numberConstraints]
 // 1 if sensitive 0 else
 m_NonZeroPartialDerivatives	= NonZeroPartialDerivatives;
 m_partialDerivsAreSet						= m_NonZeroPartialDerivatives.rows();

 int i,j;
// assert(m_NonZeroPartialDerivatives.rows() == m_numberVariables);
// assert(m_NonZeroPartialDerivatives.cols() == m_numberConstraints);
 
 m_numberOfNonZeroJacobians = 0;
 for ( i = 0 ; i < m_numberVariables; i++ )
 {
  for ( j = 0 ; j < m_numberConstraints; j++ )
  {

	   if ( m_NonZeroPartialDerivatives.cols() != 0 )
	   {
		   if ( m_NonZeroPartialDerivatives[i][j] != 0.0 )
		   {
			m_numberOfNonZeroJacobians++;
		   }
	   }
  }
 }

 m_inverseJacobianMap.resize(m_numberVariables,m_numberConstraints);
 for ( i = 0 ; i < m_numberVariables; i++ )
 {
  for ( j = 0 ; j < m_numberConstraints; j++ )
  {
	  m_inverseJacobianMap[i][j] = -1.0;
  }
 }
 
 m_nonZeroJacobianMapRows.resize(m_numberOfNonZeroJacobians);
 m_nonZeroJacobianMapCols.resize(m_numberOfNonZeroJacobians);
 m_partialDerivatives.resize(m_numberOfNonZeroJacobians);

 int k = 0;
 for ( i = 0 ; i < m_numberVariables; i++ )
 {
  for ( j = 0 ; j < m_numberConstraints; j++ )
  {
	if ( m_NonZeroPartialDerivatives.cols() != 0 )
	{
	   if ( m_NonZeroPartialDerivatives[i][j] != 0)
	   {
		m_nonZeroJacobianMapRows[k] = i;
		m_nonZeroJacobianMapCols[k] = j;
		m_inverseJacobianMap[i][j]  = k;
		k++;
	   }
	}
  }
 }


  /* set feasibility tolerance epinit  */
   lsgrg_setparameter(m_info,"InitialTolerance", 0l, m_InitialTolerance);
  /* set feasibility tolerance epnewt  */
   lsgrg_setparameter(m_info,"epnewt", 0l, m_finalTolerance);
//   lsgrg_setparameter(info,"useph0", 00, 00);
  // big = lsgrg_get_plinfy(_info);  /* get default value of +infinity */
  /* set optimality  tolerance epstop  */
   lsgrg_setparameter(m_info,"epstop", 0l, m_StoppingTolerance);

/*  outputFlag=1; 
  FILE *grgoutfile;
  grgoutfile = fopen("grgoutfile.txt", "w");
  lsgrg_set_ioout(m_info,grgoutfile);

  lsgrg_set_inprnt(m_info,(int) outputFlag);
  lsgrg_set_otprnt(m_info,(int) outputFlag);
  lsgrg_set_printlevel(m_info,(int) outputFlag);
sos
*/

}


/***************************************************************
**	Class   : LSGRGSolver 
**	Function: setPartialDerivative
**	Returns : nothing
**	Comment : function sets partial derivatives into private data
****************************************************************/

void LSGRGSolver::setPartialDerivative(double setVal,int ixVar,int iConstraint)
{
  assert ( m_inverseJacobianMap[ixVar][iConstraint] != -1 );
 m_partialDerivatives[m_inverseJacobianMap[ixVar][iConstraint]] = setVal;
}


/***************************************************************
**	Class   : LSGRGSolver 
**	Function: setPartialDerivative
**	Returns : nothing
**	Comment : function retrieves partial derivatives from private data
****************************************************************/


double LSGRGSolver::getPartialDerivative(int ixVar,int iConstraint)
{
 if ( m_inverseJacobianMap[ixVar][iConstraint] == -1 )
  return 0.0;
 double z = m_partialDerivatives[m_inverseJacobianMap[ixVar][iConstraint]];
 return z;
}


/***************************************************************
**	Class   : LSGRGSolver 
**	Function: solve
**	Returns : nothing
**	Comment : function populates LSGRGSolverSolution
****************************************************************/

void LSGRGSolver::solve( CVector& x, int maximizeFlag, int& returnCode)
{
	 m_xInitialGuess	=	x;

	 int retVal,n;
	 int nfunctionsIn = 1+m_numberConstraints;//number o constraints+ objective fcn
	 int nobjIn = nfunctionsIn;//sosm_numberConstraints+1;

	 int nlinearVar = 0;
	 for ( n = 0 ; n < m_linearVariableIndex.getsize(); n++ )
	 {
		if ( m_linearVariableIndex[n] != 0 ){
			nlinearVar++;
		}
	 }

//	 m_linearVariableIndex.resize(0);
	 
	 iVector linearVars(nlinearVar+1);
	 linearVars[0] = nlinearVar;
	 for ( n = 0 ; n < m_linearVariableIndex.getsize(); n++ )
	 {
		if ( m_linearVariableIndex[n] != 0 )
		{
		  linearVars[n+1] = n+1;
		}
	 }
 
	 CVector blVar(m_numberVariables);
	 CVector buVar(m_numberVariables);
	 for ( n = 0 ; n < m_numberVariables; n++ )
	 {
	  blVar[n] = m_xValBounds[n][0]; 
	  buVar[n] = m_xValBounds[n][1]; 
	 }
	 CVector blCon(m_numberConstraints);
	 CVector buCon(m_numberConstraints);
	 for ( n = 0 ; n < m_numberConstraints; n++ )
	 {
	  blCon[n] = m_BoundsOnConstraints[n][0]; 
	  buCon[n] = m_BoundsOnConstraints[n][1]; 
	 }

	 m_LSGRGSolverSolution.m_finalFunctionValues.resize(nfunctionsIn+1);
	 m_LSGRGSolverSolution.m_finalLagrangeMultiplyers.resize(nfunctionsIn+1);
	 m_LSGRGSolverSolution.m_finalReducedGradients.resize(m_numberVariables+1);
	 
	 long nbind;//              nbr of binding constraints
	 long nnonb;//              nbr of nonbasic variables
	 iVector indexOfBindingConstraints(nfunctionsIn+1);//m_numberConstraints);so
	 iVector indexOfNonBasicVars(m_numberVariables);

	 //void (*partialderivs)(double* x, int n, double* paij, int* iprow, int* ipcol, int *nnz,void* vp);
	 
	 P_PARSH pfnparsh = NULL;
	 if(m_partialDerivsAreSet == 0)
	 {
		//partialderivs = NULL ;
		 pfnparsh = NULL;
	 }
	 else
	 {
		//partialderivs = partialDerivs;
		pfnparsh = LSGRGSolver::partialDerivs;
	 }

	 m_LSGRGSolverSolution.m_finalxValues = m_xInitialGuess;

	 /*
	 retVal =  grgsub(m_info,m_numberVariables, nfunctionsIn,nobjIn,maximizeFlag,
		 linearVars.getPtr() ,
		 blVar.getPtr()-1,buVar.getPtr()-1,blCon.getPtr()-1,buCon.getPtr()-1,
		 m_LSGRGSolverSolution.m_finalxValues.getPtr()-1,m_LSGRGSolverSolution.m_finalFunctionValues.getPtr()-1,m_LSGRGSolverSolution.m_finalLagrangeMultiplyers.getPtr()-1,
		 indexOfNonBasicVars.getPtr()-1,
		 m_LSGRGSolverSolution.m_finalReducedGradients.getPtr()-1,indexOfBindingConstraints.getPtr()-1,&nbind,&nnonb,
		 p_gcomp,partialDerivs, m_numberOfNonZeroJacobians,this);
	*/

	//
	// A.G. We need to do this in order to call the member function
	//
	P_GCOMP pfncomp = LSGRGSolver::p_gcomp;	// just to be explicit
	m_info->m_pOwner = this;

	retVal = grgsub(m_info, 
					m_numberVariables, 
					nfunctionsIn, 
					nobjIn,maximizeFlag,
					(long*)linearVars.getPtr(),
					blVar.getPtr()-1, 
					buVar.getPtr()-1, 
					blCon.getPtr()-1, 
					buCon.getPtr()-1,
					m_LSGRGSolverSolution.m_finalxValues.getPtr()-1, 
					m_LSGRGSolverSolution.m_finalFunctionValues.getPtr()-1,
					m_LSGRGSolverSolution.m_finalLagrangeMultiplyers.getPtr()-1,
					(long*)indexOfNonBasicVars.getPtr()-1,
					m_LSGRGSolverSolution.m_finalReducedGradients.getPtr()-1,
					(long*)indexOfBindingConstraints.getPtr()-1,
					&nbind,
					&nnonb,
					pfncomp,
					pfnparsh, 
					m_numberOfNonZeroJacobians);

	 x = m_LSGRGSolverSolution.m_finalxValues;
 
	 returnCode = retVal;
}

static const double EPS = 1.0e-07;

void integratorFcn(double x, double y[], double dydx[],void *vp)
{
  RungeKuttaIntegrator* rkutta = (RungeKuttaIntegrator*)vp;
  rkutta->integrands(rkutta->m_res,x);
  dydx = rkutta->m_res.getPtr();
}


void RungeKuttaIntegrator::integrate(CVector& integral)
{ 		
  int nok,nbad,kount,n,k;
  double dxsav = 0.0;
 	
  int nvar = integral.getsize(); 	
  CVector ystart(nvar+1);

  for ( n = 0 ; n <= nvar; n++ )
  	ystart[n] = 0.0;
  	
  for ( n = 0 ; n < nvar; n++ )
	integral[n] = 0.0;
	  
  int nlimits = m_IntegrationLimits.getsize()-1;
  if ( m_IntegrationLimits[0] >= m_IntegrationLimits[nlimits] )
  {			
  	  return;
  }		
  for ( n = 0 ; n < nlimits; n++ )
  {		
	  for ( k = 0 ; k <= nvar; k++ )
		ystart[k] = 0.0;
		
  	  MlEqMaths::odeint(ystart,nvar,m_IntegrationLimits[n],m_IntegrationLimits[n]*(1.0+m_bump+EPS),1e-8,m_h_init,m_h_min,&nok,&nbad,integratorFcn,MlEqMaths::rkqs,NULL,NULL,kount,(void *) &integratorFcn,0,dxsav,m_scale);
	  for ( k = 0 ; k < nvar; k++ )
		integral[k] += ystart[k+1];


	  for ( k = 1 ; k <= nvar; k++ )
		ystart[k] = 0.0;

  	  MlEqMaths::odeint(ystart,nvar,m_IntegrationLimits[n+1]*(1.0-m_bump+EPS),m_IntegrationLimits[n+1],1e-8,m_h_init,m_h_min,&nok,&nbad,integratorFcn,MlEqMaths::rkqs,NULL,NULL,kount,(void *) &integratorFcn,0,dxsav,m_scale);
	  for ( k = 0 ; k < nvar; k++ )
		integral[k] += ystart[k+1];

	  for ( k = 1 ; k <= nvar; k++ )
		ystart[k] = 0.0;

  	  MlEqMaths::odeint(ystart,nvar,m_IntegrationLimits[n]*(1.0+m_bump+EPS),m_IntegrationLimits[n+1]*(1.0-m_bump+EPS),m_eps,m_h_init,m_h_min,&nok,&nbad,integratorFcn,MlEqMaths::rkqs,NULL,NULL,kount,(void *) &integratorFcn,0,dxsav,m_scale);
	  for ( k = 0 ; k < nvar; k++ )
		integral[k] += ystart[k+1];
  }		
} 		


void RungeKuttaIntegrator::initialize(int nfunctions,CVector& integrationLimits,double eps,double hInit,double hMin,double bump)
{
	m_res.resize(nfunctions);

	m_h_min				=	hMin;
	m_eps				=	eps;
	m_bump				=	bump;
	m_h_init			=	hInit;
	m_scale				=	10.0;
	m_IntegrationLimits =	integrationLimits;
}


double EPSILON = 1e-20;
	
void cSparseLinSolve::solve(CVector& solution,int& niter,double& xerror)
{	
	linbcg(b);
	niter = iter;
	xerror = error;
		
	if ( solution.getsize() != ndim )
		throw(" dimesion mismatch in CVector solution "); 
	
	for ( int n = 0 ; n < ndim; n++ )
		solution[n] = x[n+1];
}	
	
	
void cSparseLinSolve::solve(CMatrix& solution,const CMatrix& rhs,CVector& x_guess,int& niter,double& error)
{	
	int neq = rhs.cols();
	niter;error; // crem
	int i,j;
	b.resize(rhs.rows()+1);
	x.resize(ndim+1);


	if ( x_guess.getsize() != 0 )
	   for ( i = 0 ; i < ndim; i++ )
		  x[i+1] = x_guess[i];
	
	for ( i = 0 ; i < rhs.cols(); i++ )
	{
		for ( j = 0 ; j < ndim; j++ )
			b[j+1] = rhs[j][i];
		
		linbcg(b);		
		for ( j = 0 ; j < ndim; j++ )
			solution[j][i] = x[j+1];		
	}	
}		
		
		
		
cSparseLinSolve::cSparseLinSolve( const CMatrix& i_a,int i_itol,double i_tol,int i_itmax,double i_thresh)
{	
	a = &i_a;
	ndim = i_a.rows(); 

	
	if ( i_a.cols() != ndim )
		throw("input matrix must be a square matrix ");
	
	b.resize(ndim+1);
	x.resize(ndim+1);
	
	itmax = i_itmax;
	thresh = i_thresh;
	itol = i_itol;
	tol = i_tol;
	sprsin();
}	
	
	
cSparseLinSolve::cSparseLinSolve(const CMatrix& i_a, CVector &i_b,CVector& x_guess,int i_itol,double i_tol,int i_itmax,double i_thresh)
{	
	
	a = &i_a;
	ndim = i_a.rows();
	
	if ( i_a.cols() != ndim )
		throw("input matrix must be a square matrix ");
	
	b.resize(ndim+1);
	x.resize(ndim+1);
	
	for (int n = 0 ; n < ndim; n++ )
	{
		b[n+1] = i_b[n];
		x[n+1] = x_guess[n];	
	}
	
	itmax = i_itmax;

	thresh = i_thresh;
	itol = i_itol;
	tol = i_tol;
	sprsin();
}	
	
	
	
	
void cSparseLinSolve::linbcg(CVector& b)
{	
	int j;
	double ak,akden,bk,bkden = 0.0,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm;
	int n = ndim;
	int m = n+1;
	CVector p(m),pp(m),r(m),rr(m),z(m),zz(m);
	
	iter = 0;
	atimes(n,x,r,0);
	
	for ( j = 1 ; j <= n ; j++ )
	{
		r[j] = b[j] - r[j];
		rr[j] = r[j];
	}
	// atimes(n,r,rr,0);
	
	znrm = 1.0;
	if ( itol == 1 )
		bnrm = snrm(n,b,itol);
	else if ( itol == 2 )
	{
		asolve(n,b,z,0);
	
		bnrm = snrm(n,z,itol);
	}
	else if ( itol == 3 || itol == 4 )
	{
		asolve(n,b,z,0);
		bnrm = snrm(n,z,itol);
		asolve(n,r,z,0);
		znrm = snrm(n,z,itol);
	}
	else
		throw("itol switch not found");
	
	asolve(n,r,z,0);
	while ( iter <= itmax )
	{	
		++iter;
		zm1nrm = znrm;
		asolve(n,rr,zz,1);
		bknum;//crem
		for ( bknum = 0.0, j = 1 ; j <= n ; j++ )
			bknum += z[j]*rr[j];
		
		if ( iter == 1 )
		{
			for ( j = 1 ; j <= n ; j++ )
			{
				p[j] = z[j];
				pp[j] = zz[j];
			}
		}
		else
		{
			bk = bknum/bkden;
			for ( j = 1 ; j <= n ; j++ )
			{
				p[j] = bk*p[j] + z[j];
				pp[j] = bk*pp[j] + zz[j];
			}
		}

		bkden = bknum;
		atimes(n,p,z,0);
		for ( akden = 0.0, j = 1; j <= n; j++ )
			akden += z[j]*pp[j];

		ak = bknum/akden;
		atimes(n,pp,zz,1);
		for ( j = 1 ; j <= n ; j++ )
		{
			x[j] += ak*p[j];
			r[j] -= ak*z[j];
			rr[j] -= ak*zz[j];
		}

		asolve(n,r,z,0);
		if ( itol == 1 || itol == 2)
		{
			znrm = 1.0;
			error = snrm(n,r,itol)/bnrm;
		}
		else if ( itol == 3 || itol == 4)
		{
			znrm = snrm(n,z,itol);
			if ( fabs(zm1nrm-znrm) > EPSILON*znrm )
			{

				dxnrm = fabs(ak)*snrm(n,p,itol);
				error = znrm/fabs(zm1nrm-znrm)*dxnrm;
			}
			else
			{
				error = znrm/bnrm;
				continue;
			}
			xnrm = snrm(n,x,itol);
			if ( error <=0.5*xnrm )
				error /= xnrm;
			else
			{
				error = znrm/bnrm;
				continue;
			}
		}
		if ( error <= tol ) break;
	}
	if (iter > itmax )
		throw(" maximum number of iterations exceeded in linbcg ");	
}		
		
		
double cSparseLinSolve::snrm(int n,CVector & sx,int itol)
{
	int  i, isamax;
	double ans;
	if ( itol <= 3 )
	{
		ans = 0.0;
		for ( i = 1 ; i <= n ; i++ )
			ans += sx[i]*sx[i];
		return sqrt(ans);
	}
	else
	{
		isamax = 1;
		for ( i = 1 ; i <= n ; i++ )
		{
			if ( fabs(sx[i] ) > fabs(sx[isamax]))
				isamax = i ;
		}
		return fabs(sx[isamax]);
	}
}	
	
	
	
void cSparseLinSolve::atimes(int n,CVector & x, CVector & r, int itrnsp)
{		
	if ( itrnsp )
		dsprstx(sa,ija,x,r,n);
	else
		dsprsax(sa,ija,x,r,n);
}	
	
	
void cSparseLinSolve::asolve(int n,CVector& b,CVector& x,int itrnsp)
{	
  itrnsp; // crem
	int i;
	for ( i = 1 ; i <= n; i++ )
		x[i] = (sa[i] != 0.0 ? b[i]/sa[i] : b[i] );
}	
	
	
void cSparseLinSolve::dsprstx(CVector& sa,iVector& ija,CVector & x,CVector& b,int n)
{	
	int i,k,j;
	if ( ija[1] != n+2 )
		throw(" dsprstx: mismatch vector and matrix ");
	for ( i = 1 ; i <= n ; i++ )
		b[i] = sa[i]*x[i];
	
	for ( i = 1 ; i <= n ; i++ )
	{
		for ( k = ija[i]; k<= ija[i+1]-1 ; k++ )
		{
			j = ija[k];
			b[j] += sa[k]*x[i];
		}
	}
}	
	
	
	
void cSparseLinSolve::sprsin()
{	
	int i,j,k;
	int n = ndim;
	nmax = 0;
	for ( i = 0 ; i < n ; i++ )
		for ( j = 0 ; j < n ; j++ )
			if( ( fabs((*a)[i][j]) > EPSILON )&& ( i != j ) )
				nmax++;
	
	nmax += n + 2;
	sa.resize(nmax);
	ija.resize(nmax);
	for ( j = 1; j <= n; j++ )
		sa[j] = (*a)[j-1][j-1];
	
	ija[1] = n+2;
	k = n + 1;		
	for ( i = 1 ; i <= n; i++ )
	{
		for ( j = 1 ; j <=n; j++)
		{
			if ( fabs((*a)[i-1][j-1]) >= thresh && i!=j )
			{
				if (++k >nmax)
					throw(" nmax not big enough");
				sa[k]= (*a)[i-1][j-1];
				ija[k] = j;
			}
		}
		ija[i+1] = k+1;
	}
}	
	
	
	
	
void cSparseLinSolve::dsprsax(CVector& sa, iVector& ja,CVector& x,CVector& b,int n)
{	
	int i,k;
	ja;//crem
	if ( ija[1] != n + 2 )
		throw("mismatched vector and matrix in dsprax");
	
	for ( i = 1 ; i <= n ; i++ )
	{
		b[i] = sa[i]*x[i];
		for ( k = ija[i] ; k <= ija[i+1]-1; k++ )
			b[i] += sa[k]*x[ija[k]];
	}
}	



/****************************************************************************/



class underdetermined_fn_info
{
public:
	int n_iterations;
	int n_restarts;
	int n_functions;
	CVector tolerances;
	p_dev_underdet_func fn;
	void* user;
};

/* wrapper for the user's underdetermined rootfind function */


static int underdetermined_fn (const CVector& x, underdetermined_fn_info* mine, 
											 int calc_for_grad, CVector& f)
{	
	int i, status;
	
	if (calc_for_grad & eUNDERDETERMINED_SETUP_GRADIENT)
	{
		mine->n_restarts--;
		if (mine->n_restarts < 0)
		{
			throw("exceeded maximum allowed number of gradient evaluations");
			return ERROR;
		}
	}
	else if (calc_for_grad == eUNDERDETERMINED_CALC_BEST)
	{
		mine->n_iterations--;
		if (mine->n_iterations < 0)
		{
		  throw(" exceeded maximum allowed number of function evaluations");
		  return ERROR;
		}
	}
	
	status = mine->fn(x, mine->user, calc_for_grad, f);
	if (status != SUCCESS)
	{
		throw(" error during underdetermined root search");
		return ERROR;
	}
	
	/* we are only concerned with the function value in units of the tolerance */
	for (i = 0; i < mine->n_functions; i++)
		f[i] /= mine->tolerances[i];
	return SUCCESS;
}	

	
	

static int rootfind_underdetermined_backtrack(underdetermined_fn_info* info, CVector& current_x, 
											  CVector step_x, const CVector& stored_f, 
											  CVector& working_f, double overshoot_fraction, 
											  double overshoot_tolerance, double restart_tolerance, 
											  int j_is_approximate)
{
	double old_dot_old, new_dot_new, old_dot_new, change_squared;
	double this_length, min_length = 0.0, max_length = 1.0;
	int flipped = FALSE;
	double min_frac = -1.0;//-restart_tolerance;
	double max_frac = overshoot_fraction / sqrt(1.0 - MlEqMaths::dsqr(restart_tolerance) + MlEqMaths::dsqr(overshoot_fraction));
	/* this gives a sensible first guess for the step length to try; see its next use */
	int n_variables = current_x.getsize();
	int n_functions = stored_f.getsize();
	int i, status;

	CVector temp_x(n_variables);
	
	while (fabs(overshoot_fraction) > overshoot_tolerance)
	{
		this_length = min_length - (max_length - min_length) * min_frac/(max_frac-min_frac);

//			MlEqMaths::Max(0.2, MlEqMaths::Min(0.8, -min_frac / (max_frac - min_frac)));

		for (i = 0; i < n_variables; i++)
			temp_x[i] = current_x[i] + this_length * step_x[i];
		status = underdetermined_fn(temp_x, info, eUNDERDETERMINED_CALC_BEST, working_f);
		if (status != SUCCESS)
		{
			return status;
		}

		for (i = 0, old_dot_old = 0.0, old_dot_new = 0.0, new_dot_new = 0.0; i < n_functions; i++)
		{
			old_dot_old += MlEqMaths::dsqr(stored_f[i]);
			new_dot_new += MlEqMaths::dsqr(working_f[i]);
			old_dot_new += stored_f[i] * working_f[i];
		}
		change_squared = old_dot_old + new_dot_new - 2.0 * old_dot_new;
		if (MlEqMaths::deqz(change_squared))	
			overshoot_fraction = 0.0;
		else
			overshoot_fraction = (new_dot_new - old_dot_new) / change_squared;
	
		if (overshoot_fraction > 0.0)
		{
			/* we are still overshooting */
			max_length = this_length;
			max_frac = overshoot_fraction / ( 1.0 - MlEqMaths::dsqr(restart_tolerance) + MlEqMaths::dsqr(overshoot_fraction));
		}
		else
		{
			/* we are now undershooting */
			if (this_length > 1.0 - overshoot_tolerance && flipped)
				/* improving in the opposite direction; get new Jacobian */
				overshoot_fraction = 0.0;  /* this breaks out of the while loop */

			min_length = this_length;
			min_frac = overshoot_fraction / MlEqMaths::Min(1.0, 1.0 - MlEqMaths::dsqr(restart_tolerance) + MlEqMaths::dsqr(overshoot_fraction));
		}

		if (this_length < overshoot_tolerance && overshoot_fraction > 1.0)
		{
			/* we have tried a short step, but we are still going uphill */
			if (flipped)
			{
			
				if (j_is_approximate)
				{
					/* break out and get a new Jacobian with current_x unchanged */
	//				more efficient way ? n_functions;
					working_f = stored_f;

					return SUCCESS;
				}
				else
				{
					//assert(!j_is_approximate);
//					throw("Root not found");
					return ERROR;
				}
			}

			/* (else) we can't get a better Jacobian; linesearch backwards */
			for (i = 0; i < n_variables; i++)
				step_x[i] = -step_x[i];
	
			min_length = 0.0;
			max_length = 1.0;
			min_frac = -1.0;
			max_frac = 1.0;
			flipped = TRUE;
		}
	
/*		if (min_length > 0.0 && max_length - min_length < overshoot_tolerance)
			// close enough; break out of the loop //
			overshoot_fraction = 0.0;
*/
	}
	/* now we have completed our linesearch */
	

//	more efficient way ? n_variables;

	current_x = temp_x;

	return SUCCESS;
}	
	
/****************************************************************************/
	
	
void rootfind_underdetermined_qp_step(const CMatrix& jacobian, CMatrix& j_transpose, CMatrix& j_effective, 
									 CMatrix& j_decomposed, const CVector& stored_f, 
									 CVector& working_f, const CMatrix& weights, CVector& step_x,
									 cSparseLinSolve&  sparse)
{	
	int n_variables, n_functions;
	int i;
	n_functions = jacobian.rows();
	n_variables = jacobian.cols();
		
	/* negate f */
	for (i = 0; i < n_functions; i++){
		working_f[i] = -stored_f[i];
	}
	
	int n,m,k;
	for ( n = 0 ; n < jacobian.rows() ; n++ ){
		for (  m = 0 ; m < jacobian.cols() ; m++ ){
			j_transpose[m][n] = jacobian[n][m];
		}
	}
	
	int niter;
	double error;
	CVector x_guess;


	if ( ( weights.rows() != 0 ) || (weights.cols() != 0 ) ){

		CMatrix temp(j_transpose.rows(),j_transpose.cols());
		sparse.solve(temp,j_transpose,x_guess,niter,error);
		j_transpose = temp;
//		sparse.solve(j_transpose,j_transpose,x_guess,niter,error);

	}

	
//	at this point j_transpose contains W^{-1} J^T	

	for ( n = 0 ; n < j_effective.rows() ; n++ ){
		for ( m = 0 ; m < j_effective.cols() ; m++ ){
				j_effective[n][m] = 0.0;
	}}
		
	for ( n = 0 ; n < jacobian.rows() ; n++ ){
		for ( m = 0 ; m < j_transpose.cols() ; m++ ){
			for ( k = 0 ; k < jacobian.cols() ; k++ ){
				j_effective[n][m] += jacobian[n][k]*j_transpose[k][m];
	}}}
	
		
	CVector b(j_effective.cols());
	CVector result;  
	CVector residuals;
	double tol = sparse.tol;
	
	CMatrix genInverse;


	GVector<int> indx;
	double d=-1e30;

	MlEqMaths::ludcmp(j_effective,indx,d);
	MlEqMaths::lubksb(j_effective, indx, working_f);

//	lubksb(CMatrix& a, GVector<int>& indx, CVector& b);

//	invertMatrix(genInverse,j_effective);
    
//	CMatrix temp;
//	temp = j_transpose;

//	j_transpose = j_transpose*genInverse;

/*	int j;

	for (  i = 0 ; i < j_transpose.rows() ; i++ ){
		for ( j = 0 ; j < j_transpose.cols() ; j++ ){
			j_transpose[i][j] = 0.0;
		}
	}

	for (  i = 0 ; i < j_transpose.rows() ; i++ ){
		for ( j = 0 ; j < j_transpose.cols() ; j++ ){
			for ( k = 0 ; k < j_transpose.cols(); k++){
				j_transpose[i][j] += temp[i][k]*genInverse[k][j]; 
		}
	}}
*/

	for (  n = 0 ; n < j_transpose.rows() ; n++ ){
		step_x[n]=0.0;
	}

	for (  n = 0 ; n < j_transpose.rows() ; n++ ){
		for ( m = 0 ; m < j_transpose.cols() ; m++ ){
			step_x[n] += j_transpose[n][m]*working_f[m];
	}}


	/* now step_x contains the QP step */ 

}		

	
static void rootfind_underdetermined_solve_ex(const CVector& initial_x,const  CVector& bump_x, 
											 underdetermined_fn_info* info, const CMatrix& weights, 
											 CMatrix& jacobian, CMatrix& j_transpose, CMatrix& j_effective, 
											 CMatrix& j_decomposed, CVector& current_x, 
											 CVector& step_x, CVector stored_f, CVector& working_f, 
											 CVector& last_good_f, CVector& found_x)
{		
	static const double overshoot_tolerance = 0.15;
	static const double restart_tolerance = 0.45;//1.6;//sos0.45;
	int i, j, status;
	int n_variables = initial_x.getsize();
	int n_functions = info->n_functions;
	int j_is_approximate, restart_now = FALSE, finished;
	double new_dot_new, old_dot_new, old_dot_old;
	double step_size, change_squared, overshoot_fraction;

	int itol = 1;
	double tol = 1e-6;
	int itmax = 300;//50;
	double thresh = EPSILON;

	cSparseLinSolve sparse(weights,itol,tol,itmax,thresh);

	int i_iter = 0;
	for ( ; ; )
	{	
		/* set up the Jacobian */
		status = underdetermined_fn(current_x, info, 
									eUNDERDETERMINED_CALC_GRADIENT|eUNDERDETERMINED_SETUP_GRADIENT, 
									stored_f);
		
		for (i = 0; i < n_variables; i++)
		{
			const double base_x_i = current_x[i];
			double bump_x_i;
			current_x[i] += bump_x[i];
			underdetermined_fn(current_x, info, eUNDERDETERMINED_CALC_GRADIENT|i, working_f);
			
			bump_x_i = current_x[i] - base_x_i;
			for (j = 0; j < n_functions; j++)
				jacobian[j][i] = (working_f[j] - stored_f[j]) / bump_x_i;

			current_x[i] = base_x_i;
		}

		/* get the true function value */
		if (!restart_now)
		{
			/* we are here for the first time and don't know f */
			underdetermined_fn(current_x, info, eUNDERDETERMINED_CALC_BEST, stored_f);
			
			/* test convergence */
			finished = TRUE;
			for (i = 0; i < n_functions && finished; i++)
				/* recall that working_f[i] has already been normalized to tolerances[i] */
				if (fabs(stored_f[i]) > 1.0)
					finished = FALSE;
			if (finished)
			{
				found_x = current_x;
				return ;
			}
		}
		else
//			more efficient way ? n_functions;
			stored_f = last_good_f;

		j_is_approximate = FALSE;
		restart_now = FALSE;

		/* iterate QP step */
		do
		{

			i_iter++;

			rootfind_underdetermined_qp_step(jacobian, j_transpose, j_effective, j_decomposed,
													  stored_f, working_f, weights, step_x,sparse);


			for (i = 0, step_size = 0.0; i < n_variables; i++)
			{
				step_size += step_x[i] * step_x[i];
				current_x[i] += step_x[i];
			}
			status = underdetermined_fn(current_x, info, eUNDERDETERMINED_CALC_BEST, working_f);

			/* monitor to prevent overshoot */
			for (i = 0, old_dot_old = 0.0, old_dot_new = 0.0, new_dot_new = 0.0; i < n_functions; i++)
			{
				old_dot_old += MlEqMaths::dsqr(stored_f[i]);
				new_dot_new += MlEqMaths::dsqr(working_f[i]);
				old_dot_new += stored_f[i] * working_f[i];
			}
			change_squared = old_dot_old + new_dot_new - 2.0 * old_dot_new;
			if (!MlEqMaths::deqz(change_squared))	
			{
				overshoot_fraction = (new_dot_new - old_dot_new) / change_squared;
				if (overshoot_fraction > restart_tolerance)
				{
					/* we are going backwards; linesearch for an acceptable step */
					for (i = 0; i < n_variables; i++)
						current_x[i] -= step_x[i];
					rootfind_underdetermined_backtrack(info, current_x, step_x, stored_f, 
																working_f, overshoot_fraction, 
																overshoot_tolerance, restart_tolerance, 
																j_is_approximate);
					/* this updates current_x and working_f */
//					more efficient way ? n_functions;
					last_good_f = working_f;

					restart_now = TRUE;
				/* end of case where restart is necessary */
				}
				else if (overshoot_fraction > overshoot_tolerance)
				{
					double rollback = overshoot_fraction * MlEqMaths::Min(overshoot_fraction, overshoot_tolerance) / overshoot_tolerance;
					double temp_i;
					for (i = 0; i < n_variables; i++)
					{
						temp_i = rollback * step_x[i];
						step_x[i] -= temp_i;
						current_x[i] -= temp_i;
					}
					step_size *= MlEqMaths::dsqr(1.0 - rollback);
					underdetermined_fn(current_x, info, FALSE, working_f);
					
				}
			}
			else
			{
				/* the function value is completely unchanged */
				if (j_is_approximate)
				{
//					more efficient way ? n_functions;
					last_good_f = working_f;

					restart_now = TRUE;
				}
				else
				{
					throw("Root not found");
				}
				
			}

			/* test convergence */
			finished = TRUE;
			for (i = 0; i < n_functions && finished; i++)
				/* recall that working_f[i] has already been normalized to tolerances[i] */
				if (fabs(working_f[i]) > 1.0)
					finished = FALSE;
			if (finished)
			{
				found_x = current_x;
				return;
			}

			if (!restart_now)
			{
				double swap;
				/* Broyden's update to Jacobian */
				for (i = 0; i < n_functions; i++)
				{
					swap = working_f[i];
					working_f[i] -= stored_f[i];
					stored_f[i] = swap;
					/* now stored_f contains the latest point, working_f contains the recent change */
					for (j = 0; j < n_variables; j++)
						working_f[i] -= jacobian[i][j] * step_x[j];
				}
				/* now working_f contains y-Js; here comes Broyden */
				for (i = 0; i < n_functions; i++)
					for (j = 0; j < n_variables; j++)
						jacobian[i][j] += working_f[i] * step_x[j] / step_size;
				j_is_approximate = TRUE;
			}
		}
		while (!restart_now);
	}
	throw("Root not found");
}	
	
	
void rootfind_underdetermined_solve(const CVector& initial_x, const CVector& bump_x, const CVector& tolerances, 
								   p_dev_underdet_func fn, void* vp, int max_tries, 
								   int max_restarts, const CMatrix& weights, CVector& found_x)
{	
	CMatrix jacobian;
	CMatrix j_effective;
	CVector current_x;
	CVector step_x;
	CVector stored_f;
	CVector working_f;
	CVector last_good_f;
	CMatrix j_transpose;
	CMatrix  j_decomposed;

	found_x.resize(initial_x.getsize());
	
	int n_variables = initial_x.getsize();
	int n_functions = tolerances.getsize();;
	underdetermined_fn_info info;
		
	
	jacobian.resize(n_functions, n_variables);
	j_transpose.resize(n_variables, n_functions);
	j_effective.resize(n_functions, n_functions);
	j_decomposed.resize(n_functions,n_functions);
	
	current_x.resize(initial_x.getsize());
	current_x = initial_x;
	
	step_x.resize(n_variables);
	stored_f.resize(n_functions);
	working_f.resize(n_functions);
	last_good_f.resize(n_functions);
	info.tolerances.resize(tolerances.getsize());
	(info.tolerances) = tolerances;
		
	info.n_iterations = max_tries;
	
	info.n_restarts = max_restarts;
	info.n_functions = n_functions;
	info.fn = fn;
	info.user = vp;
	
	rootfind_underdetermined_solve_ex(initial_x, bump_x, &info, weights, jacobian, j_transpose, 
											   j_effective, j_decomposed, current_x, step_x, stored_f, 
											   working_f, last_good_f, found_x);
}	
	
	



