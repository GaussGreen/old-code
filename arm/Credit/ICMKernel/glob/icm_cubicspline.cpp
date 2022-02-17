#include "ARMKernel\glob\firsttoinc.h" 
#include <nag.h>
#include "ICMKernel\glob\icm_enums.h"
#include "ICMKernel\util\icm_macro.h"
#include "ICMKernel\glob\icm_cubicspline.h"
#include "ICMKernel\util\icm_brentsolver.h"
#include "ICMKernel\util\icm_utils.h"


void CCubicspline::build(const double * xvalues,const double * yvalues,int size,eEMMODE emmode,
						 eSPLINE_EMMODE condition, double leftSlope,double rightSlope,
						 vector<double>& alpha,bool calibrate)
{
	// destruction de this
	destroy();

	m_emmode	= emmode;
	m_alpha		= alpha;
	m_calibrate = calibrate;
	m_leftSlope = leftSlope;
	m_rightSlope = rightSlope;

	if((m_size = size) == 0) return;

	m_xvalues	= new double [size];
	m_yvalues	= new double [size];

	for(int k = 0; k < m_size; m_xvalues[k] = xvalues[k], m_yvalues[k] = yvalues[k++]);

	m_Condition = condition;

	buildSpline();
}

double CCubicspline::EvaluateDerivSecond(const double& x)
{
	vector<double> alpha = m_alpha;
	//alpha.resize(m_cal_size+1);
	alpha[m_cal_size]=x;

	CCubicspline c(m_xvalues,m_yvalues,m_size,m_emmode,m_Condition,m_leftSlope,m_rightSlope,
					alpha,false);
	double value = 0.;
	value = c.m_bk[m_cal_size];

	return (value);
}

void CCubicspline::buildSpline()
{
	if(m_size < 3)
	{
		destroy();

		return;
	}

	m_needtorebuild = false;

	switch (m_Condition)
	{
		case _emC1:
			buildSplineWithC1Condition();
			return;
		case _emC2:
			buildSplineWithC2Condition();
			return;
		case _emC1_alpha:
			buildSplineOptimized();
			return;
		case _emOptimized:
			buildSplineOptimizedMulti();
			return;
		case _emOrder4:
			buildSplineOrder4();
			return;
		case _emOrder4Convex:
			buildSplineOrder4CondDerivs();
			return;
		case _emOptimizedC4:
			buildSplineOptimizedMultiC4();
			return;
		case _emShumaker:
			buildSplineShumaker();
			return;
		case _emNone:
		default:;
	}

	double * fk = new double [m_size], * hk = new double [m_size];

	for(int k = 0; k < m_size - 1; k++) hk[k] = m_xvalues[k+1] - m_xvalues[k];
	for(k = 1; k < m_size - 1; k++) fk[k] = 6.*((m_yvalues[k+1] - m_yvalues[k]) / hk[k] - (m_yvalues[k] - m_yvalues[k-1]) / hk[k-1]);

	fk[1]			-= m_leftSlope * hk[0];
	fk[m_size-2]	-= m_rightSlope * hk[m_size-2];

	findMk(fk,hk);

	if(m_ak != 0) delete [] m_ak;
	if(m_bk != 0) delete [] m_bk;
	m_ak = new double [m_size];
	m_bk = new double [m_size];

	for(k = 0;k < m_size-1; k++) m_ak[k] = (m_yvalues[k+1] - m_yvalues[k]) / hk[k] - hk[k] * (m_mk[k+1] - m_mk[k]) / 6.;
	for(k = 0;k < m_size; k++) m_bk[k] = m_yvalues[k] - m_mk[k] * hk[k] * hk[k] / 6.;

	delete [] hk;
	delete [] fk;
}

void CCubicspline::findMk(double * fk, double * hk)
{
	if(m_mk != 0) delete [] m_mk;
	int k, n = m_size-1;

	m_mk = new double [m_size];

	// résolution par pivot de gauss
	double * beta = new double [m_size] , * gprime = new double [m_size];

	beta[n-1]	= 2. * (hk[n-1] - hk[n-2]);
	gprime[n-1] = fk[n-1];

	for(k = n-2; k > 0; k--)
	{
		beta[k]		= 2. * (hk[k] + hk[k-1]) - hk[k] * hk[k-1] / beta[k+1];
		gprime[k]	= fk[k] - hk[k] * gprime[k+1] / beta[k+1];		
	}

	m_mk[1] = gprime[1] / beta[1];
	for(k = 2; k < n; k++) m_mk[k] = (gprime[k] - hk[k-1] * m_mk[k-1]) / beta[k];

	m_mk[0] = 0.;
	m_mk[n] = 0.;

	delete [] beta;
	delete [] gprime;
}

void CCubicspline::buildSplineOptimized()
{
	delete [] m_mk;
	delete [] m_ak;
	delete [] m_bk;
	delete [] m_ck;
	delete [] m_dk;

	m_mk = 0;
	
	int k, n = m_size - 1;

	m_ak	= new double [n];
	m_bk	= new double [n];
	m_ck	= new double [n+1];
	m_dk	= new double [n];

	double * delta = new double [n];

	// les intervalles
	for(k = 0; k < n; k++) delta[k] = m_xvalues[k+1] - m_xvalues[k];

	// la valeur des polynomes aux points (xk)
	for(k = 0; k < n; m_dk[k] = m_yvalues[k++]);

	// les pentes
	double epsilon = 0.5;
	double mid = 0.5;
	double _inf = mid -epsilon;
	double _sup = mid +epsilon;
	int nbiter = 100;

	m_ck[0] = 0.;
	for(k = 1; k < n; k++)
	{
		if (m_alpha.size()>0){
			m_cal_size=k;

			//if (m_calibrate){
			//	m_alpha[k] = RootFinder1D(ff1::mem_call(&CCubicspline::EvaluateDerivSecond,(*this))).Dichotomy(_inf,_sup,100,1.E-5,1.E-2,false);}
			//	m_alpha[k] = RootFinder1D(ff1::mem_call(&CCubicspline::EvaluateDerivSecond,(*this))).SimpleNewtonRaphson(_inf,_sup,nbiter,100,1.E-5,1.E-2,false);}
			//	if (m_alpha[k]==CREDIT_DEFAULT_VALUE)
			//	{m_alpha[k]=mid;}
			m_ck[k] = (m_alpha[k]*(m_yvalues[k] - m_yvalues[k-1]) / delta[k-1] + (1.-m_alpha[k])*(m_yvalues[k+1] - m_yvalues[k]) / delta[k]);
		}
		else
			m_ck[k] = 0.5*((m_yvalues[k] - m_yvalues[k-1]) / delta[k-1] + (m_yvalues[k+1] - m_yvalues[k]) / delta[k]);
	}
	m_ck[n] = 0.;

	// les coefficients
	double * alpha = new double [n];
	double * beta = new double [n];

	for(k = 0; k < n; k++)
	{
		alpha[k]	= m_yvalues[k+1] - m_dk[k] - m_ck[k] * delta[k];
		beta[k]		= m_ck[k+1] - m_ck[k];
	}

	for(k = 0; k < n; k++)
	{
		m_ak[k]	= (- 2. * alpha[k] + delta[k] * beta[k]) / (delta[k] * delta[k] * delta[k]);
		m_bk[k]	= (3. * alpha[k] / delta[k] - beta[k]) / delta[k];
	}

	delete [] delta;
	delete [] alpha;
	delete [] beta;
}

void CCubicspline::buildSplineWithC2Condition()
{
	delete [] m_mk;
	delete [] m_ak;
	delete [] m_bk;
	delete [] m_ck;
	delete [] m_dk;

	m_mk = 0;
	
	int k, n = m_size - 1;

	m_ak	= new double [n];
	m_bk	= new double [n];
	m_ck	= new double [n+1];
	m_dk	= new double [n];

	double * delta = new double [n];

	// les intervalles
	for(k = 0; k < n; k++) delta[k] = m_xvalues[k+1] - m_xvalues[k];

	// la valeur des polynomes aux points (xk)
	for(k = 0; k < n; m_dk[k] = m_yvalues[k++]);

	// les pentes
	m_ck[0] = 0.;
	for(k = 1; k < n; k++)
	{
		m_ck[k] = -2.*m_ck[k-1] + 0.5* delta[k-1]*delta[k-1] + 3.*(m_yvalues[k] - m_yvalues[k-1]) / delta[k];
	}
	m_ck[n] = 0.;

	// les coefficients
	for(k = 0; k < n; k++)
	{
		m_bk[k]	= -0.5 * delta[k];
		m_ak[k]	= (m_ck[k+1]-m_ck[k]-2.*m_bk[k]*delta[k])/(3.*delta[k]*delta[k]);
	}

	delete [] delta;
}


void CCubicspline::buildSplineWithC1Condition()
{
	delete [] m_mk;
	delete [] m_ak;
	delete [] m_bk;
	delete [] m_ck;
	delete [] m_dk;

	m_mk = 0;
	
	int k, n = m_size - 1;

	m_ak	= new double [n];
	m_bk	= new double [n];
	m_ck	= new double [n+1];
	m_dk	= new double [n];

	double * delta = new double [n];

	// les intervalles
	for(k = 0; k < n; k++) delta[k] = m_xvalues[k+1] - m_xvalues[k];

	// la valeur des polynomes aux points (xk)
	for(k = 0; k < n; m_dk[k] = m_yvalues[k++]);

	// les pentes
	m_ck[0] = 0.;
	for(k = 1; k < n; k++)
	{
		m_ck[k] = 0.5 * ((m_yvalues[k] - m_yvalues[k-1]) / delta[k-1] + (m_yvalues[k+1] - m_yvalues[k]) / delta[k]);
	}
	m_ck[n] = 0.;

	// les coefficients
	double * alpha = new double [n];
	double * beta = new double [n];

	for(k = 0; k < n; k++)
	{
		alpha[k]	= m_yvalues[k+1] - m_dk[k] - m_ck[k] * delta[k];
		beta[k]		= m_ck[k+1] - m_ck[k];
	}

	for(k = 0; k < n; k++)
	{
		m_ak[k]	= (- 2. * alpha[k] + delta[k] * beta[k]) / (delta[k] * delta[k] * delta[k]);
		m_bk[k]	= (3. * alpha[k] / delta[k] - beta[k]) / delta[k];
	}

	delete [] delta;
	delete [] alpha;
	delete [] beta;
}


static void __stdcall objfun_EvaluateDerivSecond(Integer n, double x[], double *objf,
                   double g[], Nag_Comm *comm)
{
	int i=0;
	vector<double> X;
	for (i=0;i<n;i++)
	{
	/*	if (x[i]<=1.e-20)
		{	*objf = 10.E20;
			return;	} */
		X.push_back(x[i]);
	}

    _spline_context* p = (_spline_context*)comm->p ;

	CCubicspline c(p->m_xvalues,p->m_yvalues,p->m_size,p->m_emmode,p->m_Condition,p->m_leftSlope,
				p->m_rightSlope,X,false);
	
	double result = 0.,delta=0.,deltb=0.;
	for (i=1;i<n;i++)
		{
/*		delta=c.m_xvalues[i] - c.m_xvalues[i-1];
		deltb=c.m_bk[i] - c.m_bk[i-1];
		result+= fabs(c.m_ak[i]+(1./(3.*delta))*deltb+c.m_rightSlope);
*/
		result+= fabs(c.m_bk[i]);
		;}

	*objf = result;
}

void CCubicspline::buildSplineOptimizedMulti()
{
	using namespace OptimTools;

	double _inf = 0.;
	double _sup = 1.;

	vector<double> X;X.resize(m_size);
	vector<double> bound_inf_X;bound_inf_X.resize(m_size);
	vector<double> bound_sup_X;bound_sup_X.resize(m_size);

	_spline_context ctxt;
	ctxt.m_xvalues=m_xvalues;
	ctxt.m_yvalues=m_yvalues;
	ctxt.m_size=m_size;
	ctxt.m_emmode=m_emmode;
	ctxt.m_Condition=_emC1_alpha;
	ctxt.m_leftSlope=m_leftSlope;
	ctxt.m_rightSlope=m_rightSlope;


	for (int i=0;i<m_size;i++)
	{
		bound_inf_X[i]=_inf;
		bound_sup_X[i]=_sup;
		X[i] = (_inf+_sup)/2.;
	}

	double result = OptimTools::NagMultiOptimisator((void*) &ctxt,
												objfun_EvaluateDerivSecond,
												X,
												bound_inf_X,
												bound_sup_X,1.E-3,100);

	m_alpha = X;

	m_Condition=_emC1_alpha;
	buildSplineOptimized();
}

void CCubicspline::buildSplineOrder4()
{
	if (m_mk) delete m_mk;
	if (m_ak) delete m_ak;
	if (m_bk) delete m_bk;
	if (m_ck) delete m_ck;
	if (m_dk) delete m_dk;
	if (m_ek) delete m_ek;

	m_mk = 0;
	
	//int k, n = m_size - 1;
	int k, n = m_size;

	m_ak	= new double [n];
	m_bk	= new double [n];
	m_ck	= new double [n+1];
	m_dk	= new double [n];
	m_ek	= new double [n];

	double * delta = new double [n];

	// les intervalles
	for(k = 0; k < n; k++) 
	{delta[k] = m_xvalues[k+1] - m_xvalues[k];}

	// la valeur des polynomes aux points (xk)
	for(k = 0; k < n; m_dk[k] = m_yvalues[k++]);

	// les pentes
	if (m_alpha.size()>0)
	{m_ck[0] = m_alpha[0]*(m_yvalues[1] - m_yvalues[0])/delta[0];}
	//m_ck[0] = (m_yvalues[1] - m_yvalues[0])/delta[0];
	m_ck[n] = 0.;

	for(k = 1; k < n; k++)
	{
		if (m_alpha.size()>0)
		{m_ck[k] = (m_alpha[k]*(m_yvalues[k] - m_yvalues[k-1]) / delta[k-1] + (1.-m_alpha[k])*(m_yvalues[k+1] - m_yvalues[k]) / delta[k]);}
		else
		{m_ck[k] = 0.5*((m_yvalues[k] - m_yvalues[k-1]) / delta[k-1] + (m_yvalues[k+1] - m_yvalues[k]) / delta[k]);}

		if ( (m_ck[k]<0.) && (((m_yvalues[k] - m_yvalues[k-1]) / delta[k-1])>0.) )
		{m_ck[k] = m_yvalues[k] - m_yvalues[k-1] / delta[k-1];}
		//else
		//{m_ck[k] = 0.;}

		if (m_ck[k]>m_ck[k-1])
			{m_ck[k]=m_ck[k];}

	}
	
	for(k = 0; k < n; k++) 
	{
		m_bk[k] = -m_rightSlope -(6./(2.*delta[k]))*(m_ck[k+1]+m_ck[k]) +(6./(delta[k]*delta[k]))*(m_dk[k+1]-m_dk[k]);
	}

	// les coefficients
	double * alpha = new double [n];
	double * beta = new double [n];

	for(k = 0; k < n; k++)
	{
		alpha[k]	= (m_dk[k+1] - m_dk[k]) - m_ck[k] * delta[k] - m_bk[k] * delta[k] * delta[k];
		beta[k]		= m_ck[k+1] - m_ck[k] - 2.* m_bk[k] * delta[k];
	}

	for(k = 0; k < n; k++)
	{
		m_ak[k]	= (4. * alpha[k] - delta[k] * beta[k]) / (delta[k] * delta[k] * delta[k]);
		m_ek[k]	= (delta[k] * beta[k] - 3.* alpha[k]) / (delta[k] * delta[k] * delta[k] * delta[k]);
	}

	delete [] delta;
	delete [] alpha;
	delete [] beta;
}


void CCubicspline::buildSplineOrder4CondDerivs()
{
	int i=0;
	int size = 0;
	double epsilon = m_leftSlope;

	buildSplineOrder4();

	vector<double> X;
	vector<double> Y;
	vector<double> AK;
	vector<double> BK;
	vector<double> CK;
	vector<double> DK;
	vector<double> EK;

	vector<double> X0;
	vector<double> Y0;
	vector<double> AK0;
	vector<double> BK0;
	vector<double> CK0;
	vector<double> DK0;
	vector<double> EK0;

	vector<double> X1;

	for (i=0;i<m_size;i++) 
	{
		X0.push_back(m_xvalues[i]);
		Y0.push_back(m_yvalues[i]);
		AK0.push_back(m_ak[i]);
		BK0.push_back(m_bk[i]);
		CK0.push_back(m_ck[i]);
		DK0.push_back(m_dk[i]);
		EK0.push_back(m_ek[i]);
	}

	size = (int)((m_xvalues[m_size-1]-m_xvalues[0])/epsilon);

	for (i=0;i<size;i++) 
	{	X1.push_back(m_xvalues[0]+i*epsilon);}

	MergeVector(X1,X0,X);	

	for (i=0;i<X.size();i++) 
	{	
		Y.push_back((*this)(X[i]));
		double interpol = LinearVectorInterpol(X0,CK0,X[i]);
		CK.push_back(interpol);
		interpol = FlatVectorInterpol(X0,AK0,X[i]);
		AK.push_back(interpol);
		interpol = FlatVectorInterpol(X0,BK0,X[i]);
		BK.push_back(interpol);
		interpol = FlatVectorInterpol(X0,DK0,X[i]);
		DK.push_back(interpol);
		interpol = FlatVectorInterpol(X0,EK0,X[i]);
		EK.push_back(interpol);
	}

	int l=1;
	for (i=1;i<X.size();i++) 
	{
		double delta = (X[i]-X[i-1]);
		Y[i] = Y[i-1] + delta*delta*delta/4.*AK[i-1] + delta/4.*(CK[i]+3*CK[i-1]) + delta*delta/2.*BK[i-1];
		
/*		if (CHECK_EQUAL(X[i],X0[l]))
		{
			Y[i]=Y0[l];
			l++;
		}
*/
	}

	if (m_xvalues) delete m_xvalues;
	if (m_yvalues) delete m_yvalues;
	if (m_mk) delete m_mk;
	if (m_ak) delete m_ak;
	if (m_bk) delete m_bk;
	if (m_ck) delete m_ck;
	if (m_dk) delete m_dk;
	if (m_ek) delete m_ek;

	m_xvalues=NULL;
	m_xvalues=NULL;
	m_mk=NULL;
	m_ak=NULL;
	m_bk=NULL;
	m_ck=NULL;
	m_dk=NULL;
	m_ek=NULL;

	m_size = X.size();

	m_xvalues	= new double [m_size];
	m_yvalues	= new double [m_size];

	for(int k = 0; k < m_size;k++)
	{
		m_xvalues[k] = X[k];
		m_yvalues[k] = Y[k];
	}

	m_mk = 0;
	
	//int k, n = m_size - 1;
	int n = m_size;

	m_ak	= new double [n];
	m_bk	= new double [n];
	m_ck	= new double [n+1];
	m_dk	= new double [n];
	m_ek	= new double [n];

	double * delta = new double [n];

	// les intervalles
	for(k = 0; k < n; k++) 
	{delta[k] = m_xvalues[k+1] - m_xvalues[k];}

	// la valeur des polynomes aux points (xk)
	for(k = 0; k < n; m_dk[k] = m_yvalues[k++]);

	// les pentes
	m_ck[0] = 0.;m_ck[n] = 0.;
	for(k = 1; k < n; k++)
	{
		m_ck[k] = CK[k];

		//if (k>0)
		if (m_ck[k]>m_ck[k-1])
			{m_ck[k]=m_ck[k];}

	}

	for(k = 0; k < n; k++) 
	{
		m_bk[k] = -m_rightSlope -(6./(2.*delta[k]))*(m_ck[k+1]+m_ck[k]) +(6./(delta[k]*delta[k]))*(m_dk[k+1]-m_dk[k]);
	}

	// les coefficients
	double * alpha = new double [n];
	double * beta = new double [n];

	for(k = 0; k < n; k++)
	{
		alpha[k]	= (m_dk[k+1] - m_dk[k]) - m_ck[k] * delta[k] - m_bk[k] * delta[k] * delta[k];
		beta[k]		= m_ck[k+1] - m_ck[k] - 2.* m_bk[k] * delta[k];
	}

	for(k = 0; k < n; k++)
	{
		m_ak[k]	= (4. * alpha[k] - delta[k] * beta[k]) / (delta[k] * delta[k] * delta[k]);
		m_ek[k]	= (delta[k] * beta[k] - 3.* alpha[k]) / (delta[k] * delta[k] * delta[k] * delta[k]);
	}

	delete [] delta;
	delete [] alpha;
	delete [] beta;
}


static void __stdcall objfun_EvaluateC4cond(Integer n, double x[], double *objf,
                   double g[], Nag_Comm *comm)
{
	int i=0;
	vector<double> X;
	for (i=0;i<n;i++)
	{
		X.push_back(x[i]);
	}

    _spline_context* p = (_spline_context*)comm->p ;

	CCubicspline c(p->m_xvalues,p->m_yvalues,p->m_size,p->m_emmode,p->m_Condition,p->m_leftSlope,
				p->m_rightSlope,X,false);

	vector<double> X0;
	vector<double> Y0;

	//norme 2
	double result = 0.,delta=0.,deltb=0.;
	for (i=0;i<n;i++)
	{
		delta= (c(p->m_xvalues[i])-p->m_yvalues[i]);
		result+= delta*delta;
	}

	*objf = result;
}

void CCubicspline::buildSplineOptimizedMultiC4()
{
	using namespace OptimTools;

	double _inf = 0.;
	double _sup = 1.;

	vector<double> X;X.resize(m_size);
	vector<double> bound_inf_X;bound_inf_X.resize(m_size);
	vector<double> bound_sup_X;bound_sup_X.resize(m_size);

	_spline_context ctxt;
	ctxt.m_xvalues=m_xvalues;
	ctxt.m_yvalues=m_yvalues;
	ctxt.m_size=m_size;
	ctxt.m_emmode=m_emmode;
	ctxt.m_Condition=_emOrder4Convex;
	ctxt.m_leftSlope=m_leftSlope;
	ctxt.m_rightSlope=m_rightSlope;


	for (int i=0;i<m_size;i++)
	{
		bound_inf_X[i]=_inf;
		bound_sup_X[i]=_sup;
		X[i] = (_inf+_sup)/2.;
	}

	double result = OptimTools::NagMultiOptimisator((void*) &ctxt,
												objfun_EvaluateC4cond,
												X,
												bound_inf_X,
												bound_sup_X,1.E-6,200);

	m_alpha = X;

	m_Condition=_emOrder4Convex;
	buildSplineOrder4CondDerivs();
}

void CCubicspline::buildSplineShumaker()
{
	if (m_mk) delete m_mk;
	if (m_ak) delete m_ak;
	if (m_bk) delete m_bk;
	if (m_ck) delete m_ck;
	if (m_dk) delete m_dk;
	if (m_ek) delete m_ek;
	if (m_fk) delete m_fk;

	int k=0, j=0, n = m_size;

	m_ak	= new double [n];
	m_bk	= new double [n];
	m_ck	= new double [n];
	m_dk	= new double [n];
	m_ek	= new double [n];
	m_fk	= new double [n];
	m_mk	= new double [n];

	double * delta = new double [n];
	double * deltax = new double [n];
	double * deltay = new double [n];
	double * d = new double [n];
	double * L = new double [n];
	double * LT = new double [n];
	double * dbar = new double [n];

	double a=0, b=0, sg=0, sd=0, alpha=0, beta=0, w=0;

	//Init pentes et cordes
	for(k = 1; k < n; k++) 
	{
		deltax[k] = m_xvalues[k] - m_xvalues[k-1];
		deltay[k] = m_yvalues[k] - m_yvalues[k-1];
		LT[k] = sqrt(deltax[k]*deltax[k] + deltay[k]*deltay[k]);
		L[k] = LT[k];
		delta[k] = deltay[k] / deltax[k];
	}

	for(k = 1; k < n; k++) 
	{
		j=k;
		sd = sg = 0.;
		while ((delta[j] == delta[j+1]) && (j<n) )
		{
			sd = sd + LT[j+1];
			j++;
		}
		j=k;
		while ((delta[j] == delta[j-1]) && (j>1) )
		{
			sd = sd + LT[j-1];
			j--;
		}
		LT[k] = L[k] + sg + sd;
	}

	for(k = 0; k < n; k++) 
		L[k] = LT[k];

	//Estimation pente pondérées
	for(k = 1; k < n-1; k++)
		d[k] = (L[k]*delta[k]+L[k+1]*delta[k+1]) / (L[k]+L[k+1]);
	d[0]= (3*delta[1]-d[1])/2;
	d[n-1]=0.;

	//Détermination points et coef
	for(k = 1; k < n; k++) 
	{
		if (m_dk[k]+m_dk[k-1]==2*delta[k])
		{
			m_ak[k] = m_yvalues[k-1];
			m_bk[k] = d[k-1];
			m_ck[k] = (d[k]-d[k-1])/(2*deltax[k]);
			//Pas de pts supp :  on init au pt de départ
			m_mk[k] = m_xvalues[k-1];
			m_dk[k] = m_yvalues[k-1];
			m_ek[k] = d[k-1];
			m_fk[k] = (d[k]-d[k-1])/(2*deltax[k]);
		}
		else
		{
			//Points
			a = d[k] - delta[k];
			b = d[k-1] - delta[k];
			if (a*b>=0)
				m_mk[k] = (m_xvalues[k]+m_xvalues[k-1])/2;
			else
			{
				if (fabs(a)>fabs(b))	
					m_mk[k] = m_xvalues[k] - (-b * deltax[k])/(d[k]-d[k-1]);
				else
					m_mk[k] = m_xvalues[k-1] + (a * deltax[k])/(d[k]-d[k-1]);
			}
			//Coef
			alpha = m_mk[k]-m_xvalues[k-1];
			beta = m_xvalues[k] - m_mk[k];
			w = alpha / deltax[k];

			dbar[k] = 2*delta[k] - (w*d[k-1]+(1-w)*d[k]);

			//1st part
			m_ak[k] = m_yvalues[k-1];
			m_bk[k] = d[k-1];
			m_ck[k] = (dbar[k]-d[k-1])/(2*alpha);
			
			//2nd part
			m_dk[k] = m_yvalues[k-1] + d[k-1]*alpha + (dbar[k]-d[k-1])*alpha/2;
			m_ek[k] = dbar[k];
			m_fk[k] = (d[k]-dbar[k])/(2*beta);
		}
	}
	
	//Libération mémoire
	delete [] delta;
	delete [] deltax;
	delete [] deltay;
	delete [] d;
	delete [] L;
	delete [] LT;
	delete [] dbar;
}