#include "ARMKernel\glob\firsttoinc.h"
#include "ICMKernel\glob\icm_calibrator.h"
#include "ICMKernel\inst\icm_mez.h"
#include "ICMKernel\util\icm_matrix.h"
#include "ICMKernel\util\icm_RootFinderND.h"
#include <nag.h>
#include "ICMKernel\glob\icm_maths.h"
#include "ICMKernel\pricer\icm_pricer_adviser.h"

//----------------------------------------------------------------------------
// Global evaluation functions
//----------------------------------------------------------------------------
double ICM_Calibrator::EvaluateFunci(int K, double* x)
{
	int i=0,size=0,j=0,k=0;	
	double result = 0.,val=0.;
	size = itsPricers.size();

	for (i=0;i<itsParamsVector.size();i++)
	{
		ARM_Vector* param = itsCalParameters->GetColVect(itsParamsVector[i].c_str());
		for (j=0;j<itsParamsVectorTS[i];j++)
		{	param->Elt(j)=x[k];
			k++;}
	}

	(itsPricers[K])->ResetPricer();
	(itsPricers[K])->SetParameters(*itsCalParameters);
	(itsPricers[K])->RefreshParameters();

	val = (itsPricers[K])->ComputePrice(itsPricingType) ;

	if ((val<itsPriceVectorBid[K])||(val>itsPriceVectorAsk[K])) 
		val -= (itsPriceVectorBid[K]+itsPriceVectorAsk[K])/2.;
	else
		val =0.;

	if ((itsNorm) && (itsPriceVectorBid[K]+itsPriceVectorAsk[K]))
		{val /= (itsPriceVectorBid[K]+itsPriceVectorAsk[K])/2.;}

	result = ABS(val); 
	return result;
}

double ICM_Calibrator::EvaluateSumFuncs(double* x)
{
	int i=0,size=0,j=0,k=0;	
	double result = 0.,val=0.;
	size = itsPricers.size();

	for (i=0;i<itsParamsVector.size();i++)
	{
		ARM_Vector* param = itsCalParameters->GetColVect(itsParamsVector[i].c_str());
		for (j=0;j<itsParamsVectorTS[i];j++)
		{	param->Elt(j)=x[k];
			k++;}
	}

	for (i=0;i<size;i++)
	{
		(itsPricers[i])->ResetPricer();
		(itsPricers[i])->SetParameters(*itsCalParameters);
		(itsPricers[i])->RefreshParameters();
		val = (itsPricers[i])->ComputePrice(itsPricingType);

		if ((val<itsPriceVectorBid[i])||(val>itsPriceVectorAsk[i])) 
			val -= (itsPriceVectorBid[i]+itsPriceVectorAsk[i])/2.;
		else
			val =0.;

		if ((itsNorm) && (itsPriceVectorBid[i]+itsPriceVectorAsk[i]))
		{val /= (itsPriceVectorBid[i]+itsPriceVectorAsk[i])/2.;}

		result += ABS(val); 
	}

	//result = sqrt(result);
	return result;
}

void ICM_Calibrator::EvaluateFuncs(int nbfuncs,double* x,double* funcsvalues)
{
	int i=0,size=0,j=0,k=0;	
	double result = 0.,val=0.;

	for (i=0;i<nbfuncs;i++)
	{funcsvalues[i]=EvaluateFunci(i,x);}
}

//----------------------------------------------------------------------------
// Ganso Call
//----------------------------------------------------------------------------
static void objfun_Optimize_Ganso(int* n, double* x, double* objf)
{
	_context_cal* p = (_context_cal*)OptimTools::GansoLibContext ;
	*objf = p->m_Calibrator->EvaluateSumFuncs(x);
}

//----------------------------------------------------------------------------
// PSwarm Call
//----------------------------------------------------------------------------
extern "C" struct Stats stats;

double objfun_Optimize_PSWarm(int n, double *y, double *lb, double *ub)
{
	double fx;
	vector<double> x;x.resize(n);
	_context_cal* p = (_context_cal*)OptimTools::GansoLibContext ;
	double* scal = OptimTools::scale_ = p->m_scale;

	stats.objfunctions++;

	for(int i=0;i<n;i++)
	{
	  if(y[i]<lb[i] || y[i]>ub[i])
	  {return 1e20;}
	  x[i]=y[i]*scal[i];
	}

// FIXMEFRED: mig.vc8 (28/05/2007 10:22:16):cast
	fx = p->m_Calibrator->EvaluateSumFuncs(&(*x.begin()));

	return fx;
}

//----------------------------------------------------------------------------
// Nag Call
//----------------------------------------------------------------------------
static void __stdcall objfun_Optimize(Integer n, double x[], double *objf,
                   double g[], Nag_Comm *comm)
{
	_context_cal* p = (_context_cal*)comm->p ;
	*objf = p->m_Calibrator->EvaluateSumFuncs(x);
}

static void __stdcall objfun_Optimize_BFGS(long M, long N, double X[],
                                    double F[], double fjac[],
                                    long tdfjac, Nag_Comm* comm)
{
	_context_cal* p = (_context_cal*)comm->p ;
    p->m_Calibrator->funcNagCall(2,M, N, X, F, fjac);
}

void ICM_Calibrator::NAG_ComputeGradient(int m, int n,double* params, double* GRAD_FUNC)
{
    int i, j;

	#define ICM_DEF_TOL         1.0e-4 
	#define ICM_GRAD_EPS        1.0e-5

	vector<double> X_PLUS_H;X_PLUS_H.resize(n);
	vector<double> FX_PLUS_H;FX_PLUS_H.resize(m);
	vector<double> X_MOINS_H;X_MOINS_H.resize(n);
	vector<double> FX_MOINS_H;FX_MOINS_H.resize(m);
	vector<double> H;H.resize(n);

    for (j = 0; j < n; j++)
    {
// FIXMEFRED: mig.vc8 (28/05/2007 10:22:38):cast
		MEMCPY(&(*X_PLUS_H.begin()), params, n*sizeof(double));
       MEMCPY(&(*X_MOINS_H.begin()), params, n*sizeof(double));

       H[j] = ICM_GRAD_EPS*MAX(fabs(params[j]), 1.0);

       // F(X+H)
       X_PLUS_H[j] = params[j]+H[j];

       EvaluateFuncs(m,&(*X_PLUS_H.begin()),&(*FX_PLUS_H.begin()));
 
       // F(X-H)
       X_MOINS_H[j] = params[j]-H[j];

       EvaluateFuncs(m,&(*X_MOINS_H.begin()),&(*FX_MOINS_H.begin()));

       for (i = 0; i < m; i++)
       { double funcDeriv = GRAD_FUNC[i*n+j] = (FX_PLUS_H[i]-FX_MOINS_H[i])/(2.0*H[j]); }
    }
}


void ICM_Calibrator::funcNagCall(int MODE,long M,long N,double* X,double* F,double* FJAC)
{
    switch(MODE)
    {
        case 0 :
        {EvaluateFuncs(M,X, F);};
        break;
        case 1 :
        {if (itsGradientComp)
		{NAG_ComputeGradient(M, N, X, FJAC);}};
        break;
        case 2 :
        {EvaluateFuncs(M,X, F);
         if (itsGradientComp)
		 {NAG_ComputeGradient(M, N, X, FJAC);}};
        break;
    };
}

//----------------------------------------------------------------------------
// Levenberg - Marquardt
//----------------------------------------------------------------------------
void LMFunction_(double *x, double *objf, int n, int m, void *adata)
{
	_context_cal* p = (_context_cal*)OptimTools::GansoLibContext ;
	p->m_Calibrator->EvaluateFuncs(m,x,objf);
}



ICM_Parameters* ICM_Calibrator::Optimize()
{
	int i=0,j=0,k=0;

	for (i=0;i<itsPricers.size();i++)
		{if (itsPricers[i]) delete itsPricers[i];
		itsPricers[i]=NULL;}

	itsPricers.resize(itsSecurityVector.size());

	ICM_Pricer_Advisor advisor;

	for (i=0;i<itsPricers.size();i++)
	{itsPricers[i]=	advisor.GeneratePricer(itsSecurityVector[i],itsModel,itsPricerType,0,itsParameters,itsModel->GetStartDate());}

	itsnbParams=0;

	for (i=0;i<itsParamsVector.size();i++)
	{itsnbParams += (int) itsParamsVectorTS[i];}

	if (itsCalParameters)
		delete itsCalParameters;
	itsCalParameters = (ICM_Parameters*) itsParameters->Clone();

	vector<double> X;X.resize(itsnbParams);
	vector<double> bound_inf_X;bound_inf_X.resize(itsnbParams);
	vector<double> bound_sup_X;bound_sup_X.resize(itsnbParams);

	for (i=0;i<itsParamsVector.size();i++)
	{
		ARM_Vector* param = itsParameters->GetColVect(itsParamsVector[i].c_str());
		ARM_Vector* param_inf = itsParameters_inf->GetColVect(itsParamsVector[i].c_str());
		ARM_Vector* param_sup = itsParameters_sup->GetColVect(itsParamsVector[i].c_str());

		for (j=0;j<itsParamsVectorTS[i];j++)
		{
			X[k] = param->Elt(j);
			bound_inf_X[k] = param_inf->Elt(j);
			bound_sup_X[k] = param_sup->Elt(j);
			k++;
		}
		
	}

	_context_cal ctxt;
	ctxt.m_Calibrator = this;
	ctxt.m_Parameters = itsCalParameters;

	int maxiter = 100;
	double tol = 1.E-3;
	double stepmax = 1.E-3;
	double norm = 0.;
	double alpha =0.5;
	double MAXSIM = 0.;
	double MAXTIME = 0.;
	double TRACE = 0.;

	if (itsOptimParameters)
	{
		ARM_Vector* Ptol = itsOptimParameters->GetColVect("TOL");
		if (Ptol) {tol=(double)Ptol->Elt(0);}
		ARM_Vector* Pmaxiter = itsOptimParameters->GetColVect("MAXITER");
		if (Pmaxiter) {maxiter=(int)Pmaxiter->Elt(0);}
		ARM_Vector* Pnorm = itsOptimParameters->GetColVect("NORM");
		if (Pnorm) {itsNorm=norm=(double)Pnorm->Elt(0);}
		ARM_Vector* Pstepmax = itsOptimParameters->GetColVect("STEPMAX");
		if (Pstepmax) {stepmax=(double)Pstepmax->Elt(0);}
		ARM_Vector* Plineserach = itsOptimParameters->GetColVect("LINESEARCH");
		if (Plineserach) {alpha=(double)Plineserach->Elt(0);}
		ARM_Vector* POptType = itsOptimParameters->GetColVect("OPTTYPE");
		if (POptType) {itsOptType=(qOPTIMIZE_TYPE)(int)POptType->Elt(0);}
		ARM_Vector* pMAXTIME = itsOptimParameters->GetColVect("MAXTIME");
		if (pMAXTIME) {MAXTIME= (double)pMAXTIME->Elt(0);}
		ARM_Vector* pMAXSIM = itsOptimParameters->GetColVect("MAXSIM");
		if (pMAXSIM) {MAXSIM= (double)pMAXSIM->Elt(0);}
		ARM_Vector* pTRACE = itsOptimParameters->GetColVect("TRACE");
		if (pTRACE) {TRACE= (double)pTRACE->Elt(0);}
	}

	OptimTools::GansoLibContext = &ctxt;
	((_context_cal*)OptimTools::GansoLibContext)->m_Calibrator = this;

	double error = 0.;
	void* func = NULL;

	switch (itsOptType)
	{
	case qOPT_LEVMAR :
	func = (void*)LMFunction_;
	break;
	case qOPT_GLOBAL_UNCON :
	func = (void*)objfun_Optimize;
	break;
	case qOPT_GANSO :
	func = (void*)objfun_Optimize_Ganso;
	break;
	case qOPT_NAG_NOLIN_BC_BFGS :
	func = (void*)objfun_Optimize_BFGS;
	break;
	case qOPT_PSWARM :
	func = (void*)objfun_Optimize_PSWarm;
	break;
	case qOPT_NAG_NOLIN_BC :
	default :
	func = (void*)objfun_Optimize;
	break;
	}

	if (MAXSIM)
	{
		error = OptimTools::AleaOptimSearch(itsOptType,(void*) &ctxt,func,X,bound_inf_X,bound_sup_X,itsPricers.size(),
							  MAXSIM,MAXTIME,tol,maxiter,stepmax,alpha,TRACE);
	}
	else
	{
		error = OptimTools::MultiOptimisator(itsOptType,(void*) &ctxt,func,X,bound_inf_X,bound_sup_X,itsPricers.size(),
							  tol,maxiter,stepmax,alpha);
	}

	ARM_Vector* Cal = new ARM_Vector(itsParameters->GetColVect(itsParamsVector[0].c_str())->size(),0.);
	Cal->Elt(0)= error;
	itsCalParameters->Push(Cal,"ERROR");

	k=0;
	for (i=0;i<itsParamsVector.size();i++)
	{
		ARM_Vector* param = itsCalParameters->GetColVect(itsParamsVector[i].c_str());

		for (j=0;j<itsParamsVectorTS[i];j++)
		{param->Elt(j)=X[k];k++;}
	}

	return itsCalParameters;
}

