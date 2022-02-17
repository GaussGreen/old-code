/********************************************************************************/
/*! \file icm_calibrator.h
 * 
 *  \brief Describes a Correlation Object
 *  \author Pouponneau Damien
 *	\version 1.0
 *	\date   February 2007 */
/*
 *********************************************************************************/

#ifndef _ICM_CALIBRATOR_H_
#define _ICM_CALIBRATOR_H_

/*********************************************************************************/
/*! \class  ICM_Correlation icm_calibrator.h "icm_calibrator.h"
 *  \author D Pouponneau
 *	\version 1.0
 *	\date   February 2007
 *	\brief 
/***********************************************************************************/

#include "ARMKernel\glob\firsttoinc.h"
#include "ARMKernel\glob\linalg.h"
#include "ARMKernel\inst\security.h"
#include "ICMKernel\glob\icm_enums.h"
#include "ICMKernel\util\icm_qmatrix.h"
#include "ICMKernel\pricer\icm_pricer.h"
#include <nag.h>

class ICM_Security;
class ICM_Parameters;

static void __stdcall objfun_Optimize(Integer n, double x[], double *objf,
                   double g[], Nag_Comm *comm);

static void objfun_Optimize_Ganso(int* n, double* x, double* objf);


class ICM_Calibrator : public ARM_Object
{
public:
vector<ICM_Security*>	itsSecurityVector;
vector<double>			itsPriceVectorBid;
vector<double>			itsPriceVectorAsk;
vector<string>			itsParamsVector;
vector<double>			itsParamsVectorTS;
ARM_Model*				itsModel;
ICM_Parameters*			itsOptimParameters;
ICM_Parameters*			itsParameters;
ICM_Parameters*			itsParameters_inf;
ICM_Parameters*			itsParameters_sup;
ARM_CLASS_NAME			itsPricerType;

ICM_Parameters*			itsCalParameters;

vector<ICM_Pricer*>		itsPricers;
qCMPMETH				itsPricingType;

int						itsnbParams;
double					itsNorm;
bool					itsGradientComp;
qOPTIMIZE_TYPE			itsOptType;

public:
	
	void Init(void)
	{
		itsModel=NULL;	
		itsParameters=NULL;
		itsParameters_inf=NULL;
		itsParameters_sup=NULL;
		itsnbParams = 0;
		itsCalParameters=NULL;
		itsOptimParameters=NULL;
		itsNorm=0.;
		itsGradientComp=true;
		itsOptType=qOPT_NAG_NOLIN_BC;
	}

	void Set(vector<ICM_Security*>&	SecurityVector,
			vector<double>&			PriceVectorBid,
			vector<double>&			PriceVectorAsk,
			vector<string>&			ParamsVector,
			vector<double>&			ParamsVectorTS,
			ARM_Model*				Model,
			ARM_CLASS_NAME			PricerType,
			ICM_Parameters*			Parameters,
			ICM_Parameters*			Parameters_inf,
			ICM_Parameters*			Parameters_sup,
			qCMPMETH				PricingType,
			ICM_Parameters*			OptimParameters)
	{
	itsSecurityVector = SecurityVector;
	itsPriceVectorBid= PriceVectorBid;
	itsPriceVectorAsk= PriceVectorAsk;
	itsParamsVector= ParamsVector;
	itsParamsVectorTS= ParamsVectorTS;
	itsModel=Model;	
	itsPricerType=PricerType;
	itsParameters=Parameters;
	itsParameters_inf=Parameters_inf;
	itsParameters_sup=Parameters_sup;
	itsPricingType = PricingType;
	itsOptimParameters=OptimParameters;
	}

	ICM_Calibrator() {Init();}		

	~ICM_Calibrator() 
	{
		for (int i=0;i<itsPricers.size();i++)
		{if (itsPricers[i]) delete itsPricers[i];
		itsPricers[i]=NULL;}
	}		

	ICM_Calibrator(vector<ICM_Security*>&	SecurityVector,
			vector<double>&			PriceVectorBid,
			vector<double>&			PriceVectorAsk,
			vector<string>&			ParamsVector,
			vector<double>&			ParamsVectorTS,
			ARM_Model*				Model,
			ARM_CLASS_NAME			PricerType,
			ICM_Parameters*			Parameters,
			ICM_Parameters*			Parameters_inf,
			ICM_Parameters*			Parameters_sup,
			qCMPMETH				PricingType,
			ICM_Parameters*			OptimParameters) 
	{
		Init();
		Set(SecurityVector,PriceVectorBid,PriceVectorAsk,ParamsVector,ParamsVectorTS,Model,PricerType,Parameters,Parameters_inf,Parameters_sup,PricingType,OptimParameters);
	}

public :
	
	ICM_Parameters* Optimize();

	void funcNagCall(int MODE,long M,long N,double* X,double* F,double* FJAC);
	double EvaluateFunci(int i, double* x);
	void EvaluateFuncs(int nbfuncs,double* x,double* funcsvalues);
	double EvaluateSumFuncs(double* x);
	void NAG_ComputeGradient(int m, int n,double* params, double* GRAD_FUNC);

};


class _context_cal
{
public:
	ICM_Parameters*			m_Parameters;
	ICM_Calibrator*			m_Calibrator;
	void*					m_evrythg;
	double*					m_scale;

	_context_cal(const _context_cal& in)
	{
		m_Parameters = in.m_Parameters;
		m_Calibrator = in.m_Calibrator;
		m_scale = in.m_scale;
		m_evrythg = in.m_evrythg;
	}

	void reset()
	{
		m_Parameters = NULL;
		m_Calibrator = NULL;
		m_scale = NULL;
		m_evrythg=NULL;
	}

	_context_cal()
	{reset();}
};


#endif
