/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  Partial Differential Equation Solving framework 
 *
 *	\file PDE_Solver.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date December 2005
 */
 
#ifndef _GP_CF_PDE_SOLVER_H
#define _GP_CF_PDE_SOLVER_H

#include "firsttoinc.h"

#include <vector>

#include "gpbase/port.h"
#include "gpbase/curve.h"
#include "gpbase/curvetypedef.h"
#include "gpbase/gplinalgconvert.h"


using std::vector;

CC_BEGIN_NAMESPACE(ARM)


//////////////////////////////////////////////////////////////////////
//
//
//
//
//  In order to use the schemas, the PDE should be deriving from one of the classes below:
//
//		either  PDE_JTX_Explicit_Schema , PDE_JTX_Theta_Schema or PDE_JTX_MidPoint_Schema 
//
//		 if it has (J)umps (T)ime dependent and (X) space dependent coefficients
//
//
//
/////////////////////////////////////////////////////////////////////


class PDE_JTX
{
	public:
		double f0;
		vector<double>* fgrid;
		double maturity;
		vector<double>* tgrid;
		virtual	double localdrift(int&,int&)=0;
		virtual	double localvol(int&,int&)=0;
		virtual	double localjump_probability(int&)=0;
		virtual	double localdriftjump(int&,int&)=0;
		virtual	double localvoljump(int&,int&)=0;
		virtual double localriskneutralrate(int&)=0;
		virtual double terminalvalue(int&)=0;
		virtual double jumpupintegral(int& k,int& i,double& eta)=0;
		virtual double jumpdownintegral(int& k,int& i,double& eta)=0;
		virtual double NormalizationService(double fn,int& k, int& i)=0;
		int payofftype;
};



class PDE_TX
{
	public:
		double f0;
		vector<double>* fgrid;
		double maturity;
		vector<double>* tgrid;
		virtual	double localdrift(int&,int&)=0;
		virtual	double localvol(int&,int&)=0;
		virtual double localriskneutralrate(int&)=0;
		virtual double terminalvalue(int&)=0;
		virtual double NormalizationService(double fn,int& k, int& i)=0;
		int payofftype;
};

class PDE_T
{
	public:
		double f0;
		vector<double>* fgrid;
		double maturity;
		vector<double>* tgrid;
		virtual	double localdrift(int&)=0;
		virtual	double localvol(int&)=0;
		virtual double localriskneutralrate(int&)=0;
		virtual double terminalvalue(int&)=0;
			virtual double NormalizationService(double fn,int& k)=0;
		int payofftype;
};


class PDE_X
{
	public:
		double f0;
		vector<double>* fgrid;
		double maturity;
		vector<double>* tgrid;
		virtual	double localdrift(int&)=0;
		virtual	double localvol(int&)=0;
		virtual double localriskneutralrate()=0;
		virtual double terminalvalue(int&)=0;
			virtual double NormalizationService(double fn, int& i)=0;
		int payofftype;
};

////////////////////////////////////////////////////////////////////////
///
///
///           Explicit Schema
///
///
////////////////////////////////////////////////////////////////////////

class PDE_JTX_Explicit_Schema :public PDE_JTX
{
public:
	void Set_Theta(double theta0){}
	void backwardcolumn(vector<double>& pn);
	
};

class PDE_TX_Explicit_Schema :public PDE_TX
{
public:
	void Set_Theta(double theta0){}
	void backwardcolumn(vector<double>& pn);
};


class PDE_T_Explicit_Schema :public PDE_T
{
public:
	void Set_Theta(double theta0){}
	void backwardcolumn(vector<double>& pn);
	
};

class PDE_X_Explicit_Schema :public PDE_X
{
public:
	void Set_Theta(double theta0){}
	void backwardcolumn(vector<double>& pn);
};


////////////////////////////////////////////////////////////////////////
///
///
///           Theta Schema
///
///
////////////////////////////////////////////////////////////////////////


class PDE_JTX_Theta_Schema :public PDE_JTX
{
public:
	double  theta;
	void Set_Theta(double theta0)
	{
		theta=theta0;
	}
	void backwardcolumn(vector<double>& pn);
};

class PDE_TX_Theta_Schema :public PDE_TX
{
public:
	double  theta;
	void Set_Theta(double theta0)
	{
		theta=theta0;
	}
	void backwardcolumn(vector<double>& pn);
};

class PDE_T_Theta_Schema :public PDE_T
{
public:
	double  theta;
	void Set_Theta(double theta0)
	{
		theta=theta0;
	}
	void backwardcolumn(vector<double>& pn);
};

class PDE_X_Theta_Schema :public PDE_X
{
public:
	double  theta;
	void Set_Theta(double theta0)
	{
		theta=theta0;
	}
	void backwardcolumn(vector<double>& pn);
};




////////////////////////////////////////////////////////////////////////
///
///
///           MidPoint Second Order  Schema
///
///
////////////////////////////////////////////////////////////////////////

class PDE_JTX_MidPoint_Schema :public PDE_JTX
{
public:
	void Set_Theta(double theta0){}
	void backwardcolumn(vector<double>& pn);
};

class PDE_TX_MidPoint_Schema :public PDE_TX
{
public:
	void Set_Theta(double theta0){}
	void backwardcolumn(vector<double>& pn);
};

class PDE_T_MidPoint_Schema :public PDE_T
{
public:
	void Set_Theta(double theta0){}
	void backwardcolumn(vector<double>& pn);
};

class PDE_X_MidPoint_Schema :public PDE_X
{
public:
	void Set_Theta(double theta0){}
	void backwardcolumn(vector<double>& pn);
};




////////////////////////////////////////////////////////////////////////
///
///
///           Generic Solving Function 
///
///
////////////////////////////////////////////////////////////////////////

/// return the all collumn
template <typename PDEType> void SolveColumn(PDEType *C,vector<double>& pn)
{
	return C->backwardcolumn(pn);
}

/// compute the price
template <typename PDEType> double Solve(PDEType *C)
{
	int L=C->fgrid->size();
	vector<double> pn(L);
	C->backwardcolumn(pn);
	std::vector<double>* fk_GP=CreateARMGPVectorFromVECTOR( *(C->fgrid));
	std::vector<double>* pn_GP=CreateARMGPVectorFromVECTOR( pn);
	ARM_Curve cc(*fk_GP,*pn_GP,new ARM_LinInterpCstExtrapolDble);
	double result=cc.Interpolate(C->f0);
	return result; /// autrefois c'etait pn[(L-1)/2];
}


/// compute the price and the delta
template <typename PDEType> void DeltaSolve(PDEType *C,double& result,double& delta,double relativeShift)
{
	int L=C->fgrid->size();
	vector<double> pn(L);
	C->backwardcolumn(pn);
	std::vector<double>* fk_GP=CreateARMGPVectorFromVECTOR( *(C->fgrid));
	std::vector<double>* pn_GP=CreateARMGPVectorFromVECTOR( pn);
	//ARM_Curve cc(*fk_GP,*pn_GP,new ARM_LinInterpCstExtrapolDble);
	//result=cc.Interpolate(C->f0);
	//delta = (cc.Interpolate(C->f0*(1.+relativeShift))-cc.Interpolate(C->f0*(1.-relativeShift)))/(2.*relativeShift);

}



CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

