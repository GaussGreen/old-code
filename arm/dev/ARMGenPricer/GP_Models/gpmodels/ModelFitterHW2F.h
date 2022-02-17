/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 * 
 * \file ModelFitterHW2F.h
 *  \brief Calib Method dedicated to HW2F
 * 
 *	\author  J. Messines
 *	\version 1.0
 *	\date June 2006
 */

#ifndef _INGPMODELS_ModelFitterHW2F_H
#define _INGPMODELS_ModelFitterHW2F_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpcalib/modelfitter.h"
#include "gpnumlib/typedef.h"
#include "gpinfra/pricingstates.h"
#include "gpcalib/modelfitterdes.h"
#include "gpnumlib/numfunction.h"
#include "gpnumlib/levmarq.h"


CC_BEGIN_NAMESPACE( ARM )


class ARM_ModelFitterHW2F : public ARM_ModelFitter
{
private:
	size_t	itsNbEquations;
	bool	itsWithOptim;
public:
	
	/// constructor / destructor ) operator =
	ARM_ModelFitterHW2F(ARM_PricingModel* model, 
        const ARM_StdPortfolioPtr portfolio, 
        const ARM_ModelParamVector& calibParam,
		size_t nbEquations,
        ARM_ModelFitterPtr& linkedModelFitter	= ARM_ModelFitterPtr(NULL ));

    ARM_ModelFitterHW2F(const ARM_ModelFitterHW2F& rhs);
	virtual ~ARM_ModelFitterHW2F();

	ARM_ModelFitterHW2F& operator = (const ARM_ModelFitterHW2F& rhs);

    /// the core of the calibration
	virtual void Calibrate();
	virtual void StoreDatas();
	virtual void SolveEquations();
	virtual void CheckCoherence();
	
	/// standard ARM Object support
    virtual ARM_Object* Clone() const { return new ARM_ModelFitterHW2F(*this); }
	virtual string ModelFitterName() const;
    virtual string toString(const string& indent="", const string& nextIndent="") const;
	inline size_t GetNbEquations() const {return itsNbEquations;};
	double xi(double mrs,double t1,double t2) const {return (mrs<K_NEW_DOUBLE_TOL&&-mrs<K_NEW_DOUBLE_TOL)?t2-t1:(exp(mrs*t2)-exp(mrs*t1))/mrs;};
};

class ccalib  
{
protected:
	double m_vol;
	double m_ratio;
	double m_correl;

	double m_minVol;
	double m_minRatio;
	double m_minCorrel;
	double m_maxCorrel;

	ARM_GP_Vector m_a;
	ARM_GP_Vector m_b;
	ARM_GP_Vector m_c;
	ARM_GP_Vector m_d;
	ARM_GP_Vector m_sv;

	ARM_GP_Vector m_weight;
	bool m_switchToOptim;

public:
	ccalib(const ARM_GP_Vector& a, const ARM_GP_Vector& b, const ARM_GP_Vector& c, const ARM_GP_Vector& d, ARM_GP_Vector& sv)
		:m_a(a),m_b(b),m_c(c),m_d(d),m_sv(sv),m_minVol(1e-4),m_minRatio(0.1),m_minCorrel(-0.9999),m_maxCorrel(0.999),m_switchToOptim(false)
	{
		m_weight.resize(a.size());
		for (size_t k=0;k<a.size();k++)
			m_weight[k]=1.;
	};
	virtual ~ccalib(){};

	double getVol(){return m_vol;};
	double getRatio(){return m_ratio;};
	double getCorrel(){return m_correl;};

	void setVol(double vol){m_vol=CC_Max<double>(vol,m_minVol);};
	void setRatio(double ratio){m_ratio=CC_Max<double>(ratio,m_minRatio);};
	void setCorrel(double correl){m_correl=CC_Min<double>(CC_Max<double>(correl,m_minCorrel),m_maxCorrel);};
	
	virtual void solve(){};

	void keepInit();
	void switchTo1eqVol();
	void switchTo1eqRatio();
	void switchTo2eqCorrel();

	void setOptim(bool optim){m_switchToOptim=optim;};
	virtual void solveOptim();
	virtual void initOptim(ARM_GP_Vector& fguess){};
	virtual void setParam(double p[]){};

	double getWeight(int k) { return m_weight[k];};
	double getErr(int k)	{ return m_a[k]*m_vol*m_vol+m_b[k]*m_vol*m_vol*m_ratio*m_ratio+m_c[k]*m_vol*m_vol*m_ratio*m_correl-m_d[k];}

	class OptimCalib : public ARM_LEVMARQFunc
	{
	private:
		ccalib *			_this;
		
	public:
		OptimCalib(ccalib * calib)
		{
			_this			= calib;
		}
		~OptimCalib()
		{
		}

	public:
		void operator()(double p[], double hx[], int m, int n, void * adata = NULL) const
		{
			_this->setParam(p);
			for(int k = 0; k < n; k++)
				hx[k]	= _this->getWeight(k)*_this->getErr(k)*10000;
		}
	};
};

class ccalib3eq : public ccalib  
{
public:
	ccalib3eq(const ARM_GP_Vector& a, const ARM_GP_Vector& b, const ARM_GP_Vector& c, const ARM_GP_Vector& d, ARM_GP_Vector& sv)
		:ccalib(a,b,c,d,sv){};
	virtual ~ccalib3eq(){};
	virtual void solve();

	virtual void initOptim(ARM_GP_Vector& fguess);
	virtual void setParam(double p[]);

	double det()	{ return m_a[0]*(m_b[1]*m_c[2]-m_b[2]*m_c[1])-m_a[1]*(m_b[0]*m_c[2]-m_b[2]*m_c[0])+m_a[2]*(m_b[0]*m_c[1]-m_b[1]*m_c[0]);};
	double det1()	{ return m_d[0]*(m_b[1]*m_c[2]-m_b[2]*m_c[1])-m_d[1]*(m_b[0]*m_c[2]-m_b[2]*m_c[0])+m_d[2]*(m_b[0]*m_c[1]-m_b[1]*m_c[0]);};
	double det2()	{ return m_a[0]*(m_d[1]*m_c[2]-m_d[2]*m_c[1])-m_a[1]*(m_d[0]*m_c[2]-m_d[2]*m_c[0])+m_a[2]*(m_d[0]*m_c[1]-m_d[1]*m_c[0]);};
	double det3()	{ return m_a[0]*(m_b[1]*m_d[2]-m_b[2]*m_d[1])-m_a[1]*(m_b[0]*m_d[2]-m_b[2]*m_d[0])+m_a[2]*(m_b[0]*m_d[1]-m_b[1]*m_d[0]);};
};

class ccalib2eq : public ccalib  
{
protected:
	int m_i;
	int m_j;
public:
	ccalib2eq(const ARM_GP_Vector& a, const ARM_GP_Vector& b, const ARM_GP_Vector& c, const ARM_GP_Vector& d, ARM_GP_Vector& sv, int i, int j)
		:ccalib(a,b,c,d,sv),m_i(i),m_j(j)
	{
		for (size_t k=0;k<m_a.size();k++)
		{
			if ((m_i!=k) && (m_j!=k))
				m_weight[k]=0.;
		}
	};
	virtual ~ccalib2eq(){};
};

class ccalib2eqVol : public ccalib2eq  
{
private:
	double m_sigma;
public:
	ccalib2eqVol(const ARM_GP_Vector& a, const ARM_GP_Vector& b, const ARM_GP_Vector& c, const ARM_GP_Vector& d, ARM_GP_Vector& sv, int i, int j,double sigma)
		:ccalib2eq(a,b,c,d,sv,i,j),m_sigma(sigma){};
	virtual ~ccalib2eqVol(){};
	void solve();

	virtual void initOptim(ARM_GP_Vector& fguess);
	virtual void setParam(double p[]);

	double det() {return m_sv[0]*m_sv[0]*m_sv[0]*m_sv[0]*(m_b[m_i]*m_c[m_j]-m_b[m_j]*m_c[m_i]);}
	double det1(){return m_b[m_i]*m_sv[0]*m_sv[0]*(m_d[m_j]-m_a[m_j]*m_sv[0]*m_sv[0])-m_b[m_j]*m_sv[0]*m_sv[0]*(m_d[m_i]-m_a[m_i]*m_sv[0]*m_sv[0]);}
	double det2(){return (m_d[m_i]-m_a[m_i]*m_sv[0]*m_sv[0])*m_c[m_j]*m_sv[0]*m_sv[0]-(m_d[m_j]-m_a[m_j]*m_sv[0]*m_sv[0])*m_c[m_i]*m_sv[0]*m_sv[0];}
};

class ccalib2eqRho : public ccalib2eq  
{
private:
	double m_rho;
public:
	ccalib2eqRho(const ARM_GP_Vector& a, const ARM_GP_Vector& b, const ARM_GP_Vector& c, const ARM_GP_Vector& d, ARM_GP_Vector& sv, int i, int j,double rho)
		:ccalib2eq(a,b,c,d,sv,i,j),m_rho(rho){};
	virtual ~ccalib2eqRho(){};
	virtual void solve();

	virtual void initOptim(ARM_GP_Vector& fguess);
	virtual void setParam(double p[]);

};

class ccalib2eqRatio : public ccalib2eq  
{
private:
	double m_scale;
public:
	ccalib2eqRatio(const ARM_GP_Vector& a, const ARM_GP_Vector& b, const ARM_GP_Vector& c, const ARM_GP_Vector& d, ARM_GP_Vector& sv, int i, int j,double scale)
		:ccalib2eq(a,b,c,d,sv,i,j),m_scale(scale){};
		virtual ~ccalib2eqRatio(){};
	virtual void solve();

	virtual void initOptim(ARM_GP_Vector& fguess);
	virtual void setParam(double p[]);

	double det() {return (m_a[m_i]+m_b[m_i]*m_sv[1]*m_sv[1])*m_c[m_j]*m_sv[1]-(m_a[m_j]+m_b[m_j]*m_sv[1]*m_sv[1])*m_c[m_i]*m_sv[1];}
	double det1(){return m_d[m_i]*m_c[m_j]*m_sv[1]-m_d[m_j]*m_c[m_i]*m_sv[1];}
	double det2(){return (m_a[m_i]+m_b[m_i]*m_sv[1]*m_sv[1])*m_d[m_j]-(m_a[m_j]+m_b[m_j]*m_sv[1]*m_sv[1])*m_d[m_i];}
};

class ccalib1eq : public ccalib  
{
protected:
	int m_i;
public:
	ccalib1eq(const ARM_GP_Vector& a, const ARM_GP_Vector& b, const ARM_GP_Vector& c, const ARM_GP_Vector& d, ARM_GP_Vector& sv, int i)
		: ccalib(a,b,c,d,sv),m_i(i)
	{
		for (size_t k=0;k<m_a.size();k++)
		{
			if (m_i!=k)
				m_weight[k]=0.;
		}
	};
	virtual ~ccalib1eq(){};
};

class ccalib1eqRatio : public ccalib1eq  
{
public:
	ccalib1eqRatio(const ARM_GP_Vector& a, const ARM_GP_Vector& b, const ARM_GP_Vector& c, const ARM_GP_Vector& d, ARM_GP_Vector& sv, int i)
		: ccalib1eq(a,b,c,d,sv,i){};
	virtual ~ccalib1eqRatio(){};
	virtual void solve();

	virtual void initOptim(ARM_GP_Vector& fguess);
	virtual void setParam(double p[]);

};

class ccalib1eqRho : public ccalib1eq 
{
public:
	ccalib1eqRho(const ARM_GP_Vector& a, const ARM_GP_Vector& b, const ARM_GP_Vector& c, const ARM_GP_Vector& d, ARM_GP_Vector& sv, int i)
		: ccalib1eq(a,b,c,d,sv,i){};
	virtual ~ccalib1eqRho(){};
	virtual void solve();

	virtual void initOptim(ARM_GP_Vector& fguess);
	virtual void setParam(double p[]);

};

class ccalib1eqVol : public ccalib1eq  
{
public:
	ccalib1eqVol(const ARM_GP_Vector& a, const ARM_GP_Vector& b, const ARM_GP_Vector& c, const ARM_GP_Vector& d, ARM_GP_Vector& sv, int i)
		: ccalib1eq(a,b,c,d,sv,i){};
	virtual ~ccalib1eqVol(){};
	virtual void solve();

	virtual void initOptim(ARM_GP_Vector& fguess);
	virtual void setParam(double p[]);

};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
