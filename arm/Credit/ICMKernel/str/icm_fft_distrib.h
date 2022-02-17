#error no longer part of the project
// AbstractCopulaCalculator.h: interface for the AbstractCopulaCalculator class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FFTDISTRIB_H__4CF84D2D_3364_11D7_8820_0002A5630CEB__INCLUDED_)
#define AFX_FFTDISTRIB_H__4CF84D2D_3364_11D7_8820_0002A5630CEB__INCLUDED_

#pragma warning (disable : 4786 )

#include <string>
#include <vector>
#include <map>

#include "icm_MathSrv.h"
#include "ICMKernel\util\icm_fft.h"
#include "ICMKernel\util\icm_ComplexHermiteIntegration.h"
#include "ARMKernel\glob\linalg.h"


class ICM_FFTDISTRIB  
{
//-----  Constructors/destructors
public:
	ICM_FFTDISTRIB(int size,double* pdef,double beta);
	ICM_FFTDISTRIB(int size,double* pdef,std::vector<double>& beta_vector);
	ICM_FFTDISTRIB(int size,double* pdef,ARM_Vector* beta_vector);
	ICM_FFTDISTRIB(int size,std::vector<double> pdef, std::vector<double>& beta_vector);

	virtual ~ICM_FFTDISTRIB();

//----- Assessors
public:
	int							getSize()		const {return m_size;}
	const std::vector<double>&	getDensity()	const;

	
	void computebarrier(); //clone
	void setpdefperturb(double* pdef_perturb);
	void setpdef(double* pdef);



//----- Utilities

	virtual void init();
	bool isComputed(const long& index) const;
	bool isComputed() const {return m_flag;}

	virtual double* getValue(double t)	const;

	virtual void	computeDensity();
	// ----- For Computing Delta
	virtual double* getValuePerturbed(double t)	const;

	virtual void	computeDensityPerturbed(const long& index);

	const std::vector<double>& getDensityPerturbed(const long& index) const;

	void	addDensityPerturbed(const long& index, const std::vector<double>&);

	double compute_expectedlosstranche(double tranche_up, 
									   double tranche_down, 
									   double lossunit);


	double compute_expectedlosstranche_perturb(double tranche_up, 
											   double tranche_down, 
											   double lossunit,
											   long index);

protected:
	void	clear();

//----- Tests
 
//-----  Attributes
protected:
	int		m_size;
	int		m_nbPoints;
	bool	m_flag;

	double*	m_cosx; // =cos(2*pi*x) for discretisation = (k/2^n) for computing fct caracteristic
	double*	m_sinx; // =sin(2*pi*x) for discretisation = (k/2^n) for computing fct caracteristic
	double* m_datas; // for fft
	double	m_beta; // beta
	std::vector<double> m_beta_vector; // for assessor, m_dens[k] = P(N(T)=k)

	double*	m_b; // barrier
	double*	m_pdef_perturb; // barrier
	double*	m_pdef; // barrier

	mutable int	m_ind_x; // discretisation for computing each Phi(x)
	mutable int	m_ind_t; // discretisation for computing Phi(x) by an integral of 20 points

	mutable double* m_complexe;

	mutable double* m_g; // values of g = product
	mutable double* m_f; // values of f = composante of g

	mutable int m_index;
	mutable double m_b_perturb;

	std::vector<double> m_dens; // for assessor, m_dens[k] = P(N(T)=k)
	std::map<long, std::vector<double> > m_dens_perturb;

	// For normalisation
	double			m_normNbPoints;
	static double	m_normPi;
	static double	m_2Pi;
	static	double	m_sqrt2;

};
#endif 