// AbstractCopulaCalculator.cpp: implementation of the AbstractCopulaCalculator class.
//
//////////////////////////////////////////////////////////////////////
#include "ARMKernel\glob\firsttoinc.h"

#include "ICMKernel/glob/icm_maths.h" 
#include "icm_fft_distrib.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

double
ICM_FFTDISTRIB::m_normPi = 1./sqrt(PI);
//----------------------------------------------------------------------------------------------------------------------
double
ICM_FFTDISTRIB::m_2Pi = 2.*PI;
//----------------------------------------------------------------------------------------------------------------------
double
ICM_FFTDISTRIB::m_sqrt2 = sqrt(2.);
//----------------------------------------------------------------------------------------------------------------------
ICM_FFTDISTRIB::ICM_FFTDISTRIB(int size,double* pdef,double beta):
	m_size	(size),
	m_flag	(false),
	m_complexe (NULL), 
	m_cosx	(NULL), 
	m_sinx	(NULL), 
	m_datas	(NULL), 
	m_b		(NULL),
	m_g		(NULL),
	m_beta  (beta),
	m_f		(NULL),
	m_pdef_perturb   (NULL),
	m_pdef   (NULL)

{
	init();
	setpdef(pdef);
	computebarrier();

	m_beta_vector.resize(m_size);

	for(int k=0;k<m_size;k++)	{m_beta_vector[k]=beta;}

}

ICM_FFTDISTRIB::ICM_FFTDISTRIB(int size,double* pdef,std::vector<double>& beta_vector):
	m_size	(size),
	m_flag	(false),
	m_complexe (NULL), 
	m_cosx	(NULL), 
	m_sinx	(NULL), 
	m_datas	(NULL), 
	m_b		(NULL),
	m_g		(NULL),
	m_f		(NULL),
	m_pdef_perturb   (NULL),
	m_pdef   (NULL)

{
	init();
	setpdef(pdef);

	m_beta_vector.resize(m_size);
	for(int k=0;k<m_size;k++) 	{m_beta_vector[k]=beta_vector[k];}

	computebarrier();
}


ICM_FFTDISTRIB::ICM_FFTDISTRIB(int size,double* pdef,ARM_Vector* beta_vector):
	m_size	(size),
	m_flag	(false),
	m_complexe (NULL), 
	m_cosx	(NULL), 
	m_sinx	(NULL), 
	m_datas	(NULL), 
	m_b		(NULL),
	m_g		(NULL),
	m_f		(NULL),
	m_pdef_perturb   (NULL),
	m_pdef   (NULL)

{
	init();
	setpdef(pdef);

	m_beta_vector.resize(m_size);
	for(int k=0;k<m_size;k++) 	{m_beta_vector[k]=beta_vector->Elt(k);}

	computebarrier();
}



ICM_FFTDISTRIB::ICM_FFTDISTRIB(int size,std::vector<double> pdef,std::vector<double>& beta_vector):
	m_size	(size),
	m_flag	(false),
	m_complexe (NULL), 
	m_cosx	(NULL), 
	m_sinx	(NULL), 
	m_datas	(NULL), 
	m_b		(NULL),
	m_g		(NULL),
	m_f		(NULL),
	m_pdef_perturb   (NULL),
	m_pdef   (NULL)

{ 
	init();
	int pdef_size = pdef.size();
	double* p_def = new double[pdef_size];
	for (int i=0;i<pdef_size;i++)
		p_def[i] = pdef[i];

	setpdef(p_def);

	m_beta_vector.resize(m_size);
	for(int k=0;k<m_size;k++) 	{m_beta_vector[k]=beta_vector[k];}

	computebarrier();
}
//----------------------------------------------------------------------------------------------------------------------
ICM_FFTDISTRIB::~ICM_FFTDISTRIB()
{
	if(m_complexe !=NULL) {delete [] m_complexe; m_complexe = NULL;}

	clear();
}
//----------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------
void
ICM_FFTDISTRIB::clear()
{
	if(m_cosx		!=NULL) {delete [] m_cosx;	m_cosx = NULL;}
	if(m_sinx		!=NULL) {delete [] m_sinx;	m_sinx = NULL;}
	if(m_datas		!=NULL) {delete [] m_datas; m_datas = NULL;}
	if(m_b			!=NULL) {delete [] m_b;		m_b = NULL;}
	if(m_g			!=NULL) {delete [] m_g;		m_g = NULL;}
	if(m_f			!=NULL) {delete [] m_f;		m_f = NULL;}
	if(m_pdef_perturb	!=NULL) {delete [] m_pdef_perturb; m_pdef_perturb = NULL;}
	if(m_pdef	!=NULL) {delete [] m_pdef; m_pdef = NULL;}


}
//----------------------------------------------------------------------------------------------------------------------
void
ICM_FFTDISTRIB::init()
{
	clear();
	m_dens.clear();
	m_dens_perturb.clear();
	
	m_b	= new double [m_size];
	m_pdef_perturb	= new double [m_size];

	m_pdef	= new double [m_size];
	m_f	= new double [m_size*20];

	if(m_complexe == NULL) m_complexe = new double [2];

	int k = 1;
	while(pow(2.,k)<m_size+1) k++;
	m_nbPoints = (int)pow(2.,k-1); // Resolution is 2*m_nbpoints, we use symetry
	m_normNbPoints = 1./((double) 2.*m_nbPoints);

	m_sinx		= new double [m_nbPoints+1];
	m_cosx		= new double [m_nbPoints+1];
	m_datas		= new double [1+m_nbPoints*4];
	m_g	= new double [(m_nbPoints+1)*2*20];
		
	// Discretisation initialization
	for(k=0;k<=m_nbPoints;k++) {
		double x = m_2Pi * k * m_normNbPoints;
		m_sinx[k] = sin(x);
		m_cosx[k] = cos(x);
	}



}
//----------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------
void
ICM_FFTDISTRIB::computebarrier()
{
	m_flag = false;
	for(int k=0;k<m_size;k++) 
	{
		if (m_pdef[k] == 0)  {m_b[k] = -10;}
		else {m_b[k] = ep::MathSrv::invCumNorm(m_pdef[k]);}
	}
}


void
ICM_FFTDISTRIB::setpdef(double* pdef)
{
	m_flag = false;
	for(int k=0;k<m_size;k++) 
	{
		if (pdef[k]<0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
            "ERROR : Proba de defaut négative ");
		}
		else
			m_pdef[k] = pdef[k];
	}
}

void
ICM_FFTDISTRIB::setpdefperturb(double* pdef_perturb)
{
	for(int k=0;k<m_size;k++) m_pdef_perturb[k] = pdef_perturb[k];
}

//----------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------
bool
ICM_FFTDISTRIB::isComputed(const long& index) const
{
	std::map<long, std::vector<double> >::const_iterator it = m_dens_perturb.find(index);
	if(it==m_dens_perturb.end()) return false;
	else return true;
}
//----------------------------------------------------------------------------------------------------------------------
void
ICM_FFTDISTRIB::computeDensity()
{	
	for(int k=0;k<=m_nbPoints;k++) 
	{
		m_ind_x = k;    // variable avec laquelle est calculée la fctn caractéristique
		m_ind_t = -1;   // variable d'hermitte (de la gausienne)
		ComplexHermiteIntegration(ff1::const_mem_call(&ICM_FFTDISTRIB::getValue, (*this) )).Integrate(m_datas[2*k+1], m_datas[2*k+2]); 
	}
		
	for(k=1;k<m_nbPoints;k++)
	{
			m_datas[2*(k+m_nbPoints)+1] = m_datas[2*(m_nbPoints-k)+1];
			m_datas[2*(k+m_nbPoints)+2] = -m_datas[2*(m_nbPoints-k)+2];
	}

	// Compute the fft inverse:
	fourier::fft(m_datas, 2*m_nbPoints, -1);

	// Fill density vector
	m_dens.resize(m_size+1);
	for( k=0;k<=m_size;k++) m_dens[k] = m_normNbPoints*m_datas[1+2*k];

	m_flag = true;
}
//----------------------------------------------------------------------------------------------------------------------
const std::vector<double>&
ICM_FFTDISTRIB::getDensity()	const
{
	if(!m_flag) ((ICM_FFTDISTRIB&)*this).computeDensity();
	return m_dens;
}
//----------------------------------------------------------------------------------------------------------------------
void
ICM_FFTDISTRIB::computeDensityPerturbed(const long& index)
{
	if(!m_flag) computeDensity();

	// Rajout olivier le 7/05 (le if)
	if (m_pdef_perturb[index] == 0) {m_b_perturb = -10;}
	else  {m_b_perturb = ep::MathSrv::invCumNorm(m_pdef_perturb[index]);}


	m_index = index;
			for(int k=0;k<=m_nbPoints;k++) {
				m_ind_x = k;
				m_ind_t = -1;
				ComplexHermiteIntegration(ff1::const_mem_call(&ICM_FFTDISTRIB::getValuePerturbed, (*this) )).Integrate(m_datas[2*k+1], m_datas[2*k+2]); 
			}
			for(k=1;k<m_nbPoints;k++) {
				m_datas[2*(k+m_nbPoints)+1] = m_datas[2*(m_nbPoints-k)+1];
				m_datas[2*(k+m_nbPoints)+2] = -m_datas[2*(m_nbPoints-k)+2];
			}
	// Compute the fft inverse:

	fourier::fft(m_datas, 2*m_nbPoints, -1);

		// Fill density vector
	std::vector<double>& vec = m_dens_perturb[index];
	vec.resize(m_size+1);
	for( k=0;k<=m_size;k++) vec[k] = m_normNbPoints*m_datas[1+2*k];

}
//----------------------------------------------------------------------------------------------------------------------
const std::vector<double>&
ICM_FFTDISTRIB::getDensityPerturbed(const long& index) const
{
	std::map<long, std::vector<double> >::const_iterator it = m_dens_perturb.find(index);
	return m_dens_perturb.find(index)->second;
}
//----------------------------------------------------------------------------------------------------------------------
void
ICM_FFTDISTRIB::addDensityPerturbed(const long& index, const std::vector<double>& ref)
{
	m_dens_perturb[index] = ref;
}
//----------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
double*
ICM_FFTDISTRIB::getValue(double t) const
{
	m_ind_t +=1;
	int ind=0;

	m_complexe[0] = m_normPi;
	m_complexe[1] = 0.;
	double seuil=0.;

	for(int k=0;k<m_size;k++) {

		 seuil	= m_b[k]-m_beta_vector[k]*m_sqrt2*t;

		seuil /= sqrt(1.-m_beta_vector[k]*m_beta_vector[k]);

	//	double fk= ep::MathSrv::cumNorm(seuil);
		double fk = NAG_cumul_normal(seuil);
		
		 m_f[m_ind_t*m_size+k] = fk;

		double re =1.+(m_cosx[m_ind_x]-1.)*fk; //double re = 1+(1.-m_cosx[m_ind_x])*fk;
		double im = m_sinx[m_ind_x]*fk; //double im =-m_sinx[m_ind_x]*(fk);
		
		double tmp = m_complexe[0];
		// Update Real part
		m_complexe[0] *= re;
		m_complexe[0] -= im*m_complexe[1];
		 
		// Update Imaginary part
		m_complexe[1] *= re;
		m_complexe[1] += im*tmp;
	}

	 ind = m_ind_x*40 + m_ind_t*2;

	m_g[ind]	= m_complexe[0];
	m_g[ind+1]	= m_complexe[1];

	return m_complexe;
}

double*
ICM_FFTDISTRIB::getValuePerturbed(double t)	const
{
	m_ind_t += 1;
	double seuil=0.,fk=0.;

			 
		fk = m_f[m_ind_t*m_size+m_index];
		double ratio = 1.+(2.*fk*(m_cosx[m_ind_x] - 1.))+(2*fk*fk*(1.-m_cosx[m_ind_x]));
		seuil=0.;

		if(ratio==0.) { // |g_k(x)|=0 , low computing

			m_complexe[0] = m_normPi;
			m_complexe[1] = 0.;
			
			for(int k=0;k<m_size;k++) {

				if(k!=m_index) fk = m_f[m_ind_t*m_size+k];
				else {
					seuil	= m_b_perturb-m_beta_vector[k]*m_sqrt2*t;
					seuil /= sqrt(1.-m_beta_vector[k]*m_beta_vector[k]);
					//fk= ep::MathSrv::cumNorm(seuil);
					fk = NAG_cumul_normal(seuil);
				}

				double re =1.+(m_cosx[m_ind_x]-1.)*fk; //double re = 1+(1.-m_cosx[m_ind_x])*fk;
				double im = m_sinx[m_ind_x]*fk; //double im =-m_sinx[m_ind_x]*(fk);
				double tmp = m_complexe[0];
				// Update Real part
				m_complexe[0] *= re;
				m_complexe[0] -= im*m_complexe[1];
				// Update Imaginary part
				m_complexe[1] *= re;
				m_complexe[1] += im*tmp;
			}

		} 
		else { // fast Computing
			ratio = 1./ratio;
			int ind = m_ind_x*40 + m_ind_t*2;
			m_complexe[0] = m_g[ind];
			m_complexe[1] = m_g[ind + 1];

			seuil = m_b_perturb-m_beta_vector[m_index]*m_sqrt2*t;
			seuil /= sqrt(1.-m_beta_vector[m_index]*m_beta_vector[m_index]);

		//	double fk_prime= ep::MathSrv::cumNorm(seuil);
			double fk_prime = NAG_cumul_normal(seuil);

		//	double re = (2.*fk*fk_prime-fk-fk_prime)*(1.-m_cosx[m_ind_x]) + 1.;

			double re = (2.*fk*fk_prime-fk-fk_prime)*(1.-m_cosx[m_ind_x]) + 1.;
			double im = -(fk-fk_prime) * m_sinx[m_ind_x];

			// Update g value
			double tmp = m_complexe[0];
			// Update Real part
			m_complexe[0] *= re;
			m_complexe[0] -= im*m_complexe[1];
			m_complexe[0] *= ratio;
			// Update Imaginary part
			m_complexe[1] *= re;
			m_complexe[1] += im* tmp;
			m_complexe[1] *= ratio;
		}

	return m_complexe;
}

double
ICM_FFTDISTRIB::compute_expectedlosstranche(double tranche_up, 
									    double tranche_down, 
									    double lossunit)
{

	//	lossunit = its_LossRate[0];

	double exp_loss_tranche=0.;
	int lup=(int)floor(tranche_up/lossunit);
	int ldown=(int)floor(tranche_down/lossunit);
	//int i=0;	

	std::vector<double> distrib;
	distrib.resize(m_size+1);

	distrib=getDensity();
	double tmp=0.;
	
	// Cas ou la lossunit est trop petite et telle que T_up ne soit jamais atteinte
	int indice_max = lup;
	if (lup > m_size) indice_max = m_size;

	for(int l=0;l<=indice_max;l++) 	{	tmp+=distrib[l];	}

	double taildistrib=1.-tmp;

	for( l=ldown+1;l<=indice_max;l++) 
		{exp_loss_tranche+=distrib[l]*(l*lossunit-tranche_down);}

	exp_loss_tranche+=(tranche_up-tranche_down)*taildistrib;


return exp_loss_tranche;
}


double
ICM_FFTDISTRIB::compute_expectedlosstranche_perturb(double tranche_up, 
											    double tranche_down, 
											    double lossunit,
											    long index)
{
	double exp_loss_tranche_perturb=0.;
	int lup=(int)floor(tranche_up/lossunit);
	int ldown=(int)floor(tranche_down/lossunit);
	

	std::vector<double> distrib;
	distrib.resize(m_size+1);

		computeDensityPerturbed(index);
	distrib=getDensityPerturbed(index);

	double tmp=0.;

	int indice_max = lup;
	if (lup > m_size) indice_max = m_size;

	for(int l=0;l<=indice_max;l++) 	{	tmp+=distrib[l];	}

	double taildistrib=1.-tmp;

	for( l=ldown+1;l<=lup;l++) 	{exp_loss_tranche_perturb+=distrib[l]*(l*lossunit-tranche_down);}

	exp_loss_tranche_perturb+=(tranche_up-tranche_down)*taildistrib;


return exp_loss_tranche_perturb;
}
