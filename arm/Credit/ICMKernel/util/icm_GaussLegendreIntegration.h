
#ifndef _F1_GAUSSLEGENDREINTEGRATION_H_
#define _F1_GAUSSLEGENDREINTEGRATION_H_

#include "icm_functors.h"
//#include <double.h>
#include <math.h>

// ----------------------------------------------------------------------------------
//	class		ICM_GaussLegendreIntegration
//  author		L. Jacquel
//	version		1.0
//	date		September 2004
//	file		ICM_GaussLegendreIntegration.h

//	brief		Definition of a Gauss Legendre Integrator
// ----------------------------------------------------------------------------------

// Compute integral
//     /- x2
//	  / 
//   /    f(x) dx
//  /
//-/x1
// with a Gauss-Legendre approximation at degree n


#include "icm_gaussian.h"



template <class F> class GaussLegendreIntegration_t
{
	public:
		GaussLegendreIntegration_t(F f, double x_min = 0.0, double x_max = 1.0):m_f(f){xmin = x_min; xmax = x_max;}
		virtual ~GaussLegendreIntegration_t() {}

		double	Integrate_11();
		double	Integrate_21();
		double	Integrate_51();
		double	Integrate_101();

		double	xmin;
		double	xmax;
	private :
		F		m_f; 
		static double GaussLegendreCoeff_11[11][2];
		static double GaussLegendreCoeff_21[21][2];
		static double GaussLegendreCoeff_51[51][2];
		static double GaussLegendreCoeff_101[101][2];
};

//----------------------------------------------------------------------------
//	Helper. 
template <class F> inline GaussLegendreIntegration_t<F> GaussLegendreIntegration(F f, double x_min, double x_max)
{ return GaussLegendreIntegration_t<F>(f, x_min, x_max) ; }
//----------------------------------------------------------------------------

template <class F> double GaussLegendreIntegration_t<F>::GaussLegendreCoeff_11[11][2] = { 
{0.010885670927, 0.027834283558},
{0.056468700116, 0.062790184732},
{0.134923997213, 0.093145105463},
{0.240451935397, 0.116596882296},
{0.365228422024, 0.131402272255},
{0.500000000000, 0.136462543389},
{0.634771577976, 0.131402272255},
{0.759548064603, 0.116596882296},
{0.865076002787, 0.093145105463},
{0.943531299884, 0.062790184732},
{0.989114329073, 0.027834283558}
};

template <class F> double GaussLegendreIntegration_t<F>::GaussLegendreCoeff_21[21][2] = {
	{0.003123945690,	0.008008614129},
	{0.016386580717,	0.018476894870},
	{0.039950332925,	0.028567212713},
	{0.073318317708,	0.038050056814},
	{0.115780018262,	0.046722211728},
	{0.166430597901,	0.054398649584},
	{0.224190582056,	0.060915708027},
	{0.287828939896,	0.066134469317},
	{0.355989341599,	0.069943697396},
	{0.427219072920,	0.072262201995},
	{0.500000000000,	0.073040566825},
	{0.572780927080,	0.072262201995},
	{0.644010658401,	0.069943697396},
	{0.712171060104,	0.066134469317},
	{0.775809417944,	0.060915708027},
	{0.833569402099,	0.054398649584},
	{0.884219981738,	0.046722211728},
	{0.926681682292,	0.038050056814},
	{0.960049667075,	0.028567212713},
	{0.983613419283,	0.018476894870},
	{0.996876085310,	0.008008614129}
};


template <class F> double GaussLegendreIntegration_t<F>::GaussLegendreCoeff_51[51][2] = { 
{0.000545004576, 0.001398403586},
{0.002869369782, 0.003250168889},
{0.007042004132, 0.005092595649},
{0.013048315990, 0.006916317003},
{0.020866075693, 0.008714357362},
{0.030466227999, 0.010479994201},
{0.041813068845, 0.012206650287},
{0.054864390985, 0.013887899297},
{0.069571644409, 0.015517485645},
{0.085880118088, 0.017089346602},
{0.103729143950, 0.018597634462},
{0.123052322757, 0.020036738143},
{0.143777771211, 0.021401303999},
{0.165828389412, 0.022686255704},
{0.189122147700, 0.023886813120},
{0.213572391824, 0.024998510075},
{0.239088165317, 0.026017210968},
{0.265574547857, 0.026939126157},
{0.292933008387, 0.027760826048},
{0.321061771656, 0.028479253860},
{0.349856196832, 0.029091736991},
{0.379209166776, 0.029595996961},
{0.409011486521, 0.029990157888},
{0.439152289491, 0.030272753467},
{0.469519449925, 0.030442732422},
{0.500000000000, 0.030499462421},
{0.530480550075, 0.030442732422},
{0.560847710509, 0.030272753467},
{0.590988513479, 0.029990157888},
{0.620790833224, 0.029595996961},
{0.650143803168, 0.029091736991},
{0.678938228344, 0.028479253860},
{0.707066991613, 0.027760826048},
{0.734425452143, 0.026939126157},
{0.760911834683, 0.026017210968},
{0.786427608176, 0.024998510075},
{0.810877852300, 0.023886813120},
{0.834171610588, 0.022686255704},
{0.856222228789, 0.021401303999},
{0.876947677243, 0.020036738143},
{0.896270856050, 0.018597634462},
{0.914119881912, 0.017089346602},
{0.930428355591, 0.015517485645},
{0.945135609015, 0.013887899297},
{0.958186931155, 0.012206650287},
{0.969533772001, 0.010479994201},
{0.979133924307, 0.008714357362},
{0.986951684010, 0.006916317003},
{0.992957995868, 0.005092595649},
{0.997130630218, 0.003250168889},
{0.999454995424, 0.001398403586},
};

template <class F> double GaussLegendreIntegration_t<F>::GaussLegendreCoeff_101[101][2] = { 
{0.000140330235, 0.000360115853},
{0.000739244005, 0.000837946379},
{0.001816133605, 0.001315682437},
{0.003370206906, 0.001792219699},
{0.005400004562, 0.002267052505},
{0.007903589185, 0.002739717342},
{0.010878564857, 0.003209759135},
{0.014322082721, 0.003676726837},
{0.018230844583, 0.004140172813},
{0.022601106364, 0.004599652970},
{0.027428681798, 0.005054727090},
{0.032708946491, 0.005504959204},
{0.038436842377, 0.005949918001},
{0.044606882576, 0.006389177230},
{0.051213156654, 0.006822316103},
{0.058249336287, 0.007248919699},
{0.065708681333, 0.007668579360},
{0.073584046281, 0.008080893077},
{0.081867887102, 0.008485465885},
{0.090552268475, 0.008881910228},
{0.099628871390, 0.009269846341},
{0.109089001115, 0.009648902609},
{0.118923595525, 0.010018715922},
{0.129123233787, 0.010378932025},
{0.139678145379, 0.010729205858},
{0.150578219456, 0.011069201881},
{0.161813014529, 0.011398594404},
{0.173371768470, 0.011717067892},
{0.185243408824, 0.012024317270},
{0.197416563410, 0.012320048215},
{0.209879571220, 0.012603977438},
{0.222620493588, 0.012875832955},
{0.235627125624, 0.013135354347},
{0.248887007907, 0.013382293013},
{0.262387438422, 0.013616412402},
{0.276115484726, 0.013837488246},
{0.290057996333, 0.014045308768},
{0.304201617320, 0.014239674892},
{0.318532799111, 0.014420400429},
{0.333037813462, 0.014587312257},
{0.347702765611, 0.014740250487},
{0.362513607585, 0.014879068614},
{0.377456151660, 0.015003633660},
{0.392516083951, 0.015113826302},
{0.407678978123, 0.015209540983},
{0.422930309212, 0.015290686014},
{0.438255467535, 0.015357183665},
{0.453639772691, 0.015408970235},
{0.469068487619, 0.015445996118},
{0.484526832718, 0.015468225844},
{0.500000000000, 0.015475638120},
{0.515473167282, 0.015468225844},
{0.530931512381, 0.015445996118},
{0.546360227309, 0.015408970235},
{0.561744532465, 0.015357183665},
{0.577069690788, 0.015290686014},
{0.592321021877, 0.015209540983},
{0.607483916049, 0.015113826302},
{0.622543848340, 0.015003633660},
{0.637486392415, 0.014879068614},
{0.652297234389, 0.014740250487},
{0.666962186538, 0.014587312257},
{0.681467200889, 0.014420400429},
{0.695798382680, 0.014239674892},
{0.709942003667, 0.014045308768},
{0.723884515274, 0.013837488246},
{0.737612561578, 0.013616412402},
{0.751112992093, 0.013382293013},
{0.764372874376, 0.013135354347},
{0.777379506412, 0.012875832955},
{0.790120428780, 0.012603977438},
{0.802583436590, 0.012320048215},
{0.814756591176, 0.012024317270},
{0.826628231530, 0.011717067892},
{0.838186985471, 0.011398594404},
{0.849421780544, 0.011069201881},
{0.860321854621, 0.010729205858},
{0.870876766213, 0.010378932025},
{0.881076404475, 0.010018715922},
{0.890910998885, 0.009648902609},
{0.900371128610, 0.009269846341},
{0.909447731525, 0.008881910228},
{0.918132112898, 0.008485465885},
{0.926415953719, 0.008080893077},
{0.934291318667, 0.007668579360},
{0.941750663713, 0.007248919699},
{0.948786843346, 0.006822316103},
{0.955393117424, 0.006389177230},
{0.961563157623, 0.005949918001},
{0.967291053509, 0.005504959204},
{0.972571318202, 0.005054727090},
{0.977398893636, 0.004599652970},
{0.981769155417, 0.004140172813},
{0.985677917279, 0.003676726837},
{0.989121435143, 0.003209759135},
{0.992096410815, 0.002739717342},
{0.994599995438, 0.002267052505},
{0.996629793094, 0.001792219699},
{0.998183866395, 0.001315682437},
{0.999260755995, 0.000837946379},
{0.999859669765, 0.000360115853},
};

//----------------------------------------------------------------------------


template <class F>
double
GaussLegendreIntegration_t<F>::Integrate_11()
{
	double ret(0.);
	double xi;
	double	size	=	xmax - xmin;

	for (int i=0;i<11;i++)
	{
		xi	=	xmin + size * GaussLegendreCoeff_11[i][0];
		ret += m_f(xi) * GaussLegendreCoeff_11[i][1] * size * exp(- 0.5 * xi * xi);
	}
	return ret;
}

template <class F>
double
GaussLegendreIntegration_t<F>::Integrate_21()
{
	double ret(0.);
	double xi;
	double	size	=	xmax - xmin;
	for (int i=0;i<21;i++)
	{
		xi	=	xmin + size * GaussLegendreCoeff_21[i][0];
		ret += m_f(xi) * GaussLegendreCoeff_21[i][1] * size * exp(- 0.5 * xi * xi);
	}
	return ret;
}


template <class F>
double
GaussLegendreIntegration_t<F>::Integrate_51()
{
	double ret(0.);
	double xi;
	double	size	=	xmax - xmin;
	
	for (int i=0;i<51;i++)
	{
		xi	=	xmin + size * GaussLegendreCoeff_51[i][0];
		ret += m_f(xi) * GaussLegendreCoeff_51[i][1] * size * exp(- 0.5 * xi * xi);
	}
	return ret;
}

template <class F>
double
GaussLegendreIntegration_t<F>::Integrate_101()
{
	double ret(0.);
	double xi;
	double	size	=	xmax - xmin;

	for (int i=0;i<101;i++)
	{
		xi	=	xmin + size * GaussLegendreCoeff_101[i][0];
		ret += m_f(xi) * GaussLegendreCoeff_101[i][1] * size * exp(- 0.5 * xi * xi);
	}
	return ret;
}

//----------------------------------------------------------------------------

#endif




/*
#define EPS_GL 3.0e-11

template <class F> class GaussLegendreIntegration_t
{
public:

	GaussLegendreIntegration_t(F f, int n_step, double x1_bound, double x2_bound):m_f(f){Init();n = n_step; x1 = x1_bound; x2 = x2_bound;}
	
	virtual ~GaussLegendreIntegration_t() {Reset();}
	double Integrate();

	void	ComputeAbscissasAndWeights(double x1, double x2, double* x_out, double* w_out, int n);
	void	ComputeAbscissasAndWeights(double x1, double x2, int n);

	void	SetGaussLegendreDegree(const int n_deg) {n = n_deg;}
	int		GetGaussLegendreDegree() {return n;}

	void	GetGaussLegendreAbscissasAndWeights(double* x_out, double* w_out);
	void	GetGaussLegendreAbscissasAndWeights(double* x_out, double* w_out, int& n_out);

	void	ResetArrays();
	void	ResizeArrays(int	n_deg);

	double test(/*std::string& errStr*///);
/*
private:
	void	Init();
	void	Reset();

private :
	int		n;
	double	x1;
	double	x2;

	double*	x;
	double*	w;

	F		m_f;
};

template <class F>
void
GaussLegendreIntegration_t<F>::Init()
{
	n = 21;
	x1 = 0.0;
	x2 = 1.0;
	x	=	NULL;
	w	=	NULL;

	double	the_val	=	test();
}

template <class F>
void
GaussLegendreIntegration_t<F>::Reset()
{
	ResetArrays();
}

//----------------------------------------------------------------------------
//	Helper. 
//template <class F> inline GaussLegendreIntegration_t<F> GaussLegendreIntegration()
//{ return GaussLegendreIntegration_t<B> ; }

//	Helper. 
template <class F> inline GaussLegendreIntegration_t<F> GaussLegendreIntegration(F f, int n_step, double x1_bound, double x2_bound)
{ return GaussLegendreIntegration_t<F>(f, n_step, x1_bound, x2_bound) ; }

//----------------------------------------------------------------------------

template <class F>
double
GaussLegendreIntegration_t<F>::Integrate()
{
	double ret(0.);

	ComputeAbscissasAndWeights(x1, x2, n);

	for(int i=0;i<n;i++) ret += m_f(x[i])*w[i];
	return ret;
}
//----------------------------------------------------------------------------

template <class F>
void
GaussLegendreIntegration_t<F>::ResizeArrays(int	n_deg)
{
	n	=	n_deg;
	if (x)
		delete[] x;
	x	=	new	double[n];

	if (w)
		delete[] w;
	w	=	new	double[n];

}

//----------------------------------------------------------------------------

template <class F>
void
GaussLegendreIntegration_t<F>::ResetArrays()
{
	if (x)
		delete[] x;

	if (w)
		delete[] w;
}

//----------------------------------------------------------------------------

template <class F>
double
GaussLegendreIntegration_t<F>::test(/*std::string& errStr*///)
//{
//	return 0.0;
//	bool ret(true);
	
//	errStr = "Testing GaussLegendreIntegration";

/*
	int NGL	=	101;

	double	x1	=	0.0;
	double	x2	=	1.0;

	double*	the_x	=	NULL;
	double*	the_w	=	NULL;

	ComputeAbscissasAndWeights(x1, x2, the_x, the_w, NGL);

	// ----------------------------------------------------
	
	// Integration limits
	vmin = -6;   vmax =  6;

	int i;
	double	xi, f_xi;
	double	interval_size	=	(x2 - x1);

	double	The_Integral	=	0.0;

	for (i=0;i<NGL;i++)
	{
		xi		=	interval_size * the_x[i] + x1;
		//	the function as of xi
		f_xi	=	1;
		coef	=	StdNorDensity(xi) * the_w[i] * interval_size;

		The_Integral	+=	coef * f_xi;
	}
	// ----------------------------------------------------

	return	The_Integral;
*/
/*	B b;
	double val = GaussLegendreIntegration_t(ff1::mem_call(&B::f, b), 51, 0.0, 1.0).Integrate();	
	return val;
//	SDMTEST((val>1.77245385)&&(val<1.77245386))
*/
//	return ret;
//}
/*
template <class F>
void
GaussLegendreIntegration_t<F>::ComputeAbscissasAndWeights(double x1, double x2, double* x_out, double* w_out, int n)
{
	ComputeAbscissasAndWeights(x1, x2, n);
	GetGaussLegendreAbscissasAndWeights(x_out, w_out);
}

template <class F>
void
GaussLegendreIntegration_t<F>::ComputeAbscissasAndWeights(double x1, double x2, int n)
//Given the lower and upper limits of integration x1 and x2, and given n, this routine returns
//arrays x[0..n-1] and w[0..n-1] of length n, containing the abscissas and weights of the Gauss-
//Legendre n-point quadrature formula.
{
	int m,j,i;
	double z1,z,xm,xl,pp,p3,p2,p1;

	if (n <= 0) return;
	ResizeArrays(n);

	//	High precision is a good idea for this routine.
	
	m=(n+1)/2;
	//	The roots are symmetric in the interval, so we only have to .nd half of them.
	xm=0.5*(x2+x1);
	xl=0.5*(x2-x1);
	for (i=1;i<=m;i++)
	{
		//	Loop over the desired roots.
		z=cos(3.141592654*(i-0.25)/(n+0.5));
		//	Starting with the above approximation to the ith root, we enter the main loop of
		//	refinement by Newton’s method.
		do
		{
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++)
			{
				//	Loop up the recurrence relation to get the Legendre polynomial evaluated at z.
				p3=p2;
				p2=p1;
				p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
			}

			//	p1 is now the desired Legendre polynomial. We next compute pp, its derivative,
			//	by a standard relation involving also p2, the polynomial of one lower order.

			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp;
			//	Newton’s method.
		} while (fabs(z-z1) > EPS_GL);

		x[i-1]=xm-xl*z;
		//	Scale the root to the desired interval,
		x[n-i]=xm+xl*z;
		//	and put in its symmetric counterpart.
		w[i-1]=2.0*xl/((1.0-z*z)*pp*pp);
		//	Compute the weight
		w[n-i]=w[i-1];
		//	and its symmetric counterpart.
	}
}


template <class F>
void
GaussLegendreIntegration_t<F>::GetGaussLegendreAbscissasAndWeights(double* x_out, double* w_out)
{
	x_out	=	x;
	w_out	=	w;
}

template <class F>
void
GaussLegendreIntegration_t<F>::GetGaussLegendreAbscissasAndWeights(double* x_out, double* w_out, int& n_out)
{
	n_out	=	n;
	x_out	=	x;
	w_out	=	w;
}

#endif	// _F1_GAUSSLEGENDREINTEGRATION_H_

/*
//----------------------------------------------------------------------------
// For test purpose
class B {
	public: 
	B(){};
	~B(){};
	double f(double x) {return StdNorDensity(x);}
};
//----------------------------------------------------------------------------
*/