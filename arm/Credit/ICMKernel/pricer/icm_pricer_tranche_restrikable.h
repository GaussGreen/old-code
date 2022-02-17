#ifndef _ICM_PRICER_RESTRIKABLE_CDO_H
#define _ICM_PRICER_RESTRIKABLE_CDO_H

#include <nag.h>                      
#include <nags.h>                    
#include <nagg01.h>                 
#include <nagg05.h> 


#include "ICMKernel\pricer\icm_pricer.h"
//#include "ICMKernel\util\icm_matrix.h"


/*********************************************************************************/
/*  \class  ICM_Pricer_Tranche_Restrikable and all usefull 
/*  \classes/utils for MC pricing of this product
/*  \author N.ALAYA
/*	\version 1.0
/*	\date   April 2007
/***********************************************************************************/

class ARM_Model; 
class ARM_Security; 

class DiffStochIntensityHetero;
class GaussHermite;

class ICM_Pricer_Tranche_Restrikable : public ICM_Pricer
{

private:

	int		itsNbSimul ;
	double	itsTimeStep ;

	double	itsSigma ;		
	double	itsMRS ;	

	double	itsBeta0;		
	double	itsBeta1;	
	double	itsBeta2;	
	double	itsBeta3;

	int		itsIsHetero;	

	double	itsRecov;				
	double	itsAvgSpreadFreq ;

	ARM_Vector* itsDatesYF;		
	ARM_Vector* itsprobas;	

	double	itsInitialSpread; 
	bool	itsInitialSpreadFlag ;

	double itsELCoeff ;
	double itsLossCoeff ;
	int itsUsePeriodIsEL;

	double	itsFeeLegUnit ; // Feeleg = feelegUnit * initspread 
	bool	itsFeeLegFlg ;
	
public:

	void Init(void);

	void Set(ARM_Security* sec , ARM_Model* mod , const ICM_Parameters& params , const ARM_Date& asof );


	ICM_Pricer_Tranche_Restrikable (void) 
	{
		Init();
	}

	~ICM_Pricer_Tranche_Restrikable(void);


	ICM_Pricer_Tranche_Restrikable(ARM_Security* sec , ARM_Model* mod , const ICM_Parameters& params , const ARM_Date& asof)
	{
		Init();
		Set(sec,mod,params,asof);
	}

	virtual double FeeLegPV () ;
	virtual double DefLegPV () ;
	virtual double ComputeSpread() ;
	
	virtual double Accrued() ;
	virtual double ComputeDuration(void) ;
	virtual double ComputeSpread(const double& MtM = 0.) ;
 	virtual	double ComputeImpliedVol(const double& Price) ;



	void ResetPricer(void)  ;
	void View(char* id = NULL, FILE* ficOut = NULL) ;

protected:

	virtual double ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
		const std::string& plot,
		const std::string& label,
		double epsilon, double epsilonGamma =0 ) ; 

	//virtual void DoPriceVector(const std::string& measure) ;
	//virtual double DoPrice(const std::string& measure); 

private:
	void PriceTranche(double &fl, double &dl,double &K1, double&K2, double &t , double &xt ,double &beta, DiffStochIntensityHetero &diffusion , GaussHermite &gaussHermite ,int &nbNoms, double &lgd, double** probloss_v,int numK, int adj =0 );
	double GetPrice(int legtype, double SpreadInit);
	double GetInitialSpread();// First to calculate 
	void Compute() ;
	double CptDuration() ;
	double GetFeeLegUnit() ;
	void ComputeInitialProbas();
	void ComputeBetas();
	void SetSigmaFromModel();

}; 

#endif

/* __________________________________________________________________________________________________________*/
 

double Uniform(double &seed);

double Normal(double &seed);

void GenerateVar(double lambda, double dt, double &seed, double &W1, double& W2, double& dY);

double GenerateDefaut(double lambda, double dt, double &seed);
double GenerateDefaut2(double prodadef, double &seed);


void CptProbaLossHomog(double ** &probaloss_v, double * &proba_pf_v, double &prob_noms, double* &vsyst, 
				  int &degreeHerm, double &beta, int &nbNoms, int nbDefmax);

void CptProbaLossHeterog(double ** &probaloss_v, double ** &proba_noms_v, double *&prob_noms, double* &vsyst, 
				  int &degreeHerm, double &beta, int &nbNoms, int nbDefmax);


/* __________________________________________________________________________________________________________*/
 




/* __________________________________________________________________________________________________________*/

//////////////////////////////
//
//		Quelques constantes mathématiques et autres ...
//
//////////////////////////////

#define _last		-1001

//#define PI			3.1415926535897932			// ... PI

#define SQRT2		1.4142135623730951			// = sqrt(2.)

#define SQRT2PI		2.5066282746310002			// = sqrt(2 * PI)

#define INVSQRT2PI	0.3989422804014327			// = 1 / sqrt(2 * PI)

#define INVSQRTPI	0.56418958354775628			// = 1 / sqrt(PI)

#define INFTY		1.e+20						// bonne approximation de l'infini

#define DEPS		1.e-12						// précision sur les doubles pour l'opérateur == 

#define BPToPct		0.0015811388301				// = sqrt(250)/10000 (passage des vols de bp en %)

static const double	DNOVALUE(-999.);

static const int	INOVALUE(-999);

/* __________________________________________________________________________________________________________*/


class R2RFunction
{
public:
	virtual double operator()(double x) const = 0;

	virtual double derivative(double x,double eps = 1.e-4) const
	{
		return (this->operator ()(x + eps) - this->operator ()(x - eps)) / (2. * eps);
	}
};



/////////////////////////////////////
//
//		définition de quelques fonctions mathématiques simples 
//
/////////////////////////////////////


#define	DMIN(a,b)		 a > b ? b : a

inline double	DMAX(double a,double b)				{return a > b ? a : b;};



inline int		IMAX(int p,int q)					{return p > q ? p : q;};

inline int		IMIN(int p,int q)					{return p > q ? q : p;};

//inline double	SQR(double x)						{return x * x;};

inline double	POW(double x,double y)				{return x == 0. && y == 0. ? DNOVALUE : x == 1. || y == 0. ? 1. : x == 0. ? 0. : exp(y * log( fabs(x) ) );};

//inline int		SIGN(double x)						{return x >= 0. ? 1 : -1;};

inline int		SIGN_(int p)							{return p >= 0 ? 1 : -1;};

inline double	SIGN_(double a,double b)				{return b >= 0. ? fabs(a) : - fabs(a);};

//inline double	DIV(double x,double y)				{return y == 0. ? SIGN(x) * INFTY : x / y;};

//inline int		ROUND(double x)						{return x - int(x) > 0.5 ? int(x)+1 : int(x);};

inline int		MOD(int p,int q)					{return q == 0 ? 0 : p % q;};

inline bool		ODD(int p)							{return ( p % 2 == 0 );};

inline long		FAC(long p)							{return p == 0 ? 1 : p * FAC(p - 1);};

inline int		PGCD(int p,int q)					{return q == 0 ? p : PGCD(q,p % q);};

inline double	PYTHAG(double a,double b)			{return sqrt(a * a + b * b);};

/////////////
// précision de la machine
/////////////
inline double	machprec()
{
	double z,tsmall = 1.;
	do
	{
		tsmall *= 0.5;
		z		= tsmall + 1.;
	} while (z > 1.);

	return 4.*tsmall;
}

///////////////
//	fonctions simples avec un nombre d'arguments optionel
///////////////

#define iEndArgList		-9999



/* __________________________________________________________________________________________________________*/


enum ePOLYNOMIALS	{_hermite12,_hermite20,_hermite40,_hermite80,_nohermite,
						_legendre8,_legendre12,_legendre20,_legendre40,_legendre80,_nolegendre};

enum eNUMINTEGRAL	{_simpson,_legendre,_hermite,_romberg,_nonumint};
class GaussHermite
{
private :
	ePOLYNOMIALS _enumDegree;
	eNUMINTEGRAL _typeInt;
	int _degree;
	double *_x;
	double *_w;
public :
	GaussHermite(int degree,eNUMINTEGRAL typeInt = _hermite)
	{
		_x = _w = 0;
		_degree = degree;
		_typeInt = typeInt;
		switch (_typeInt)
		{
			case _hermite:
				_enumDegree = getHermitePOLYNOMIAL(degree);
				break;
			case _legendre :
				_enumDegree = getLegendrePOLYNOMIAL(degree);
				break;
			
			default:
				_enumDegree = _nohermite;
		}
		
		setPolynomialCoeff(degree, _enumDegree);

	}

	double *getx(void) {return _x;}
	double *getw(void) {return _w;}
	int		getDeg(void) {return _degree;}

	ePOLYNOMIALS getHermitePOLYNOMIAL(int degree)
	{
		switch(degree)
		{
		case 12:
			return _hermite12;

		case 20:
			return _hermite20;

		case 40:
			return _hermite40;

		case 80:
			return _hermite80;

		default:
			return _nohermite;
		}
	}
		ePOLYNOMIALS getLegendrePOLYNOMIAL(int degree)
	{
		switch(degree)
		{
		case 8:
			return _legendre8;

		case 12:
			return _legendre12;

		case 20:
			return _legendre20;

		case 40:
			return _legendre40;

		case 80:
			return _legendre80;

		default:
			return _nolegendre;
		}
	}
			double GaussIntegration(double* func, double a = -1000.0, double b = 1000.0)
	{
		switch (_typeInt)
		{
			case _hermite:
				return GaussHermiteInt(func);
			case _legendre :
				return GaussLegendreInt(func, a, b);
			
			default:
				return 0.0;
		}
	}


	double GaussHermiteInt(double* func)
	{
		
		if(_x == 0)
		{
			if(ODD(_degree) == false) return DNOVALUE;

			hermiteCoeff(_x,_w,_degree);
		}

		double integ = 0.;

		for(int k = 0; k < _degree; k++)
		{
			integ += _w[k]*func[k];
		}

		return integ;
	}
	double GaussLegendreInt(double* func, double a, double b)
	{
		if(a >= b) return 0.;
		
		if(_x == 0)
		{
			if(ODD(_degree) == false) return DNOVALUE;

			//legendreCoeff(_x,_w,_degree);
			return 0.0;
		}

		double  scale = b - a;
		double integ = 0.;

		for(int k = 0; k < _degree; k++)
		{
			//x = a + scale*_x[k];
			integ += func[k]*_w[k];
		}
		integ *= scale;

		return integ;
	}



	void setPolynomialCoeff(int& degree, ePOLYNOMIALS polynome)
	{
		if (_x)
			delete _x;
		if (_w)
			delete _w;
		
		int k;

		switch(polynome)
		{
		case _legendre8:
			degree = 8;
			_x = new double [degree];
			_w = new double [degree];
			k = 0;
			_x[k++] = 0.019855071751232;
			_x[k++] = 0.10166676129319;
			_x[k++] = 0.23723379504184;
			_x[k++] = 0.40828267875218;
			_x[k++] = 0.59171732124782;
			_x[k++] = 0.76276620495816;
			_x[k++] = 0.89833323870681;
			_x[k++] = 0.98014492824877;

			k = 0;
			_w[k++] = 0.050614268145188;
			_w[k++] = 0.11119051722669;
			_w[k++] = 0.15685332293894;
			_w[k++] = 0.18134189168918;
			_w[k++] = 0.18134189168918;
			_w[k++] = 0.15685332293894;
			_w[k++] = 0.11119051722669;
			_w[k++] = 0.050614268145188;
			return;

		case _legendre12:
			degree = 12;
			_x = new double [degree];
			_w = new double [degree];
			k = 0;
			_x[k++] = 0.00921968287664038;
			_x[k++] = 0.0479413718147626;
			_x[k++] = 0.115048662902848;
			_x[k++] = 0.206341022856691;
			_x[k++] = 0.31608425050091;
			_x[k++] = 0.437383295744266;
			_x[k++] = 0.562616704255735;
			_x[k++] = 0.68391574949909;
			_x[k++] = 0.793658977143309;
			_x[k++] = 0.884951337097152;
			_x[k++] = 0.952058628185237;
			_x[k++] = 0.99078031712336;
			
			k = 0;
			_w[k++] = 0.0235876681932542;
			_w[k++] = 0.0534696629976591;
			_w[k++] = 0.0800391642716732;
			_w[k++] = 0.101583713361533;
			_w[k++] = 0.116746268269177;
			_w[k++] = 0.124573522906701;
			_w[k++] = 0.124573522906701;
			_w[k++] = 0.116746268269177;
			_w[k++] = 0.101583713361533;
			_w[k++] = 0.0800391642716732;
			_w[k++] = 0.0534696629976591;
			_w[k++] = 0.0235876681932542;
			return;

		case _legendre20:
			degree = 20;
			_x = new double [degree];
			_w = new double [degree];
			k = 0;
			_x[k++] = 0.00343570040745256;
			_x[k++] = 0.0180140363610431;
			_x[k++] = 0.043882785874337;
			_x[k++] = 0.0804415140888906;
			_x[k++] = 0.126834046769925;
			_x[k++] = 0.181973159636742;
			_x[k++] = 0.244566499024586;
			_x[k++] = 0.31314695564229;
			_x[k++] = 0.386107074429177;
			_x[k++] = 0.461736739433251;
			_x[k++] = 0.538263260566749;
			_x[k++] = 0.613892925570823;
			_x[k++] = 0.68685304435771;
			_x[k++] = 0.755433500975414;
			_x[k++] = 0.818026840363258;
			_x[k++] = 0.873165953230075;
			_x[k++] = 0.919558485911109;
			_x[k++] = 0.956117214125663;
			_x[k++] = 0.981985963638957;
			_x[k++] = 0.996564299592547;

			k = 0;
			_w[k++] = 0.0088070035695753;
			_w[k++] = 0.0203007149001935;
			_w[k++] = 0.0313360241670545;
			_w[k++] = 0.0416383707883524;
			_w[k++] = 0.0509650599086202;
			_w[k++] = 0.0590972659807588;
			_w[k++] = 0.0658443192245883;
			_w[k++] = 0.071048054659191;
			_w[k++] = 0.0745864932363018;
			_w[k++] = 0.0763766935653629;
			_w[k++] = 0.0763766935653629;
			_w[k++] = 0.0745864932363018;
			_w[k++] = 0.071048054659191;
			_w[k++] = 0.0658443192245883;
			_w[k++] = 0.0590972659807588;
			_w[k++] = 0.0509650599086202;
			_w[k++] = 0.0416383707883524;
			_w[k++] = 0.0313360241670545;
			_w[k++] = 0.0203007149001935;
			_w[k++] = 0.0088070035695753;
			return;

		case _legendre40:
			degree = 40;
			_x = new double [degree];
			_w = new double [degree];
			k = 0;
			_x[k++] = 0.000881145;
			_x[k++] = 0.004636881;
			_x[k++] = 0.011370025;
			_x[k++] = 0.02104159;
			_x[k++] = 0.033593596;
			_x[k++] = 0.048950597;
			_x[k++] = 0.067020248;
			_x[k++] = 0.087693885;
			_x[k++] = 0.110847174;
			_x[k++] = 0.136340872;
			_x[k++] = 0.164021658;
			_x[k++] = 0.193723055;
			_x[k++] = 0.225266437;
			_x[k++] = 0.258462099;
			_x[k++] = 0.293110398;
			_x[k++] = 0.329002955;
			_x[k++] = 0.365923907;
			_x[k++] = 0.40365121;
			_x[k++] = 0.441957965;
			_x[k++] = 0.480613791;
			_x[k++] = 0.519386209;
			_x[k++] = 0.558042035;
			_x[k++] = 0.59634879;
			_x[k++] = 0.634076093;
			_x[k++] = 0.670997045;
			_x[k++] = 0.706889602;
			_x[k++] = 0.741537901;
			_x[k++] = 0.774733563;
			_x[k++] = 0.806276945;
			_x[k++] = 0.835978342;
			_x[k++] = 0.863659128;
			_x[k++] = 0.889152826;
			_x[k++] = 0.912306115;
			_x[k++] = 0.932979752;
			_x[k++] = 0.951049403;
			_x[k++] = 0.966406404;
			_x[k++] = 0.97895841;
			_x[k++] = 0.988629975;
			_x[k++] = 0.995363119;
			_x[k++] = 0.999118855;

			k = 0;
			_w[k++] = 0.002260639;
			_w[k++] = 0.005249142;
			_w[k++] = 0.008210529;
			_w[k++] = 0.011122925;
			_w[k++] = 0.013968503;
			_w[k++] = 0.016730098;
			_w[k++] = 0.019391084;
			_w[k++] = 0.021935454;
			_w[k++] = 0.024347904;
			_w[k++] = 0.026613923;
			_w[k++] = 0.028719885;
			_w[k++] = 0.030653121;
			_w[k++] = 0.032402007;
			_w[k++] = 0.033956023;
			_w[k++] = 0.035305824;
			_w[k++] = 0.036443291;
			_w[k++] = 0.037361585;
			_w[k++] = 0.038055181;
			_w[k++] = 0.038519909;
			_w[k++] = 0.038752974;
			_w[k++] = 0.038752974;
			_w[k++] = 0.038519909;
			_w[k++] = 0.038055181;
			_w[k++] = 0.037361585;
			_w[k++] = 0.036443291;
			_w[k++] = 0.035305824;
			_w[k++] = 0.033956023;
			_w[k++] = 0.032402007;
			_w[k++] = 0.030653121;
			_w[k++] = 0.028719885;
			_w[k++] = 0.026613923;
			_w[k++] = 0.024347904;
			_w[k++] = 0.021935454;
			_w[k++] = 0.019391084;
			_w[k++] = 0.016730098;
			_w[k++] = 0.013968503;
			_w[k++] = 0.011122925;
			_w[k++] = 0.008210529;
			_w[k++] = 0.005249142;
			_w[k++] = 0.002260639;
			return;

		case _legendre80:
			degree = 80;
			_x = new double [degree];
			_w = new double [degree];
			k = 0;
			_x[k++] = 0.00022308867418469;
			_x[k++] = 0.00117506780088116;
			_x[k++] = 0.00288622951715583;
			_x[k++] = 0.00535434875012225;
			_x[k++] = 0.00857571363068548;
			_x[k++] = 0.0125454297071361;
			_x[k++] = 0.0172574554781003;
			_x[k++] = 0.0227046168281825;
			_x[k++] = 0.0288786193450636;
			_x[k++] = 0.035770061413777;
			_x[k++] = 0.0433684487141211;
			_x[k++] = 0.0516622102806146;
			_x[k++] = 0.0606387161608931;
			_x[k++] = 0.0702842966684444;
			_x[k++] = 0.0805842632098723;
			_x[k++] = 0.0915229306592682;
			_x[k++] = 0.103083641247697;
			_x[k++] = 0.115248789932479;
			_x[k++] = 0.127999851208201;
			_x[k++] = 0.14131740731895;
			_x[k++] = 0.155181177828986;
			_x[k++] = 0.16957005050694;
			_x[k++] = 0.184462113476563;
			_x[k++] = 0.199834688585124;
			_x[k++] = 0.215664365938645;
			_x[k++] = 0.231927039551434;
			_x[k++] = 0.248597944055607;
			_x[k++] = 0.265651692414727;
			_x[k++] = 0.283062314584121;
			_x[k++] = 0.300803297059015;
			_x[k++] = 0.318847623250256;
			_x[k++] = 0.337167814626149;
			_x[k++] = 0.355735972557744;
			_x[k++] = 0.374523820803863;
			_x[k++] = 0.393502748571166;
			_x[k++] = 0.412643854083676;
			_x[k++] = 0.431917988595428;
			_x[k++] = 0.451295800779207;
			_x[k++] = 0.470747781423789;
			_x[k++] = 0.490244308371602;
			_x[k++] = 0.509755691628396;
			_x[k++] = 0.52925221857621;
			_x[k++] = 0.548704199220792;
			_x[k++] = 0.568082011404571;
			_x[k++] = 0.587356145916323;
			_x[k++] = 0.606497251428833;
			_x[k++] = 0.625476179196136;
			_x[k++] = 0.644264027442255;
			_x[k++] = 0.662832185373851;
			_x[k++] = 0.681152376749743;
			_x[k++] = 0.699196702940984;
			_x[k++] = 0.716937685415877;
			_x[k++] = 0.734348307585272;
			_x[k++] = 0.751402055944392;
			_x[k++] = 0.768072960448565;
			_x[k++] = 0.784335634061354;
			_x[k++] = 0.800165311414875;
			_x[k++] = 0.815537886523435;
			_x[k++] = 0.830429949493059;
			_x[k++] = 0.844818822171013;
			_x[k++] = 0.858682592681049;
			_x[k++] = 0.872000148791798;
			_x[k++] = 0.88475121006752;
			_x[k++] = 0.896916358752302;
			_x[k++] = 0.908477069340731;
			_x[k++] = 0.919415736790127;
			_x[k++] = 0.929715703331555;
			_x[k++] = 0.939361283839106;
			_x[k++] = 0.948337789719385;
			_x[k++] = 0.956631551285878;
			_x[k++] = 0.964229938586222;
			_x[k++] = 0.971121380654936;
			_x[k++] = 0.977295383171817;
			_x[k++] = 0.982742544521899;
			_x[k++] = 0.987454570292863;
			_x[k++] = 0.991424286369314;
			_x[k++] = 0.994645651249877;
			_x[k++] = 0.997113770482844;
			_x[k++] = 0.998824932199118;
			_x[k++] = 0.999776911325815;

			k = 0;
			_w[k++] = 0.000572475001593489;
			_w[k++] = 0.00133176679475636;
			_w[k++] = 0.0020901565623475;
			_w[k++] = 0.0028454612256952;
			_w[k++] = 0.00359645238405746;
			_w[k++] = 0.00434197263463013;
			_w[k++] = 0.00508088302055144;
			_w[k++] = 0.00581205706039887;
			_w[k++] = 0.00653438079620066;
			_w[k++] = 0.00724675402025453;
			_w[k++] = 0.00794809179186284;
			_w[k++] = 0.00863732602813464;
			_w[k++] = 0.00931340710414951;
			_w[k++] = 0.009975305439071;
			_w[k++] = 0.0106220130578909;
			_w[k++] = 0.0112525451231662;
			_w[k++] = 0.011865941432965;
			_w[k++] = 0.0124612678820577;
			_w[k++] = 0.0130376178837825;
			_w[k++] = 0.0135941137502431;
			_w[k++] = 0.0141299080286384;
			_w[k++] = 0.0146441847916339;
			_w[k++] = 0.0151361608797789;
			_w[k++] = 0.0156050870940573;
			_w[k++] = 0.0160502493367438;
			_w[k++] = 0.0164709696988226;
			_w[k++] = 0.0168666074923057;
			_w[k++] = 0.0172365602258769;
			_w[k++] = 0.0175802645223737;
			_w[k++] = 0.017897196976708;
			_w[k++] = 0.0181868749529179;
			_w[k++] = 0.018448857319138;
			_w[k++] = 0.0186827451193652;
			_w[k++] = 0.0188881821810007;
			_w[k++] = 0.0190648556572388;
			_w[k++] = 0.0192124965034797;
			_w[k++] = 0.0193308798870382;
			_w[k++] = 0.0194198255295259;
			_w[k++] = 0.0194791979813847;
			_w[k++] = 0.0195089068281533;
			_w[k++] = 0.0195089068281533;
			_w[k++] = 0.0194791979813847;
			_w[k++] = 0.0194198255295259;
			_w[k++] = 0.0193308798870382;
			_w[k++] = 0.0192124965034797;
			_w[k++] = 0.0190648556572388;
			_w[k++] = 0.0188881821810007;
			_w[k++] = 0.0186827451193652;
			_w[k++] = 0.018448857319138;
			_w[k++] = 0.0181868749529179;
			_w[k++] = 0.017897196976708;
			_w[k++] = 0.0175802645223737;
			_w[k++] = 0.0172365602258769;
			_w[k++] = 0.0168666074923057;
			_w[k++] = 0.0164709696988226;
			_w[k++] = 0.0160502493367438;
			_w[k++] = 0.0156050870940573;
			_w[k++] = 0.0151361608797789;
			_w[k++] = 0.0146441847916339;
			_w[k++] = 0.0141299080286384;
			_w[k++] = 0.0135941137502431;
			_w[k++] = 0.0130376178837825;
			_w[k++] = 0.0124612678820577;
			_w[k++] = 0.011865941432965;
			_w[k++] = 0.0112525451231662;
			_w[k++] = 0.0106220130578909;
			_w[k++] = 0.009975305439071;
			_w[k++] = 0.00931340710414951;
			_w[k++] = 0.00863732602813464;
			_w[k++] = 0.00794809179186284;
			_w[k++] = 0.00724675402025453;
			_w[k++] = 0.00653438079620066;
			_w[k++] = 0.00581205706039887;
			_w[k++] = 0.00508088302055144;
			_w[k++] = 0.00434197263463013;
			_w[k++] = 0.00359645238405746;
			_w[k++] = 0.0028454612256952;
			_w[k++] = 0.0020901565623475;
			_w[k++] = 0.00133176679475636;
			_w[k++] = 0.000572475001593489;
			
			return;

		case _hermite12:
			degree = 12;
			_x = new double [degree];
			_w = new double [degree];
			k = 0;
			_x[k++] = 5.50090170446774;
			_x[k++] = 4.27182584793228;
			_x[k++] = 3.22370982877009;
			_x[k++] = 2.25946445100079;
			_x[k++] = 1.34037519715161;
			_x[k++] = 0.444403001944138;
			_x[k++] = -0.444403001944138;
			_x[k++] = -1.34037519715161;
			_x[k++] = -2.25946445100079;
			_x[k++] = -3.22370982877009;
			_x[k++] = -4.27182584793228;
			_x[k++] = -5.50090170446774;

			k = 0;
			_w[k++] = 0.000000149992716763473;
			_w[k++] = 0.0000483718492258277;
			_w[k++] = 0.00220338068752962;
			_w[k++] = 0.0291166879123168;
			_w[k++] = 0.146967048045091;
			_w[k++] = 0.321664361512307;
			_w[k++] = 0.321664361512307;
			_w[k++] = 0.146967048045091;
			_w[k++] = 0.0291166879123168;
			_w[k++] = 0.00220338068752962;
			_w[k++] = 0.0000483718492258277;
			_w[k++] = 0.000000149992716763473;
			return;

		case _hermite20:
			degree = 20;
			_x = new double [degree];
			_w = new double [degree];
			k = 0;
			_x[k++] = 7.61904854167976;
			_x[k++] = 6.51059015701366;
			_x[k++] = 5.5787388058932;
			_x[k++] = 4.73458133404606;
			_x[k++] = 3.94396735065732;
			_x[k++] = 3.18901481655339;
			_x[k++] = 2.45866361117237;
			_x[k++] = 1.74524732081413;
			_x[k++] = 1.04294534880275;
			_x[k++] = 0.346964157081356;
			_x[k++] = -0.346964157081356;
			_x[k++] = -1.04294534880275;
			_x[k++] = -1.74524732081413;
			_x[k++] = -2.45866361117237;
			_x[k++] = -3.18901481655339;
			_x[k++] = -3.94396735065732;
			_x[k++] = -4.73458133404606;
			_x[k++] = -5.5787388058932;
			_x[k++] = -6.51059015701366;
			_x[k++] = -7.61904854167976;
			
			k = 0;
			_w[k++] = 1.25780067243793e-013;
			_w[k++] = 2.48206236231518e-010;
			_w[k++] = 6.12749025998295e-008;
			_w[k++] = 4.40212109023085e-006;
			_w[k++] = 0.000128826279961929;
			_w[k++] = 0.00183010313108049;
			_w[k++] = 0.013997837447101;
			_w[k++] = 0.0615063720639759;
			_w[k++] = 0.161739333984;
			_w[k++] = 0.260793063449555;
			_w[k++] = 0.260793063449555;
			_w[k++] = 0.161739333984;
			_w[k++] = 0.0615063720639759;
			_w[k++] = 0.013997837447101;
			_w[k++] = 0.00183010313108049;
			_w[k++] = 0.000128826279961929;
			_w[k++] = 4.40212109023085e-006;
			_w[k++] = 6.12749025998295e-008;
			_w[k++] = 2.48206236231518e-010;
			_w[k++] = 1.25780067243793e-013;
			return;

		case _hermite40:
			degree = 40;
			_x = new double [degree];
			_w = new double [degree];
			k = 0;
			_x[k++] = 11.453377841549;
			_x[k++] = 10.481560534674;
			_x[k++] = 9.6735563669340;
			_x[k++] = 8.9495045438556;
			_x[k++] = 8.2789406236595;
			_x[k++] = 7.6461637645415;
			_x[k++] = 7.0417384064538;
			_x[k++] = 6.4594233775838;
			_x[k++] = 5.8948056753720;
			_x[k++] = 5.3446054457201;
			_x[k++] = 4.8062871920939;
			_x[k++] = 4.2778261563627;
			_x[k++] = 3.7575597761690;
			_x[k++] = 3.2440887329999;
			_x[k++] = 2.7362083404654;
			_x[k++] = 2.2328592186349;
			_x[k++] = 1.7330905906317;
			_x[k++] = 1.2360320047992;
			_x[k++] = 0.74087072528593;
			_x[k++] = 0.24683289602272;
			_x[k++] = -0.24683289602272;
			_x[k++] = -0.74087072528593;
			_x[k++] = -1.2360320047992;
			_x[k++] = -1.7330905906317;
			_x[k++] = -2.2328592186349;
			_x[k++] = -2.7362083404654;
			_x[k++] = -3.2440887329999;
			_x[k++] = -3.7575597761690;
			_x[k++] = -4.2778261563627;
			_x[k++] = -4.8062871920939;
			_x[k++] = -5.3446054457201;
			_x[k++] = -5.8948056753720;
			_x[k++] = -6.4594233775838;
			_x[k++] = -7.0417384064538;
			_x[k++] = -7.6461637645415;
			_x[k++] = -8.2789406236595;
			_x[k++] = -8.9495045438556;
			_x[k++] = -9.6735563669340;
			_x[k++] = -10.481560534674;
			_x[k++] = -11.453377841549;

			k = 0;
			_w[k++] = 1.4618398738670e-029;
			_w[k++] = 4.8204679401929e-025;
			_w[k++] = 1.4486094315492e-021;
			_w[k++] = 1.1222752068253e-018;
			_w[k++] = 3.3898534432428e-016;
			_w[k++] = 4.9680885291897e-014;
			_w[k++] = 4.0376385816887e-012;
			_w[k++] = 1.9891185260245e-010;
			_w[k++] = 6.3258971885386e-009;
			_w[k++] = 1.3603424215727e-007;
			_w[k++] = 2.0488974360781e-006;
			_w[k++] = 2.2211771432440e-005;
			_w[k++] = 0.00017707292879895;
			_w[k++] = 0.0010558790169001;
			_w[k++] = 0.0047735448818156;
			_w[k++] = 0.016537844142543;
			_w[k++] = 0.044274555202205;
			_w[k++] = 0.092176579169911;
			_w[k++] = 0.14992111176333;
			_w[k++] = 0.19105900966168;
			_w[k++] = 0.19105900966168;
			_w[k++] =  0.14992111176333;
			_w[k++] = 0.092176579169911;
			_w[k++] = 0.044274555202205;
			_w[k++] = 0.016537844142543;
			_w[k++] = 0.0047735448818156;
			_w[k++] = 0.0010558790169001;
			_w[k++] = 0.00017707292879895;
			_w[k++] = 2.2211771432440e-005;
			_w[k++] = 2.0488974360781e-006;
			_w[k++] = 1.3603424215727e-007;
			_w[k++] = 6.3258971885386e-009;
			_w[k++] = 1.9891185260245e-010;
			_w[k++] = 4.0376385816887e-012;
			_w[k++] = 4.9680885291897e-014;
			_w[k++] = 3.3898534432428e-016;
			_w[k++] = 1.1222752068253e-018;
			_w[k++] = 1.4486094315492e-021;
			_w[k++] = 4.8204679401929e-025;
			_w[k++] = 1.4618398738670e-029;
			return;

		case _hermite80:
			degree = 80;
			_x = new double [degree];
			_w = new double [degree];
			k = 0;
			_x[k++] = 16.811977874859;
			_x[k++] = 15.954724894673;
			_x[k++] = 15.246348584550;
			_x[k++] = 14.615177388338;
			_x[k++] = 14.033856524353;
			_x[k++] = 13.488299320645;
			_x[k++] = 12.970060141888;
			_x[k++] = 12.473575934119;
			_x[k++] = 11.994938265320;
			_x[k++] = 11.531268507680;
			_x[k++] = 11.080368463136;
			_x[k++] = 10.640510727646;
			_x[k++] = 10.210305800835;
			_x[k++] = 9.7886140496506;
			_x[k++] = 9.3744852364933;
			_x[k++] = 8.9671157030768;
			_x[k++] = 8.5658172635359;
			_x[k++] = 8.1699940967438;
			_x[k++] = 7.7791252448249;
			_x[k++] = 7.3927511292267;
			_x[k++] = 7.0104630027663;
			_x[k++] = 6.6318945846959;
			_x[k++] = 6.2567153440974;
			_x[k++] = 5.8846250450935;
			_x[k++] = 5.5153492699408;
			_x[k++] = 5.1486357083411;
			_x[k++] = 4.7842510530526;
			_x[k++] = 4.4219783794553;
			_x[k++] = 4.0616149143831;
			_x[k++] = 3.7029701201253;
			_x[k++] = 3.3458640349965;
			_x[k++] = 2.9901258236544;
			_x[k++] = 2.6355924993643;
			_x[k++] = 2.2821077873744;
			_x[k++] = 1.9295211039619;
			_x[k++] = 1.5776866299194;
			_x[k++] = 1.2264624605194;
			_x[k++] = 0.87570981654072;
			_x[k++] = 0.52529230289626;
			_x[k++] = 0.17507520287869;
			_x[k++] = -0.17507520287869;
			_x[k++] = -0.52529230289626;
			_x[k++] = -0.87570981654072;
			_x[k++] = -1.2264624605194;
			_x[k++] = -1.5776866299194;
			_x[k++] = -1.9295211039619;
			_x[k++] = -2.2821077873744;
			_x[k++] = -2.6355924993643;
			_x[k++] = -2.9901258236544;
			_x[k++] = -3.3458640349965;
			_x[k++] = -3.7029701201253;
			_x[k++] = -4.0616149143831;
			_x[k++] = -4.4219783794553;
			_x[k++] = -4.7842510530526;
			_x[k++] = -5.1486357083411;
			_x[k++] = -5.5153492699408;
			_x[k++] = -5.8846250450935;
			_x[k++] = -6.2567153440974;
			_x[k++] = -6.6318945846959;
			_x[k++] = -7.0104630027663;
			_x[k++] = -7.3927511292267;
			_x[k++] = -7.7791252448249;
			_x[k++] = -8.1699940967438;
			_x[k++] = -8.5658172635359;
			_x[k++] = -8.9671157030768;
			_x[k++] = -9.3744852364933;
			_x[k++] = -9.7886140496506;
			_x[k++] = -10.210305800835;
			_x[k++] = -10.640510727646;
			_x[k++] = -11.080368463136;
			_x[k++] = -11.531268507680;
			_x[k++] = -11.994938265320;
			_x[k++] = -12.473575934119;
			_x[k++] = -12.970060141888;
			_x[k++] = -13.488299320645;
			_x[k++] = -14.033856524353;
			_x[k++] = -14.615177388338;
			_x[k++] = -15.246348584550;
			_x[k++] = -15.954724894673;
			_x[k++] = -16.811977874859;

			k = 0;
			_w[k++] = 1.6676172424931e-062;
			_w[k++] = 1.6148147294588e-056;
			_w[k++] = 8.8417336957561e-052;
			_w[k++] = 9.9534648162192e-048;
			_w[k++] = 3.8322199546539e-044;
			_w[k++] = 6.5965480524107e-041;
			_w[k++] = 5.9776715330824e-038;
			_w[k++] = 3.1804454466265e-035;
			_w[k++] = 1.0736474797675e-032;
			_w[k++] = 2.4360675685908e-030;
			_w[k++] = 3.8835747628103e-028;
			_w[k++] = 4.5051594398809e-026;
			_w[k++] = 3.9121581732439e-024;
			_w[k++] = 2.6028723242580e-022;
			_w[k++] = 1.3528365065127e-020;
			_w[k++] = 5.5835185692130e-019;
			_w[k++] = 1.8557109471024e-017;
			_w[k++] = 5.0266251141603e-016;
			_w[k++] = 1.1213385520192e-014;
			_w[k++] = 2.0789538854682e-013;
			_w[k++] = 3.2290148009037e-012;
			_w[k++] = 4.2312638553690e-011;
			_w[k++] = 4.7070652484468e-010;
			_w[k++] = 4.4700307909793e-009;
			_w[k++] = 3.6415407638465e-008;
			_w[k++] = 2.5560809741474e-007;
			_w[k++] = 1.5519301819855e-006;
			_w[k++] = 8.1787428111035e-006;
			_w[k++] = 3.7528431210811e-005;
			_w[k++] = 0.00015034404791796;
			_w[k++] = 0.00052713230670707;
			_w[k++] = 0.0016210281196417;
			_w[k++] = 0.0043803901587737;
			_w[k++] = 0.010418184453696;
			_w[k++] = 0.021838988865634;
			_w[k++] = 0.040396396021903;
			_w[k++] = 0.065999311488051;
			_w[k++] = 0.095313042531020;
			_w[k++] = 0.12173797707836;
			_w[k++] = 0.13756964881328;
			_w[k++] = 0.13756964881328;
			_w[k++] = 0.12173797707836;
			_w[k++] = 0.095313042531020;
			_w[k++] = 0.065999311488051;
			_w[k++] = 0.040396396021903;
			_w[k++] = 0.021838988865634;
			_w[k++] = 0.010418184453696;
			_w[k++] = 0.0043803901587737;
			_w[k++] = 0.0016210281196417;
			_w[k++] = 0.00052713230670707;
			_w[k++] = 0.00015034404791796;
			_w[k++] = 3.7528431210811e-005;
			_w[k++] = 8.1787428111035e-006;
			_w[k++] = 1.5519301819855e-006;
			_w[k++] = 2.5560809741474e-007;
			_w[k++] = 3.6415407638465e-008;
			_w[k++] = 4.4700307909793e-009;
			_w[k++] = 4.7070652484468e-010;
			_w[k++] = 4.2312638553690e-011;
			_w[k++] = 3.2290148009037e-012;
			_w[k++] = 2.0789538854682e-013;
			_w[k++] = 1.1213385520192e-014;
			_w[k++] = 5.0266251141603e-016;
			_w[k++] = 1.8557109471024e-017;
			_w[k++] = 5.5835185692130e-019;
			_w[k++] = 1.3528365065127e-020;
			_w[k++] = 2.6028723242580e-022;
			_w[k++] = 3.9121581732439e-024;
			_w[k++] = 4.5051594398809e-026;
			_w[k++] = 3.8835747628103e-028;
			_w[k++] = 2.4360675685908e-030;
			_w[k++] = 1.0736474797675e-032;
			_w[k++] = 3.1804454466265e-035;
			_w[k++] = 5.9776715330824e-038;
			_w[k++] = 6.5965480524107e-041;
			_w[k++] = 3.8322199546539e-044;
			_w[k++] = 9.9534648162192e-048;
			_w[k++] = 8.8417336957561e-052;
			_w[k++] = 1.6148147294588e-056;
			_w[k++] = 1.6676172424931e-062;
			return;
		}
	}

	
	void hermiteCoeff(double*& x,double*& w,int n)
	{
		if (_x)
			delete  _x;
		if (_w)
			delete  _w;

		_x = new double [n];
		_w = new double [n];

		double sqrt2 = sqrt(2.);
		double sqrtpi = sqrt(PI);

		gauher(x,w,n);

		int k;
		for(k = 0; k < n; _x[k++] *= sqrt2);
		for(k = 0; k < n; _w[k++] /= sqrtpi);
	}

	void gauher(double* x, double* w,int n)
	{
		int error = 0;
		int i,its,j,m;
		double p1,p2,p3,pp,z,z1;

		m=(n+1)/2;

		for(i = 1; i <= m; i++)
		{
			if (i==1)
			{
				z = sqrt((double)(2*n+1))-1.85575*pow((double)(2*n+1),-0.16667);
			}
			else if (i==2)
			{
				z -= 1.14*pow((double)n,0.426)/z;
			}
			else if (i == 3)
			{
				z = 1.86*z-0.86*x[1];
			}
			else if (i == 4)
			{
				z = 1.91*z-0.91*x[2];
			}
			else
			{
				z = 2.0*z-x[i-2];
			}
			for(its = 1; its < 20; its++)
			{
				p1=0.7511255444649425;
				p2= 0.0;
				for(j=1;j<=n;j++)
				{
					p3=p2;
					p2=p1;
					p1=z*sqrt(2.0/j)*p2-sqrt(((double)(j-1))/j)*p3;
				}
				pp=sqrt((double)2*n)*p2;
				z1=z;
				z=z1-p1/pp;
				if (fabs(z-z1) <= 3.0e-14) break;
			}

			if (its > 20) error = 1;

			x[i-1] = z;
			x[n-i] = -z;
			w[i-1] = 2.0/(pp*pp);
			w[n-i] = w[i-1];
		}
	}

};

/*___________________________________________________________________________________________________________________________________________________*/

static const double gaussBound = 7.;

static const double gaussA[5] = {
1.161110663653770e-002,3.951404679838207e-001,2.846603853776254e+001,
1.887426188426510e+002,3.209377589138469e+003
};

static const double gaussB[5] = {
1.767766952966369e-001,8.344316438579620e+000,1.725514762600375e+002,
1.813893686502485e+003,8.044716608901563e+003
};

static const double gaussC[9] = {
2.15311535474403846e-8,5.64188496988670089e-1,8.88314979438837594e00,
6.61191906371416295e01,2.98635138197400131e02,8.81952221241769090e02,
1.71204761263407058e03,2.05107837782607147e03,1.23033935479799725E03
};

static const double gaussD[9] = {
1.00000000000000000e00,1.57449261107098347e01,1.17693950891312499e02,
5.37181101862009858e02,1.62138957456669019e03,3.29079923573345963e03,
4.36261909014324716e03,3.43936767414372164e03,1.23033935480374942e03
};

static const double gaussP[6] = {
1.63153871373020978e-2,3.05326634961232344e-1,3.60344899949804439e-1,
1.25781726111229246e-1,1.60837851487422766e-2,6.58749161529837803e-4
};

static const double gaussQ[6] = {
1.00000000000000000e00,2.56852019228982242e00,1.87295284992346047e00,
5.27905102951428412e-1,6.05183413124413191e-2,2.33520497626869185e-3
};

inline double Gauss(double x, double bound)
{
	double y, z;

	//y = fabs(x);
	y = (x > 0. ? x : -x);

	if(y >= bound) return x > 0. ? 1 : 0.;

	if(y <= 0.46875*SQRT2) 
	{
		// evaluate erf() for |u| <= sqrt(2)*0.46875 
		z = y * y;
		y = x * ((((gaussA[0] * z + gaussA[1]) * z + gaussA[2]) * z + gaussA[3]) * z + gaussA[4])
			/((((gaussB[0] * z + gaussB[1]) * z + gaussB[2]) * z + gaussB[3]) * z + gaussB[4]);
		return 0.5 + y;
	}

	z = exp(-y * y / 2.) / 2.;
	if (y <= 4.0) 
	{
		// evaluate erfc() for sqrt(2)*0.46875 <= |u| <= sqrt(2)*4.0 
		y = y / SQRT2;
		y =
		((((((((gaussC[0] * y + gaussC[1]) * y + gaussC[2]) * y + gaussC[3]) * y + gaussC[4]) * y + gaussC[5]) * y + gaussC[6]) * y + gaussC[7]) * y + gaussC[8])
		/((((((((gaussD[0] * y + gaussD[1])* y + gaussD[2])* y + gaussD[3]) * y + gaussD[4]) * y + gaussD[5]) * y + gaussD[6]) * y + gaussD[7]) * y + gaussD[8]);

		y =  z * y;
	} 
	else
	{
		// evaluate erfc() for |u| > sqrt(2)*4.0 
		z = z * SQRT2 / y;
		y = 2. / (y * y);
		y = y * (((((gaussP[0] * y + gaussP[1]) * y + gaussP[2]) * y + gaussP[3]) * y + gaussP[4]) * y + gaussP[5])
		/(((((gaussQ[0] * y + gaussQ[1]) * y + gaussQ[2]) * y + gaussQ[3]) * y + gaussQ[4]) * y + gaussQ[5]);
		y = z * (INVSQRTPI - y);
	}

	return (x < 0.0 ? y : 1. - y);
}

inline double Gauss(double x)
{
	return Gauss(x,gaussBound);
}

inline double Gauss(double u,double mean,double var, double bound = gaussBound)
{
	return Gauss(mean + sqrt(var) * mean, bound);
}


static const double oG_p = 0.2316419;
static const double oG_b1 = 0.31938153;
static const double oG_b2 = -0.356563782;
static const double oG_b3 = 1.781477937;
static const double oG_b4 = -1.821255978;
static const double oG_b5 = 1.330274429;

inline double oGauss(double x,double bound)
{
	

	double xp,t,Ln,b;
	
	//xp = fabs(x);
	xp = (x > 0. ? x : -x);

	if(xp >= bound) return x > 0. ? 1. : 0.;

	t	= 1 / (1 + oG_p * xp);

	b	= (oG_b2 + t * (oG_b3 + t * (oG_b4 + t * oG_b5)));
	Ln	= 1 - exp(-xp * xp / 2.) * INVSQRT2PI * t * (oG_b1 + t * b);

	if(x < 0.) Ln = 1. - Ln;
	
	return Ln;
}

inline double oGauss(double x)
{
	return oGauss(x,gaussBound);
}

inline double Gauss(double x,double mean,double var)
{
	return oGauss(mean + sqrt(var) * x);
}

inline double GaussDensity(double x,double mean,double var,double bound = gaussBound)
{
	return var < DEPS ? 0. : fabs((x - mean) / sqrt(var)) > bound ? 0. : exp(-0.5 * (x - mean) * (x - mean) / var) * INVSQRT2PI / sqrt(var);
}

inline double GaussDensity(double x, double bound)
{
	return fabs(x) > bound ? 0. : exp(- 0.5 * x * x) * INVSQRT2PI;
}

inline double GaussDensity(double x)
{
	return GaussDensity(x,gaussBound);
}

static const double Abinorm[4] = {
0.3253030,0.4211071,0.1334425,0.006374323
};

static const double Bbinorm[4] = {
0.1337764,0.6243247,1.3425378,2.2626645
};

inline double BiGauss(double x,double y,double rho)
{
	if(rho == 1.) return Gauss(DMIN(x,y));

	// les polynomes

	if(x < 0. && y < 0. && rho < 0)
	{
		double tmp,yal = 0.;
		double rho2 = 1. - rho * rho;
		double a = x / sqrt(2. * rho2);
		double b = y / sqrt(2. * rho2);

		int i,k;
		for(i = 0; i < 4; i++)
		{
			for(k = 0; k < 4; k++)
			{
				tmp = x * (2. * Bbinorm[i] - a) + b * (2. * Bbinorm[k] - b) + 2. * rho * (a - Bbinorm[i]) * (b - Bbinorm[k]);
				yal += Abinorm[i] * Abinorm[k] * exp(yal);
			}
		}

		return yal;
	}
	else if(x * y * rho > 0.)
	{
		int sigx = SIGN_(x);
		int sigy = SIGN_(y);

		double xy = x * x - 2. * rho * x * y + y * y;
		double r1 = sigx * (rho * x - y) / sqrt(xy);
		double r2 = sigy * (rho * y - x) / sqrt(xy);
		double dd = .25 * (1. - sigx * sigy);

		return BiGauss(x,0.,r1) + BiGauss(y,0.,r2) + dd;
	}
	else if(x <= 0.)
	{
		return Gauss(x) - BiGauss(x,-y,-rho);
	}
	else if(y <= 0.)
	{
		return Gauss(y) - BiGauss(-x,y,-rho);
	}
	else if(rho <= 0.)
	{
		return Gauss(x) + Gauss(y) - 1. + BiGauss(-x,-y,rho);
	}

	return 0.;
}

inline double BiGaussDensity(double x, double y, double rho)
{
	// check le cas dégénéré
	if(rho == 1.) return GaussDensity(DMIN(x,y));

	double rho2 = 1. - rho * rho;
	double xy = - (x * x - 2. * rho * x * y + y * y) / (2. * rho2);

	return exp(xy)/(2. * PI * sqrt(rho2));
}

	// les polynomes : 
static const double Acndev[4] = {
2.50662823884,-18.61500062529,41.39119773534,-25.44106049637
};

static const double Bcndev[4] = {
-8.47351093090,23.08336743743,-21.06224101826,3.13082909833
};

static double Ccndev[9] = {
0.3374754822726147,0.9761690190917186,0.1607979714918209,
0.0276438810333863,0.0038405729373609,0.0003951896411919,
0.0000321767881768,0.0000002888167364,0.0000003960315187
};

inline double cndev(double u)
{
	// petit rappel : 
	// si x suit une loi F() alors F(x) suit une loi uniforme et si u est uniforme alors F^(-1)(u) suit une loi F()

	// u est uniform sur [0,1] : on renvoie la gaussienne

	// les cas dégénérés : 
	if(u > 1 - DEPS) return gaussBound;
	if(u < DEPS) return -gaussBound;

	double x,r;
	x = u - 0.5;
	if(fabs(x) < 0.42) 
	{
		r = x*x;
		r = x*(((Acndev[3]*r + Acndev[2])*r + Acndev[1])*r + Acndev[0])/
			((((Bcndev[3]*r + Bcndev[2])*r + Bcndev[1])*r + Bcndev[0])*r + 1.0);
		return r;
	}

	r = u;
	if(x > 0.) r = 1. - u;
	r = log(-log(r));
	r = Ccndev[0] + r*(Ccndev[1] + r*(Ccndev[2] + r*(Ccndev[3]+r*(Ccndev[4] + r*(Ccndev[5] + r*(Ccndev[6] + r*(Ccndev[7] + r*Ccndev[8])))))));
	if(x < 0.) r = -r;

	if(r < -gaussBound) r = -gaussBound;
	if(r > gaussBound) r = gaussBound;
	
	return r;	
}

static const double norminvA[6] = {
-3.969683028665376e+01,  2.209460984245205e+02,
-2.759285104469687e+02,  1.383577518672690e+02,
-3.066479806614716e+01,  2.506628277459239e+00
};

static const double norminvB[5] = {
-5.447609879822406e+01,  1.615858368580409e+02,
-1.556989798598866e+02,  6.680131188771972e+01,
-1.328068155288572e+01
};

static const double norminvC[6] = {
-7.784894002430293e-03, -3.223964580411365e-01,
-2.400758277161838e+00, -2.549732539343734e+00,
4.374664141464968e+00,  2.938163982698783e+00
};

static const double norminvD[4] = {
7.784695709041462e-03,  3.224671290700398e-01,
2.445134137142996e+00,  3.754408661907416e+00
};

inline double NormInv(double p)
{
	double q, t, u;

	if(p > 1 - DEPS) return gaussBound;
	if(p < DEPS) return -gaussBound;

	q = DMIN(p,1. - p);

	if(q > 0.02425) 
	{
		// Rational approximation for central region. 
		u = q - 0.5;
		t = u * u;
		u = u * (((((norminvA[0] * t+ norminvA[1]) * t + norminvA[2]) * t + norminvA[3]) * t + norminvA[4]) * t + norminvA[5])
		/(((((norminvB[0] * t + norminvB[1]) * t + norminvB[2]) * t + norminvB[3]) * t + norminvB[4]) * t + 1.);
	}
	else
	{
		// Rational approximation for tail region.
		t = sqrt(-2. * log(q));
		u = (((((norminvC[0] * t + norminvC[1]) * t + norminvC[2]) * t + norminvC[3]) * t + norminvC[4]) * t + norminvC[5])
		/((((norminvD[0] * t + norminvD[1]) * t + norminvD[2]) * t + norminvD[3]) * t + 1.);
	}

	// The relative error of the approximation has absolute value less
	// than 1.15e-9.  One iteration of Halley's rational method (third
	// order) gives full machine precision... 

	//t = Gauss(u) - q;					// error 
	//t = t * SQRT2PI * exp(u * u / 2.);	// f(u)/df(u) 
	//u = u - t / (1. + u * t / 2.);		// Halley's method 

	return (p > 0.5 ? - u : u);
}


inline double NormInvNag(double p)
{
  NagError fail; 

  INIT_FAIL(fail) ;

  double ret = g01fac(Nag_LowerTail, p,&fail);

  return ret ;
}



inline double probaSachV( double &p, double &v, double &beta, double &beta2)
{
	double pv = (p - beta*v)/sqrt(1-beta2);

	return oGauss(pv,gaussBound);
	//return Gauss(pv,gaussBound);
}



class DiffStochIntensity  
{
protected :
	double	itsMRS; // meanreversion
	double	itssigma; // vol de lambda
	//double  itsIR;
	ARM_ZeroCurve* itsZeroCurve ; // 0-1 : agregation
	int		itssizeProbaInit;  // size vect proba initiales
	ARM_Vector* itsdateProba;
	ARM_Vector* itsprobaInitPF; // P(to > t) en 0 moyenne du PF

	double	itsrecovery;	
	
	double	itsbaseCDS; //daycount

	int		itsDegreeHerm;
	GaussHermite	itsGaussHermite;
	
	int		itsDegreeLegendre;
	GaussHermite	itsGaussLegendre;

	double	itsdt;
	double	itsDealMatu;
	double	itsDiffmatu;
	double	itsdtReset;

	double	itsBetas[4];
	double	itsStrikes[4];

	double	***itsCorrEPV; // ARM_Matrix
	double	*itsCorrET[4];// ARM_Matrix

	

public :
	DiffStochIntensity(double mrs, double sigma, double recovery, ARM_ZeroCurve* zc, ARM_Vector* dateProba,
						ARM_Vector* probaInit, int degreeHermite , double dt , double dealMatu, double diffmatu,
						double dtReset, double beta0 , double beta1, double beta2, double beta3 , double K0, double K1, double K2, double K3,
						int degreeLegendre)
		: itsMRS(mrs), itssigma(sigma), itsrecovery(recovery),  itsdateProba(dateProba), itsprobaInitPF(probaInit) ,
		itsDegreeHerm(degreeHermite), itsDegreeLegendre(degreeLegendre),itsdt(dt), itsDealMatu(dealMatu), itsDiffmatu(diffmatu),
		itsdtReset(dtReset),itsGaussHermite(degreeHermite,_hermite),itsGaussLegendre(degreeLegendre,_legendre)
	{
		itsbaseCDS = 0.25;
		itssizeProbaInit	= probaInit ->GetSize();
		itsCorrEPV = NULL;
		itsBetas[0] = beta0;
		itsBetas[1] = beta1;
		itsBetas[2] = beta2;
		itsBetas[3] = beta3;

		itsStrikes[0] = K0;
		itsStrikes[1] = K1;
		itsStrikes[2] = K2;
		itsStrikes[3] = K3;
		itsZeroCurve = NULL;


		itsZeroCurve= (ARM_ZeroCurve*) zc->Clone();


	}
	
	double getDt(void) { return itsdt;}
	double getDtReset(void){ return itsdtReset;}

	double ProbaSurvPF(double &t, double &T, double &xt);
	double ProbaSurvInitPF(double &t);
	double MgE(double &t, double &T, double &xt);

	// spread CDS en t sachant xt
	double SpreadPF(double &t, double &T, double &xt, double &duration);
	double GetCorrET(int i, int j){ return itsCorrET[i][j] ;} 
	
	~DiffStochIntensity(void)
	{
		/*if(itsdateProba)
			delete[] itsdateProba;
		itsdateProba = NULL;

		if ( itsprobaInitPF)
			delete[] itsprobaInitPF;
		itsprobaInitPF = NULL;*/
		if (itsZeroCurve)	
			delete itsZeroCurve;
		itsZeroCurve = NULL;

	}

	double lambdaY(double &t);

	void CptProbaLossHomog(double ** &probaloss_v, double * &proba_pf_v, double &prob_noms, double* &vsyst, 
				  int &degreeHerm, double &beta, int &nbNoms, int nbDefmax, double t, double T, int adj = 0);


	void cptCorrEpv(double sigma);	

};



class DiffStochIntensityHetero : public DiffStochIntensity
{
protected :
	
	int		itsnbNoms; // nb de nom du PF	
	ARM_Vector* itsprobaInitNom; 

	double	itsLgd;
	double ** itsprobloss_v;

	

public :
	DiffStochIntensityHetero(double mrs, double sigma, double recovery, ARM_ZeroCurve* zc,int nbNoms, ARM_Vector* dateProba,
						ARM_Vector* probaInit, int degreeHermite , double dt , double dealMatu, double diffmatu,
						double dtReset, double beta0 , double beta1, double beta2, double beta3 , double K0, double K1, double K2, double K3,
						int degreeLegendre)
		: DiffStochIntensity(mrs, sigma, recovery , zc, dateProba, probaInit,degreeHermite,dt,dealMatu,
			diffmatu,dtReset,beta0,beta1,beta2,beta3,K0,K1,K2,K3,degreeLegendre)
	{
				
		itsnbNoms = nbNoms;
		itsLgd = 1.0/itsnbNoms * ( 1.0-itsrecovery);

		if (itsnbNoms > 0)
		{
			itsprobaInitNom = new ARM_Vector(itsnbNoms);
			for (int i = 0; i < probaInit->GetSize()-1; i++)
				itsprobaInitNom->Elt(i) = probaInit->Elt(i+1);
		}
		
		itsprobloss_v = new double*[itsnbNoms+1];
		for ( int k=0;k<=itsnbNoms;k++)
		{
			itsprobloss_v[k] = new double[itsDegreeHerm];
		}
		cptCorrET();
	}
	
	double ProbaSurvNom(int &nom, double &t, double &T, double &xt);

	// proba de survi vu en 0 
	// P(to > t)
	double ProbaSurvInitNom(int &nom, double &t);

	
	void PrixTranche(double &fl, double &dl, double &K1, double &K2, double &t, double &T, 
				 double &xt, double &beta, int adj);

	void PrixTranche0(double &fl, double &dl, double &K1, double &K2, double &t, double &T, 
				 double &beta, int adj);

	void cptCorrET(void);

	
	~DiffStochIntensityHetero(void)
	{
		if ( itsprobaInitNom)
			delete[] itsprobaInitNom;
		itsprobaInitNom = NULL;

	}
	
};


