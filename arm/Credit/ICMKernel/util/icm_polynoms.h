
#if !defined(ICM_POLYNOM_CLASS)
#define ICM_POLYNOM_CLASS

#include <string>
#include <vector>
#include "ARMKernel\glob\linalg.h"
#include "ICMKernel/util/icm_macro.h"


/*********************************************************************************/
/*! \class  ICM_Poly icm_polynoms.h "icm_polynoms.h"
 *  \author D. Pouponneau
 *	\version 1.0
 *	\date   June 2004
 *	\file   icm_polynoms.h
 *	\brief for polynoms computing 
**********************************************************************************/

#define POLY_DEFORDER	10
class ICM_Poly
{
private:
	double* itsPoly;					// aggregated or not  , might point to the internal buffer
	int		itsOrder;					// -1 means null. 
private:
	//	internal preallocated buffers
	double itsPoly_[POLY_DEFORDER+1]; 
public:
	ICM_Poly()
	{
		Init(); 
	}
	ICM_Poly(const ICM_Poly& poly) 
	{
		Init(); 
		if (poly.itsOrder>POLY_DEFORDER) 
		{
			itsPoly = new double[poly.itsOrder+1];
		}
		else 
		{
			itsPoly=itsPoly_; 
		}
		itsOrder=poly.itsOrder ; 
		if (itsOrder!=-1)		memcpy(itsPoly,poly.itsPoly,sizeof(double)*(itsOrder+1));
	}	
private:
	void Init()
	{
		itsPoly = 0;
		itsOrder = -1;
	}
public:
	~ICM_Poly()
	{
		if (itsOrder>POLY_DEFORDER) 
		{
			if (itsPoly)			delete[] itsPoly;
		}
		itsPoly=0; 
	}
	double& operator()(int pos) 
	{
		if (pos>itsOrder || pos<=-1 ) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Poly::() "<<pos) ;
		return itsPoly[pos]; 
	}
	int order() const { return itsOrder; }
	void ResetOrder(int order)
	{
		if (itsOrder>POLY_DEFORDER) 
		{
			if (itsPoly) delete [] itsPoly; 
		}
		if (order==-1) { Init(); return; }
		if (order>POLY_DEFORDER) itsPoly=new double[order+1]; 
		else itsPoly=itsPoly_; 
		itsOrder=order; 
		memset(itsPoly,'\0',(itsOrder+1)*sizeof(double)); 
	}
	// ------------------------------
	//	Constructeur Polynome constant
	// ------------------------------
	inline ICM_Poly& operator= (const ICM_Poly& poly) 
	{
		if (this!=&poly)
		{
			this->~ICM_Poly(); 
			new(this)ICM_Poly(poly); 
		}
		return *this; 
	}

	double Evaluate(double x) const 
	{
		double value =0,x2=1.;
		if (itsPoly == NULL) return -999.;
		for (int i =0; i<=itsOrder; i++)
		{
			value += itsPoly[i]*x2;
			x2 *= x;
		}
		return (value);
	}
	void Integrate(void)
	{
		if (itsOrder>POLY_DEFORDER-1) 
		{
			// un peu à la hussarde
			ICM_Poly temp; temp.ResetOrder(itsOrder+1) ; 
			for (int i =1; i<itsOrder+2; i++)
				temp.itsPoly[i]= itsPoly[i-1]/((double)i);
			temp.itsPoly[0]=0; 
			(*this)=temp; 
		}
		else 
		{
			// no reallocation needed within the poly buffer. 
			for (int i =1; i<itsOrder+2; i++)
				itsPoly[i]= itsPoly[i-1]/((double)i);
			itsPoly[0]=0; 
			//itsDerivatePoly=0; 
			itsOrder++; 
		}
	}		

	double EvaluateDerivatives(double x)
	{
		if (itsOrder==-1) return -999;		// should throw here.. 
		if (itsOrder==0) return 0 ;			// derivative of constant is 0
		double value=0,x2=1; 
		for(int i=1;i<=itsOrder;i++) 
		{
			value += i*itsPoly[i]*x2; 
			x2*=x; 
		}
		return (value); 
	}
	double EvaluateIntegraleExp(double a,double t1,double t2) const 
	{
		double value =0.,t1_2 =t1,t2_2=t2,ui=0.,uiprec=0.;
		double ex1 = exp(a*t1);
		double ex2 = exp(a*t2);

		int i = 0;

		if (itsOrder==0) return (itsPoly[0]/a)*(ex2-ex1);

		uiprec = (ex2-ex1)*(1/a);
		value = uiprec *itsPoly[0];

		for (i =1; i<=itsOrder; i++)
		{
			ui = (t2_2*ex2-t1_2*ex1)*(1/a) - (((double)i)/a) * uiprec;
			uiprec = ui;
			t1_2 *= t1;t2_2 *= t2;
			value += ui*itsPoly[i];
		}
 
		return (value);
	}

	void somme(const ICM_Poly& p1,const ICM_Poly& p2) 
	{
		int order1 = p1.itsOrder;
		int order2 = p2.itsOrder;
		int Order = order1;
		int i=0;

		if (order1 > order2) Order = order1;
		else 
		if (order2 > order1) Order = order2;

		ICM_Poly temp ; // required : this often is p1 or p2
		temp.ResetOrder(Order); 

		for	(i=0; i<=order1; i++)
			temp.itsPoly[i] += p1.itsPoly[i];

		for	(i=0; i<=order2; i++)
			temp.itsPoly[i] += p2.itsPoly[i];
		*this=temp; 
	}

	void multiplication (const ICM_Poly& p1,const ICM_Poly& p2) 
	{
		int order1 = p1.itsOrder;
		int order2 = p2.itsOrder;

		ICM_Poly temp ;		// required : this often is p1 or p2
		temp.ResetOrder(order1+order2); 
		for (int i=0; i<=order1;i++)
			for (int j=0; j<=order2;j++)
			{
				temp.itsPoly[i+j] += p1.itsPoly[i]*p2.itsPoly[j];
			}
		*this=temp; 

	}
	bool operator==(const ICM_Poly&ref) const
	{
		if (itsOrder!=ref.itsOrder) return false; 
		if (itsOrder==-1) return true; 
		for(int i=0;i<=itsOrder;i++) 
			if (itsPoly[i]!=ref.itsPoly[i]) return false; 
		return true; 
	}
	// ---------------------------------------------------------------
	// Full Taylor Expansion for expression Exp(a0+a1*X+a2*X²) on [t1,t2] returns alpha & Px 
	// for Exp(alpha*X)*Px
	// ---------------------------------------------------------------
	static inline void TaylorExpansion(const double& a0,
								const double& a1,
								const double& a2,
								const double& t1,
								const double& t2,
								double& alpha,
								ICM_Poly& Poly_x)
	{
		double Px = 0.;
		// double* coefs = NULL;
		int degree = 1;
		double epsilon = 1.E-6;
		int maxloop = 10;
		int deg = 0;
		double tmp = 0.;
		double tc = (t1+t2)/2.;

		//definition du Polynome Q d'ordre 2
		
		ICM_Poly Q ; 
		Q.ResetOrder(2); 
		Q(0) = a0;
		Q(1) = a1;
		Q(2) = a2;

		alpha = Q.EvaluateDerivatives(tc);

		//definition du Polynome P d'ordre 1
		double	Qtc	=	Q.Evaluate(tc);
		tmp = Qtc-tc*alpha;
		Px = exp(tmp);
		ICM_Poly P ; 
		P.ResetOrder(0); 
		P(0)= Px; 
			
		ICM_Poly Q2 ; 
		Q2.ResetOrder(2); 
		Q2(0)=Q(0) - tmp;
		Q2(1)=Q(1) - alpha;
		Q2(2)=a2;

		ICM_Poly R ; 
		R.ResetOrder(0) ;
		R(0)=1; 

		// JLA ICM_Poly Qd;
		// JLA Qd = R;
		ICM_Poly Qd(R); 
		deg = 1;

 		double exp_q2_t1 = exp(Q2.Evaluate(t1)); 
		double exp_q2_t2 = exp(Q2.Evaluate(t2)); 
		while (((fabs((exp_q2_t1  -  R.Evaluate(t1))/exp_q2_t1)>epsilon) ||
			 (fabs((exp_q2_t2 -  R.Evaluate(t2))/exp_q2_t2)>epsilon)) &&
			 degree <maxloop)

		{
			Qd.multiplication(Qd,Q2);
			deg *= degree;

			ICM_Poly Qd2(Qd); 

			for (int i =0; i<=Qd.order(); i++)
					// Qd2.itsPoly[i] /= (double)deg;
					Qd2(i) /= (double)deg;

			R.somme(R,Qd2);

			degree++;
		}

		
		Poly_x.multiplication(R,P); 

		return;
	}
};


#endif // !defined(ICM_POLYNOM_CLASS)




