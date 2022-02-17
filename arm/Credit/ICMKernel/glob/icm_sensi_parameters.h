
#ifndef _ICM_SENSI_PARAMETERS_H_
#define _ICM_SENSI_PARAMETERS_H_

#include "ICMKernel/glob/icm_enums.h"

// ***********************************************************************************************
// *							DEFINITION DES SENSIS PARAMETERS											 *	
// ***********************************************************************************************

class ICM_BASIS_BUMP
{
	public :
	string			  itsLabel;
	vector<string>    itsTenors;
	double			  itsEpsilon;

	void Reset(void)
	{
		itsTenors.clear();
		itsLabel = "NONE";
		itsEpsilon = CREDIT_DEFAULT_VALUE;
	}

	ICM_BASIS_BUMP(const vector<string>& Tenors = vector<string>(),
				   const string& label = "NONE",
				   const double& eps = CREDIT_DEFAULT_VALUE ) 
	{
		Reset();
		Set(Tenors,label,eps);
	}

	inline void Set(const vector<string>& Tenors = vector<string>(),
					const string& label = "NONE",
					const double& eps = CREDIT_DEFAULT_VALUE ) 	
	{
		itsLabel = label;
		itsTenors = Tenors;
		itsEpsilon = eps ;

	}

	inline void GetTenors(ARM_CRV_TERMS tab,ARM_Vector& vector)
	{
		memset(tab,'\0',sizeof(char)*ARM_NB_TERMS*ARM_NB_MAX_CHAR_TERMS);
		vector.Resize(itsTenors.size());

		for (int i=0;i<itsTenors.size();i++)
		{
			strcpy(tab[i],itsTenors[i].c_str());
			vector.Elt(i)=itsEpsilon;
		}
	}

	inline bool IsParallelShift(void)
	{
		if (itsTenors.size()==0)
			return true;
		return false;
	} 
};

template <class T> class ICM_BUMP
{
	public :
	vector<T>		itsBumps;
	qSENSITIVITY_TYPE itsSensiType;
	double			itsEpsilon;

	ICM_BUMP(void)
	{
		Reset();
	}

	ICM_BUMP(const qSENSITIVITY_TYPE& SensiType,
			const double& eps = CREDIT_DEFAULT_VALUE) 
	{
		Reset();
		Set(SensiType,eps);
	}
	
	inline void Set(const qSENSITIVITY_TYPE& SensiType = ICMSPREAD_TYPE,
					const double& eps = CREDIT_DEFAULT_VALUE) 	
	{
		itsSensiType = SensiType;
		itsEpsilon = FindDefaultEps(SensiType,eps);
	}

	void Reset(void)
	{
		itsSensiType = ICMSPREAD_TYPE;
		itsBumps.clear();
		itsEpsilon = CREDIT_DEFAULT_VALUE;
	}

	inline void Push(T& bump)
	{
		int sz=itsBumps.size();
		itsBumps.resize(sz+1);
		//bump.itsSensiType = itsSensiType;
		itsBumps[sz] = bump;
	}

	inline double FindDefaultEps(const qSENSITIVITY_TYPE& SensiType,const double& eps)
	{
		if (eps != CREDIT_DEFAULT_VALUE)
			return eps;

		switch (SensiType)
		{			
		case ICMSPREAD_TYPE:
			return 0.01; //Bump sur le spread fixé à +1 bp.
		case ICMIRCURVE_TYPE:
			return 0.01;	//Bump sur les taux fixé à +1bp	
		case ICMRECOVERY_TYPE:
			return 0.1;  //Bump fixe sur le recovery de +0.1
		case ICMCORRELATION_TYPE:
			return 0.1;  //Bump fixe sur la correl de +0.1
		case ICMBETA_TYPE:
			return 0.1;  //Bump fixe sur la correl de +0.1
		case ICMSPRELSHIFT_TYPE:
			return 0.1 ;//Bump relatif des spreads de +0.1
		case ICMCORREL_STRIKE_DOWN_TYPE:
			return -1. ;//Bump relatif des spreads de -1
		case ICMCORREL_STRIKE_UP_TYPE:
			return 1. ;//Bump relatif des spreads de +1
		case ICM_GREEK_VEGA_TYPE:
			return 0.01 ;//Bump sur la volatilité fixé à 0.01
		case ICMRECOVERY_BAR_TYPE:
		case ICM_RECOVERY_DEFAULT_TYPE:
		case ICM_ISSUER_DEFAULT:
// 17783 		case ICM_FAST_SPREAD_TYPE:
		case ICM_SAMECORRELATION_TYPE:
		case ICM_SAMEBETA_TYPE:
		case ICMIRCURVE_WITHOUT_DEFCURVE_TYPE:
		case ICM_IRCURVE_WITH_CPN:
		case ICMBETA_WITH_SPREAD_SHIFT:
		case ICM_THETA_Type:
// 17783 		case ICM_FAST_RECOVERY_TYPE:
		case ICM_DTR_TYPE:
			return CREDIT_DEFAULT_VALUE ;
		default :
			return CREDIT_DEFAULT_VALUE ;
		;
		}
	}

};

#endif  //_ICM_ENUMS_H_
