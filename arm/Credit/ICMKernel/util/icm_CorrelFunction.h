#ifndef _ICM_CorrelFunction_H_
#define _ICM_CorrelFunction_H_

#include <vector>
#include "ARMKernel/glob/linalg.h"
#include "ICMKernel/glob/icm_enums.h"



// 17783 using namespace std;
/*********************************************************************************/
/*! \class  ICM_CorrelFunction
 *  \author CN
 *	\version 1.0
 *	\date   June 2005
 *	\file   ICM_CorrelFunction.h
 *		Class describing the Correlation function in Random Factors Model
/*****************************************************************************************/

class ICM_CorrelFunction : public ARM_Object
{

private :
	
	int						itsFactorSize;			//Factor Size (default = 1)
	qCopula_TYPE			itsCopulaType;			//Copula type
	qCAL_INDEX_CORR_TYPE	itsCorrelFunctionType;	//Function Type (PWC, PWL)
	int						itsNbLevels;			//Nb of Different level (linear function)
	std::vector<double>		itsLevels;				//Correlation Levels 
	std::vector<double>		itsThreshold;			//Correlation Threshold (nb threshold = nb levels - 1)
	std::vector<double>		itsLinCoef;				//Linear Coef (nb = nb level)

	//Va Intermediare
	double itsM;
	double itsNU;
	double itsVolatility;
	double itsMarginalProb;

public : 

	//Constructor

	ICM_CorrelFunction() {Init();}
	
	ICM_CorrelFunction(int FactorSize , qCopula_TYPE copule, qCAL_INDEX_CORR_TYPE type, vector<double> &param,double Volatility);

	//Destructor
	~ICM_CorrelFunction()
	{
		itsLevels.clear();
		itsThreshold.clear();
		itsLinCoef.clear();
	}	

	//Init
	void Init(void)
	{
		itsFactorSize = 1;
		itsCopulaType = qGAUSSIAN;
		itsCorrelFunctionType = qCAL_PWC_CORREL;

		itsNbLevels = 1;
		itsLevels.clear();
		itsThreshold.clear();
		itsLinCoef.clear();
		itsM=0.;
		itsNU=0.;
		itsMarginalProb=0.;
		itsVolatility=0.;
	}
	
	void InitVAInter(void);
	void Set(int FactorSize, qCopula_TYPE copule, qCAL_INDEX_CORR_TYPE type, vector<double> &param,double Volatility);

	inline void SetFactorSize(int size) { itsFactorSize = size;}
	inline void SetNbLevels(int size) { itsNbLevels = size;}
	inline void SetCorrelFunctionType(qCAL_INDEX_CORR_TYPE type) {itsCorrelFunctionType = type;}
	inline void SetCopulaType(qCopula_TYPE copule) {itsCopulaType = copule;}
	inline void SetMarginalProb(double prob) { itsMarginalProb = prob;}

	//Set Parameters
	inline void SetLevel(int index, double level) 
	{ 
		if (index < 0 || index >= itsNbLevels)
	        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "SetLevel(CorrelFunction) : Incompatible arguments");
		else
			itsLevels[index] = level;
	}
	
	inline void SetThreshold(int index, double thres) 
	{ 
		if (index < 0 || index >= itsNbLevels-1)
	        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "SetLevel(CorrelFunction) : Incompatible arguments");
		else
			itsThreshold[index] = thres;
	}

	inline int					GetFactorSize(void) {return itsFactorSize;}
	inline int					GetNbLevels(void) {return itsNbLevels;}
	inline qCAL_INDEX_CORR_TYPE	GetCorrelFunctionType(void) {return itsCorrelFunctionType;}
	inline qCopula_TYPE			GetCopulaType(void) {return itsCopulaType;}
	inline double				GetM(void) {return itsM;}
	inline double				GetNU(void) {return itsNU;}
	inline double				GetMarginalProb(void) {return itsMarginalProb;}
	
	inline double GetLevel(int index) 
	{ 
		double res = 0.;
		if (index < 0 || index >= itsNbLevels)
	        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "GetLevel(CorrelFunction) : Incompatible arguments");
		else
			res = itsLevels[index];
		return res;
	}
	
	inline double GetThreshold(int index) 
	{ 
		double res = 0.;
		if (index < 0 || index >= itsNbLevels-1)
	        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "GetLevel(CorrelFunction) : Incompatible arguments");
		else
			res = itsThreshold[index];
		return res;
	}
	
	// ----------------------------
	//	Copy of members data
	// ----------------------------

	void BitwiseCopy(const ARM_Object* src)
	{
		ARM_Object::BitwiseCopy(src);

		ICM_CorrelFunction * CorrelFunction = (ICM_CorrelFunction *) src;
		SetFactorSize(CorrelFunction->GetFactorSize());
		SetCorrelFunctionType(CorrelFunction->GetCorrelFunctionType());
		SetCopulaType(CorrelFunction->GetCopulaType());
		
		itsNbLevels = CorrelFunction->itsNbLevels;
		itsLevels = CorrelFunction->itsLevels;
		itsThreshold = CorrelFunction->itsThreshold;
		itsLinCoef = CorrelFunction->itsLinCoef;	
		itsM=CorrelFunction->itsM;
		itsNU=CorrelFunction->itsNU;
		itsMarginalProb=CorrelFunction->itsMarginalProb;
	}


	// -------------
	//	Copy Method 
	// -------------
	
	void Copy(const ARM_Object* src)
	{
		ARM_Object::Copy(src) ;
		BitwiseCopy(src) ;
	}

	// --------------
	//	Clone Method
	// --------------

	ARM_Object * Clone(void)
	{
		ICM_CorrelFunction * theClone = new ICM_CorrelFunction() ;
		theClone->Copy(this) ;
		return(theClone) ;
	}

	// Return Beta coef given the factor
	double BetaCond(const double &factor);

	// Return the Marginal Default Probabilitie
	double DiffMarginalDefprob(const double &seuil);

	// Find the Normal threshold given marginal prob
	double EstimeSeuil(const double &MarginalDefprob);

};


#endif 