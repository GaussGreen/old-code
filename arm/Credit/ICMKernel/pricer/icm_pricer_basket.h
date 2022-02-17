
#ifndef _ICM_Pricer_Basket_H
#define _ICM_Pricer_Basket_H

#include "ICMKernel\pricer\icm_pricer_security.h"


/*********************************************************************************/
/*! \class  ICM_Pricer_Basket icm_pricer_basket.h "icm_pricer_basket.h"
 *  \author Damien Pouponneau 
 *	\version 1.0
 *	\date   March 2004
 *	\file   ICM_Pricer_Basket.h
 *	\brief  Pricer for Basket Security */
/***********************************************************************************/

class ICM_Pricer_Basket : public ICM_Pricer_Security
{
private:

	ARM_Vector * itsBetas;

	double			itsEqCorrelUp;		//TRANCHES -correl results
	bool			itsEqCorrelUpFlg;
	double			itsEqCorrelDown;
	bool			itsEqCorrelDownFlg;
	
public :
	ICM_Pricer_Basket(void) { Init();}

	virtual ~ICM_Pricer_Basket() ;

	void Init() ;

	void Set(ARM_Security *sec, ARM_Model *mod,const ICM_Parameters&params,const ARM_Date&asof);

	virtual void ComputeBeta();

 
	

	ARM_Vector * GetBetas()
	{
		if (!itsBetas)
			ComputeBeta() ;

		return itsBetas ;
	}


	//JLA unused 	ARM_Vector * GetBetasPointer() { return itsBetas ; }		

	void SetBetasPointer(ARM_Vector * betas) 
	{ 
		if (itsBetas)
			delete itsBetas;
		itsBetas = betas; 
	}
	
	double GetBeta(const std::string& issuer) ;


	
protected:
	virtual double ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
		const std::string& plot,
		const std::string& label,
		double epsilon, double epsilonGamma =0 ); 
	public:

/**
	virtual double ComputeCorrelSmile(int typesmile, double MktData, double seed, double UpfrontPay, int DataType);
**/ 
	virtual void BeforePrice(const std::string& label,qSENSITIVITY_TYPE  ) {}
	virtual void Reset(void) ;
	virtual void ResetLU(void) {} ;

	
	//Correl Up
	virtual void ComputeEqCorrelUp();
	inline virtual void setEqCorrelUp(const double & CorrelUp) { itsEqCorrelUp = CorrelUp; itsEqCorrelUpFlg = true;}
	inline virtual double getEqCorrelUp(){ 
										if (itsEqCorrelUpFlg) return itsEqCorrelUp; 
										else ICMTHROW(ERR_INVALID_DATA," CorrelUp = 0");
										return 0.;
										}
	//Correl Down
	virtual void ComputeEqCorrelDown();
	inline virtual void setEqCorrelDown(const double & CorrelDown) { itsEqCorrelDown = CorrelDown; itsEqCorrelDownFlg = true;}
	inline virtual double getEqCorrelDown(){ 
										if (itsEqCorrelDownFlg) return itsEqCorrelDown; 
										else ICMTHROW(ERR_INVALID_DATA," CorrelDown = 0");
										return 0.;
	}

/**
	virtual double Get_MixCopula_Factor (int SmileType, double BC1, double BC2, double Seed1, 
										 double Seed2,double Accuracy, int FactorType);

	virtual double Get_ReducedMixCopula_Factor(int SmileType, double BC1, double BC2, double SeedIndep,
											   double SeedBeta, int FactorType) ;
**/ 

 



	void View(char* id, FILE* ficOut);
	virtual double DoPrice(qCMPMETH measure);
	virtual double ComputeFlatCorrel(){ICMTHROW(ERR_INVALID_OPERATION,"Not Implemented"); }

};

#endif

