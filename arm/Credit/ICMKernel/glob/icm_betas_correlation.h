/********************************************************************************/
/*! \file icm_corrmatrix.h
 * 
 *  \brief Describes a Correlation Matrix Object
 *  \author 
 *	\version 1.0
 *	\date   April 2004 */
/*
 *********************************************************************************/

#ifndef _ICM_BETAS_CORRELATION_H_
#define _ICM_BETAS_CORRELATION_H_

#include "ICMKernel/glob/icm_correlation.h"

/*********************************************************************************/
/*! \class  ICM_Beta_Correlation icm_betas_correlation.h "icm_betas_correlation.h"
 *  \author D Pouponneau
 *	\version 1.0
 *	\date   April 2004
 *	\file   icm_betas_correlation.h
 *	\brief Creates a beta correlation object */
/***********************************************************************************/

class ICM_Beta_Correlation : public ICM_Correlation
{
private:
	
	ARM_Vector* itsBetas;
	double itsFixedBeta;

protected:
	void Init(void) ;
	

public: 
	ARM_Vector* GetBetas(void) const {return itsBetas;}
	void Set(const ARM_Date& AsOf,
		const std::string& name,
			 const ARM_Vector& betas, 
			 const std::vector<std::string>& labels,
			 const ARM_IRIndex*index1,const ARM_IRIndex*index2) ;


	void SetBetas(const ARM_Vector& betas) ;
	void SetBetas(const ARM_Vector& betas,const std::vector<std::string>&labels) ;
	ARM_Vector* GetBetas(const std::vector<std::string>&labels,int size) const ;


	ICM_Beta_Correlation() {Init();}		
	ICM_Beta_Correlation(const ARM_Date& AsOf,
						const string& name, 
						const ARM_Vector& betas, 
						const std::vector<std::string>&,
						const ARM_IRIndex*index1,
						const ARM_IRIndex*index2) ;	
	
	ICM_Beta_Correlation(const ARM_Date& AsOf,
						const double& beta_fixed,
						const string& name,
						const ARM_IRIndex*index1,
						const ARM_IRIndex*index2) ;


	virtual double GetCorrelation(const std::string&  issuer1,
								  const std::string&  issuer2,
								  double maturity = CREDIT_DEFAULT_VALUE,
								  double strike = CREDIT_DEFAULT_VALUE,
								  double actualYF = CREDIT_DEFAULT_VALUE) ;


	virtual double GetBeta(const std::string& issuer,
						   double maturity = CREDIT_DEFAULT_VALUE,
						   double strike = CREDIT_DEFAULT_VALUE,
						   double actualYF = CREDIT_DEFAULT_VALUE) ;

	virtual ICM_QMatrix<double> ComputeCholeskyMatrix() ;


	~ICM_Beta_Correlation() ;

	// ----------------------------
	//	Copy of members data
	// ----------------------------
	void BitwiseCopy(const ARM_Object* src) ;


	double GetFixedBeta() { return itsFixedBeta;}
	void SetFixedBeta(double beta) {itsFixedBeta = beta;}

	// -------------
	//	Copy Method 
	// -------------
	void Copy(const ARM_Object* src) ;
	// --------------
	//	Clone Method
	// --------------
	ARM_Object* Clone(void) ;


public:

	void View(char* id, FILE* ficOut);

	virtual void ResetBetas() {}
	virtual ICM_Correlation* GenerateShiftCorrel(const std::string& label,
												qSENSITIVITY_TYPE typesensi,
												double epsilon /** = CREDIT_DEFAULT_VALUE **/ ) ; 


	virtual ICM_Correlation* GenerateShiftBetas(const std::string& label,
												qSENSITIVITY_TYPE typesensi,
												double epsilon = CREDIT_DEFAULT_VALUE) ;

};



#endif