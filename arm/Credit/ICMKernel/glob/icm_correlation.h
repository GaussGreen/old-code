/********************************************************************************/
/*! \file icm_correlation.h
 * 
 *  \brief Describes a Correlation Object
 *  \author 
 *	\version 1.0
 *	\date   October 2004 */
/*
 *********************************************************************************/

#ifndef _ICM_CORRELATION_H_
#define _ICM_CORRELATION_H_

/*********************************************************************************/
/*! \class  ICM_Correlation icm_corrmatrix.h "icm_corrmatrix.h"
 *  \author D Pouponneau
 *	\version 1.0
 *	\date   April 2004
 *	\brief Creates a Basis correlation object 
/***********************************************************************************/

#include "ARMKernel\glob\firsttoinc.h"
#include "ARMKernel\glob\linalg.h"
#include "ARMKernel\inst\security.h"
#include "ICMKernel\glob\icm_enums.h"
#include "ICMKernel\util\icm_qmatrix.h"


class ICM_Pricer;
class ARM_IRIndex; 

class ICM_Correlation : public ARM_Object
{
private:
	ARM_Date itsAsOf;
	std::string  itsStructName;
	std::vector<std::string> itsLabels; 
	ARM_IRIndex* itsIndex1;						//	optionnal, aggregated
	ARM_IRIndex* itsIndex2;						//	optionnal, aggregated
public:
	
	const ARM_Date& GetAsOfDate(void) const {return itsAsOf;}
protected:
	void Init(void) ; 
	

	void Set(const ARM_Date& AsOf,
		const std::vector<std::string>& labels,
		const std::string& name, 
		const ARM_IRIndex* index1,
		const ARM_IRIndex* index2); 
	ICM_Correlation() 
	{
		Init();
	}		
	ICM_Correlation(const ARM_Date& AsOf,
					const std::vector<std::string>& labels,
					const std::string& name, 
					const ARM_IRIndex* index1,
					const ARM_IRIndex* index2) ;

	void SetIndex1(const ARM_IRIndex* index);
	void SetIndex2(const ARM_IRIndex* index);
	void SetAsOf(const ARM_Date& aDate);
	void SetLabels(const std::vector<std::string>& labels);

public:
	virtual ARM_CLASS_NAME GetRootName(void) ;

	virtual void SetProportionsInfos(const std::string&indexname,
							 const double& proportion,
							 const double& forcedstrikelow = CREDIT_DEFAULT_VALUE,
							 const double& forcedstrikehigh = CREDIT_DEFAULT_VALUE){};
	

	virtual double GetCorrelation(const std::string& issuer1,
								  const std::string& issuer2,
								  double maturity = CREDIT_DEFAULT_VALUE,
								  double strike = CREDIT_DEFAULT_VALUE,
								  double actualYF = CREDIT_DEFAULT_VALUE)
	{
		return CREDIT_DEFAULT_VALUE;
	}
	virtual double GetCorrelation(double yf) const 
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Correlation::GetCorrelation(yf) not implemented"); 
	}


	virtual void ResetBetas()
	{
		ICMTHROW(ERR_INVALID_MODEL,typeid(*this).name()<<": not implemented ResetBetas"); 
	}

	virtual double GetBeta(const std::string& issuer,
						   double maturity = CREDIT_DEFAULT_VALUE,
						   double strike = CREDIT_DEFAULT_VALUE,
						   double actualYF = CREDIT_DEFAULT_VALUE) 
	{
		return CREDIT_DEFAULT_VALUE;
	}

	ARM_Vector GetBetaVector(const std::vector<std::string>& labels , // char** labels, 
									  int size,
									  double maturity = CREDIT_DEFAULT_VALUE,
									  double strike = CREDIT_DEFAULT_VALUE);
	
	virtual ARM_Vector* ComputeBetas(int nbissuers,
									const std::vector<std::string>& labels,
									const ARM_Vector&nominals, 
									 const ARM_Date & Maturity,
 									 ARM_Model* model = NULL)
	{ 
		return new ARM_Vector(GetBetaVector(labels,nbissuers));
	}

	virtual ICM_QMatrix<double> ComputeCholeskyMatrix() {ICM_QMatrix<double> test; return test;}

	virtual void ComputeStrikesEq(ICM_Pricer* pricer,const qRescalType& rescal/*=qRescal_Std_Maturity*/){}
	virtual void ComputeStrikesEq_Equity(ICM_Pricer* pricer){}
	virtual void ComputeStrikesEq_ELoss(ICM_Pricer* pricer){}


	virtual ~ICM_Correlation() ;

	void BitwiseCopy(const ARM_Object* src) ;
	virtual void Copy(const ARM_Object* src) ;
	virtual ARM_Object* Clone(void) ; 

	inline int GetSize() { return itsLabels.size();}
	const std::vector<std::string>& GetLabels() const { return itsLabels; }

	const std::string& GetLabel(unsigned int i) const 
	{
		if (i>=itsLabels.size()) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"GetLabel: out of bounds "<<i); 
		return itsLabels[i] ;
	}

	int GetLabelNo(const std::string& ) const ;


	virtual ICM_Correlation* GenerateShiftCorrel(const std::string& label,
												 qSENSITIVITY_TYPE typesensi,
												 double epsilon = CREDIT_DEFAULT_VALUE)
	{
		ICMTHROW(ERR_INVALID_MODEL,typeid(*this).name()<<":Unimplemented GenerateShiftCorrel");
	}

	virtual ICM_Correlation* GenerateShiftBetas(const std::string& label,
												qSENSITIVITY_TYPE typesensi,
												double epsilon /**= CREDIT_DEFAULT_VALUE **/ )
	{
		ICMTHROW(ERR_INVALID_MODEL,typeid(*this).name()<<":Unimplemented GenerateShiftBetas");
	}

	void View(char* id, FILE* ficOut) ;

	virtual void SetForcedStrikeType(qCorrel_By_Strike type){return;}


	//Equivalent Strike Down
	virtual double GetEqStrikeDown(const std::string&indexname) {return 0.;}

	//Equivalent Strike Up
	virtual double GetEqStrikeUp(const std::string& indexname) {return 0.;}

	//Equivalent Correl Strike Down
	virtual double GetCorrelStrikeDown(double maturity) {return 0.;}

	//Equivalent Correl Strike Up
	virtual double GetCorrelStrikeUp(double maturity) {return 0.;}

	//Equivalent Smile Strike Down
	virtual void GetSmileStrikeDown(const std::string& indexname, vector<double>& vMatu, vector<double>& vStrikes){return;};

	//Equivalent Smile Strike Up
	virtual void GetSmileStrikeUp(const std::string& indexname, vector<double>& vMatu, vector<double>& vStrikes){return;};

	void SetStructName(const std::string& StructName) {itsStructName =  StructName;}
	const std::string GetStructName(void) const { return itsStructName;}

	virtual void GetCorrelationTerms(ARM_Vector& vector) {};

	virtual int GetIndexSize() { return 0;}
	virtual void SetInterpType(const int& interpType) {}
	virtual int GetInterpType() {return 0;}

	int		Get_IssuersRankFromLabel(const std::string&Label) const ;

public:
	const ARM_IRIndex*  GetIndex1() const { return itsIndex1; } 
	const ARM_IRIndex*  GetIndex2() const { return itsIndex2; } 
};


#endif
