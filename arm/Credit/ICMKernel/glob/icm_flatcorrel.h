#ifndef FLAT_CORREL_H_
#define FLAT_CORREL_H_

#include "ICMKernel/glob/icm_correlation.h"

class ICM_FlatCorrel : public ICM_Correlation
{
	double itsCorrel; 
public:
	ICM_FlatCorrel(const ARM_Date& AsOf,const std::string& structName,
		const ARM_IRIndex*ix1,const ARM_IRIndex*ix2,double correlValue); 
	void Init(void);
	ICM_FlatCorrel();
// FIXMEFRED: mig.vc8 (28/05/2007 10:23:53):missing return type
	void Set(const ARM_Date& AsOf,const std::string& structName,
		const ARM_IRIndex*ix1,const ARM_IRIndex*ix2,double correlValue); 
	ICM_FlatCorrel(const ICM_FlatCorrel&) ; 
	ICM_FlatCorrel& operator=(const ICM_FlatCorrel&); 
	~ICM_FlatCorrel() ;
	virtual ARM_Object* Clone(void) const ; 
	virtual ARM_Object* Clone() ; 
	void BitwiseCopy(const ARM_Object* src);
	virtual void Copy(const ARM_Object* src);
	virtual double GetCorrelation(const std::string& issuer1,
								  const std::string& issuer2,
								  double maturity ,
								  double strike ,
								  double actualYF) ;
	virtual double GetCorrelation(double yf) const ; 
	virtual ICM_Correlation* GenerateShiftCorrel(const std::string& label,
												 qSENSITIVITY_TYPE typesensi,
												 double epsilon ) ;
	virtual ICM_Correlation* GenerateShiftBetas(const std::string& label,
												qSENSITIVITY_TYPE typesensi,
												double epsilon  ) ;
	void View(char* id, FILE* ficOut) ;

	virtual double GetBeta(const std::string& issuer,
						   double maturity = CREDIT_DEFAULT_VALUE,
						   double strike = CREDIT_DEFAULT_VALUE,
						   double actualYF = CREDIT_DEFAULT_VALUE) 
	{
		return sqrt(itsCorrel);
	}
};

// -------------------------------------------------------------------------------------------------
inline double 
ICM_FlatCorrel::GetCorrelation(const std::string& issuer1,
								  const std::string& issuer2,
								  double maturity ,
								  double strike ,
								  double actualYF)
{
	return itsCorrel ; 
}
// -------------------------------------------------------------------------------------------------
inline double 
ICM_FlatCorrel::GetCorrelation(double yf) const
{
	return itsCorrel ; 
}
// -------------------------------------------------------------------------------------------------

#endif	// FLAT_CORREL_H_
