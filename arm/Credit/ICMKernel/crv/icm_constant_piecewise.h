#ifndef _ICM_CONSTANT_PIECEWISE_CURVE_H
#define _ICM_CONSTANT_PIECEWISE_CURVE_H

#include "ICMKernel\crv\icm_defaultcurve.h"


class ICM_Constant_Piecewise : public ICM_DefaultCurve
{
public:
	ICM_Constant_Piecewise(const ARM_Date& asOf,
				const std::vector<std::string>& terms,
					 ARM_Vector* rates,
					 double Recovery,
					 ARM_ZeroCurve* zc,
					 int intRule,
					 int adjStartDate,
					 qCDS_ADJ adj ,  	
					 const std::string&  ccy  ,
					 const std::string& label  ,
					 bool issummitcurve ,
					 const ARM_VolCurve* VolCurve ,
					 long PayFreq ,
					 // bool	isBrentCalib ,
					 qDEFCURVE_CALIB_ALGO	calibAlgo,
					 const std::string& calibData,
					 int LagStartDate,
					 const ICM_Parameters& params,
					 ARM_Vector* UpFronts=NULL,
					 vector<ARM_ReferenceValue*>& AmortRefValues = vector<ARM_ReferenceValue*>());

	void Set (const ARM_Date& asOf,
		const std::vector<std::string>& terms,
                   ARM_Vector* rates,
				   double Recovery,
				   ARM_ZeroCurve* zc,
					 int intRule,
					 int adjStartDate,
				   qCDS_ADJ adj ,  	
                   const std::string&  ccy ,
				   const std::string& label ,
				   bool issummitcurve ,
				   const ARM_VolCurve* VolCurve /*=NULL*/,
				   long PayFreq /*=K_QUARTERLY*/,
					qDEFCURVE_CALIB_ALGO	calibAlgo,
					const std::string& calibData,
				   int LagStartDate,
				   const ICM_Parameters&params,
				   ARM_Vector* UpFronts=NULL,
				   vector<ARM_ReferenceValue*>& AmortRefValues = vector<ARM_ReferenceValue*>());

	ICM_Constant_Piecewise(const ARM_Date& asOf,
					 const ARM_Vector& Dates,
					 const ARM_Vector& Inputs,
					 int Intensity_or_Survivalprobas,
					 double Recovery,
					 ARM_ZeroCurve* zc,  
					 int intRule,
					 int adjStartDate,
					 qCDS_ADJ adj ,  
					 const std::string& ccy,
					 const std::string& label,
					 const ARM_VolCurve* VolCurve ,
					 qDEFCURVE_CALIB_ALGO	calibAlgo,
					 const std::string& calibData,
					 int LagStartDate,
					 const ICM_Parameters&); 

	ICM_Constant_Piecewise(const ARM_Date& asOf,
					 double* YearFractions,
					 ARM_Vector* Inputs,
					 int Intensity_or_Survivalprobas,
					 double Recovery,
					 ARM_ZeroCurve* zc,  
					int intRule,
					int adjStartDate,
					qCDS_ADJ adj ,  
					const std::string&  ccy ,
					const std::string& label ,
					 const ARM_VolCurve* VolCurve ,
					 qDEFCURVE_CALIB_ALGO	calibAlgo,
					 const std::string& calibData,
					 int LagStartDate);

	ICM_Constant_Piecewise(const ARM_Date& asOf,
					 const ARM_Vector& YForDates,
					 const ARM_Vector& rates,
					 double Recovery,
					 ARM_ZeroCurve* zc,
					 int intRule,
					 int adjStartDate,
					 qCDS_ADJ adj ,
					 const std::string&  ccy ,
					 const std::string& label ,
					 bool issummitcurve ,
					 const ARM_VolCurve* VolCurve ,
					 const bool& isdatesininput ,
					 long PayFreq ,
					 qDEFCURVE_CALIB_ALGO	calibAlgo,
					 const std::string& calibData,
					 int LagStartDate);

	void Set (const ARM_Date& asOf,
                   const ARM_Vector& YForDates,
				   const ARM_Vector& rates,
				   double Recovery,
				   ARM_ZeroCurve* zc,
					 int intRule,
					 int adjStartDate,
					 qCDS_ADJ adj ,
					 const std::string&  ccy ,
					 const std::string& label ,
				   bool issummitcurve ,
				   const ARM_VolCurve* VolCurve ,
				   const bool& isdatesininput ,
				   long PayFreq ,
					 qDEFCURVE_CALIB_ALGO	calibAlgo,
					 const std::string& calibData,
				   int LagStartDate);
public :
    ICM_Constant_Piecewise(void) { Init();}

	ICM_Constant_Piecewise(const ICM_Constant_Piecewise&ref); 
private:
	void Init();
public:
	virtual ~ICM_Constant_Piecewise(void){}

protected:
	virtual double SurvivalFunction(const double& yearterm) const ;
	virtual void CptTermsSurvivalProba(void);

protected:
	virtual void ResetLambda(const int& indice,const double& lambda);
public:

	virtual void View(char* id = NULL, FILE* ficOut = NULL );

	virtual ICM_DefaultCurve* GenerateShiftCurve(const std::vector<std::string>&  pTerm, 
											  const ARM_Vector& epsilon ,
											  qSENSITIVITY_TYPE mode ) const ;


public :
	virtual double DefProbInverse(const double& DefaultProba) const ;

public:
 
    virtual ARM_Object* Clone(void);
 
private:
	ARM_Vector* CptLambdas(const ARM_Vector& InputSP, const ARM_Vector& Dates) const ;
	ARM_Vector* CptLambdas(const ARM_Vector& InputSP, double* InputYF) const;
public:
	virtual ICM_DefaultCurve* GenDefCurve(const ARM_Vector& dates,
										  const ARM_Vector& spreads,
										  const double& recovery,
										  const string& label,
										  const bool& isdates,
										  ARM_ZeroCurve* ircurve) const ;
protected:
	virtual void Calibrate() 
	{
		ICM_DefaultCurve::Calibrate(); 
	}
private: 
	
	ICM_Constant_Piecewise& operator=(const ICM_Constant_Piecewise&); //NA 
};


#endif /*---- End of file ----*/
