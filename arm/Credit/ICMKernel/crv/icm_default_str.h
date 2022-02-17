#ifndef _ICM_STR_DEFCURVE_H
#define _ICM_STR_DEFCURVE_H

#include "ARMKernel\util\refvalue.h"
#include "ARMKernel\crv\zeroflat.h"
#include "ICMKernel\glob\icm_enums.h"
#include "ICMKernel\crv\ICM_Constant_Piecewise.h"
#include "ICMKernel\str\icm_util_str.h"


class ICM_DefCurvStr;
class ICM_DefaultCurve;
 

class ICM_DefCurvStr : public ICM_DefaultCurve
{
private :
	double itsPente;			// Pente annuelle 
	double itsSpreadDate;
public :
	// --------------------  CONSTRUCTEURS ------------------------------
	// Constructeur par défaut
	ICM_DefCurvStr() { Init();}
	ICM_DefCurvStr(const ICM_DefCurvStr&ref); 
	
	// Constructeur : Calibration
	ICM_DefCurvStr(ICM_DefaultCurve* dc);
	ICM_DefCurvStr( const ARM_Vector* Spread,
					const ARM_Vector* YF_plot,
					double Recovery,
					ARM_ZeroCurve* ZCurv,
					const std::string& label);
	
	void Set(const ARM_Vector* Spread,
			 const ARM_Vector* YF_plot,
			 double Recovery,
			 ARM_ZeroCurve* ZCurv,
			 const std::string& label);
	
	
	// Constructeur & Set avec un vecteur de Date or YF
	ICM_DefCurvStr(const vector<double>& Spread,
				   const vector<double>& YForDates,
				   const double& Recovery,
				   ARM_ZeroCurve* ZCurv,
				   const std::string& Label,
				   const bool& IsDates);

	void Set(const vector<double>& Spread,
			 const vector<double>& YForDates,
			 const double& Recovery,
			 ARM_ZeroCurve* ZCurv,
			 const std::string& Label,
			 const bool& IsDates);
 

	// Get et Set
	double GetPente()	{ return itsPente; }
	void SetPente(double value) {itsPente = value;}

	double GetSpreadDate() { return itsSpreadDate;}
	void SetSpreadDate(double value) {itsSpreadDate = value;}
	// -------- Methodes de bases ---------------
	void Init();

// private:void BitwiseCopy(const ARM_Object* src);
public:
   //  virtual void Copy(const ARM_Object* srcZc);
    virtual ARM_Object* Clone(void);

	virtual ~ICM_DefCurvStr(void) {}
	virtual void View(char* id, FILE* ficOut);
	// ------------------ Implémentation spécifique à la classe ------------------------------
	virtual double SurvivalFunction(const double& yearterm) const ;
	ICM_Constant_Piecewise* ConvertStrToARM(void);
	virtual ICM_DefaultCurve* GenerateShiftCurve(const std::vector<std::string>&pTerm, 
												 const ARM_Vector&epsilon ,
												 qSENSITIVITY_TYPE mode  ) const ;
/** JLA Removed 
private:
	virtual ICM_DefaultCurve* xGenerateShiftCurve(double epsilon, qSENSITIVITY_TYPE mode) const ;
	**/ 
public:

	virtual void Calibrate(/** ICM_DefCurvStr* defcurve = NULL**/ );
	void CptSProbaPlot();
	virtual double FwdSpread(const ARM_Date& t1, const ARM_Date& t2) const ;
	// virtual double FwdSpread(const double& t1, const double& t2);

	virtual ICM_DefaultCurve* GenDefCurve(const ARM_Vector& dates,
										  const ARM_Vector& spreads,
										  const double& recovery,
										  const std::string& label,
										  const bool& isdates,
										  ARM_ZeroCurve* ircurve) const ; 

private:
	virtual double DefProbInverse(const double& PriceIn) const 
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_DefCurvStr::DefProbInverse: not implemented"); 
	}

private:
	
	ICM_DefCurvStr& operator=(const ICM_DefCurvStr&); //NA 
};

#endif
