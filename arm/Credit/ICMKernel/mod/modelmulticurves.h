
#ifndef _MODEL_MULTI_CURVES_H
#define _MODEL_MULTI_CURVES_H




#include "ICMKernel\mod\icm_defcurvemodel.h"
#include "ICMKernel\util\icm_qmatrix.h"

namespace ARM 
{
	class ARM_InfCurv; 
}; 

class ICM_Credit_Index;
class ICM_MktDataMng;
class ICM_Correlation; 


/*********************************************************************************/
/*! \class  ICM_ModelMultiCurves modelmulticurves.h "modelmulticurves.h"
 *  \author Damien Pouponneau 
 *	\version 1.0
 *	\date   June 2003
 *	\brief  Model Multi Curves 
 *	Cloning Enabled
 *		Some inputs are cloned and the Model takes ownership of the inputs
 *	Cloning Disabled
 *		Inputs are not clonded, the Model is not owner of the inputs
 *		BUT it remains owner of the container structure (asif a vector<>) */
/***********************************************************************************/

class ICM_ModelMultiCurves : public ICM_DefaultCurveModel        
{        
    protected:

		// int				   itsNbDefCurves;		//Number of default curves
		// const ICM_DefaultCurve** itsDefaultCurves;	//Vector of default curves
		std::vector<const ICM_DefaultCurve*> itsDefaultCurves; 
		ARM_Vector		itsRecoveryRates; 
		ARM_CorrelManager* itsCorrelManager;	//Correlation manager for IR products

		double			   itsIndependantPart;
		double			   itsFullCorrelPart;
	private : 
		ICM_MktDataMng*		itsMktDataMng;			// conteneur of Inflation Curve and interest rates curve et autres courbes
		ICM_Correlation*   itsCorrelation;		//0..1, Correlation object
    public:

		

        ICM_ModelMultiCurves(void)  {Init();}

		ICM_ModelMultiCurves(/** int NbDefCurves,
							const ICM_DefaultCurve** DefaultCurves, **/ 
							const std::vector<const ICM_DefaultCurve*> defCurves,
							ARM_ZeroCurve* DiscountCurve,
							const ARM_Vector& RecoveryRates,
							ICM_Correlation* Correlation = NULL,
							ARM_VolCurve* volcurve = NULL,
							bool FlgClone = true,
							ARM::ARM_InfCurv* infcurve = NULL,
							ARM_ZeroCurve* CpnIRCurve = NULL,
							ARM_CorrelManager* correlmng = NULL):ICM_DefaultCurveModel(defCurves[0],DiscountCurve,volcurve,FlgClone)
		{
			Init();
			Set(defCurves,/** NbDefCurves,DefaultCurves,**/  DiscountCurve,RecoveryRates,Correlation,infcurve,CpnIRCurve);
		}
		/**  JLA: non const version now useless.

		ICM_ModelMultiCurves(int NbDefCurves,
							ICM_DefaultCurve** DefaultCurves,  
							const std::vector<const ICM_DefaultCurve*> defCurves,
							ARM_ZeroCurve* DiscountCurve,
							const ARM_Vector& RecoveryRates,
							ICM_Correlation* Correlation = NULL,
							ARM_VolCurve* volcurve = NULL,
							bool FlgClone = true,
							ARM::ARM_InfCurv* infcurve = NULL,
							ARM_ZeroCurve* CpnIRCurve = NULL,
							ARM_CorrelManager* correlmng = NULL):ICM_DefaultCurveModel(DefaultCurves[0],DiscountCurve,volcurve,FlgClone)
		{
			Init();
			Set(NbDefCurves,(const ICM_DefaultCurve**)DefaultCurves,DiscountCurve,RecoveryRates,Correlation,infcurve,CpnIRCurve);
		}
		**/ 
		/** JLA TEMP 
		ICM_ModelMultiCurves(int NbDefCurves,
							ICM_DefaultCurve** DefaultCurves,
							ARM_ZeroCurve* DiscountCurve,
							const ARM_Vector& RecoveryRates,
							const ICM_QMatrix<double>&Correlation ,
							ARM_VolCurve* volcurve=NULL,
							bool FlgClone = true,
							ARM::ARM_InfCurv* infcurve = NULL,
							ARM_ZeroCurve* CpnIRCurve = NULL,
							ARM_CorrelManager* correlmng = NULL):ICM_DefaultCurveModel(DefaultCurves[0],DiscountCurve,volcurve,FlgClone)
		{
			Init();
			Set(NbDefCurves,DefaultCurves,DiscountCurve,RecoveryRates,Correlation,infcurve,CpnIRCurve);
			}
			**/ 

		ICM_ModelMultiCurves(/** int NbDefCurves,
							 const ICM_DefaultCurve** DefaultCurves, **/ 
							 const std::vector<const ICM_DefaultCurve*> defCurves,
							ARM_ZeroCurve* DiscountCurve, // useless
							const ARM_Vector& RecoveryRates,
							ICM_Correlation* Correlation,
							ICM_MktDataMng*		pMktDataMng,
							ARM_VolCurve* volcurve=NULL,
							bool FlgClone = true,
							ARM_CorrelManager* correlmng = NULL):ICM_DefaultCurveModel(defCurves[0],DiscountCurve,volcurve,FlgClone)
		{
			Init();
			Set(defCurves,/** NbDefCurves,DefaultCurves,**/ DiscountCurve,RecoveryRates,Correlation, pMktDataMng);
		}
		/** JLA: non const version now useless.
		ICM_ModelMultiCurves( int NbDefCurves,
							 ICM_DefaultCurve** DefaultCurves,  
							 const std::vector<const ICM_DefaultCurve*> defCurves,
							ARM_ZeroCurve* DiscountCurve, // useless
							const ARM_Vector& RecoveryRates,
							ICM_Correlation* Correlation,
							ICM_MktDataMng*		pMktDataMng,
							ARM_VolCurve* volcurve=NULL,
							bool FlgClone = true,
							ARM_CorrelManager* correlmng = NULL):ICM_DefaultCurveModel(DefaultCurves[0],DiscountCurve,volcurve,FlgClone)
		{
			Init();
			Set(defCurves, NbDefCurves,(const ICM_DefaultCurve**)DefaultCurves, DiscountCurve,RecoveryRates,Correlation, pMktDataMng);
		} **/ 

		void SetCorrelation(ICM_Correlation* correl) ; 
		void SetCorrelationPtr(ICM_Correlation* correl) { itsCorrelation = correl;}

		ICM_Correlation* GetCorrelation(void) const { return itsCorrelation;}

		// ICM_ModelMultiCurves(ARM_ZeroCurve* DiscountCurve);
		ICM_ModelMultiCurves(const ICM_ModelMultiCurves& multi);
        virtual ~ICM_ModelMultiCurves(void);

//         void BitwiseCopy(const ARM_Object* src);

//         virtual void Copy(const ARM_Object* src);  
        virtual ARM_Object* Clone(void); 

		// virtual bool IsCreditModel(void) { return true; };

		void View(char* id, FILE* ficOut);

		void Set(/** int NbDefCurves,
				 const ICM_DefaultCurve** DefaultCurves, **/ 
				 const std::vector<const ICM_DefaultCurve*>& defCurves,
				 ARM_ZeroCurve* DiscountCurve,
				 const ARM_Vector& RecoveryRates,
				 ICM_Correlation* Correlation = NULL,
				 ARM::ARM_InfCurv* infcurve = NULL,
				 ARM_ZeroCurve* CpnIRCurve = NULL,
				 ARM_CorrelManager* correlmng = NULL);

 
		void Set(/** int NbDefCurves,
				const ICM_DefaultCurve** DefaultCurves, **/ 
				const std::vector<const ICM_DefaultCurve*>& defCurves,
				ARM_ZeroCurve* DiscountCurve,
				const ARM_Vector& RecoveryRates,
				ICM_Correlation* Correlation,
				ICM_MktDataMng * pMktDataMng,
				ARM_CorrelManager* correlmng = NULL);
		
		unsigned int GetNbDefCurves() const { return itsDefaultCurves.size(); }
 
public:
		const ICM_DefaultCurve*	GetDefaultCurves(unsigned int i) const ; 
 		const ARM_Vector& GetRecoveryRates() const { return itsRecoveryRates; } 

 
		const ICM_DefaultCurve* GetDefaultCurve(const std::string& name) const 
		{
 			for (int i = 0; i<itsDefaultCurves.size(); i++) 
			{
				if (itsDefaultCurves[i]->GetLabel()==name) return itsDefaultCurves[i];
			}
			ICMTHROW(ERR_INVALID_ARGUMENT,"GetDefaultCurve: Can't get "<<name); 
			return NULL; // useless 
		}
		const ICM_DefaultCurve* GetDefaultCurve(int Num) const 
		{
			if (Num<itsDefaultCurves.size() && Num>=0 ) return itsDefaultCurves[Num] ;
			ICMTHROW(ERR_INVALID_ARGUMENT,"GetDefaultCurve: Can't get curve i="<<Num); 
			return NULL; //useless
		}

		// LJ -> CC
		void SetDefaultCurve(int Num, const ICM_DefaultCurve* data);

		double GetRecoveryRate(const std::string& label) const ;

		double GetIndependantPart() const { return itsIndependantPart;}
		double GetFullCorrelPart() const { return itsFullCorrelPart;}

		void SetIndependantPart(double value) { itsIndependantPart = value;}
		void SetFullCorrelPart(double value) { itsFullCorrelPart = value;}

		// last Flag in order to allow Parallel Shift computation
		// char** GetUnionTenorsForCollateral(const std::vector<std::string>& labels,int& sizeout, bool ParallelFlag = false);
		void GetUnionYFForCollateral(const std::vector<std::string>&labels,ARM_Vector& Output) const ;

		void CptNoCurve (const std::string& ,int& outNoCurve) const ;


		ICM_ModelMultiCurves* GenerateShiftModel(qSENSITIVITY_TYPE typesensi,
			const std::string& plot, 
			const std::string& label,
			int& nocurve,
			double epsilon );

		// void   SetDefaultCurves(const ICM_DefaultCurve** curves);
		void   SetDefaultCurves(const std::vector<const ICM_DefaultCurve*>& curves);
		// void   SetDefaultCurves(ICM_DefaultCurve** curves) ;

		void SetCpnInfCurve(ARM::ARM_InfCurv* infcurve) ;
		ARM::ARM_InfCurv* GetCpnInfCurve(void) const ;

		void SetCpnIRCurve(ARM_ZeroCurve* ircurve) ;
		ARM_ZeroCurve* GetCpnIRCurve(void) const ;

		void SetCorrelManager(ARM_CorrelManager* CorrelManager) ; 

		const ARM_CorrelManager* GetCorrelManager(void) const {return itsCorrelManager;}

		ICM_MktDataMng* GetMktDataMng() const { return itsMktDataMng;}

		void SetMktDataMng(ICM_MktDataMng* pMktDataMng) ;
protected:
	void Init(); // protected because customized credit multicurve model ?
private:
	
	ICM_ModelMultiCurves& operator=(const ICM_ModelMultiCurves&ref); //NA 
};


inline const ICM_DefaultCurve*	
ICM_ModelMultiCurves::GetDefaultCurves(unsigned int i) const
{
	if (i>=itsDefaultCurves.size()) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_ModelMultiCurves::GetDefaultCurves: out of bonds "<<i); 
	return itsDefaultCurves[i]; 
}
#endif
