/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		ICM_Correlation_Fitting.H
	PROJECT:	UTIL
	
	DESCRIPTION:	Fit a Correlation Structure

  -----------------------------------------------------------------

 	CREATION:	February, 2006

	LAST MODIF:	February, 2006
  -----------------------------------------------------------------
   
	ICM Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */

#ifndef __ICM_Correlation_Fitting_H__
#define __ICM_Correlation_Fitting_H__

#include "ICMKernel\glob\icm_corrmatrix.h"

class ICM_Correlation_Fitting : public ARM_Object
{
	public:
	
		int		its_dim;

		int		its_NbMaxIter;
		int		its_NbIter;

		// -------------------------------------------------
		// -------------------------------------------------

		// -------------------------------------------------
		// 2F
		vector<int>		its_Sector_Membership;
		// -------------------------------------------------

		vector<double>	its_Steps;
		vector<double>	its_correl_parameters;

		int		its_CorrelationSize;
		ICM_CorrMatrix*			its_CorrelationMatrix;			// Correlation Matrix Object

		qCORRELATION_FIT_TYPE	its_Correlation_Fit_Type;

		// -------------------------------------------------
		// -------------------------------------------------

		// OUTPUTS

		double			its_1F_Single_Value;
		vector<double>	its_1F_Betas;

		double			its_2F_Inter_Sector_Value;
		double			its_2F_Intra_Sector_Value;
		vector<double>	its_2F_Betas;
		vector<double>	its_2F_Lambdas;

		// -------------------------------------------------
		// -------------------------------------------------

	public:

		// FUNCTIONS
	inline void Init(void)
	{
		its_CorrelationSize	=	0;
		its_CorrelationMatrix	=	NULL;

		its_NbMaxIter	=	100;

		its_Sector_Membership.clear();

		its_correl_parameters.clear();
		its_Steps.clear();

		its_Correlation_Fit_Type	=	qCORREL_1F_SINGLE_VALUE;

		its_1F_Single_Value	=	0.0;
		its_1F_Betas.clear();

		its_2F_Inter_Sector_Value	=	0.0;
		its_2F_Intra_Sector_Value	=	0.0;
		its_2F_Betas.clear();
		its_2F_Lambdas.clear();
	}

	public: 

	ICM_Correlation_Fitting() {Init();}

	void Set()
	{
	}

	~ICM_Correlation_Fitting() 
	{

		its_Sector_Membership.clear();

		its_correl_parameters.clear();
		its_Steps.clear();

		if (its_CorrelationMatrix)
			delete its_CorrelationMatrix;
		its_CorrelationMatrix	=	NULL;
		
		its_1F_Betas.clear();

		its_2F_Betas.clear();
		its_2F_Lambdas.clear();
	}

	// ----------------------------
	//	Copy of members data
	// ----------------------------
	void BitwiseCopy(const ARM_Object* src)
	{
	    ICM_Correlation_Fitting* Correl = (ICM_Correlation_Fitting *) src;
		
		its_Sector_Membership	=	Correl->its_Sector_Membership;

		its_CorrelationSize	=	Correl->its_CorrelationSize;
		its_correl_parameters	=	Correl->its_correl_parameters;
		its_Steps		=	Correl->its_Steps;
		its_dim			=	Correl->its_dim;

		its_NbMaxIter	=	Correl->its_NbMaxIter;
		its_NbIter		=	Correl->its_NbIter;

		its_Correlation_Fit_Type		=	Correl->its_Correlation_Fit_Type;

		its_1F_Single_Value	=	Correl->its_1F_Single_Value;
		its_1F_Betas		=	Correl->its_1F_Betas;

		its_2F_Inter_Sector_Value	=	Correl->its_2F_Inter_Sector_Value;
		its_2F_Intra_Sector_Value	=	Correl->its_2F_Intra_Sector_Value;
		its_2F_Betas	=	Correl->its_2F_Betas;
		its_2F_Lambdas	=	Correl->its_2F_Lambdas;
	}

	// -------------
	//	Copy Method 
	// -------------
	void Copy(const ARM_Object* src)
	{
		ARM_Object::Copy(src);
		BitwiseCopy(src);
	}

	// --------------
	//	Clone Method
	// --------------
	ARM_Object* Clone(void)
	{
		 ICM_Correlation_Fitting* theClone = new ICM_Correlation_Fitting();
		 theClone->Copy(this);
		 return(theClone);
	}

	void View(char* id, FILE* ficOut);

	public:

		ICM_CorrMatrix*	Get_CorrelationMatrix(void) { return its_CorrelationMatrix;}

		void	Set_CorrelationMatrix(ICM_CorrMatrix* correlmatrix)
		{
			if (its_CorrelationMatrix)
				delete its_CorrelationMatrix;
			its_CorrelationMatrix = (ICM_CorrMatrix*) correlmatrix->Clone();
		}

		void	Set_Correlation_Fit_Type(const qCORRELATION_FIT_TYPE& value){its_Correlation_Fit_Type = value;}
		void	Get_Correlation_Fit_Type(qCORRELATION_FIT_TYPE& value){value = its_Correlation_Fit_Type;}

		void	Set_Sector_Membership(vector<int>& value){its_Sector_Membership = value;}
		void	Get_Sector_Membership(vector<int>& value){value = its_Sector_Membership;}
		int		GetSectorId(int Num){return its_Sector_Membership[Num];}

	public:

		// COMPUTATION
		void	Compute_BestFit_CorrelationMatrix();

		// GET RESULTS
		void	Get_BestFit_Results(vector<double>& data);

	protected:

		double	TheBetaBestFitFunction(double* X);
		void	TheBetaBestFitGradientFunction(double* X, double* Gradient);
		double	The2FSameInterSameIntraBestFitFunction(double* X);
		void	The2FSameInterSameIntraBestFitGradientFunction(double* X, double* Gradient);

		void	ExtractRootsForBestFit(double* X); //, vector<double>& Correl_Parameters)

};


#endif