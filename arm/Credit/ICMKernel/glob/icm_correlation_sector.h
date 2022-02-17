/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		ICM_Correlation_Sector.H
	PROJECT:	MOD
	
	DESCRIPTION:	Containor for a Sector Correlation

  -----------------------------------------------------------------

 	CREATION:	February, 2006

	LAST MODIF:	February, 2006
  -----------------------------------------------------------------
   
	ICM Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */

#ifndef __ICM_Correlation_Sector_H__
#define __ICM_Correlation_Sector_H__


#include "ICMKernel\glob\icm_enums.h"
#include "ICMKernel/glob/icm_correlation.h"


class ICM_Correlation_Sector : public ICM_Correlation
{
	public:
	
		// -------------------------------------------------
		// SECTORIAL APPROACH
		int			its_Nb_Sectors;
		// int			its_Nb_Names;
		
		vector<int>		its_Sector_Membership;

		//double			its_Single_intra_sector_correlation;
		//double			its_Single_inter_sector_correlation;

		vector<double>	its_Betas;
		vector<double>	its_Lambdas;

		vector<double>	its_Betas_Down;
		vector<double>	its_Lambdas_Down;

		// -------------------------------------------------

		qTWO_FACTORS_CORRELATION_TYPE	its_CorrelationType;

	inline void Init(void)
	{
		ICM_Correlation::Init(); 
		SetName(ICM_CORRELATION_SECTOR);
		
		its_Sector_Membership.clear();

		its_Betas.clear();
		its_Lambdas.clear();
		its_Betas_Down.clear();
		its_Lambdas_Down.clear();
	}

public: 

	ICM_Correlation_Sector() {Init();}

	//Constructeur Single Intra & Inter value

	ICM_Correlation_Sector(const ARM_Date& AsOf,
				const std::string& structName,
				qTWO_FACTORS_CORRELATION_TYPE  The_CorrelationType,
				const std::vector<std::string>& labels,
				const int&	The_Nb_Sectors,
				const vector<int>&	The_Sector_Membership,
				double		The_Single_intra_sector_correlation,
				double		The_Single_inter_sector_correlation)
	{
		Init();

		Set(AsOf,structName,The_CorrelationType, labels, The_Nb_Sectors, The_Sector_Membership,
			The_Single_intra_sector_correlation, The_Single_inter_sector_correlation);
	}	

	void Set(
				const ARM_Date& AsOf,
				const std::string& structName,
				qTWO_FACTORS_CORRELATION_TYPE  The_CorrelationType,
				const std::vector<std::string>& labels,
				const int&	The_Nb_Sectors,
				const vector<int>	The_Sector_Membership,
				double		The_Single_intra_sector_correlation,
				double		The_Single_inter_sector_correlation)
	{
		ICM_Correlation::Set(AsOf,labels,structName,(ARM_IRIndex*)0,(ARM_IRIndex*)0);


		its_CorrelationType	=	The_CorrelationType;

		its_Nb_Sectors	=	The_Nb_Sectors;
		

		its_Sector_Membership	=	The_Sector_Membership;

		//Coefficients 
		its_Betas.resize(labels.size());	
		its_Lambdas.resize(labels.size());	
		for (int i=0; i<labels.size(); i++)
		{
			its_Betas[i]	=	sqrt(The_Single_intra_sector_correlation);
			its_Lambdas[i]	=	sqrt(The_Single_inter_sector_correlation) / sqrt(The_Single_intra_sector_correlation);	
		}
		//Copy correl down
		its_Betas_Down	=	its_Betas;
		its_Lambdas_Down	=	its_Lambdas;
	}

	//Constructeur different Intra & same inter
	ICM_Correlation_Sector(
				const ARM_Date& AsOf,
				const std::string& structName,
				qTWO_FACTORS_CORRELATION_TYPE  The_CorrelationType,
				const std::vector<std::string>& labels,
				const int&	The_Nb_Sectors,
				const vector<int>	The_Sector_Membership,
				double		The_Single_intra_sector_correlation,		//Not relevant, for compatibility
				double		The_Single_inter_sector_correlation,		//Not relevant, for compatibility
				const vector<double>&	The_Betas,
				const vector<double>&	The_Lambdas)
	{
		Init();

		Set(AsOf,structName,The_CorrelationType, labels, The_Nb_Sectors, The_Sector_Membership, The_Betas, The_Lambdas);
	}	

	ICM_Correlation_Sector(
				const ARM_Date& AsOf,
				const std::string& structName,
				qTWO_FACTORS_CORRELATION_TYPE  The_CorrelationType,
				const std::vector<std::string>& labels,
				const int&	The_Nb_Sectors,
				const vector<int>	The_Sector_Membership,
				const vector<double>&	The_Betas,
				const vector<double>&	The_Lambdas)
	{
		Init();

		Set(AsOf,structName,The_CorrelationType, labels, The_Nb_Sectors, The_Sector_Membership, The_Betas, The_Lambdas);
	}	

	void Set(
				const ARM_Date& AsOf,
				const std::string& structName,
				qTWO_FACTORS_CORRELATION_TYPE  The_CorrelationType,
				const std::vector<std::string>& labels,
				const int&	The_Nb_Sectors,
				const vector<int>	The_Sector_Membership,
				//double&		The_Single_intra_sector_correlation,
				//double&		The_Single_inter_sector_correlation,
				const vector<double>&	The_Betas,
				const vector<double>&	The_Lambdas)
	{
		ICM_Correlation::Set(AsOf,labels, structName,(ARM_IRIndex*)0,(ARM_IRIndex*)0);

		its_CorrelationType	=	  The_CorrelationType;

		its_Nb_Sectors	=	The_Nb_Sectors;
		
		
		its_Sector_Membership	=	The_Sector_Membership;

		its_Betas	=	The_Betas;
		its_Lambdas	=	The_Lambdas;

		//Copy correl down
		its_Betas_Down	=	its_Betas;
		its_Lambdas_Down	=	its_Lambdas;

	}


	//Constructeur different Intra & same inter by Strike
	ICM_Correlation_Sector(
				const ARM_Date& AsOf,
				const std::string& structName,
				qTWO_FACTORS_CORRELATION_TYPE  The_CorrelationType,
				const std::vector<std::string>& labels,
				const int&	The_Nb_Sectors,
				const vector<int>	The_Sector_Membership,
				double		The_Single_intra_sector_correlation,				//Not relevant, for compatibility
				double		The_Single_inter_sector_correlation,				//Not relevant, for compatibility
				const vector<double>&	The_Betas,
				const vector<double>&	The_Lambdas,
				const vector<double>&	The_Betas_Down,
				const vector<double>&	The_Lambdas_Down)
	{
		Init();

		Set(AsOf,structName,The_CorrelationType, labels, The_Nb_Sectors,  The_Sector_Membership, The_Betas, The_Lambdas,The_Betas_Down,The_Lambdas_Down);
	}	

	ICM_Correlation_Sector(
				const ARM_Date& AsOf,
				const std::string& structName,
				qTWO_FACTORS_CORRELATION_TYPE  The_CorrelationType,
				const std::vector<std::string>& labels,
				const int&	The_Nb_Sectors,
				const vector<int>&	The_Sector_Membership,
				const vector<double>&	The_Betas,
				const vector<double>&	The_Lambdas,
				const vector<double>&	The_Betas_Down,
				const vector<double>&	The_Lambdas_Down)
	{
		Init();

		Set(AsOf,structName,The_CorrelationType, labels, The_Nb_Sectors, The_Sector_Membership, The_Betas, The_Lambdas,The_Betas_Down,The_Lambdas_Down);
	}	

	void Set(
				const ARM_Date& AsOf,
				const std::string& structName,
				qTWO_FACTORS_CORRELATION_TYPE  The_CorrelationType,
				const std::vector<std::string>& labels,
				const int&	The_Nb_Sectors,
				const vector<int>&	The_Sector_Membership,
				const vector<double>&	The_Betas,
				const vector<double>&	The_Lambdas,
				const vector<double>&	The_Betas_Down,
				const vector<double>&	The_Lambdas_Down)
	{
		ICM_Correlation::Set(AsOf,labels,structName, (ARM_IRIndex*)0,(ARM_IRIndex*)0);

		its_CorrelationType	=	 The_CorrelationType;

		its_Nb_Sectors	=	The_Nb_Sectors;
		// its_Nb_Names	=	The_Nb_Names;
		
		its_Sector_Membership	=	The_Sector_Membership;

		its_Betas	=	The_Betas;
		its_Lambdas	=	The_Lambdas;

		its_Betas_Down	=	The_Betas_Down;
		its_Lambdas_Down	=	The_Lambdas_Down;
	}



	~ICM_Correlation_Sector() 
	{
		 
		its_Sector_Membership.clear(); 
		its_Betas.clear();
		its_Lambdas.clear();
		its_Betas_Down.clear();
		its_Lambdas_Down.clear();
	};

	// ----------------------------
	//	Copy of members data
	// ----------------------------
	void BitwiseCopy(const ARM_Object* src)
	{
	    ICM_Correlation_Sector* Correl = (ICM_Correlation_Sector *) src;

		its_Nb_Sectors	=	Correl->its_Nb_Sectors;
		// its_Nb_Names	=	Correl->its_Nb_Names;

		//its_Single_intra_sector_correlation	=	Correl->its_Single_intra_sector_correlation;
		//its_Single_inter_sector_correlation	=	Correl->its_Single_inter_sector_correlation;

		its_Sector_Membership	=	Correl->its_Sector_Membership;
		its_Betas		=	Correl->its_Betas;
		its_Lambdas		=	Correl->its_Lambdas;
		its_Betas_Down		=	Correl->its_Betas_Down;
		its_Lambdas_Down		=	Correl->its_Lambdas_Down;
	
		its_CorrelationType	=	Correl->its_CorrelationType;

	}

	// -------------
	//	Copy Method 
	// -------------
	void Copy(const ARM_Object* src)
	{
		ICM_Correlation::Copy(src);
		BitwiseCopy(src);
	}

	// --------------
	//	Clone Method
	// --------------
	ARM_Object* Clone(void)
	{
		 ICM_Correlation_Sector* theClone = new ICM_Correlation_Sector();
		 theClone->Copy(this);
		 return(theClone);
	}

	void View(char* id, FILE* ficOut);

//----- Assessors

	public:

		inline void SetNbSectors(const int& value){its_Nb_Sectors = value;}
		inline void GetNbSectors(int& value)const{value = its_Nb_Sectors;}


		inline void SetSingle_Intra_Sector_Correlation(const double& value)
		{
			for (int i=0; i<GetLabels().size();i++)
			{
				its_Betas[i] = sqrt(value);
			}
			its_Betas_Down=its_Betas;
		}
		
		inline void SetSingle_Inter_Sector_Correlation(const double& value)
		{
			for (int i=0; i<GetLabels().size();i++)
			{
				its_Lambdas[i] = sqrt(value) / its_Betas[i];
			}
			its_Betas_Down=its_Betas;
		}
	
		inline void SetCorrelationType(const qTWO_FACTORS_CORRELATION_TYPE& value){its_CorrelationType = value;}
		inline void GetCorrelationType(qTWO_FACTORS_CORRELATION_TYPE& value){value = its_CorrelationType;}

		void	Set_Sector_Membership(const vector<int>& value){its_Sector_Membership = value;}
		void	Get_Sector_Membership(vector<int>& value) const {value = its_Sector_Membership;}
		const vector<int> &Get_Sector_Membership() const {return  its_Sector_Membership;}

		int		GetSectorId(int Num)const {return its_Sector_Membership[Num];}

		void	Set_Betas(const vector<double>& value){its_Betas = value;}
		void	Get_Betas(vector<double>& value)const {value = its_Betas;}

		void	Set_Lambdas(const vector<double>& value){its_Lambdas = value;}
		void	Get_Lambdas(vector<double>& value)const{value = its_Lambdas;}

		void	Set_Betas_Down(const vector<double>& value){its_Betas_Down = value;}
		void	Get_Betas_Down(vector<double>& value)const {value = its_Betas_Down;}

		void	Set_Lambdas_Down(const vector<double>& value){its_Lambdas_Down = value;}
		void	Get_Lambdas_Down(vector<double>& value)const {value = its_Lambdas_Down;}


};



#endif