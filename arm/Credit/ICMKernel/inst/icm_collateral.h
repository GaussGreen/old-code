

#if !defined(_ICM_COLLATERAL_H_)
#define _ICM_COLLATERAL_H_


/*********************************************************************************/
/*! \class  ICM_Collateral icm_collateral.h "icm_collateral.h"
 *  \author D. Pouponneau
 *	\version 1.0
 *	\date   June 2003
 *	\brief  <B> Description of a given collateral </B> */
/***********************************************************************************/

/*
	Modif CN : ajout vecteur de signe pour Long Short CDO
	Date : 15/03/05
	Auteur : C. Nicklaus
	Notional > 0  : usual case (PorS given by he tranche description)
	Notional < 0  : short case
	itsNbIssuersShort : Nb Short Issuer in the Collateral
	itsShortIndices Vector of short indices
	itsLongIndices  Vector of Long indices
*/

/*
	Modif CN : ajout d'un coefficient lineaire de recovery
	Date : 25/10/05
	Auteur : C. Nicklaus
	itsRecovCoef : coefficient linéaire appliqué au recovery sur une tranche
	pour l'instant un seul coefficient pour le collateral
*/

/**
		ICM_Collateral 

		This class holds the collateral definition :
		- basket of names (issuer names)
		- associated with notionals values that are reference values
		
		Notionals can be defined at several dates, however it is assumer
		that their sign is not changing 

**/ 

#include "ARMKernel/glob/armglob.h"
#include "ICMKernel/glob/icm_enums.h"
#include "ICMKernel/util/icm_vector.h"
#include "ARMKernel/util/refvalue.h"
#include <string>

class ICM_Collateral : public ARM_Object  
{
private: 
	int								itsNbIssuersShort;				//Number of Short issuers (if = 0 Std)
	std::vector<std::string>		itsIssuersLabels;		//Vector of Labels
	ICM_Vector<ARM_ReferenceValue>	itsIssuersNotionals; 
	double							itsRecovCoef;			//Recovery coefficient
	bool							itsFlgFullHomogeneous;	//for full homogeneous case
	ICM_Collateral*					itsCollateralInDefault; //CollateralForDefaultIssuers
public:

	ICM_Collateral() ;
	ICM_Collateral(const ICM_Collateral&ref) ; 
	virtual ARM_Object* Clone(void) ;
	// ICM_Collateral(const std::vector<std::string>& Labels,
	// 			   const ARM_Vector&  Notionals,
	// 			   const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& 	IsInDefault,
	// 			   const double& RecovCoef = CREDIT_DEFAULT_VALUE) ;
	 
	ICM_Collateral(const std::vector<std::string>& labels,
					const ARM_Vector & Notionals) ;

	// void Set(		const std::vector<std::string>& Labels,
	// 			   const ARM_Vector & Notionals,
	// 			   const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/&	IsInDefault,
	// 			   const double& RecovCoef = CREDIT_DEFAULT_VALUE);
	void Set(const std::vector<std::string>& labels,
					const ARM_Vector&  Notionals ); 
 
	void Set(const std::vector<std::string>& labels,
			const ICM_Vector<ARM_ReferenceValue>& notionals) ;

	virtual ~ICM_Collateral() ;
	
	void Init();

	bool isNotionalConstant() const ; 
	const std::vector<std::string>& GetIssuersLabels() const { return itsIssuersLabels; }
	const std::string& GetIssuersLabels(unsigned int i) const ; 
	const ICM_Vector<ARM_ReferenceValue>& GetIssuersNotionals() const { return itsIssuersNotionals; }
	
	void GetConstantNotionals(ARM_Vector&notionals) const ; 
	void GetNotionals(const ARM_Date&date,ARM_Vector&notionals) const; 
	
	int getIssuerPosition(const std::string&name) const 
	{
		for (int i=0; i<itsIssuersLabels.size(); i++)
			if (itsIssuersLabels[i]==name) return i; 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Collateral::getIssuerPosition: not found "<<name); 
	}
	void GetNotionalDates(ICM_Vector<ARM_Date>&) const; 
	double GetIssuersNotional(unsigned int i) const 
	{
	 	if (!itsIssuersNotionals[i].IsConstant()) 
	 		ICMTHROW(ERR_INVALID_ARGUMENT,"GetIssuersNotional: issuer "<<i<<" is not constant"); 
	 	return itsIssuersNotionals[i].CptReferenceValue(0); 
	}
	double GetIssuersNotional(unsigned int i,const ARM_Date&date) const 
	{
	 	return itsIssuersNotionals[i].CptReferenceValue(date); 
	}
	int GetNbIssuers() const { return itsIssuersLabels.size() ;}
	int GetNbIssuersShort(void) const { return itsNbIssuersShort;}

	void ExcludeIssuer(const std::string&label);

	void SetIssuersInfos(const std::vector<std::string>& IssuersLabels, const ARM_Vector& IssuersNotionals);

	void View(char* id, FILE* ficOut);

	double SumNotionals(const ARM_Date& ) const ;

	bool ContainIssuer(const std::string&label) const ;
	void TransferIssuerInDefaultCollateral(const std::string&label);
						   
	
	bool IsLongShort(void) const
	{ 
		return ((itsIssuersLabels.size()!=itsNbIssuersShort)&&(itsNbIssuersShort));
	}
	
	bool IsFullHomog(void) const  
	{
		return itsFlgFullHomogeneous;
	}

	void GetIssuersInDefault(vector<string>& issuers,const bool& indef = true);

	//Recovery coefficient
	bool GetRecovCoefFlg(void) const 
	{
		return itsRecovCoef !=CREDIT_DEFAULT_VALUE ;
	}
	double GetRecovCoef(void) const				{ return itsRecovCoef;}
	void SetRecovCoef(const double& RecovCoef)	{itsRecovCoef = RecovCoef;}

private:
	
	void CptFullHomog(void) ; 
	void PushIssuer(const std::string& _IssuersLabel,
							const ARM_ReferenceValue& _notional
							);
	void CptInternals(); 
private:
	ICM_Collateral& operator=(const ICM_Collateral&) ;	//NA 
	// void BitwiseCopy(const ARM_Object* src) ;			//NA 
	// void Copy(const ARM_Object* src) ;			//NA 
};

#endif 
 