//////////////////////////////////////////////////////////////////////
// ICM_RangeFactor.h: interface for the ICM_RangeFactor class.
//
//////////////////////////////////////////////////////////////////////


#if !defined(_ICM_RANGE_FACTOR_H_)
#define _ICM_RANGE_FACTOR_H_

/********************************************************************************/
/*! \file icm_RangeFactor.h
 * 
 *  \brief Describes a Range Factor object
 *  \author 
 *	\version 1.0
 *	\date   April 2005 */
/*
 *********************************************************************************/
// A chaque intervalle définit par une borne sup et une borne inf, on fait correspondre une valeur
// il faut que les vecteurs Inf, Max et Value soient de même dimension

#include "ARMKernel\glob\linalg.h"
#include "ICMKernel/util./icm_macro.h"
#include <vector>
// 17783 using namespace std;

class ICM_RangeFactor: public ARM_Object
{
private:
	unsigned int itsnbRow;			// nb de ligne
	unsigned int itsIndice;			// Situe Intervalle dans lequel on se trouve
	std::vector<double> itsVInf;			// Vector des bornes Inférieurs des intervalles	
	std::vector<double> itsVMax;			// Vector des bornes Superieurs des intervalles
	std::vector<double> itsVValueMin;	// Valeurs Min correspondantes
	std::vector<double> itsVValueMax;	// Valeurs Max correspondantes
	std::vector<double> itsVValueMid;	// Valeurs Mid correspondantes

public :
	void Init();
	
	ICM_RangeFactor() { Init(); }
	
	ICM_RangeFactor(const std::vector<double>& inf,
					const std::vector<double>& max,
					const std::vector<double>& valuemin,
					const std::vector<double>& valuemax)
	{
		Init();
		Set(inf,max,valuemin,valuemax);
	}

	void Set(const std::vector<double>& inf,
			 const std::vector<double>& max,
			 const std::vector<double>& valuemin,
			 const std::vector<double>& valuemax);

	~ICM_RangeFactor()
	{
		itsVInf.clear();
		itsVMax.clear();
		itsVValueMin.clear();
		itsVValueMax.clear();
		itsVValueMid.clear();
	};

	// ----------------------------
	//	Get et Set
	// ----------------------------
	int GetNbRow()		{ return itsnbRow;}
	int GetIndice()		{ return itsIndice;}
	std::vector<double> GetInf()		 { return itsVInf;}
	std::vector<double> GetMax()		 { return itsVMax;}
	std::vector<double> GetValueMin() { return itsVValueMin;}
	std::vector<double> GetValueMax() { return itsVValueMax;}
	std::vector<double> GetValueMid() { return itsVValueMid;}
	
	void SetInf(const std::vector<double>& inf) 	{ itsVInf = inf;}
	void SetMax(const std::vector<double>& max)	{ itsVMax = max;}
	void SetValueMin(const std::vector<double>& ValueMin) 	{ itsVValueMin = ValueMin;}
	void SetValueMax(const std::vector<double>& ValueMax) 	{ itsVValueMax = ValueMax;}
	void SetIndice(const unsigned int& indice)
	{
		if (indice >= itsnbRow)
		{
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_RangeFactor Set Indice: depassement d'indice");	
		}
		else
			itsIndice = indice;
	}
	
	// ----------------------------
	//	Copy of members data
	// ----------------------------
	void BitwiseCopy(const ARM_Object* srcRange)
	{
			ICM_RangeFactor* RangeFactor = (ICM_RangeFactor*) srcRange; 
	
			itsnbRow = RangeFactor->itsnbRow;
			itsIndice = RangeFactor->itsIndice;
			itsVInf = RangeFactor->itsVInf;
			itsVMax = RangeFactor->itsVMax;
			itsVValueMin = RangeFactor->itsVValueMin;
			itsVValueMax = RangeFactor->itsVValueMax;
			itsVValueMid = RangeFactor->itsVValueMid;
	}

	// -------------
	//	Copy Method 
	// -------------
	void Copy(const ARM_Object* src)
	{
		ARM_Object::Copy(src) ;
 
		BitwiseCopy(src) ;
	}

	// --------------
	//	Clone Method
	// --------------
	ARM_Object* Clone(void)
	{
	   ICM_RangeFactor* theClone = new ICM_RangeFactor();

	   theClone->Copy(this);

	   return(theClone);
	}

	// --------------
	//	View Method
	// --------------
	void View(char* id, FILE* ficOut);

	// --------------
	//	Other Method
	// --------------
	void computeMid(void);
	int FindIndex(double x);
};
#endif