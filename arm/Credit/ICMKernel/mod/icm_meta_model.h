
#if !defined(_ICM_META_MODEL_H_)
#define _ICM_META_MODEL_H_

#include "ARMKernel\mod\model.h"
#include "ICMKernel\util\icm_matrix.h"


/*********************************************************************************/
/*! \class  ICM_Meta_Model icm_meta_model.h "icm_meta_model.h"
 *  \author D. Pouponneau
 *	\version 1.0
 *	\date   Febuary 2005
 *	\brief  Description of a vector of models */
/***********************************************************************************/


class ICM_Meta_Model : public ARM_Model
{
private: 

	int					    itsNbModel;
	vector<ARM_CLASS_NAME>	itsPricerType;
	ARM_Model**				itsModels;  

private: 
	
	void Init();

public:


	ICM_Meta_Model()
	{
		Init();
	}	


	ICM_Meta_Model(ARM_Model** Models,
				  const int& nbofmodels,
				  const vector<ARM_CLASS_NAME>& PricerType) : ARM_Model()
	{
		Init();

		Set(Models, nbofmodels,PricerType);
	}

	void Set(ARM_Model** Models,
				  const int& nbofmodels,
				  const vector<ARM_CLASS_NAME>& PricerType);


	virtual ~ICM_Meta_Model()
	{

		if (itsNbModel)
		{
			delete[] itsModels;
			itsModels = NULL;
		}

		itsPricerType.clear();

		itsNbModel = 0;
	}	

	void BitwiseCopy(const ARM_Object* srcleg);

	void Copy(const ARM_Object* srcleg);

	void SetPricerType(const vector<ARM_CLASS_NAME>& pricertype)
	{
		itsPricerType.clear();

		itsPricerType = pricertype;
	}

	void GetPricerType(vector<ARM_CLASS_NAME>& pricertype) { pricertype = itsPricerType;}

	ARM_Object* Clone(void);

	ARM_Model** GetModels()
	{
		return itsModels;
	}

	ARM_Model* GetModel(const int& n)
	{
		if (n>=0 && n<itsNbModel)
			return itsModels[n];
		else
			return NULL;
	}

	int GetNbModel() {return itsNbModel;};

	void View(char* id, FILE* ficOut)
	{
		FILE* fOut;
		char  fOutName[200];

		if ( ficOut == NULL )
		{
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);

		(void) unlink(fOutName);

		fOut = fopen(fOutName, "w"); 
		}
		else
		{
		fOut = ficOut;
		} 

		for (int i=0; i<itsNbModel ; i++)
		{
		// Affichage de la matrice
		fprintf(fOut, "\t\t\t ----------------- Model N°%ld ----------------- \n",i);

		itsModels[i]->View(id, fOut);
		}

		if ( ficOut == NULL )
		{
		fclose(fOut);
		}
	}

};


#endif 
