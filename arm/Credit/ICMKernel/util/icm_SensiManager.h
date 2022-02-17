#error no longer part of the project

#ifndef _ICM_SENSIMANAGER_H_
#define _ICM_SENSIMANAGER_H_

/*********************************************************************************/
/*! \class  ICM_Sensi2D icm_SensiManager.h "icm_SensiManager.h"
 *  \author D Pouponneau
 *	\version 1.0
 *	\date   Jully 2004
 *	\file   icm_SensiManager.h
/***********************************************************************************/

#include "ICMKernel/glob/icm_enums.h"
#include "ICMKernel/util/icm_qmatrix.h"

class ICM_Sensi2D : public ARM_Object
{
private :

qSENSITIVITY_TYPE	itsType;
std::vector<std::string>	itsMktData;
char**			 	itsBumpProfile;

// int					itsMktDataSize;
int					itsBumpProfileSize;

ICM_QMatrix<double>* itsMatrix;
void*				itsObjectMatrix;

public :

ICM_Sensi2D()
{
	Init();
}

~ICM_Sensi2D()
{
	int i=0;
	int j=0;

	// for (i=0; i<itsMktDataSize; i++)
	// {
	// 	if (itsMktData[i])
	// 		delete[] itsMktData[i];
	// 	itsMktData[i] = NULL;
	// }

	// if (itsMktData)
	// 	delete[] itsMktData;
	// itsMktData = NULL;

	for (j=0; j<itsBumpProfileSize; j++)
	{
		if (itsBumpProfile[j])
			delete[] itsBumpProfile[j];
		itsBumpProfile[j] = NULL;
	}

	if (itsBumpProfile)
		delete[] itsBumpProfile;
	itsBumpProfile = NULL;

	if (itsMatrix)
		delete itsMatrix;
	itsMatrix = NULL;
}


ICM_Sensi2D(qSENSITIVITY_TYPE st,const std::vector<std::string>& MktData,char**  BumpProfile,int BumpProfileSize,void* ObjectMatrix)
{
	Init();

	Set(st,MktData,BumpProfile,BumpProfileSize,ObjectMatrix);
}

void Init()
{
	itsType = ICMSPREAD_TYPE;
	// itsMktData = NULL;
	itsBumpProfile = NULL;
	// itsMktDataSize = -999;
	itsBumpProfileSize=-999;
	itsMatrix = NULL;
	itsObjectMatrix = NULL;
}

void Set(qSENSITIVITY_TYPE st,const std::vector<std::string>& MktData,char**  BumpProfile,int BumpProfileSize,void* ObjectMatrix)
{
	int i=0;
	itsType = st;
	// itsMktDataSize =MktDataSize;
	itsBumpProfileSize=BumpProfileSize;

	if (itsMatrix) delete itsMatrix;
	itsMatrix = new ICM_QMatrix<double>(MktData.size(),BumpProfileSize,0.);	

	//if (itsMktData)
	//{
	//for (i=0; i<MktDataSize;i++)
	//	if (itsMktData[i]) delete[] itsMktData[i];
	//}
	//
	//itsMktData = new char*[MktDataSize];

	//for (i=0; i<MktDataSize;i++)
	//{
	//	itsMktData[i] = new char[_size_zclabel_];
	//	strcpy(itsMktData[i],MktData[i]);
	//}
	itsMktData=MktData; 
	if (itsBumpProfile)
	{
	for (i=0; i<BumpProfileSize;i++)
		if (itsBumpProfile[i]) delete[] itsBumpProfile[i];
	}

	itsBumpProfile = new char*[BumpProfileSize];

	for (i=0; i<BumpProfileSize;i++)
	{
		itsBumpProfile[i] = new char[_size_zclabel_];
		strcpy(itsBumpProfile[i],BumpProfile[i]);
	}
	
	if (itsMatrix) delete itsMatrix;
	itsMatrix = new ICM_QMatrix<double>(MktData.size(),BumpProfileSize);

	itsObjectMatrix = ObjectMatrix;
}

qSENSITIVITY_TYPE	GetType() {return itsType;}
const std::vector<std::string>&	GetMktData() {return itsMktData;}
const std::string& GetMktData(unsigned int i) { return itsMktData[i];}
char**			 	GetBumpProfile() {return itsBumpProfile;}
int					GetMktDataSize() { return itsMktData.size();}
int					GetBumpProfileSize() { return itsBumpProfileSize;}
ICM_QMatrix<double>* GetMatrix() {return itsMatrix;}
void*				GetObjectMatrix() {return itsObjectMatrix;}

void				SetObjectMatrix(void* obj) {itsObjectMatrix = obj;}
	// ----------------------------
	//	Copy of members data
	// ----------------------------
	void BitwiseCopy(const ARM_Object* src)
	{
		int i=0;

	    ICM_Sensi2D* sensi2d = (ICM_Sensi2D *) src;
		itsType = sensi2d->itsType;

		// if (sensi2d->itsMktData)
		// {
		// if (itsMktData)
		// {
		// for (i=0; i<itsMktDataSize;i++)
		// 	if (itsMktData[i]) delete[] itsMktData[i];
		// }
		// 
		// itsMktData = new char*[sensi2d->itsMktDataSize];
		// 
		// for (i=0; i<sensi2d->itsMktDataSize;i++)
		// {
		// 	itsMktData[i] = new char[_size_zclabel_];
		// 	strcpy(itsMktData[i],sensi2d->itsMktData[i]);
		// }
		// }
		itsMktData=sensi2d->itsMktData ;

		if (sensi2d->itsBumpProfile)
		{
		if (itsBumpProfile)
		{
		for (i=0; i<itsBumpProfileSize;i++)
			if (itsBumpProfile[i]) delete[] itsBumpProfile[i];
		}

		itsBumpProfile = new char*[sensi2d->itsBumpProfileSize];

		for (i=0; i<sensi2d->itsBumpProfileSize;i++)
		{
			itsBumpProfile[i] = new char[_size_zclabel_];
			strcpy(itsBumpProfile[i],sensi2d->itsBumpProfile[i]);
		}
		}

		// itsMktDataSize=sensi2d->itsMktDataSize;
		itsBumpProfileSize=sensi2d->itsBumpProfileSize;

		if (itsMatrix)
			delete itsMatrix;

		itsMatrix = (ICM_QMatrix<double>*)sensi2d->itsMatrix->Clone();;

		if (sensi2d->itsObjectMatrix)
			itsObjectMatrix = sensi2d->itsObjectMatrix;
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
     ICM_Sensi2D* theClone = new ICM_Sensi2D();

     theClone->Copy(this);
 
     return(theClone);
	}

	double GetValue(const std::string& MktData, const std::string & BumpProfile)
	{
		int i=0,j=0;
		bool status = true;

		for (i=0;i<itsMktData.size();i++)
		{
			if (!strcmp(itsMktData[i].c_str(),MktData.c_str())) 
			{
				status = false;
				break;
			}
		}

		if (status) return NULL;
		status = true;
		for (j=0;j<itsBumpProfileSize;j++)
		{
			if (!strcmp(itsBumpProfile[j],BumpProfile.c_str()))
			{
				status = false;
				break;
			}
		}

		if (status) 
			throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,"ERROR : Unable to find value");

		return itsMatrix->Getvalue(i,j);
	}

	int GetMktDataIndice(const char* MktData)
	{
		int i=0;
		bool status = true;

		for (i=0;i<itsMktData.size();i++)
		{
			if (!strcmp(itsMktData[i].c_str(),MktData)) 
			{
				status = false;
				break;
			}
		}

		if (status) return -1;

		return (i);
	}

	void SetValue(const char* MktData, char* BumpProfile,double value)
	{
		int i=0,j=0;
		bool status = true;

		for (i=0;i<itsMktData.size();i++)
			if (!strcmp(itsMktData[i].c_str(),MktData))
			{
				status = false;
				break;
			}

		for (j=0;j<itsBumpProfileSize;j++)
			if (!strcmp(itsBumpProfile[j],BumpProfile))
			{
				status = false;
				break;
			}

		if (status) 
			throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,"ERROR : Unable to set value");

		itsMatrix->SetValue(i,j,value);
	}

};


/*********************************************************************************/
/*! \class  ICM_SensiManager icm_SensiManager.h "icm_SensiManager.h"
 *  \author D Pouponneau
 *	\version 1.0
 *	\date   Jully 2004
 *	\file   icm_SensiManager.h
/***********************************************************************************/


class ICM_SensiManager : public ARM_Object
{
	private :

	ICM_Sensi2D** itsSensiVector;
	int				itsSensiSize;	

	public :

	ICM_SensiManager()
	{
		Init();
	}

	void Init()
	{
		itsSensiVector = NULL;
		itsSensiSize = 0;
	}	

	~ICM_SensiManager()
	{
		for (int i=0; i<itsSensiSize; i++)
		{
			delete itsSensiVector[i];
			itsSensiVector[i] = NULL;
		}

		if (itsSensiVector)
			delete itsSensiVector;
		itsSensiVector = NULL;
	}	

	ICM_Sensi2D** GetSensiVector() {return itsSensiVector;}

	void push_back(ICM_Sensi2D* sensi2d)
	{
		int i =0;
		bool status = true;

		for (i=0; i<itsSensiSize; i++) //on check si le type n'existe pas déjà
		{
			if (itsSensiVector[i]->GetType() == sensi2d->GetType()) 
			{
				delete itsSensiVector[i];
				itsSensiVector[i] = sensi2d;
				status = false;
				break;
			}
		}
		
		if (status)
		{
		ICM_Sensi2D** SensiVector = new ICM_Sensi2D*[itsSensiSize+1];

		for (i=0; i<itsSensiSize;i++)
			SensiVector[i] = itsSensiVector[i];

		SensiVector[itsSensiSize] = sensi2d;

		if (itsSensiVector)
			delete[] itsSensiVector;
		itsSensiVector = SensiVector;
		}

		itsSensiSize++;
	}

	// ----------------------------
	//	Copy of members data
	// ----------------------------
	void BitwiseCopy(const ARM_Object* src)
	{
		int i=0;

	    ICM_SensiManager* sensiM = (ICM_SensiManager *) src;

		if (sensiM->itsSensiVector)
		{
		if (itsSensiVector)
		{
		for (i=0; i<itsSensiSize;i++)
			if (itsSensiVector[i]) delete itsSensiVector[i];
		}

		itsSensiVector = new ICM_Sensi2D*[sensiM->itsSensiSize];

		for (i=0; i<sensiM->itsSensiSize;i++)
			itsSensiVector[i] = (ICM_Sensi2D*) sensiM->itsSensiVector[i]->Clone();
		}

		itsSensiSize = sensiM->itsSensiSize;
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
     ICM_SensiManager* theClone = new ICM_SensiManager();

     theClone->Copy(this);
 
     return(theClone);
	}

	double GetValue(qSENSITIVITY_TYPE type, const char* MktData, char* BumpProfile)
	{
		int i=0;
		bool status = true;

		for (i=0; i<itsSensiSize; i++)
			if (itsSensiVector[i]->GetType() == type)
			{
				status = false;
				break;
			}

		if (status) 
			throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,"ERROR : Unable to find value");

		itsSensiVector[i]->GetValue(MktData,BumpProfile);	
	}

	void SetValue(qSENSITIVITY_TYPE type, const char* MktData, char* BumpProfile, double value)
	{
		int i=0;
		bool status = true;

		for (i=0; i<itsSensiSize; i++)
			if (itsSensiVector[i]->GetType() == type) 
			{
				status = false;
				break;
			}

		if (status) 
			throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,"ERROR : Unable to set value");

		itsSensiVector[i]->SetValue(MktData,BumpProfile,value);	
	}

	ICM_Sensi2D* GetSensiMatrix(qSENSITIVITY_TYPE type)
	{
		int i=0;
		bool status = true;

		for (i=0; i<itsSensiSize; i++)
		{		
			if (itsSensiVector[i]->GetType() == type)
			{		
				status = false;
				break;
			}
		}

		if (status) return NULL;

		return (itsSensiVector[i]);	
	}

};

#endif