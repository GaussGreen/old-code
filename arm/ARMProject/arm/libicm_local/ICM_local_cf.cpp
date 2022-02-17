
#pragma warning(disable :4786)

#include "firstToBeIncluded.h"
#include <ARM\libarm_local\ARM_local_persistent.h>
#include <ARM\libarm_local\ARM_local_glob.h>

#include <ARMKernel\inst\security.h>
#include <ICMKernel\inst\icm_gen.h>
#include <ICMKernel\util\icm_matrix.h>

#include <ARM\local_xlarm\ARM_local_interglob.h>




long ICMLOCAL_CashFlows (const VECTOR<CCString>& matrice,
						 int nbrows,
						 int nbcolumns,
						 ARM_result& result,
						 long objId)
{
	long genId = 0;

	ICM_GenCF* prev = NULL;
	ICM_GenCF* newp = NULL;
	ICM_Matrix<ARM_Vector>* ICMMATRIX = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	char* sDate = new char[11];
	const char* tmp = NULL;
	char* tmp2 = NULL;
	char* tmp3 = NULL;
	int j =0,active = 0;
	double value = 0.;

	try
	{

		ICMMATRIX = new ICM_Matrix<ARM_Vector>();
		
		for (int i = 0; i<nbcolumns; i++)
		{
			CCString ColName = (const char*) matrice[i]; 
			ColName.trim_right();
			ColName.toUpper();

			tmp = (const char*) ColName; 
			tmp2 = (char*) ColName; 
			tmp2[2]='\0';


			ARM_Vector* vector = new ARM_Vector(nbrows-1,0.);

			for (j = 1; j<nbrows; j++) 
			{
				active = i + j*nbcolumns;

				if ((strcmp(tmp2,"D_") == NULL) || (strcmp(tmp,"SCHEDULE") == NULL))
				{
				tmp3 = (char*) matrice[active]; 
				value = atof(tmp3);
				Local_XLDATE2ARMDATE(value,sDate); 
				ARM_Date date = (ARM_Date) sDate;
				vector->InitElt(j-1,date.GetJulian());
				if (tmp3)
					delete[] tmp3;
				}
				else
				{
				tmp3 = (char*) matrice[active]; 
				value = atof(tmp3);
				vector->InitElt(j-1,value);
				if (tmp3)
					delete[] tmp3;
				}

			}

			ICMMATRIX->Push(vector, (char*)tmp);
		}

		newp = new ICM_GenCF(ICMMATRIX);

		if (ICMMATRIX)
			delete ICMMATRIX;
		ICMMATRIX = NULL;

		if (sDate)
			delete sDate;
		sDate = NULL;


		if (newp == NULL)
		{
			result.setMsg ("ARM_ERR: Object Generic Cash Flows is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			genId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ICM_GenCF*)newp);

			if (genId == RET_KO)
			{
				if (newp)
					delete newp;
				newp = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(genId);

			return ARM_OK;
		}
		else
		{
			ARM_Object* CF2 = LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			prev = (ICM_GenCF*) (CF2);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(CF2, ICM_GENCF) == 1)
			{
				if (prev)
				{
					delete prev;
					prev = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ICM_GenCF*)(newp), objId);

				return ARM_OK;
			}
			else
			{
				if (newp)
					delete newp;
				newp = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (newp)
			delete newp;
		newp = NULL;

		if (sDate)
			delete sDate;
		sDate = NULL;

		ARM_RESULT();
	}
}


long ICMLOCAL_Parameters (const VECTOR<CCString>& matrice,
						  int nbrows,
						  int nbcolumns,
						  ARM_result& result,
						  long objId)
{
	long genId = 0;

	ICM_Parameters* prev = NULL;
	ICM_Parameters* newp = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	char* sDate = new char[11];
	const char* tmp = NULL;
	char* tmp2 = NULL;
	char* tmp3 = NULL;
	int j =0,active = 0;

	try
	{

		newp = new ICM_Parameters();
		
		for (int i = 0; i<nbcolumns; i++)
		{
			CCString ColName = (const char*) matrice[i]; 
			ColName.trim_right();
			ColName.toUpper();

			tmp = (const char*) ColName; 
			tmp2 = (char*) ColName; 
			tmp2[4]='\0';

			if ((strcmp(tmp2,"STR_") != NULL))
			{	
			ARM_Vector* vector = new ARM_Vector(nbrows-1,0.);

			for (j = 1; j<nbrows; j++) 
			{
				double value = 0.;
				active = i + j*nbcolumns;
				tmp3 = (char*) matrice[active]; 
				value = atof(tmp3);
				vector->InitElt(j-1,value);
				if (tmp3) delete[] tmp3;
			}

			newp->Push(vector, (char*)tmp);
			}
			else
			{	
			ICM_QMatrix<string>* vector = new ICM_QMatrix<string>(nbrows-1,1);

			for (j = 1; j<nbrows; j++) 
			{
				string value;
				active = i + j*nbcolumns;
				tmp3 = (char*) matrice[active]; 
				value = (string)(const char*)(tmp3);
				(*vector)(j-1,0)=value;
				if (tmp3) delete[] tmp3;
			}

			newp->Push(vector, (char*)tmp);
			}

		}

		if (newp == NULL)
		{
			result.setMsg ("ARM_ERR: Object Generic Cash Flows is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			genId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ICM_Parameters*)newp);

			if (genId == RET_KO)
			{
				if (newp)
					delete newp;
				newp = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(genId);

			return ARM_OK;
		}
		else
		{
			ARM_Object* CF2 = LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			prev = (ICM_Parameters*) (CF2);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(CF2, ICM_PARAMETERS) == 1)
			{
				if (prev)
				{
					delete prev;
					prev = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ICM_GenCF*)(newp), objId);

				return ARM_OK;
			}
			else
			{
				if (newp)
					delete newp;
				newp = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (newp)
			delete newp;
		newp = NULL;

		ARM_RESULT();
	}
}
