#pragma warning(disable : 4541)

#include "firstToBeIncluded.h"
#include "ARM_local_persistent.h"
#include "ARM_local_glob.h"

#include <ARM\local_xlarm\ARM_local_interglob.h>

#include <util\exercise.h>




long ARMLOCAL_BERMUDANXSTYLE(VECTOR<double>& xDates,
                             VECTOR<double>& expiryDates,
							 ARM_result& result,
							 long objId)
{
	long xstyleId;

	long size = xDates.size();


	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 )
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");

		return(ARM_KO);
	}

	ARM_ExerciseStyle* xStyle = NULL;
	ARM_Vector*        vDates = NULL;
   
	VECTOR<CCString> date_str;

	char  sDate[30];


	CCString msg("");


	try
	{
		vDates = new ARM_Vector(size);

        ARM_Vector expiryJulDates;


		for (int i = 0; i < size; i++)
		{
			Local_XLDATE2ARMDATE(xDates[i], sDate);
			
            vDates->Elt(i) = ((ARM_Date) sDate).GetJulian();
		}

        int expirySize = expiryDates.size();

        if ( expirySize != 0 )
        {
           // In fact : expirySize == size
           expiryJulDates.Resize(expirySize);

           for (int i = 0; i < expirySize; i++)
           {
			   Local_XLDATE2ARMDATE(expiryDates[i], sDate);
			
               expiryJulDates.Elt(i) = ((ARM_Date) sDate).GetJulian();
           }
        }

		if ( objId == -1 )
		{
			xStyle = new ARM_ExerciseStyle(vDates);

			xStyle->SetExerciseType(K_BERMUDAN);

            if ( expirySize != 0 )
            {
               xStyle->SetExerciseEndDates(&expiryJulDates);
			   
               xStyle->SetExerciseDatesSize(expiryJulDates.GetSize());
            }

			if (vDates)
			   delete vDates;
			vDates = NULL;

			if ( xStyle == NULL )
			{
			   result.setMsg("ARM_ERR: XStyle is null!");
				
               return(ARM_KO);
			}

			CREATE_GLOBAL_OBJECT();
 
			xstyleId = LOCAL_PERSISTENT_OBJECTS->SetPersistent(xStyle);

			if ( xstyleId == RET_KO )
			{
				if (xStyle)
				   delete xStyle;
				xStyle = NULL;

				result.setMsg("ARM_ERR: Pb with inserting object");				
				
                return(ARM_KO);
			}
	
            result.setLong(xstyleId);

			return(ARM_OK);
		}
		else
		{
			xStyle = (ARM_ExerciseStyle *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);
 
			if ( LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(xStyle, ARM_EXERCISE_STYLE) == 1 )
			{
				xStyle->Set(vDates);
			
                if ( expirySize != 0 )
                {
                   xStyle->SetExerciseEndDates(&expiryJulDates);
			   
                   xStyle->SetExerciseDatesSize(expiryJulDates.GetSize());
                }

				xStyle->SetExerciseType(K_BERMUDAN);

				if (vDates)
				   delete vDates;
				vDates = NULL;

				result.setLong(objId);

				return(ARM_OK);
			}
			else
			{
				result.setMsg("ARM_ERR: previous object is not of a good type");
				
                return(ARM_KO);
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (vDates)
		   delete vDates;
		vDates = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_EUROPEANXSTYLE(double xdate,
							 ARM_result& result,
							 long objId)
{
	long xstyleId;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	ARM_ExerciseStyle* xStyle = NULL;
   
	char* sDate = new char[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE (xdate,sDate);

		if (objId == -1)
		{
			xStyle = new ARM_ExerciseStyle((ARM_Date) sDate);
 
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			if (xStyle == NULL)
			{
				result.setMsg ("ARM_ERR: XStyle is null");
				return ARM_KO;
			}

			CREATE_GLOBAL_OBJECT();
 
			xstyleId = LOCAL_PERSISTENT_OBJECTS->SetPersistent(xStyle);

			if (xstyleId == RET_KO)
			{
				if (xStyle)
					delete xStyle;
				xStyle = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}
			result.setLong(xstyleId);

			return ARM_OK;
		}
		else
		{
			xStyle = (ARM_ExerciseStyle *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);
 
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(xStyle, ARM_EXERCISE_STYLE) == 1)
			{
				xStyle->Set((ARM_Date) sDate);
				
				if (sDate)
					delete [] sDate;
				sDate = NULL;

				result.setLong(objId);

				return ARM_OK;
			}
			else
			{
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (sDate)
			delete [] sDate;
		sDate = NULL;

		ARM_RESULT();
	}
}

long ARMLOCAL_AMERICANXSTYLE (double xStartDate,
							  double xEndDate,
							  ARM_result& result,
							  long objId)
{
	long xstyleId;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	ARM_ExerciseStyle* xStyle = NULL;
   
	char* sStartDate = new char[11];
	char* sEndDate = new char[11];

	CCString msg ("");

	try
	{
		Local_XLDATE2ARMDATE (xStartDate,sStartDate);
		Local_XLDATE2ARMDATE (xEndDate,sEndDate);

		if (objId == -1)
		{
			xStyle = new ARM_ExerciseStyle((ARM_Date) sStartDate, (ARM_Date) sEndDate);
 
			if (sStartDate)
				delete [] sStartDate;
			sStartDate = NULL;

			if (sEndDate)
				delete [] sEndDate;
			sEndDate = NULL;

			if (xStyle == NULL)
			{
				result.setMsg ("ARM_ERR: XStyle is null");
				return ARM_KO;
			}

			CREATE_GLOBAL_OBJECT();
 
			xstyleId = LOCAL_PERSISTENT_OBJECTS->SetPersistent(xStyle);

			if (xstyleId == RET_KO)
			{
				if (xStyle)
					delete xStyle;
				xStyle = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}
			result.setLong(xstyleId);

			return ARM_OK;
		}
		else
		{
			xStyle = (ARM_ExerciseStyle *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);
 
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(xStyle, ARM_EXERCISE_STYLE) == 1)
			{
				xStyle->Set((ARM_Date) sStartDate, (ARM_Date) sEndDate);

				if (sStartDate)
					delete [] sStartDate;
				sStartDate = NULL;

				if (sEndDate)
					delete [] sEndDate;
				sEndDate = NULL;

				result.setLong(objId);

				return ARM_OK;
			}
			else
			{
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (sStartDate)
			delete [] sStartDate;
		sStartDate = NULL;

		if (sEndDate)
			delete [] sEndDate;
		sEndDate = NULL;

		ARM_RESULT();
	}
}

