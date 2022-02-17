#ifndef _ARM_LOCAL_PERSISTENT_H
#define _ARM_LOCAL_PERSISTENT_H


#include <glob\armglob.h>
#include <glob\expt.h>
#include <ICMKernel/util/icm_macro.h>
#include <time.h>

//class ARM_Object;

#define LOCAL_MAX_OBJECTS_NUMBER    10000000

#define ARM_OK					0
#define ARM_KO					-1


class LocalPersistent ;
extern LocalPersistent* LOCAL_PERSISTENT_OBJECTS;

class LocalPersistent
{
public:
	//	JLA
	//
	//		Singleton - like pattern
	//		Use get to access the instance
	//
	//		The only difference is that the instance is not static member of the class,
	//		but static pointer 
	// 
	static LocalPersistent& get() 
	{
		if (!LOCAL_PERSISTENT_OBJECTS) 
			LOCAL_PERSISTENT_OBJECTS=new LocalPersistent(); 
		if (!LOCAL_PERSISTENT_OBJECTS) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"LocalPersistent::LocalPersistent: can't create !"); 
		return *LOCAL_PERSISTENT_OBJECTS; 
	}
public:

	// Array of objects
 
	ARM_Object**    LOCAL_OBJECTS_ARRAY;            
    int             index;
 

	LocalPersistent(void)
	{
		LOCAL_OBJECTS_ARRAY = new ARM_Object*[LOCAL_MAX_OBJECTS_NUMBER];

		index = -1;

		for ( int i = 0; i < LOCAL_MAX_OBJECTS_NUMBER; i++)
		{
			LOCAL_OBJECTS_ARRAY[i] = (ARM_Object *) NULL;
		}

	}

	~LocalPersistent(void)
	{
		try
		{
			for ( int i = 0; i <= index; i++)
			{
				if (LOCAL_OBJECTS_ARRAY[i])
				   delete LOCAL_OBJECTS_ARRAY[i];
        
				LOCAL_OBJECTS_ARRAY[i] = NULL;
			}
			//JLA Added: 
			delete [] LOCAL_OBJECTS_ARRAY ;

			LOCAL_OBJECTS_ARRAY = NULL;
		}
		catch(...)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			 "Persistency Pb : Pb in destructor : Unknown Object");
		}
	}

	void FreeAllObjects(void)
	{
		try
		{
			for ( int i = 0; i <= index; i++)
			{
				if (LOCAL_OBJECTS_ARRAY[i])
				   delete LOCAL_OBJECTS_ARRAY[i];
        
				LOCAL_OBJECTS_ARRAY[i] = NULL;
			}

			index = -1;
		}
		catch(...)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			 "Persistency Pb : Pb in FreeAllObjects : Unknown Object");
		}
	}

	long GetObjectIdFromPointer(ARM_Object* obj)
	{
		try
		{
			for ( int i = 0; i <= index; i++)
			{
				if ( LOCAL_OBJECTS_ARRAY[i] == obj )
				{
				   return(i);
				}
			}

			return(ARM_NULL_OBJECT);
		}
		catch(...)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			 "Persistency Pb : Unknown Object Ptr");

			return NULL;
		}
	}

	/// basically rename the function to avoid conflict with SDK function GetObject!
	inline ARM_Object* GetARMObject(long objectId)
	{
		return GetObject(objectId);
	}

	ARM_Object* GetObject(long objectId)
	{
		try
		{
			if (( objectId > index) 
				||
				( objectId < 0 )
				||
				( objectId >= LOCAL_MAX_OBJECTS_NUMBER )
			   )
			{
			   return((ARM_Object *) NULL);
			}

			return(LOCAL_OBJECTS_ARRAY[objectId]);
		}
		catch(...)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			 "Persistency Pb : Unknown Object");

			return NULL;
		}
	}


	//	JLA -
	//
	//		This methods insert the object into the cache; 
	//		id=-9999	:	this is new object
	//		id<0		:	this is an error
	//		id>0		:	previous object has beel already deleted...
	//
	//	It seems that null objects can be inserted, however it might be
	//	that null object are filtered before. 
	//
	int SetPersistent(ARM_Object* obj, long id=-9999)
	{
		if ( id == -9999 )
		{
			if ( index >= LOCAL_MAX_OBJECTS_NUMBER-1 )
			{
				return(RET_KO);
			}

			++index;

			LOCAL_OBJECTS_ARRAY[index] = (ARM_Object *) obj;

			return(index);
		}

		if ( id < 0 )
		{
		   return(RET_KO);
		}

		if ( id > index )
		{
			return(RET_KO);
		}
		else
		{
			LOCAL_OBJECTS_ARRAY[id] = (ARM_Object *) obj;

			return(id);
		}
	}

	//	JLA
	//
	//	If this function is called with a NEGATIVE ID,
	//		then it will returns the new id for the object.
	//		
	//	If it is call with POSITIVE ID,
	//		then it will remove the previous object and
	//		inserts the new one. 
	//
	//	Advantages; 
	//		Function is not static, hence can only be called to an allocated cache
	//		Function automatically remove previous object : caller code is simplified 
	//
	long adopt(ARM_Object* obj, long id)
	{
		if (!obj) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"LocalPersistent::adopt: can't insert null object") ; 
		if (id>index) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"LocalPersistent::adopt: can't find #"<<id<<" max:"<<index); 
		if (id<0)
		{
			if ( index >= LOCAL_MAX_OBJECTS_NUMBER-1 )
				ICMTHROW(ERR_INVALID_ARGUMENT,"LocalPersistent::adopt: can't add any more object, max="<<LOCAL_MAX_OBJECTS_NUMBER); 
			++index ;
			LOCAL_OBJECTS_ARRAY[index] = obj; 
			return(index);
		}
		if (LOCAL_OBJECTS_ARRAY[id]) delete  LOCAL_OBJECTS_ARRAY[id] ;
		LOCAL_OBJECTS_ARRAY[id] = obj;
		return id; 
	}
	// JLA 
	//
	//	This function attemps to convert an ARM_Object to 
	//		a concrete type. In any case of failure an exception
	//		will be returned
	//
	template <class T>  void convert(long id,T*& ret)
	{
		LocalPersistent::get(); 
		if (id==-1L) 
		{ 
			ret=0; return ; 
		}
		if ( (id<0) || (id>index)) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"LocalPersistent::convert: wrong id "<<id) ;
		ARM_Object* item = LOCAL_OBJECTS_ARRAY[id] ; 
		if (!item) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"LocalPersistent::convert: null object "<<id) ;
		T* tmp = dynamic_cast<T*>(item); 
		if (!tmp) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"LocalPersistent::convert: can't convert "
				<<typeid(item).name()<<" to " << typeid(T).name() ) ; 
		ret=tmp; 
	}
	//	JLA
	//
	//	This function will parse the string an return the id
	//	Equivalent to LocalGetNumObjectId, but with strings. 
	// 
	long getObjectId (const std::string& stringId)
	{
		/// we suppose that any ARM Object start by L
		/// and is followed by 4 capital letters or numbers characters followed by '_'
		/// if this is not the case, return an error
		if (stringId.size()<6) return ARM_KO; 
		if ( stringId[0] != 'L' ) return ARM_KO; 
		for ( int i = 1; i < 5; ++i)
		{
			if ((( stringId[i] < 'A' ) || ( stringId[i] > 'Z' ))
				&&  
				(( stringId[i] < '0' ) || ( stringId[i] > '9' ))
			   ) return ARM_KO;
		}
		if ( stringId[5] != '_' ) return ARM_KO ;

		// now we extract the desired word... 
		std::string tmp = stringId.substr(6,stringId.find(' ',6)); 
		if (tmp.empty()) return ARM_KO; 
		return atol(tmp.c_str());
	}
	std::string getStringId(long objId,const std::string objClass)
	{
		if (objId<0) throw std::runtime_error("can't create stringId"); 
		// TMP: objId a tester
		char buf[50];
		char sHour[4], sMin[4], sSec[4];

		struct tm *newtime;
		time_t long_time;

		time( &long_time );                /* Get time as long integer. */
		newtime = localtime( &long_time ); /* Convert to local time. */
		
		if (newtime->tm_hour < 10)
			sprintf (sHour, "0%1d", newtime->tm_hour);
		else
			sprintf (sHour, "%2d", newtime->tm_hour);
		
		if (newtime->tm_min < 10)
			sprintf (sMin, "0%1d", newtime->tm_min);
		else
			sprintf (sMin, "%2d", newtime->tm_min);
		
		if (newtime->tm_sec < 10)
			sprintf (sSec, "0%1d", newtime->tm_sec);
		else
			sprintf (sSec, "%2d", newtime->tm_sec);

		sprintf (buf, "%ld %s:%s:%s", objId,sHour,sMin,sSec);
		std::string ret = objClass + "_" + buf ;
		return (ret);
	}

	void FreeObject(long id)
	{
		if (( id > index )	||
			( id < 0  ))
		{
			return;
		}

		if (LOCAL_OBJECTS_ARRAY[id])
		{
			delete LOCAL_OBJECTS_ARRAY[id];

			LOCAL_OBJECTS_ARRAY[id] = (ARM_Object *) NULL;
		}
	}

	static int LOCAL_IS_OBJECT_CLASS_OK(ARM_Object* obj, ARM_CLASS_NAME oClassName)
	{
		if ( obj == NULL )
		{
			return 0;
		} 
		if ( obj->GetName() != oClassName )
		{
			return 0;
		}
		return 1;
	}

	static int LOCAL_IS_OBJECT_ROOT_CLASS_OK(ARM_Object* obj, ARM_CLASS_NAME oClassName)
	{
		if ( obj == NULL )
		{ 
			return 0;
		}
		if ( obj->GetRootName() != oClassName )
		{
			return 0;
		}
		return 1;
	}

	static int LOCAL_IS_SWAPLEG_OK(ARM_Object* swapLeg)
	{
		if ( swapLeg == NULL )
		{
		   return 0;
		}

		switch(swapLeg->GetName())
		{
			case ARM_INFLEG:
			case ARM_FIXZC:
			case ARM_SWAPLEG:
			case ARM_FIXEDLEG:
			case ARM_SPREADLEG:
			case ARM_CMSLEG:
			case ARM_CMTLEG:
			case ARM_TMLEG:
			case ARM_CORRIDORLEG:
			{
				return 1;
				break;
			}
			default:
			{
				return 0;
			}
		}
	}


	// For Credit ************************************************

	/** 
	static int LOCAL_IS_CREDIT_OK(ARM_Object* swapLeg)
	{
		if ( swapLeg == NULL )
		{
		   return 0;
		}

		switch(swapLeg->GetName())
		{
			case ICM_CDS:
			case ICM_BOND:
			case ICM_LEG:
			case ICM_FRN:
			case ICM_CLN:
			case ICM_ZCCLN:
			{
				return 1;
				break;
			}
			default:
			{
				return 0;
			}
		}
	}

	
	static int LOCAL_IS_CREDIT_CDSLIKE_OK(ARM_Object* swapLeg)
	{
		if ( swapLeg == NULL )
		{
		   return 0;
		}

		switch(swapLeg->GetName())
		{
			case ICM_CDS:
			case ICM_FTD:
			case ICM_MEZ:
			case ICM_NTD:
			{
				return 1;
				break;
			}
			default:
			{
				return 0;
			}
		}
	}
	 
	static int LOCAL_IS_CREDIT(ARM_Object* swapLeg)
	{
		if ( swapLeg == NULL )
		{
		   return 0;
		}

		switch(swapLeg->GetName())
		{
			case ICM_CDS:
			case ICM_LEG:
			case ICM_BOND:
			{
				return 1;
				break;
			}
			default:
			{
				return 0;
			}
		}
	}
	**/ 
	// *******************************************************************

	static int LOCAL_IS_CURVE_CLASS_OK(ARM_Object* obj)
	{
		if ( obj == NULL )
		{ 
			return 0;
		}
		if ( ( obj->GetRootName() == ARM_ZERO_CURVE )
			|| ( obj->GetRootName() == ARM_VOL_CURVE ) )
		{
			return 1;
		}
		return 0;
	}

};

inline int CHECK_GLOBAL_OBJECT(LocalPersistent* myPersistentObjects)
{
	int result (1);

	if ( myPersistentObjects == (LocalPersistent *) NULL )
	{
		result = 0;
	}
	return result;
}

#define CREATE_GLOBAL_OBJECT() {  LocalPersistent::get(); }

/** #define CREATE_GLOBAL_OBJECT() \
{ \
	if ( LOCAL_PERSISTENT_OBJECTS == (LocalPersistent *) NULL ) \
	{ \
		LOCAL_PERSISTENT_OBJECTS = new LocalPersistent(); \
	} \
}**/ 




#endif // _ARM_LOCAL_PERSISTENT_H
