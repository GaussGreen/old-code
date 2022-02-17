#ifndef ARM_RESULT_H
#define ARM_RESULT_H

#ifndef _LOCAL_ARM
	#include "req.h"
#endif


#include <CCString.h>
#include <CCmessage.h>


#ifndef _LOCAL_ARM
	#define ARM_TYPE_DOUBLE		1
	#define ARM_TYPE_STRING		3
	#define ARM_TYPE_UNKNOWN	-1
#endif


class ARM_result
{
public:

	ARM_result ():	objShortName("LANYC"),
					objId(-9999),
					dVal(-9999.),
					sVal("")
					{};

#ifndef _LOCAL_ARM
	ARM_result (ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req);	
#endif	

	~ARM_result () {};


	long getRetCode () const 
	{
		return retCode;
	}
	void setRetCode (long n_retCode)
	{
		retCode = n_retCode;
	}

	CCString getMsg () const
	{
		return msg;
	}
	void setMsg (const CCString& n_msg)
	{
		msg = n_msg;
	}
	void setMsgString(const std::string& n_msg)
	{
		msg = n_msg.c_str(); 
	}
	void setArray(double pVal,long i)
	{
		dArray[i] = pVal;
	}

	double getArray(long i)
	{
		return dArray[i];
	}

#ifndef _LOCAL_ARM
	long getObjectId () const
	{
		return getLong ();
	}

	long getLong () const
	{
		if(value._d () == ARM_CorbaRequestModule::ARM_CORBA_LONG)
		{
			return value.longVal ();
		}

		return 0;
	}

	double getDouble () const
	{
		if(value._d () == ARM_CorbaRequestModule::ARM_CORBA_DOUBLE)
		{
			return value.doubleVal ();
		}

		return 0.0;
	}

	void setDouble (double n_value)
	{
		value.doubleVal (n_value);
	}

	void setString (CCString pVal)
	{
		value.stringVal (pVal.GetStr());
	}

	CCString getString () const
	{
		CCString CCStr;

		if(value._d () == ARM_CorbaRequestModule::ARM_CORBA_STRING)
		{
			CCStr.Set (value.stringVal ());
		}

		return (CCStr);
	}

	int getType () const
	{
		switch(value._d ())
		{
			case ARM_CorbaRequestModule::ARM_CORBA_STRING:
				return ARM_TYPE_STRING;
			case ARM_CorbaRequestModule::ARM_CORBA_DOUBLE:
				return ARM_TYPE_DOUBLE;
			default:
				return ARM_TYPE_UNKNOWN;
		}
	}

	void set (ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req);
#else

	long getLong () const
	{
		return objId;
	}

	void setLong (long id)
	{
		objId = id;
	}

	CCString getShortName() const
	{
		return objShortName;
	}

	void setShortName(const char* shortName )
	{
		objShortName = shortName;
	}


	double getDouble () const
	{
		return dVal;
	}

	void setDouble (double pVal)
	{
		dVal = pVal;
	}

	void setString (CCString pVal)
	{
		sVal = pVal;
	}

	CCString getString () const
	{
		return sVal;
	}

	void setStringInVect (CCString pVal)
	{
		svVal.push_back(pVal);
	}

	const VECTOR<CCString>& getStringVector () const
	{
		return svVal;
	}

#endif

private:
	long retCode;

	double		dArray[1000];

#ifndef _LOCAL_ARM
	ARM_CorbaRequestModule::ARM_CORBA_VALUE value; 
#else
	long		objId;
	double		dVal;
	CCString	sVal;
	VECTOR<CCString> svVal;
#endif

	CCString msg;
	CCString objShortName;
};










#endif	// ARM_RESULT_H

// EOF %M%
