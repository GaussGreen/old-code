#ifndef _ARM_SCALAR_DATA
#define _ARM_SCALAR_DATA

//#include "abstractmarketclass.h"
#include "armval.h"



class ARM_ScalarData : public ARM_AbstractMarketClass
{
	private :
		ARM_Val	itsValue;

	public :
		ARM_ScalarData(void)
		{
			Init();
		}

		ARM_ScalarData(double aData)
		{
			Init();
			SetValue(aData);
		}

		ARM_ScalarData(ARM_ScalarData& aData) : ARM_AbstractMarketClass(aData)
        {
            Init();
            BitwiseCopy(&aData);
        }

		~ARM_ScalarData(void)
		{
		}

        ARM_ScalarData& operator = (ARM_ScalarData& srcVal)
        {
            (*this).ARM_AbstractMarketClass::operator = (srcVal);
 
            BitwiseCopy(&srcVal);
 
            return	(*this);
        }

        int operator == (ARM_ScalarData& srcVal);

        int operator < (ARM_ScalarData& srcVal);

		void Init(void)
		{
			SetName(ARM_SCALAR_DATA);
			itsValue.SetDoubleVal(0.);
		}

		void BitwiseCopy(const ARM_Object* aData)
        {
			ARM_ScalarData*	vData = (ARM_ScalarData*)aData;
			itsValue = vData->itsValue;
		}

        void Copy(const ARM_Object* srcVal)
        {
            ARM_AbstractMarketClass::Copy(srcVal);
 
            BitwiseCopy(srcVal);
        }
 
        ARM_Object* Clone(void)
        {
            ARM_ScalarData* theClone = new ARM_ScalarData();
 
            theClone->Copy(this);
 
            return	theClone;
        }
	
		double GetValue(void)
		{
			return	itsValue.GetValue().doubleVal;
		}

		void SetValue(double aData)
		{
			itsValue.SetDoubleVal(aData);
		}

		void View(char* id = NULL, FILE* ficOut = NULL);
};


#endif