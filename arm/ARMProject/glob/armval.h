/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : armval.h                                                     */
/*                                                                            */
/* DESCRIPTION : global utilities                                             */
/*                                                                            */
/* DATE        : Thu Feb 1998                                                 */
/*                                                                            */
/*----------------------------------------------------------------------------*/
#ifndef _ARMVAL_H
#define _ARMVAL_H


#include <string.h>

#include "armglob.h"

#include "retcode.h"



class ARM_Val : public ARM_Object
{
    private:

             ARM_VALUE_TYPE valType;

             ARM_VALUE      val;

    public:


        ARM_Val(void)
        {
            SetName(ARM_VAL);

            valType = ARM_DOUBLE;

            val.doubleVal = 0.0;
        }

        ARM_Val(int v)
        {
            SetName(ARM_VAL);

            SetIntVal(v);
        }

        ARM_Val(double v)
        {
            SetName(ARM_VAL);
 
            SetDoubleVal(v);
        }

        ARM_Val(char* v)
        {
            SetName(ARM_VAL);
 
            SetStringVal(v);
        }

        ARM_Val(ARM_Val& srcVal) : ARM_Object(srcVal)
        {
            SetName(ARM_VAL);
 
            BitwiseCopy(&srcVal);
        }

       ~ARM_Val(void)
        {
        }

        ARM_Val& operator = (ARM_Val& srcVal)
        {
            (*this).ARM_Object::operator = (srcVal);
 
            BitwiseCopy(&srcVal);
 
            return(*this);
        }

        int operator == (ARM_Val& srcVal);

        int operator < (ARM_Val& srcVal);

        void BitwiseCopy(const ARM_Object* srcVal)
        {
            ARM_Val* value = (ARM_Val *) srcVal;


            valType = value->valType;

            switch(value->valType)
            {
                case ARM_INT :
                {
                    val.intVal = value->val.intVal; 
                };
                break;

                case ARM_DOUBLE :
                {
                    val.doubleVal = value->val.doubleVal; 
                };
                break;
            
                case ARM_STRING:
                {
                    strcpy(val.stringVal, value->val.stringVal);
                };
                break;

                default : 
                break;
            };
        }

        void Copy(const ARM_Object* srcVal)
        {
            ARM_Object::Copy(srcVal);
 
            BitwiseCopy(srcVal);
        }
 
        ARM_Object* Clone(void)
        {
            ARM_Val* theClone = new ARM_Val();
 
 
            theClone->Copy(this);
 
            return(theClone);
        }

        ARM_VALUE GetValue(void)
        {
            return(val);
        }

        void SetValue(ARM_VALUE value)
        {
            val = value;
        }

        ARM_VALUE_TYPE GetValueType(void)
        {
            return(valType);
        }

        void SetIntVal(int v)
        {
            valType = ARM_INT;

            val.intVal = v;
        }

        void SetDoubleVal(double v)
        {
            valType = ARM_DOUBLE;

            val.doubleVal = v;
        }

        void SetStringVal(char* str)
        {
            valType = ARM_STRING;

            strcpy(val.stringVal, str);
        }

		void View(char* id = NULL, FILE* ficOut = NULL);
};









#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
