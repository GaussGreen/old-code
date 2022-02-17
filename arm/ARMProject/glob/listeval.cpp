/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : listeval.cpp                                                 */
/*                                                                            */
/* DESCRIPTION : Methods of listeval Object                                   */
/*                                                                            */
/* DATE        :     Sep 30 1996                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/


#include "listeval.h"


/*----------------------------------------------------------------------------*/
/*   Constructors and destructor                                              */
/*----------------------------------------------------------------------------*/


ARM_Listeval::ARM_Listeval(void)
{

	SetName(ARM_LISTEVAL);
	itsSize=0;
	itsBeginPeriod = NULL;
	itsEndPeriod = NULL;
	itsValue = NULL;
}



ARM_Listeval::ARM_Listeval(int size, ARM_Vector* jdatdeb, ARM_Vector* jdatfin)
{

    SetName(ARM_LISTEVAL);

	itsSize = size;
	itsBeginPeriod = jdatdeb;
	itsEndPeriod = jdatfin;
	itsValue = NULL;

}



ARM_Listeval::ARM_Listeval(int size, ARM_Vector* jdatdeb, 
                           ARM_Vector* jdatfin, ARM_Vector* value)
{

    SetName(ARM_LISTEVAL);

	itsSize = size;
	itsBeginPeriod = jdatdeb;
	itsEndPeriod =  jdatfin;
	itsValue =  value;

}



ARM_Listeval::ARM_Listeval(double value)
{

	SetName(ARM_LISTEVAL);
 
    itsSize = 1;
	itsValue = new ARM_Vector(1, value);
    itsBeginPeriod = NULL;
    itsEndPeriod = NULL;


}



void ARM_Listeval::Set(int size, ARM_Vector* jdatdeb, ARM_Vector* jdatfin, 
						ARM_Vector* value)
{

	itsSize = size;
	itsBeginPeriod = jdatdeb;
	itsEndPeriod = jdatfin;
	itsValue = value;

}



ARM_Listeval::ARM_Listeval(const ARM_Listeval& l)
{

	SetName(ARM_LISTEVAL);

	itsSize = l.itsSize;
	itsBeginPeriod = l.itsBeginPeriod;
	itsEndPeriod =  l.itsEndPeriod;
	itsValue =  l.itsValue;

}


ARM_Listeval& ARM_Listeval::operator = (const ARM_Listeval& l)
{
	(*this).ARM_Object::operator =(l);
    itsSize = l.itsSize;
    itsBeginPeriod = l.itsBeginPeriod;
    itsEndPeriod =  l.itsEndPeriod;
    itsValue =  l.itsValue;

	return(*this);

}




/*----------------------------------------------------------------------------*/
/* Methods                                                                    */
/*----------------------------------------------------------------------------*/

int ARM_Listeval::GetFlagBelong(double jdatePar)
{
	int i, flag = 0;

	for (i=0; i< itsSize; i++)
	{
    	if (jdatePar >= (*itsBeginPeriod)[i]
     	    && jdatePar <= (*itsEndPeriod)[i])
        {
		   flag = 1;
        }
	}

	return(flag);
}



double ARM_Listeval::GetValue(double jdatePar)
{
	int i;
	double res=0.0;



	if (itsSize = 1)
    {
       res = (*itsValue)[0];
       return(res);
    }

	for (i=0; i< itsSize; i++)
	{
       if (jdatePar >= (*itsBeginPeriod)[i] 
                     && jdatePar <= (*itsEndPeriod)[i])
       {
	      res = (*itsValue)[i];
       }
	}

	return(res);
}



void ARM_Listeval::SetValue(double aValue)
{
	int i;

	for (i=0; i< itsSize; i++)
		(*itsValue)[i] = aValue;
}



/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
