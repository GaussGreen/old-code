/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : listeval.h                                                   */
/*                                                                            */
/* DESCRIPTION : Object dealing with a list of dates and value                */
/*                                                                            */
/* DATE        :     Sep 30 1996                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/
 
 
#ifndef _LISTEVAL_H
#define _LISTEVAL_H


#include "armglob.h"
#include "linalg.h"




class ARM_Listeval : public ARM_Object
{
    private:

		int        itsSize ;
        ARM_Vector* itsBeginPeriod ;      // Borne comprise
        ARM_Vector* itsEndPeriod ;        // Borne comprise
        ARM_Vector* itsValue ;            // Borne comprise

    public:

        ARM_Listeval (void);

        ARM_Listeval (int size, ARM_Vector* jDateDeb, ARM_Vector* jDateFin);
        // Other Formats should be provided in further versions

        ARM_Listeval (int size, ARM_Vector* jDateDeb, 
                      ARM_Vector* jDateFin, ARM_Vector* value);

        ARM_Listeval (double value);

        void Set(int size, ARM_Vector* jdeb, ARM_Vector* jfin, ARM_Vector* value);

        ~ARM_Listeval(void) 
        {
            if (itsBeginPeriod) 
               delete itsBeginPeriod;

            if (itsEndPeriod) 
               delete itsEndPeriod;

            if (itsValue) 
               delete itsValue;
        }

        int GetSize (void) { return itsSize;}

		ARM_Listeval(const ARM_Listeval &);

		ARM_Listeval& operator = (const ARM_Listeval &);

        int GetFlagBelong(double jDatePar);

		double GetValue(double jDatePar);

		void SetValue(double aValue);
};

#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
