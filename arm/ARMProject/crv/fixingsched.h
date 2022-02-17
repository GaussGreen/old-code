/*----------------------------------------------------------------------*


     File: fixingsched.h

     Class: ARM_FixingSched, a class to manage fixings


 *----------------------------------------------------------------------*/
#ifndef _FIXINGSCHED_H
#define _FIXINGSCHED_H


#include <firsttoinc.h>

#include "dates.h"

#include "gpbase/curve.h"
#include "gpbase/curvetypedef.h"


#include <string>
#include <map>
#include <vector>

using std::map;
using std::string;
using ARM::ARM_Curve;



class ARM_FixingSched : public ARM_Object
{
       private:

              ARM_Date itsAsOfDate;

              map< string , ARM_Curve* > itsRatesFixings; // Libor fixings

              map< string , ARM_Curve* > itsFxFixings;    // Fx fixings

              void CleanUp(void);

       public:

              ARM_FixingSched(void);
			  ARM_FixingSched(const ARM_FixingSched& rhs); // to be implemented
			  ARM_FixingSched& operator=( const ARM_FixingSched& rhs ); //to be implemented

             ~ARM_FixingSched(void);
			 
			 ARM_FixingSched(const ARM_Date& asOfDate);


              void BitwiseCopy(const ARM_Object* src);
	          void Copy(const ARM_Object* src);
	          ARM_Object* Clone(void);
	          void View(char* id = NULL, FILE* ficOut = NULL);

              ARM_Date GetAsOfDate(void)
              {
                  return(itsAsOfDate);
              }

              void SetAsOfDate(const ARM_Date& asOfDate)
              {
                   itsAsOfDate = asOfDate;
              }

			  int GetNbFixings(void);
              ARM_Curve* GetLiborFixing(const string& liborKey); // ex: USD_LIBOR_6M                                                                                                                                  
                                                                 //     USD_LIBOR_10Y

              double GetLiborFixing(const ARM_Date& resetDate, 
                                    const string& ccy,
                                    const string& index,
                                    const string& matu);


              ARM_Curve* GetFxFixing(const string& fxKey);// ex: USD_JPY
                                                          //     EUR_JPY

              double GetFxFixing(const ARM_Date& resetDate, 
                                 const string& ccy1,
                                 const string& ccy2);

              // Append fixings

              void AddLiborFixing(const string& liborKey, ARM_Curve* fixings);

              void AddFxFixing(const string& fxKey, ARM_Curve* fixings);                                                               
};


















#endif
/*---------------------------------------------------------------------------*/
/*---- End Of File */
