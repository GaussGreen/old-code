#ifndef _ARM_ABSTRACTMARKETCLASS_H
#define _ARM_ABSTRACTMARKETCLASS_H


                      /*******************************************/
                      /*                                         */
                      /*       A Market data abstract class:     */
                      /*   Useful for market data inheritance    */
                      /*                                         */
                      /*******************************************/


#include "armglob.h"
#include <string>
#include "dates.h"


class ARM_AbstractMarketClass : public ARM_Object
{
    private:
	
	    string itsStrCurrency;
	
        string itsStrType;			// "ZERO", "BASIS USD", "VOL SWOPT HIST", ...

        string itsExternalIndex;	// EURIB, LIBOR, AUBB ...

        string itsExternalCrvId;	// MO, ...: The real curve Id!!!


	protected :

		ARM_Date itsAsOfDate;


        void Init(void)
        {
           string defCcy(ARM_DEFAULT_COUNTRY);


           // Initialize

           SetName(ARM_ABSTRACTMARKETCLASS);


           itsStrCurrency   = defCcy;

           itsStrType       = string("");

           itsExternalIndex = string("");

           itsExternalCrvId = string("");
        }

    public:

   	    ARM_AbstractMarketClass(void)
	    {
            Init();
	    }
	    
        ARM_AbstractMarketClass(const ARM_AbstractMarketClass& mktObj) : ARM_Object(mktObj)
	    {
            Init();

            BitwiseCopy(&mktObj);
	    }

		virtual ~ARM_AbstractMarketClass(void)
        {
	    }

		virtual int IsMarketData(void)
        {
            return(1);
        }   

        void BitwiseCopy(const ARM_Object* srcObj)
        {
            ARM_AbstractMarketClass* src = (ARM_AbstractMarketClass *) srcObj;


            itsStrCurrency   = src->itsStrCurrency;

            itsStrType       = src->itsStrType;

            itsExternalIndex = src->itsExternalIndex;

            itsExternalCrvId = src->itsExternalCrvId;

			itsAsOfDate = src->itsAsOfDate;
        }
	    
		virtual void Copy(const ARM_Object* srcMkt)
        {
		    ARM_Object::Copy(srcMkt);
		    
		    BitwiseCopy(srcMkt);
	    }
	    
		virtual ARM_Object* Clone(void)
	    {
		    ARM_AbstractMarketClass* theClone = new ARM_AbstractMarketClass();
		    
		    theClone->Copy(this);
		    
		    return(theClone);
	    }

		virtual ARM_AbstractMarketClass& operator = (const ARM_AbstractMarketClass& srcMkt)
	    {
		    (*this).ARM_Object::operator = (srcMkt);

            BitwiseCopy(&srcMkt);

		    return(*this);
	    }

        string GetStrCurrency(void) const
        {
             return(itsStrCurrency);
        }

        void SetStrCurrency(string& strCcy)
        {
             itsStrCurrency = strCcy;
        }

        void SetStrCurrency(char* aCcy)
        {
             string ccy(aCcy);

             itsStrCurrency = ccy;
        }

        string GetStrType(void) const
        {
             return(itsStrType);
        }

        void SetStrType(string& aType)
        {
            itsStrType = aType;
        }

        string GetExternalIndex(void) const
        {
            return(itsExternalIndex);
        }

        void SetExternalIndex(string& index)
        {
             itsExternalIndex = index;   
        }

        string GetExternalCrvId(void) const
        {
            return(itsExternalCrvId);
        }

        void SetExternalCrvId(string& cvId)
        {
             itsExternalCrvId = cvId;   
        }

		ARM_Date GetAsOfDate(void) const
		{
			return	itsAsOfDate;
		}

		void SetMktExternalCharacteristics(string& aIndex,
										   string& aCurrency,
										   string& aCvName,
										   string& aType)
		{  
			SetExternalIndex(aIndex);
			SetStrCurrency(aCurrency);
			SetExternalCrvId(aCvName);
			SetStrType(aType);
		}

		//virtuelle car par exemple pour la volCurve, on set aussi itsLastKnownDate
		virtual void SetAsOfDate(ARM_Date& aDate)
		{
			itsAsOfDate = aDate;
		}

        void View(char* id = NULL, FILE* ficOut = NULL)
        {
		    FILE* fOut;
		    char fOutName[200];
		    
		    
		    if ( ficOut == NULL )
		    {
			    ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
			    
			    fOut = fopen(fOutName, "w");
		    }
		    else
		    {
			    fOut = ficOut;
		    }
		    
		    fprintf(fOut, "\n\n >>>>>>> MARKET DATA characteristics : \n");
		    
            fprintf(fOut, "\n\n CURRENCY : %s \n", (const char *) itsStrCurrency.c_str());
            fprintf(fOut, "\n\n     TYPE : %s \n", (const char *) itsStrType.c_str());
            fprintf(fOut, "\n\n    INDEX : %s \n", (const char *) itsExternalIndex.c_str());
            fprintf(fOut, "\n\n     CVID : %s \n", (const char *) itsExternalCrvId.c_str());
			fprintf(fOut, "\n\n AsOfDate : %s \n", (const char *) itsAsOfDate.toString('.').c_str());

            fprintf(fOut, "\n\n <<<<<<< End Of MARKET DATA characteristics : \n");

            if ( ficOut == NULL )
            {
		        fclose(fOut);
            }
        }	
};


#endif
