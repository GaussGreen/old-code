/*----------------------------------------------------------------------*


     File: fixingsched.cpp

     Class: ARM_FixingSched, a class to manage fixings

 *----------------------------------------------------------------------*/






/*----------------------------------------------------------------------*/

#include "fixingsched.h"

#include "gpbase/gpvector.h"


using namespace std;



CC_USING_NS(std,pair)


ARM_FixingSched::ARM_FixingSched(void)
{
	SetName(ARM_FIXING_SCHED);
}



ARM_FixingSched::ARM_FixingSched(const ARM_Date& AsOf)
{
    itsAsOfDate = AsOf;
	SetName(ARM_FIXING_SCHED);
}


void ARM_FixingSched::CleanUp(void)
{
    if (!itsRatesFixings.empty())
	{
	   map< string , ARM_Curve* >::iterator iter = itsRatesFixings.begin();

	   for ( ; iter != itsRatesFixings.end() ; iter++)
       {	
		   if (iter->second)
           {
			  delete iter->second;
				
              iter->second = NULL;
           }
		}
	}

    if (!itsFxFixings.empty() )
	{
		map< string , ARM_Curve* >::iterator iter = itsFxFixings.begin();

		for ( ; iter != itsFxFixings.end() ; iter++ )
		{	
			if (iter->second)
			{
			   delete iter->second;
				
               iter->second = NULL;
			}
		}
	}  
}



void ARM_FixingSched::BitwiseCopy(const ARM_Object* src)
{

	ARM_FixingSched* fixings = dynamic_cast<ARM_FixingSched*>(const_cast<ARM_Object*>(src));


    CleanUp();

    itsAsOfDate = fixings->itsAsOfDate;

	if (fixings)
	{
	   map< string , ARM_Curve* >::iterator iter = fixings->itsRatesFixings.begin();
		
       for ( ; iter != fixings->itsRatesFixings.end() ; iter++ )
		   itsRatesFixings[iter->first] = (ARM_Curve *)(iter->second)->Clone();


       // FX fixings

	   map<string , ARM_Curve*>::iterator fxIter = fixings->itsFxFixings.begin();

	   for ( ; fxIter != fixings->itsFxFixings.end() ; fxIter++ )
		   itsFxFixings[fxIter->first] = (ARM_Curve *)(fxIter->second)->Clone();
	}
}



void ARM_FixingSched::Copy(const ARM_Object* src)
{
	ARM_Object::Copy(src);

	BitwiseCopy(src);
}



ARM_Object* ARM_FixingSched::Clone()
{
	ARM_FixingSched* theClone = new ARM_FixingSched();

	theClone->Copy(this);

	return(theClone);
}



void ARM_FixingSched::View(char* id, FILE* ficOut)
{
    FILE* fOut;
    char fOutName[200];

    if ( ficOut == NULL )
    {
       ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
       (void) unlink(fOutName);
       fOut = fopen(fOutName, "w");
    }
    else
       fOut = ficOut;

    fprintf(fOut, "\n       >>>>>>>>>>> Libor and FX fixings <<<<<<<<<<<\n");

    fprintf(fOut, "\n ================> Libor fixings:\n");

    map< string , ARM_Curve* >::const_iterator liborIter    = itsRatesFixings.begin();

	map< string , ARM_Curve* >::const_iterator liborIterEnd	= itsRatesFixings.end();

	
	for ( ; liborIter != liborIterEnd; ++liborIter )
	{
		fprintf(fOut, "\n\nIndex:\t%s\n", ((*liborIter).first).c_str());

		ARM_Curve* liborFixings= ((*liborIter).second);

	    liborFixings->View(id, fOut);
	}

    fprintf(fOut, "\n ================> FX fixings:\n");

    map< string , ARM_Curve* >::const_iterator fxIter    = itsFxFixings.begin();

	map< string , ARM_Curve* >::const_iterator fxIterEnd = itsFxFixings.end();

	
	for ( ; fxIter != fxIterEnd; ++fxIter)
	{
        fprintf(fOut, "\n\n FX:\t%s\n", ((*fxIter).first).c_str());

		ARM_Curve* fxFixings = ((*fxIter).second);

	    fxFixings->View(id, fOut);
	}

    if ( ficOut == NULL )
       fclose(fOut);
}



ARM_FixingSched::~ARM_FixingSched(void)
{
    CleanUp();
}



int ARM_FixingSched::GetNbFixings(void)
{
	return itsRatesFixings.size();
}


ARM_Curve* ARM_FixingSched::GetLiborFixing(const string& liborKey) // ex: USD_LIBOR_6M
{
    ARM_Curve* liborFixing = NULL;

  
    map< string , ARM_Curve* >::const_iterator iter = itsRatesFixings.find(liborKey);
	
    if ( iter == itsRatesFixings.end() )
	{
	   return(NULL);   
	}

	liborFixing = (ARM_Curve *) iter->second;

    return(liborFixing);
}



void ARM_FixingSched::AddLiborFixing(const string& liborKey, ARM_Curve* fixings)
{
     map< string , ARM_Curve* >::iterator iter = itsRatesFixings.find(liborKey);

     if ( iter != itsRatesFixings.end() )
		delete iter->second;

     itsRatesFixings[liborKey] = (ARM_Curve *) fixings->Clone(); 
}



double ARM_FixingSched::GetLiborFixing(const ARM_Date& resetDate,
                                       const string& ccy,
                                       const string& index,
                                       const string& matu)
{
    double liborFixing = -1.0;

    string liborKey = ccy+string("_")+index+string("_")+matu;

    ARM_Curve* fixings = GetLiborFixing(liborKey);

    if ( fixings == NULL )
    {
       return(liborFixing);
    }

    double resetLag = (resetDate-itsAsOfDate);

    liborFixing = fixings->ExactSearch(resetLag);

    return(liborFixing);
}

ARM_Curve* ARM_FixingSched::GetFxFixing(const string& fxKey) // ex: USD_JPY 
{
	string fxName(fxKey);
	if(fxName.size() == 6) //in case EURJPY
		fxName = fxKey.substr(0,3)+"_"+fxKey.substr(3,3);
  
    map< string , ARM_Curve* >::const_iterator iter = itsFxFixings.find(fxName);
	
	ARM_Curve* fxFixing = NULL;
    if ( iter != itsFxFixings.end() )
	   	fxFixing = (ARM_Curve *) iter->second;

    return(fxFixing);
}

double ARM_FixingSched::GetFxFixing(const ARM_Date& resetDate, 
                                    const string& ccy1,
                                    const string& ccy2)
{
    double fxFixing = -1.0;

    string fxKey    = ccy1+"_"+ccy2;

    ARM_Curve* fixings = GetFxFixing(fxKey);

    if ( fixings == NULL )
    {
       return(fxFixing);
    }

    double resetLag = (resetDate-itsAsOfDate);

    fxFixing = fixings->ExactSearch(resetLag);

    return(fxFixing);
}

void ARM_FixingSched::AddFxFixing(const string& fxKey, ARM_Curve* fixings)
{
     map< string , ARM_Curve* >::iterator iter = itsFxFixings.find(fxKey);

     if ( iter != itsFxFixings.end() )
		delete iter->second;

     itsFxFixings[fxKey] = (ARM_Curve *) fixings->Clone(); 
}














/*---------------------------------------------------------------------------*/
/*---- End Of File */