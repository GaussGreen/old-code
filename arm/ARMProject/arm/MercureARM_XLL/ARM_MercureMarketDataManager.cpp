#include "ARM_MercureMarketDataManager.h"
#include "MarketDataDictionary.h"
#include "expt.h"
#include <string>


void ARM_MercureMarketDataManager::View(char* id, FILE* ficOut)
{	
    FILE*	vOutputFile;
    char	vOutputFileName[40];

    if ( ficOut == NULL )
    {
        ARM_GetViewFile(ARM_VIEW_FILE, id, vOutputFileName);        
        vOutputFile = fopen(vOutputFileName, "w");
    }
    else
    {
        vOutputFile = ficOut;
    }

    try
    {
		map<string, MarketDataDictionary*>&	vDictionaries = GetMarketDataManager()->GetMktDictionaries();
		map<string, MarketDataDictionary*>::iterator	vDicoIter = vDictionaries.find(itsAsOfDate);

		if( vDicoIter != vDictionaries.end() )
		{
			MarketDataDictionary*	vMktDataDico = vDicoIter->second;
			map<string, MetaMarketData*>::iterator	vDataIter = vMktDataDico->GetDictionary().begin();
			map<string, MetaMarketData*>::iterator	vDataEnd = vMktDataDico->GetDictionary().end();

			for(; vDataIter != vDataEnd; vDataIter++)
			{
				string	vKey = vDataIter->first;
				string	vStars("*");
				int	vKeyLenght = vKey.length();
				for(int i=1; i<vKeyLenght; i++)
					vStars += "*";

				fprintf(vOutputFile, "\n\n%s", vStars.c_str());
				fprintf(vOutputFile, "\n%s\n", vKey.c_str());
				fprintf(vOutputFile, "%s\n", vStars.c_str());

				MetaMarketData*	vMktData = vDataIter->second;
				vMktData->GetMarketData().View(id, vOutputFile);
			}
		}
    }

	catch(Exception& )
    {
        throw Exception(__LINE__, __FILE__, ERR_CALC_DET_PB, "Error reading result string");
    }

    if ( ficOut == NULL )
    {
       fclose(vOutputFile);
    }
}
