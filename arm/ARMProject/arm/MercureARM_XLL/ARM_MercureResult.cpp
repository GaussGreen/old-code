#include "ARM_MercureResult.h"
#include "expt.h"
#include "CCString.h"
#include "gpbase\stringmanip.h"
#include "Global.h"
#include "Meta2IRFXModel.h"
#include "MetaBSGenModel.h"
#include "MetaBSModel.h"
#include "MetaBSSmiled.h"
#include "MetaMarkovTree.h"
#include "MetaMixtureModel.h"
#include "MetaModel3F.h"
#include "MetaModelFXOption3F.h"
#include "MetaQModel.h"
#include "MetaCRAModel.h"

using ARM::stringToUpper;
using ARM::stringGetUpper;
using namespace mercure;

void ARM_MercureResult::View(char* id, FILE* ficOut)
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
		int	vLen = strlen(vOutputFileName);
		if(vLen > 3)
		{
			vector<string>*	vOutput;
			char*	vExtension = vOutputFileName + vLen - 3;
			CCString	cExtension(vExtension);

			if(cExtension == "xml")
				vOutput = itsXMLOutput;
			else
				vOutput = itsTXTOutput;

			vector<string>::iterator	it = vOutput->begin();
			for (it; it != vOutput->end(); it++) 
			{
				fprintf(vOutputFile, "%s\n\n", it->c_str());
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

void ARM_MercureResult::ViewSensi(string& aSensiType)
{	
	string*	vOutputString = itsXMLOutput->begin();

	const string	vStringBegin( string("<Market_Data_Dim1><Class>") + aSensiType );
	const string	vStringEnd("</Hedge_Ratio_Class>");

	int	vIndexBegin = vOutputString->find(vStringBegin);
	int	vIndexEnd = vOutputString->find(vStringEnd, vIndexBegin);

	int		vSensiStringLenght = vIndexEnd - vIndexBegin;
	string	vSensiString = vOutputString->substr( vIndexBegin, vSensiStringLenght);


    FILE*	vOutputFile;
    char	vOutputFileName[40];

    CCString	id( CCString("123") + CCString("alekrafi") );
	ARM_GetViewFile(ARM_VIEW_FILE, id, vOutputFileName);        
	strcat(vOutputFileName, ".xml");
    vOutputFile = fopen(vOutputFileName, "w");

    try
    {
		fprintf(vOutputFile, "<Vega>\n");
		fprintf(vOutputFile, "%s\n", vSensiString.c_str());
		fprintf(vOutputFile, "</Vega>\n");
    }

	catch(Exception& )
    {
        throw Exception(__LINE__, __FILE__, ERR_CALC_DET_PB, "Error reading XML result string");
    }

	fclose(vOutputFile);
}

/// Get the string located at "tagName" in an XML source
/// Research starts at "startPoint"
string GetValue(char* tagName, string& source, int startPoint)
{
	char tag[50];
	sprintf(tag, "<%s>", tagName);
	string tagBegin = string(tag);

	sprintf(tag, "</%s>", tagName);
	string tagEnd = string(tag);
	
	int valuePos1 = source.find(tagBegin, startPoint);
	if (valuePos1 == string::npos) // not found !
		return string();

	int valuePos2 = source.find(tagEnd, startPoint);

	string value = source.substr(valuePos1 + tagBegin.size(), valuePos2 - valuePos1 - tagBegin.size());

	return value;
}

/// Display value computed for "aSensiName" and plot "aPlot"
/// this value must have been computed in Mercure postprocessing
double ARM_MercureResult::GetPostProcessedData(string& aSensiName, string& aPlot)
{
	double sensiValue = 0.0;
	bool sensiFound = false;
	vector<string>*	vOutput = itsXMLOutput;

	vector<string>::iterator it = vOutput->begin();
	for (it; it != vOutput->end(); it++) 
	{
		int position = 0;
		stringToUpper(it[0]);
		/// PV case :
		/// 1st arg (sensiName) must be "PV"
		/// 2nd arg (plot) is the postprocessed data (if 2nd arg is empty, get global PV)
		if (aSensiName == "PV")
		{
			if (aPlot != "")
			{
				while (it[0].find("<HEDGE_RATIO>", position) != string::npos)
				{
					int ppPos = it[0].find("<POSTPROCESSED>", position);
					int ppNamePos = it[0].find(aPlot, position);
					if ( (ppPos != string::npos) && ( ppNamePos != string::npos))
					{
						string value = GetValue("VALUE", it[0], ppNamePos);
						sensiValue = atof(value.c_str());
						sensiFound = true;
					}
					position = it[0].find("</HEDGE_RATIO>", position);
				}
			}
			else
			{
				// extract global post processed value
				int gppPos = it[0].find("<GLOBALPOSTPROCESSED>", position);
				if (gppPos != string::npos)
				{
					string value = GetValue("VALUE", it[0], gppPos);
					sensiValue = atof(value.c_str());
					sensiFound = true;
				}
				position = it[0].find("</GLOBALPOSTPROCESSED>", position);
			}
		}
		/// 1st order sensi case :
		/// 1st arg is the market data
		/// 2nd arg the plot
		else
		{
			string dataClass, dataCcy, dataIndex;
			int slashPos1 = aSensiName.find("/");
			dataClass = aSensiName.substr(0, slashPos1);
			int slashPos2 = aSensiName.find("/", slashPos1+1);
			dataCcy = aSensiName.substr(slashPos1+1, slashPos2-slashPos1-1);
			dataIndex = aSensiName.substr(slashPos2+1, aSensiName.size()-slashPos2-1);

			while (it[0].find("<HEDGE_RATIO_CLASS>", position) != string::npos)
			{
				int mdPos = it[0].find("<MARKET_DATA_DIM1>", position);
				string scenarioClass = GetValue("CLASS", it[0], mdPos);
				string scenarioCcy = GetValue("CCY", it[0], mdPos);
				string scenarioIndex = GetValue("INDEX", it[0], mdPos);
				if ( (scenarioClass == dataClass) && (scenarioCcy == dataCcy) && (scenarioIndex == dataIndex) )
				{
					// extract global post processed value
					int gppPos = it[0].find("<GLOBALPOSTPROCESSED>", position);
					int plotPos = it[0].find(aPlot, gppPos);
					if ( (gppPos != string::npos) && (plotPos != string::npos) )
					{
						string value = GetValue("VALUE", it[0], plotPos);
						sensiValue = atof(value.c_str());
						sensiFound = true;
					}
				}
				position = it[0].find("</HEDGE_RATIO_CLASS>", ++position);
			}
		}
	}

	if (!sensiFound)
	{
		char msg[300];
		if (it[0].find("ERROR"))
			sprintf(msg, "ARM_MercureResult::GetPostProcessedData : an error occured during hedge process\n");
		else
			sprintf(msg, "ARM_MercureResult::GetPostProcessedData : '%s' has not been specified in config file\n", aSensiName.c_str());
		
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg);						
	}

	return sensiValue;
}

void ARM_MercureHelp::View(char* id, FILE* ficOut)
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

	MetaModel* model = NULL;

	fprintf(fOut, "\n\nEXPECTED EXTRA MODEL PARAMETERS\n\n");
	fprintf(fOut, "NB 1 : Optional parameters have default values\n");
	fprintf(fOut, "If no default value is specified, the parameter must be in the HedgeRatio file\n\n");
	fprintf(fOut, "NB 2 : 'YES' and 'NO' values can be replaced by string 'Y', 'N', or integer 1, 0\n\n");

	if ( (stringGetUpper(itsMercureObjectName) == stringGetUpper(Model3FType)) || (itsMercureObjectName == "ALL") )
	{
		model = new MetaModel3F();
		model->View(id, fOut);
		delete model;
	}
	if ( (stringGetUpper(itsMercureObjectName) == stringGetUpper(QModelType)) || (itsMercureObjectName == "ALL") )
	{
		model = new MetaQModel();
		model->View(id, fOut);
		delete model;
	}
	if ( (stringGetUpper(itsMercureObjectName) == stringGetUpper(BSType)) || (itsMercureObjectName == "ALL") )
	{
		model = new MetaBSModel();
		model->View(id, fOut);
		delete model;
	}
	if ( (stringGetUpper(itsMercureObjectName) == stringGetUpper(FXOption3FType)) || (itsMercureObjectName == "ALL") )
	{
		model = new MetaModelFXOption3F();
		model->View(id, fOut);
		delete model;
	}
	if ( (stringGetUpper(itsMercureObjectName) == stringGetUpper(BSSmiledType)) || (itsMercureObjectName == "ALL") )
	{
		model = new MetaBSSmiled();
		model->View(id, fOut);
		delete model;
	}
	if ( (stringGetUpper(itsMercureObjectName) == stringGetUpper(BSGenType)) || (itsMercureObjectName == "ALL") )
	{
		model = new MetaBSGenModel();
		model->View(id, fOut);
		delete model;
	}
	if ( (stringGetUpper(itsMercureObjectName) == stringGetUpper(CRAModelType)) || (itsMercureObjectName == "ALL") )
	{
		model = new MetaCRAModel();
		model->View(id, fOut);
		delete model;
	}
	if ( (stringGetUpper(itsMercureObjectName) == stringGetUpper(MarkovTreeType)) || (itsMercureObjectName == "ALL") )
	{
		model = new MetaMarkovTree();
		model->View(id, fOut);
		delete model;
	}
	if ( (stringGetUpper(itsMercureObjectName) == stringGetUpper(TwoIRFXModelType)) || (itsMercureObjectName == "ALL") )
	{
		model = new Meta2IRFXModel();
		model->View(id, fOut);
		delete model;
	}
	if ( (stringGetUpper(itsMercureObjectName) == stringGetUpper(MixtureType)) || (itsMercureObjectName == "ALL") )
	{
		model = new MetaMixtureModel();
		model->View(id, fOut);
		delete model;
	}

    if ( ficOut == NULL )
        fclose(fOut);
}