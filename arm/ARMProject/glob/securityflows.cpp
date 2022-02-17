#include "firsttoinc.h"
#include "securityflows.h"
#include "linalg.h"


ARM_SecurityFlows::ARM_SecurityFlows()
{
	SetName(ARM_SECURITY_FLOWS);
}


ARM_SecurityFlows::ARM_SecurityFlows(vector<string>& aLabels, vector<ARM_Vector*>& aValues)
{
	SetName(ARM_SECURITY_FLOWS);

	int	vSize = aLabels.size();
	if(vSize != aValues.size())
	{
		throw "Vector of labels and vector of values should have the same size";
	}

	itsLabels = aLabels;
	itsValues.resize(vSize);
	int	i=0;
	for(vector<ARM_Vector*>::iterator vIter = aValues.begin(); vIter != aValues.end(); vIter++)
	{
		ARM_Vector*	vValue = (ARM_Vector*)(*vIter)->Clone();
		itsValues.at(i) = vValue;
		itsData[itsLabels.at(i)] = vValue;
		i++;
	}
}


ARM_SecurityFlows::~ARM_SecurityFlows()
{
	for(vector<ARM_Vector*>::iterator vIter = itsValues.begin(); vIter != itsValues.end(); vIter++)
		delete	*vIter;
}


ARM_Vector*	ARM_SecurityFlows::GetValues(string& aLabel)
{
	map<string, ARM_Vector*>::iterator	vIter = itsData.find(aLabel);
	
	if(vIter == itsData.end())
		throw	"Check key label";

	return	vIter->second;
}


void ARM_SecurityFlows::SetValues(string& aLabel, ARM_Vector* aValues)
{
	map<string, ARM_Vector*>::iterator	vIter = itsData.find(aLabel);
	
	if(vIter != itsData.end())
	{
		ARM_Vector*	vCurrentValue = vIter->second;

		if(vCurrentValue != aValues)
			delete	vCurrentValue;
	}

	itsData[aLabel] = aValues;
}


void ARM_SecurityFlows::View(char* id, FILE* ficOut)
{
    FILE*	fOut;
    char	fOutName[200];


    if( ficOut == NULL )
    {
       ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);

       (void) unlink(fOutName);

       fOut = fopen(fOutName, "w");
    }
    else
    {
       fOut = ficOut;
    }

    fprintf(fOut, "\n\n >>>>>>>>>>>>>>>>>>>>>> Security Flows <<<<<<<<<<<<<<<<<<<<<<<\n\n");

	for(vector<string>::iterator vIter = itsLabels.begin(); vIter != itsLabels.end(); vIter++)
	{
		fprintf(fOut, "%s \t", (*vIter).c_str());
	}

	vector<ARM_Vector*>::iterator vIter2 = itsValues.begin();
	int	vSize = (*vIter2)->size();
	for(int i=0; i<vSize; i++)
	{
		fprintf(fOut, "\n");

		for(vIter2 = itsValues.begin(); vIter2 != itsValues.end(); vIter2++)
		{
			fprintf(fOut, "%lf \t", (*vIter2)->Elt(i));
		}
	}
    
	fprintf(fOut, "\n\n <<<<<<<<<<<<<<<<<<<<<< Security Flows >>>>>>>>>>>>>>>>>>>>>>>\n\n");

    if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}