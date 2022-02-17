#ifndef ARM_MERCURE_RESULT_H
#define ARM_MERCURE_RESULT_H

#include "armglob.h"
#include <string>
#include <vector>

using namespace std;


class ARM_MercureResult : public ARM_Object
{
	private:
		vector<string>*	itsXMLOutput;
		vector<string>*	itsTXTOutput;

	public:
		ARM_MercureResult(void)
		{
			itsXMLOutput = NULL;
			itsTXTOutput = NULL;
			SetName(ARM_MERCURE_RESULT);
		}

		ARM_MercureResult(vector<string>* aXMLOutput, vector<string>* aTXTOutput)
		{
			itsXMLOutput = aXMLOutput;
			itsTXTOutput = aTXTOutput;
			SetName(ARM_MERCURE_RESULT);
		}

		~ARM_MercureResult(void)
		{
			delete	itsXMLOutput;
			delete	itsTXTOutput;
		}

		virtual void View(char* id = NULL, FILE* ficOut = NULL);

		void ViewSensi(string& aSensiType);

		double GetPostProcessedData(string& aSensiType, string& aPlot);
};


class ARM_MercureHelp : public ARM_Object
{
	private:
		string itsMercureObjectName;

	public:
		ARM_MercureHelp(void)
		{
			SetName(ARM_MERCURE_HELP);
		}
		
		ARM_MercureHelp(string objectName)
		{
			itsMercureObjectName = objectName;
			SetName(ARM_MERCURE_HELP);
		}

		virtual void View(char* id = NULL, FILE* ficOut = NULL);
};

#endif