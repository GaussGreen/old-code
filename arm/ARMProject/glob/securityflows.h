#ifndef _SECURITYFLOWS_H
#define _SECURITYFLOWS_H

#include "armglob.h"
#include <string>
#include <vector>
#include <map>

class ARM_Vector;
using namespace std;


class ARM_SecurityFlows : public ARM_Object
{
	private :
		vector<string>				itsLabels;
		vector<ARM_Vector*>			itsValues;
		map<string, ARM_Vector*>	itsData;

	public :
		ARM_SecurityFlows();
		ARM_SecurityFlows(vector<string>& aLabels, vector<ARM_Vector*>& aValues);
		~ARM_SecurityFlows();

		ARM_Vector*	GetValues(string& aLabel);
		void SetValues(string& aLabel, ARM_Vector* aValues);
		void View(char* id, FILE* ficOut);
};


#endif