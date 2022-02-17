#ifndef _SECURITYFLOWS_H
#define _SECURITYFLOWS_H

#include <string>
#include <vector>
#include <map>

class SecurityFlows
{
	private :
		vector<string>				itsLabels;
		vector<ARM_Vector*>			itsValues;
		map<string, ARM_Vector*>	itsData;

	public :
		SecurityFlows();
		SecurityFlows(vector<string>& aLabels, vector<ARM_Vector*>& aValues);
		~SecurityFlows();

		ARM_Vector*	GetValues(string& aLabel);
		void SetValues(string& aLabel, ARM_Vector* aValues);
}


#endif