#ifndef _INDEXINDEXCORRELCUBE_H
#define _INDEXINDEXCORRELCUBE_H


#include "firsttoinc.h"

#include "hypercube.h"
#include <string>
#include <map>



class ARM_IndexIndexCorrelCube : public ARM_HyperCube
{
	public :

		ARM_IndexIndexCorrelCube(void);

		ARM_IndexIndexCorrelCube(vector<ARM_VolCurve*>& aCorrelList, 
								 vector<string>& aTenor1List, vector<string>& aTenor2List);

		virtual ~ARM_IndexIndexCorrelCube(void);


		void Init(void);

		virtual ARM_Object* Clone(void);

		void Copy(const ARM_Object* aCurveIn);

		void BitwiseCopy(const ARM_Object* srcObject);

		bool isCorrelCurve(string& aCorrelStr) const;


		void View(char* id = NULL, FILE* ficOut = NULL);

		virtual double ComputeIndexIndexCorrel(	string aTenor1, double aExpiry1,
												string aTenor2, double aExpiry2) const;

		virtual double ComputeIndexIndexCorrel(	double aTenor1, double aExpiry1,
												double aTenor2, double aExpiry2) const;

		void OrderTenorsAndExpiries(int& aTenor1, int& aTenor2, 
									string& aTenor1Str, string& aTenor2Str,
									double& aExpiry1, double& aExpiry2) const;

		double ComputeCorrel(int aTenor1, int aTenor2, 
							 string aTenor1Str, string aTenor2Str,
							 double aExpiry1, double aExpiry2) const;

		ARM_VolCurve* GetCorrelDiag(string& aTenor1, vector<string>& aTenor2) const;

		ARM_VolCurve* GetCorrelCurve(string& aCorrelStr) const;

		void BumpVolatility(double value, int nthLine = 0, int nthCol = 0,
                            int isCumul = K_NO, int isAbsolute = K_YES);
};


#endif
