/*---------------------------------------------------------------------------------------*/
/*                                                                                       */
/* ARM_HyperCube class : a cube of cubes useful for managing correlations                */
/*                                                                                       */
/*---------------------------------------------------------------------------------------*/
#ifndef _HYPERCUBE_H
#define _HYPERCUBE_H



#include "firsttoinc.h"

#include "volcube.h"
#include <string>
#include <map>
#include <vector>

using namespace std;

class ARM_HyperCube : public ARM_VolCube
{

	protected :
	
		map<string, ARM_VolCurve*>	itsVolCurves;

		bool	itsIntersurfaceInterpol;

		// Liste des tenors en nb de mois
		vector<int>	itsTenorList;

		// A mettre ailleurs en template
		void GetBounds(const vector<int>& aVector, int aTenor, int& aLowerBound, int& aUpperBound) const;


	public :

		ARM_HyperCube(void);

		ARM_HyperCube(vector<ARM_VolCurve*>& aVolCubeList, vector<string>& aTenorList);

		virtual ~ARM_HyperCube(void);

		void SetIntersurfaceInterpol(bool aValue)	{itsIntersurfaceInterpol = aValue;}


		void Init(void);

		virtual ARM_Object* Clone(void);

		void Copy(const ARM_Object* aCurveIn);

		void BitwiseCopy(const ARM_Object* srcObject);


		void View(char* id = NULL, FILE* ficOut = NULL);

		virtual double ComputeVolatility(double aOptMat, double aTenor1,
                                         double aTenor2, double aStrike = 0.0);

		virtual double ComputeCorrelByExpiry(double aExpiry, double aTenor1, double aTenor2);

		virtual double ComputeHyperCorrel(double aOptMat, double aTenor1,
										  double aTenor2, double aStrike = 0.0);

		virtual ARM_VolCurve*	GetVolCurve(string& aTenor);

		vector<int>	GetTenorList(void);

		ARM_VolCube*	CreateCorrelCubeByExpiry(vector<double>& aTenor1List, 
												 vector<double>& aTenor2List, 
												 vector<double>& aExpiryList);

		void BumpVolatility(double value, int nthLine = 0, int nthCol = 0,
                            int isCumul = K_NO, int isAbsolute = K_YES);

        virtual void BumpSmile(double value, double aCubeTenor = 0., double aSmileTenor = 0.,
							   int nthLine = 0, int nthCol = 0,
							   int isCumul = K_NO, int isAbsolute = K_YES);
};


#endif
/*---------------------------------------------------------------------------------------*/
/*---- End Of File ----*/