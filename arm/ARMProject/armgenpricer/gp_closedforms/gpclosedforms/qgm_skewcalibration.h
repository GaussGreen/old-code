/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file sabr_calibration.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_QGM_SKEWCALIBRATION_H
#define _GP_CF_QGM_SKEWCALIBRATION_H

#include "gpbase/port.h"
#include "gpbase/typedef.h"
#include "gpbase/gpvector.h"


CC_BEGIN_NAMESPACE(ARM)
class Optimization_Result_Set;

class QGM_ParameterSet
{
private:
		std::vector<double>* _X;
		double _objective;

public:
		std::vector<double>* get_X() {return _X;}
		double get_objective() {return _objective;}

		QGM_ParameterSet()
		:_X(),_objective(0)
		{}

		QGM_ParameterSet(Optimization_Result_Set* r);
		~QGM_ParameterSet() {delete _X;}
};

QGM_ParameterSet* QGM_CalibrateSkew(std::vector<double>* tagetVect,
	  std::vector<double>* weights,
      std::vector<double>* precisions,
	  std::vector<double>* InitVector,
	  std::vector<double>* LBoundVector,
	  std::vector<double>* UBoundVector,
	  int algorithm
	  );

CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

