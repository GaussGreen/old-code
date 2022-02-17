/*===========================================================
  Name    : dk_wrapper_dkf_tree.cpp
  Owner   : DK
  Created : DK  
  Comment : 3F tree using the dk model.
  Time    : 16:38 12/09/2002
=============================================================*/
#include <iostream>
#include <strstream>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifndef WIN32
#include <sys/param.h>
#endif 
#define ACC 1.0e-04
#define JMAX 5000

#include "dk_utils.h"

#include "DKMaille.h"
#include "DKMaille2D.h"

#include "dk_wrapper_dk3f_tree.h"
#include "dk_multifactor_skewtree.h"

// creates an array of times, timeArray[], passing through each of the inputDates 
// the array starts at inputDates[1] and ends at inputDates[howManyDates] */
// howManySteps is updated to return the actual number of steps in the tree */
void datesForSteps(DKMaille<double> inputDates, 
                   int howManyDates, 
                   int *howManySteps, 
                   DKMaille<double> &timeArray)
{
	double startdate;
	double stepsize;
	double newstepsize;
	int stepgap;
	int counter;
  register int j;
  register int i;
  timeArray[0]=0.;
	counter=0;
 
	startdate=inputDates.at(1);
	stepsize=(inputDates.at(howManyDates)-startdate)/(*howManySteps);
 
	for(i=2;i<=howManyDates;i++)
	{
		// Check dates are not equal
		if(inputDates.at(i)!=inputDates.at(i-1))
		{
      if((inputDates.at(i)-inputDates.at(i-1))/stepsize<1.) {
        newstepsize=(inputDates.at(i)-inputDates.at(i-1));
        stepgap=1;
      }
      else 
      {
      stepgap=(int)(1.*(inputDates.at(i)-inputDates.at(i-1))/stepsize); 
      if((double)stepgap<(inputDates.at(i)-inputDates.at(i-1))/stepsize-0.5) stepgap+=1;
      newstepsize=(inputDates.at(i)-inputDates.at(i-1))/(double)stepgap;		
      }
      for (j=1;j<=stepgap;j++)
			{
				counter+=1;
				timeArray.at(counter)=timeArray.at(counter-1)+newstepsize/365.;
			}
		}
	}
	(*howManySteps)=counter; 
} // datesForSteps(...)

