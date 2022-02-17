//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : LocalVolGrid.cpp
//
//   Description : local vol grid support interpolations
//
//   Author      : Ning Shen
//
//   Date        : 10 Jan 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/TRACE.hpp"
#include "edginc/Nrfns.hpp"
#include "edginc/Maths.hpp"
#include "edginc/LocalVolGrid.hpp"

DRLIB_BEGIN_NAMESPACE

#define TRACE_VIEW(X)                                   \
    do {                                                \
        if (TRACE_enabled) {                            \
            TRACE_BLOCK(#X);                            \
            (X).outputWrite("", "", TRACE_stream());    \
        }                                               \
    }                                                   \
    while (0)


LocalVolGrid::LocalVolGrid(const DoubleArray& fwds, 
                           CVolProcessedDVFConstSP vol, 
                           const DateTimeArray& futurePathDates,
                           const DateTimeArray& cumDates,
                           int numVolGrid, 
                           double stdevGridRange):
gridDates(futurePathDates), dateOffset(0) {

    static const string method = "LocalVolGrid::LocalVolGrid";

    try {
        TRACE_METHOD;

    // debug testing
    //#define PRINT_DEBUG_GRID
    #ifdef PRINT_DEBUG_GRID
        FILE * file = fopen("C:\\temp\\MCLV-grid.txt","w");
    #endif
    // debug testing

        const double MAX_DRIFT_RAND_NUMBER = 3.0; /* what bounds to use when
                                                    calculating max drift */
        // for spline
        const double y1 = 2e30; // default to natural spline
        const double yn = 2e30;

        static const string routine("LocalVolGrid::generateLVGrid");

        numVolGrid = (numVolGrid/2)*2 + 1; // make it an odd number to have a grid point at middle

        // lay out grids in stdev units
    //    const double lvStdevRange[21]  = {-6.0, -4.0, -3.0, -2.5, -2.0, -1.5, -1.0, -0.75, -0.5, -0.25,
    //                                                  0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0};
    /*    const double lvStdevRange[41] = {-6.0, -5.0, -4.0, -3.5, -3.0, -2.75, -2.5, -2.0, -1.75, -1.5, -1.25,
        -1.0, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,
        0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0,
        1.25, 1.5, 1.75, 2.0, 2.5, 2.75, 3.0, 3.5, 4.0, 5.0, 6.0};
    */
        vector<double> lvStdevRange(numVolGrid);
        for (int k=0; k<numVolGrid; k++)
        {
            lvStdevRange[k] = -stdevGridRange + 2.0*k*stdevGridRange/(numVolGrid-1);
        }

        int numPathSteps = futurePathDates.size();
        int iStep, j;

        double totalVar = 0.0; // used for max drift calc
        double sqrTotalVar = 0.0;

        TimeMetricConstSP metric = vol->GetTimeMetric();

        DateTimeArray t(2);
        t[1] = futurePathDates[0];

        // local vol*sqrt(dt) for steps=0 to T-1
        spotGrid.resize(numPathSteps); // will remove last one
        lvGrid.resize(numPathSteps);
        lvY2.resize(numPathSteps);

        double theDrift = 1.0;

        // compute up to T first for averaging
        for (iStep =0; iStep < numPathSteps; iStep++)
        {
            if (iStep == numPathSteps -1)
            {
                int days = t[1].daysDiff(t[0]);
                t[0] = t[1];
                t[1] = t[0].rollDate(days);
            }
            else
            {
                t[0] = t[1];
                t[1] = futurePathDates[iStep+1];

                CSliceDouble strike(const_cast<double*>(&fwds[iStep]), 1);
                double fwdVar = 0.0;
                CSliceDouble var(&fwdVar, 1);
                vol->CalcLocVar(&strike, t, &var, false); // no intradayinterp
                
                double drift = fwds[iStep+1] / fwds[iStep] * exp(-0.5*fwdVar);

                theDrift *= Maths::max(1.0, drift);
                totalVar += fwdVar;
            }

            sqrTotalVar = sqrt(totalVar);

            spotGrid[iStep].resize(numVolGrid);
            lvGrid[iStep].resize(numVolGrid);
            lvY2[iStep].resize(numVolGrid);
            // set up spot grid for lv interp
            for (j = 0; j<numVolGrid; j++)
                spotGrid[iStep][j] = fwds[iStep]*exp(-0.5*totalVar + lvStdevRange[j]*sqrTotalVar);
            //spotGrid[iStep][j] = fwds[0]*exp(lvStdevRange[j]*sqrTotalVar);

            CSliceDouble spots(&spotGrid[iStep][0], numVolGrid);
            CSliceDouble locVol(&lvGrid[iStep][0], numVolGrid);

            vol->CalcLocVol(&spots, t, &locVol, false); // no intradayinterp
        }

        // scale by sqrDt and prepare cubic spline coefficients 
	    // if step is not across div date, average the local vol
	    int divIdx=0, isDiv;
        for (iStep =0; iStep < numPathSteps-1; iStep++)
        {
            // estimate end point
            //double vol_dt = vol1 * simSqrDt[step];
            //double s2 = spot * drift * exp(vol_dt * rand*0 - 0.5*vol_dt*vol_dt);
    		
            double sqrDt = sqrt(metric->yearFrac(futurePathDates[iStep], futurePathDates[iStep+1]));

		    while( divIdx < cumDates.size() &&
			    futurePathDates[iStep] > cumDates[divIdx] )
		    {
				    divIdx++;	// this algorithm handles multiple div on the same date
		    }

		    // no local vol averaging if this or next period is across div date, 
		    isDiv = false;
		    if( divIdx < cumDates.size() )
		    {
			    if( futurePathDates[iStep] == cumDates[divIdx] ||
				    futurePathDates[iStep+1] > cumDates[divIdx] )
				    isDiv = true; // this period across div date
			    else if ( iStep < (numPathSteps-1) && ( 
				    futurePathDates[iStep+1] == cumDates[divIdx] || 
				    futurePathDates[iStep+2] > cumDates[divIdx] ))
				    isDiv = true; // next period across div date
		    }

		    if ( isDiv )
		    {
			    for(j =0; j<numVolGrid; j++)
				    lvGrid[iStep][j] *= sqrDt;
		    }
		    else
		    {
			    // prepare spline for next step for interpolation
			    spline(&*spotGrid[iStep+1].begin()-1, &*lvGrid[iStep+1].begin()-1, 
				    numVolGrid, y1, yn, &*lvY2[iStep+1].begin()-1);

			    // averaging step and step+1 vols for same grid points
			    double vol_next;
			    for(j =0; j<numVolGrid; j++)
			    {
				    // interpolate vol at current grid for next step
				    splint(&*spotGrid[iStep+1].begin()-1,
					    &*lvGrid[iStep+1].begin()-1, &*lvY2[iStep+1].begin()-1, 
					    spotGrid[iStep+1].size(), spotGrid[iStep][j], &vol_next);
				    // taking average variance and multiply by sqrt(dt)
				    lvGrid[iStep][j] = sqrDt*sqrt((lvGrid[iStep][j]*lvGrid[iStep][j]+vol_next*vol_next+lvGrid[iStep][j]*vol_next)/3.0);
			    }
            }
            // re spline after averaging and put dt in
            spline(&*spotGrid[iStep].begin()-1, &*lvGrid[iStep].begin()-1, 
                numVolGrid, y1, yn, &*lvY2[iStep].begin()-1);

    // debug testing
    #ifdef PRINT_DEBUG_GRID
            int j0;
            if (iStep==0)
            {
                fprintf(file, "Step, dt, spot_voldt, ");
                for (j0=-numVolGrid/2; j0<=numVolGrid/2; j0++)
                    fprintf(file, "%d,", j0);
            }

            fprintf(file, "\n%d, %f, spot,", iStep, sqrDt*sqrDt);
            for (j0=0; j0<numVolGrid; j0++)
                fprintf(file, "%f,", spotGrid[iStep][j0]);

            fprintf(file, "\n, , vol_dt, ");
            for (j0=0; j0<numVolGrid; j0++)
                fprintf(file, "%f,", lvGrid[iStep][j0]);
    #endif
    // debug testing
        }
    // debug testing
    #ifdef PRINT_DEBUG_GRID
        fclose(file);
    #endif

        // remove the last unused point
        spotGrid.resize(numPathSteps-1);
        lvGrid.resize(numPathSteps-1);
        lvY2.resize(numPathSteps-1);

        // here we're simulating what happens in generatePath
        // want theDrift to represent [reasonable] worst
        // case drift. A value of 2 for MAX_DRIFT_RAND_NUMBER
        // means that we catch >96% of paths
        _maxDrift = exp(sqrTotalVar*MAX_DRIFT_RAND_NUMBER);
    	
	    gridStart.resize(numPathSteps-1);
	    gridSize.resize(numPathSteps-1);
        gridDates.resize(numPathSteps-1);
        for(iStep =0; iStep < numPathSteps-1; iStep++)
	    {
		    gridStart[iStep] = log(spotGrid[iStep][0] / fwds[iStep]);
            // gridStart[iStep] = log(spotGrid[iStep][0]);
		    gridSize[iStep] = log(spotGrid[iStep][1] / spotGrid[iStep][0]);
	    }
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}


// generate LV grids for each asset vol
ILocalVolGridSP LocalVolGrid::createLocalVolGrid(const DoubleArray& fwds, 
                                                 CVolProcessedDVFConstSP vol, 
                                                 const DateTimeArray& futurePathDates,
                                                 const DateTimeArray& cumDates,
                                                 int numVolGrid, 
                                                 double stdevGridRange) {
    LocalVolGridSP grid(new LocalVolGrid(fwds, vol, futurePathDates, cumDates, numVolGrid, stdevGridRange));
    return grid;
}


void LocalVolGrid::CalcLocVol(const CSliceDouble& strikes, int iDate, CSliceDouble& locVol) {
    static const string method = "LocalVolGrid::CalcLocVol";
    
    try {
        for(int i = 0; i < strikes.size(); i++) {
            locVol[i] = interpLocVar(iDate, log(strikes[i]));
        }
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}


double LocalVolGrid::interpLocVar(int step, double logSpot) const {
    // Offset step index
    step += dateOffset;
    double lo = (logSpot - gridStart[step])/gridSize[step];
    int klo = (int) lo;
    if ( klo<0 ) 
        return lvGrid[step][0];
    else if (klo >= (spotGrid[step].size() - 1) )
        return lvGrid[step][spotGrid[step].size()-1];

    vector<double>::const_iterator ya = lvGrid[step].begin() + klo;
    return (*ya) + (lo - klo)*(*(ya+1) - *ya);
}


double LocalVolGrid::maxDrift() const {
    return _maxDrift;
}


void LocalVolGrid::rollTime(const DateTime& valueDate) {
    // Offset the date index
    int index = valueDate.numPastDates(gridDates);
    dateOffset = index;
}


DRLIB_END_NAMESPACE
