#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include "error2.h"
#include "payoff.h"
#include "proba_utils.h"
#include "gaussian.h"
#include "student.h"
#include "dependence.h"


typedef void(*PAYOFF_FUNCTION)(CREDIT_PORTFOLIO *,
                               double*,
                               long*,
                               double*,
                               double*);
typedef struct
{
    double time;
    long index;
} TimeStruct;

/* --------------------------------------------------------------------------
// Compare
// returns 0,>0 or <0 given the comparison of TimeStruct
// this function is a parameter of qsort in TranchePayoff
*/ 
int static Compare( const void *arg1,
                    const void *arg2 )
{
    return (int) (((TimeStruct *) arg1)->time - ((TimeStruct *) arg2)->time);
}

/* --------------------------------------------------------------------------
// TranchePayoff
// returns the payoff of a tranche
// nbTimes is the number of steps of the payoff
// time is an Array[nbTime] of step time
// loss is an Array[nbTime] of the corresponding loss
*/
int TranchePayoff(
    CREDIT_PORTFOLIO *port,/* (I) portfolio spec */
    double *defaultTime,   /* (I) scenario default times */
    long   *nbTime,        /* (O) nb step in payoff function */
    double *time,          /* (O) step times */
    double *loss)          /* (O) payoff at steps */
{
    static char routine[] = "TranchePayoff";
    int status = FAILURE;
    int i = 0;
    int j = 0;
    long index;
    long nbNames =port->nbNames;
    double trancheLoss = 0.0;
    double previousTrancheLoss =0.0;
    double totalLoss = 0.0;
    double totalNotional = 0.0;
    TimeStruct *defaultTimeStruct = malloc(nbNames*sizeof(TimeStruct));
    memset(defaultTimeStruct,0,sizeof(TimeStruct));
    if(defaultTimeStruct==NULL) goto RETURN;

    for(i=0;i<nbNames;i++)
    {
        defaultTimeStruct[i].time = defaultTime[i];
        defaultTimeStruct[i].index = i;
    }
    /* sort defaultTimeStruct */
    qsort(defaultTimeStruct, (size_t) nbNames, sizeof(TimeStruct), Compare);

    /* computes total notional */
    for(i=0;i<nbNames;i++)
    {
        totalNotional += port->notional[i];
    }

    for(i=0;i<nbNames;i++)
    {
        index = defaultTimeStruct[i].index;
        previousTrancheLoss = trancheLoss;
        if(defaultTimeStruct[i].time<=
            port->nameMaturity[index])
        {
            totalLoss += port->notional[index]*(1-port->recovery[index]);
            trancheLoss = min(
                max(totalLoss - port->strike1*totalNotional, 0),
                (port->strike2 - port->strike1)*totalNotional);
            if(trancheLoss != previousTrancheLoss)
            {
                time[j] = defaultTimeStruct[i].time;
                loss[j] = trancheLoss;
                j++;
            }
        }
   }
   *nbTime = j;
   status = SUCCESS;
RETURN:
   if(defaultTimeStruct) free(defaultTimeStruct);
   if(status == FAILURE)
   {
         DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);   
   }
   return status;
}

/* ---------------------------------------------------------------------------
// ExpectedPayoff
//
*/
double *ExpectedPayoff(CREDIT_PORTFOLIO *port,
                       SCENARIO_SIM *scenario,
                       const double *samplingTime, 
                       long nbSample)
{
    static PAYOFF_FUNCTION payoffFunction = &TranchePayoff;
    long i,j,k;
    long nbPaths = scenario->nbPaths;
    long nbNames = port->nbNames;
    long nbTimes = 0;
    double weight;
    double *time = malloc(nbNames*sizeof(double));
    double *loss = malloc(nbNames*sizeof(double));
    double *expectedLoss = calloc(nbSample,sizeof(double)); 
    long samplingIdxA = 0;
    long samplingIdxB = 0;
    for(j=0;j<nbPaths;j++)
    {
        weight = scenario->weight[j];
        payoffFunction( port,
                        &scenario->defaultTime[j*nbNames],
                        &nbTimes,
                        time,
                        loss);

        /* find the first loss samplingIdx */
        for (k=0;k<nbSample;k++)
        {
                if(time[0]>samplingTime[k] && time[0]<samplingTime[k+1])
                {
                    samplingIdxA = k;
                    break;
                }
                if(k == nbSample - 1)
                {
                    samplingIdxA = k;
                }
        }
        
        for(i=0;i<nbTimes;i++)
        {
            if(i == nbTimes-1)
            {
                for(k=samplingIdxA;k<nbSample;k++)
                {
                    expectedLoss[k] += loss[i]*weight;
                }
                break;
            }

            for (k=samplingIdxA;k<nbSample;k++)
            {
                if(k == nbSample-1)
                {
                    samplingIdxB = k + 1;
                    break;
                }
                if(time[i+1]>=samplingTime[k] && time[i+1]<samplingTime[k+1])
                {
                    samplingIdxB = k;
                    break;
                }
            }

            for(k=samplingIdxA;k<samplingIdxB;k++)
            {
                if(k<nbSample)
                {
                    expectedLoss[k] += loss[i]*weight;
                }
            }
            samplingIdxA = samplingIdxB;
        }
    }        
    return expectedLoss;
}

/* ---------------------------------------------------------------------------
// UniformSimCreate
//
*/
UNIFORM_SIM UniformSimCreate(long nbPaths,
                             long nbNames,
                             int copulaType,
                             const double *beta,
                             long freedomDegree)
{
    UNIFORM_SIM us;
    long j;
    double weight = 1./ nbPaths;
    us.nbNames = nbNames;
    us.nbPaths = nbPaths;
    us.uniformDeviate = malloc(nbPaths*nbNames*sizeof(double));
    us.weight = malloc(nbPaths*sizeof(double));
    if(copulaType == 0)
    {
        GaussianCopulatedUniformDeviates(   us.uniformDeviate,
                                            nbNames,
                                            nbPaths,
                                            beta);
    }
    else if(copulaType == 1)
    {
        StudentCopulatedUniformDeviates(us.uniformDeviate,
                                        nbNames,
                                        nbPaths,
                                        beta,
                                        freedomDegree);
    }
    else
    {
        DependenceCopulatedUniformDeviates(us.uniformDeviate,
                                        nbNames,
                                        nbPaths);
    }

    for(j=0;j<nbPaths;j++)
    {
        us.weight[j] = weight;
    }

    return us;
}

/* ---------------------------------------------------------------------------
// UniformSimFree
//
*/
void UniformSimFree(UNIFORM_SIM *us)
{
    if(us->uniformDeviate) free(us->uniformDeviate);
    if(us->weight) free(us->weight);
    us->uniformDeviate = NULL;
    us->weight = NULL;
}


/* ---------------------------------------------------------------------------
// SpreadCurveCreate
//
*/
SPREAD_CURVE SpreadCurveCreate(double *time, double *spread, long nbTimes)
{
    static char routine[] = "SpreadCurveCreate";
    int status = FAILURE;
    SPREAD_CURVE sc;
    sc.spread = malloc(nbTimes*sizeof(double));
    if(sc.spread==NULL) goto RETURN;
    sc.time = malloc(nbTimes*sizeof(double));
    if(sc.time==NULL) goto RETURN;
    sc.nbTimes = nbTimes;
    sc.spread = spread;
    sc.time = time;
    status = SUCCESS;
RETURN:
    if(status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);    
    }
    return sc;
}

/* ---------------------------------------------------------------------------
// SpreadCurveFree
//
*/
void SpreadCurveFree(SPREAD_CURVE *sc)
{
    if(sc->spread) free(sc->spread);
    if(sc->time) free(sc->time);
    sc->spread = NULL;
    sc->time = NULL;
}

/* ---------------------------------------------------------------------------
// ScenarioSimCreate
//
*/
SCENARIO_SIM ScenarioSimCreate(UNIFORM_SIM *us,
                               SPREAD_CURVE *sc)
{
    SCENARIO_SIM ss;
    long nbPaths = us->nbPaths;
    long nbNames = us->nbNames;
    long i,j;
    double *time;
    double *spread;
    long nbTimes;
    ss.nbNames = nbNames;
    ss.nbPaths = nbPaths;
    ss.defaultTime = malloc(nbPaths*nbNames*sizeof(double));
    ss.weight = malloc(nbPaths*sizeof(double));
    for(i=0;i<nbNames;i++)
    {
        time = sc[i].time;
        spread= sc[i].spread;
        nbTimes = sc[i].nbTimes;
        for(j=0;j<nbPaths;j++)
        {
           ss.defaultTime[i+j*nbNames] =
                DefaultTimeFromUniform( time,
                                        spread,
                                        nbTimes,
                                        us->uniformDeviate[i+j*nbNames]);
        }
    }

    memcpy(ss.weight,us->weight, nbPaths*sizeof(double));
    return ss;
}

/* ---------------------------------------------------------------------------
// ScenarioSimCreate2
//
*/
SCENARIO_SIM ScenarioSimCreate2(UNIFORM_SIM *us1,
                                UNIFORM_SIM *us2,
                                SPREAD_CURVE *sc,
                                const double *alpha /* [nbNames] */)
{
    SCENARIO_SIM ss;
    long nbPaths = us1->nbPaths;
    long nbNames = us1->nbNames;
    long i,j,l;
    double *time;
    double *spread;
    long nbTimes;

    ss.nbNames = nbNames;
    ss.nbPaths = nbPaths;
    ss.defaultTime = malloc(nbPaths*nbNames*sizeof(double));
    ss.weight = malloc(nbPaths*sizeof(double));
    
    for(i=0;i<nbNames;i++)
    {
        time = sc[i].time;
        nbTimes = sc[i].nbTimes;

        for(j=0;j<nbPaths;j++)
        {
                // 1st copula
                spread = sc[i].spread;
                for(l=0;l<nbTimes;l++)
                {
                    spread[l] *= alpha[i];
                }

                ss.defaultTime[i+j*nbNames] =
                    DefaultTimeFromUniform( time,
                                            spread,
                                            nbTimes,
                                            us1->uniformDeviate[i+j*nbNames]);

                // 2nd copula
                spread = sc[i].spread;
                for(l=0;l<nbTimes;l++)
                {
                    spread[l] *= 1.- alpha[i];
                }

                ss.defaultTime[i+j*nbNames] = min(
                    ss.defaultTime[i+j*nbNames],
                    DefaultTimeFromUniform( time,
                                            spread,
                                            nbTimes,
                                            us2->uniformDeviate[i+j*nbNames]));
        }
    }
    
    memcpy(ss.weight,us1->weight,nbPaths*sizeof(double));
    return ss;
}

/* ---------------------------------------------------------------------------
// ScenarioSimFree
//
*/
void ScenarioSimFree(SCENARIO_SIM *ss)
{
    if(ss->defaultTime) free(ss->defaultTime);
    if(ss->weight) free(ss->weight);
    ss->defaultTime = NULL;
    ss->weight = NULL;
}

/* --------------------------------------------------------------------------
//  CumulatedSpreadToFlatForward
//
*/
double *CumulatedSpreadToFlatForward(double *time, double *cumSpread, long nbTimes)
{
    long i = 0;
    double *flatForward = malloc(nbTimes*sizeof(double));
    flatForward[0] = cumSpread[0];
    for(i=1;i<nbTimes;i++)
    {
        flatForward[i]  = time[i] * cumSpread[i] -
                                time[i-1] * cumSpread[i-1];
        flatForward[i] /= (time[i]-time[i-1]);
    }
    return flatForward;
}

/* --------------------------------------------------------------------------
//  CopulaProduct
//
*/
UNIFORM_SIM CopulaProduct(UNIFORM_SIM u_s1,
                          UNIFORM_SIM u_s2,
                          const double *alpha)
{
    UNIFORM_SIM u_s;
    long nbNames = u_s1.nbNames;
    long nbPaths = u_s1.nbPaths;
    double a,b,alp;
    long i = 0;
    long j = 0;
    u_s.uniformDeviate = malloc(nbNames*nbPaths*sizeof(double));
    u_s.weight = malloc(nbPaths*sizeof(double));
    for(i=0;i<nbNames;i++)
    {
        alp = alpha[i];
        a = 1./alp;
        b = 1./(1-alp);
        for(j=0;j<nbPaths;j++)
        {
            u_s.uniformDeviate[i+j*nbNames] = max(
                    pow(u_s1.uniformDeviate[i+j*nbNames],a) ,
                    pow(u_s2.uniformDeviate[i+j*nbNames],b));
        }
    }

    for(j=0;j<nbPaths;j++)
    {
        //// weights?
    }
  
    return u_s;
}

/* -----------------------------------------------------------------------------
// DefaultTimeFromUniform
// this function returns a default time
// given a triggerLevel in [0,1],
// a flat forward spread 
// spread[i] = constant spread between time[i]  and time[i+1]
*/
double DefaultTimeFromUniform(const double *time,
                              const double *spread,
                              long nbTimes,
                              double triggerLevel)
{
    double sum = 0.0;
    double logTriggerLevel = log(triggerLevel);
    long i = 0;
    if(logTriggerLevel==0.0)
    {
        return 0.;
    }

    if((sum+=spread[0]*(time[0]))>-logTriggerLevel)
    {
        return -(logTriggerLevel+sum - spread[0]*(time[i]))/
                                                        spread[i];        
    }
    
    for(i=0;i<nbTimes-1;i++)
    {
        if((sum+=spread[i]*(time[i]-time[i-1]))>-logTriggerLevel)
        {
            return -(logTriggerLevel+sum - spread[i]*(time[i]-time[i-1]))/
                                                        spread[i] + time[i-1];        
        }
    }
    return -(logTriggerLevel+sum)/spread[nbTimes-1] + time[nbTimes-1];
}
