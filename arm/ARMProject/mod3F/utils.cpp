
#include "math.h"
#include "utils.h"



// Fonction de répartition de la  loi normale
double CumNormDist(double Z)
{

    double P0, P1, P2, P3, P4, P5, P6, Q0, Q1, Q2, Q3, Q4, Q5, Q6, Q7, P, EXPNTL, CUTOFF, ROOTPI, ZABS;

    P0 = 220.2068679123761;
    P1 = 221.2135961699311;
    P2 = 112.0792914978709;
    P3 = 33.91286607838300;
    P4 = 6.373962203531650;
    P5 = .7003830644436881;
    P6 = .03526249659989109;

    Q0 = 440.4137358247522;
    Q1 = 793.8265125199484;
    Q2 = 637.3336333788311;
    Q3 = 296.5642487796737;
    Q4 = 86.78073220294608;
    Q5 = 16.06417757920695;
    Q6 = 1.755667163182642;
    Q7 = .08838834764831844;


    ROOTPI = 2.506628274631001;
    CUTOFF = 7.071067811865475;

    ZABS = fabs(Z);

    if (ZABS > 37.)
    {
        P = 0.;
    }
    else
    {

        EXPNTL = exp( -ZABS*ZABS/2.);
        if (ZABS < CUTOFF)
        {
            P = EXPNTL*((((((P6*ZABS + P5)*ZABS + P4)*ZABS + P3)*ZABS + P2)*ZABS + P1)*ZABS + P0)
                / (((((((Q7*ZABS + Q6)*ZABS + Q5)*ZABS + Q4)*ZABS + Q3)*ZABS + Q2)*ZABS + Q1)*ZABS + Q0);
        }
        else
        {
            P = EXPNTL/(ZABS + 1./(ZABS + 2./(ZABS + 3./(ZABS + 4./(ZABS + 0.65)))))/ROOTPI;
        }
    }
    if ( Z > 0 ) P = 1 - P;
    return P;

}

double NormDist(double x, double x0, double Sigma)
{
    return 1./(sqrt(2.*Pi)*Sigma)*exp(-0.5*pow((x-x0)/Sigma,2));
}


double Integration(double (*func)(double x, double *Y), double *Y, double LB, double HB, int NbSteps)
{
    double integrale=0.;

    double f1,f2;
    f1=(*func)(LB,Y);
    double dx=1.*(HB-LB)/NbSteps;

    for(int i=0;i<(int) NbSteps;i++)
    {
        f2=(*func)(LB+(i+1.)/NbSteps*(HB-LB),Y);
        integrale+=(f1+f2)/2.*dx;
        f1=f2;
    }

    return integrale;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Create_Strip_Analytics(DKMaille<double> &dStrip,
                            DKMaille<double> &dStripIntervalles,
                            DKMaille<double> NoticeDates,
                            int N1,
                            int N2,
                            int N3,
                            double dOptimalDate,
                            double LastDate,
                            int* NbPasTotal,
                            int* NbPasBeforeLastNotice)
{
    int ui, uiSlices;
    dStrip.resize(1);
    dStrip.at(0) = 0.;

    // determination du nombre de pas avant la premiere notice
    uiSlices = (int) ( N1 * NoticeDates.at(0) );
    if( uiSlices < 4 || uiSlices < N1 )
    {
        uiSlices = FMAX( 4, N1 );
    }

    // on insert uiSlices avant la premiere notice
    for( ui = 1; ui <= uiSlices; ui++ )
    {
        dStrip.insert( (double)(ui) / (double)(uiSlices) * NoticeDates.at(0) );
    }


    // determination de la NoticeDate la plus proche de l'OptimalDate
    int Optimal=-1;
    for( ui = 0; ui < NoticeDates.entries(); ui++ )
    {
        if( NoticeDates.at(ui) > dOptimalDate )
        {
            if( ui == 0 )
            {
                Optimal = 0;
            }
            else
            {
                if( ( NoticeDates.at(ui) - dOptimalDate ) > ( dOptimalDate - NoticeDates.at(ui-1) ) )
                {
                    Optimal = ui-1;
                }
                else
                {
                    Optimal = ui;
                }
            }
            ui = NoticeDates.entries();
        }

    }

    // si Optimal est encore a -1  c que dOptimalDate est superieure a la derniere notice
    if( Optimal == -1 )
    {
        Optimal = NoticeDates.entries()-1;
    }



    // On insert les slices jusqu'a la derniere notice
    for( ui = 0; ui < Optimal; ui++ )
    {
        uiSlices = (int) ( (double) (N2) * ( NoticeDates.at(ui+1) - NoticeDates.at(ui) ) + 0.7 );

        for( int j = 1; j <= uiSlices; j++ )
        {
            dStrip.insert( NoticeDates.at(ui) + (double)(j) / (double)(uiSlices)
                           * ( NoticeDates.at(ui+1) - NoticeDates.at(ui) ) );
        }
    }

    for( ui = Optimal; ui < NoticeDates.entries()-1; ui++ )
    {
        uiSlices = (int) ( (double) (N3) * ( NoticeDates.at(ui+1) - NoticeDates.at(ui) ) + 0.7 );

        for( int j = 1; j <= uiSlices; j++ )
        {
            dStrip.insert( NoticeDates.at(ui) + (double)(j) / (double)(uiSlices)
                           * ( NoticeDates.at(ui+1) - NoticeDates.at(ui) ) );
        }
    }

    // dStrip est remplit jusqu'a la derniere notice
    *NbPasBeforeLastNotice = dStrip.entries() - 1;

    // On remplit jusqu'a la LastDate
    if(Optimal!=NoticeDates.entries()-1)
    {
        uiSlices = (int) ( N3 * ( LastDate - NoticeDates.at(NoticeDates.entries()-1) ) + 0.7 );
    }
    else
    {
        uiSlices = (int) ( N2 * ( LastDate - NoticeDates.at(NoticeDates.entries()-1) ) + 0.7 );
    }

    for( ui = 1; ui <= uiSlices; ui++ )
    {
        dStrip.insert( NoticeDates.at(NoticeDates.entries()-1)
                       + (double)(ui) / (double)(uiSlices) * (LastDate-NoticeDates.at(NoticeDates.entries()-1)) );
    }

    // dStrip est maintenant remplit jusqu'a la fin
    *NbPasTotal = dStrip.entries() - 1;



    dStripIntervalles.resize( *NbPasTotal );
    for( ui = 1; ui < dStrip.entries(); ui++ )
    {
        dStripIntervalles.at(ui-1) = dStrip.at(ui) - dStrip.at(ui-1);
    }

}
// Create_Strip_Analytics



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double ZC_interpole(double T1,
                    double T2,
                    double ZC1,
                    double ZC2,
                    double T)
{
    if( T < T1 || T > T2 )
        throw("ZC_interpole : Interpolation au dela des bornes");
    double R1, R2, R;
    R1 = -log(ZC1)/T1;
    R2 = -log(ZC2)/T2;

    R = R1 + (R2 - R1) * (T - T1) / (T2 - T1);

    return exp(- R * T );
}
// ZC_interpole

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double FromRateToX(double r, double q, double K)
{
    if(q != 0.)
    {
        return log(1. + q * ( r / K - 1.) ) / q;
    }
    else if( K != 0 )
    {
        return r / K - 1.;
    }
    else
    {
        return r;
    }
    return 0.;
}
// FromRateToX

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double FromXToRate(double X, double q, double K)
{
    if(q != 0.)
    {
        double r = K * ( 1. + ( exp( q * X ) - 1.) / q );
        if( r > 5.)
        {
            return 5.;
        }
        else
        {
            return r;
        }
    }
    else if( K != 0. )
    {
        double r = K * (1 + X);
        if( r > 5.)
        {
            return 5.;
        }
        else
        {
            return r;
        }
    }
    else
    {
        return X;
    }
    return 0.;
}
// FromXToRate


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double FromXToRate_deriv(double X, double q, double K)
{
    if(q != 0.)
    {
        double Xmax = log(q * (5./K - 1.) + 1) / q;
        if(X <= Xmax)
        {
            return K *  exp( q * X );
        }
        else
        {
            return 0.;
        }
    }
    else if( K != 0.)
    {
        double Xmax = (5./K - 1.);
        if(X <= Xmax)
        {
            return K;
        }
        else
        {
            return 0.;
        }
    }
    else
    {
        return 1.;
    }
    return 0.;
}
// FromXToRate_deriv


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double BondPrice(   double A,
                    double MeanReversion,
                    double t,
                    double T,
                    double r)
{
    if( T < t )
        throw ("BondPrice : T must be superior than t");

    double B;
    B = (1.-exp(-MeanReversion*(T-t)))/MeanReversion;

    return A * exp( - B * r);

}
// BondPrice

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

