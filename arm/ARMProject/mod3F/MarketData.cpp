#include <math.h>
#include "MarketData.h"



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


CcyMarketData::CcyMarketData()
{

}

CcyMarketData::~CcyMarketData()
{

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

MarketData::MarketData()
{

}

MarketData::~MarketData()
{

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


double CcyMarketData::ZC(double T)
{
    if( T == 0 )
    {
        return 1.;
    }

    double R;
    if( T <= ZC_YearFraction.at(1) )
    {
        R = -log(ZC_Rates.at(1)) / ZC_YearFraction.at(1);
        return exp(-R*T);
    }

    double R1, R2;

    int i = 1;

    while( i < (ZC_YearFraction.entries()-2) && ZC_YearFraction.at(i+1) <= T )
    {
        i++;
    }

    R1 = -log(ZC_Rates.at(i)) / ZC_YearFraction.at(i);
    R2 = -log(ZC_Rates.at(i+1)) / ZC_YearFraction.at(i+1);
    R = R1+(T-ZC_YearFraction.at(i))*(R2-R1) / (ZC_YearFraction.at(i+1)-ZC_YearFraction.at(i));

    return exp(-R*T);
}
// CcyMarketData::ZC


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double CcyMarketData::TauxZC(double T)
{
    if( T <= ZC_YearFraction.at(1) )
    {
        return -log(ZC_Rates.at(1))/ZC_YearFraction.at(1);
    }

    double R, R1, R2;
    int i = 1;

    while( i < (ZC_YearFraction.entries()-2) && ZC_YearFraction.at(i+1) <= T )
    {
        i++;
    }

    R1 = -log(ZC_Rates.at(i)) / ZC_YearFraction.at(i);
    R2 = -log(ZC_Rates.at(i+1)) / ZC_YearFraction.at(i+1);
    R = R1+(T-ZC_YearFraction.at(i))*(R2-R1) / (ZC_YearFraction.at(i+1)-ZC_YearFraction.at(i));

    return R;
}
// CcyMarketData:TauxZC


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


double CcyMarketData::BasisZC(double T)
{
    if( T == 0 )
    {
        return 1.;
    }

    double R;

    if( T <= Basis_ZC_YearFraction.at(1) )
    {
        R=-log(Basis_ZC_Rates.at(1))/Basis_ZC_YearFraction.at(1);
        return exp(-R*T);
    }

    double R1, R2;

    int i=1;

    while( i < (Basis_ZC_YearFraction.entries()-2) && Basis_ZC_YearFraction.at(i+1) <= T )
    {
        i++;
    }

    R1 = -log(Basis_ZC_Rates.at(i)) / Basis_ZC_YearFraction.at(i);

    R2 = -log(Basis_ZC_Rates.at(i+1)) / Basis_ZC_YearFraction.at(i+1);

    R = R1+(T-Basis_ZC_YearFraction.at(i))*(R2-R1) / (Basis_ZC_YearFraction.at(i+1)-Basis_ZC_YearFraction.at(i));

    return exp(-R*T);
}
// CcyMarketData::BasisZC


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double CcyMarketData::BasisTauxZC(double T)
{
    if( T <= Basis_ZC_YearFraction.at(1) )
    {
        return -log(Basis_ZC_Rates.at(1)) / Basis_ZC_YearFraction.at(1);
    }

    double R, R1, R2;
    int i=1;

    while( i < (Basis_ZC_YearFraction.entries()-2) && Basis_ZC_YearFraction.at(i+1) <= T )
    {
        i++;
    }

    R1 = -log(Basis_ZC_Rates.at(i)) / Basis_ZC_YearFraction.at(i);

    R2 = -log(Basis_ZC_Rates.at(i+1)) / Basis_ZC_YearFraction.at(i+1);

    R = R1+(T-Basis_ZC_YearFraction.at(i))*(R2-R1) / (Basis_ZC_YearFraction.at(i+1)-Basis_ZC_YearFraction.at(i));

    return R;
}
// CcyMarketData::BasisTauxZC



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CcyMarketData::Init(double SpotDate,
                    DKMaille<double> dDates,
                    DKMaille<double> dRates,
                    DKMaille<double> dDatesNoBasis,
                    DKMaille<double> dRatesNoBasis)
{
    int i;

    ZC_Dates=dDatesNoBasis;
    ZC_YearFraction.resize(ZC_Dates.entries());
    ZC_Rates=dRatesNoBasis;
    Basis_ZC_Dates=dDates;
    Basis_ZC_YearFraction.resize(Basis_ZC_Dates.entries());
    Basis_ZC_Rates=dRates;



    for( i = 0; i < ZC_YearFraction.entries(); i++ )
    {
        ZC_YearFraction.at(i)=(ZC_Dates.at(i)-SpotDate)/365.;
    }

    for( i = 0; i < Basis_ZC_YearFraction.entries(); i++ )
    {
        Basis_ZC_YearFraction.at(i)=(Basis_ZC_Dates.at(i)-SpotDate)/365.;
    }
}
// CcyMarketData::Init


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CcyMarketData::GetConstantSigma()
{
    int i;

    ConstSigma.resize(Sigma.entries());

    for( i = 0; i < Sigma.entries(); i++)
    {
        ConstSigma.at(i) = 1.;
    }
}
// CcyMarketData::GetConstantSigma


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CcyMarketData::SetForwardRates(DKMaille<double> T)
{
    ForwardRates.resize(T.entries());
    ForwardDates.resize(T.entries());

    for( int i = 0 ; i < T.entries() - 1; i++ )
    {
        ForwardDates.at(i) = T.at(i);
        ForwardRates.at(i) = ( BasisZC(T.at(i)) / BasisZC(T.at(i+1)) - 1.) / (T.at(i+1) - T.at(i));
    }

    ForwardDates.at(ForwardDates.entries()-1) = T.at(T.entries()-1);
    ForwardRates.at(ForwardRates.entries()-1) = ForwardRates.at(ForwardRates.entries()-2);
}
// CcyMarketData::GetForwardRates


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double CcyMarketData::GetBforBondPrice( double MeanReversion,
                                        double t,
                                        double T,
                                        DKMaille<double> q,
                                        DKMaille<double> K)
{
    if( T < t || T > ForwardDates.at(ForwardDates.entries()-1) )
        throw("CcyMarketData.BondPrice : Probleme avec les dates");

    // On determine la position de t et T au milieu des ForwardDates
    int iMin, iMax;

    iMin = 1;
    while(iMin < ForwardDates.entries()-1 &&  t > ForwardDates.at(iMin) )
    {
        iMin++;
    }
    if( iMin ==  ForwardDates.entries()-1 )
        throw("CcyMarketData.BondPrice : Probleme dans la determination de iMin");

    iMax = iMin -1;

    while( iMax < ForwardDates.entries()-1 && T > ForwardDates.at(iMax +1) )
    {
        iMax++;
    }

    double B=0;

    if( iMax == iMin - 1)
    {
        B = ( q.at(iMax) * (ForwardRates.at(iMax) - K.at(iMax)) + K.at(iMax) ) * ( 1.- exp( -MeanReversion * (T-t) ) ) / MeanReversion;
    }
    else
    {
        B +=  ( q.at(iMin-1) * (ForwardRates.at(iMin-1) - K.at(iMin-1)) + K.at(iMin-1) ) * ( 1.- exp( -MeanReversion * (ForwardDates.at(iMin)-t) ) );

        for( int i = iMin; i < iMax; i++)
        {
            B += ( q.at(i) * (ForwardRates.at(i) - K.at(i)) + K.at(i) ) * ( exp( -MeanReversion * (ForwardDates.at(i)-t) ) - exp( -MeanReversion * (ForwardDates.at(i+1)-t) ) );
        }

        B += ( q.at(i) * (ForwardRates.at(iMax) - K.at(iMax)) + K.at(iMax) ) * ( exp( -MeanReversion * (ForwardDates.at(iMax)-t) ) - exp( -MeanReversion * (T-t) ) );

        B /= MeanReversion;
    }

    return B;
}
// CcyMarketData::BondPrice


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double CcyMarketData::GetInstantaneousRate(double T)
{
    int i=1;

    while( i < ForwardDates.entries() && T >= ForwardDates.at(i+1) )
    {
        i++;
    }

    return ForwardRates.at(i);
}
// CcyMarketData::GetInstantaneousRate


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CcyMarketData::Print(FILE* file)
{
    ZC_Dates.print(file);
    ZC_YearFraction.print(file);
    ZC_Rates.print(file);

    Basis_ZC_Dates.print(file);
    Basis_ZC_YearFraction.print(file);
    Basis_ZC_Rates.print(file);

    ForwardRates.print(file);
    ForwardDates.print(file);

    SigmaDates.print(file);
    Sigma.print(file);
    ConstSigma.print(file);
}
// CcyMarketData::Print

