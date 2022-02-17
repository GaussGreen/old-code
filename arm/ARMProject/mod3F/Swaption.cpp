
#include "Swaption.h"
#include "hw_vfdk_analytics.h"
#include "hw_vfdk_LDHD_lattice.h"
#include "utils.h"
#include "math.h"
#include "dk_utils.h"


FILE* fic;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


Swaption::Swaption()
{
    CallPut = -1.;
    PriceZCinTree = false;
    QModel = false;
}

Swaption::~Swaption()
{

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


DKMaille<double>& SwaptionPrice( double dSpotDate,
                                 DKMaille<double> &dDates,
                                 DKMaille<double> &dRates,
                                 DKMaille<double> &dRatesNoBasis,
                                 DKMaille<double> &dStdDevX,
                                 DKMaille<double> &dStdDevY,
                                 DKMaille2D<double> &dStdDevZ,
                                 double dMeanReversion,
                                 double dSmileParameter1,
                                 double dSmileParameter2,
                                 double dIsSwaptionCalibrationWithBasis,
                                 DKMaille2D<double> &dBoosterData,
                                 double dNumTimeLinesBeforeFirstNotice,
                                 double dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
                                 double dNumTimeLinesPerYear,
                                 double dOptimal)
{
    //fic = fopen("C:\\Documents and Settings\\spannetier\\desktop\\DBGFIC.txt", "w");

    //Creation de l'objet Swaption
    Swaption swaption;

    int NbRows=dBoosterData.rows();

    DKMaille<double> dNoticeDates;
    DKMaille<double> dFloatResetDates(NbRows);
    DKMaille<double> dFixedPaymentDates(NbRows);
    DKMaille<double> dFloatPaymentDates(NbRows);
    DKMaille<double> dIndexStartDates(NbRows);
    DKMaille<double> dIndexEndDates(NbRows);
    DKMaille<double> dFixedAccrualBasis(NbRows);
    DKMaille<double> dFloatAccrualBasis(NbRows);
    DKMaille<double> dFixedCoupon(NbRows);
    DKMaille<double> dBasisSpread(NbRows);

    for(int i=0;i<NbRows;i++)
    {
        if((dBoosterData.at(i,0)-dSpotDate)/365.>0)
            dNoticeDates.insert((dBoosterData.at(i,0)-dSpotDate)/365.);
        dFloatResetDates.at(i)=(dBoosterData.at(i,13)-dSpotDate)/365.;
        dFixedPaymentDates.at(i)=(dBoosterData.at(i,2)-dSpotDate)/365.;
        dFloatPaymentDates.at(i)=(dBoosterData.at(i,5)-dSpotDate)/365.;
        dIndexStartDates.at(i)=(dBoosterData.at(i,3)-dSpotDate)/365.;
        dIndexEndDates.at(i)=(dBoosterData.at(i,4)-dSpotDate)/365.;
        dFixedAccrualBasis.at(i)=dBoosterData.at(i,10);
        dFloatAccrualBasis.at(i)=dBoosterData.at(i,11);
        dFixedCoupon.at(i)=dBoosterData.at(i,12);
        dBasisSpread.at(i)=dBoosterData.at(i,6);
    }

    double LastMaturity=(FMAX(dBoosterData.at(dBoosterData.rows()-1,2),dBoosterData.at(dBoosterData.rows()-1,5))-dSpotDate)/365.;



    swaption.ccyMarketData.Init(dSpotDate,
                                dDates,
                                dRates,
                                dDates,
                                dRatesNoBasis);


    // Bootstrappting du Sigma de t
    DKMaille<double> dSigmaDates;
    DKMaille<double> dSigma;




    GetBootstrappedVolForQModelInTree(  dSigmaDates,
                                        dSigma,
                                        dNoticeDates,
                                        swaption.ccyMarketData.ZC_Dates,
                                        swaption.ccyMarketData.ZC_Rates,
                                        swaption.ccyMarketData.Basis_ZC_Rates,
                                        dStdDevX,
                                        dStdDevY,
                                        dStdDevZ,
                                        dSpotDate,
                                        LastMaturity,
                                        dNumTimeLinesBeforeFirstNotice,
                                        dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
                                        dNumTimeLinesPerYear,
                                        dOptimal,
                                        dMeanReversion,
                                        dSmileParameter1);





    swaption.Init(  dSpotDate,
                    dNoticeDates,
                    dNumTimeLinesBeforeFirstNotice,
                    dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
                    dNumTimeLinesPerYear,
                    dOptimal,
                    LastMaturity,
                    dMeanReversion,
                    dSmileParameter1,
                    dFloatResetDates,
                    dFixedPaymentDates,
                    dFloatPaymentDates,
                    dIndexStartDates,
                    dIndexEndDates,
                    dFixedAccrualBasis,
                    dFloatAccrualBasis,
                    dFixedCoupon,
                    dBasisSpread);

    swaption.SetVol(dSigmaDates,
                    dSigma);

    swaption.CreateTree();

    static DKMaille<double> Price;

    Price = swaption.Price((int) (dSmileParameter2));

    //fclose(fic);

    return Price;
}
// SwaptionPrice

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Swaption::Init( double dSpotDate,
                DKMaille<double> dNoticeDates,
                double dNumTimeLinesBeforeFirstNotice,
                double dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
                double dNumTimeLinesPerYear,
                double dOptimal,
                double LastMaturity,
                double dMeanReversion,
                double smileParameter,
                DKMaille<double> dFloatResetDates,
                DKMaille<double> dFixedPaymentDates,
                DKMaille<double> dFloatPaymentDates,
                DKMaille<double> dIndexStartDates,
                DKMaille<double> dIndexEndDates,
                DKMaille<double> dFixedAccrualBasis,
                DKMaille<double> dFloatAccrualBasis,
                DKMaille<double> dFixedCoupon,
                DKMaille<double> dBasisSpread)
{
    SpotDate = dSpotDate;
    MeanReversion = dMeanReversion;
    NbNotices=dNoticeDates.entries();
    NbCoupons=dFloatPaymentDates.entries();


    NoticeDates = dNoticeDates;
    FloatResetDates = dFloatResetDates;
    FixedPaymentDates = dFixedPaymentDates;
    FloatPaymentDates = dFloatPaymentDates;
    IndexStartDates = dIndexStartDates;
    IndexEndDates = dIndexEndDates;
    FixedAccrualBasis = dFixedAccrualBasis;
    FloatAccrualBasis = dFloatAccrualBasis;
    FixedCoupon = dFixedCoupon;
    BasisSpread= dBasisSpread;


    Create_Strip_Analytics(	T,
                            dT,
                            dNoticeDates,
                            dNumTimeLinesBeforeFirstNotice,
                            dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
                            dNumTimeLinesPerYear,
                            dOptimal,
                            LastMaturity,
                            &NbPasTotal,
                            &NbPasBeforeLastNotice);


    // determination des parametres de smile
    if( smileParameter != 0.)
    {
        QModel = true;
        ccyMarketData.SetForwardRates(T);

    }

    Setq(smileParameter);
    SetK();
}
// Swaption::Init

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Swaption::ModifyBoosterData(DKMaille<double> dFloatResetDates,
                            DKMaille<double> dFixedPaymentDates,
                            DKMaille<double> dFloatPaymentDates,
                            DKMaille<double> dIndexStartDates,
                            DKMaille<double> dIndexEndDates,
                            DKMaille<double> dFixedAccrualBasis,
                            DKMaille<double> dFloatAccrualBasis,
                            DKMaille<double> dFixedCoupon,
                            DKMaille<double> dBasisSpread)
{
    FloatResetDates.clear();
    FixedPaymentDates.clear();
    FloatPaymentDates.clear();
    IndexStartDates.clear();
    IndexEndDates.clear();
    FixedAccrualBasis.clear();
    FloatAccrualBasis.clear();
    FixedCoupon.clear();
    BasisSpread.clear();

    NbCoupons=dFloatPaymentDates.entries();


    FloatResetDates = dFloatResetDates;
    FixedPaymentDates = dFixedPaymentDates;
    FloatPaymentDates = dFloatPaymentDates;
    IndexStartDates = dIndexStartDates;
    IndexEndDates = dIndexEndDates;
    FixedAccrualBasis = dFixedAccrualBasis;
    FloatAccrualBasis = dFloatAccrualBasis;
    FixedCoupon = dFixedCoupon;
    BasisSpread= dBasisSpread;
}
// Swaption::ModifyBoosterData


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Swaption::SetVol(   DKMaille<double> dSigmaDates,
                    DKMaille<double> dSigma)
{
    // On interpole les Sigma aux dates de l'arbre
    TransformVolatilitiesSingle(dSigmaDates,T,dSigma);

    ccyMarketData.SigmaDates.clear();
    ccyMarketData.Sigma.clear();

    ccyMarketData.SigmaDates = T;
    ccyMarketData.Sigma = dSigma;
}
// SetVol

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Swaption::CreateTree()
{
    double X0 = 0;

    // Determination de la volatilite moyenne pour l'arbre
    ccyMarketData.GetConstantSigma();

    // Determination de la MeanReversion ajustee pour se placer en Vol constante
    SetMeanReversionForTree();

    // si PriceZCinTree est false on ne construit l'arbre que jusqu'a la derniere notice
    // sinon on le construit jusqu'a la derniere date
    if( PriceZCinTree )
    {
        tree.Init(X0, NbPasTotal, T, ccyMarketData.ConstSigma, MeanReversionInTree);
    }
    else
    {
        tree.Init(X0, NbPasBeforeLastNotice, T, ccyMarketData.ConstSigma, MeanReversionInTree);
    }

    tree.SetMultiplicateurInConstVolCase(ccyMarketData.Sigma);


    tree.SetArrowDebreu(FromXToRate, FromXToRate_deriv, ccyMarketData, q, K);

    if( PriceZCinTree )
    {
        tree.SetZCinTree();
    }
    else
    {
        CalibrateZCinTree();
    }
}
// Swaption::CreateTree


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Swaption::CleanTree()
{
    ccyMarketData.ConstSigma.clear();

    MeanReversionInTree.clear();

    tree.Clean();
}
// Swaption::CleanTree


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Swaption::CalibrateZCinTree()
{
    int i, j, FirstResetDateAfterNoticeIndex;
    double NoticeDate, r0;
    Node* nodes;

    A_FixedPayment.resize(NbNotices, NbCoupons);
    A_FloatPayment.resize(NbNotices, NbCoupons);
    A_IndexStart.resize(NbNotices, NbCoupons);
    A_IndexEnd.resize(NbNotices, NbCoupons);

    if(QModel)
    {
        B_FixedPayment.resize(NbNotices, NbCoupons);
        B_FloatPayment.resize(NbNotices, NbCoupons);
        B_IndexStart.resize(NbNotices, NbCoupons);
        B_IndexEnd.resize(NbNotices, NbCoupons);
    }

    int uiNotice=0;
    int uiSlice=0;
    while ( uiNotice < NbNotices )
    {
        if( fabs(T.at(uiSlice)-NoticeDates.at(uiNotice))<1.e-7  )
        {
            tree.slices[uiSlice].isNotice=true;
            NoticeDate = NoticeDates.at(uiNotice);
            FirstResetDateAfterNoticeIndex = GetFirstResetDateAfterNotice(uiNotice);
            nodes=tree.slices[uiSlice].nodes;

            for( i = FirstResetDateAfterNoticeIndex; i < NbCoupons; i++ )
            {
                A_FixedPayment.at(uiNotice,i) = 0;
                A_FloatPayment.at(uiNotice,i) = 0;
                A_IndexStart.at(uiNotice,i) = 0;
                A_IndexEnd.at(uiNotice,i) = 0;

                if(QModel)
                {

                    B_FixedPayment.at(uiNotice,i) = ccyMarketData.GetBforBondPrice( MeanReversion,
                                                    NoticeDate,
                                                    FixedPaymentDates.at(i),
                                                    q,
                                                    K);
                    B_FloatPayment.at(uiNotice,i) = ccyMarketData.GetBforBondPrice( MeanReversion,
                                                    NoticeDate,
                                                    FloatPaymentDates.at(i),
                                                    q,
                                                    K);
                    B_IndexStart.at(uiNotice,i) = ccyMarketData.GetBforBondPrice(   MeanReversion,
                                                  NoticeDate,
                                                  IndexStartDates.at(i),
                                                  q,
                                                  K);
                    B_IndexEnd.at(uiNotice,i) = ccyMarketData.GetBforBondPrice( MeanReversion,
                                                NoticeDate,
                                                IndexEndDates.at(i),
                                                q,
                                                K);
                }
            }


            for( j = 0; j < tree.slices[uiSlice].nbNodes; j++ )
            {
                for( i = FirstResetDateAfterNoticeIndex; i < NbCoupons; i++ )
                {
                    if( QModel )
                    {
                        r0 = q.at(uiSlice) * (ccyMarketData.GetInstantaneousRate(T.at(uiSlice)) - K.at(uiSlice)) + K.at(uiSlice);

                        A_FixedPayment.at(uiNotice,i) += nodes[j].ArrowDebreu * exp( - B_FixedPayment.at(uiNotice,i) * nodes[j].r / r0);

                        A_FloatPayment.at(uiNotice,i) += nodes[j].ArrowDebreu * exp( - B_FloatPayment.at(uiNotice,i) * nodes[j].r / r0);

                        A_IndexStart.at(uiNotice,i) += nodes[j].ArrowDebreu * exp( - B_IndexStart.at(uiNotice,i) * nodes[j].r / r0);

                        A_IndexEnd.at(uiNotice,i) += nodes[j].ArrowDebreu * exp( - B_IndexEnd.at(uiNotice,i) * nodes[j].r / r0);
                    }
                    else
                    {
                        A_FixedPayment.at(uiNotice,i) += nodes[j].ArrowDebreu*
                                                         BondPrice(1.,MeanReversion, NoticeDate, FixedPaymentDates.at(i), nodes[j].r);
                        A_FloatPayment.at(uiNotice,i) += nodes[j].ArrowDebreu*
                                                         BondPrice(1.,MeanReversion, NoticeDate, FloatPaymentDates.at(i), nodes[j].r);
                        A_IndexStart.at(uiNotice,i) += nodes[j].ArrowDebreu*
                                                       BondPrice(1.,MeanReversion, NoticeDate, IndexStartDates.at(i), nodes[j].r);
                        A_IndexEnd.at(uiNotice,i) += nodes[j].ArrowDebreu*
                                                     BondPrice(1.,MeanReversion, NoticeDate, IndexEndDates.at(i), nodes[j].r);
                    }
                }
            }

            for( i = FirstResetDateAfterNoticeIndex; i < NbCoupons; i++ )
            {
                A_FixedPayment.at(uiNotice,i) = ccyMarketData.BasisZC(FixedPaymentDates.at(i)) / A_FixedPayment.at(uiNotice,i);
                A_FloatPayment.at(uiNotice,i) = ccyMarketData.BasisZC(FloatPaymentDates.at(i)) / A_FloatPayment.at(uiNotice,i);
                A_IndexStart.at(uiNotice,i) = ccyMarketData.BasisZC(IndexStartDates.at(i)) / A_IndexStart.at(uiNotice,i);
                A_IndexEnd.at(uiNotice,i) = ccyMarketData.BasisZC(IndexEndDates.at(i)) / A_IndexEnd.at(uiNotice,i);
            }

            // on est tombe sur une noticedate donc on incremente
            uiNotice++;
        }

        // on incremente l'indice des slices
        uiSlice++;
    }
}
// Swaption::CalibrateZCinTree

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

DKMaille<double>& Swaption::Price(int PriceOnWhatNotice)
{
    static DKMaille<double> Price;

    if(PriceZCinTree)
    {
        Price = PriceInTree(PriceOnWhatNotice);
    }
    else
    {
        Price = PriceInTreeWithAnalytics(PriceOnWhatNotice);
    }
    return Price;
}
//Price

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

DKMaille<double>& Swaption::StraddlePrice(int PriceOnWhatNotice)
{
    static DKMaille<double> Price;

    DKMaille<double> Price1;
    DKMaille<double> Price2;


    if(PriceZCinTree)
    {
        Price1 = PriceInTree(PriceOnWhatNotice);
        CallPut *= -1;
        Price2 = PriceInTree(PriceOnWhatNotice);
        CallPut *= -1;
    }
    else
    {
        Price1 = PriceInTreeWithAnalytics(PriceOnWhatNotice);
        CallPut *= -1;
        Price2 = PriceInTreeWithAnalytics(PriceOnWhatNotice);
        CallPut *= -1;
    }

    Price.resize(Price1.entries());

    for( int i = 0 ; i < Price.entries(); i++ )
    {
        Price.at(i) = Price1.at(i) + Price2.at(i);
    }

    return Price;
}
//StraddlePrice

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


DKMaille<double>& Swaption::PriceInTreeWithAnalytics(int PriceOnWhatNotice)
{
    if(PriceOnWhatNotice > NbNotices )
    {
        throw("Swaption.PriceInTreeWithAnalytics : PriceOnWhatNotice est superieur au nombre de notice");
    }

    int i, j, FirstResetDateAfterNoticeIndex, iNoticeRestantes;
    double NoticeDate, r0;
    Node* nodes;


    double SWAP, SWO, SWOi, SWOB, SWOBi, discount1, discount2, discount3;
    int uiNotice, uiSlice;
    SWO = 0.;
    SWAP = 0.;
    SWOB = 0.;


    // On calcule le prix de l'europeenne sur la slice PriceOnWhatNotice si PriceOnWhatNotice est > 0
    if( PriceOnWhatNotice > 0)
    {

        // on determine la slice de la notice numero PriceOnWhatNotice
        uiNotice = -1;
        uiSlice = 0;

        while( uiNotice + 1 < PriceOnWhatNotice && uiSlice < tree.nbSlices - 1)
        {
            uiSlice++;
            if( tree.slices[uiSlice].isNotice )
            {
                uiNotice++;
            }
        }

        NoticeDate = NoticeDates.at(uiNotice);
        FirstResetDateAfterNoticeIndex = GetFirstResetDateAfterNotice(uiNotice);
        nodes = tree.slices[uiSlice].nodes;

        if( fabs(NoticeDate-T.at(uiSlice)) > 1.e-7)
            throw ("Swaption.PriceInTreeWithAnalytics : on calcule un prix sur une slice qui n'est pas une notice date");


        for( j = 0; j < tree.slices[uiSlice].nbNodes; j++ )
        {
            SWOi = 0.;
            for( i = FirstResetDateAfterNoticeIndex; i < NbCoupons; i++ )
            {
                // FixedCoupon
                if( QModel )
                {
                    r0 = q.at(uiSlice) * (ccyMarketData.GetInstantaneousRate(T.at(uiSlice)) - K.at(uiSlice) ) + K.at(uiSlice);

                    discount3 = A_FixedPayment.at(uiNotice,i) * exp( - B_FixedPayment.at(uiNotice,i) * nodes[j].r / r0 );
                }
                else
                {
                    discount3 = BondPrice(  A_FixedPayment.at(uiNotice,i),
                                            MeanReversion,
                                            NoticeDate,
                                            FixedPaymentDates.at(i),
                                            nodes[j].r);
                }

                SWOi -= FixedCoupon.at(i)*FixedAccrualBasis.at(i)*discount3;

                // Float Coupon
                if( QModel )
                {
                    r0 = q.at(uiSlice) * (ccyMarketData.GetInstantaneousRate(T.at(uiSlice)) - K.at(uiSlice) ) + K.at(uiSlice);

                    discount3 = A_FloatPayment.at(uiNotice,i) * exp( - B_FloatPayment.at(uiNotice,i) * nodes[j].r / r0 );

                    discount1 = A_IndexStart.at(uiNotice,i) * exp( - B_IndexStart.at(uiNotice,i) * nodes[j].r / r0);

                    discount2 = A_IndexEnd.at(uiNotice,i) * exp( - B_IndexEnd.at(uiNotice,i) * nodes[j].r / r0);
                }
                else
                {
                    discount3 = BondPrice(  A_FloatPayment.at(uiNotice,i),
                                            MeanReversion,
                                            NoticeDate,
                                            FloatPaymentDates.at(i),
                                            nodes[j].r);
                    discount1 = BondPrice(  A_IndexStart.at(uiNotice,i),
                                            MeanReversion,
                                            NoticeDate,
                                            IndexStartDates.at(i),
                                            nodes[j].r);
                    discount2 = BondPrice(  A_IndexEnd.at(uiNotice,i),
                                            MeanReversion,
                                            NoticeDate,
                                            IndexEndDates.at(i),
                                            nodes[j].r);
                }

                SWOi += (discount1/discount2 - 1.)*discount3 + BasisSpread.at(i)*FloatAccrualBasis.at(i)*discount3;
            }
            SWAP += SWOi*nodes[j].ArrowDebreu;
            SWO += FMAX(CallPut*SWOi,0)*nodes[j].ArrowDebreu;
        }
    }

    // si PriceOnWhatNotice est egal a 0 on calcule le prix de la bermudeenne
    else
    {
        uiSlice = tree.nbSlices - 1;
        iNoticeRestantes = NbNotices - 1;


        while( iNoticeRestantes >= 0 && uiSlice >= 0 )
        {

            SWAP = 0.;
            SWO =  0.;
            SWOB = 0.;

            // Si on est sur une Notice on calcule le prix de l'europenne
            if( tree.slices[uiSlice].isNotice )
            {
                NoticeDate = NoticeDates.at(iNoticeRestantes);
                FirstResetDateAfterNoticeIndex = GetFirstResetDateAfterNotice(iNoticeRestantes);
                nodes = tree.slices[uiSlice].nodes;

                if( fabs(NoticeDate-T.at(uiSlice)) > 1.e-7)
                    throw ("Swaption.PriceInTreeWithAnalytics : on calcule un prix sur une slice qui n'est pas une notice date");


                for( j = 0; j < tree.slices[uiSlice].nbNodes; j++ )
                {
                    SWOi = 0.;
                    SWOBi = 0.;
                    for( i = FirstResetDateAfterNoticeIndex; i < NbCoupons; i++ )
                    {
                        // FixedCoupon
                        if( QModel )
                        {
                            r0 = q.at(uiSlice) * (ccyMarketData.GetInstantaneousRate(T.at(uiSlice)) - K.at(uiSlice) ) + K.at(uiSlice);

                            discount3 = A_FixedPayment.at(iNoticeRestantes,i) * exp( - B_FixedPayment.at(iNoticeRestantes,i) * nodes[j].r / r0 );
                        }
                        else
                        {
                            discount3 = BondPrice(  A_FixedPayment.at(iNoticeRestantes,i),
                                                    MeanReversion,
                                                    NoticeDate,
                                                    FixedPaymentDates.at(i),
                                                    nodes[j].r);
                        }

                        SWOi -= FixedCoupon.at(i)*FixedAccrualBasis.at(i)*discount3;

                        // Float Coupon
                        if( QModel )
                        {
                            r0 = q.at(uiSlice) * (ccyMarketData.GetInstantaneousRate(T.at(uiSlice)) - K.at(uiSlice) ) + K.at(uiSlice);

                            discount3 = A_FloatPayment.at(iNoticeRestantes,i) * exp( - B_FloatPayment.at(iNoticeRestantes,i) * nodes[j].r / r0 );

                            discount1 = A_IndexStart.at(iNoticeRestantes,i) * exp( - B_IndexStart.at(iNoticeRestantes,i) * nodes[j].r / r0);

                            discount2 = A_IndexEnd.at(iNoticeRestantes,i) * exp( - B_IndexEnd.at(iNoticeRestantes,i) * nodes[j].r / r0);
                        }
                        else
                        {
                            discount3 = BondPrice(  A_FloatPayment.at(iNoticeRestantes,i),
                                                    MeanReversion,
                                                    NoticeDate,
                                                    FloatPaymentDates.at(i),
                                                    nodes[j].r);
                            discount1 = BondPrice(  A_IndexStart.at(iNoticeRestantes,i),
                                                    MeanReversion,
                                                    NoticeDate,
                                                    IndexStartDates.at(i),
                                                    nodes[j].r);
                            discount2 = BondPrice(  A_IndexEnd.at(iNoticeRestantes,i),
                                                    MeanReversion,
                                                    NoticeDate,
                                                    IndexEndDates.at(i),
                                                    nodes[j].r);
                        }

                        SWOi += (discount1/discount2 - 1.)*discount3 + BasisSpread.at(i)*FloatAccrualBasis.at(i)*discount3;
                    }

                    SWAP += SWOi*nodes[j].ArrowDebreu;
                    SWOi = FMAX(CallPut*SWOi,0);
                    SWO += SWOi*nodes[j].ArrowDebreu;


                    if( uiSlice < tree.nbSlices - 1 )
                    {
                        SWOBi = tree.slices[uiSlice + 1].nodes[nodes[j].k - tree.slices[uiSlice + 1].nodeMin + 1].SWOi * nodes[j].P[1]
                                + tree.slices[uiSlice + 1].nodes[nodes[j].k - tree.slices[uiSlice + 1].nodeMin ].SWOi * nodes[j].P[2]
                                + tree.slices[uiSlice + 1].nodes[nodes[j].k - tree.slices[uiSlice + 1].nodeMin - 1].SWOi * nodes[j].P[3];
                        SWOBi *= exp(- nodes[j].r * dT.at(uiSlice));
                    }

                    nodes[j].SWOi = FMAX(SWOi, SWOBi);
                    SWOB += nodes[j].SWOi * nodes[j].ArrowDebreu;
                }

                // On decremente le nombre de Notice restant a pricer
                iNoticeRestantes--;

            }

            // si on n'est pas sur une notice on ne fait qu'actualiser les SWOi
            else if( uiSlice < tree.nbSlices - 1)
            {
                nodes = tree.slices[uiSlice].nodes;
                for( j = 0; j < tree.slices[uiSlice].nbNodes; j++ )
                {
                    SWOBi = tree.slices[uiSlice + 1].nodes[nodes[j].k - tree.slices[uiSlice + 1].nodeMin + 1].SWOi * nodes[j].P[1]
                            + tree.slices[uiSlice + 1].nodes[nodes[j].k - tree.slices[uiSlice + 1].nodeMin ].SWOi * nodes[j].P[2]
                            + tree.slices[uiSlice + 1].nodes[nodes[j].k - tree.slices[uiSlice + 1].nodeMin - 1].SWOi * nodes[j].P[3];

                    SWOBi *= exp(- nodes[j].r * dT.at(uiSlice));
                    nodes[j].SWOi = FMAX(SWOi, SWOBi);
                }
            }

            // on decrement l'indice de slice a tous les coups (Notice ou pas)
            uiSlice--;
        }

        // si iNoticeRestantes est encore positif c que l'on n'a pas price les europeenes sur toutes les notices
        if( iNoticeRestantes >= 0 )
        {
            throw ("Swaption.PriceInTreeWithAnalytics : on ne detecte pas toutes les notices en pricant la bermudeenne");
        }
    }



    static DKMaille<double> Price(3);

    Price.at(0) = SWAP * 10000.;
    Price.at(1) = SWOB * 10000.;
    Price.at(2) = SWO * 10000.;

    return Price;
}
// Swaption::PriceInTreeWithAnalytics


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

DKMaille<double>& Swaption::PriceInTree(int PriceOnWhatNotice)
{
    if(PriceOnWhatNotice > NbNotices )
    {
        throw("Swaption.PriceInTreeWithAnalytics : PriceOnWhatNotice est superieur au nombre de notice");
    }

    int i, j, iNoticeRestantes, FirstResetDateAfterNoticeIndex;
    double NoticeDate;

    Node* nodes;

    // on set le parametre isNotice des slices
    // a faire peut-etre au moment de la creation de l'arbre
    int uiNotice=0;
    int uiSlice=0;
    while ( uiNotice < NbNotices  && uiSlice < tree.nbSlices )
    {
        if( fabs(T.at(uiSlice)-NoticeDates.at(uiNotice))<1.e-7  )
        {
            tree.slices[uiSlice].isNotice=true;
            uiNotice++;
        }
        uiSlice++;
    }

    if( uiNotice < NbNotices)
        throw("Swaption.PriceInTree : Tree mal construit, il n'y a pas de slice pour toutes les Notice Dates");

    // On calcule le prix de l'europeenne sur la slice PriceOnWhatNotice
    double SWAP, SWO, SWOi, SWOB, SWOBi, discount1, discount2, discount3;
    SWO=0.;
    SWAP = 0.;
    SWOB = 0.;

    // On calcule le prix de l'europeenne sur la slice PriceOnWhatNotice si PriceOnWhatNotice est > 0
    if( PriceOnWhatNotice > 0)
    {

        // on determine la slice de la notice numero PriceOnWhatNotice
        uiNotice=-1;
        uiSlice=0;
        while( uiNotice + 1 < PriceOnWhatNotice && uiSlice < tree.nbSlices - 1)
        {
            uiSlice++;
            if( tree.slices[uiSlice].isNotice )
            {
                uiNotice++;
            }
        }

        NoticeDate = NoticeDates.at(uiNotice);
        FirstResetDateAfterNoticeIndex = GetFirstResetDateAfterNotice(uiNotice);
        nodes = tree.slices[uiSlice].nodes;

        if( fabs(NoticeDate-T.at(uiSlice)) > 1.e-7)
            throw ("PriceInTree : on calcule un prix sur une slice qui n'est pas une notice date");



        for( j = 0; j < tree.slices[uiSlice].nbNodes; j++ )
        {
            SWOi = 0.;

            for( i = FirstResetDateAfterNoticeIndex; i < NbCoupons; i++ )
            {
                discount3 = tree.CalculateZCinTree(uiSlice,j,FixedPaymentDates.at(i));

                SWOi -= FixedCoupon.at(i)*FixedAccrualBasis.at(i)*discount3;

                discount3 = tree.CalculateZCinTree(uiSlice,j,FloatPaymentDates.at(i));
                discount1 = tree.CalculateZCinTree(uiSlice,j,IndexStartDates.at(i));
                discount2 = tree.CalculateZCinTree(uiSlice,j,IndexEndDates.at(i));

                SWOi += (discount1/discount2 - 1.)*discount3 + BasisSpread.at(i)*FloatAccrualBasis.at(i)*discount3;
            }

            SWAP += SWOi*nodes[j].ArrowDebreu;
            SWO += FMAX(CallPut*SWOi,0)*nodes[j].ArrowDebreu;
        }
    }

    // si PriceOnWhatNotice est egal a 0 on calcule le prix de la bermudeenne
    else
    {
        uiSlice = tree.nbSlices - 1;
        iNoticeRestantes = NbNotices - 1;


        while( iNoticeRestantes >= 0 && uiSlice >= 0 )
        {

            SWAP = 0.;
            SWO =  0.;
            SWOB = 0.;

            // Si on est sur une Notice on calcule le prix de l'europenne
            if( tree.slices[uiSlice].isNotice )
            {
                NoticeDate = NoticeDates.at(iNoticeRestantes);
                FirstResetDateAfterNoticeIndex = GetFirstResetDateAfterNotice(iNoticeRestantes);
                nodes = tree.slices[uiSlice].nodes;

                if( fabs(NoticeDate-T.at(uiSlice)) > 1.e-7)
                    throw ("Swaption.PriceInTreeWithAnalytics : on calcule un prix sur une slice qui n'est pas une notice date");


                for( j = 0; j < tree.slices[uiSlice].nbNodes; j++ )
                {
                    SWOi = 0.;
                    SWOBi = 0.;
                    for( i = FirstResetDateAfterNoticeIndex; i < NbCoupons; i++ )
                    {
                        // FixedCoupon
                        discount3 = tree.CalculateZCinTree(uiSlice,j,FloatPaymentDates.at(i));

                        SWOi -= FixedCoupon.at(i)*FixedAccrualBasis.at(i)*discount3;

                        // Float Coupon
                        discount3 = tree.CalculateZCinTree(uiSlice,j,FloatPaymentDates.at(i));
                        discount1 = tree.CalculateZCinTree(uiSlice,j,IndexStartDates.at(i));
                        discount2 = tree.CalculateZCinTree(uiSlice,j,IndexEndDates.at(i));

                        SWOi += (discount1/discount2 - 1.)*discount3 + BasisSpread.at(i)*FloatAccrualBasis.at(i)*discount3;
                    }

                    SWAP += SWOi*nodes[j].ArrowDebreu;
                    SWOi = FMAX(CallPut*SWOi,0);
                    SWO += SWOi*nodes[j].ArrowDebreu;


                    if( uiSlice < tree.nbSlices - 1 )
                    {
                        SWOBi = tree.slices[uiSlice + 1].nodes[nodes[j].k - tree.slices[uiSlice + 1].nodeMin + 1].SWOi * nodes[j].P[1]
                                + tree.slices[uiSlice + 1].nodes[nodes[j].k - tree.slices[uiSlice + 1].nodeMin ].SWOi * nodes[j].P[2]
                                + tree.slices[uiSlice + 1].nodes[nodes[j].k - tree.slices[uiSlice + 1].nodeMin - 1].SWOi * nodes[j].P[3];
                        SWOBi *= exp(- nodes[j].r * dT.at(uiSlice));
                    }

                    nodes[j].SWOi = FMAX(SWOi, SWOBi);
                    SWOB += nodes[j].SWOi * nodes[j].ArrowDebreu;
                }

                // On decremente le nombre de Notice restant a pricer
                iNoticeRestantes--;

            }

            // si on n'est pas sur une notice on ne fait qu'actualiser les SWOi
            else if( uiSlice < tree.nbSlices - 1)
            {
                nodes = tree.slices[uiSlice].nodes;
                for( j = 0; j < tree.slices[uiSlice].nbNodes; j++ )
                {
                    SWOBi = tree.slices[uiSlice + 1].nodes[nodes[j].k - tree.slices[uiSlice + 1].nodeMin + 1].SWOi * nodes[j].P[1]
                            + tree.slices[uiSlice + 1].nodes[nodes[j].k - tree.slices[uiSlice + 1].nodeMin ].SWOi * nodes[j].P[2]
                            + tree.slices[uiSlice + 1].nodes[nodes[j].k - tree.slices[uiSlice + 1].nodeMin - 1].SWOi * nodes[j].P[3];

                    SWOBi *= exp(- nodes[j].r * dT.at(uiSlice));
                    nodes[j].SWOi = FMAX(SWOi, SWOBi);
                }
            }

            // on decrement l'indice de slice a tous les coups (Notice ou pas)
            uiSlice--;
        }

        // si iNoticeRestantes est encore positif c que l'on n'a pas price les europeenes sur toutes les notices
        if( iNoticeRestantes >= 0 )
        {
            throw ("Swaption.PriceInTreeWithAnalytics : on ne detecte pas toutes les notices en pricant la bermudeenne");
        }
    }

    static DKMaille<double> Price(3);
    Price.at(0) = SWAP * 10000.;
    Price.at(1) = SWOB * 10000.;
    Price.at(2) = SWO * 10000.;
    return Price;
}
// Swaption::PriceInTree


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Swaption::SetMeanReversionForTree()
{
    MeanReversionInTree.resize(ccyMarketData.Sigma.entries());

    for( int i = 0; i < MeanReversionInTree.entries()-1; i++)
    {
        MeanReversionInTree.at(i) = MeanReversion
                                    + (ccyMarketData.Sigma.at(i+1)-ccyMarketData.Sigma.at(i))/(T.at(i+1)-T.at(i))/ccyMarketData.Sigma.at(i);
    }

    MeanReversionInTree.at(MeanReversionInTree.entries()-1)=MeanReversionInTree.at(MeanReversionInTree.entries()-2);
}
// Swaption:SetMeanReversionForTree


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int Swaption::GetFirstResetDateAfterNotice(int uiNotice)
{
    int i = 0.;

    while ( i < FloatResetDates.entries() && FloatResetDates.at(i) < NoticeDates.at(uiNotice) )
    {
        i++;
    }

    if( i == FloatResetDates.entries() )
    {
        throw("Swaption.GetFirstResetDateAfterNotice : Pas de ResetDate apres la NoticeDate");
    }

    return i;
}
// Swaption:GetFirstResetDateAfterNotice


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void GetBootstrappedVolForQModelInTree( DKMaille<double>& SigmaDates,
                                        DKMaille<double>& Sigma,
                                        DKMaille<double> dNoticeDates,
                                        DKMaille<double> ZC_Dates,
                                        DKMaille<double> ZC_Rates,
                                        DKMaille<double> Basis_ZC_Rates,
                                        DKMaille<double> dStdDevX,
                                        DKMaille<double> dStdDevY,
                                        DKMaille2D<double> dStdDevZ,
                                        double dSpotDate,
                                        double dFinalMaturity,
                                        double dNumTimeLinesBeforeFirstNotice,
                                        double dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
                                        double dNumTimeLinesPerYear,
                                        double dOptimal,
                                        double dMeanReversion,
                                        double dSmileParameter)
{
    // Sort Notice Dates for bootstrapping
    DKMaille<double> dNoticeDatesSelected;
    SelectBootstrappingDates(dNoticeDatesSelected,dNoticeDates,dStdDevX);

    dNoticeDates.resize(dNoticeDatesSelected.entries()+2);
    DKMaille<double> dSwapStartDates(dNoticeDatesSelected.entries());
    DKMaille<double> dSwapEndDates(dNoticeDatesSelected.entries());

    dNoticeDates.at(0)=dSpotDate;

    for( unsigned int uiD = 0; uiD < dNoticeDatesSelected.entries(); uiD++ )
    {
        dNoticeDates.at(uiD+1) = dSpotDate+dNoticeDatesSelected.at(uiD) * 365.;

        // dNoticePeriod is set to 0 days for cash settlement in the bootstrapping procedure. It will be questionable to change it
        dSwapStartDates.at(uiD) = dNoticeDates.at(uiD+1) + 0.;

        // Only final maturity diagonal is allowed for calibration purposes for the moment
        // Change the below line to make it accept any date
        dSwapEndDates.at(uiD) = dFinalMaturity *365 + dSpotDate;
    }

    dNoticeDates.at(dNoticeDates.entries()-1) = dFinalMaturity *365 + dSpotDate;;

    DKMaille<double> dModelParameters(2);
    dModelParameters.at(0) = 1.;
    dModelParameters.at(1) = dMeanReversion;

    DKMaille<double> VolStrip;

    VolStrip = Bootstrapping_DK3F_Numerical_SP( dStdDevZ,
               dStdDevX,
               dStdDevY,
               ZC_Dates,
               ZC_Rates,
               Basis_ZC_Rates,
               dNoticeDates,
               dSwapStartDates,
               dSwapEndDates,
               dModelParameters,
               dSpotDate,
               dNumTimeLinesBeforeFirstNotice,
               dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
               dNumTimeLinesPerYear,
               dOptimal,
               dSmileParameter);

    Sigma.resize(VolStrip.entries()+2);
    SigmaDates.resize(VolStrip.entries()+2);

    for( uiD = 0; uiD < SigmaDates.entries(); uiD++ )
    {
        SigmaDates.at(uiD) = (dNoticeDates.at(uiD)-dSpotDate)/365.;
    }

    for( uiD = 0; uiD < VolStrip.entries(); uiD++ )
    {
        Sigma.at(uiD) = VolStrip.at(uiD);
    }
    Sigma.at(Sigma.entries() - 2) = Sigma.at(Sigma.entries()-3);
    Sigma.at(Sigma.entries() - 1) = Sigma.at(Sigma.entries()-3);

}
// GetBootstrappedVolForQModelInTree


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

DKMaille<double> Bootstrapping_DK3F_Numerical_SP(  DKMaille2D<double>& dStdDevZ,
        DKMaille<double>& dStdDevX,
        DKMaille<double>& dStdDevY,
        DKMaille<double>& ZC_Dates,
        DKMaille<double>& ZC_Rates,
        DKMaille<double>& Basis_ZC_Rates,
        DKMaille<double>& dNoticeDates,
        DKMaille<double>& dSwapStartDates,
        DKMaille<double>& dSwapEndDates,
        DKMaille<double>& dModelParameters,
        double dJulianObservationDate,
        double dNumTimeLinesBeforeFirstNotice,
        double dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
        double dNumTimeLinesPerYear,
        double dOptimal,
        double dSmileParameter)
{
    if( dSwapStartDates.entries() != dNoticeDates.entries() - 2 || dSwapEndDates.entries() != dNoticeDates.entries() - 2 )
    {
        throw("GetBootstrappedVolForQModelInTree: Wrong number of entries in dSwapStartDate or dSwapEndDate");
    }

    for( int i = 0; i < dSwapStartDates.entries(); i++ )
    {
        if( dSwapStartDates.at(i)  !=  dNoticeDates.at(i+1) )
        {
            throw("GetBootstrappedVolForQModelInTree: Les dSwapStartDate doivent etre egale au dNoticeDates");
        }

        if( dSwapEndDates.at(i)  !=  dSwapEndDates.at(0) )
        {
            throw("GetBootstrappedVolForQModelInTree: Les dSwapEndDate doivent toutes etre egales a dFinalMaturity");
        }
    }

    if( dModelParameters.at(0) != 1.)
    {
        throw("GetBootstrappedVolForQModelInTree: Le bootstrapping n'a ete developpe que en un facteur");
    }

    if( dJulianObservationDate != dNoticeDates.at(0) )
    {
        throw("GetBootstrappedVolForQModelInTree: La date d'observation doit etre egale a la spotdate");
    }

    double dNoticePeriodInDays = 0.;
    double dSpotDate = dNoticeDates.at(0);
    double dFinalMaturity = ( dSwapEndDates.at(0) - dSpotDate ) /365.;
    double dMeanReversion = dModelParameters.at(1);


    DKMaille<double> dNoticeDatesSelected;
    dNoticeDatesSelected.resize(dNoticeDates.entries()-2);

    for( i = 0; i < dNoticeDatesSelected.entries(); i++ )
    {
        dNoticeDatesSelected.at(i) = ( dNoticeDates.at(i+1) - dNoticeDates.at(0) ) /365.;
    }



    // On genere les dates du swap sous-jacent.
    DKMaille<double> dSwaptionDates;
    DKMaille<double> dFixedAccrualBasisBis;
    DKMaille<double> dFloatResetDates;
    DKMaille<double> dFixedPaymentDates;
    DKMaille<double> dFloatPaymentDates;
    DKMaille<double> dIndexStartDates;
    DKMaille<double> dIndexEndDates;
    DKMaille<double> dFixedAccrualBasis;
    DKMaille<double> dFloatAccrualBasis;

    DKMaille<double> dFixedCoupon;
    DKMaille<double> dBasisSpread;

    CreateVanillaSwap(  dSpotDate,
                        (dFinalMaturity - dNoticeDatesSelected.at(0)),
                        dNoticeDatesSelected.at(0),
                        dSwaptionDates,
                        dFixedAccrualBasisBis,
                        dFloatPaymentDates,
                        dFloatAccrualBasis,
                        dBasisSpread,
                        dNoticePeriodInDays,
                        ZC_Dates,
                        ZC_Rates,
                        Basis_ZC_Rates);


    // On met toutes les dates au format boosterdata
    int nbRows = dFloatPaymentDates.entries();
    dFloatResetDates.resize(nbRows);
    dFixedPaymentDates.resize(nbRows);
    dFloatPaymentDates.resize(nbRows);
    dIndexStartDates.resize(nbRows);
    dIndexEndDates.resize(nbRows);
    dFixedAccrualBasis.resize(nbRows);
    dFloatAccrualBasis.resize(nbRows);
    dFixedCoupon.resize(nbRows);

    for( i = 0; i < nbRows; i++)
    {
        dFloatResetDates.at(i) = (dSwaptionDates.at(i) - dSpotDate)/ 365.;
        dFixedPaymentDates.at(i) = (dSwaptionDates.at(i + 1) - dSpotDate)/ 365.;
        dFloatPaymentDates.at(i) = (dFloatPaymentDates.at(i) - dSpotDate)/ 365.;
        dIndexStartDates.at(i) = (dSwaptionDates.at(i) - dSpotDate)/ 365.;
        dIndexEndDates.at(i) = dFixedPaymentDates.at(i);
        dFixedAccrualBasis.at(i) = dFixedAccrualBasisBis.at(i+1);
    }


    // on genere une swaption qui permettra de calibrer la vol
    Swaption swaption;

    swaption.ccyMarketData.Init(dSpotDate,
                                ZC_Dates,
                                Basis_ZC_Rates,
                                ZC_Dates,
                                ZC_Rates);

    swaption.Init(  dSpotDate,
                    dNoticeDatesSelected,
                    dNumTimeLinesBeforeFirstNotice,
                    dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
                    dNumTimeLinesPerYear,
                    dOptimal,
                    dFinalMaturity,
                    dMeanReversion,
                    dSmileParameter,
                    dFloatResetDates,
                    dFixedPaymentDates,
                    dFloatPaymentDates,
                    dIndexStartDates,
                    dIndexEndDates,
                    dFixedAccrualBasis,
                    dFloatAccrualBasis,
                    dFixedCoupon,
                    dBasisSpread);



    // on determine le vevteurs des vol a fitter
    DKMaille<double> InputVol(dNoticeDatesSelected.entries());
    for( i = 0; i < InputVol.entries(); i++ )
    {
        InputVol.at(i) = InterpolateMatrix( dStdDevZ,
                                            dNoticeDatesSelected.at(i),
                                            dFinalMaturity-dNoticeDatesSelected.at(i),
                                            dStdDevX,
                                            dStdDevY);
    }

    static DKMaille<double> Sigma;
    swaption.BootstrapVolForQModelInTree(   Sigma,
                                            InputVol,
                                            dNoticePeriodInDays);

    return Sigma;

}
// GetBootstrappedVolForQModelInTree


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Swaption::BootstrapVolForQModelInTree(  DKMaille<double>& OutputVol,
                                        DKMaille<double> InputVol,
                                        double dNoticePeriodInDays)
{
    int i, j, iCounter, nbRows;
    double dAdjSwapRate, dAdjAnnuity, dSwapRate, dMarketPrice, Price1, Price2;
    double dGuess, dGuessMin, dGuessMax, dGuessPrime, PriceMin, PriceMax;

    DKMaille<double> dFloatResetDates;
    DKMaille<double> dFixedPaymentDates;
    DKMaille<double> dFloatPaymentDates;
    DKMaille<double> dIndexStartDates;
    DKMaille<double> dIndexEndDates;
    DKMaille<double> dFixedAccrualBasis;
    DKMaille<double> dFloatAccrualBasis;
    DKMaille<double> dFixedCoupon;
    DKMaille<double> dBasisSpread;

    DKMaille<double> swaptionDates;
    DKMaille<double> accrualPeriods;
    DKMaille<double> basisSwaptionDates;
    DKMaille<double> basis;
    DKMaille<double> basisAccrualPeriods;

    // On initialise le vecteur de vol avec la premiere vol a fitter
    DKMaille<double> SigmaDates;
    DKMaille<double> Sigma;
    Sigma.resize(InputVol.entries()+2);
    SigmaDates.resize(InputVol.entries()+2);

    SigmaDates.at(0) = 0.;
    Sigma.at(0) = InputVol.at(0);

    for( i = 1; i <= InputVol.entries(); i++)
    {
        SigmaDates.at(i) = NoticeDates.at(i - 1);
        Sigma.at(i) = InputVol.at(0);
    }

    SigmaDates.at(Sigma.entries()-1) = T.at(NbPasTotal);
    Sigma.at(Sigma.entries()-1) = InputVol.at(0);

    for( i = 0; i < InputVol.entries(); i++)
    {
        dGuess = Sigma.at(i);
        dGuessMin = 0.0000001;
        dGuessMax = 10.;
        PriceMin = 0.;
        PriceMax = 1.e10;

        //On genere les vecteurs necessaires au calcul du MarketPrice
        CreateVanillaSwap(  SpotDate,
                            (T.at(NbPasTotal) - NoticeDates.at(i)),
                            NoticeDates.at(i),
                            swaptionDates,
                            accrualPeriods,
                            basisSwaptionDates,
                            basisAccrualPeriods,
                            basis,
                            dNoticePeriodInDays,
                            ccyMarketData.ZC_Dates,
                            ccyMarketData.ZC_Rates,
                            ccyMarketData.Basis_ZC_Rates);



        // Get AdjSwap Rate
        dAdjSwapRate=GetAdjSwapRate(SpotDate,
                                    swaptionDates,
                                    accrualPeriods,
                                    basisSwaptionDates,
                                    basis,
                                    basisAccrualPeriods,
                                    ccyMarketData.ZC_Dates,
                                    ccyMarketData.ZC_Rates,
                                    ccyMarketData.Basis_ZC_Rates);

        // Get AdjSwap Annuity
        dAdjAnnuity=GetAdjAnnuity(  SpotDate,
                                    swaptionDates,
                                    accrualPeriods,
                                    ccyMarketData.ZC_Dates,
                                    ccyMarketData.ZC_Rates,
                                    ccyMarketData.Basis_ZC_Rates,
                                    SpotDate);

        // Get Swap Rate
        dSwapRate=GetSwapRate(  SpotDate,
                                swaptionDates,
                                accrualPeriods,
                                ccyMarketData.ZC_Dates,
                                ccyMarketData.ZC_Rates,
                                ccyMarketData.Basis_ZC_Rates);


        // Get straddle price from AbsoluteVol
        dMarketPrice=2.*GetPriceFromAbsVol( SpotDate,
                                            SpotDate,
                                            swaptionDates,
                                            accrualPeriods,
                                            dNoticePeriodInDays,
                                            InputVol.at(i),
                                            ccyMarketData.ZC_Dates,
                                            ccyMarketData.ZC_Rates,
                                            ccyMarketData.Basis_ZC_Rates);

        if( i != 0 )
        {
            // On met toutes les dates au format boosterdata
            nbRows = basisSwaptionDates.entries();

            dFloatResetDates.resize(nbRows);
            dFixedPaymentDates.resize(nbRows);
            dFloatPaymentDates.resize(nbRows);
            dIndexStartDates.resize(nbRows);
            dIndexEndDates.resize(nbRows);
            dFixedAccrualBasis.resize(nbRows);
            dFloatAccrualBasis.resize(nbRows);
            dFixedCoupon.resize(nbRows);
            dBasisSpread.resize(nbRows);

            for(int ui = 0; ui < nbRows; ui++)
            {
                dFloatResetDates.at(ui) = (swaptionDates.at(ui) - SpotDate)/ 365.;
                dFixedPaymentDates.at(ui) = (swaptionDates.at(ui + 1) - SpotDate)/ 365.;
                dFloatPaymentDates.at(ui) = (basisSwaptionDates.at(ui) - SpotDate)/ 365.;
                dIndexStartDates.at(ui) = (swaptionDates.at(ui) - SpotDate)/ 365.;
                dIndexEndDates.at(ui) = dFixedPaymentDates.at(ui);
                dFixedAccrualBasis.at(ui) = accrualPeriods.at(ui+1);
                dFloatAccrualBasis.at(ui) = basisAccrualPeriods.at(ui);
                dBasisSpread.at(ui) = basis.at(ui);
            }


            ModifyBoosterData(  dFloatResetDates,
                                dFixedPaymentDates,
                                dFloatPaymentDates,
                                dIndexStartDates,
                                dIndexEndDates,
                                dFixedAccrualBasis,
                                dFloatAccrualBasis,
                                dFixedCoupon,
                                dBasisSpread);
        }

        // On set le taux de la swaption
        for( j = 0; j < FixedCoupon.entries(); j++)
        {
            FixedCoupon.at(j) = dSwapRate;
        }

        CleanTree();

        SetVol(SigmaDates,Sigma);

        CreateTree();

        Price1 = StraddlePrice(i+1).at(2);

        iCounter = 0;

        while( fabs(Price1 - dMarketPrice) > 0.01 && iCounter < 100 && fabs(dGuessMax-dGuessMin) > 1.e-4 )
        {
            iCounter++;

            if( Price1 > dMarketPrice )
            {
                dGuessMax = dGuess;
                PriceMax = Price1;
            }
            else
            {
                dGuessMin = dGuess;
                PriceMin = Price1;
            }

            dGuessPrime = dGuess * 1.001;

            for( j = i; j < Sigma.entries() ;j++ )
            {
                Sigma.at(j) = dGuessPrime;

            }

            CleanTree();

            SetVol(SigmaDates,Sigma);

            CreateTree();

            Price2 = StraddlePrice(i+1).at(2);

            dGuessPrime = (1. + (dMarketPrice - Price1)/(Price2 - Price1) * 0.001) * dGuess;

            if( dGuessPrime < dGuessMin )
            {
                if(dGuess == dGuessMin)
                {
                    dGuessPrime = FMIN( ( dGuess + dGuessMax ) /2. , dGuess * 2.);
                }
                else
                {
                    dGuessPrime = (dGuessMin + dGuess) / 2.;
                }
            }
            if(dGuessPrime > dGuessMax)
            {
                dGuessPrime = (dGuessMax + dGuess) / 2.;
            }
            dGuess = dGuessPrime;

            for( j = i; j < Sigma.entries() ;j++ )
            {
                Sigma.at(j) = dGuess;
            }

            CleanTree();

            SetVol(SigmaDates,Sigma);

            CreateTree();

            Price1 = StraddlePrice(i+1).at(2);
        }
        if(fabs(dGuessMin - dGuessMax) > 0.0001 && fabs(Price1 - dMarketPrice) > 0.01 )
        {
            throw("Swaption.BootstrapVolForQModelInTree : Echec dans la calibration du Sigma ");
        }
    }

    OutputVol.resize(Sigma.entries()-2);
    for( i = 0; i < OutputVol.entries(); i++ )
    {
        OutputVol.at(i) = Sigma.at(i);
    }
}
// Swaption::BootstrapVolForQModelInTree

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Swaption::Print(FILE* file)
{

    T.print(file);
    dT.print(file);

    MeanReversionInTree.print(file);

    K.print(file);
    NoticeDates.print(file);
    FloatResetDates.print(file);
    FixedPaymentDates.print(file);
    FloatPaymentDates.print(file);
    IndexStartDates.print(file);
    IndexEndDates.print(file);
    FixedAccrualBasis.print(file);
    FloatAccrualBasis.print(file);
    FixedCoupon.print(file);
    BasisSpread.print(file);
    A_FixedPayment.print(file);
    A_FloatPayment.print(file);
    A_IndexStart.print(file);
    A_IndexEnd.print(file);
    B_FixedPayment.print(file);
    B_FloatPayment.print(file);
    B_IndexStart.print(file);
    B_IndexEnd.print(file);

    ccyMarketData.Print(fic);
}
// Swaption::Print


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Swaption::SetK()
{
    K.resize(T.entries());

    for( int i = 0; i < K.entries(); i++ )
    {
        if( QModel )
        {
            K.at(i) = FMAX( ccyMarketData.ForwardRates.at(i), 0.002);
        }
        else
        {
            K.at(i) = 0.;
        }
    }
}
// Swaption::SetK


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Swaption::Setq(double SmileParameter)
{
    q.resize(T.entries());

    for( int i = 0; i < q.entries(); i++ )
    {
        if( !QModel || ccyMarketData.ForwardRates.at(i) < 0.002 )
        {
            q.at(i) = 0.;
        }
        else
        {
            q.at(i) = SmileParameter;
        }
    }
}
// Swaption::Setq
