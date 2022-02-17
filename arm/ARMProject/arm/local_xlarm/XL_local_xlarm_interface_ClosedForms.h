       
        {
        		" Local_VanillaOption_Normal",		/// name of the C++ function
           
                " RRRRRR",					    /// 6 parametres = 5 d'entree + 1 parametre de retour 
				
				" ARM_CF_Norm_VanillaOption",
				" undelying, volatility, strike, maturity, callOrPut",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the close formula of E(S-K)+ or E(S-K)- where S is Normal",
				" underling value (real value: 0.05 for 5%)",
				" normal volatility value(real value: 0.05 for 5%)",
				" strike value (real value: 0.01 for 1%)",
				" maturity (real value: 0.1 for 10%)",
				" callOrPut C(ALL) or P(UT)",
        },
		{
        		" Local_VanillaOption_Normal_Der",		/// name of the C++ function
               
                " RRRRRRR",					    /// 7 parametres = 6 d'entree + 1 parametre de retour 				
				" ARM_CF_Norm_VanillaOption_Der",				
				" Index,undelying, volatility, strike, maturity, callOrPut",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
              
                " gives the close formula of first derivative of E(S-K)+ or E(S-K)- where S is Normal",
				" index of the dimension: Delta(0=fwd), Vega(1=vol), Digital=(2=K), Theta(3=t)",
				" underling value (real value: 0.05 for 5%)",
				" normal volatility value(real value: 0.05 for 5%)",
				" strike value (real value: 0.01 for 1%)",
				" maturity (real value: 0.1 for 10%)",
				" callOrPut C(ALL) or P(UT)",
        },
		{
        		" Local_VanillaOption_Normal_Der2",		/// name of the C++ function
                " RRRRRRRR",					    /// 8 parametres = 7 d'entree + 1 parametre de retour 
				" ARM_CF_Norm_VanillaOption_Der2",
				" Index1,Index2,undelying, volatility, strike, maturity, callOrPut",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the merge between close formulas and numerical formulas of second derivatives of E(S-K)+ or E(S-K)- where S is Normal",
				" index of the dimension 1: 0=fwd, 1=vol, 2=k, 3=t",
				" index of the dimension 2: 0=fwd, 1=vol, 2=k, 3=t",
				" underling value (real value: 0.05 for 5%)",
				" normal volatility value(real value: 0.05 for 5%)",
				" strike value (real value: 0.01 for 1%)",
				" maturity (real value: 0.1 for 10%)",
				" callOrPut C(ALL) or P(UT)",
        },
		{
        		" Local_DoubleDigital_Normal",		/// name of the C++ function
                " RRRRRRRRRRRRRRR",					    /// 8 parametres = 7 d'entree + 1 parametre de retour 
				" ARM_CF_Norm_DoubleDigital",
				" Fwd1,Fwd2,Maturity,Strike1,Spread1,Strike2,Spread2,Vol1Plus,Vol1Minus,Vol2Plus,Vol2Minus,Correl,CallOrPut1,CallOrPut2",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " Normal Formula for double digital with call spread replication",
				" First forward",
				" Second forward",
				" Maturity",
				" First Strike",
				" First Strike Spread",
				" Second Strike",
				" Second Strike Spread",
				" Normal Vol for First Strike + Spread",
				" Normal Vol for First Strike - Spread",
				" Normal Vol for Second Strike + Spread",
				" Normal Vol for Second Strike - Spread",
				" Correlation",
				" Call or Put for First Digital",
				" Call or Put for Second Digital",
        },
        {
        		" Local_LogNormal_SpreadOption",		/// name of the C++ function
                " RRRRRRRRRRR",					/// 11 parametres = 10 d'entree + 1 parametre de retour 
				" ARM_CF_LN_SpreadOpt",
				" index1,index2,sigma1,sigma2,correlation,strike,maturity,callorput,optiontype,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the FORWARD VALUE of a spread option (2Lognormals)",
				" fwd value (real value: 0.05 for 5%)",
				" fwd value (real value: 0.05 for 5%)",
				" LN vol (real value: 0.1 for 10%)",
				" LN vol (real value: 0.1 for 10%)",
				" correlation (real value: 0.1 for 10%)",
				" real value: 0.05 for 5%",
				" year fraction (1 for one year)",
				" C(ALL) or P(UT)",
				" DIGITALOPTION,SPREADOPTION,PAYFIRST,PAYSECOND",
				" number of points in the numerical integration (default=64)"
        },
		{
        		" Local_LogNormal_SpreadOption_Calibrate_Correlation",		/// name of the C++ function
                " RRRRRRRRRRR",					/// 11 parametres = 10 d'entree + 1 parametre de retour 
				" ARM_CF_LN_SpreadOptCorr",
				" index1,index2,sigma1,sigma2,optprice,strike,maturity,callorput,optiontype,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the implicit correlation of a spread option price (2Lognormals)",
				" fwd value (real value: 0.05 for 5%)",
				" fwd value (real value: 0.05 for 5%)",
				" LN vol (real value: 0.1 for 10%)",
				" LN vol (real value: 0.1 for 10%)",
				" option price",
				" real value: 0.05 for 5%",
				" year fraction (1 for one year)",
				" C(ALL) or P(UT)",
				" DIGITALOPTION,SPREADOPTION,PAYFIRST,PAYSECOND",
					" number of points in the numerical integration (default=64)"
        },
		{
        		" Local_LogNormal_SpreadOption_Der",		/// name of the C++ function
                " RRRRRRRRRRRR",					/// 12 parametres = 11 d'entree + 1 parametre de retour 
                " ARM_CF_LN_SpreadOptDer",
				" derivative index, index1 value,index2 value,sigma1,sigma2,correlation,strike,maturity,callorput,optiontype,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives a first derivative of the FORWARD VALUE of a spread option (2Lognormals)",
				" index of the dimension: 0=fwd1, 1=fwd2, 2=sig1, 3=sig2, 4=rho, 5=k, 6=t",
				" fwd value (real value: 0.05 for 5%)",
				" fwd value (real value: 0.05 for 5%)",
				" LN vol (real value: 0.1 for 10%)",
				" LN vol (real value: 0.1 for 10%)",
				" correlation (real value: 0.1 for 10%)",
				" real value: 0.05 for 5%",
				" year fraction (1 for one year)",
				" C(ALL) or P(UT)",
				" DIGITALOPTION,SPREADOPTION,PAYFIRST,PAYSECOND",
				" number of points in the numerical integration (default=64)"
        },
		{
        		" Local_LogNormal_SpreadOption_Der2",		/// name of the C++ function
                " RRRRRRRRRRRRR",					/// 13 parametres = 12 d'entree + 1 parametre de retour 
                " ARM_CF_LN_SpreadOptDer2",
				" derivative index1,derivative index2, index1 value,index2 value,sigma1,sigma2,correlation,strike,maturity,callorput,optiontype,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives a second derivative of the FORWARD VALUE of a spread option (2Lognormals)",
				" index of the dimension: 0=fwd1, 1=fwd2, 2=sig1, 3=sig2, 4=rho, 5=k, 6=t",
				" index of the dimension: 0=fwd1, 1=fwd2, 2=sig1, 3=sig2, 4=rho, 5=k, 6=t",
				" fwd value (real value: 0.05 for 5%)",
				" fwd value (real value: 0.05 for 5%)",
				" LN vol (real value: 0.1 for 10%)",
				" LN vol (real value: 0.1 for 10%)",
				" correlation (real value: 0.1 for 10%)",
				" real value: 0.05 for 5%",
				" year fraction (1 for one year)",
				" C(ALL) or P(UT)",
				" DIGITALOPTION,SPREADOPTION,PAYFIRST,PAYSECOND",
				" number of points in the numerical integration (default=64)"
        },
		{
        		" Local_Smiled_LogNormal_SpreadOption",		/// name of the C++ function
                " RRRRRRRRRRRRR",					/// 11 parametres = 10 d'entree + 1 parametre de retour 
				" ARM_CF_LNSmiled_SpreadOpt",
				" index1,index2,sigma1,sigma2,correlation,strike,maturity,slope1,slope2,callorput,optiontype,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the FORWARD VALUE of a spread option with smile slope (2Lognormals)",
				" fwd value (real value: 0.05 for 5%)",
				" fwd value (real value: 0.05 for 5%)",
				" LN vol (real value: 0.1 for 10%)",
				" LN vol (real value: 0.1 for 10%)",
				" correlation (real value: 0.1 for 10%)",
				" real value: 0.05 for 5%",
				" year fraction (1 for one year)",
				" smile slope for the index 1",
				" smile slope for the index 2",
				" C(ALL) or P(UT)",
				" DIGITALOPTION,SPREADOPTION,PAYFIRST,PAYSECOND",
				" number of points in the numerical integration (default=64)"
        },
		{
        		" Local_Smiled_LogNormal_SpreadOption_Calibrate_Correlation",		/// name of the C++ function
                " RRRRRRRRRRRRR",					/// 13 parametres = 12 d'entree + 1 parametre de retour 
				" ARM_CF_LNSmiled_SpreadOptCorr",
				" index1,index2,sigma1,sigma2,optprice,strike,maturity,slope1,slope2,callorput,optiontype,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the implicit correlation of a spread option price  with smile slope(2Lognormals)",
				" fwd value (real value: 0.05 for 5%)",
				" fwd value (real value: 0.05 for 5%)",
				" LN vol (real value: 0.1 for 10%)",
				" LN vol (real value: 0.1 for 10%)",
				" option price",
				" real value: 0.05 for 5%",
				" year fraction (1 for one year)",
				" smile slope for the index 1",
				" smile slope for the index 2",
				" C(ALL) or P(UT)",
				" DIGITALOPTION,SPREADOPTION,PAYFIRST,PAYSECOND",
					" number of points in the numerical integration (default=64)"
        },
		{
        		" Local_Smiled_LogNormal_SpreadOption_Der",		/// name of the C++ function
                " RRRRRRRRRRRRRR",					/// 14 parametres = 13 d'entree + 1 parametre de retour 
                " ARM_CF_LNSmiled_SpreadOptDer",
				" derivative index, index1 value,index2 value,sigma1,sigma2,correlation,strike,maturity,slope1,slope2,callorput,optiontype,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives a first derivative of the FORWARD VALUE of a spread option  with smile slope (2Lognormals)",
				" index of the dimension: 0=fwd1, 1=fwd2, 2=sig1, 3=sig2, 4=rho, 5=k, 6=t",
				" fwd value (real value: 0.05 for 5%)",
				" fwd value (real value: 0.05 for 5%)",
				" LN vol (real value: 0.1 for 10%)",
				" LN vol (real value: 0.1 for 10%)",
				" correlation (real value: 0.1 for 10%)",
				" real value: 0.05 for 5%",
				" year fraction (1 for one year)",
				" smile slope for the index 1",
				" smile slope for the index 2",
				" C(ALL) or P(UT)",
				" DIGITALOPTION,SPREADOPTION,PAYFIRST,PAYSECOND",
				" number of points in the numerical integration (default=64)"
        },
		{
        		" Local_Smiled_LogNormal_SpreadOption_Der2",		/// name of the C++ function
                " RRRRRRRRRRRRRRR",					/// 15 parametres = 14 d'entree + 1 parametre de retour 
                " ARM_CF_LNSmiled_SpreadOptDer2",
				" derivative index1,derivative index2, index1 value,index2 value,sigma1,sigma2,correlation,strike,maturity,slope1,slope2,callorput,optiontype,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives a second derivative of the FORWARD VALUE of a spread option  with smile slope(2Lognormals)",
				" index of the dimension: 0=fwd1, 1=fwd2, 2=sig1, 3=sig2, 4=rho, 5=k, 6=t",
				" index of the dimension: 0=fwd1, 1=fwd2, 2=sig1, 3=sig2, 4=rho, 5=k, 6=t",
				" fwd value (real value: 0.05 for 5%)",
				" fwd value (real value: 0.05 for 5%)",
				" LN vol (real value: 0.1 for 10%)",
				" LN vol (real value: 0.1 for 10%)",
				" correlation (real value: 0.1 for 10%)",
				" real value: 0.05 for 5%",
				" year fraction (1 for one year)",
				" smile slope for the index 1",
				" smile slope for the index 2",
				" C(ALL) or P(UT)",
				" DIGITALOPTION,SPREADOPTION,PAYFIRST,PAYSECOND",
				" number of points in the numerical integration (default=64)"
        },
		{
        		" Local_Normal_SpreadOption",		/// name of the C++ function
                " RRRRRRRRRR",					/// 10 parametres = 9 d'entree + 1 parametre de retour 
				" ARM_CF_Norm_SpreadOpt",
				" index1,index2,sigma1,sigma2,correlation,strike,maturity,callorput,optiontype",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the FORWARD VALUE of a spread option (2normals)",
				" fwd value (real value: 0.05 for 5%)",
				" fwd value (real value: 0.05 for 5%)",
				" LN vol (real value: 0.1 for 10%)",
				" LN vol (real value: 0.1 for 10%)",
				" correlation (real value: 0.1 for 10%)",
				" real value: 0.05 for 5%",
				" year fraction (1 for one year)",
				" C(ALL) or P(UT)",
				" DIGITALOPTION,SPREADOPTION,PAYFIRST,PAYSECOND"
		},
		{
        		" Local_Normal_SpreadOption_Der",		/// name of the C++ function
                " RRRRRRRRRRR",					/// 11 parametres = 10 d'entree + 1 parametre de retour 
                " ARM_CF_Norm_SpreadOptDer",
				" derivative index, index1 value,index2 value,sigma1,sigma2,correlation,strike,maturity,callorput,optiontype",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives a first derivative of the FORWARD VALUE of a spread option (2normals)",
				" index of the dimension: 0=fwd1, 1=fwd2, 2=sig1, 3=sig2, 4=rho, 5=k, 6=t",
				" fwd value (real value: 0.05 for 5%)",
				" fwd value (real value: 0.05 for 5%)",
				" LN vol (real value: 0.1 for 10%)",
				" LN vol (real value: 0.1 for 10%)",
				" correlation (real value: 0.1 for 10%)",
				" real value: 0.05 for 5%",
				" year fraction (1 for one year)",
				" C(ALL) or P(UT)",
				" DIGITALOPTION,SPREADOPTION,PAYFIRST,PAYSECOND"
        },
		{
        		" Local_Normal_SpreadOption_Der2",		/// name of the C++ function
                " RRRRRRRRRRRR",					/// 12 parametres = 12 d'entree + 1 parametre de retour 
                " ARM_CF_Norm_SpreadOptDer2",
				" derivative index1,derivative index2, index1 value,index2 value,sigma1,sigma2,correlation,strike,maturity,callorput,optiontype",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives a second derivative of the FORWARD VALUE of a spread option (2normals)",
				" index of the dimension: 0=fwd1, 1=fwd2, 2=sig1, 3=sig2, 4=rho, 5=k, 6=t",
				" index of the dimension: 0=fwd1, 1=fwd2, 2=sig1, 3=sig2, 4=rho, 5=k, 6=t",
				" fwd value (real value: 0.05 for 5%)",
				" fwd value (real value: 0.05 for 5%)",
				" LN vol (real value: 0.1 for 10%)",
				" LN vol (real value: 0.1 for 10%)",
				" correlation (real value: 0.1 for 10%)",
				" real value: 0.05 for 5%",
				" year fraction (1 for one year)",
				" C(ALL) or P(UT)",
				" DIGITALOPTION,SPREADOPTION,PAYFIRST,PAYSECOND"
        },
		{
        		" Local_Gaussian_SABR_Power_SpreadOption",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRRR",					/// 20 parametres = 19 d'entree + 1 parametre de retour 
                " ARM_CF_SABRGauss_SpreadOpt",
				" index1,index2 ,alpha1,beta1,rho1,nu1,alpha2,beta2,rho2,nu2,correlation,maturity,flag,a1,b1,k1,a2,b2,k2,params",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the  FORWARD VALUE of a power spread option (SABR+Gaussian)",
				" fwd value (real value: 0.05 for 5%)",
				" fwd value (real value: 0.05 for 5%)",
				" alpha index1",
				" beta index1",
				" rho index1",
				" nu index1",
				" alpha index2",
				" beta index2",
				" rho index2",
				" nu index2",
				" copula correlation: 0.5 for 50%",
				" year fraction (1 for one year)",
				" Type of SABR formulas : SABR_IMPLNVOL(2),SABR_A,SABR_G,DIRECTEXACT,DIRECTGEOMETRIC,DIRECTARITHMETIC,NORMALEXACT,NORMALGEOMETRIC,NORMALARITHMETIC,ANALYTICZP0,ANALYTICZP2",
				" a1: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" b1: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" k1: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" a2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" b2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" k2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" Nb(=120),alpha_exp(=0.1),alpha_tanh(=1.5),kb_tanh=(=0.02)"
        },
		{
        		" Local_Gaussian_SABR_Power_Digital_SpreadOption",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRR",					/// 18 parametres = 17 d'entree + 1 parametre de retour 
                " ARM_CF_SABRGauss_DigSpreadOpt",
				" index1,index2 ,alpha1,beta1,rho1,nu1,alpha2,beta2,rho2,nu2,correlation,maturity,flag,k1,a2,b2,k2,params",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the  FORWARD VALUE of a power spread option (SABR+Gaussian)",
				" fwd value (real value: 0.05 for 5%)",
				" fwd value (real value: 0.05 for 5%)",
				" alpha index1",
				" beta index1",
				" rho index1",
				" nu index1",
				" alpha index2",
				" beta index2",
				" rho index2",
				" nu index2",
				" copula correlation: 0.5 for 50%",
				" year fraction (1 for one year)",
				" Type of SABR formulas : SABR_IMPLNVOL(2),SABR_A,SABR_G,DIRECTEXACT,DIRECTGEOMETRIC,DIRECTARITHMETIC,NORMALEXACT,NORMALGEOMETRIC,NORMALARITHMETIC,ANALYTICZP0,ANALYTICZP2",
				" k1: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" a2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" b2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" k2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" Nb(=120),alpha_exp(=0.1),alpha_tanh(=1.5),kb_tanh=(=0.02)"
        },
		{
        		" Local_Gaussian_SABR_Power_SpreadOption_Der",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRRRR",					/// 21 parametres = 20 d'entree + 1 parametre de retour 
                " ARM_CF_SABRGauss_SpreadOptDer",
				" derivative rank, index1,index2 ,alpha1,beta1,rho1,nu1,alpha2,beta2,rho2,nu2,correlation,maturity,flag,a1,b1,k1,a2,b2,k2,params",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives a first derivative of the  FORWARD VALUE of a power spread option (SABR+Gaussian)",
				" index of the dimension: 0=fwd1, 1=fwd2, 2=alpha1, 3=beta1, 4=rho1, 5=nu1, 6=alpha2, 7=beta2, 8=rho2, 9=nu2, 10= corr, 11=t",
				" fwd value (real value: 0.05 for 5%)",
				" fwd value (real value: 0.05 for 5%)",
				" alpha index1",
				" beta index1",
				" rho index1",
				" nu index1",
				" alpha index2",
				" beta index2",
				" rho index2",
				" nu index2",
				" copula correlation: 0.5 for 50%",
				" year fraction (1 for one year)",
				" Type of SABR formulas : SABR_IMPLNVOL(2),SABR_A,SABR_G,DIRECTEXACT,DIRECTGEOMETRIC,DIRECTARITHMETIC,NORMALEXACT,NORMALGEOMETRIC,NORMALARITHMETIC,ANALYTICZP0,ANALYTICZP2",
				" a1: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" b1: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" k1: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" a2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" b2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" k2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
		//		" number of points in the numerical integration (default=64)"
			
        },
		{
        		" Local_Gaussian_SABR_Power_Digital_SpreadOption_Der",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRR",					/// 19 parametres = 18 d'entree + 1 parametre de retour 	
				" ARM_CF_SABRGauss_DigSpreadOptDer",
				" derivative rank, index1,index2 ,alpha1,beta1,rho1,nu1,alpha2,beta2,rho2,nu2,correlation,maturity,flag,k1,a2,b2,k2,params",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives a first derivative of the  FORWARD VALUE of a power spread option (SABR+Gaussian)",
				" index of the dimension: 0=fwd1, 1=fwd2, 2=alpha1, 3=beta1, 4=rho1, 5=nu1, 6=alpha2, 7=beta2, 8=rho2, 9=nu2, 10= corr, 11=t",
				" fwd value (real value: 0.05 for 5%)",
				" fwd value (real value: 0.05 for 5%)",
				" alpha index1",
				" beta index1",
				" rho index1",
				" nu index1",
				" alpha index2",
				" beta index2",
				" rho index2",
				" nu index2",
				" copula correlation: 0.5 for 50%",
				" year fraction (1 for one year)",
				" Type of SABR formulas : SABR_IMPLNVOL(2),SABR_A,SABR_G,DIRECTEXACT,DIRECTGEOMETRIC,DIRECTARITHMETIC,NORMALEXACT,NORMALGEOMETRIC,NORMALARITHMETIC,ANALYTICZP0,ANALYTICZP2",
				" k1: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" a2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" b2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" k2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
		//		" number of points in the numerical integration (default=64)"
			
        },
		{
        		" Local_Gaussian_SABR_Power_SpreadOption_Der2",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRRRRR" ,					/// 22 parametres = 21 d'entree + 1 parametre de retour 
                " ARM_CF_SABRGauss_SpreadOptDer2",
				" derivative rank 1,derivative rank 2 , index1,index2 ,alpha1,beta1,rho1,nu1,alpha2,beta2,rho2,nu2,correlation,maturity,flag,a1,b1,k1,a2,b2,k2,params",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives a second derivative of the  FORWARD VALUE of a power spread option (SABR+Gaussian)",
				" index of the dimension: 0=fwd1, 1=fwd2, 2=alpha1, 3=beta1, 4=rho1, 5=nu1, 6=alpha2, 7=beta2, 8=rho2, 9=nu2, 10= corr, 11=t",
				" index of the dimension: 0=fwd1, 1=fwd2, 2=alpha1, 3=beta1, 4=rho1, 5=nu1, 6=alpha2, 7=beta2, 8=rho2, 9=nu2, 10= corr, 11=t",
				" fwd value (real value: 0.05 for 5%)",
				" fwd value (real value: 0.05 for 5%)",
				" alpha index1",
				" beta index1",
				" rho index1",
				" nu index1",
				" alpha index2",
				" beta index2",
				" rho index2",
				" nu index2",
				" copula correlation: 0.5 for 50%",
				" year fraction (1 for one year)",
				" Type of SABR formulas : SABR_IMPLNVOL(2),SABR_A,SABR_G,DIRECTEXACT,DIRECTGEOMETRIC,DIRECTARITHMETIC,NORMALEXACT,NORMALGEOMETRIC,NORMALARITHMETIC,ANALYTICZP0,ANALYTICZP2",
				" a1: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" b1: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" k1: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" a2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" b2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
//				" k2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
//				" number of points in the numerical integration (default=64)"
        },
			{
        		" Local_Gaussian_SABR_Power_Digital_SpreadOption_Der2",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRRR" ,					/// 22 parametres = 21 d'entree + 1 parametre de retour 		
				" ARM_CF_SABRGauss_DigSpreadOptDer2",
				" derivative rank 1,derivative rank 2 , index1,index2 ,alpha1,beta1,rho1,nu1,alpha2,beta2,rho2,nu2,correlation,maturity,flag,k1,a2,b2,k2,params",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives a second derivative of the  FORWARD VALUE of a power spread option (SABR+Gaussian)",
				" index of the dimension: 0=fwd1, 1=fwd2, 2=alpha1, 3=beta1, 4=rho1, 5=nu1, 6=alpha2, 7=beta2, 8=rho2, 9=nu2, 10= corr, 11=t",
				" index of the dimension: 0=fwd1, 1=fwd2, 2=alpha1, 3=beta1, 4=rho1, 5=nu1, 6=alpha2, 7=beta2, 8=rho2, 9=nu2, 10= corr, 11=t",
				" fwd value (real value: 0.05 for 5%)",
				" fwd value (real value: 0.05 for 5%)",
				" alpha index1",
				" beta index1",
				" rho index1",
				" nu index1",
				" alpha index2",
				" beta index2",
				" rho index2",
				" nu index2",
				" copula correlation: 0.5 for 50%",
				" year fraction (1 for one year)",
				" Type of SABR formulas : SABR_IMPLNVOL(2),SABR_A,SABR_G,DIRECTEXACT,DIRECTGEOMETRIC,DIRECTARITHMETIC,NORMALEXACT,NORMALGEOMETRIC,NORMALARITHMETIC,ANALYTICZP0,ANALYTICZP2",
				" k1: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" a2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" b2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" k2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" Nb(=120),alpha_exp(=0.1),alpha_tanh(=1.5),kb_tanh=(=0.02)"
        },
		{
        		" Local_Gaussian_SABR_Power_SpreadOption_Certitude",		/// name of the C++ function
                " RRRR",					/// 4 parametres = 3 d'entree + 1 parametre de retour 
                " ARM_CF_SABR_SpreadOptCertitude",
				" correlation,maturity,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the probability of the domain taken into accout in the integration",
				" copula correlation: 0.5 for 50%",
				" year fraction (1 for one year)",
				" number of points in the numerical integration (default=64)"
        },
		{
        		" Local_CF_LN_VanillaOption",		/// name of the C++ function
                " RRRRRRR",					/// 7 parametres = 6 d'entree + 1 parametre de retour 
                " ARM_CF_LN_VanillaOption",
				" forward,volatility,bondprice,strike,maturity,CallPut",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the present value of the black and Scholes option",
				" forward value: 0.5 for 50%",
				" total volatility of the forward value: 0.5 for 50% (equal to v*sqrt[t])",
				" discount factor : 0.85 for 85%",
				" strike: 0.5 for 50%",
				" maturity: in years",
				" (C)all=1 or (P)ut=-1 "
        },
		{
        		" Local_CF_BlackSholes",		/// name of the C++ function
                " RRRRRR",					/// 6 parametres = 5 d'entree + 1 parametre de retour 
                " ARM_CF_LN_Option",
				" forward,totalvolatility,bondprice,strike,CallPut",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the present value of the black and Scholes option",
				" forward value: 0.5 for 50%",
				" total volatility of the forward value: 0.5 for 50% (equal to v*sqrt[t])",
				" discount factor : 0.85 for 85%",
				" strike: 0.5 for 50%",
				" (C)all=1 or (P)ut=-1 "
        },
		{
        		" Local_CF_BlackSholes_Der",		/// name of the C++ function
                " RRRRRRR",					/// 7 parametres = 6 d'entree + 1 parametre de retour 
                " ARM_CF_LN_OptionDer",
				" index,forward,totalvolatility,bondprice,strike,CallPut",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the present value of the black and Scholes option",
				" index to derive : 0=forward, 1=totalvolatility, 2=discount, 3=strike",
				" forward value: 0.5 for 50%",
				" total volatility of the forward value: 0.5 for 50% (equal to v*sqrt[t])",
				" discount factor : 0.85 for 85%",
				" strike: 0.5 for 50%",
				" (C)all=1 or (P)ut=-1 "
        },
		{
        		" Local_CF_BlackSholes_Der2",		/// name of the C++ function
                " RRRRRRRR",					/// 8 parametres = 7 d'entree + 1 parametre de retour 
                " ARM_CF_LN_OptionDer2",
				" index1,index2,forward,totalvolatility,bondprice,strike,CallPut",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the present value of the black and Scholes option",
				" index to derive : 0=forward, 1=totalvolatility, 2=discount, 3=strike",
				" index to derive : 0=forward, 1=totalvolatility, 2=discount, 3=strike",
				" forward value: 0.5 for 50%",
				" total volatility of the forward value: 0.5 for 50% (equal to v*sqrt[t])",
				" discount factor : 0.85 for 85%",
				" strike: 0.5 for 50%",
				" (C)all=1 or (P)ut=-1 "
        },
		{
        		" Local_CF_BlackSholes_ImplicitVolatility",		/// name of the C++ function
                " RRRRRRRR",					/// 8 parametres = 7 d'entree + 1 parametre de retour 
                " ARM_CF_LN_ImplicitVol",
				" forward,bondprice,strike,maturity,CallPut,optprice,algo",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the present value of the black and Scholes option",
				" forward value: 0.5 for 50%",
				" discount factor : 0.85 for 85%",
				" strike: 0.5 for 50%",
				" maturity (in years) ",
				" (C)all=1 or (P)ut=-1 ",
				" option present price ",
				" Algorithm (default=1)"
        },
		{
        		" Local_CF_BlackSholes_ImplicitVol",		/// name of the C++ function
                " RRRRRRR",					/// 7 parametres = 6 d'entree + 1 parametre de retour 
                " ARM_CF_LN_ImpVol",
				" forward,bondprice,strike,CallPut,optprice,algo",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives total volatility value of the black and Scholes option (sqrt(T)*sigma)",
				" forward value: 0.5 for 50%",
				" discount factor : 0.85 for 85%",
				" strike: 0.5 for 50%",
				" (C)all=1 or (P)ut=-1 ",
				" option present price ",
				" Algorithm (default=1"
        },
		{
        		" Local_Gaussian_ShiftedLN_Power_SpreadOption",		/// name of the C++ function
                " RRRRRRRRRRRRRRRR",					/// 16 parametres = 15 d'entree + 1 parametre de retour 
                " ARM_CF_LNShifted_SpreadOpt",
				" index1,index2,sigma1,alpha1,sigma2,alpha2,correlation,maturity,a1,b1,k1,a2,b2,k2,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the  FORWARD VALUE of a power spread option (LNShifted+Gaussian)",
				" fwd value (real value: 0.05 for 5%)",
				" fwd value (real value: 0.05 for 5%)",
				" sigma index1",
				" alpha index1",
				" sigma index2",
				" alpha index2",
				" copula correlation: 0.5 for 50%",
				" year fraction (1 for one year)",
				" a1: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" b1: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" k1: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" a2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" b2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" k2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" number of points in the numerical integration (default=64)"
        },
		{
        		" Local_Gaussian_ShiftedLN_Power_SpreadOption_Der",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRR",					/// 17 parametres = 16 d'entree + 1 parametre de retour 
                " ARM_CF_LNShifted_SpreadOptDer",
				" derivative rank, index1,index2,sigma1 ,alpha1,sigma2,alpha2,correlation,maturity,a1,b1,k1,a2,b2,k2,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives a first derivative of the  FORWARD VALUE of a power spread option (ShiftedLN+Gaussian)",
				" index of the dimension: 0=fwd1, 1=fwd2, 2=sigma1, 3=alpha1, 4=sigma2, 5=alpha2, 6= corr, 7=t",
				" fwd value (real value: 0.05 for 5%)",
				" fwd value (real value: 0.05 for 5%)",
				" sigma index1",
				" alpha index1",
				" sigma index2",
				" alpha index2",
				" copula correlation: 0.5 for 50%",
				" year fraction (1 for one year)",
				" a1: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" b1: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" k1: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" a2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" b2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" k2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" number of points in the numerical integration (default=64)"
        },
		{
        		" Local_Gaussian_ShiftedLN_Power_SpreadOption_Der2",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRR" ,					/// 18 parametres = 17 d'entree + 1 parametre de retour 
                " ARM_CF_LNShifted_SpreadOptDer2",
				" derivative rank 1,derivative rank 2 , index1,index2 ,sigma1 ,alpha1,sigma2,alpha2,correlation,maturity,a1,b1,k1,a2,b2,k2,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives a second derivative of the  FORWARD VALUE of a power spread option (ShiftedLN+Gaussian)",
				" index of the dimension: 0=fwd1, 1=fwd2, 2=sigma1, 3=alpha1, 4=sigma2, 5=alpha2, 6= corr, 7=t",
				" index of the dimension: 0=fwd1, 1=fwd2, 2=sigma1, 3=alpha1, 4=sigma2, 5=alpha2, 6= corr, 7=t",
				" fwd value (real value: 0.05 for 5%)",
				" fwd value (real value: 0.05 for 5%)",
				" sigma index1",
				" alpha index1",
				" sigma index2",
				" alpha index2",
				" copula correlation: 0.5 for 50%",
				" year fraction (1 for one year)",
				" a1: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" b1: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" k1: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" a2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" b2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" k2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" number of points in the numerical integration (default=64)"
        },
		{
        		" Local_Merton_JumpDiffusion",		/// name of the C++ function
                " RRRRRRRRRR" ,					/// 10 parametres = 9 d'entree + 1 parametre de retour 
                " ARM_CF_MertonJmpDiff_Option",
				" index,strike ,timetomaturity ,sigma,lambda,muJ,sigmaJ,callput,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the value of an option under the Merton Jump Diffusion Process",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" time to maturity (years)",
				" volatility (0.15 for 15% vol)",
				" probability of jump (0.1 for 10% yearly)",
				" expected value of the Log of the jump",
				" volatility of the jump (0.2 for 20% yearly)",
				" Call (CALL) or Put (PUT) ",
				" nb max of jumps taken into account (10 by default) "
        },
		{
        		" Local_Merton_JumpDiffusion_Der",		/// name of the C++ function
                " RRRRRRRRRRR" ,					/// 11 parametres = 10 d'entree + 1 parametre de retour 
                " ARM_CF_MertonJmpDiff_OptionDer",
				" index,forward,strike ,timetomaturity ,sigma,lambda,muJ,sigmaJ,callput,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the derivative of an option under the Merton Jump Diffusion Process",
				" index of the dimension: 0=fwd, 1=k, 2=t, 3=sigma, 4=lambda, 5=muJ, 6= sigmaJ",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" time to maturity (years)",
				" volatility (0.15 for 15% vol)",
				" probability of jump (0.1 for 10% yearly)",
				" expected value of the Log of the jump",
				" volatility of the jump (0.2 for 20% yearly)",
				" Call (CALL) or Put (PUT) ",
				" nb max of jumps taken into account (10 by default) "
        },
		{
        		" Local_Merton_JumpDiffusion_Der2",		/// name of the C++ function
                " RRRRRRRRRRRR" ,					/// 12 parametres = 10 d'entree + 1 parametre de retour 
                " ARM_CF_MertonJmpDiff_OptionDer2",
				" index,forward,strike ,timetomaturity ,sigma,lambda,muJ,sigmaJ,callput,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the second derivative of an option under the Merton Jump Diffusion Process",
				" index of the first dimension: 0=fwd, 1=k, 2=t, 3=sigma, 4=lambda, 5=muJ, 6= sigmaJ",
				" index of the second dimension: 0=fwd, 1=k, 2=t, 3=sigma, 4=lambda, 5=muJ, 6= sigmaJ",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" time to maturity (years)",
				" volatility (0.15 for 15% vol)",
				" probability of jump (0.1 for 10% yearly)",
				" expected value of the Log of the jump",
				" volatility of the jump (0.2 for 20% yearly)",
				" Call (CALL) or Put (PUT) ",
				" nb max of jumps taken into account (10 by default) "
        },
		{
        		" Local_SABR_ImplicitVol",		/// name of the C++ function
                " RRRRRRRRRRRRR" ,					/// 13 parametres = 12 d'entree + 1 parametre de retour 
                " ARM_CF_SABR_ImplicitVol",
				" forward,strike ,timetomaturity ,alpha,beta,rho,nu,flag,[n],[Alpha-Exp],[Alpha-Tanh],[Kb-Tanh]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the SABR implicit vol of a call",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" time to maturity (years)",
				" alpha (0.15 for 15% vol)",
				" beta",
				" rho (-0.5 for -50%)",
				" nu (0.15 for 15% vol)",
				" Type of SABR formula : SABR_IMPLNVOL or SABR_DIRECT_SC1",
				" number of Laguerre points for the integration (if necessary)",
				" Alpha for StrikeCuter Exp",
				" Alpha for StrikeCuter Tanh",
				" Kb for StrikeCuter Tanh"
        },
		{
        		" Local_SABR_ImplicitVol_Der",		/// name of the C++ function
                " RRRRRRRRRRRRRR" ,					/// 14 parametres = 13 d'entree + 1 parametre de retour 
                " ARM_CF_SABR_ImplicitVolDer",
				" index,forward,strike ,timetomaturity ,alpha,beta,rho,nu,flag,[n],[Alpha-Exp],[Alpha-Tanh],[Kb-Tanh]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the derivative of an SABR implicit vol of a call",
				" index of the dimension: 0=fwd, 1=k, 2=t, 3=alpha, 4=beta, 5=rho, 6= nu",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" time to maturity (years)",
				" alpha (0.15 for 15% vol)",
				" beta",
				" rho (-0.5 for -50%)",
				" nu (0.15 for 15% vol)",
				" Type of SABR formula : SABR_IMPLNVOL(2),SABR_A,SABR_G,DIRECTEXACT,DIRECTGEOMETRIC,DIRECTARITHMETIC,NORMALEXACT,NORMALGEOMETRIC,NORMALARITHMETIC,ANALYTICZP0,ANALYTICZP2",
				" number of Laguerre points for the integration (if necessary)",
				" Alpha for StrikeCuter Exp",
				" Alpha for StrikeCuter Tanh",
				" Kb for StrikeCuter Tanh"
        },
		{
        		" Local_SABR_ImplicitVol_Der2",		/// name of the C++ function
                " RRRRRRRRRRRRRRR" ,					/// 15 parametres = 14 d'entree + 1 parametre de retour 
                " ARM_CF_SABR_ImplicitVolDer2",
				" index1,index2 ,forward,strike ,timetomaturity ,alpha,beta,rho,nu,flag,[n],[Alpha-Exp],[Alpha-Tanh],[Kb-Tanh]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the second derivative of an SABR implicit vol of a call",
				" index of the first dimension: 0=fwd, 1=k, 2=t, 3=alpha, 4=beta, 5=rho, 6= nu",
				" index of the second dimension: 0=fwd, 1=k, 2=t, 3=alpha, 4=beta, 5=rho, 6= nu",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" time to maturity (years)",
				" alpha (0.15 for 15% vol)",
				" beta",
				" rho (-0.5 for -50%)",
				" nu (0.15 for 15% vol)",
				" Type of SABR formula : SABR_IMPLNVOL(2),SABR_A,SABR_G,DIRECTEXACT,DIRECTGEOMETRIC,DIRECTARITHMETIC,NORMALEXACT,NORMALGEOMETRIC,NORMALARITHMETIC,ANALYTICZP0,ANALYTICZP2",
				" number of Laguerre points for the integration (if necessary)",
				" Alpha for StrikeCuter Exp",
				" Alpha for StrikeCuter Tanh",
				" Kb for StrikeCuter Tanh"
        },
		{
        		" Local_SABR_FromGaussianToDistribution",		/// name of the C++ function
                " RRRRRRRRRRRRR" ,					/// 10 parametres = 9 d'entree + 1 parametre de retour 
                " ARM_CF_SABR_FromGaussianToDistribution",
				" forward,strike ,timetomaturity ,alpha,beta,rho,nu,flag,[n],[a_exp],[a_tanh],[k_tanh]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the SABR value of the forward value of a vanilla option",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" time to maturity (years)",
				" alpha (0.15 for 15% vol)",
				" beta",
				" rho (-0.5 for -50%)",
				" nu (0.15 for 15% vol)",
				" Type of SABR formula : SABR_IMPLNVOL(2),SABR_A,SABR_G,DIRECTEXACT,DIRECTGEOMETRIC,DIRECTARITHMETIC,NORMALEXACT,NORMALGEOMETRIC,NORMALARITHMETIC,ANALYTICZP0,ANALYTICZP2",
				" number of Laguerre points for the integration (if necessary)"
				" Strike Cuter coef exp ",
				" Strike Cuter coef tanh ",
				" Strike Cuter coef cut off KB "
        },
		{
        		" Local_SABR_FromDistributionToGaussian",		/// name of the C++ function
                " RRRRRRRRRRRRR" ,					/// 10 parametres = 19 d'entree + 1 parametre de retour 
                " ARM_CF_SABR_FromDistributionToGaussian",
				" forward,strike ,timetomaturity ,alpha,beta,rho,nu,flag,[n],[a_exp],[a_tanh],[k_tanh]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the SABR value of the forward value of a vanilla option",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" time to maturity (years)",
				" alpha (0.15 for 15% vol)",
				" beta",
				" rho (-0.5 for -50%)",
				" nu (0.15 for 15% vol)",
				" Type of SABR formula : SABR_IMPLNVOL(2),SABR_A,SABR_G,DIRECTEXACT,DIRECTGEOMETRIC,DIRECTARITHMETIC,NORMALEXACT,NORMALGEOMETRIC,NORMALARITHMETIC,ANALYTICZP0,ANALYTICZP2",
				" number of Laguerre points for the integration (if necessary)",
				" Strike Cuter coef exp ",
				" Strike Cuter coef tanh ",
				" Strike Cuter coef cut off KB "
        },
		{
        		" Local_SABR_VanillaOption",		/// name of the C++ function
                " RRRRRRRRRRRRRR" ,					/// 14 parametres = 13 d'entree + 1 parametre de retour 
                " ARM_CF_SABR_Option",
				" forward,strike ,timetomaturity ,alpha,beta,rho,nu,callput,flag,[n],[Alpha-Exp],[Alpha-Tanh],[Kb-Tanh]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the SABR value of the forward value of a vanilla option",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" time to maturity (years)",
				" alpha (0.15 for 15% vol)",
				" beta",
				" rho (-0.5 for -50%)",
				" nu (0.15 for 15% vol)",
				" CALL or PUT",
				" Type of SABR formula : SABR_IMPLNVOL(2),SABR_A,SABR_G,DIRECTEXACT,DIRECTGEOMETRIC,DIRECTARITHMETIC,NORMALEXACT,NORMALGEOMETRIC,NORMALARITHMETIC,ANALYTICZP0,ANALYTICZP2",
				" number of Laguerre points for the integration (if necessary)",
				" Alpha for StrikeCuter Exp",
				" Alpha for StrikeCuter Tanh",
				" Kb for StrikeCuter Tanh"
        },
		{
        		" Local_SABR_VanillaOption_Der",		/// name of the C++ function
                " RRRRRRRRRRRRRRR" ,					/// 15 parametres = 14 d'entree + 1 parametre de retour 
                " ARM_CF_SABR_OptionDer",
				" index,forward,strike ,timetomaturity ,alpha,beta,rho,nu,callput,flag,[n],[Alpha-Exp],[Alpha-Tanh],[Kb-Tanh]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the derivative of  the forward value of a vanilla option",
				" index of the dimension: 0=fwd, 2=k, 2=t, 3=alpha, 4=beta, 5=rho, 6= nu",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" time to maturity (years)",
				" alpha (0.15 for 15% vol)",
				" beta",
				" rho (-0.5 for -50%)",
				" nu (0.15 for 15% vol)",
				" CALL or PUT",
				" Type of SABR formula : SABR_IMPLNVOL(2),SABR_A,SABR_G,DIRECTEXACT,DIRECTGEOMETRIC,DIRECTARITHMETIC,NORMALEXACT,NORMALGEOMETRIC,NORMALARITHMETIC,ANALYTICZP0,ANALYTICZP2",
				" number of Laguerre points for the integration (if necessary)",
				" Alpha for StrikeCuter Exp",
				" Alpha for StrikeCuter Tanh",
				" Kb for StrikeCuter Tanh"
        },
		{
        		" Local_SABR_VanillaOption_Der2",		/// name of the C++ function
                " RRRRRRRRRRRRRRRR" ,					/// 16 parametres = 15 d'entree + 1 parametre de retour 
                " ARM_CF_SABR_OptionDer2",
				" index1,index2 ,forward,strike ,timetomaturity ,alpha,beta,rho,nu,callput,flag,[n],[Alpha-Exp],[Alpha-Tanh],[Kb-Tanh]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the second derivative of  the forward value of a vanilla option",
				" index of the first dimension: 0=fwd, 2=k, 2=t, 3=alpha, 4=beta, 5=rho, 6= nu",
				" index of the second dimension: 0=fwd, 2=k, 2=t, 3=alpha, 4=beta, 5=rho, 6= nu",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" time to maturity (years)",
				" alpha (0.15 for 15% vol)",
				" beta",
				" rho (-0.5 for -50%)",
				" nu (0.15 for 15% vol)",
				" CALL or PUT",
				" Type of SABR formula : SABR_IMPLNVOL(2),SABR_A,SABR_G,DIRECTEXACT,DIRECTGEOMETRIC,DIRECTARITHMETIC,NORMALEXACT,NORMALGEOMETRIC,NORMALARITHMETIC,ANALYTICZP0,ANALYTICZP2",
				" number of Laguerre points for the integration (if necessary)",
				" Alpha for StrikeCuter Exp",
				" Alpha for StrikeCuter Tanh",
				" Kb for StrikeCuter Tanh"
        },
		{
        		" Local_SABR_From_Sigma_To_Alpha",		/// name of the C++ function
                " RRRRRRRRRR" ,					/// 10 parametres = 9 d'entree + 1 parametre de retour 
                " ARM_CF_SABR_FromSigmaToAlpha",
				" forward,strike ,timetomaturity ,sigma,beta,rho,nu,flag,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the alpha associated with a given implicit volatility",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" time to maturity (years)",
				" sigma (0.15 for 15% vol)",
				" beta",
				" rho (-0.5 for -50%)",
				" nu (0.15 for 15% vol)",
				" Type of SABR formula : SABR_IMPLNVOL(2),SABR_A,SABR_G,DIRECTEXACT,DIRECTGEOMETRIC,DIRECTARITHMETIC,NORMALEXACT,NORMALGEOMETRIC,NORMALARITHMETIC,ANALYTICZP0,ANALYTICZP2",
				" number of terms for the complex N[x] (if necessary)"
        },
		{
        		" Local_BS_EuroBarriere",		/// name of the C++ function
                " RRRRRRRRRRR" ,					/// 11 parametres = 10 d'entree + 1 parametre de retour 
                " ARM_CF_LN_BarrierOption",
				" forward, strike, barrier, rebate, volatility, expiry, discount,callput, in or out, up or down",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the forward value of a barrier option with a logNormal model",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" barriere (real value: 0.05 for 5%)",
				" rebate (real value: 0.05 for 5%)",
				" volatility (0.15 for 15% vol)",
				" expiry (years: 0.5 for 6m)",
				" discount (0.95 for exemple)",
				" CALL or PUT",
				" IN or OUT Barrier ",
				" Up or DOWN Barrier"
        },
		{
        		" Local_BS_EuroBarriere_Der",		/// name of the C++ function
                " RRRRRRRRRRRR" ,					/// 12 parametres = 11 d'entree + 1 parametre de retour 
                " ARM_CF_LN_BarrierOptionDer",
				" index,forward, strike, barrier, rebate, volatility, expiry,discount, callput, in or out, up or down",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the derivative ofthe forward value of a barrier option with a LogNormal model",
				" index of the dimension: 0=fwd, 1=k, 2=b, 3=r, 4=v, 5=t",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" barriere (real value: 0.05 for 5%)",
				" rebate (real value: 0.05 for 5%)",
				" volatility (0.15 for 15% vol)",
				" expiry (years: 0.5 for 6m)",
				" discount (0.95 for exemple)",
				" CALL or PUT",
				" IN or OUT Barrier ",
				" Up or DOWN Barrier"
        },
		{
        		" Local_BS_EuroBarriere_Der2",		/// name of the C++ function
                " RRRRRRRRRRRRR" ,					/// 13 parametres = 12 d'entree + 1 parametre de retour 
                " ARM_CF_LN_BarrierOptionDer2",
				" index1,index2 ,forward, strike, barrier, rebate, volatility, expiry,discount, callput, in or out, up or down",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the second derivative of  the forward value of a barrier option with a LogNormal model",
				" index1 of the first dimension: 0=fwd, 1=k, 2=b, 3=r, 4=v, 5=t",
				" index2 of the second dimension: 0=fwd, 1=k, 2=b, 3=r, 4=v, 5=t",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" barriere (real value: 0.05 for 5%)",
				" rebate (real value: 0.05 for 5%)",
				" volatility (0.15 for 15% vol)",
				" expiry (years: 0.5 for 6m)",
				" discount (0.95 for exemple)",
				" CALL or PUT",
				" IN or OUT Barrier ",
				" Up or DOWN Barrier"
        },
		{
        		" Local_BS_EuroDoubleBarriere",		/// name of the C++ function
                " RRRRRRRRRR" ,					/// 10 parametres = 9 d'entree + 1 parametre de retour 
                " ARM_CF_LN_DoubleBarrierOption",
				" forward,strike ,barrierUp ,barrierDown, volatility, expiry,interestrate,dividendyield, callput",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the forward value of a option barrier option with a LogNormal model",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" barriere up (real value: 0.05 for 5%)",
				" barriere down (real value: 0.05 for 5%)",
				" volatility (0.15 for 15% vol)",
				" expiry (years: 0.5 for 6m)",
				" IR (0.05 for 05%)",
				" Div (0.01 for 1%)",
				" CALL or PUT"
        },
		{
        		" Local_BS_EuroDoubleBarriere_Der",		/// name of the C++ function
                " RRRRRRRRR" ,					/// 9 parametres = 8 d'entree + 1 parametre de retour 
                " ARM_CF_LN_DoubleBarrierOptionDer",
				" index,forward,strike ,barrierUp ,barrierDown,volatility,expiry,callput",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the derivative of the forward value of a double barrier option with a LogNormal model",
				" index of the dimension: 0=fwd, 1=k, 2=bup, 3=bdown, 4=v, 5=t",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" barriere up (real value: 0.05 for 5%)",
				" barriere down (real value: 0.05 for 5%)",
				" volatility (0.15 for 15% vol)",
				" expiry (years: 0.5 for 6m)",
				" CALL or PUT"
        },
		{
        		" Local_BS_EuroDoubleBarriere_Der2",		/// name of the C++ function
                " RRRRRRRRRR" ,					/// 10 parametres = 8 d'entree + 1 parametre de retour 
                " ARM_CF_LN_DoubleBarrierOptionDer2",
				" index1,index2 ,forward,strike ,barrierUp ,barrierDown,volatility,expiry,callput",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the second derivative of the forward value of a double barrier option with a LogNormal model",
				" index1 of the first dimension: 0=fwd, 1=k, 2=bup, 3=bdown, 4=v, 5=t",
				" index2 of the second dimension: 0=fwd, 1=k, 2=bup, 3=bdown, 4=v, 5=t",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" barrier up (real value: 0.05 for 5%)",
				" barrier down (real value: 0.05 for 5%)",
				" volatility (0.15 for 15% vol)",
				" expiry (years: 0.5 for 6m)",
				" CALL or PUT"
        },
		{
        		" Local_BS_PartialTime_Start_SingleBarrier",		/// name of the C++ function
                " RRRRRRRRRR" ,					/// 10 parametres = 9 d'entree + 1 parametre de retour 
                " ARM_CF_LN_PartialTime_Start_SingleBarrier",
				" forward,strike,barriere,rebate,volatility,T_end_barrier,T_maturity,callput,option_type",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the forward value of a single european barriere option, whose barrier end before maturity with a Black and Sholes model",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" barriere (real value: 0.05 for 5%)",
				" rebate (real value: 0.05 for 5%)",
				" volatility (0.15 for 15% vol)",
				" end of the barrier (years)",
				" time to maturity (years)",
				" CALL or PUT",
				" Option Type: DOWN_AND_IN,UP_AND_IN,DOWN_AND_OUT,UP_AND_OUT"
        },
		{
        		" Local_BS_PartialTime_Start_SingleBarrier_Der",		/// name of the C++ function
                " RRRRRRRRRRR" ,					/// 11 parametres = 10 d'entree + 1 parametre de retour 
                " ARM_CF_LN_PartialTime_Start_SingleBarrierDer",
				" index,forward,strike,barriere,rebate,volatility,T_end_barrier,T_maturity,callput,option_type",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the derivative of the forward value of a single european barriere option, whose barrier end before maturity with a Black and Sholes model",
				" index: 0=fwd,1=strike,2=barrier,3=rebate,4=volatility,5=TendBarr,6=Tmaturity",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" barriere (real value: 0.05 for 5%)",
				" rebate (real value: 0.05 for 5%)",
				" volatility (0.15 for 15% vol)",
				" end of the barrier (years)",
				" time to maturity (years)",
				" CALL or PUT",
				" Option Type: DOWN_AND_IN,UP_AND_IN,DOWN_AND_OUT,UP_AND_OUT"
        },
		{
        		" Local_BS_PartialTime_Start_SingleBarrier_Der2",		/// name of the C++ function
                " RRRRRRRRRRRR" ,					/// 12 parametres = 11 d'entree + 1 parametre de retour 
                " ARM_CF_LN_PartialTime_Start_SingleBarrierDer2",
				" index,index2,forward,strike,barriere,rebate,volatility,T_end_barrier,T_maturity,callput,option_type",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the 2nd derivative of the forward value of a single european barriere option, whose barrier end before maturity with a Black and Sholes model",
				" index1: 0=fwd,1=strike,2=barrier,3=rebate,4=volatility,5=TendBarr,6=Tmaturity",
				" index2: 0=fwd,1=strike,2=barrier,3=rebate,4=volatility,5=TendBarr,6=Tmaturity",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" barriere (real value: 0.05 for 5%)",
				" rebate (real value: 0.05 for 5%)",
				" volatility (0.15 for 15% vol)",
				" end of the barrier (years)",
				" time to maturity (years)",
				" CALL or PUT",
				" Option Type: DOWN_AND_IN,UP_AND_IN,DOWN_AND_OUT,UP_AND_OUT"
        },
		{
        		" Local_BS_PartialTime_End_SingleBarrier",		/// name of the C++ function
                " RRRRRRRRRR" ,					/// 10 parametres = 9 d'entree + 1 parametre de retour 
                " ARM_CF_LN_PartialTime_End_SingleBarrier",
				" forward,strike,barriere,rebate,volatility,T_start_barrier,T_maturity,callput,option_type",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the forward value of a single european barriere option, whose barrier end before maturity with a Black and Sholes model",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" barriere (real value: 0.05 for 5%)",
				" rebate (real value: 0.05 for 5%)",
				" volatility (0.15 for 15% vol)",
				" start of the barrier (years)",
				" time to maturity (years)",
				" CALL or PUT",
				" Option Type: CROSS_AND_IN,CROSS_AND_OUT,INSIDE_UP_AND_IN,INSIDE_UP_AND_OUT,INSIDE_DOWN_AND_IN,INSIDE_DOWN_AND_OUT"
        },
		{
        		" Local_BS_PartialTime_End_SingleBarrier_Der",		/// name of the C++ function
                " RRRRRRRRRRR" ,					/// 11 parametres = 10 d'entree + 1 parametre de retour 
                " ARM_CF_LN_PartialTime_End_SingleBarrierDer",
				" index,forward,strike,barriere,rebate,volatility,T_start_barrier,T_maturity,callput,option_type",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the forward value of a single european barriere option, whose barrier end before maturity with a Black and Sholes model",
				" index: 0=fwd,1=strike,2=barrier,3=rebate,4=volatility,5=TendBarr,6=Tmaturity",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" barriere (real value: 0.05 for 5%)",
				" rebate (real value: 0.05 for 5%)",
				" volatility (0.15 for 15% vol)",
				" start of the barrier (years)",
				" time to maturity (years)",
				" CALL or PUT",
				" Option Type: CROSS_AND_IN,CROSS_AND_OUT,INSIDE_UP_AND_IN,INSIDE_UP_AND_OUT,INSIDE_DOWN_AND_IN,INSIDE_DOWN_AND_OUT"
        },
		{
        		" Local_BS_PartialTime_End_SingleBarrier_Der2",		/// name of the C++ function
                " RRRRRRRRRRRR" ,					/// 12 parametres = 11 d'entree + 1 parametre de retour 
                " ARM_CF_LN_PartialTime_End_SingleBarrierDer2",
				" index,forward,strike,barriere,rebate,volatility,T_start_barrier,T_maturity,callput,option_type",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the forward value of a single european barriere option, whose barrier end before maturity with a Black and Sholes model",
				" index1: 0=fwd,1=strike,2=barrier,3=rebate,4=volatility,5=TendBarr,6=Tmaturity",
				" index2: 0=fwd,1=strike,2=barrier,3=rebate,4=volatility,5=TendBarr,6=Tmaturity",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" barriere (real value: 0.05 for 5%)",
				" rebate (real value: 0.05 for 5%)",
				" volatility (0.15 for 15% vol)",
				" start of the barrier (years)",
				" time to maturity (years)",
				" CALL or PUT",
				" Option Type: CROSS_AND_IN,CROSS_AND_OUT,INSIDE_UP_AND_IN,INSIDE_UP_AND_OUT,INSIDE_DOWN_AND_IN,INSIDE_DOWN_AND_OUT"
        },
		{
        		" Local_BS_SingleBarrier_2Asset",		/// name of the C++ function
                " RRRRRRRRRRR" ,					/// 11 parametres = 10 d'entree + 1 parametre de retour 
                " ARM_CF_LN_SingleBarrier_2Assets",
				" forward1,strike1,forward2,strike2,volatility1,volatility2,correlation,T_maturity,callput,option_type",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the forward value of a single european barriere option, paid with and option on another asset, with a Black and Sholes model",
				" fwd1 value (real value: 0.05 for 5%)",
				" strike1 (real value: 0.05 for 5%)",
				" fwd2 value (real value: 0.05 for 5%)",
				" strike2 (real value: 0.05 for 5%)",
				" volatility1 (0.15 for 15% vol)",
				" volatility2 (0.15 for 15% vol)",
				" end of the barrier (years)",
				" time to maturity (years)",
				" CALL or PUT",
				" Option Type: DOWN_AND_IN,UP_AND_IN,DOWN_AND_OUT,UP_AND_OUT"
        },
		{
        		" Local_BS_SingleBarrier_2Asset_Der",		/// name of the C++ function
                " RRRRRRRRRRRR" ,					/// 12 parametres = 11 d'entree + 1 parametre de retour 
                " ARM_CF_LN_SingleBarrier_2AssetsDer",
				" Index,forward1,strike1,forward2,strike2,volatility1,volatility2,correlation,T_maturity,callput,option_type",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the forward value of a single european barriere option, paid with and option on another asset, with a Black and Sholes model",
				"index: 0=fwd1,1=strike1,2=fwd2,3=strike2,4=vol1,5=vol2,6=corr,7=Tmaturity",
				" fwd1 value (real value: 0.05 for 5%)",
				" strike1 (real value: 0.05 for 5%)",
				" fwd2 value (real value: 0.05 for 5%)",
				" strike2 (real value: 0.05 for 5%)",
				" volatility1 (0.15 for 15% vol)",
				" volatility2 (0.15 for 15% vol)",
				" end of the barrier (years)",
				" time to maturity (years)",
				" CALL or PUT",
				" Option Type: DOWN_AND_IN,UP_AND_IN,DOWN_AND_OUT,UP_AND_OUT"
        },
		{
        		" Local_BS_SingleBarrier_2Asset_Der2",		/// name of the C++ function
                " RRRRRRRRRRRRR" ,					/// 13 parametres = 12 d'entree + 1 parametre de retour 
                " ARM_CF_LN_SingleBarrier_2AssetsDer2",
				" Index,Index2,forward1,strike1,forward2,strike2,volatility1,volatility2,correlation,T_maturity,callput,option_type",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the forward value of a single european barriere option, paid with and option on another asset, with a Black and Sholes model",
				" index1: 0=fwd1,1=strike1,2=fwd2,3=strike2,4=vol1,5=vol2,6=corr,7=Tmaturity",
				" index2: 0=fwd1,1=strike1,2=fwd2,3=strike2,4=vol1,5=vol2,6=corr,7=Tmaturity",
				" fwd1 value (real value: 0.05 for 5%)",
				" strike1 (real value: 0.05 for 5%)",
				" fwd2 value (real value: 0.05 for 5%)",
				" strike2 (real value: 0.05 for 5%)",
				" volatility1 (0.15 for 15% vol)",
				" volatility2 (0.15 for 15% vol)",
				" end of the barrier (years)",
				" time to maturity (years)",
				" CALL or PUT",
				" Option Type: DOWN_AND_IN,UP_AND_IN,DOWN_AND_OUT,UP_AND_OUT"
        },
		{
        		" Local_Bivariate",		/// name of the C++ function
                " RRRRRR" ,					/// 4 parametres = 3 d'entree + 1 parametre de retour 
                " ARM_CF_Util_Bivariate",
				" x,y,correlation,[p],[q]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives E[ X^p Y^q 1(X<x) 1(Y<y) ], and therefore the bivariate cdf if p = q = 0",
				" x",
				" y",
				" correlation (between -1 and 1 )",
				" p",
				" q"
        },
		{
        		" Local_Gamma",		/// name of the C++ function
                " RR" ,					/// 2 parametres = 1 d'entree + 1 parametre de retour 
                " ARM_CF_Util_Gamma",
				" x",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the gamma function",
				" x"
        },
		{
        		" Local_ImcompleteBeta",		/// name of the C++ function
                " RRRR" ,					/// 4 parametres = 3 d'entree + 1 parametre de retour 
                " ARM_CF_Util_IncompleteBeta",
				" a , b, x",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the IncompleteBeta function",
				" a>0",
				" b>0",
				" x should be between 0 and 1"
        },
		{
        		" Local_InverseImcompleteBeta",		/// name of the C++ function
                " RRRRR" ,					/// 4 parametres = 3 d'entree + 1 parametre de retour 
                " ARM_CF_Util_InverseIncompleteBeta",
				" a , b, z0, x",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the IncompleteBeta function",
				" a>0",
				" b>0",
				"z0",
				" x should be between 0 and 1"
        },
		{
        		" Local_Hypergeometric2F1",		/// name of the C++ function
                " RRRRR" ,					/// 4 parametres = 3 d'entree + 1 parametre de retour 
                " ARM_CF_Util_Hypergeometric2F1",
				" a , b, c, x",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the Hypergeometric2F1 function",
				" a",
				" b",
				" c>0",
				" x should be between -1 and 1"
        },
		{
        		" Local_ImaginaryPart_ComplexErf",		/// name of the C++ function
                " RRRR" ,					/// 4 parametres = 3 d'entree + 1 parametre de retour 
                " ARM_CF_Util_ImCmpxErf",
				" real part,imaginary part ,nb terms",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the imaginary part of the erf function extended to a complex argument",
				" real part of the argument",
				" imaginary part of the argument",
				" nb terms (usually between 5 and 10)"
        },
		{
        		" Local_RealPart_ComplexErf",		/// name of the C++ function
                " RRRR" ,					/// 4 parametres = 3 d'entree + 1 parametre de retour 
                " ARM_CF_Util_ReCmpxErf",
				" real part,imaginary part ,nb terms",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the real part of the erf function extended to a complex argument",
				" real part of the argument",
				" imaginary part of the argument",
				" nb terms (usually between 5 and 10)"
        },
		{
        		" Local_NormalCummulative_Inverse",		/// name of the C++ function
                " RR" ,					/// 2 parametres = 1 d'entree + 1 parametre de retour 
                " ARM_CF_Util_NormalInverse",
				" x",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the inverse of the normal function",
				" argument",
        },
		{
        		" Local_NormalCummulative",		/// name of the C++ function
                " RR" ,					/// 2 parametres = 1 d'entree + 1 parametre de retour 
                " ARM_CF_Util_Normal",
				" x",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the normal function",
				" argument",
        },
		{
        		" Local_GHeston_VanillaOption",		/// name of the C++ function
                " RRRRRRRRRRRRRR" ,					/// 14 parametres = 12 d'entree + 1 parametre de retour 
                " ARM_CF_GHeston_Option",
				" forward,strike ,Vol initiale,timetomaturity ,long term limit,mean reveting speed,vol of vol,correlation,jump probability, jump size,jump volatility, callput,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the Generalized Heston value of the forward value of a vanilla option",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" vol initiale (real value: 0.05 for 5%)",
				" time to maturity (years)",
				" long terme limit (0.15 for 15% vol)",
				" mean reverting speed",
				" volatility of volatility",
				" correlation between the volatility and the underlying (-0.5 for -50%)",
				" jump probability (0.15 for 15%)",
				" jump size (expectation of the return)",
				" jump volatility ",
				" CALL or PUT",
				" number of Legendre points for the integration "
        },
		{
        		" Local_GHeston_VanillaOption_Der",		/// name of the C++ function
                " RRRRRRRRRRRRRRR" ,					/// 15 parametres = 13 d'entree + 1 parametre de retour 
                " ARM_CF_GHeston_OptionDer",
				" index,forward,strike ,Vol initiale,timetomaturity ,long term limit,mean reveting speed,vol of vol,correlation,jump probability, jump size,jump volatility, callput,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the derivative of the Generalized Heston forward value of a vanilla option",
				" index of the dimension: 0=fwd, 1=k,2=v0, 3=t, 4=longterm, 5=speed,6=volvol, 7=rho, 8= lambda,9=muJ,10=sigmaJ",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" vol initiale (real value: 0.05 for 5%)",
				" time to maturity (years)",
				" long terme limit (0.15 for 15% vol)",
				" mean reverting speed",
				" volatility of volatility",
				" correlation between the volatility and the underlying (-0.5 for -50%)",
				" jump probability (0.15 for 15%)",
				" jump size (expectation of the return)",
				" jump volatility ",
				" CALL or PUT",
				" number of Legendre points for the integration "
        },
		{
        		" Local_GHeston_VanillaOption_Der2",		/// name of the C++ function
                " RRRRRRRRRRRRRRRR" ,					/// 16 parametres = 14 d'entree + 1 parametre de retour 
                " ARM_CF_GHeston_OptionDer2",
				" index1,index2 ,forward,strike ,Vol initiale,timetomaturity ,long term limit,mean reveting speed,vol of vol,correlation,jump probability, jump size,jump volatility, callput,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the second derivative of the Generalized Heston forward value of a vanilla option",
				" index of the dimension: 0=fwd, 1=k,2=v0, 3=t, 4=longterm, 5=speed,6=volvol, 7=rho, 8= lambda,9=muJ,10=sigmaJ",
				" index of the dimension: 0=fwd, 1=k,2=v0, 3=t, 4=longterm, 5=speed,6=volvol, 7=rho, 8= lambda,9=muJ,10=sigmaJ",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" vol initiale (real value: 0.05 for 5%)",
				" time to maturity (years)",
				" long terme limit (0.15 for 15% vol)",
				" mean reverting speed",
				" volatility of volatility",
				" correlation between the volatility and the underlying (-0.5 for -50%)",
				" jump probability (0.15 for 15%)",
				" jump size (expectation of the return)",
				" jump volatility ",
				" CALL or PUT",
				" number of Legendre points for the integration "
        },
		{
        		" Local_CEV_VanillaOption",		/// name of the C++ function
                " RRRRRRRRR" ,					/// 9 parametres = 8 d'entree + 1 parametre de retour 
                " ARM_CF_CEV_VanillaOption",
				" forward,strike,time_to_maturity,drift,volatility,beta,callput,[nbterms]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the forward value of a european vanilla option with a CEV model",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" time to reset of the option",
				" drift (necessary, put 0.00001 if 0)",
				" volatility (0.15 for 15% vol)",
				" beta (act on the slope of the smile)",
				" CALL or PUT",
				" nb of bessel terms (typically 40)"
        },
		{
        		" Local_CEV_VanillaOption_Der",		/// name of the C++ function
                " RRRRRRRRRR" ,					/// 9 parametres = 8 d'entree + 1 parametre de retour 
                " ARM_CF_CEV_VanillaOptionDer",
				" index,forward,strike,time_to_maturity,drift,volatility,beta,callput,[nbterms]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives a derivative of a forward value of a european vanilla option with a CEV model",
				" index of the dimension: 0=t, 1=fwd, 2=strike, 3=sig, 4=beta, 5=mu",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" time to reset of the option",
				" volatility (0.15 for 15% vol)",
				" beta (act on the slope of the smile)",
				" CALL or PUT",
				" nb of bessel terms (typically 40)"
        },
		{
        		" Local_CEV_VanillaOption_Der2",		/// name of the C++ function
                " RRRRRRRRRRR" ,					/// 9 parametres = 8 d'entree + 1 parametre de retour 
                " ARM_CF_CEV_VanillaOptionDer2",
				" index1,index2,forward,strike,time_to_maturity,drift,volatility,beta,callput,[nbterms]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives a second derivative of a forward value of a european vanilla option with a CEV model",
				" index of the dimension1: 0=t, 1=fwd, 2=strike, 3=sig, 4=beta, 5=mu",
				" index of the dimension2: 0=t, 1=fwd, 2=strike, 3=sig, 4=beta, 5=mu",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" time to reset of the option",
				" volatility (0.15 for 15% vol)",
				" beta (act on the slope of the smile)",
				" CALL or PUT",
				" nb of bessel terms (typically 40)"
        },
		{
        		" Local_CEV_DoubleBarrierOption",		/// name of the C++ function
                " RRRRRRRRRRR" ,					/// 11 parametres = 10 d'entree + 1 parametre de retour 
                " ARM_CF_CEV_DoubleBarrierOption",
				" forward,strike,time_to_maturity,Barrierdown,Barrierup,drift,volatility,beta,callput,[nbterms]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the forward value of a  vanilla DoubleBarrier Option with a CEV model",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" time to reset of the option",
				" Barrier Down :",
				" Barrier up :",
				" Drift ",
				" volatility (0.15 for 15% vol)",
				" beta (act on the slope of the smile)",
				" CALL or PUT",
				" nb of spectral terms (typically 6)"
        },
		{
        		" Local_CEV_DoubleBarrierOption_Der",		/// name of the C++ function
                " RRRRRRRRRRRR" ,					/// 11 parametres = 10 d'entree + 1 parametre de retour 
                " ARM_CF_CEV_DoubleBarrierOptionDer",
				" index,forward,strike,time_to_maturity,Barrierdown,Barrierup,drift,volatility,beta,callput,[nbterms]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives a derivative of the forward value of a  vanilla DoubleBarrier Option with a CEV model",
				" index of the dimension1: 0=t, 1=fwd, 2=strike, 3=sig, 4=beta, 5=mu, 6=Barrierdown, 7=Barrierup",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" time to reset of the option",
				" Barrier Down :",
				" Barrier up :",
				" Drift ",
				" volatility (0.15 for 15% vol)",
				" beta (act on the slope of the smile)",
				" CALL or PUT",
				" nb of spectral terms (typically 6)"
        },
		{
        		" Local_CEV_DoubleBarrierOption_Der2",		/// name of the C++ function
                " RRRRRRRRRRRRR" ,					/// 11 parametres = 10 d'entree + 1 parametre de retour 
                " ARM_CF_CEV_DoubleBarrierOptionDer2",
				" index,index2,forward,strike,time_to_maturity,Barrierdown,Barrierup,drift,volatility,beta,callput,[nbterms]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives a second derivative of the forward value of a  vanilla DoubleBarrier Option with a CEV model",
				" index of the dimension1: 0=t, 1=fwd, 2=strike, 3=sig, 4=beta, 5=mu, 6=Barrierdown, 7=Barrierup",
				" index of the dimension2: 0=t, 1=fwd, 2=strike, 3=sig, 4=beta, 5=mu, 6=Barrierdown, 7=Barrierup",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" time to reset of the option",
				" Barrier Down :",
				" Barrier up :",
				" Drift ",
				" volatility (0.15 for 15% vol)",
				" beta (act on the slope of the smile)",
				" CALL or PUT",
				" nb of spectral terms (typically 6)"
        },
		{
        		" Local_CEV_BarrierOption",		/// name of the C++ function
                " RRRRRRRRRRR" ,					/// 12 parametres = 11 d'entree + 1 parametre de retour 
                " ARM_CF_CEV_BarrierOption",
				" forward,strike,time_to_maturity,Barrier,drift,volatility,beta,optype,callput,[nbterms]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the forward value of a  Barrier Option with a CEV model",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" time to reset of the option",
				" Barrier ",
				" Drift ",
				" volatility (0.15 for 15% vol)",
				" beta (act on the slope of the smile)",
				" optype should : DOWN_AND_IN,UP_AND_IN,DOWN_AND_OUT,UP_AND_OUT",
				" CALL or PUT",
				" nb of spectral terms (typically 50)"
        },
		{
        		" Local_CEV_BarrierOption_Der",		/// name of the C++ function
                " RRRRRRRRRRRR" ,					/// 13 parametres = 12 d'entree + 1 parametre de retour 
                " ARM_CF_CEV_BarrierOptionDer",
				" index,forward,strike,time_to_maturity,Barrier,drift,volatility,beta,optype,callput,[nbterms]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the forward value of a  Barrier Option with a CEV model",
				" index of the dimension1: 0=t, 1=fwd, 2=strike, 3=sig, 4=beta, 5=mu, 6=Barrier",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" time to reset of the option",
				" Barrier ",
				" Drift ",
				" volatility (0.15 for 15% vol)",
				" beta (act on the slope of the smile)",
				" optype should : DOWN_AND_IN,UP_AND_IN,DOWN_AND_OUT,UP_AND_OUT",
				" CALL or PUT",
				" nb of spectral terms (typically 50)"
        },
		{
        		" Local_CEV_BarrierOption_Der2",		/// name of the C++ function
                " RRRRRRRRRRRRRR" ,					/// 14 parametres = 13 d'entree + 1 parametre de retour 
                " ARM_CF_CEV_BarrierOptionDer2",
				" index,index2,forward,strike,time_to_maturity,Barrier,drift,volatility,beta,optype,callput,[nbterms]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the forward value of a  Barrier Option with a CEV model",
				" index of the dimension1: 0=t, 1=fwd, 2=strike, 3=sig, 4=beta, 5=mu, 6=Barrier",
				" index of the dimension2: 0=t, 1=fwd, 2=strike, 3=sig, 4=beta, 5=mu, 6=Barrier",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" time to reset of the option",
				" Barrier ",
				" Drift ",
				" volatility (0.15 for 15% vol)",
				" beta (act on the slope of the smile)",
				" optype should : DOWN_AND_IN,UP_AND_IN,DOWN_AND_OUT,UP_AND_OUT",
				" CALL or PUT",
				" nb of spectral terms (typically 50)"
        },
		{
        		" Local_LogNormal_TriSpreadOption",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRR",					/// 19 parametres = 18 d'entree + 1 parametre de retour 
                " ARM_CF_LN_TriSpreadOpt",
				" index1 value,index2 value,index3 value,sigma1,sigma2,sigma3,correlation12,correlation13,correlation23,mu1,mu2,mu3,a0,a1,a2,a3,maturity,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives  the FORWARD VALUE of a tri index  spread option (3Lognormals)",
				" fwd value 1 (real value: 0.05 for 5%)",
				" fwd value 2 (real value: 0.05 for 5%)",
				" fwd value 3 (real value: 0.05 for 5%)",
				" LN vol 1 (real value: 0.1 for 10%)",
				" LN vol 2 (real value: 0.1 for 10%)",
				" LN vol 3 (real value: 0.1 for 10%)",
				" correlation 12 (real value: 0.1 for 10%)",
				" correlation 13 (real value: 0.1 for 10%)",
				" correlation 23 (real value: 0.1 for 10%)",
				" mu1 : drift de l'index 1",
				" mu2 : drift de l'index 2",
				" mu3 : drift de l'index 3",
				" coefficient 0 : (moins strike): 0.05 for 5%",
				" coefficient 1 : ",
				" coefficient 2 : ",
				" coefficient 3 : ",
				" year fraction (1 for one year)",
				" number of points in the numerical integration (default=64)"
        },
		{
        		" Local_LogNormal_TriSpreadOption_Der",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRR",					/// 20 parametres = 19 d'entree + 1 parametre de retour 
                " ARM_CF_LN_TriSpreadOptDer",
				" derivative index1, index1 value,index2 value,index3 value,sigma1,sigma2,sigma3,correlation12,correlation13,correlation23,mu1,mu2,mu3,a0,a1,a2,a3,maturity,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives a derivative of the FORWARD VALUE of a tri index  spread option (3Lognormals)",
				" index of the dimension: 0=fwd1, 1=fwd2, 2=fwd3, 3=sig1, 4=sig2, 5=sig3, 6=rho12, 7=rho13,8=rho23,9=mu1,10=mu2,11=mu3,12=a0,13=a1,14=a2,15=a3,16=t",
				" fwd value 1 (real value: 0.05 for 5%)",
				" fwd value 2 (real value: 0.05 for 5%)",
				" fwd value 3 (real value: 0.05 for 5%)",
				" LN vol 1 (real value: 0.1 for 10%)",
				" LN vol 2 (real value: 0.1 for 10%)",
				" LN vol 3 (real value: 0.1 for 10%)",
				" correlation 12 (real value: 0.1 for 10%)",
				" correlation 13 (real value: 0.1 for 10%)",
				" correlation 23 (real value: 0.1 for 10%)",
				" mu1 : drift de l'index 1",
				" mu2 : drift de l'index 2",
				" mu3 : drift de l'index 3",
				" coefficient 0 : (moins strike): 0.05 for 5%",
				" coefficient 1 : ",
				" coefficient 2 : ",
				" coefficient 3 : ",
				" year fraction (1 for one year)",
				" number of points in the numerical integration (default=64)"
        },
		{
        		" Local_LogNormal_TriSpreadOption_Der2",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRRR",					/// 21 parametres = 20 parametre d'entree + 1 parametre de retour 
                " ARM_CF_LN_TriSpreadOptDer2",
				" derivative index1,derivative index2, fwd1 value,fwd2 value,fwd3 value,sigma1,sigma2,sigma3,correlation12,correlation13,correlation23,mu1,mu2,mu3,a0,a1,a2,a3,maturity,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives a second derivative of the FORWARD VALUE of a tri index  spread option (3Lognormals)",
				" index of the dimension: 0=fwd1, 1=fwd2, 2=fwd3, 3=sig1, 4=sig2, 5=sig3, 6=rho12, 7=rho13,8=rho23,9=mu1,10=mu2,11=mu3,12=a0,13=a1,14=a2,15=a3,16=t",
				" index of the dimension 2: 0=fwd1, 1=fwd2, 2=fwd3, 3=sig1, 4=sig2, 5=sig3, 6=rho12, 7=rho13,8=rho23,9=mu1,10=mu2,11=mu3,12=a0,13=a1,14=a2,15=a3,16=t",
				" fwd value 1 (real value: 0.05 for 5%)",
				" fwd value 2 (real value: 0.05 for 5%)",
				" fwd value 3 (real value: 0.05 for 5%)",
				" LN vol 1 (real value: 0.1 for 10%)",
				" LN vol 2 (real value: 0.1 for 10%)",
				" LN vol 3 (real value: 0.1 for 10%)",
				" correlation 12 (real value: 0.1 for 10%)",
				" correlation 13 (real value: 0.1 for 10%)",
				" correlation 23 (real value: 0.1 for 10%)",
				" mu1 : drift de l'index 1",
				" mu2 : drift de l'index 2",
				" mu3 : drift de l'index 3",
				" coefficient 0 : (moins strike): 0.05 for 5%",
				" coefficient 1 ",
				" coefficient 2 ",
				" coefficient 3 ",
				" year fraction (1 for one year)",
				" number of points in the numerical integration (default=64)"
        },
		{
        		" Local_LogNormal_TriSpreadDigitalOption",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRR",					/// 20 parametres = 19 d'entree + 1 parametre de retour 
                " ARM_CF_LN_TriSpreadDigiOpt",
				" fwd1,fwd2,fwd3,sigma1,sigma2,sigma3,correlation12,correlation13,correlation23,mu1,mu2,mu3,a0,a1,a2,a3,maturity,CallPut,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives a  FORWARD VALUE of a tri index digital spread option (3Lognormals)",
				" fwd value 1 (real value: 0.05 for 5%)",
				" fwd value 2 (real value: 0.05 for 5%)",
				" fwd value 3 (real value: 0.05 for 5%)",
				" LN vol 1 (real value: 0.1 for 10%)",
				" LN vol 2 (real value: 0.1 for 10%)",
				" LN vol 3 (real value: 0.1 for 10%)",
				" correlation 12 (real value: 0.1 for 10%)",
				" correlation 13 (real value: 0.1 for 10%)",
				" correlation 23 (real value: 0.1 for 10%)",
				" mu1 : drift de l'index 1",
				" mu2 : drift de l'index 2",
				" mu3 : drift de l'index 3",
				" coefficient 0 : (moins strike): 0.05 for 5%",
				" coefficient 1 : ",
				" coefficient 2 : ",
				" coefficient 3 : ",
				" year fraction (1 for one year)",
				" Call or Put",
				" number of points in the numerical integration (default=64)"
        },
		{
        		" Local_LogNormal_TriSpreadDigitalOption_Der",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRRR",					/// 21 parametres = 20 d'entree + 1 parametre de retour 
                " ARM_CF_LN_TriSpreadDigiOptDer",
				" Index Derivative,fwd1,fwd2,fwd3,sigma1,sigma2,sigma3,correlation12,correlation13,correlation23,mu1,mu2,mu3,a0,a1,a2,a3,maturity,CallPut,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives a  FORWARD VALUE of a tri index digital spread option (3Lognormals)",
				" index of the dimension: 0=fwd1, 1=fwd2, 2=fwd3, 3=sig1, 4=sig2, 5=sig3, 6=rho12, 7=rho13,8=rho23,9=mu1,10=mu2,11=mu3,12=a0,13=a1,14=a2,15=a3,16=t",
				" fwd value 1 (real value: 0.05 for 5%)",
				" fwd value 2 (real value: 0.05 for 5%)",
				" fwd value 3 (real value: 0.05 for 5%)",
				" LN vol 1 (real value: 0.1 for 10%)",
				" LN vol 2 (real value: 0.1 for 10%)",
				" LN vol 3 (real value: 0.1 for 10%)",
				" correlation 12 (real value: 0.1 for 10%)",
				" correlation 13 (real value: 0.1 for 10%)",
				" correlation 23 (real value: 0.1 for 10%)",
				" mu1 : drift de l'index 1",
				" mu2 : drift de l'index 2",
				" mu3 : drift de l'index 3",
				" coefficient 0 : (moins strike): 0.05 for 5%",
				" coefficient 1 : ",
				" coefficient 2 : ",
				" coefficient 3 : ",
				" year fraction (1 for one year)",
				" Call or Put",
				" number of points in the numerical integration (default=64)"
        },
		{
        		" Local_LogNormal_TriSpreadDigitalOption_Der2",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRRRR",					/// 22 parametres = 21 d'entree + 1 parametre de retour 
                " ARM_CF_LN_TriSpreadDigiOptDer2",
				" Index Derivative,Index Derivative 2,fwd1,fwd2,fwd3,sigma1,sigma2,sigma3,correlation12,correlation13,correlation23,mu1,mu2,mu3,a0,a1,a2,a3,maturity,CallPut,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives a  FORWARD VALUE of a tri index digital spread option (3Lognormals)",
				" index of the dimension: 0=fwd1, 1=fwd2, 2=fwd3, 3=sig1, 4=sig2, 5=sig3, 6=rho12, 7=rho13,8=rho23,9=mu1,10=mu2,11=mu3,12=a0,13=a1,14=a2,15=a3,16=t",
				" index of the dimension 2: 0=fwd1, 1=fwd2, 2=fwd3, 3=sig1, 4=sig2, 5=sig3, 6=rho12, 7=rho13,8=rho23,9=mu1,10=mu2,11=mu3,12=a0,13=a1,14=a2,15=a3,16=t",
				" fwd value 1 (real value: 0.05 for 5%)",
				" fwd value 2 (real value: 0.05 for 5%)",
				" fwd value 3 (real value: 0.05 for 5%)",
				" LN vol 1 (real value: 0.1 for 10%)",
				" LN vol 2 (real value: 0.1 for 10%)",
				" LN vol 3 (real value: 0.1 for 10%)",
				" correlation 12 (real value: 0.1 for 10%)",
				" correlation 13 (real value: 0.1 for 10%)",
				" correlation 23 (real value: 0.1 for 10%)",
				" mu1 : drift de l'index 1",
				" mu2 : drift de l'index 2",
				" mu3 : drift de l'index 3",
				" coefficient 0 : (moins strike): 0.05 for 5%",
				" coefficient 1 : ",
				" coefficient 2 : ",
				" coefficient 3 : ",
				" year fraction (1 for one year)",
				" Call or Put",

        },
		{
        		" Local_LogNormal_TriSpreadDigitalOption2",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRR",					/// 20 parametres = 19 d'entree + 1 parametre de retour 
                " ARM_CF_LN_TriSpreadDigiOpt2",
				" fwd1,fwd2,fwd3,sigma1,sigma2,sigma3,correlation12,correlation13,correlation23,mu1,mu2,mu3,a0,a2,a3,b0,b1,maturity,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives a  FORWARD VALUE of a digital spread option that pays 1 if (b0+b1*S1)> and a0+a2*S2+a3*S3>0",
				" fwd value 1 (real value: 0.05 for 5%)",
				" fwd value 2 (real value: 0.05 for 5%)",
				" fwd value 3 (real value: 0.05 for 5%)",
				" LN vol 1 (real value: 0.1 for 10%)",
				" LN vol 2 (real value: 0.1 for 10%)",
				" LN vol 3 (real value: 0.1 for 10%)",
				" correlation 12 (real value: 0.1 for 10%)",
				" correlation 13 (real value: 0.1 for 10%)",
				" correlation 23 (real value: 0.1 for 10%)",
				" mu1 : drift de l'index 1",
				" mu2 : drift de l'index 2",
				" mu3 : drift de l'index 3",
				" coefficient a0 : (moins strike): 0.05 for 5%",
				" coefficient a2 : ",
				" coefficient a3 : ",
				" coefficient b0 : ",
				" coefficient b1 : ",
				" year fraction (1 for one year)",
				" number of points in the numerical integration (default=64)"
        },
		{
        		" Local_LogNormal_TriSpreadDigitalOption2_Der",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRRR",					/// 21 parametres = 20 d'entree + 1 parametre de retour 
                " ARM_CF_LN_TriSpreadDigiOpt2Der",
				" Index Derivative,fwd1,fwd2,fwd3,sigma1,sigma2,sigma3,correlation12,correlation13,correlation23,mu1,mu2,mu3,a0,a2,a3,b0,b1,maturity,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives Derivative of a  FORWARD VALUE of a digital spread option that pays (b0+b1*S1)+ if a0+a2*S2+a3*S3>0",
				" index of the dimension: 0=fwd1, 1=fwd2, 2=fwd3, 3=sig1, 4=sig2, 5=sig3, 6=rho12, 7=rho13,8=rho23,9=mu1,10=mu2,11=mu3,12=a0,13=a2,14=a3,15=b0,16=b1,17=t",
				" fwd value 1 (real value: 0.05 for 5%)",
				" fwd value 2 (real value: 0.05 for 5%)",
				" fwd value 3 (real value: 0.05 for 5%)",
				" LN vol 1 (real value: 0.1 for 10%)",
				" LN vol 2 (real value: 0.1 for 10%)",
				" LN vol 3 (real value: 0.1 for 10%)",
				" correlation 12 (real value: 0.1 for 10%)",
				" correlation 13 (real value: 0.1 for 10%)",
				" correlation 23 (real value: 0.1 for 10%)",
				" mu1 : drift de l'index 1",
				" mu2 : drift de l'index 2",
				" mu3 : drift de l'index 3",
				" coefficient a0 : (moins strike): 0.05 for 5%",
				" coefficient a2 : ",
				" coefficient a3 : ",
				" coefficient b0 : ",
				" coefficient b1 : ",
				" year fraction (1 for one year)",
				" number of points in the numerical integration (default=64)"
        },
		{
        		" Local_LogNormal_TriSpreadDigitalOption2_Der2",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRRRR",					/// 22 parametres = 21 d'entree + 1 parametre de retour 
                " ARM_CF_LN_TriSpreadDigiOpt2Der2",
				" Index Derivative,Index Derivative 2,fwd1,fwd2,fwd3,sigma1,sigma2,sigma3,correlation12,correlation13,correlation23,mu1,mu2,mu3,a0,a2,a3,b0,b1,maturity,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives second derivatives of a FORWARD VALUE of a digital spread option that pays (b0+b1*S1)+ if a0+a2*S2+a3*S3>0",
				" index of the dimension: 0=fwd1, 1=fwd2, 2=fwd3, 3=sig1, 4=sig2, 5=sig3, 6=rho12, 7=rho13,8=rho23,9=mu1,10=mu2,11=mu3,12=a0,13=a2,14=a3,15=b0,16=b1,17=t",
				" index of the dimension: 0=fwd1, 1=fwd2, 2=fwd3, 3=sig1, 4=sig2, 5=sig3, 6=rho12, 7=rho13,8=rho23,9=mu1,10=mu2,11=mu3,12=a0,13=a2,14=a3,15=b0,16=b1,17=t",
				" fwd value 1 (real value: 0.05 for 5%)",
				" fwd value 2 (real value: 0.05 for 5%)",
				" fwd value 3 (real value: 0.05 for 5%)",
				" LN vol 1 (real value: 0.1 for 10%)",
				" LN vol 2 (real value: 0.1 for 10%)",
				" LN vol 3 (real value: 0.1 for 10%)",
				" correlation 12 (real value: 0.1 for 10%)",
				" correlation 13 (real value: 0.1 for 10%)",
				" correlation 23 (real value: 0.1 for 10%)",
				" mu1 : drift de l'index 1",
				" mu2 : drift de l'index 2",
				" mu3 : drift de l'index 3",
				" coefficient a0 : (moins strike): 0.05 for 5%",
				" coefficient a2 : ",
				" coefficient a3 : ",
				" coefficient b0 : ",
				" coefficient b1 : ",
				" year fraction (1 for one year)",

        },
		{
        		" Local_Lognormal_Asian_VanillaOption",		/// name of the C++ function
                " RRRRRRRRR",					/// 9 parametres = 8 d'entree + 1 parametre de retour 
                " ARM_CF_LN_Asian_VanillaOption",
				" Spot, strike, maturity, interest_rate, volatility, xabs, CallorPut, nb ",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " Spot (real value: 0.05 for 5%)",
				" Strike (real value: 0.05 for 5%)",
				" year fraction (1 for one year)",
				" Interest Rate (real value: 0.05 for 5%)",
				" Vol  (real value: 0.1 for 10%)",
				" alpha for Integration abcisse  (real value: 0.1 for 10%)",
				" Call or Put",
				" nb for the integration",
        },
		{
        		" Local_Lognormal_Asian_VanillaOption_Der",		/// name of the C++ function
                " RRRRRRRRRR",					/// 10 parametres = 9 d'entree + 1 parametre de retour 
                " ARM_CF_LN_Asian_VanillaOptionDer",
				" Index, Spot, strike, maturity, interest_rate, volatility, xabs, CallorPut, nb ",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
				" index of the dimension: 0=spot,1=strike,2=T,3=IR, 4=vol", 
                " Spot (real value: 0.05 for 5%)",
				" Strike (real value: 0.05 for 5%)",
				" year fraction (1 for one year)",
				" Interest Rate (real value: 0.05 for 5%)",
				" Vol  (real value: 0.1 for 10%)",
				" alpha for Integration abcisse  (real value: 0.1 for 10%)",
				" Call or Put",
				" nb for the integration",
        },
		{
        		" Local_Lognormal_Asian_VanillaOption_Der2",		/// name of the C++ function
                " RRRRRRRRRRR",					/// 11 parametres = 10 d'entree + 1 parametre de retour 
                " ARM_CF_LN_Asian_VanillaOptionDer2",
				" Index, Index2,Spot, strike, maturity, interest_rate, volatility, xabs, CallorPut, nb ",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
				" index of the dimension: 0=spot,1=strike,2=T,3=IR, 4=vol", 
				" index of the dimension2: 0=spot,1=strike,2=T,3=IR, 4=vol", 
                " Spot (real value: 0.05 for 5%)",
				" Strike (real value: 0.05 for 5%)",
				" year fraction (1 for one year)",
				" Interest Rate (real value: 0.05 for 5%)",
				" Vol  (real value: 0.1 for 10%)",
				" alpha for Integration abcisse  (real value: 0.1 for 10%)",
				" Call or Put",
				" nb for the integration",
        },
		{
        		" Local_GaussianIntegrals_Legendre_Coeffs",		/// name of the C++ function
                " RRRR",					/// 4 parametres = 3 d'entree + 1 parametre de retour 
                " ARM_CF_Util_GaussianIntegrals_Legendre_Coeffs",
				" a,b,nb",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
				" return the roots and weigthts for an integral between a and b",
				 
				" start of the integral",
				" end of the integral",
				" np points"
        },
		{
        		" Local_GaussianIntegrals_Hermite_Coeffs",		/// name of the C++ function
                " RR",											/// 2 parametres = 1 d'entree + 1 parametre de retour 
                " ARM_CF_Util_GaussianIntegrals_Hermite_Coeffs",
				" nb ",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
				" np points", 
        },
		{
        		" Local_GaussianIntegrals_Laguerre_Coeffs",		/// name of the C++ function
                " RRR",											/// 3 parametres = 2 d'entree + 1 parametre de retour 
                " ARM_CF_Util_GaussianIntegrals_Laguerre_Coeffs",
				" alpha,nb ",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
				" computes Gaussian Laguerre integral coefficients",
				" alpha",
				" np points", 
        },
		{
        		" Local_ConvertFXOptionToStrike",		/// name of the C++ function
                " RRRRRR",								/// 6 parametres = 5 d'entree + 1 parametre de retour 
                " ARM_CF_ConvertFXOptionToStrike",
				" TargetDelta, Fwd, TotalVol, CallPut,initValue",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
				" Converts fx option to strike",
				" double (TargetDelta)",
				" double (Fwd)",
				" double (TotalVol)",
				" double (CallPut)",
				" double (initValue)",
        },
		{
        		" Local_GHeston_VanillaOption_ModelVector",		/// name of the C++ function
                " RRRRRRRRRRRRRRRR" ,					/// 14 parametres = 12 d'entree + 1 parametre de retour 
                " ARM_CF_GHeston_Option_ModelVector",
				" strike ,timetomaturity,callput,interpolation,maturities,forward,Vol initiale ,long term limit,mean reveting speed,vol of vol,correlation,jump probability, jump size,jump volatility,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the Generalized Heston value of the forward value of a vanilla option",
				" strike (real value: 0.05 for 5%)",
				" time to maturity (years)",
				" CALL or PUT",
				" Interpolation : LINEAR or PARABOLIC",
				" maturity vector",
				" fwd value (real value: 0.05 for 5%) vector",
				" vol initiale (real value: 0.05 for 5%) vector",
				" long terme limit (0.15 for 15% vol) vector",
				" mean reverting speed vector",
				" volatility of volatility vector",
				" correlation between the volatility and the underlying (-0.5 for -50%) vector",
				" jump probability (0.15 for 15%) vector",
				" jump size (expectation of the return) vector",
				" jump volatility vector",
				" number of Legendre points for the integration "
        },
		{
        		" Local_GHeston_Implicit_Volatility_ModelVector",		/// name of the C++ function
                " RRRRRRRRRRRRRRR" ,					/// 14 parametres = 12 d'entree + 1 parametre de retour 
                " ARM_CF_GHeston_ImplicitVol_ModelVector",
				" strike ,timetomaturity,interpolation,maturities,forward,Vol initiale ,long term limit,mean reveting speed,vol of vol,correlation,jump probability, jump size,jump volatility,[n]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the Generalized Heston value of the forward value of a vanilla option",
				" strike (real value: 0.05 for 5%)",
				" time to maturity (years)",
				" Interpolation : LINEAR or PARABOLIC",
				" maturity vector",
				" fwd value (real value: 0.05 for 5%) vector",
				" vol initiale (real value: 0.05 for 5%) vector",
				" long terme limit (0.15 for 15% vol) vector",
				" mean reverting speed vector",
				" volatility of volatility vector",
				" correlation between the volatility and the underlying (-0.5 for -50%) vector",
				" jump probability (0.15 for 15%) vector",
				" jump size (expectation of the return) vector",
				" jump volatility vector",
				" number of Legendre points for the integration "
        },
		{
        		" Local_SABR_Model_Calibrate",		/// name of the C++ function
                " RRRRRRRRRRRRR" ,					/// 12 parametres = 11 d'entree + 1 parametre de retour 
                " ARM_CF_SABR_Smile_Calibrate",
				" strikes, volatilities , forward, maturity, SABR Type, nbsteps, algorithm Type ,alpha, beta, rho ,nu,weigths",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
				" ",
                " gives the calibrated values of SABR Parameters",
				" strikes ector (real value: 0.05 for 5%)",
				" implied volatilities vector ( 0.5 for 50%)",
				" Underlying value (0.05 for 5%)",
				" maturity (in years)",
				" Type of SABR formula: SABR_IMPLNVOL(2),SABR_A,SABR_G,DIRECTEXACT,DIRECTGEOMETRIC,DIRECTARITHMETIC,NORMALEXACT,NORMALGEOMETRIC,NORMALARITHMETIC,ANALYTICZP0,ANALYTICZP2",	
				" nb of steps in gaussian integration of analytical formula (<=175) ",
				" Algorithm Type: NLIN_LSQ,LSQ_DERIV,..",
				" inital value for alpha  >0" ,
				" initial value for beta (between 0.5 and 0.9999",
				" initial value for rho between -0.9999 and +0.9999",
				" initial value for nu >0",
				" vector of weigths"
		},
		{
        		" Local_SABR_Model_Calibrate_WithForward",		/// name of the C++ function
                " RRRRRRRRRRRRR" ,					/// 12 parametres = 11 d'entree + 1 parametre de retour 
                " ARM_CF_SABR_Smile_Calibrate_WithForward",
				" strike vector, price vector,maturity,type of SABR,nbsteps,algorithm,alpha0,beta0,rho0,nu0,f0,weigths",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the calibrated values of SABR Parameters",
				" strike vectors (real value: 0.05 for 5%)",
				" price vector",
				" maturity (in years)",
				" Type of SABR formula : SABR_IMPLNVOL(2),SABR_A,SABR_G,DIRECTEXACT,DIRECTGEOMETRIC,DIRECTARITHMETIC,NORMALEXACT,NORMALGEOMETRIC,NORMALARITHMETIC,ANALYTICZP0,ANALYTICZP2",	
				" nb of steps in gaussian integration of analytical formula (<=175) ",
				" Algorithm: Values= NLIN_LSQ,LSQ_DERIV,..",
				" inital value for alpha  >0" ,
				" initial beta (between 0.5 and 0.9999",
				" initial value for rho between -0.9999 and +0.9999",
				" initial value for nu >0",
				" forward price (0.05 for 5%)",
				" vector of weigths"
		},
		{
        		" Local_SABR_Model_Calibrate_BetaFixed",		/// name of the C++ function
                " RRRRRRRRRRRRRRRR" ,					/// 16 parametres = 15 d'entree + 1 parametre de retour 
                " ARM_CF_SABR_Smile_Calibrate_BetaFixed",
				" strikes, impliedVol, forward, beta, maturity, SABR Type, nbsteps, algorithm Type ,alpha, rho ,nu,weigths,a_exp,a_tanh,kt_tanh",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the calibrated values of SABR Parameters",
				" strike vectors (real value: 0.05 for 5%)",
				" price vector",
				" maturity (in years)",
				" f forward price ",
				" beta (between 0.5 and 0.9999",
				" Type of SABR formula : SABR_IMPLNVOL(2),SABR_A,SABR_G,DIRECTEXACT,DIRECTGEOMETRIC,DIRECTARITHMETIC,NORMALEXACT,NORMALGEOMETRIC,NORMALARITHMETIC,ANALYTICZP0,ANALYTICZP2",	
				" nb of steps in gaussian integration of analytical formula (<=175) ",
				" Algorithm: Values= NLIN_LSQ,LSQ_DERIV,..",
				" inital value for alpha  >0" ,
				" initial value for rho between -0.9999 and +0.9999",
				" initial value for nu >0",
				" vector of weigths",
				" strike cuter coef exp (def=0.01)",
				" strike cuter coef tanh (def=1.5)",
				" strike cuter K cutoff tanh (def=0.3)"
		},
		{
        		" Local_SABR_Model_Calibrate_BetaFixed_Linked",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRR" ,					/// 18 parametres = 17 d'entree + 1 parametre de retour 
                " ARM_CF_SABR_Smile_Calibrate_BetaFixed_linked",
				" strike vector, price vector,forward,beta,maturity,type of SABR,nbsteps,algorithm,alpha0,rho0,nu0,alphap,rhop,nup,rweight_alpha,rweight_rho,rweight_nu,weigths",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the calibrated values of SABR Parameters",
				" strike vectors (real value: 0.05 for 5%)",
				" price vector",
				" maturity (in years)",
				" f forward price ",
				" beta (between 0.5 and 0.9999",
				" Type of SABR formula : SABR_IMPLNVOL(2),SABR_A,SABR_G,DIRECTEXACT,NORMALEXACT,NORMALGEOMETRIC,NORMALARITHMETIC,ANALYTICZP0"	
				" nb of steps in gaussian integration of analytical formula (<=175) ",
				" Algorithm: Values= NLIN_LSQ,LSQ_DERIV,..",
				" inital value for alpha  >0" ,
				" initial value for rho between -0.9999 and +0.9999",
				" initial value for nu >0",
				" alphap: previous alpha (linked)",
				" rhop: previous rho (linked)",
				" nup : previous nu (linked)",
				" rweight_alpha : weight for the alpha linked part",
				" rweight_rho : weight for the rho linked part",
				" rweight_nu : weight for the nu linked part",
				" vector of weigths"
		},
		{
				" Local_SABR_Calibrate_BetaFixedToOne",
				" RRRRRRR",
				" ARM_CF_SABR_Calibrate_BetaFixedToOne",
				" strike vector, imp vol vector, forward, maturity, [atmvol], [weight]",
				" 1",
				XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
				" smile calibration with sabr beta = 1",
				" ",
				" ",
				" ",
				" ",
				" if force atm vol calibration",
				" vector of weights"
		},
		{
        		" Local_GHeston_Model_Calibrate_Total",		/// name of the C++ function
                " RRRRRRRRRRRRRRR" ,					/// 15 parametres = 14 d'entree + 1 parametre de retour 
                " ARM_CF_GHeston_Smile_Calibrate_Total",
				" strike vector, price vector,forward,maturity,nbsteps,algorithm,V0,omega, theta, ksi,rho,muJ,sigmaJ,lambda",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the calibrated values of SABR Parameters",
				" strike vectors (real value: 0.05 for 5%)",
				" price vector",
				" f forward price ",
				" maturity (in years)",
				" nb of steps in gaussian integration of analytical formula (<=175) ",
				" Algorithm: Values= NLIN_LSQ,LSQ_DERIV,..",
				" inital value for V0  >0" ,
				" inital value for omega  >0" ,
				" inital value for theta  >0" ,
				" inital value for ksi  >0" ,
				" initial value for rho between -0.9999 and +0.9999",
				" inital value for muJ  " ,
				" inital value for sigmaJ  >0" ,
				" initial value for lambda >0"
		},
		{
        		" Local_GHeston_Model_Calibrate_NoJump",		/// name of the C++ function
                " RRRRRRRRRRRRRRR" ,					/// 15 parametres = 14 d'entree + 1 parametre de retour 
                " ARM_CF_GHeston_Smile_Calibrate_NoJump",
				" strike vector, price vector,forward,maturity,muJ,sigmaJ,lambda,nbsteps,algorithm,V0,omega, theta, ksi,rho",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the calibrated values of SABR Parameters",
				" strike vectors (real value: 0.05 for 5%)",
				" price vector",
				" f forward price ",
				" maturity (in years)",
				" muJ ",
				" sigmaJ >0",
				" lambda >0",
				" nb of steps in gaussian integration of analytical formula (<=175) ",
				" Algorithm: Values= NLIN_LSQ,LSQ_DERIV,..",
				" inital value for V0  >0" ,
				" inital value for omega  >0" ,
				" inital value for theta  >0" ,
				" inital value for ksi  >0" ,
				" initial value for rho between -0.9999 and +0.9999"	
		},

		{
        		" Local_Student_SABR_Power_SpreadOption",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRRRR",					/// 22 parametres = 21 d'entree + 1 parametre de retour 
                " ARM_CF_SABRStudent_SpreadOpt",
				" index1,index2 ,alpha1,beta1,rho1,nu1,alpha2,beta2,rho2,nu2,copulaStruct,maturity,flag,a1,b1,k1,a2,b2,k2,Params",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the  FORWARD VALUE of a power spread option (SABR+Gaussian)",
				" fwd value (real value: 0.05 for 5%)",
				" fwd value (real value: 0.05 for 5%)",
				" alpha index1",
				" beta index1",
				" rho index1",
				" nu index1",
				" alpha index2",
				" beta index2",
				" rho index2",
				" nu index2",
				" correlation, student degre",
				" year fraction (1 for one year)",
				" Type of SABR formulas : SABR_IMPLNVOL(2),SABR_A,SABR_G,DIRECTEXACT,DIRECTGEOMETRIC,DIRECTARITHMETIC,NORMALEXACT,NORMALGEOMETRIC,NORMALARITHMETIC,ANALYTICZP0,ANALYTICZP2",
				" a1: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" b1: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" k1: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" a2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" b2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" k2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" Nb(=120),alpha_exp(=0.1),alpha_tanh(=1.5),kb_tanh=(=0.02)"
        },
		{
        		" Local_Mepi_VanillaOption_STOBS",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRRR",					/// 20 parametres = 19 d'entree + 1 parametre de retour 
                " ARM_CF_Mepi_VanillaOption_STOBS",
				" Today,P0,K,T,zero_curve,b,YearlyFees,cashspread,SABRModel,minExp,maxExp,riskFac,InitProtect,FinalProtect,AveragingPeriodNb,Asianreset,CallOrPut,TotalPeriodNb,nz,nbh",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the  present VALUE of a Mepi Option in the stochastic volatility Black and Sholes Setting",
				" Date of today",
				" Present Value of the initial Portfolio",
				" Strike",
				" Maturity in Years",
				" zero curve",
				" borrowing spread",
				" YearlyFees ",
				" CashSpread ",
				" Model SABR",
				" Minimal Exposition (0.2 for 20 %)",
				" Maximal Exposition (0.2 for 20 %)",
				" Risk Factor",
				" Initial Garantee, in case of linear increasing garanty",
				" Garantied fraction at maturity (1 means 100%)",
				" Averaging Periods Number ",
				" Already Computed Average",
				" Call Or Put",
				" Number of periods",
				" size of the pathstate container",
				" nb of Legendre points for the volatility Integration",

        },
		{
        		" Local_Mepi_VanillaOption_STOBS_delta",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRRR",					/// 20 parametres = 19 d'entree + 1 parametre de retour 
                " ARM_CF_Mepi_VanillaOption_STOBS_delta",
				" Today,P0,K,T,zero_curve,b,YearlyFees,cashspread,SABRModel,minExp,maxExp,riskFac,InitProtect,FinalProtect,AveragingPeriodNb,Asianreset,CallOrPut,TotalPeriodNb,nz,nbh",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the  present VALUE of a Mepi Option in the stochastic volatility Black and Sholes Setting",
				" Date of today",
				" Present Value of the initial Portfolio",
				" Strike",
				" Maturity in Years",
				" zero curve",
				" borrowing spread",
				" YearlyFees ",
				" CashSpread ",
				" Model SABR",
				" Minimal Exposition (0.2 for 20 %)",
				" Maximal Exposition (0.2 for 20 %)",
				" Risk Factor",
				" Initial Garantee, in case of linear increasing garanty",
				" Garantied fraction at maturity (1 means 100%)",
				" Averaging Periods Number ",
				" Already Computed Average",
				" Call Or Put",
				" Number of periods",
				" size of the pathstate container",
				" nb of Legendre points for the volatility Integration",

        },
		{
        		" Local_Mepi_VanillaOption_SABR",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRRR",					/// 23 parametres = 22 d'entree + 1 parametre de retour 
                " ARM_CF_Mepi_VanillaOption_SABR",
				" Today,P0,K,Maturity,zero_curve,b,YearlyFees,cashspread,SABRModel,spot_name,minExp,maxExp,riskFac,InitProtect,FinalProtect,AveragingPeriodNb,Asianreset,CallOrPut,ResetFreq,NbMC2",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the  present VALUE of a Mepi Option in the stochastic volatility Black and Sholes Setting",
				" Date of today",
				" Present Value of the initial Portfolio",
				" Strike",
				" Maturity Date",
				" zero curve",
				" borrowing spread",
				" YearlyFees ",
				" CashSpread ",
				" Model SABR",
				" spot name",
				" Minimal Exposition (0.2 for 20 %)",
				" Maximal Exposition (0.2 for 20 %)",
				" Risk Factor",
				" Initial Garantee, in case of linear increasing garanty",
				" Garantied fraction at maturity (1 means 100%)",
				" Averaging Periods Number ",
				" Already Computed Average",
				" Call Or Put",
				" Reset Frequency (12 for monthly)",
				" Number of Path MC2",
        },
		{
        		" Local_Mepi_VanillaOption_SABR_delta",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRRR",					/// 22 parametres = 21 d'entree + 1 parametre de retour 
                " ARM_CF_Mepi_VanillaOption_SABR_delta",
				" Today,P0,K,Maturity,zero_curve,b,YearlyFees,cashspread,SABRModel,spotname,minExp,maxExp,riskFac,InitProtect,FinalProtect,AveragingPeriodNb,Asianreset,CallOrPut,ResetFreq,NbMC2",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the  present VALUE of a Mepi Option in the stochastic volatility Black and Sholes Setting",
				" Date of today",
				" Present Value of the initial Portfolio",
				" Strike",
				" Maturity Date",
				" zero curve",
				" borrowing spread",
				" YearlyFees ",
				" CashSpread ",
				" Model SABR",
				" Spot Name",
				" Minimal Exposition (0.2 for 20 %)",
				" Maximal Exposition (0.2 for 20 %)",
				" Risk Factor",
				" Initial Garantee, in case of linear increasing garanty",
				" Garantied fraction at maturity (1 means 100%)",
				" Averaging Periods Number ",
				" Already Computed Average",
				" Call Or Put",
				" Reset Frequency (12 for monthly)",
				" Number of Path MC2",

        },
				{
        		" Local_StochasticVol_LN_VanillaOption",		/// name of the C++ function
                " RRRRRRRRRRRR" ,					/// 12 parametres = 11 d'entree + 1 parametre de retour 
                " ARM_CF_StochasticVol_LN_Geo_VanillaOption",
				" spot,strike,time_to_maturity,r,volatility,VolDrift,VolVol,Averaging,Reset,callput,[nbterms]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the spot value of a european vanilla option with a BS model with stochastic volatility",
				" spot value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" time to reset of the option",
				" interest rate (0.05 for 5% vol)",
				" volatility (0.15 for 15% vol)",
				" VolDrift: drift of the volatility",
				" VolVol: Volatility of the volatility",
				" Averaging preriod of time",
				" Reset Value for the average of Nav up to now",
				" CALL or PUT",
				" nb of Legendre terms (typically 20)"
        },
			{
        		" Local_StochasticVol_LN_VanillaOption_Der",		/// name of the C++ function
                " RRRRRRRRRRRRR" ,					/// 13 parametres = 12 d'entree + 1 parametre de retour 
                " ARM_CF_StochasticVol_LN_Geo_VanillaOptionDer",
				" index,spot,strike,time_to_maturity,r,volatility,VolDrift,VolVol,Averaging,Reset,callput,[nbterms]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the spot value of a european vanilla option with a BS model with stochastic volatility",
				" index : 0=spot,1=strike,2=maturity,3=r,4=vol,5=VolDrift,6=VolVol,7=Averaging,8=Reset",
				" spot value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" time to reset of the option",
				" interest rate (0.05 for 5% vol)",
				" volatility (0.15 for 15% vol)",
				" VolDrift: drift of the volatility",
				" VolVol: Volatility of the volatility",
				" Averaging preriod of time",
				" Reset Value for the average of Nav up to now",
				" CALL or PUT",
				" nb of Legendre terms (typically 20)"
        },
		{
        		" Local_StochasticVol_LN_VanillaOption_Der2",		/// name of the C++ function
                " RRRRRRRRRRRRRR" ,					/// 14 parametres = 13 d'entree + 1 parametre de retour 
                " ARM_CF_StochasticVol_LN_Geo_VanillaOptionDer2",
				" index,index2,spot,strike,time_to_maturity,r,volatility,VolDrift,VolVol,Averaging,Reset,callput,[nbterms]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the spot value of a european vanilla option with a BS model with stochastic volatility",
				" index : 0=spot,1=strike,2=maturity,3=r,4=vol,5=VolDrift,6=VolVol,7=Averaging,8=Reset",
				" index : 0=spot,1=strike,2=maturity,3=r,4=vol,5=VolDrift,6=VolVol,7=Averaging,8=Reset",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" time to reset of the option",
				" interest rate (0.05 for 5% vol)",
				" volatility (0.15 for 15% vol)",
				" VolDrift: drift of the volatility",
				" VolVol: Volatility of the volatility",
				" Averaging preriod of time",
				" Reset Value for the average of Nav up to now",
				" CALL or PUT",
				" nb of Legendre terms (typically 20)"
        },
		{
        		" Local_StochasticVol_LN_Ari_VanillaOption",		/// name of the C++ function
                " RRRRRRRRRRRR" ,					/// 12 parametres = 11 d'entree + 1 parametre de retour 
                " ARM_CF_StochasticVol_LN_Ari_VanillaOption",
				" spot,strike,time_to_maturity,r,volatility,VolDrift,VolVol,Averaging,Reset,callput,[nbterms]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the spot value of a european vanilla option with a BS model with stochastic volatility",
				" spot value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" time to reset of the option",
				" interest rate (0.05 for 5% vol)",
				" volatility (0.15 for 15% vol)",
				" VolDrift: drift of the volatility",
				" VolVol: Volatility of the volatility",
				" Averaging preriod of time",
				" Reset Value for the average of Nav up to now",
				" CALL or PUT",
				" nb of Legendre terms (typically 20)"
        },
			{
        		" Local_StochasticVol_LN_Ari_VanillaOption_Der",		/// name of the C++ function
                " RRRRRRRRRRRRR" ,					/// 13 parametres = 12 d'entree + 1 parametre de retour 
                " ARM_CF_StochasticVol_LN_Ari_VanillaOptionDer",
				" index,spot,strike,time_to_maturity,r,volatility,VolDrift,VolVol,Averaging,Reset,callput,[nbterms]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the spot value of a european vanilla option with a BS model with stochastic volatility",
				" index : 0=spot,1=strike,2=maturity,3=r,4=vol,5=VolDrift,6=VolVol,7=Averaging,8=Reset",
				" spot value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" time to reset of the option",
				" interest rate (0.05 for 5% vol)",
				" volatility (0.15 for 15% vol)",
				" VolDrift: drift of the volatility",
				" VolVol: Volatility of the volatility",
				" Averaging preriod of time",
				" Reset Value for the average of Nav up to now",
				" CALL or PUT",
				" nb of Legendre terms (typically 20)"
        },
		{
        		" Local_StochasticVol_LN_Ari_VanillaOption_Der2",		/// name of the C++ function
                " RRRRRRRRRRRRRR" ,					/// 14 parametres = 13 d'entree + 1 parametre de retour 
                " ARM_CF_StochasticVol_LN_Ari_VanillaOptionDer2",
				" index,index2,spot,strike,time_to_maturity,r,volatility,VolDrift,VolVol,Averaging,Reset,callput,[nbterms]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the spot value of a european vanilla option with a BS model with stochastic volatility",
				" index : 0=spot,1=strike,2=maturity,3=r,4=vol,5=VolDrift,6=VolVol,7=Averaging,8=Reset",
				" index : 0=spot,1=strike,2=maturity,3=r,4=vol,5=VolDrift,6=VolVol,7=Averaging,8=Reset",
				" spot value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" time to reset of the option",
				" interest rate (0.05 for 5% vol)",
				" volatility (0.15 for 15% vol)",
				" VolDrift: drift of the volatility",
				" VolVol: Volatility of the volatility",
				" Averaging preriod of time",
				" Reset Value for the average of Nav up to now",
				" CALL or PUT",
				" nb of Legendre terms (typically 20)"
        },
	{
        	" Local_OptimizeSkewVector",		/// name of the C++ function
            " RRRRRRRR",						/// 7 parametres = 6 d'entree + 1 parametre de retour 
            " ARM_GP_OptimizeSkewVector",
			" target,weight,precision,Init, LBound, UBound, aldo",
            " 1",							/// visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
	},
    {
        	" Local_PXL_OptimizeSkewVector",	/// name of the C++ function
            " RRRRRRRR",						/// 7 parametres = 6 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_OptimizeSkewVector",
			" target,weight,precision,Init, LBound, UBound, aldo",
            " 0",							/// not visible in excel
            XLLOCALARM_SEC_GROUP,
            " ",
            " ",
    },
	{
        		" Local_StochasticVol_N_VanillaOption",		/// name of the C++ function
                " RRRRRRRRRR" ,					/// 10 parametres = 9 d'entree + 1 parametre de retour 
                " ARM_CF_StochasticVol_N_VanillaOption",
				" forward,strike,time_to_maturity,r,StdDev,VolDrift,VolVol,callput,[nbterms]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the forward value of a european vanilla option with a Normal model with stochastic volatility",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" time to reset of the option",
				" interest rate (0.05 for 5% vol)",
				" Standard deviation (0.15 for 15% )",
				" VolDrift: drift of the volatility",
				" VolVol: Volatility of the volatility",
				" CALL or PUT",
				" nb of Legendre terms (typically 20)"
        },
		{
        		" Local_GLambda_From_SABR_Calibrate",		/// name of the C++ function
                " RRRRRRRRRRRRRRRR" ,					/// 16 parametres = 14 d'entree + 1 parametre de retour 
                " ARM_CF_GLambda_From_SABR_Calibrate",
				" forward,alpha,beta,rho,nu,T,SabrType,scope,Initial_L1,Initial_L2,Initial_L3,Initial_L4,Initial_L5,Initial_L6,algorithm",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the calibrated GLambda Parameters from SABR Parameters",
				" Forward (ofthe SABR model) (real value: 0.05 for 5%)",
				" alpha (ofthe SABR model) (real value: 0.05 for 5%)",
				" beta (ofthe SABR model)",
				" rho (ofthe SABR model)",
				" nu (ofthe SABR model) (real value: 0.05 for 5%)",
				" maturity (in years)",
				" Type of SABR formula : SABR_IMPLNVOL(2),SABR_A,SABR_G,DIRECTEXACT,DIRECTGEOMETRIC,DIRECTARITHMETIC,NORMALEXACT,NORMALGEOMETRIC,NORMALARITHMETIC,ANALYTICZP0,ANALYTICZP2",	
				" 0.2 means strikes from 0.8 to 1.2 are used ",
				" inital value for l1 " ,
				" inital value for l2  >0",
				" inital value for l3  >0  ",
				" inital value for l4  >0",
				" inital value for l5  >0  ",
				" inital value for l6  >0  ",
				" Algorithm: Values= NAG_OPT_NLIN_LSQ,NAG_OPT_LSQ_DERIV,.."
				
		},
		{
        		" Local_Student_GLambda_Power_SpreadOption",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRRR",					/// 21 parametres = 21 d'entree + 1 parametre de retour 
                " ARM_CF_GlambdaStudent_GSpreadOpt",
				" l1a,l2a,l3a,l4a,l5a,l6a,l1b,l2b,l3b,l4b,l5b,l6b,correlation,degre,b1,k1,a2,b2,k2,n",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the  FORWARD VALUE of a power spread option (GLambda+Student)",
				" underlying 1 : L1",
				" underlying 1 : L2",
				" underlying 1 : L3",
				" underlying 1 : L4",
				" underlying 1 : L5",
				" underlying 1 : L6",
				" underlying 2 : L1",
				" underlying 2 : L2",
				" underlying 2 : L3",
				" underlying 2 : L4",
				" underlying 2 : L5",
				" underlying 2 : L6",
				" copula correlation: 0.5 for 50%",
				" copula degre:",
				" b1: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" k1: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" a2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" b2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" k2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" number of points in the numerical gauss integration"
        },
		{
        		" Local_Lambert_Function",		/// name of the C++ function
                " RR" ,					/// 4 parametres = 3 d'entree + 1 parametre de retour 
                " ARM_CF_Util_Lambert_Function",
				" z",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the main solution of w exp(w) = z",
				" z should be > -1/e "
        },
		{
        		" Local_CF_BlackSholesTimeValue",		/// name of the C++ function
                " RRRRRR",					/// 6 parametres = 5 d'entree + 1 parametre de retour 
                " ARM_CF_LN_Option_TimeValue",
				" forward,totalvolatility,bondprice,strike,CallPut",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the present value of the time value of a black and Scholes option",
				" forward value: 0.5 for 50%",
				" total volatility of the forward value: 0.5 for 50% (equal to v*sqrt[t])",
				" discount factor : 0.85 for 85%",
				" strike: 0.5 for 50%",
				" (C)all=1 or (P)ut=-1 "
        },
		{
        		" Local_CF_BlackSholesTimeValue_ImplicitVol",		/// name of the C++ function
                " RRRRRR",					/// 6 parametres = 5 d'entree + 1 parametre de retour 
                " ARM_CF_LN_Option_TimeValue_ImpVol",
				" forward,bondprice,strike,CallPut,optprice",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the present value of the black and Scholes option",
				" forward value: 0.5 for 50%",
				" discount factor : 0.85 for 85%",
				" strike: 0.5 for 50%",
				" (C)all=1 or (P)ut=-1 ",
				" option present price "
        },
		{
        		" Local_ShiftedLogNormal_Quantile",		/// name of the C++ function
                " RRRRRR" ,					/// 4 parametres = 3 d'entree + 1 parametre de retour 
                " ARM_CF_Util_ShiftedLogNormal_Quantile",
				" f,p,t,sigma,m",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " returns the quantile of the shifted lognormal distribution",
				" f ",
				" probability (parameter) ",
				" t ",
				" sigma ",
				" shift "
        },
		{
        		" Local_ShiftedLogNormal_Distribution",		/// name of the C++ function
                " RRRRRR" ,					/// 4 parametres = 3 d'entree + 1 parametre de retour 
                " ARM_CF_Util_ShiftedLogNormal_Distribution",
				" f,x,t,sigma,m",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " returns the quantile of the shifted lognormal distribution",
				" f ",
				" shifted Log values (parameter) ",
				" t ",
				" sigma ",
				" shift "
        },
		{
        		" Local_SABR_Quantile",		/// name of the C++ function
                " RRRRRRRRRRRRR" ,					/// 13 parametres = 12 d'entree + 1 parametre de retour 
                " ARM_CF_Util_SABR_Quantile",
				" f,p,T,alpha, beta, rho ,nu,SABR Type, nbsteps,alpha_exp,alpha_tanh,kb_tanh",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " returns the quantile of the SABR distribution",
				" forward",
				" probability (parameter)",
				" maturity",
				" alpha ",
				" beta ",
				" rho ",
				" nu ",
				" Type of SABR formula: SABR_IMPLNVOL(2),SABR_DIRECT_SC1(2),ANALYTICZP0",	
				" nb of steps in gaussian integration of analytical formula (def=120) ",
				" Strike Cuter Exp : alpha (def=0.01)",
				" Strike Cuter Tanh : alpha (def=1.5) ",
				" Strike Cuter Tanh : Kb (def=0.03)"
			
        },
		{
        		" Local_SABR_Distribution",		/// name of the C++ function
                " RRRRRRRRRRRRR" ,					/// 13 parametres = 12 d'entree + 1 parametre de retour 
                " ARM_CF_Util_SABR_Distribution",
				" f,x,T,alpha, beta, rho ,nu,SABR Type, nbsteps,alpha_exp,alpha_tanh,kb_tanh",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " returns the quantile of the SABR distribution",
				" forward",
				" SABR Variable value (parameter)",
				" maturity",
				" alpha ",
				" beta ",
				" rho ",
				" nu ",
				" Type of SABR formula: SABR_IMPLNVOL(2),SABR_DIRECT_SC1(2),ANALYTICZP0",	
				" nb of steps in gaussian integration of analytical formula (def=120) ",
				" Strike Cuter Exp : alpha (def=0.01)",
				" Strike Cuter Tanh : alpha (def=1.5) ",
				" Strike Cuter Tanh : Kb (def=0.03)"
			
        },
		{
        		" Local_Student_GLambda_Power_Digital_SpreadOption",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRR",					/// 19 parametres = 18 d'entree + 1 parametre de retour 
                " ARM_CF_GlambdaStudent_Digital_SpreadOpt",
				" l1a,l2a,l3a,l4a,l5a,l6a,l1b,l2b,l3b,l4b,l5b,l6b,correlation,degre,a2,b2,k2,n",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the  FORWARD VALUE of a power spread option (GLambda+Student)",
				" underlying 1 : L1",
				" underlying 1 : L2",
				" underlying 1 : L3",
				" underlying 1 : L4",
				" underlying 1 : L5",
				" underlying 1 : L6",
				" underlying 2 : L1",
				" underlying 2 : L2",
				" underlying 2 : L3",
				" underlying 2 : L4",
				" underlying 2 : L5",
				" underlying 2 : L6",
				" copula correlation: 0.5 for 50%",
				" copula degre:",
				" a2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" b2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" k2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" number of points in the numerical gauss integration"
        },
		{
        		" Local_Student_GLambda_Power_Index2Digital_SpreadOption",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRR",					/// 19 parametres = 18 d'entree + 1 parametre de retour 
                " ARM_CF_GlambdaStudent_Index2Digital_SpreadOpt",
				" l1a,l2a,l3a,l4a,l5a,l6a,l1b,l2b,l3b,l4b,l5b,l6b,correlation,degre,a2,b2,k2,n",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the  FORWARD VALUE of a power spread option (GLambda+Student)",
				" underlying 1 : L1",
				" underlying 1 : L2",
				" underlying 1 : L3",
				" underlying 1 : L4",
				" underlying 1 : L5",
				" underlying 1 : L6",
				" underlying 2 : L1",
				" underlying 2 : L2",
				" underlying 2 : L3",
				" underlying 2 : L4",
				" underlying 2 : L5",
				" underlying 2 : L6",
				" copula correlation: 0.5 for 50%",
				" copula degre:",
				" a2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" b2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" k2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" number of points in the numerical gauss integration"
        },
		{
        		" Local_Student_GLambda_Power_Index1Digital_SpreadOption",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRR",					/// 19 parametres = 18 d'entree + 1 parametre de retour 
                " ARM_CF_GlambdaStudent_Index1Digital_SpreadOpt",
				" l1a,l2a,l3a,l4a,l5a,l6a,l1b,l2b,l3b,l4b,l5b,l6b,correlation,degre,a2,b2,k2,n",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the  FORWARD VALUE of a power spread option (GLambda+Student)",
				" underlying 1 : L1",
				" underlying 1 : L2",
				" underlying 1 : L3",
				" underlying 1 : L4",
				" underlying 1 : L5",
				" underlying 1 : L6",
				" underlying 2 : L1",
				" underlying 2 : L2",
				" underlying 2 : L3",
				" underlying 2 : L4",
				" underlying 2 : L5",
				" underlying 2 : L6",
				" copula correlation: 0.5 for 50%",
				" copula degre:",
				" a2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" b2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" k2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" number of points in the numerical gauss integration"
        },
		{
        		" Local_Student_GLambda_SpreadOption",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRR",					/// 19 parametres = 18 d'entree + 1 parametre de retour 
                " ARM_CF_GlambdaStudent_SpreadOpt",
				" l1a,l2a,l3a,l4a,l5a,l6a,l1b,l2b,l3b,l4b,l5b,l6b,correlation,degre,k2,n",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the  FORWARD VALUE of a power spread option (GLambda+Student)",
				" underlying 1 : L1",
				" underlying 1 : L2",
				" underlying 1 : L3",
				" underlying 1 : L4",
				" underlying 1 : L5",
				" underlying 1 : L6",
				" underlying 2 : L1",
				" underlying 2 : L2",
				" underlying 2 : L3",
				" underlying 2 : L4",
				" underlying 2 : L5",
				" underlying 2 : L6",
				" copula correlation: 0.5 for 50%",
				" copula degre:",
				" k2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" number of points in the numerical gauss integration"
        },
		{
        		" Local_GLambda_Distribution",		/// name of the C++ function
                " RRRRRRRR" ,					/// 8 parametres = 7 d'entree + 1 parametre de retour 
                " ARM_CF_Util_GLambda_Distribution",
				" l1,l2,l3,l4,l5,l6,x",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " returns the Distribution of the Glambda distribution",
				" l1 ",
				" l2 ",
				" l3 ",
				" l4 ",
				" l5 ",
				" l6 ",
				" x "
        },
		{
        		" Local_GLambda_Quantile",		/// name of the C++ function
                " RRRRRRRR" ,					/// 8 parametres = 7 d'entree + 1 parametre de retour 
                " ARM_CF_Util_GLambda_Quantile",
				" l1,l2,l3,l4,l5,l6,x",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " returns the quantile of the shifted Glambda distribution",
				" l1 ",
				" l2 ",
				" l3 ",
				" l4 ",
				" l5 ",
				" l6 ",
				" x "
        },
		{
        		" Local_Student_Quantile",		/// name of the C++ function
                " RRR" ,					/// 4 parametres = 3 d'entree + 1 parametre de retour 
                " ARM_CF_Util_Student_Quantile",
				" rank,x",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " returns the quantile of the student distribution",
				" degre ",
				" x "
        },
		{
        		" Local_Student_Distribution",		/// name of the C++ function
                " RRR" ,					/// 4 parametres = 3 d'entree + 1 parametre de retour 
                " ARM_CF_Util_Student_Distribution",
				" rank,x",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " returns the Distribution of the student distribution",
				" degre ",
				" x "
        },
		{
				" Local_ImcompleteBeta_Inverse",		/// name of the C++ function
                " RRRR" ,					/// 4 parametres = 3 d'entree + 1 parametre de retour 
                " ARM_CF_Util_ImcompleteBeta_Inverse",
				" a,b,x",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " returns the imcomplete beta regularized",
				" a ",
				" b ",
				" x "
        },
		{
				" Local_Student_QIntegral",		/// name of the C++ function
                " RRRR" ,					/// 4 parametres = 3 d'entree + 1 parametre de retour 
                " ARM_CF_Util_Student_QIntegral",
				" nu,H,x",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " returns the Student Q Integral",
				" nu : degre ",
				" H : niveau ",
				" x "
        },
		{
				" Local_Normal_ImpliedVol",		/// name of the C++ function
                " RRRRR" ,					/// 5 parametres = 4 d'entree + 1 parametre de retour 
                " ARM_CF_Norm_ImpliedVol",
				" F,K,opt,CallPut",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " returns the Normal_ImpliedVol",
				" forward ",
				" strike ",
				" option price ",
				" Call or Put"
        },
		{
				" Local_Normal_Digital_ImpliedVol",		/// name of the C++ function
                " RRRRR" ,					/// 5 parametres = 4 d'entree + 1 parametre de retour 
                " ARM_CF_Norm_Digital_ImpliedVol",
				" F,K,opt,CallPut",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " returns the Normal_Digital_ImpliedVol",
				" forward ",
				" strike ",
				" option price ",
				" Call or Put"
        },
		{
				" Local_Hypergeometric_Whittaker_W",		/// name of the C++ function
                " RRRR" ,					/// 4 parametres = 3 d'entree + 1 parametre de retour 
                " ARM_CF_Util_Hypergeometric_Whittaker_W",
				" a,b,x",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " returns the confluent whittaker W",
				" a : ",
				" b : ",
				" x "
        },
		{
				" Local_Hypergeometric_Whittaker_M",		/// name of the C++ function
                " RRRR" ,					/// 4 parametres = 3 d'entree + 1 parametre de retour 
                " ARM_CF_Util_Hypergeometric_Whittaker_M",
				" a,b,x",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " returns the confluent whittaker M",
				" a : ",
				" b : ",
				" x "
        },
		{
				" Local_Bessel_Y",		/// name of the C++ function
                " RRR" ,					/// 4 parametres = 3 d'entree + 1 parametre de retour 
                " ARM_CF_Util_Bessel_Y",
				" nu,x",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " returns the Bessel Y",
				" nu : ",
				" x "
        },
		{
				" Local_Bessel_I",		/// name of the C++ function
                " RRR" ,					/// 4 parametres = 3 d'entree + 1 parametre de retour 
                " ARM_CF_Util_Bessel_I",
				" nu,x",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " returns the Bessel Y",
				" nu : ",
				" x "
        },
		{
				" Local_Bessel_J",		/// name of the C++ function
                " RRR" ,					/// 4 parametres = 3 d'entree + 1 parametre de retour 
                " ARM_CF_Util_Bessel_J",
				" nu,x",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " returns the Bessel J",
				" nu : ",
				" x "
        },
		{
				" Local_Bessel_K",		/// name of the C++ function
                " RRR" ,					/// 4 parametres = 3 d'entree + 1 parametre de retour 
                " ARM_CF_Util_Bessel_K",
				" nu,x",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " returns the Bessel Y",
				" nu : ",
				" x "
        },

		{
        		" Local_Shifted_Heston_VanillaOption",		/// name of the C++ function
                " RRRRRRRRRRRRRRRR" ,					/// 12 parametres = 10 d'entree + 1 parametre de retour 
                " ARM_CF_Shifted_Heston_Option",
				" forward,strike ,Variance initiale,timetomaturity ,Long Term Variance,mean reveting speed,vol of vol,correlation, shift,callput,[n1],[n],[nS],[nO],[prec]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the Shifted Heston value of the forward value of a vanilla option",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" (real value: 0.05 for 5%)",
				" time to maturity (years)",
				" (0.15 for 15% vol)",
				" mean reverting speed",
				" volatility of volatility",
				" correlation between the volatility and the underlying (-0.5 for -50%)",
				" shift ",
				" CALL or PUT",
				" number of Legendre points for the integration (Stage 1) ",
				" number of Legendre points for the integration ",
				" number of stage (default 7)",
				" number of of ascillation/stage (default 3) ",
				" relative precision (defaul=1e-5)"
        },

		{
        		" Local_SABR_Heston_VanillaOption",		/// name of the C++ function
                " RRRRRRRRRRRRRRRR" ,					/// 12 parametres = 10 d'entree + 1 parametre de retour 
                " ARM_CF_SABR_Heston_Option",
				" forward,strike ,Vol initiale,timetomaturity ,long term limit,mean reveting speed,vol of vol,correlation, beta,callput,[n1],[n],[nS],[nO],[prec]",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the Shifted Heston value of the forward value of a vanilla option",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" vol initiale (real value: 0.05 for 5%)",
				" time to maturity (years)",
				" long terme limit (0.15 for 15% vol)",
				" mean reverting speed",
				" volatility of volatility",
				" correlation between the volatility and the underlying (-0.5 for -50%)",
				" beta ",
				" CALL or PUT",
				" number of Legendre points for the integration (Stage 1) ",
				" number of Legendre points for the integration ",
				" number of stage (default 7)",
				" number of of ascillation/stage (default 3) ",
				" relative precision (defaul=1e-5)"
        },
		{
        		" Local_GLambda_CompleteSpreadoption",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRR",					/// 20 parametres = 19 d'entree + 1 parametre de retour 
                " ARM_CF_GlambdaStudent_Complete_SpreadOpt",
				" l1a,l2a,l3a,l4a,l5a,l6a,l1b,l2b,l3b,l4b,l5b,l6b,DiscVec,correlation,degre,k2,n",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the  Present Value of a complete spread option (GLambda+Student)",
				" underlying 1 : L1 Vector",
				" underlying 1 : L2 Vector",
				" underlying 1 : L3 Vector",
				" underlying 1 : L4 Vector",
				" underlying 1 : L5 Vector",
				" underlying 1 : L6 Vector",
				" underlying 2 : L1 Vector",
				" underlying 2 : L2 Vector",
				" underlying 2 : L3 Vector",
				" underlying 2 : L4 Vector",
				" underlying 2 : L5 Vector",
				" underlying 2 : L6 Vector",
				" Discount Vector ",
				" copula correlation: 0.5 for 50%",
				" copula degre:",
				" k2: (of Max(a1*S1+b1*S2-k1,0) if a2*S1+b2*S2-k2>0)",
				" number of points in the numerical gauss integration"
        },
		{
        		" Local_JumpDiffusion_Mepi_Call",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRRR",					/// 20 parametres = 19 d'entree + 1 parametre de retour 
                " ARM_CF_JumpDiffusion_Mepi_Call",
				" P0,f0,T,K,R,Emin,Lmax,gamma0,gamma1,sig,lambda,sigJ,r,s,mu,fees,voldrift,volvol,CallPut,params",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the  Present Value of a Vanilla Mepi Option ",
				" Initial Value of the Portfolio",
				" Initial Value ofthe Cushion",
				" Maturity",
				" Strike on the cushion",
				" Risk Factor",
				" Minimum value invested",
				" Maximum factor of portfolio borrowed",
				" Initial garantied level (0.8 for 80%)",
				" Final garantied level (1.0 for 100%)",
				" volatility of the underlying index",
				" jump probability",
				" Variance of the jump size",
				" instantaneous rate",
				" borrowing spread",
				" underlying drift",
				" Annual fees (0.02 for 2%)",
				" Volatility Drift",
				" Volatility of volatility",
				" (C)all / (P)ut",
				" Params : 1:log space acceleration,2:Schema (ex,imp,mid) ",
			
        },
		{
        		" Local_Util_TrigonalSolve",		/// name of the C++ function
                " RRRRR" ,					/// 5 parametres = 4 d'entree + 1 parametre de retour 
                " ARM_CF_Util_TrigonalSolve",
				" Underdiagonal, diagonal, upperdiagonal, Constant",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives Solution of the trigonal system",
				" Vector of Under diagonals",
				" Vector of diagonals",
				" Vector of upper diagonals",
				" Vector of constant"
			
		},
		{
        		" Local_BiSABR_SpreadOption",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRR" ,					/// 18 parametres = 17 d'entree + 1 parametre de retour 
                " ARM_CF_BiSABR_SpreadOption",
				"  F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2, K, T,callput, rhos, rhov, rhoc12, rhoc21,flag",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the forward value of a BiSABR vanilla spreadoption",
				" Forward underlying 1",
				" SABR alpha 1",
				" SABR beta 1",
				" SABR rho 1",
				" SABR nu 1",
				" Forward underlying 2",
				" SABR alpha 2",
				" SABR beta 2",
				" SABR rho 2",
				" SABR nu 2",
				" strike spreadoption ",
				" maturity spreadoption ",
				" Call or Put ",
				" correlation under1 / under2",
				" correlation vol1 / vol2 ",
				" correlation under1 / vol 2",
				" correlation under2 / vol 1",
				" flag: normal=0,bilog corrected=1,complex=10"
			
		},
		{
        		" Local_Hypergeometric_Appell",		/// name of the C++ function
                " RRRRRRRRRR" ,					/// 8 parametres = 7 d'entree + 1 parametre de retour 
                " ARM_CF_Util_Hypergeometric_Appell",
				"  a,b1,b2,c,x,xim,y,yim,nb",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the hypergeometric appell F_2",
				" a",
				" b1",
				" b2",
				" c",
				" x ",
				" xim :(x^2+xim^2) should be <1",
				" y ",
				" yim :(y^2+yim^2) should be <1",
				" nb iteration"	
		},
		{
        		" Local_Util_Eigenvalues4",		/// name of the C++ function
                " RRRRRRR" ,					/// 5 parametres = 4 d'entree + 1 parametre de retour 
                " ARM_CF_Util_Eigenvalues4",
				"  rho12, rho13,  rho14, rho23, rho24, rho34",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives eigenvalues of the 4X4 correlation matrix",
				" rho12",
				" rho13",
				" rho14",
				" rho23",
				" rho24",
				" rho34"
			
		},
		{
        		" Local_Util_Eigenvalues3",		/// name of the C++ function
                " RRRR" ,					/// 5 parametres = 4 d'entree + 1 parametre de retour 
                " ARM_CF_Util_Eigenvalues3",
				"  rho12, rho13, rho23",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives eigenvalues of the 3X3 correlation matrix",
				" rho12",
				" rho13",
				" rho23"
			
		},
		{
        		" Local_BiSABR_Calibrate",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRR" ,					/// 17 parametres = 16 d'entree + 1 parametre de retour 
                " ARM_CF_BiSABR_Calibrate",
				" F1_vec, alpha1_vec,beta1_vec ,rho1_vec,nu1_vec,F2_vec, alpha2_vec,beta2_vec ,rho2_vec,nu2_vec,strike_vec,Maturity_vec,callput_vec,price_vec,weights_vec,initialparams",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " returns the calibrated correlations for a BiSABR Model",
				" F1 vector",
				" alpha1 vector ",
				" beta1 vector ",
				" rho1 vector ",
				" nu1 vector ",
				" F2 vector ",
				" alpha2 vector ",
				" beta2 vector ",
				" rho2 vector ",
				" nu2 vector ",
				" strike vector ",
				" maturity vector ",
				" call put vector ",
				" price vector ",
				" weight vector, ",
				" param initiaux ; rhos,rhov,rhoc12,rhoc21 "
		},
{
        		" Local_LN_DigitalOption",		/// name of the C++ function
                " RRRRRR" ,					/// 6 parametres = 5 d'entree + 1 parametre de retour 
                " ARM_CF_LN_DigitalOption",
				" forward, strike,maturity ,callput,volatility",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " returns the digital vanilla option in the LN Model",
				" forward",
				" strike ",
				" maturity ",
				" callput ",
				" volatility "
		},
		{
        		" Local_Util_BiSABR_CorrelationEvolution",		/// name of the C++ function
                " RRRRRRRRR" ,					/// 9 parametres = 8 d'entree + 1 parametre de retour 
                " ARM_CF_Util_BiSABR_CorrelationEvolution",
				"   rho1,rho2,rhos,rhov,rhoc12,rhoc21,newrho1,newrho2",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives acceptable values for newrhov,newrhoc12,newrhoc21",
				" rho1",
				" rho2",
				" rhos",
				" rhov",
				" rhoc12",
				" rhoc21",
				" newrho1",
				" newrho2"
			
		},

		{
        		" Local_LN_RatioOption",		/// name of the C++ function
                " RRRRRRRRRRR" ,				/// 11 parametres = 10 d'entree + 1 parametre de retour 
                " ARM_CF_LN_RatioOption",
				"   L1,Mu1,Sigma1,L2,Mu2,Sigma2,rho,K,T,Callput",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " return E[S2/S1 1(S1>K)] for the call",
				" L1",
				" Mu1",
				" Sigma1",
				" L2",
				" Mu2",
				" Sigma2",
				" rho",
				" K",
				" T",
				"CallPut"
			
		},
		{
        		" Local_LN_ProductOption",		/// name of the C++ function
                " RRRRRRRRRRR" ,				/// 11 parametres = 10 d'entree + 1 parametre de retour 
                " ARM_CF_LN_ProductOption",
				"   L1,Mu1,Sigma1,L2,Mu2,Sigma2,rho,K,T,Callput",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " return E[S1*S2 1(S1>K)] for the call",
				" L1",
				" Mu1",
				" Sigma1",
				" L2",
				" Mu2",
				" Sigma2",
				" rho",
				" K",
				" T",
				"CallPut"
			
		},
{
        		" Local_BS_EuroBarriere_ImpliedVol",		/// name of the C++ function
                " RRRRRRRRRRR" ,					/// 10 parametres = 9 d'entree + 1 parametre de retour 
                " ARM_CF_LN_BarrierOption_ImpliedVol",
				" forward, strike, barrier, rebate, OptPrice, expiry, discount,callput, in or out, up or down",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the implied LN volatility of a barrier option",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" barriere (real value: 0.05 for 5%)",
				" rebate (real value: 0.05 for 5%)",
				" option price",
				" expiry (years: 0.5 for 6m)",
				" discount (0.95 for exemple)",
				" CALL or PUT",
				" IN or OUT Barrier ",
				" Up or DOWN Barrier"
        },
		{
        		" Local_BS_EuroDoubleBarriere_ImpliedVol",		/// name of the C++ function
                " RRRRRRRRRR" ,					/// 8 parametres = 7 d'entree + 1 parametre de retour 
                " ARM_CF_LN_DoubleBarrierOption_ImpliedVol",
				" forward,strike ,barrierUp ,barrierDown, OptPrice, expiry, interest rate,dividend yield,callput",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the implied LN volatility of a option double barrier ",
				" fwd value (real value: 0.05 for 5%)",
				" strike (real value: 0.05 for 5%)",
				" barriere up (real value: 0.05 for 5%)",
				" barriere down (real value: 0.05 for 5%)",
				" option price",
				" expiry (years: 0.5 for 6m)",
				" 0.05 for 5%",
				" 0.05 for 5%",
				" CALL or PUT"
        },
		{
        		" Local_SABR_GaussianSABRDigitalCall",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRR" ,					/// 17 parametres = 16 d'entree + 1 parametre de retour 
                " ARM_CF_SABR_GaussianSABRDigitalCall",
				" f1,alpha1,beta1,rho1,nu1,SABRFlag1,f2,alpha2,beta2,rho2,nu2,SABRFlag2,rho,K,T,legendreNb,alpha_exp,alpha_tanh,kb_tanh",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the Digital S1-S2-K With SABR marginals in a gaussian copula ",
				" forward value 1",
				" alpha1",
				" beta1",
				" rho1",
				" nu1",
				" SABR_IMPLNVOL(2),SABR_GEOMETRIC,SABR_ANALYTIC",
				" forward value 2",
				" alpha2",
				" beta2",
				" rho2",
				" nu2",
				" SABR_IMPLNVOL(2),SABR_GEOMETRIC,SABR_ANALYTIC",
				" copula correlation",
				" strike ",
				" maturity",
				" nb of gauss points (<180)",
				" alpha_exp(=0.1)",
				" alpha_tanh(=1.5)",
				" kb_tanh=(=0.02)"
        },
		{
        		" Local_SABR_GaussianSABRDigitalCallPayingS1",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRR" ,					/// 17 parametres = 16 d'entree + 1 parametre de retour 
                " ARM_CF_SABR_GaussianSABRDigitalCallPayingS1",
				" f1,alpha1,beta1,rho1,nu1,SABRFlag1,f2,alpha2,beta2,rho2,nu2,SABRFlag2,rho,K,T,legendreNb,alpha_exp,alpha_tanh,kb_tanh",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the Digital S1-S2-K With SABR marginals in a gaussian copula ",
				" forward value 1",
				" alpha1",
				" beta1",
				" rho1",
				" nu1",
				" SABR_IMPLNVOL(2),SABR_GEOMETRIC,SABR_ANALYTIC",
				" forward value 2",
				" alpha2",
				" beta2",
				" rho2",
				" nu2",
				" SABR_IMPLNVOL(2),SABR_GEOMETRIC,SABR_ANALYTIC",
				" copula correlation",
				" strike ",
				" maturity",
				" nb of gauss points (<180)",
				" alpha_exp(=0.1)",
				" alpha_tanh(=1.5)",
				" kb_tanh=(=0.02)"
        },
		{
        		" Local_SABR_GaussianSABRDigitalCallPayingS2",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRR" ,					/// 17 parametres = 16 d'entree + 1 parametre de retour 
                " ARM_CF_SABR_GaussianSABRDigitalCallPayingS2",
				" f1,alpha1,beta1,rho1,nu1,SABRFlag1,f2,alpha2,beta2,rho2,nu2,SABRFlag2,rho,K,T,legendreNb,alpha_exp,alpha_tanh,kb_tanh",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the Digital S1-S2-K With SABR marginals in a gaussian copula ",
				" forward value 1",
				" alpha1",
				" beta1",
				" rho1",
				" nu1",
				" SABR_IMPLNVOL(2),SABR_GEOMETRIC,SABR_ANALYTIC",
				" forward value 2",
				" alpha2",
				" beta2",
				" rho2",
				" nu2",
				" SABR_IMPLNVOL(2),SABR_GEOMETRIC,SABR_ANALYTIC",
				" copula correlation",
				" strike ",
				" maturity",
				" nb of gauss points (<180)",
				" alpha_exp(=0.1)",
				" alpha_tanh(=1.5)",
				" kb_tanh=(=0.02)"
        },

{
        		" Local_SABR_GaussianSABRDigitalCallPayingS3",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRR" ,					/// 19 parametres = 18 d'entree + 1 parametre de retour 
                " ARM_CF_SABR_GaussianSABRDigitalCallPayingS3",
				" f1,alpha1,beta1,rho1,nu1,SABRFlag1,f2,alpha2,beta2,rho2,nu2,SABRFlag2,f3,sigma3,correlations,K,T,params",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the Digital S1-S2-K With SABR marginals in a gaussian copula ",
				" forward value 1",
				" alpha1",
				" beta1",
				" rho1",
				" nu1",
				" SABR_IMPLNVOL(2),SABR_GEOMETRIC,SABR_ANALYTIC",
				" forward value 2",
				" alpha2",
				" beta2",
				" rho2",
				" nu2",
				" SABR_IMPLNVOL(2),SABR_GEOMETRIC,SABR_ANALYTIC",
				" forward 3"
				" lognormal volatility 3",
				" rho12,rho13,rho23",
				" strike ",
				" maturity",
				" nbpoints,alpha_exp,alpha_tanh,kb_tanh"
        },
		{
				" Local_TarnProxy_Create",
                " RRRRRRRRRRRRRRRRR",						// 16 parametres = 15 d'entree + 1 parametre de retour
                " ARM_CF_TarnProxy_Create",
                " resetDates,fwds,densityFunctors, discountFactors, levPrec,lev,fix,cap,floor,fees,dcf,target,globalCap,globalFloor,correlInput,nbSimul",
				" 1",
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
				" ",
				" Creates a fake tarn pricer"
        },
		{
				" Local_PXL_TarnProxy_Create",
                " RRRRRRRRRRRRRRRRR",						// 16 parametres = 15 d'entree + 1 parametre de retour
                " PXL_ARM_CF_TarnProxy_Create",
                " resetDates,fwds,densityFunctors, discountFactors, levPrec,lev,fix,cap,floor,fees,dcf,target,globalCap,globalFloor,correlInput,nbSimul",
				" 0",
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " "
        },
		{
        	" Local_TarnProxy_GetPrice",		/// name of the C++ function
            " RRR",						/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_CF_TarnProxy_GetPrice",
			" proxy, [nbPayoff]",
            " 1",								/// visible in excel
            XLLOCALARM_CLOSEDFORME_GROUP,
            " ",
            " ",
            " Price of Tarn Proxy", 
			" "
			" "
		},
		{
				" Local_VBMinMaxProxy_Create",
                " RRRRRRRRRRRRRRRRRRRRRR",						// 16 parametres = 15 d'entree + 1 parametre de retour
                " ARM_CF_VBMinMaxProxy_Create",
                " asof, resetDates,fwdrates, TotalVol,LeftVol,RightVol,Nu,Rho,nbsimul,SabrDiff,typevb,[type1sens],[type2sens],[rate1Lev],[rate1Add],[caprateLev],[caprateAdd],[vbLev],[maxchoice],[minchoice],[minmaxfreq]",
				" 1",
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
				" ",
				" Creates a fake vol bond min max pricer"
        },
		{
				" Local_PXL_VBMinMaxProxy_Create",
                " RRRRRRRRRRRRRRRRRRRRRR",						// 16 parametres = 15 d'entree + 1 parametre de retour
                " PXL_ARM_CF_VBMinMaxProxy_Create",
                " asof, resetDates,fwdrates, TotalVol,LeftVol,RightVol,Nu,Rho,nbsimul,SabrDiff,typevb,[type1sens],[type2sens],[choiceformax],[rate1Lev],[rate1Add],[caprateLev],[caprateAdd],[vbLev],[maxchoice],[minchoice],[minmaxfreq]",
				" 0",
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
				" ",
				" Creates a fake vol bond min max pricer"
        },
		{
        	" Local_VBMinMaxProxy_GetInfo",		/// name of the C++ function
            " RRR",						/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_CF_VBMinMaxProxy_GetInfo",
			" proxy, info",
            " 1",								/// visible in excel
            XLLOCALARM_CLOSEDFORME_GROUP,
            " ",
            " ",
            " Price, MaxRate or MinRate", 
			" "
			" "
		},
		{
				" Local_Berm2DatesProxy_Create",
				" RRRRRRRRRRRRRR",
				" ARM_CF_Berm2DatesProxy_Create",
				" asof, resetDates, fwdRates, strike, Annuity, ratevols, partialratevols, vols, vvols, rho, ratecorrel, nbsimul, typediff",
				" 1",
				XLLOCALARM_CLOSEDFORME_GROUP,
				" ",
				" ",
				" Creates a fake bermuda swaption two periods"
		},
		{		
				" Local_Berm2DatesProxy_GetPrice",
				" RRR",
				" ARM_CF_Berm2DatesProxy_GetPrice",
				" proxy, [info]",
				" 1",
				XLLOCALARM_CLOSEDFORME_GROUP,
				" ",
				" ",
				" Price"
		},
		{
                " RRRRRRRRRRRRRRRRRRRRRR",						// 16 parametres = 15 d'entree + 1 parametre de retour
                " ARM_CF_VBMinMaxProxy_Create",
                " asof, resetDates,fwdrates, spotvol, leftvol, rightvol, nu, rho, nbSimul, SABRDiff, typeprice, [type1sens], [type2sens], [firstRateLeverage], [firstRateAdd], [capfirstRateLeverage], [capStrike], [capVBLeverage], [maxChoice], [minChoice], [minmaxFreq]",
				" 1",
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
				" ",
				" Creates a fake vol bond min max pricer"
        },
		{
				" Local_PXL_VBMinMaxProxy_Create",
                " RRRRRRRRRRRRRRRRRRRRRR",						// 16 parametres = 15 d'entree + 1 parametre de retour
                " PXL_ARM_CF_VBMinMaxProxy_Create",
                " asof, resetDates,fwdrates, spotvol, leftvol, rightvol, nu, rho, nbSimul, SABRDiff, typeprice, [type1sens], [type2sens], [firstRateLeverage], [firstRateAdd], [capfirstRateLeverage], [capStrike], [capVBLeverage], [maxChoice], [minChoice], [minmaxFreq]",
				" 0",
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
				" ",
				" Creates a fake vol bond min max pricer"
        },
		{
        	" Local_VBMinMaxProxy_GetInfo",		/// name of the C++ function
            " RRR",						/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_CF_VBMinMaxProxy_GetInfo",
			" proxy, info",
            " 1",								/// visible in excel
            XLLOCALARM_CLOSEDFORME_GROUP,
            " ",
            " ",
            " Price, MaxRate or MinRate", 
			" "
			" "
		},
		{
				" Local_SpreadVBProxy_Create",
                " RRRRRRRRRRRRRRRRRRRRRRRRRRR",						// 16 parametres = 15 d'entree + 1 parametre de retour
                " ARM_CF_SpreadVBProxy_Create",
                " asof, resetDates,fwdrates1, ratevols1,partialratevols1,vols1,nu1,rho1,fwdrates2, ratevols2,partialratevols2,vols2,nu2,rho2,Rate1Rate2Correl,typevb,nbsimul,[fwdLev],[fwdStrikes],[type1sens],[type2sens],[sabrdiff], [RateCorrelMeanRev], [RateCorrelVol],[FixedCpn],[Leverage]",
				" 1",
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
				" ",
				" Creates a fake vol bond min max pricer"
        },
		{
        	" Local_SpreadVBProxy_GetInfo",		/// name of the C++ function
            " RRR",						/// 3 parametres = 2 d'entree + 1 parametre de retour 
            " ARM_CF_SpreadVBProxy_GetInfo",
			" proxy, info",
            " 1",								/// visible in excel
            XLLOCALARM_CLOSEDFORME_GROUP,
            " ",
            " "
		},
		{
				" Local_Berm2DatesProxy_Create",
				" RRRRRRRRRRRRRR",
				" ARM_CF_Berm2DatesProxy_Create",
				" asof, resetDates, fwdRates, strike, Annuity, ratevols, partialratevols, vols, vvols, rho, ratecorrel, nbsimul, typediff",
				" 1",
				XLLOCALARM_CLOSEDFORME_GROUP,
				" ",
				" ",
				" Creates a fake bermuda swaption two periods"
		},
		{		
				" Local_Berm2DatesProxy_GetPrice",
				" RRR",
				" ARM_CF_Berm2DatesProxy_GetPrice",
				" proxy, [info]",
				" 1",
				XLLOCALARM_CLOSEDFORME_GROUP,
				" ",
				" ",
				" Price"
		},
		{
        		" Local_BiSABR_Digital_SpreadOption",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRR" ,					/// 18 parametres = 17 d'entree + 1 parametre de retour 
                " ARM_CF_BiSABR_Digital_SpreadOption",
				"  F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2, K, T,callput, rhos, rhov, rhoc12, rhoc21,flag",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " forward value of E[1{S1-S2-K>0}] for a call",
				" Forward underlying 1",
				" SABR alpha 1",
				" SABR beta 1",
				" SABR rho 1",
				" SABR nu 1",
				" Forward underlying 2",
				" SABR alpha 2",
				" SABR beta 2",
				" SABR rho 2",
				" SABR nu 2",
				" strike spreadoption ",
				" maturity spreadoption ",
				" Call or Put ",
				" correlation under1 / under2",
				" correlation vol1 / vol2 ",
				" correlation under1 / vol 2",
				" correlation under2 / vol 1",
				" flag: normal=0,bilog corrected=1,complex=10"
			
		},
			{
        		" Local_BiSABR_DigitalPaysS1_SpreadOption",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRR" ,					/// 18 parametres = 17 d'entree + 1 parametre de retour 
                " ARM_CF_BiSABR_DigitalPaysS1_SpreadOption",
				"  F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2, K, T,callput, rhos, rhov, rhoc12, rhoc21,flag",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " forward value of E[S1*1{S1-S2-K>0}] for a call",
				" Forward underlying 1",
				" SABR alpha 1",
				" SABR beta 1",
				" SABR rho 1",
				" SABR nu 1",
				" Forward underlying 2",
				" SABR alpha 2",
				" SABR beta 2",
				" SABR rho 2",
				" SABR nu 2",
				" strike spreadoption ",
				" maturity spreadoption ",
				" Call or Put ",
				" correlation under1 / under2",
				" correlation vol1 / vol2 ",
				" correlation under1 / vol 2",
				" correlation under2 / vol 1",
				" flag: normal=0,bilog corrected=1,complex=10"
			
		},
		{
        		" Local_BiSABR_DigitalPaysS2_SpreadOption",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRR" ,					/// 18 parametres = 17 d'entree + 1 parametre de retour 
                " ARM_CF_BiSABR_DigitalPaysS2_SpreadOption",
				"  F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2, K, T,callput, rhos, rhov, rhoc12, rhoc21,flag",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " forward value of E[S2*1{S1-S2-K>0}] for a call",
				" Forward underlying 1",
				" SABR alpha 1",
				" SABR beta 1",
				" SABR rho 1",
				" SABR nu 1",
				" Forward underlying 2",
				" SABR alpha 2",
				" SABR beta 2",
				" SABR rho 2",
				" SABR nu 2",
				" strike spreadoption ",
				" maturity spreadoption ",
				" Call or Put ",
				" correlation under1 / under2",
				" correlation vol1 / vol2 ",
				" correlation under1 / vol 2",
				" correlation under2 / vol 1",
				" flag: normal=0,bilog corrected=1,complex=10"
			
		},
		{
        		" Local_BiSABR_DigitalPaysS3_SpreadOption",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRR" ,					/// 19 parametres = 18 d'entree + 1 parametre de retour 
                " ARM_CF_BiSABR_DigitalPaysS3_SpreadOption",
				"  F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2, rhos, rhov, rhoc12, rhoc21,S3params,K, T,callput,flag",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " forward value of E[S3*1{S1-S2-K>0}] for a call",
				" Forward underlying 1",
				" SABR alpha 1",
				" SABR beta 1",
				" SABR rho 1",
				" SABR nu 1",
				" Forward underlying 2",
				" SABR alpha 2",
				" SABR beta 2",
				" SABR rho 2",
				" SABR nu 2",
				" correlation under1 / under2",
				" correlation vol1 / vol2 ",
				" correlation under1 / vol 2",
				" correlation under2 / vol 1",
				" S3 parameters:F,sig,rho13,rho23",
				" strike spreadoption ",
				" maturity spreadoption ",
				" Call or Put ",
				" flag: normal=0,bilog corrected=1,complex=10"
				
			
		},


		{
        		" Local_BiSABR_Quantile",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRR" ,					/// 18 parametres = 17 d'entree + 1 parametre de retour 
                " ARM_CF_Util_BiSABR_Quantile",
				" f1,alpha1,beta1,rho1,nu1,f2,alpha2,beta2,rho2,nu2,rhos,rhov,rhoc12,rhoc21,x,T,flag",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " returns the quantile of the BiSABR distribution",
				" forward1 ",
				" alpha1 ",
				" beta1",
				" rho1",
				" nu1",
				" forward2 ",
				" alpha2 ",
				" beta2",
				" rho2",
				" nu2",
				" rhos",
				" rhov"
				" rho cross vol1, underlying2"
				" rho cross vol2, underlying1"
				" probability (parameter) ",
				" maturity ",
				" flag= 0,1,10"
        },
		{
        		" Local_BiSABR_Distribution",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRR" ,					/// 18 parametres = 17 d'entree + 1 parametre de retour 
                " ARM_CF_Util_BiSABR_Distribution",
				" f1,alpha1,beta1,rho1,nu1,f2,alpha2,beta2,rho2,nu2,rhos,rhov,rhoc12,rhoc21,x,T,flag",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " returns the quantile of the BiSABR distribution",
				" forward1 ",
				" alpha1 ",
				" beta1",
				" rho1",
				" nu1",
				" forward2 ",
				" alpha2 ",
				" beta2",
				" rho2",
				" nu2",
				" rhos",
				" rhov"
				" rho cross vol1, underlying2"
				" rho cross vol2, underlying1"
				" probability (parameter) ",
				" maturity ",
				" flag= 0,1,10"
        },
			{
        		" Local_BetaEqualZeroSABR",		/// name of the C++ function
                " RRRRRRRRR" ,					/// 9 parametres = 8 d'entree + 1 parametre de retour 
                " ARM_CF_SABR_BetaEqualZero_Option",
				" f,K,T,mu,alpha,rho,nu,CallPut",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " returns an option SABR with beta =0",
				" forward ",
				" Strike ",
				" Maturity",
				" VolDrift",
				" alpha ",
				" rho",
				" nu",
				" CALL/PUT",
				" rhov"
        },
			{
        		" Local_Shifted2LogNormal_Distribution",		/// name of the C++ function
                " RRRRRRRRRR" ,					/// 10 parametres = 9 d'entree + 1 parametre de retour 
                " ARM_CF_Util_SpreadShifted2LogNormal_Distribution",
				" f1,sigma1,f2,sigma2,alpha,rho,T,n,x ",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " returns the Distribution of the Shifted2LogNormal distribution",
				" forward1 ",
				" Volatility1 ",
				" forward2 ",
				" Volatility2 ",
				" alpha ",
				" rho ",
				" T ",
				" nb of integration steps",
				" x value (strike)"
        },
			{
        		" Local_Shifted2LogNormal_Quantile",		/// name of the C++ function
                " RRRRRRRRRR" ,					/// 10 parametres = 9 d'entree + 1 parametre de retour 
                " ARM_CF_Util_SpreadShifted2LogNormal_Quantile",
				" f1,sigma1,f2,sigma2,alpha,rho,T,n,x",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " returns the Quantile of the Shifted2LogNormal distribution",
				" forward1 ",
				" Volatility1 ",
				" forward2 ",
				" Volatility2 ",
				" alpha ",
				" rho ",
				" T ",
				" nb of integration steps",
				" probability "
        },
			{
        		" Local_BiSABR_S3_SpreadOption",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRR" ,					/// 18 parametres = 17 d'entree + 1 parametre de retour 
                " ARM_CF_BiSABR_SABR_S3_SpreadOption",
				" S1Params,S2Params,S3Params,rhos,rhov,rhc12,rhoc21,correlation,T,A1,B1,K1,A2,B2,K2,flag,nbsteps",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " pays A1*spread+B1*S3-K1 if (A2,B2,K2)>0",
				" spread:forward1 params(f,alpha,beta,rho,nu) ",
				" spread:forward2 params(f,alpha,beta,rho,nu) ",
				" forward3 params(f,alpha,beta,rho,nu) ",
				" correlation between underl1 and underl2 ",
				" correlation between vol1 and vol2 ",
				" correlation between underl1 and vol2 ",
				" correlation between vol1 and underl2 ",
				" copula correlation ",
				" T ",
				" A1*(S1-S2)+B1*S3 -K1 if [A2*(S1-S2)+B2*S3 -K2]>0"
				" A1*(S1-S2)+B1*S3 -K1 if [A2*(S1-S2)+B2*S3 -K2]>0"
				" A1*(S1-S2)+B1*S3 -K1 if [A2*(S1-S2)+B2*S3 -K2]>0"
				" A1*(S1-S2)+B1*S3 -K1 if [A2*(S1-S2)+B2*S3 -K2]>0"
				" A1*(S1-S2)+B1*S3 -K1 if [A2*(S1-S2)+B2*S3 -K2]>0"
				" A1*(S1-S2)+B1*S3 -K1 if [A2*(S1-S2)+B2*S3 -K2]>0"
				" nb of integration steps",
        },
			{
				" Local_Heston_OptionPrice",
				" RRRRRRRRRRRRRR",
				" ARM_CF_HestonOptionPrice",
				" AsOfDate, ResetTime, Forward, Strike, CallOrPut, V0, Kappa, Rho, Theta, VVol, [shift], [Times], [Level]",
				" 1",
				XLLOCALARM_CLOSEDFORME_GROUP,
				" ",
				" ",
				" returns option price with heston model (either with time dependent or constant parameters)",
				" ",
				" ",
				" ",
				" ",
				" ",
				" Initial Variance",
				" Var Mean Reversion",
				" Forward / Variance correlation",
				" Long Term Variance",
				" Variance Volatility",
				" Relative Shift (default = 1)",
				" Time Pilars if time dependant level",
				" sigma (constant or time dependant)",
		},
		{
				" Local_SABR_To_Heston_SmileCalibration_Create",
				" RRRRRRRRRRRRRR",
				" ARM_CF_SABRToHestonSmileCalib_Create",
				" ResetTime, FwdRate, ATMVol, Alpha, Beta, Rho SABR, Nu, SABR Type, V0, Kappa, [Theta], [Rho], [Shift]",
				" 1",
				XLLOCALARM_CLOSEDFORME_GROUP,
				" ",
				" ",
		},
		{
				" Local_SABR_To_Heston_SmileCalibration_GetValue",
				" RRR",
				" ARM_CF_SABRToHestonSmileCalib_GetValue",
				" calibObj, info",
				" 1",
				XLLOCALARM_CLOSEDFORME_GROUP,
				" ",
		},
		{
				" Local_SmileCalib_Spread2Heston",
				" RRRRRRRRRRRRRRRRRRR",
				" ARM_CF_SmileCalib_Spread2Heston",
				" reset, fwd1, mktvols1, strikes1, constrvol1, constrK1, fwd2, mktvols2, strikes2, constrvol2, constrK2, mktvolspread, strikespread, constrvolSpread, constrKspread, v0, kappa, [theta]",
				" 1",
				XLLOCALARM_CLOSEDFORME_GROUP,
				" ",
		},
		{
				" Local_Spread2HestonVanilla",
				" RRRRRRRRRRRRRRRRRRR",
				" ARM_CF_Spread2HestonVanilla",
				" reset, fwd1, fwd2, strike, callput, v0, kappa, theta, nu, rho1, rho2, shift1, shift2, level1, level2, correl, [index1 leverage], [index2 leverage]",
				" 1",
				XLLOCALARM_CLOSEDFORME_GROUP,
				" ",
		},
		{
				" Local_Spread2Heston_TOTEMCalibration",
				" RRRRRRRRRRRRRRRRRRR",
				" ARM_CF_Spread2Heston_TOTEMCalibration",
				" TMat, TStrikes, TPrices, SchedReset, SchedDF, SchedFwd1, SchedFwd2, FwdCalibReset, Fwds1, Vols1, Strikes1, Fwds2, Vols2, Strikes2, v0, kappa, [theta], [ConstrCorrel]",
				" 1",
				XLLOCALARM_CLOSEDFORME_GROUP,
				" ",
				//" TOTEMMat, TOTEMStrikes, TOTEMPrices, FullSchedReset, FullSchedAnnuity, FullSchedFwd1, FullSchedFwd2, FwdCalibReset, LongFwds, LongVols, LongStrikes, ShortFwds, ShortVols, ShortStrikes, v0, kappa, [theta], [local]",
				//" TOTEMMat, TOTEMStrikes, TOTEMPrices, FullSchedReset, FullSchedAnnuity, FullSchedFwd1, FullSchedFwd2, FwdCalibReset, LongFwds, LongVols, LongStrikes, ShortFwds, ShortVols, ShortStrikes, v0, kappa, [theta], [locShifts]",
		},

	{
        	" Local_CIR_ModelParamsCreate",		/// name of the C++ function
            " RR",									/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " ARM_GP_CIR_ModelParamsCreate",
			" ModelParamsId",
            " 1",									/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Create the peacewize CIR model",
			" Volatility, MeanReversion, LongTermVol",
    },
	{
        	" Local_PXL_CIR_ModelParamsCreate",	/// name of the C++ function
            " RR",									/// 2 parametres = 1 d'entree + 1 parametre de retour 
            " PXL_ARM_GP_CIR_ModelParamsCreate",
			" ModelParamsId",
            " 0",									/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Create the peacewize CIR model",
			" Volatility, MeanReversion, LongTermVol",    
	},
	{
        	" Local_CIR_BondPrice",				
            " RRRRR",							
            " ARM_GP_CIR_BondPrice",
			" ModelParamsId,time,[r0],[x]",
            " 1",									/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Compute E( exp( -x.Int(rs,s=0..t) ) )",
			" Model Params",
			" Time from as of date (nbDays)",
			" r0 (default r0  =1.0)",
			" x  (default x  =1.0)",
    },
	{
        	" Local_PXL_CIR_BondPrice",				
            " RRRRR",								 
            " PXL_ARM_GP_CIR_BondPrice",
			" ModelParamsId,time,[r0]",
            " 0",									/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Compute E( exp( -x.Int(rs,s=0..t) ) )",
			" Model Params",
			" Time from as of date (nbDays)",
			" r0 (default r0  =1.0)", 
			" x  (default x  =1.0)",
	},
	{
        	" Local_CIR_BondDensity",				
            " RRRRRR",							
            " ARM_GP_CIR_BondDensity",
			" ModelParamsId,time,[r0],r,[Freqency]",
            " 1",									/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Compute the density using the inverse Laplace",
			" Model Params",
			" Time from as of date (nbDays)",
			" r0 (default r0 = 1.0)",
			" r ",
			" frequency(default frequency = 1.0)",
    },
	{
        	" Local_PXL_CIR_BondDensity",				
            " RRRRRRR",								 
            " PXL_ARM_GP_CIR_BondDensity",
			" ModelParamsId,time,[r0],r,[Freqency]",
            " 0",									/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Compute the density using the inverse Laplace",
			" Model Params",
			" Time from as of date (nbDays)",
			" r0 (default r0  =1.0)",
			" r ",
			" frequency(default frequency = 1.0)",
	},
	{
        	" Local_CIR_BondDistribution",				
            " RRRRRR",							
            " ARM_GP_CIR_BondDistribution",
			" ModelParamsId,time,[r0],nbDiscr,[Factor]",
            " 1",									/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Compute the cir distribution using the inverse Laplace",
			" Model Params",
			" Time from as of date (nbDays)",
			" r0 (default r0 = 1.0)",
			" Nb Discretization",
			" frequency(default factor = 1.0)",
    },
	{
        	" Local_PXL_CIR_BondDistribution",				
            " RRRRRR",							
            " PXL_ARM_GP_CIR_BondDistribution",
			" ModelParamsId,time,[r0],nbDiscr,[Factor]",
            " 0",									/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Compute the cir distribution using the inverse Laplace",
			" Model Params",
			" Time from as of date (nbDays)",
			" r0 (default r0 = 1.0)",
			" Nb Discretization",
			" frequency(default factor = 1.0)",
    },
	{
        	" Local_CIR_Bond_Weight",				
            " RR",							
            " ARM_GP_CIR_BondWeight",
			" ModelParamsId",
            " 1",									/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Compute P( Int(rs,s=0..t) < Infinity )",
			" Model Params",
    },
	{
        	" Local_CIR_Bond_Expectation",				
            " RR",							
            " ARM_GP_CIR_BondExpectation",
			" ModelParamsId",
            " 1",									/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Compute E( Int(rs,s=0..t))",
			" Model Params",
    },
	{
        	" Local_CIR_Bond_Variance",				
            " RR",							
            " ARM_GP_CIR_BondVariance",
			" ModelParamsId",
            " 1",									/// visible in excel
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Compute Var( Int(rs,s=0..t))",
			" Model Params",
    },
	{
        	" Local_EXP_RICCATI_Create",				
            " RRRRRRRRRR",								 
            " ARM_CF_EXP_RICCATI_Create",
			" alpha, beta, delta, lambda, x0, x1 ,x2 ,y0 ,t0",
            " 1",									
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " Equation dy/dt = delta*y(t)² + alpha*(x0-exp(-lambda*t))*y(t) + beta*(x1 - exp(-lambda*t) )*(x2 - exp(-lambda*t) ), et y(t0)=y0",
			" alpha",
			" beta",
			" delta",
			" lambda",
			" x0",
			" x1",
			" x2",
			" y0",
			" t0",
	},
	{
        	" Local_EXP_RICCATI_Price",				
            " RRR",								 
            " ARM_CF_EXP_RICCATI_Price",
			" object, t",
            " 1",									
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " solution at time t of the Exp Riccati",
			" Object",
			" t",
	},
	{
        	" Local_EXP_RICCATI_Int",				
            " RRRR",								 
            " ARM_CF_EXP_RICCATI_Int",
			" object, t, T",
            " 1",									
            XLLOCALARM_GENMODELS_GROUP,
            " ",
            " ",
            " solution at time t of the Exp Riccati",
			" Object",
			" t",
			" T",
	},

		{
        		" Local_Normal_Heston_VanillaCall",		/// name of the C++ function
                " RRRRRRRRRRRRRRRR" ,					/// 15 parametres = 14 d'entree + 1 parametre de retour 
                " ARM_CF_Normal_Heston_VanillaCall",
				" rho,lambdaV,thetaV,kappaV,V0,S0, k,T,lambdaB,callput,nbfirst,nb,NbStage,NbOscill,prec",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " returns the option price in a normal heston model",
				" correlation ",
				" mean reveresion ",
				" long term variance ",
				" vol of vol ",
				" initial variance ",
				" initial price ",
				" strike ",
				" maturity ",
				" regularisation parameter [default=0.45] ",
				" (C)all/(P)ut",
				" nb of integration steps for the first stage [default=40]",
				" nb of integration steps for the regular stage [default=40]",
				" nb of stage [default=-1]",
				" nb of oscillation per stage [default=0]",
				" precision to stop the integration [default=1e-8]",
		},
		{
					" Local_SABR_To_Heston_SmileCalibration_Create",
						" RRRRRRRRRRRRRR",
						" ARM_CF_SABRToHestonSmileCalib_Create",
						" ResetTime, FwdRate, ATMVol, Alpha, Beta, Rho SABR, Nu, SABR Type, V0, Kappa, [Theta], [Rho], [Shift]",
						" 1",
						XLLOCALARM_CLOSEDFORME_GROUP,
						" ",
						" "
		},
		{
						" Local_SABR_To_Heston_SmileCalibration_GetValue",
							" RRR",
							" ARM_CF_SABRToHestonSmileCalib_GetValue",
							" calibObj, info",
							" 1",
							XLLOCALARM_CLOSEDFORME_GROUP,
							" "
		},
		{
			" Local_SmileCalibration",
			" RRRRRRRRRR",
			" ARM_CF_SmileCalibration",
			" AsOfDate, ResetDates, Forwards, MktVols, Strikes, CalibParamId, [Constraint Strikes], [Constraint Vols], [Weights]",
			" 1",
			XLLOCALARM_CLOSEDFORME_GROUP,
			" ",
			" "
		},
		{
        		" Local_TRiSABR_VanillaOption",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRR" ,					/// 20 parametres = 19 d'entree + 1 parametre de retour 
                " ARM_CF_TRiSABR_VanillaOption",
				" Index1,Index2,Index3,rhos12,rhos23,rhos13,rhov12,rhov23,rhov13,rhoc12,rhoc21,rhoc23,rhoc32,rhoc13,rhoc31,K,T,CallPut,SABRtypeFlag",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " returns the option price of a tri SABR vanilla option",
				" f1,alpha1,beta1,rho1,nu1 ",
				" f2,alpha2,beta2,rho2,nu2 ",
				" f3,alpha3,beta3,rho3,nu3 ",
				" underlying correlation ",
				" underlying correlation ",
				" underlying correlation ",
				" volatility correlation ",
				" volatility correlation ",
				" volatility correlation ",
				" Cross underlying volatility correlation ",
				" Cross underlying volatility correlation ",
				" Cross underlying volatility correlation ",
				" Cross underlying volatility correlation ",
				" Cross underlying volatility correlation ",
				" Cross underlying volatility correlation ",
				" Strike",
				" Maturity",
				" Call (C) / Put (P)",
				" SABR_IMPLNVOL(2),ANALYTICZP0,SABR_A,SABR_G,.."
			},

{
        		" Local_Util_BiSABR_Eigenvalues",		/// name of the C++ function
                " RRRRRRR" ,					/// 5 parametres = 4 d'entree + 1 parametre de retour 
                " ARM_CF_Util_BiSABR_Eigenvalues",
				"  rho1,rho2,rhoS,rhoV,rhoc12,rhoc21",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives eigenvalues of the BiSABR correlation matrix",
				" rho1 SABR Correlation",
				" rho2 SABR Correlation",
				" S1-S2 correlation",
				" V1-V2 Correlation",
				" S1-V2 cross correlation",
				" V1-S2 cross correlation",
			
			
			
		},

		{
        		" Local_Util_TriSABR_Eigenvalues",		/// name of the C++ function
                " RRRRRRRRRRRRRRRR" ,					/// 16 parametres = 15 d'entree + 1 parametre de retour 
                " ARM_CF_Util_TriSABR_Eigenvalues",
				"  rho1,rho2,rho3,rhos12,rhos23,rhos13,rhov12,rhov23,rhov13,rhoc12,rhoc21,rhoc23,rhoc32,rhoc13,rhoc31",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives eigenvalues of the TriSABR correlation matrix",
				" rho1 SABR Correlation",
				" rho2 SABR Correlation",
				" rho3 SABR Correlation",
				" S1-S2 correlation",
				" S2-S3 correlation",
				" S1-S3 correlation",
				" V1-V2 Correlation",
				" V2-V3 Correlation",
				" V1-V3 Correlation",
				" S1-V2 cross correlation",
				" V1-S2 cross correlation",
				" S2-V3 cross correlation",
				" V2-S3 cross correlation",
				" S1-V3 cross correlation",
				" V1-S3 cross correlation"
			
		},

				{
        		" Local_Nonparametric_CompleteOption",		/// name of the C++ function
                " RRRRRRRRRRRRRRRRRR" ,					/// 17 parametres = 16 d'entree + 1 parametre de retour 
                " ARM_CF_Nonparametric_CompleteOption",
				"  Strike1_vec,Vol1_vec,Strike2_vec,Vol2_vec,S1params_vec,S2params_vec,correlation,maturity,a1,b1,k1,a2,b2,k2,nbsteps,algo,smiletype",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " nonparametric complete spreadoption",
				" Vector of strikes 1",
				" Vector of Volatility 1",
				" Vector of strike 2",
				" Vector of volatility 2",
				" S1 : indexBefore,indexAfter,flagBefore,flagAfter,Forward",
				" S2 : indexBefore,indexAfter,flagBefore,flagAfter,Forward",
				" Correlation",
				" Maturity",
				" a1",
				" b1",
				" k1",
				" a2",
				" b2",
				" k2",
				" nbsteps",
				" Algorithm",
				" Type=0(both positive) 1 (both real)"
			
		},
		{
        		" Local_Nonparametric_LogVolatility",		/// name of the C++ function
                " RRRRR" ,					/// 5 parametres = 4 d'entree + 1 parametre de retour 
                " ARM_CF_Nonparametric_LogVolatility",
				"  Strike1_vec,Vol1_vec,S1params_vec,k",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " nonparametric interpolation of a positive variable",
				" Vector of strikes 1",
				" Vector of Volatility 1",
				" S1 : indexBefore,indexAfter,flagBefore,flagAfter",
				" Strike"
			
		},
		{
        		" Local_Nonparametric_NormalVolatility",		/// name of the C++ function
                " RRRRR" ,					/// 5 parametres = 4 d'entree + 1 parametre de retour 
                " ARM_CF_Nonparametric_NormalVolatility",
				"  Strike1_vec,Vol1_vec,S1params_vec,k",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " nonparametric interpolation of a real variable",
				" Vector of strikes 1",
				" Vector of Volatility 1",
				" S1 : indexBefore,indexAfter,flagBefore,flagAfter",
				" Strike"
			
		},
		{
        		" Local_Nonparametric_NormalDistribution",		/// name of the C++ function
                " RRRRRRR" ,					/// 7 parametres = 6 d'entree + 1 parametre de retour 
                " ARM_CF_Nonparametric_NormalDistribution",
				"  Strike1_vec,Vol1_vec,S1params_vec,S,T,k",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " nonparametric Distribution of a  real variable",
				" Vector of strikes 1",
				" Vector of Volatility 1",
				" S1 : indexBefore,indexAfter,flagBefore,flagAfter",
				" Forward",
				" Maturity"
				" Strike"
			
		},
		{
        		" Local_Nonparametric_NormalQuantile",		/// name of the C++ function
                " RRRRRRR" ,					/// 7 parametres = 6 d'entree + 1 parametre de retour 
                " ARM_CF_Nonparametric_NormalQuantile",
				"  Strike1_vec,Vol1_vec,S1params_vec,S,T,k",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " nonparametric Quantile of a  real variable",
				" Vector of strikes 1",
				" Vector of Volatility 1",
				" S1 : indexBefore,indexAfter,flagBefore,flagAfter",
				" Forward",
				" Maturity"
				" Strike"
			
		},
				{
        		" Local_Nonparametric_LogNormalDistribution",		/// name of the C++ function
                " RRRRRRR" ,					/// 7 parametres = 6 d'entree + 1 parametre de retour 
                " ARM_CF_Nonparametric_LogNormalDistribution",
				"  Strike1_vec,Vol1_vec,S1params_vec,S,T,k",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " nonparametric Distribution of a positive  variable",
				" Vector of strikes 1",
				" Vector of Volatility 1",
				" S1 : indexBefore,indexAfter,flagBefore,flagAfter",
				" Forward",
				" Maturity"
				" Strike"
			
		},
		{
        		" Local_Nonparametric_LogNormalQuantile",		/// name of the C++ function
                " RRRRRRR" ,					/// 7 parametres = 6 d'entree + 1 parametre de retour 
                " ARM_CF_Nonparametric_LogNormalQuantile",
				"  Strike1_vec,Vol1_vec,S1params_vec,S,T,k",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " nonparametric Quantile of a positive  variable",
				" Vector of strikes 1",
				" Vector of Volatility 1",
				" S1 : indexBefore,indexAfter,flagBefore,flagAfter",
				" Forward",
				" Maturity"
				" Strike"
			
		},
		{
        		" Local_sabr2b_ImplicitVol",	/// name of the C++ function
                " RRRRRRRRRRR" ,				/// 11 parametres = 10 d'entree + 1 parametre de retour 
                " ARM_CF_SABR2B_ImplicitVol",
				" forward,strike,maturity,alpha,beta1,beta2,rho,nu,zero,lambda",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the SABR2B implicit vol of a call",
				" fwd value",
				" strike",
				" time to maturity (years)",
				" alpha",
				" beta1",
				" beta2",
				" rho",
				" nu",
				" zero",
				" lambda"
        },
		{
        		" Local_BiShiftedHeston_VanillaOption",	/// name of the C++ function
                " RRRRRRRRRRRRRRRRRRRRR" ,				/// 11 parametres = 10 d'entree + 1 parametre de retour 
                " ARM_CF_BiShiftedHeston_VanillaOption",
				" F1,V1,Vinfini1,lambda1,nu1,rho1,gamma1,F2,V2,Vinfini2,lambda2,nu2,rho2,gamma2,Correlations,k,T,CallPut,LambdaB,Flag",
                " 1",							/// visible in excel
                XLLOCALARM_CLOSEDFORME_GROUP,
                " ",
                " ",
                " gives the BiShiftedHeston VanillaOption",
				" forward 1",
				" InitialVariance 1",
				" InifiniteVariance 1",
				" Mean reversion 1",
				" Vol of Vol 1",
				" Correlation 1",
				" Shift 1",
				" forward 2",
				" InitialVariance 2",
				" InifiniteVariance 2",
				" Mean reversion 2",
				" Vol of Vol 2",
				" Correlation 2",
				" Shift 2",
				" {rhos,rhov,rho s1/v2,rho s2/v1}",
				" Strike",
				" Maturity",
				" (C)all / (P)ut",
				" LambdaB[def=0.1]",
				" Flag[=0]"	
        },
		{ 
			" Local_Merton_VanillaOption",
			" RRRRRRRRRRR",
			" ARM_CF_Merton_VanillaOption",
			" F, K, T, CallPut, Sigma, Lambda1, U1, Lambda2, U2, [N]",
			" 1",
			XLLOCALARM_CLOSEDFORME_GROUP,
			" "
		},
		{
			" Local_BSImpliedVol",
			" RRRRRRR",
			" ARM_CF_BSImpliedVol",
			" F, K, T, CallPut, Target, [BondPrice]",
			" 1",
			XLLOCALARM_CLOSEDFORME_GROUP,
			" "
		},
		{
			" Local_SuperNormal_Heston_VanillaOption",
			" RRRRRRRRRRRRRRRR",
			" ARM_CF_SuperNormalHeston_VanillaOption",
			" S0,K,T,CallPut,Rho1,Theta1,Kappa1,Nu1,V01,Rho2,Theta2,Kappa2,Nu2,V02,[Nb]",
			" 1",
			XLLOCALARM_CLOSEDFORME_GROUP,
			" ",
		},
			{
				" Local_MixteHeston_OptionPrice",
				" RRRRRRRRRRRRRRR",
				" ARM_CF_MixteHestonOptionPrice",
				" AsOfDate, ResetTime, Forward, Strike, CallOrPut, Sigma, V0, Kappa, Rho, Theta, VVol, [shift], [Times], [Level]",
				" 1",
				XLLOCALARM_CLOSEDFORME_GROUP,
				" ",
				" ",
				" returns option price with heston model (either with time dependent or constant parameters)",
				" ",
				" ",
				" ",
				" ",
				" ",
				" ",
				" Initial Variance",
				" Var Mean Reversion",
				" Forward / Variance correlation",
				" Long Term Variance",
				" Variance Volatility",
				" Relative Shift (default = 1)",
				" Time Pilars if time dependant level",
				" sigma (constant or time dependant)",
		},
			{
				" Local_Heston2B_OptionPrice",
				" RRRRRRRRRRRRRRRRRRR",
				" ARM_CF_Heston2BOptionPrice",
				" AsOfDate, ResetTime, Forward, Strike, CallOrPut, V01, Kappa1, Rho1, Theta1, VVol1, V02, Kappa2, Rho2, Theta2, VVol2,[shift], [Times], [Level]",
				" 1",
				XLLOCALARM_CLOSEDFORME_GROUP,
				" ",
				" ",
				" returns option price with heston model (either with time dependent or constant parameters)",
				" ",
				" ",
				" ",
				" ",
				" ",
				" ",
				" ",
				" Initial Variance",
				" Var Mean Reversion",
				" Forward / Variance correlation",
				" Long Term Variance",
				" Variance Volatility",
				" Relative Shift (default = 1)",
				" Time Pilars if time dependant level",
				" sigma (constant or time dependant)",
		},
