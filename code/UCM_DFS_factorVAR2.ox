#include <oxstd.h>
#include <oxdraw.h>
#include <oxprob.h>
#include <oxfloat.h>
#include <arma.h>
#import <modelbase>
#import <simula>
#import <maximize>
#import <database>
#import <maxsqp>
#include <packages/ssfpack/ssfpack_ex.h>

static decl s_cT, s_vY, s_mY, s_mX, s_cN;
static decl mT, mZ, mPhi, mOmega, mSigma, vPar;
static decl rho1, rho2, rho3, rho4, lambda1, lambda2, lambda3, lambda4, sigma2_mu, sigma2_gamma,
			sigma2_C1, sigma2_C2, sigma2_C3, sigma2_C4, sigma2_epsilon;

LoadData(const sFile, const ini_t, const length, const nA)
{
	decl mData, true_Y;
	mData = deletec(loadmat(sFile)');	// balance the panel
	true_Y = mData[0][ini_t+length:ini_t+length-1+nA];
	s_vY = mData[0][ini_t:ini_t+length-1];
	s_mX = mData[1:][ini_t:ini_t+length-1];
	s_cT = columns(s_vY);
	s_cN = rows(s_vY);
	return true_Y;
}

PrincipalCmp(const mX, const i_prop)		// mX must be demeaned first
{
	decl mc = variance(mX);					// also change to correlation matrix
	decl eigval, eigvec;
	eigensym(mc, &eigval, &eigvec);
	decl cumperc =cumulate(eigval' / sumr(eigval));
	decl NrCmp = sumc(cumperc .< i_prop);
	eigval = eigval[:NrCmp - 1];
	eigvec = eigvec[][:NrCmp - 1];
	decl pcaFactor = diag(eigval.^(-0.5)) * eigvec' * mX';
	decl pcaLambda = eigvec * diag(sqrt(eigval));
	return {pcaFactor, pcaLambda};
}

Antithetic(const nrAnti, const vDraw, const vHat, const dChi2)
{
	decl i, d, n = columns(vHat), err = vDraw, mret;

	for (mret = err, i=1; i<nrAnti; i++)
	{
		if (i==2)
		{
			d = sqrt(quanchi(1 - probchi(dChi2, n), n) / dChi2);
			err = vHat - d * (err - vHat);
		}
		else err = (2.0 * vHat) - err; 
		mret = mret | err;
	}
	return mret;
}

DFMdraw(const vR, const n_seq, const i_prop, const var_or)
{
	ranseed(5);
	decl pc_factor, pc_lambda, tmp_mX, tmp_mu, i_n, i_p, i_T;
	tmp_mX = vR | s_mX;
	tmp_mu = meanr(tmp_mX);
	tmp_mX = tmp_mX - tmp_mu * ones(1, columns(tmp_mX));
	[pc_factor, pc_lambda] = PrincipalCmp(tmp_mX', i_prop);
	i_n = rows(pc_lambda);
	i_p = columns(pc_lambda);
	
	decl mT_Doz, mZ_Doz, mPhi_Doz, mH_Doz, mG_Doz, mOmega_Doz, mSigma_Doz, mDelta_Doz;
	mZ_Doz = pc_lambda ~ zeros(i_n, (var_or - 1) * i_p);
	mG_Doz = diagonalize(variance(tmp_mX') - pc_lambda * pc_lambda');
	decl var_Y, var_X;
	var_Y = pc_factor[][var_or:];
	i_T = columns(var_Y);
	var_X = pc_factor[][0 : i_T - 1];
	for (decl i=1; i<var_or; i++)
	{
		var_X = var_X | pc_factor[][i : i_T - 1 + i];
	}
	olsr(var_Y, var_X, &mT_Doz);
	mH_Doz = diagcat(variance((var_Y - mT_Doz * var_X)'), zeros((var_or-1)*i_p, (var_or-1)*i_p));
	mOmega_Doz = diagcat(mH_Doz, mG_Doz);
	mT_Doz = mT_Doz | ( unit((var_or - 1) * i_p) ~ zeros((var_or - 1) * i_p, i_p) );
	mSigma_Doz = reshape((unit((var_or*i_p)^2) - mT_Doz**mT_Doz)^(-1) * vec(mH_Doz),
				 var_or*i_p, var_or*i_p) | zeros(1, var_or*i_p);
	mDelta_Doz = zeros(var_or * i_p, 1) | tmp_mu;
	mPhi_Doz = mT_Doz | mZ_Doz;
	
	decl bs_R = constant(M_NAN, n_seq, columns(vR));
	tmp_mX = tmp_mX + tmp_mu * ones(1, columns(tmp_mX));
	decl s_mGamma, mkf, mwgt, veps_mu, vepssim, err;
	veps_mu = SsfCondDensEx(DS_SMO, tmp_mX, mPhi_Doz, mOmega_Doz, mSigma_Doz, mDelta_Doz)[var_or * i_p][];
	parallel for (decl i=0; i<n_seq; i++)
	{
		if (var_or == 1)
		{
			bs_R[i][] = SsfCondDens(ST_SIM, tmp_mX, mPhi_Doz, mOmega_Doz, mSigma_Doz, mDelta_Doz)[i_p][];
		}
		else	// when var_or > 1, mOmega_Doz is not full p.d., when var_or = 1, the two algorithms are equivalent, but the first faster
		{
			s_mGamma = diag(zeros(1, var_or * i_p) ~ 1 ~ zeros(1, i_n - 1));
			mkf = KalmanFil(tmp_mX, mPhi_Doz, mOmega_Doz, mSigma_Doz, mDelta_Doz);
			mwgt = SimSmoWgt(s_mGamma, mkf, mPhi_Doz, mOmega_Doz, mSigma_Doz, mDelta_Doz);
			err = rann(1, columns(tmp_mX));
			vepssim = SimSmoDraw(s_mGamma, err, mwgt, mkf, mPhi_Doz, mOmega_Doz, mSigma_Doz, mDelta_Doz)[var_or * i_p][1:];
			vepssim = Antithetic(4, vepssim, veps_mu, err*err');
			bs_R[i][] = vR - vepssim;
		}
	}
	return vR | bs_R;
}

Synthesize(const Fil_vY, const VR, const n_sam, const iprop, const varor)
{
	decl BS_R = DFMdraw(VR, n_sam, iprop, varor);
	return ones(rows(BS_R), 1) * Fil_vY + BS_R;
}

SetSsf(const vP){
	rho1 = exp(vP[0]) ./ (1 + exp(vP[0]));
	rho2 = exp(vP[1]) ./ (1 + exp(vP[1]));
	rho3 = exp(vP[2]) ./ (1 + exp(vP[2]));
	lambda1 = M_2PI ./ (11 + 12 * exp(vP[3]) / (1 + exp(vP[3])));
	lambda2 = M_2PI ./ (23 + 17 * exp(vP[4]) / (1 + exp(vP[4])));
	lambda3 = M_2PI ./ (40 + 30 * exp(vP[5]) / (1 + exp(vP[5])));  // restrict periods according to the paper
	sigma2_mu = exp(vP[6]);
	sigma2_gamma = exp(vP[7]);
	sigma2_C1 = exp(vP[8]) * (1-rho1^2);
	sigma2_C2 = exp(vP[9]) * (1-rho2^2);
	sigma2_C3 = exp(vP[10]) * (1-rho3^2);
	sigma2_epsilon = exp(vP[11]);
	
	mZ = 1 ~ 1 ~ zeros(1, 10) ~ <1, 0, 1, 0, 1, 0>;
	mT = (1 ~ zeros(1, 17)) |
		 (0 ~ -1*ones(1, 11) ~ zeros(1, 6)) |
		 (zeros(10, 1) ~ unit(10) ~ zeros(10, 7)) |
		 (zeros(2, 12) ~ rho1 * ((cos(lambda1) ~ sin(lambda1)) | (-sin(lambda1) ~ cos(lambda1))) ~ zeros(2, 4)) |
		 (zeros(2, 14) ~ rho2 * ((cos(lambda2) ~ sin(lambda2)) | (-sin(lambda2) ~ cos(lambda2))) ~ zeros(2, 2)) |
		 (zeros(2, 16) ~ rho3 * ((cos(lambda3) ~ sin(lambda3)) | (-sin(lambda3) ~ cos(lambda3))));
	mPhi = mT | mZ;
	mOmega = diagcat(diag(sigma2_mu~sigma2_gamma~zeros(1, 10)~sigma2_C1~sigma2_C1~sigma2_C2~sigma2_C2~sigma2_C3~sigma2_C3), sigma2_epsilon);
	mSigma = diagcat(diag(-ones(1, 12)), diag(exp(vP[<8,8,9,9,10,10>]))) | zeros(1, 18);
}

LogLikelihood(const vP, const adLogLik, const avSco, const amHes)
{
	SetSsf(vP);
	decl dVar=1, Likp;
	SsfLikEx(&Likp, &dVar, s_vY, mPhi, mOmega, mSigma);
	adLogLik[0] = Likp / columns(s_vY); 
	return 1; 
}

MaxLik()
{
    decl ir, dloglikmn;
	MaxControlEps(1e-4, 1e-4);
    MaxControl(1000, 50, 1);
    ir = MaxBFGS(LogLikelihood, &vPar, &dloglikmn, 0, 1);
	SetSsf(vPar);
}

main()
{
	// Augmented Stock and Watson model
	// moving window estimate and forecasting evaluation
	decl dtime, smodeltype, total_T, s_file, vEN;
	dtime = 0;
	smodeltype = "UCM_DFS_factorVAR2";
	s_file = "Book1982_2016_autosim1.xls";
	vEN = deletec(loadmat(s_file)')[0][];
	total_T = columns(vEN);	// total data length

	// general model specification
	decl i_setlength, s_nA, s_nVar, i_prop, RMSE, MAFE;
	i_setlength = 275;					// size of the training set, 275 leads to around 100 RMSE's for DM test
	s_nA = 30;							// max step-ahead forecast
	s_nVar = 2;							// VAR(s_nVar) for state dynamics
	i_prop = 0.95;						// % of variation explained by PC's
	vPar = zeros(12, 1);				// initial parameters of UC model

	decl predErr = constant(M_NAN, total_T - i_setlength - s_nA + 1, s_nA);
	for (decl i_t0 = 0; i_t0 < total_T - i_setlength - s_nA + 1; i_t0++)	// the initial data point (in moving window analysis)
	{
		decl true_vY = LoadData(s_file, i_t0, i_setlength, s_nA);

		// first estimation
		MaxLik();

		// filtering / obtaining residuals
		decl minf, mKF, fil_vY, vRes;
		minf = KalmanInit(s_vY, mPhi, mOmega, mSigma);						// you are not supposed to understand this
		mKF = KalmanFilEx(minf, s_vY, mPhi, mOmega, mSigma);
		vRes = mKF[0][];
		fil_vY = s_vY - vRes;
		
		// Doz method and synthesize
		decl n_sample;
		n_sample = 199;														// number of samples of the simulation smoother				
		s_mY = Synthesize(fil_vY, vRes, n_sample, i_prop, s_nVar);

		// UCmodel estimation and forecast
		decl EN_pred;
		EN_pred = constant(M_NAN, rows(s_mY), s_nA);
		for (decl i=0; i<rows(s_mY); i++)
		{
			s_vY = s_mY[i][];
			MaxLik();
			EN_pred[i][] = SsfCondDensEx(ST_SMO, s_vY ~ constant(M_NAN, 1, s_nA),
								    	 mPhi, mOmega, mSigma)[18][s_cT : s_cT + s_nA - 1];
		}
		predErr[i_t0][] = meanc(EN_pred) - true_vY;
		println(sprint(i_t0+1)+"th window is finished.", predErr[i_t0][]);
	}
	println("step "+sprint(s_nA)+" is finished.");
//	savemat("UCM_DFS_factorVAR2_predErr of "+sprint(s_nA)+"step-ahead forecast.xlsx", predErr);
	RMSE = sqrt(meanc(predErr.^2))';
	MAFE = meanc(fabs(predErr))';

	println("\nForecasting evaluation for "+sprint(s_nA)+"step-ahead under model ", smodeltype, ":", "%cf",
	 {"%10.0f", "%10.4f", "%10.4f"}, "%c", {"steps", "RMSE", "MAFE"}, range(1, s_nA)' ~ RMSE ~ MAFE);
	 
}
