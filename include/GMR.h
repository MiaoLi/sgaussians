#ifndef __GMR_H__
#define __GMR_H__

#define MAX_GMM 30

struct DGMMStates {
	double *Mu;
	double **Sigma;
	double Prio;
};

struct DGMRStates{ 
	double **SigmaOI;
	double **SigmaIOIOInv;
	double **SigmaIIInv;
	double *muI;
	double *muO;
	double det; // for SigmaIOIOInv
};


// GMR Dynamics
class GMR
{
private:
	DGMRStates DGMRState[MAX_GMM];
	
	// for gauss
	double *g_diff, *g_diffp;

	// for regression
	double *r_diff, *r_diffa;
	double *x_s;
	double *h;
	double *Ctemp;
	double *Btemp;

	double **A,**D;
	double *B, *C;


	int nStates;
	int nVar;
	double delta_t;

	int in_dim;
	int out_dim;
	int *in_index, *out_index;	
	double *xti;	

	double *pdfx;

public:
	DGMMStates GMMState[MAX_GMM];

	GMR(int nStates, int nVar, char *fileName );
	GMR(int nStates, int nVar, char *f_mu, char *f_sigma, char *f_prio );	
	
	double GaussianPDF(int state, double *in_data);

	void initGMR(int in_dim, int out_dim, int *in_index, int *out_index);
	void regression(const Vector & indata, Vector & outdata, Matrix & derGMR);
	void regression(const Vector & indata, Vector & outdata);
};


#endif //__GMR_H__
