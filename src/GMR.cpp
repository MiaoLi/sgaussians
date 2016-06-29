#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "GMR.h"
#include "MathLib.h"

GMR::GMR(int nStates, int nVar, char *fileName)
{
	char f_mu[256], f_sigma[256], f_prio[256];

	sprintf(f_mu   , "%s%s", fileName, "_mu.txt"   );
	sprintf(f_sigma, "%s%s", fileName, "_sigma.txt");
	sprintf(f_prio , "%s%s", fileName, "_prio.txt" );

	GMR(nStates, nVar, f_mu, f_sigma, f_prio);
}

GMR::GMR(int nStates, int nVar, char *f_mu, char *f_sigma, char *f_prio )
{
	FILE *fid;
	int s, i, j;

	this->nStates = nStates;
	this->nVar    = nVar;

	//GMMState  = (DGMMStates *)malloc(nStates*sizeof(DGMMStates ) );
	//DGMRState = (DGMRStates *)malloc(nStates*sizeof(DGMRStates ) );

	for( s=0; s<nStates; s++ ){
		GMMState[s].Mu       = svector(nVar);
		GMMState[s].Sigma    = smatrix(nVar, nVar );
	}
	pdfx = svector(nStates);

	// f_mu	
	fid = fopen(f_mu, "r+" );
	for( i=0; i<nVar; i++ ){
		for( s=0; s<nStates; s++ ){
			fscanf(fid, "%lf", &(GMMState[s].Mu[i]) );
		}
	}
	fclose(fid);

	// f_sigma
    fid = fopen(f_sigma, "r+" );
	for( s=0; s<nStates; s++ ){
		for( i=0; i<nVar; i++ ){
			for( j=0; j<nVar; j++ ){
				fscanf(fid, "%lf", &(GMMState[s].Sigma[i][j]) );
			}
		}
	}
	fclose(fid);

	// f_prio
    fid = fopen(f_prio, "r+" );
	for( s=0; s<nStates; s++ ){
		fscanf(fid, "%lf", &(GMMState[s].Prio ) );
	}
	fclose(fid);

}

void GMR::initGMR(int in_dim, int out_dim, int *in_index, int *out_index)
{
	int i,j,s;
	int dim;
	double **covIOIO, **covII;
	double temp;

	dim = in_dim+out_dim;
	this->in_dim  = in_dim;
	this->out_dim  = out_dim;
	this->in_index  = (int *)malloc(in_dim*sizeof(int));
	this->out_index = (int *)malloc(out_dim*sizeof(int));
	xti       = svector(in_dim);
	covIOIO   = smatrix(dim, dim);
	covII     = smatrix(in_dim, in_dim);

	// for gauss
	g_diff  = svector(in_dim);
	g_diffp = svector(in_dim);

	// for regression
	h        = svector(nStates);
	x_s      = svector(in_dim);
	r_diff   = svector(in_dim);
	r_diffa  = svector(in_dim);
	Btemp    = svector(in_dim);
	Ctemp    = svector(in_dim);

	A = smatrix(dim, dim);
	D = smatrix(dim, dim);

	B = svector(dim);
	C = svector(dim);

	for( i=0; i<in_dim ; i++) this->in_index[i]   = in_index[i];
	for( i=0; i<out_dim; i++) this->out_index[i]  = out_index[i];

	// Initialize GMR Variables
	for( s=0; s<nStates; s++ ){
		DGMRState[s].SigmaOI = smatrix(in_dim, in_dim );
		DGMRState[s].SigmaIOIOInv = smatrix(dim, dim );
		DGMRState[s].SigmaIIInv = smatrix(in_dim, in_dim );
		DGMRState[s].muI = svector(in_dim);
		DGMRState[s].muO = svector(in_dim);
	}

	for( s=0; s<nStates; s++ ){
		for( i=0; i<in_dim; i++){
			DGMRState[s].muI[i] = GMMState[s].Mu[in_index[i] ];

			for( j=0; j<in_dim; j++){
				covII[i][j]                  = GMMState[s].Sigma[in_index[ i]][in_index[ j]];
				covIOIO[i][j]                = GMMState[s].Sigma[in_index[ i]][in_index[ j]];
			}
			for( j=0; j<out_dim; j++){
				covIOIO[i][in_dim+j]         = GMMState[s].Sigma[in_index[ i]][out_index[j]];
			}
		}
		for( i=0; i<out_dim; i++){
			DGMRState[s].muO[i] = GMMState[s].Mu[out_index[i]];

			for( j=0; j<in_dim; j++){
				DGMRState[s].SigmaOI[i][j]  = GMMState[s].Sigma[out_index[i]][in_index[ j]];	
				covIOIO[i][in_dim+j]        = GMMState[s].Sigma[out_index[i]][in_index[ j]];
			}
			for( j=0; j<out_dim; j++){
				covIOIO[in_dim+i][in_dim+j] = GMMState[s].Sigma[out_index[i]][out_index[j]];
			}
		}


		matsvdinvall( covIOIO, dim   , dim   , 0.00001, &temp               , DGMRState[s].SigmaIOIOInv );
		matsvdinvall( covII  , in_dim, in_dim, 0.00001, &(DGMRState[s].det) , DGMRState[s].SigmaIIInv   );
	}	
}


double GMR::GaussianPDF(int state, double *in_data)
{
	double p;

	matsub_c1(in_dim, in_data, DGMRState[state].muI, g_diff );
	matmul_c(in_dim, in_dim, DGMRState[state].SigmaIIInv, g_diff, g_diffp );
	p = in_prod(in_dim, g_diff, g_diffp );
	p = exp(-0.5*p) / sqrt(pow(2.0*3.14159,in_dim)*(abs(DGMRState[state].det)+1e-30));

    if(p < 1e-30){
		return 1e-30;
    }
	else{
		return p;
	}
}


// x, xdot : meter
// GMM is trained in mm
void GMR::regression(double *indata, double *outdata, double **derGMR)
{
	int s, i, j;	
	double pa;
	double temp;

    //for( i=0; i<in_dim; i++ ) x_s[i] = indata[i]*1000.0;
	for( i=0; i<in_dim; i++ ) x_s[i] = indata[i];
    
	
	for( s=0; s<nStates; s++ ){		
		pdfx[s] = GMMState[s].Prio*GaussianPDF(s, x_s);
	}
	pa = vsum(nStates, pdfx );

	matzero_c1( in_dim, outdata );
	matzero_c1( in_dim, Ctemp );
	for( s=0; s<nStates; s++ ){		
		h[s] = pdfx[s]/(pa + 1e-30 );
		matsub_c1(in_dim, x_s, DGMRState[s].muI, r_diff );
		matmul_c(in_dim, in_dim, DGMRState[s].SigmaIIInv, r_diff , r_diffa );
		matmul_c(in_dim, in_dim, DGMRState[s].SigmaOI   , r_diffa, r_diff  );

		for(i=0; i<in_dim; i++ ){
			outdata[i] += h[s]*(r_diff[i] + DGMRState[s].muO[i]);
			Ctemp[i]   += h[s]*r_diffa[i];
		}		
	}
	//matsmul_c1(in_dim, 0.001, outdata );

	for(i=0; i<in_dim; i++) for(j=0; j<in_dim; j++) derGMR[i][j] = 0.0;

	for( s=0; s<nStates; s++ ){		
		//A{k} = Sigma(out,in,k)*invSigma{k};
		matmul_cc(in_dim, in_dim, in_dim, DGMRState[s].SigmaOI, DGMRState[s].SigmaIIInv, A );
		
		//B{k} = Mu(out,k) - A{k}*Mu(in,k);
		matmul_c(in_dim, in_dim, A, DGMRState[s].muI, Btemp );
		matsub_c1(in_dim, DGMRState[s].muO, Btemp, B );
	
		//C{k} = ( -invSigma{k} * (x(:)-Mu(in,k)) + Ctemp)';
		matsub_c1(in_dim, x_s, DGMRState[s].muI, r_diff );
		matmul_c(in_dim, in_dim, DGMRState[s].SigmaIIInv, r_diff, r_diffa );
		matsub_c1(in_dim, Ctemp, r_diffa, C);

		temp = 0.0;
		for( i=0; i<in_dim; i++){
			//D[j](j,:) = ( C{k}*(A{k}(j,:)*x(:) + B{k}(j)) +A{k}(j,:));
			for( j=0; j<in_dim; j++) temp += A[i][j]*x_s[j];
			temp+= B[i];
			for( j=0; j<in_dim; j++) D[i][j] = h[s]*(C[j]*temp +A[i][j]);
		}
		for(i=0; i<in_dim; i++) for(j=0; j<in_dim; j++) derGMR[i][j] += D[i][j];
	}

}

void GMR::regression(double *indata, double *outdata)
{
	int s, i, j;	
	double pa;
	double temp;

    //for( i=0; i<in_dim; i++ ) x_s[i] = indata[i]*1000.0;
	for( i=0; i<in_dim; i++ ) x_s[i] = indata[i];
    
	
	for( s=0; s<nStates; s++ ){		
		pdfx[s] = GMMState[s].Prio*GaussianPDF(s, x_s);
	}
	pa = vsum(nStates, pdfx );

	matzero_c1( out_dim, outdata );
	for( s=0; s<nStates; s++ ){		
		h[s] = pdfx[s]/(pa + 1e-30 );
		matsub_c1(in_dim, x_s, DGMRState[s].muI, r_diff );
		matmul_c(in_dim, in_dim, DGMRState[s].SigmaIIInv, r_diff , r_diffa );
		matmul_c(out_dim, in_dim, DGMRState[s].SigmaOI  , r_diffa, r_diff  );

		for(i=0; i<out_dim; i++ ){
			outdata[i] += h[s]*(r_diff[i] + DGMRState[s].muO[i]);
		}		
	}

}
