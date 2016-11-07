#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <cmath>
#include <vector>
#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>
#include "bramauxiliary.h"


//#define NDEBUG


using namespace std;

// random number generator 
// see http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Generation.html#Random-Number-Generation 
gsl_rng_type const * T; // gnu scientific library rng type
gsl_rng *r; // gnu scientific rng 

int seed = -1;

// give the outputfile a unique name
// by using the create_filename function (see 
// "bramauxiliary.h")
string filename("iter_m");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  // output file 

// parameters
int skip = 100;
double tau=0.25;       // fraction of a generation between critical period of development and selection 
double rho=0.50;		// environmental auto-correlation
  
double xi_bar=0;		// mean of stationary autocorrelated noise
double sig_xi=0.001;		// standard deviation of stationary autocorrelated noise
    
double omega2=40;		// width of fitness function, 3=> strong stabilizing selection
double sigma_e2=1;    	// variance in residual component of phenotypic variation
double A=0;			// intercept (elevation) of evolution in reference environment
double B=2;			// slope of reaction norm in reference environment
double Wmax=1;			// maximal fitness
double Gaa=0.1;		// variance of elevation
double Gbb=0.045;		// variance of plasticity
double Gmm=0.045;        // variance of maternal effect
double omega_b_2=100;
double omega_m_2=100;
double freq=0.5;
double Ut=10;			// time at which step change of size delta occurs
double epsinit=0;      // initial value of environment
double delta=10;		// size of step change/amp of sine wave
double init_g = 0;
double init_b = 0;
double init_m = 0;
//int nrp=1000000-1;		// length of simulation
int nrp=10-1;		// length of simulation


// G-matrix
double GG[] = { 0, 0, 0,
                0, 0, 0,
                0, 0, 0 };

gsl_matrix_view GGx;

bool shift = true;

vector<double> abar(nrp,0);
vector<double> bbar(nrp,0);
vector<double> mbar(nrp,0);
vector<double> zbar(nrp,0);
vector<double> zsbar(nrp,0);
vector<double> gaz(nrp,0);
vector<double> gbz(nrp,0);
vector<double> gmzstar(nrp,0);
vector<double> sigma_z2(nrp,0);
vector<double> fitness(nrp,0);
vector<double> fitvar(nrp,0);
vector<double> fitexp(nrp,0);
vector<double> epst(nrp,0);
vector<double> epstT(nrp,0);
vector<double> sigz(nrp,0);
vector<double> xi_devel(nrp,0);
vector<double> xi_select(nrp,0);

void init(int argc, char **argv)
{
    DataFile <<  
            "t;" << 
            "abar;" << 
            "bbar;" << 
            "mbar;" << 
            "zbar;" << 
            "zsbar;" << 
            "gaz;" << 
            "gbz;" << 
            "gmzstar;" << 
            "sigma_z2;" << 
            "fitness;" << 
            "fitvar;" << 
            "fitexp;" << 
            "epst;" << 
            "epstT;" << 
            "sigz;" << 
            "xi_devel;" << 
            "xi_select;" << endl;


    // obtain a seed from current nanosecond count
	seed = get_nanoseconds();

    // set up the random number generators
    // (from the gnu gsl library)
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);

    // initialize values from command line
    Gaa = atof(argv[1]);
    Gbb = atof(argv[2]);
    Gmm = atof(argv[3]);

    omega_b_2 = atof(argv[4]);
    omega_m_2 = atof(argv[5]);

    shift = atoi(argv[6]);
    freq = atof(argv[7]);
    tau = atof(argv[8]);
    delta = atof(argv[9]);
    init_g = atof(argv[10]);
    init_b = atof(argv[11]);
    init_m = atof(argv[12]);

    //initialize trait values
    mbar[0] = init_m;
    abar[0] = init_g;
    bbar[0] = init_b;
    epstT[0] = epsinit;
    epst[0] = epsinit;

    sigma_z2[0] = (Gaa+Gbb*pow(epsinit,2.0)+Gmm*pow(epsinit,2.0)+sigma_e2);

    GG[0] = Gaa;
    GG[4] = Gbb;
    GG[8] = Gmm;

    GGx = gsl_matrix_view_array(GG, 3, 3);
}


void generate_envt()
{
    xi_devel[0] = gsl_ran_gaussian(r, sig_xi);

    double tmp, j;

    for (int k = 0; k < nrp; ++k)
    {
        if (rho == 0.0)
        {
            xi_select[k] = gsl_ran_gaussian(r, sig_xi);
            xi_devel[k] = gsl_ran_gaussian(r, sig_xi);
        }
        else
        {
            tmp=(rho*xi_devel[k]) + (sig_xi*sqrt(1.0-(pow(rho,2)))*gsl_ran_gaussian(r, 1.0));
            xi_select[k]=tmp;
            j=1.0/tau;

            if (j >= 2.0)
            {
                for (int i = 0; i < int(j-1);++i)
                {
                    tmp=(rho*tmp) + (sig_xi*sqrt(1-(pow(rho,2)))*gsl_ran_gaussian(r,1.0));
                }
                xi_devel[k+1]=(rho*tmp) + (sig_xi*sqrt(1-(pow(rho,2)))*gsl_ran_gaussian(r,1.0));
            }
        }
    }
}

int main(int argc, char **argv)
{
    init(argc, argv);

    double eps_t, eps_t1, eps_tT, gamma_z, gamma_b, gamma_m;
    double theta, tmp, dsigm, dsiga;

    gsl_vector *mat1 = gsl_vector_alloc(3);
    gsl_vector *mat2 = gsl_vector_alloc(3);
    gsl_vector *mat3 = gsl_vector_alloc(3);
    gsl_vector *beta = gsl_vector_alloc(3);
    gsl_vector *del = gsl_vector_alloc(3);

    for (int t = 0; t < nrp; ++t)
    {
        // obtain environments
        eps_t = shift ? epsinit+(t>Ut)*delta : sin(freq*t); 

        epst[t] = eps_t;
        eps_t = epst[t] + xi_select[t];
        eps_t1 = epst[t-1] + xi_select[t-1];

        eps_tT= shift ? epsinit+((t-tau)>Ut)*delta : sin(freq*(t-tau));
        epstT[t]=eps_tT;
        eps_tT=epstT[t] + xi_devel[t];


        // calculate sigma_z^2
        sigma_z2[t]=Gaa+Gbb*pow(eps_tT,2.0)+Gmm*(Gaa+pow(abar[t],2.0))+Gaa*mbar[t]*(1+mbar[t])+sigma_e2;
        sigma_z2[t]=0 > sigma_z2[t] ? 0 : sigma_z2[t];

        // calculate gamma
        gamma_z = 1.0/(omega2 + sigma_z2[t]);
        gamma_b=1.0/(omega_b_2+Gbb);
        gamma_m=1.0/(omega_m_2+Gmm);

        // calculate phenotypic change
        zbar[t] = abar[t] + bbar[t]*eps_tT +mbar[t]*abar[t];

        theta = A+B*eps_t;
        tmp=zbar[t]-theta;

        // set first part of the selection gradients
        gsl_vector_set(mat1, 0, -1.0/omega2*tmp*(1+mbar[t]));
        gsl_vector_set(mat1, 1,  -1.0/omega2*tmp*eps_tT);
        gsl_vector_set(mat1, 2,  -1.0/omega2*tmp*abar[t]);
        
        dsigm=Gaa*(1+2*mbar[t]);

        dsiga=2*Gmm*abar[t];

        gsl_vector_set(mat2, 0, 0.5*(-1.0/omega2)*dsiga);
        gsl_vector_set(mat2, 1, 0);
        gsl_vector_set(mat2, 2, 0.5*(-1.0/omega2)*dsigm);

        gsl_vector_set(mat3, 0, 0);
        gsl_vector_set(mat3, 1, -1.0/omega2 * omega2*bbar[t]/omega_b_2);
        gsl_vector_set(mat3, 2, -1.0/omega2 * omega2*mbar[t]/omega_m_2);

        gsl_vector_set_zero(beta);
        gsl_vector_set_zero(del);

        gsl_vector_add(beta, mat1);
        gsl_vector_add(beta, mat2);
        gsl_vector_add(beta, mat3);
        
        cout << gsl_vector_get(beta, 0) << "; " << gsl_vector_get(beta, 1) << ";" << gsl_vector_get(beta, 1) << endl;
   
        // see http://www.netlib.org/lapack/explore-html/d7/d15/group__double__blas__level2_gadd421a107a488d524859b4a64c1901a9.html#gadd421a107a488d524859b4a64c1901a9 
        gsl_blas_dgemv(CblasTrans, 1.0, &GGx.matrix, beta, 0.0, del);

        // calculate mean fitness
        fitness[t] = Wmax*sqrt(gamma_z*gamma_b*gamma_m*omega2*omega_b_2*omega_m_2)*exp(-(gamma_z/2)*(pow(zbar[t]-theta,2))- .5 * gamma_m * (mbar[t]*mbar[t]) - .5 * gamma_b * (bbar[t]*bbar[t]));
        fitvar[t] = sqrt(gamma_z*gamma_b*gamma_m*omega2*omega_b_2*omega_m_2);
        fitexp[t] = exp(-(gamma_z/2)*(pow(zbar[t]-theta,2.0))- .5 * gamma_m * (mbar[t]*mbar[t]) - .5 * gamma_b * (bbar[t]*bbar[t]));		 

        // update elevation & plasticity
        abar[t+1] = abar[t] + gsl_vector_get(del, 0);
        bbar[t+1] = bbar[t] + gsl_vector_get(del, 1);
        mbar[t+1] = mbar[t] + gsl_vector_get(del, 2);

        if (t % skip == 0)
        {
            DataFile << t << ";" <<
                        abar[t] << ";" << 
                        bbar[t] << ";" << 
                        mbar[t] << ";" << 
                        zbar[t] << ";" << 
                        zsbar[t] << ";" << 
                        gaz[t] << ";" << 
                        gbz[t] << ";" << 
                        gmzstar[t] << ";" << 
                        sigma_z2[t] << ";" << 
                        fitness[t] << ";" << 
                        fitvar[t] << ";" << 
                        fitexp[t] << ";" << 
                        epst[t] << ";" << 
                        epstT[t] << ";" << 
                        sigz[t] << ";" << 
                        xi_devel[t] << ";" << 
                        xi_select[t] << ";" << endl;
        }
    }
}


