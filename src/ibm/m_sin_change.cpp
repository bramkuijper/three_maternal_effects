// the evolution of maternal effects in a sinusoidal environment
// 
// Bram Kuijper & Rebecca B. Hoyle
//
// this code is published according to the GNU Public license v3
// https://www.gnu.org/licenses/gpl.html 
//
// When using/extending this simulation, please cite the corresponding paper(s):
// Kuijper, B & Hoyle, R.B. When to rely on maternal effects 
// as opposed to phenotypic plasticity? 
//
// You may also find the following paper interesting:
// Kuijper, B.; Johnstone, R. A. & Townley, S. (2014). The evolution of 
// multivariate maternal effects. PLoS Comp. Biol. 10: e1003550. 
// http://dx.doi.org/10.1371/journal.pcbi.1003550
//  


#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cassert>

// random number generation
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// various functions, such as unique filename creation
#include "bramauxiliary.h"

//#define NDEBUG
//
// the compilation sign should only be turned on when one wants
// to assess the complete distribution of phenotypes
// 
//#define DISTRIBUTION

using namespace std;

// number of generations
const int NumGen = 50000;

// population size
const int Npop = 5000; 

// number of generations to skip when outputting data
const int skip = 100;

// track number of survivors
int NSurv = 0;

// indicator variable if we are printing stats for this generation
bool do_stats = 0;

double epsilon = 0; // the value of the environment
double epsilon_sens = 0; // the perceived value of the environment


double theta = 0;  // phenotypic optimum
double omega2 = 0; // width of the selection function
double omega_b_2 = 0; // width of the selection function for plasticity
double omega_m_m_2 = 0; // width of the selection function for maternal effects
double omega_m_e_2 = 0; // width of the selection function for maternal effects
double omega_m_g_2 = 0; // width of the selection function for maternal effects
double wmin = 0.0; // minimal survival probability
double sigma_e = 1.0; // variance of developmental noise
double sigma_ksi = 0.1; // variance of the autocorrelated process
double rho_t = 0.5; // temporal autocorrelation
double mu_g 	  = 0.05;            // mutation rate
double sdmu         = 0.05;			 // standard deviation mutation size
double mu_m_m 	  = 0.05;            // mutation rate
double mu_m_g 	  = 0.05;            // mutation rate
double mu_m_e 	  = 0.05;            // mutation rate
double mu_b 	  = 0.05;            // mutation rate
double ksi = 0;			 // standard deviation mutation size
double tau = 0.0;       // developmental time lag

double rate = 0.0;
double intercept = 0.0;
double ampl = 0.0;

// initial values
double int_t0 = 0;
double rate_t0 = 0;
double ampl_t0 = 0;

// values after perturbation
double intptb = 0;
double rateptb = 0;
double amplptb = 0;

bool envt_is_changed = false; // whether the environment is changed or not

const int n_alleles_b = 2; // number of alleles underlying genetic architecture
const int n_alleles_g = 50; // number of alleles underlying genetic architecture
const int n_alleles_m = 2; // number of alleles underlying genetic architecture

// keep track of the current generation number
int generation = 0;
int t_change = 0;

// random seed
unsigned seed = 0;

// gnu scientific library random number generator initialization
// http://www.gnu.org/software/gsl/ 
gsl_rng_type const * T; // gnu scientific library rng type
gsl_rng *r; // gnu scientific rng 

// the individual struct
struct Individual
{
    double g[n_alleles_g]; 
    double b[n_alleles_b]; 
    double m_g[n_alleles_m]; // maternal genetic effect
    double m_e[n_alleles_m]; // maternal environmental effect
    double m_m[n_alleles_m]; // maternal phenotypic effect (cascading etc)
    double phen; // an individual's phenotype, z
    double phen_m_m; //m
    double phen_m_g; //m
    double phen_m_e; //m
    double phen_b; // b
    double phen_g; // g

    double envt;
};

// allocate a population and a population of survivors
typedef Individual Population[Npop];
Population Pop;
Population Survivors;

// generate a unique filename for the output file
string filename("sim_evolving_m");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  

#ifdef DISTRIBUTION
// generate a filename for the phenotype distribution file
string filename_new2(create_filename("sim_evolving_m_dist"));
ofstream distfile(filename_new2.c_str());
#endif //DISTRIBUTION

// initialize simulations from command line arguments
void initArguments(int argc, char *argv[])
{
	mu_g = atof(argv[1]);
	mu_m_g = atof(argv[2]);
	mu_m_e = atof(argv[3]);
	mu_m_m = atof(argv[4]);
	mu_b = atof(argv[5]);
	sdmu = atof(argv[6]);
	sigma_e = atof(argv[7]);
	sigma_ksi = atof(argv[8]);
	wmin = atof(argv[9]);
	rho_t = atof(argv[10]);
	omega2 = atof(argv[11]);
	omega_b_2 = atof(argv[12]);
	omega_m_m_2 = atof(argv[13]);
	omega_m_e_2 = atof(argv[14]);
	omega_m_g_2 = atof(argv[15]);
	tau = atof(argv[16]);

    intercept = int_t0 = atof(argv[17]);
    rate = rate_t0 = atof(argv[18]);
    ampl = ampl_t0 = atof(argv[19]);

    intptb = atof(argv[20]);
    rateptb = atof(argv[21]);
    amplptb = atof(argv[22]);

    t_change = atoi(argv[23]);
}


// mutation according to a continuum of alleles model
void MutateG(double &G)
{
	G += gsl_rng_uniform(r)<mu_g ? gsl_ran_gaussian(r, sdmu) : 0;
}

void MutateM_M(double &G)
{
	G += gsl_rng_uniform(r)<mu_m_m ? gsl_ran_gaussian(r,sdmu) : 0;
}

void MutateM_E(double &G)
{
	G += gsl_rng_uniform(r)<mu_m_e ? gsl_ran_gaussian(r,sdmu) : 0;
}

void MutateM_G(double &G)
{
	G += gsl_rng_uniform(r)<mu_m_g ? gsl_ran_gaussian(r,sdmu) : 0;
}

void MutateB(double &G)
{
	G += gsl_rng_uniform(r)<mu_b ? gsl_ran_gaussian(r,sdmu) : 0;
}

// write the parameters (typically at the end of the output file)
void WriteParameters()
{
	DataFile << endl
		<< endl
		<< "type;" << "three_types" << ";" << endl
        << "mu_g;" << mu_g << ";" << endl
        << "mu_m_m;" << mu_m_m << ";" << endl
        << "mu_m_e;" << mu_m_e << ";" << endl
        << "mu_m_g;" << mu_m_g << ";" << endl
        << "mu_b;" << mu_b << ";" << endl
        << "sdmu;" << sdmu << ";" << endl
        << "omega2;" << omega2 << ";" << endl
        << "omega_b_2;" << omega_b_2 << ";" << endl
        << "omega_m_m_2;" << omega_m_m_2 << ";" << endl
        << "omega_m_g_2;" << omega_m_g_2 << ";" << endl
        << "omega_m_e_2;" << omega_m_e_2 << ";" << endl
        << "wmin;" << wmin << ";" << endl
        << "int_t0;" << int_t0 << ";" << endl
        << "rate_t0;" << rate_t0 << ";" << endl
        << "ampl_t0;" << ampl_t0 << ";" << endl
        << "intptb;" << intptb << ";" << endl
        << "rateptb;" << rateptb << ";" << endl
        << "amplptb;" << amplptb << ";" << endl
        << "sigma_e;" << sigma_e << ";" << endl
        << "sigma_ksi;" << sigma_ksi << ";" << endl
        << "rho_t;" << rho_t << ";" << endl
        << "tau;" << tau << ";" << endl
		<< "seed;" << seed << ";"<< endl;
}

// initialize the simulation
// by giving all the individuals 
// genotypic values
//
// and doing some other stuff (e.g., random seed)
void Init()
{
    // get the timestamp (with nanosecs)
    // to initialize the seed
	seed = get_nanoseconds();
    
    // set up the random number generators
    // (from the gnu gsl library)
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);


	// initialize the whole populatin
	for (int i = 0; i < Npop; ++i)
	{
        Pop[i].phen = intercept;

        for (int j = 0; j < n_alleles_g; ++j)
        {
            Pop[i].g[j] = intercept == 0 ? 0 : intercept/n_alleles_g;
        }

        for (int j = 0; j < n_alleles_b; ++j)
        {
            Pop[i].b[j] = 0;
        }

        for (int j = 0; j < n_alleles_m; ++j)
        {
            Pop[i].m_m[j] = 0;
            Pop[i].m_e[j] = 0;
            Pop[i].m_g[j] = 0;
        }
	}
}

// create an offspring
void Create_Kid(int mother, int father, Individual &kid)
{
    double sum_g = 0; // sum over all the breeding values of the offspring coding for the actual phenotype
    double sum_b = 0; // sum over all the breeding values of the offspring coding for the norm of reaction
    double sum_m_m = 0; // sum over all the breeding values of the offspring coding for the maternal effect
    double sum_m_g = 0; // sum over all the breeding values of the offspring coding for the maternal effect
    double sum_m_e = 0; // sum over all the breeding values of the offspring coding for the maternal effect

    // we assume all loci are unlinked 
    for (int i = 0; i < n_alleles_g;++i)
    {
        kid.g[i] = i % 2 == 0 ? Survivors[mother].g[i + gsl_rng_uniform_int(r, 2)] : Survivors[father].g[i - 1 + gsl_rng_uniform_int(r, 2)];
        MutateG(kid.g[i]);
        sum_g += kid.g[i];
    }

    for (int i = 0; i < n_alleles_b; ++i)
    {
        kid.b[i] = i % 2 == 0 ? Survivors[mother].b[i + gsl_rng_uniform_int(r, 2)] : Survivors[father].b[i - 1 + gsl_rng_uniform_int(r, 2)];
        MutateB(kid.b[i]);
        sum_b += kid.b[i];
    }

    for (int i = 0; i < n_alleles_m; ++i)
    {
        kid.m_m[i] = i % 2 == 0 ? Survivors[mother].m_m[i + gsl_rng_uniform_int(r, 2)] : Survivors[father].m_m[i - 1 + gsl_rng_uniform_int(r, 2)];
        MutateM_M(kid.m_m[i]);
        sum_m_m += kid.m_m[i];
        
        kid.m_e[i] = i % 2 == 0 ? Survivors[mother].m_e[i + gsl_rng_uniform_int(r, 2)] : Survivors[father].m_e[i - 1 + gsl_rng_uniform_int(r, 2)];
        MutateM_E(kid.m_e[i]);
        sum_m_e += kid.m_e[i];
        
        kid.m_g[i] = i % 2 == 0 ? Survivors[mother].m_g[i + gsl_rng_uniform_int(r, 2)] : Survivors[father].m_g[i - 1 + gsl_rng_uniform_int(r, 2)];
        MutateM_G(kid.m_g[i]);
        sum_m_g += kid.m_g[i];
    }

    kid.phen_m_m = sum_m_m;
    kid.phen_m_g = sum_m_g;
    kid.phen_m_e = sum_m_e;
    kid.phen_g = sum_g;
    kid.phen_b = sum_b;

    // complete maternal control
    kid.phen = kid.phen_g // elevation
        + gsl_ran_gaussian(r,sigma_e)  // developmental noise
        + kid.phen_b * epsilon_sens  // plasticity
        + kid.phen_m_m * Survivors[mother].phen // maternal cascading effect
        + kid.phen_m_e * Survivors[mother].envt // maternal environmental effect
        + Survivors[mother].phen_m_g; // maternal genetic effect

    kid.envt = epsilon_sens;

    assert(isnan(kid.phen) == 0);
}


// Survival of juveniles to reproductive adults
void Survive()
{
    double W;

    double theta = epsilon;

    if (!envt_is_changed)
    {
        if (generation > t_change && fabs(epsilon - 0) < 0.01)
        {
            rate = rateptb;
            intercept = intptb;
            ampl = amplptb;

            envt_is_changed = true;
        }
    }


    NSurv = 0;

    for (int i = 0; i < Npop; ++i)
    {
        W = wmin + (1.0 -  wmin) * exp(-.5 * (
                                            pow((Pop[i].phen - theta),2.0)/omega2 
                                            + pow(Pop[i].phen_b,2.0)/omega_b_2
                                            + pow(Pop[i].phen_m_m,2.0)/omega_m_m_2
                                            + pow(Pop[i].phen_m_g,2.0)/omega_m_g_2
                                            + pow(Pop[i].phen_m_e,2.0)/omega_m_e_2
                                        )
                                    );


        assert(isnan(W) == 0);

        // let individual survive or not
        if (gsl_rng_uniform(r) < W)
        {
            Survivors[NSurv++] = Pop[i];
        }
    }

    if (NSurv == 0)
    {
        WriteParameters();
        exit(1);
    }

    // update the environment for the next generation
    // as an autocorrelated gaussian random variable
    ksi = rho_t*ksi + gsl_ran_gaussian(r, sqrt(1.0-rho_t*rho_t)*sigma_ksi);

    epsilon = intercept + ampl * sin(rate * (generation+1)) + ksi;

    // in the likely case there is a developmental timelag, tau,
    // update the environment for a number of 'sub' timesteps
    // to achieve a 'sensed' (rather than real) value of epsilon
    if (tau > 0)
    {
        int timesteps = rint(1.0 / tau);

        for (int time_i = 0; time_i < timesteps; ++time_i)
        {
            ksi = rho_t*ksi + gsl_ran_gaussian(r, sqrt(1.0-rho_t*rho_t)*sigma_ksi);
        }
    }

    // update the value of the sensed environment
    epsilon_sens = intercept + ampl * sin(rate * (generation-tau+1)) + ksi;

    // finally, let the survivors produce N offspring 
    for (int i = 0; i < Npop; ++i)
    {
        Individual Kid;

        Create_Kid(gsl_rng_uniform_int(r,NSurv), gsl_rng_uniform_int(r,NSurv), Kid);

        Pop[i] = Kid;
    }
}


// write down summary statistics
void WriteData()
{
    double meanphen = 0;
    double ssphen = 0;
    double meang = 0;
    double ssg = 0;
    double meanm_m = 0;
    double meanm_g = 0;
    double meanm_e = 0;
    double ssm_m = 0;
    double ssm_g = 0;
    double ssm_e = 0;
    double meanb = 0;
    double ssb = 0;

    // get stats from the population
    for (int i =  0; i < Npop; ++i)
    {
        // stats for m
        meang += Pop[i].phen_g;
        ssg += Pop[i].phen_g * Pop[i].phen_g;

        // stats for m
        meanm_m += Pop[i].phen_m_m;
        ssm_m += Pop[i].phen_m_m * Pop[i].phen_m_m;
        
        meanm_g += Pop[i].phen_m_g;
        ssm_g += Pop[i].phen_m_g * Pop[i].phen_m_g;
        
        meanm_e += Pop[i].phen_m_e;
        ssm_e += Pop[i].phen_m_e * Pop[i].phen_m_e;

        meanb += Pop[i].phen_b;
        ssb += Pop[i].phen_b * Pop[i].phen_b;

        meanphen += Pop[i].phen;

        ssphen += Pop[i].phen * Pop[i].phen;
    }

    DataFile << generation << ";" << epsilon << ";" << NSurv << ";" << ksi << ";";

    DataFile 
            << (meanphen/Npop) << ";"
            << ((ssphen/Npop) - pow(meanphen/Npop,2.0)) << ";"
            << (meang/Npop) << ";"
            << (ssg/Npop - pow(meang/Npop,2.0)) << ";"
            << (meanm_m/Npop) << ";"
            << (ssm_m/Npop - pow(meanm_m/Npop,2.0)) << ";" 
            << (meanm_g/Npop) << ";"
            << (ssm_g/Npop - pow(meanm_g/Npop,2.0)) << ";" 
            << (meanm_e/Npop) << ";"
            << (ssm_e/Npop - pow(meanm_e/Npop,2.0)) << ";" 
            << (meanb/Npop) << ";"
            << (ssb/Npop - pow(meanb/Npop,2.0)) << ";"  << endl;
}

// write the headers of a datafile
void WriteDataHeaders()
{
    DataFile << "generation;epsilon;nsurv;ksi;meanz;varz;meang;varg;meanm_m;varm_m;meanm_g;varm_g;meanm_e;varm_e;meanb;varb;" << endl;
}


// the guts of the code
int main(int argc, char ** argv)
{
	initArguments(argc, argv);
	WriteDataHeaders();
	Init();

	for (generation = 0; generation <= NumGen; ++generation)
	{
        do_stats = generation % skip == 0;

		Survive();

        if (do_stats)
		{
			WriteData();
		}
	}

	WriteParameters();
}
