/*
   A simulator with various combinations of extensions for the model of Evolving Open Systems (for the detail of the model, see: http://www.nature.com/articles/srep04082)

USAGE: ./a.out model_TYPE, #_of_interactions_per_species, SimulationTime, iseed

< model_TYPE >
model_TYPE % 2 -> 0: fixed M, 1: flat degree distribution in (1,M)
(model_TYPE/2) % 2 -> 0: Gaussian link weight, 1: flat weight distribution
Therefore, 0: standard model, 1: flat degree, 2: flat link weight, 3: flat deg. & weight
*/
#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <random>

// CONSTANT PARAMETERS
const int Maxtime = 100000000;
const int maxN = 10000;// It is set to be small, for the casual test on OACIS. One need larger maxN for the determination of transition point.
const int maxNINT = maxN*100;
const int InitialN = 100;

// Parameters for Observation
const int interval = 1000;// the interval for observation

// File Handles
char filename[99];
FILE *OACISfile;
FILE *skimfile;

std::mt19937_64 * pRnd;

// Functions
inline int min(int a, int b){return b < a ? b : a;}
double Gaussian(void);
int FindExt(int nspecies, int n_incubate, double *fitness);

// MAIN
int main(int argc, char *argv[]){
  if( argc != 5 ) {
    std::cerr << "Error! wrong number of arguments." << std::endl
      << "USAGE: ./a.out model_TYPE, #_of_interactions_per_species, SimulationTime, iseed" << std::endl
      << "  < model_TYPE >" << std::endl
      << "    model_TYPE % 2 -> 0: fixed M, 1: flat degree distribution in (1,M)" << std::endl
      << "    (model_TYPE/2) % 2 -> 0: Gaussian link weight, 1: flat weight distribution" << std::endl
      << "    Therefore, 0: standard model, 1: flat degree, 2: flat link weight, 3: flat deg. & weight" << std::endl;
    exit(1);
  }
  // Get the Number of Interactions
  int model_TYPE = atoi(argv[1]);
  int FlatDegree = model_TYPE % 2;
  int FlatWeight = (model_TYPE/2) % 2;
  int M = atoi(argv[2]);
  int SimTime = atoi(argv[3]);
  int iseed = atoi(argv[4]);
  iseed %= 1000;

  // File handles
  //// for OACIS
  if((OACISfile=fopen("_output.json","w"))==NULL){
    puts("HIST_F:Save to file failed");
    return 0;
  }
  //// main output files
  int iname = iseed;
  sprintf(filename, "eos_type%01d_m%02di%03d.dat", model_TYPE, M, iname);
  if((skimfile=fopen(filename,"w"))==NULL){
    exit(0);
  }
  fprintf(skimfile, "# Time, N, <f>, <m>, AC, ACfiti, ACfito, ACfoti, ACfoto, CC, CC/CCrandom, CCrandom\n");

  //Interaction Arrays
  int* from = new int[maxNINT];	
  int* to = new int[maxNINT];
  double* amp = new double[maxNINT];
  //Profile Arrays
  int* generation = new int[maxN];
  int* debut = new int[maxN];
  int* income = new int[maxN];
  int* outgo = new int[maxN];
  double* fitness = new double[maxN];
  //Statistics Arrays
  int* deg = new int[maxN];

  // Initial Settings
  //Giving a random number seed
  pRnd = new std::mt19937_64(iseed);
  int n = InitialN;
  int nlast = n;
  int nint;
  double sumf;

  for(int i=0; i< maxN; i++){
    debut[i] = -1;
    generation[i] = -1;		
    income[i] = -1;
    outgo[i] = -1;
    fitness[i] = 0.0;
  }
  for(int i=0; i< maxNINT; i++){
    from[i] = -1;
    to[i] = -1;
    amp[i] = 0.0;
  }

  // Give Initial Interactions
  nint = 0;
  for(int i=0; i < n; i++){
    debut[i] = 0;
    generation[i] = 1; // Generation starts from 1 (Note that the notation in the first paper was 0)
    int nbond = M;
    if(n < M){nbond = n;}
    for(int j = 0; j < nbond; j++){
      int tmpwith = i;
      while(tmpwith == i){// Get an off-diagonal connection
        std::uniform_int_distribution<int> dist(0,n-1);
        tmpwith = dist( *pRnd );
      }
      std::uniform_int_distribution<int> dist01(0,1);
      if( dist01(*pRnd) ){
        from[nint] = i;
        to[nint] = tmpwith;
      }
      else{
        from[nint] = tmpwith;
        to[nint] = i;				
      }
      amp[nint] = Gaussian();
      nint++;
    }
  }

  //// OACIS observers
  double oacis_nsample = 0;
  double oacis_sumn = 0.0;
  double oacis_sumt = 0.0;
  double oacis_sumtt = 0.0;
  double oacis_sumnt = 0.0;
  int oacis_maxn = 0;
  double oacis_avef = 0.0;
  double oacis_avem = 0.0;
  double oacis_asso = 0.0;
  double oacis_fiti = 0.0;
  double oacis_fito = 0.0;
  double oacis_foti = 0.0;
  double oacis_foto = 0.0;
  double oacis_cc = 0.0;
  double oacis_ccratio = 0.0;

  // BODY
  if(SimTime > Maxtime) SimTime = Maxtime;
  for(int t = 1; t < SimTime; t++){
    // Check the Limit
    if(n > maxN-1){
      break;
    }
    // Add a New Species with linear preference to the degree (if the system is sufficiently large)
    debut[n] = t;
    generation[n] = 1; // Generation starts from 1 (Note that the notation in the first paper was 0)
    int nbond = M;
    if(FlatDegree){
      std::uniform_int_distribution<int> dist(1, 2*M);
      nbond = dist( *pRnd );
    }
    if(n < M){nbond = n;}

    for(int j = 0; j < nbond; j++){
      std::uniform_int_distribution<int> dist(0,n-1);
      int resident = dist( *pRnd );
      std::uniform_int_distribution<int> dist01(0,1);
      if( dist01(*pRnd) ){
        from[nint+j] = n;
        to[nint+j] = resident;
        generation[resident]++;
      }
      else{
        from[nint+j] = resident;
        to[nint+j] = n;
      }
      if(FlatWeight){
        std::uniform_real_distribution<double> uni(-1.0,1.0);
        amp[nint+j] = uni(*pRnd);
      }else{
        amp[nint+j] = Gaussian();
      }
    }
    n++;
    nint += nbond;

    // Calculate the Fitness
    sumf = 0.0;
    for(int i = 0; i < n; i++){
      income[i] = 0;
      outgo[i] = 0;
      fitness[i] = 0.0;
    }
    for(int i=0; i < nint; i++){
      outgo[from[i]]++;
      income[to[i]]++;
      fitness[to[i]] += amp[i];
      sumf += amp[i];
    }
    // Find the Least Fit Species
    int here = FindExt(n, M, fitness);

    // (Succesive) Extinctions
    while(here > -1){
      int incdec = -1;
      if(here < n-1){incdec = 1;}
      for(int i = 0; i < nint; i++){if(from[i] == here){generation[to[i]] += incdec;}}
      // Pruning the ID Records
      for(int i = here; i < n; i++){
        debut[i] = debut[i+1];
        generation[i] = generation[i+1];
      }
      // Pruning the Interaction
      int cutcounter = 0;
      for(int i = 0; i < nint; i++){
        while((from[i] == here) || (to[i] == here)){
          // main body
          cutcounter++;
          for(int j = i; j < nint; j++){
            from[j] = from[j+1];
            to[j] = to[j+1];
            amp[j] = amp[j+1];
          }
        }
      }
      // Update the Species Index
      for(int i = 0; i < nint; i++){
        if(from[i] > here){from[i]--;}
        if(to[i] > here){to[i]--;}				
      }
      n--;
      nint -= cutcounter;
      // ReCalculate the Fitness
      sumf = 0.0;
      for(int i = 0; i < n; i++){
        income[i] = 0;
        outgo[i] = 0;
        fitness[i] = 0.0;
      }
      for(int i=0; i < nint; i++){
        outgo[from[i]]++;
        income[to[i]]++;
        fitness[to[i]] += amp[i];
        sumf += amp[i];
      }			
      // Find the Least Fit Species
      here = FindExt(n, M, fitness);
    }
    // Observation Part
    if(t % interval == 0){
      if(n < M+1){
        // Outout the Time Series WITHOUT NETWORK INFORMATION, because the system is under the Incubation Rule
        fprintf(skimfile, "%d %d %f %f ", t, n, sumf/n, ((double) nint)/n);
        fprintf(skimfile, "%f %f %f %f %f ", 0.0, 0.0, 0.0, 0.0, 0.0);
        fprintf(skimfile, "%f %f\n", 0.0, 0.0);
      }
      else{
        // Output
        ///// calculating assortativity and clustering coefficients
        //// Prepare the degree & connectivity array
        //!!!// CREATION OF CONNECTION ARRAY: must be deleted within this parse //!!!//
        for(int i=0; i<n; i++){deg[i] = 0;}
        for(int i=0; i<nint; i++){
          deg[from[i]]++;
          deg[to[i]]++;
        }
        int** connected = new int*[n];
        for(int i=0; i<n; i++){
          connected[i] = new int[deg[i]];
        }
        // must refresh deg[i]
        for(int i=0; i<n; i++){deg[i] = 0;}
        for(int i=0; i<nint; i++){
          int me = from[i];
          int you = to[i];
          connected[me][deg[me]] = you;
          connected[you][deg[you]] = me;
          deg[me]++;
          deg[you]++;
        }
        // calculate assortativities
        double fk = 0.0, fkk = 0.0;
        double fi = 0.0, fisq = 0.0;
        double fo = 0.0, fosq = 0.0;
        double tk = 0.0, tkk = 0.0;
        double ti = 0.0, tisq = 0.0;
        double tout = 0.0, tosq = 0.0; // namig for avoiding "to"
        double ft = 0.0;
        double fiti = 0.0, fito = 0.0, foti = 0.0, foto = 0.0;

        for(int i=0; i<nint; i++){
          int kfrom = deg[from[i]];
          int kfromin = income[from[i]];
          int kfromout = outgo[from[i]];
          int kto = deg[to[i]];
          int ktoin = income[to[i]];
          int ktoout = outgo[to[i]];

          fk += kfrom;
          fkk += kfrom*kfrom;
          fi += kfromin;
          fisq += kfromin*kfromin;
          fo += kfromout;
          fosq += kfromout*kfromout;
          tk += kto;
          tkk += kto*kto;
          ti += ktoin;
          tisq += ktoin*ktoin;
          tout += ktoout;
          tosq += ktoout*ktoout;
          ft += kfrom*kto;
          fiti += kfromin*ktoin;
          fito += kfromin*ktoout;
          foti += kfromout*ktoin;
          foto += kfromout*ktoout;
        }
        double aveft = ft/nint;
        double avedeg = 0.5*(fk + tk)/nint;
        double avedegsq = 0.5*(fkk + tkk)/nint;
        double assortativity = 0.0;
        if((avedegsq - avedeg*avedeg) > 0.0){assortativity = (aveft - avedeg*avedeg)/(avedegsq - avedeg*avedeg);}

        fi /= nint;
        fo /= nint;
        ti /= nint;
        tout /= nint;
        double devfi = sqrt(fisq/nint - fi*fi);
        double devfo = sqrt(fosq/nint - fo*fo);
        double devti = sqrt(tisq/nint - ti*ti);
        double devto = sqrt(tosq/nint - tout*tout);
        double assfiti = 0.0, assfito = 0.0, assfoti = 0.0, assfoto = 0.0;
        if(devfi*devfo*devti*devto > 0.0){
          assfiti = (fiti/nint - fi*ti)/(devfi*devti);
          assfito = (fito/nint - fi*tout)/(devfi*devto);
          assfoti = (foti/nint - fo*ti)/(devfo*devti);
          assfoto = (foto/nint - fo*tout)/(devfo*devto);
        }
        // calculate clustering
        int ntriple = 0;
        int ntriangle = 0;
        for(int me=0; me<n; me++){
          for(int i=0; i<deg[me]; i++){
            int you = connected[me][i];
            for(int j=0; j<deg[you]; j++){
              int him = connected[you][j];
              for(int k=0; k<deg[him]; k++){
                if(connected[him][k] == me){ntriangle++;}
                ntriple++;
              }
            }
          }
        }
        double cc = 0.0;
        if(ntriple > 0){cc = 3.0*ntriangle/ntriple;}
        double cctorandom = cc;
        if(nint > 1){cctorandom *= n*(n-1)/((double) nint);}

        // DELETION of CONNECTIVITY ARRAY
        for(int i=0; i<n; i++){delete[] connected[i];}
        delete[] connected;

        // Outout the Time Series
        fprintf(skimfile, "%d %d %f %f ", t, n, sumf/n, ((double) 2*nint)/n);
        fprintf(skimfile, "%f %f %f %f %f ", assortativity, assfiti, assfito, assfoti, assfoto);
        fprintf(skimfile, "%f %f\n", cc, cctorandom);

        oacis_nsample++;
        double sct = ((double) t)/((double) SimTime);
        oacis_sumn += n;
        oacis_sumt += sct;
        oacis_sumtt += sct*sct;
        oacis_sumnt += n*sct;
        if(oacis_maxn < n) oacis_maxn = n;
        oacis_avef += sumf/n;
        oacis_avem += ((double) 2*nint)/n;
        oacis_asso += assortativity;
        oacis_fiti += assfiti;
        oacis_fito += assfito;
        oacis_foti += assfoti;
        oacis_foto += assfoto;
        oacis_cc += cc;
        oacis_ccratio += cctorandom;
      }// normal/incubate situations
    }// observe
  }// time

  // The End of the Simulation
  if(oacis_nsample == 0){oacis_nsample = 1;}
  double norm = 1.0/((double) oacis_nsample);
  double aven = oacis_sumn*norm;
  double avet = oacis_sumt*norm;
  double avett = oacis_sumtt*norm;
  double avent = oacis_sumnt*norm;
  oacis_avef *= norm;
  oacis_avem *= norm;
  fprintf(OACISfile, "{\n");
  fprintf(OACISfile, "\"<N>\": %f,\n", aven);
  fprintf(OACISfile, "\"max N\": %d,\n", oacis_maxn);
  // least-square fit to the linear function
  fprintf(OACISfile, "\"Divergence Speed\": %f,\n", (avent - aven*avet)/(avett - avet*avet)/SimTime);
  fprintf(OACISfile, "\"Intercept\": %f,\n", (avett*aven - avent*avet)/(avett - avet*avet));
  fprintf(OACISfile, "\"Average Fitness\": %f,\n", oacis_avef);
  fprintf(OACISfile, "\"Average Degree\": %f,\n", oacis_avem);
  fprintf(OACISfile, "\"Assortativity\": %f,\n", oacis_asso*norm);
  fprintf(OACISfile, "\"Clustering Coefficient\": %f,\n", oacis_cc*norm);
  fprintf(OACISfile, "\"C.C./random C.C.\": %f\n", oacis_ccratio*norm);
  fprintf(OACISfile, "}\n");
  fclose(OACISfile);

  // FINALIZATION
  fclose(skimfile);
  delete[] from;
  delete[] to;
  delete[] amp;
  delete[] debut;
  delete[] generation;
  delete[] income;
  delete[] outgo;
  delete[] deg;
  delete[] fitness;
}

///// Return a Random Number from the Gaussian Distribution with <r> = 0.0, <r^2> = 1.0
double Gaussian(void) {
  std::normal_distribution<double> nd;
  return nd(*pRnd);
}

///// Return the Least Fit Species
int FindExt(int nspecies, int n_incubate, double *fitness){
  int here = -1;
  double min = 0.0;
  if(nspecies < n_incubate){min=-0.0000000001;}
  for(int i = 0; i < nspecies; i++){
    if(min >= fitness[i]){
      min = fitness[i];
      here = i;
    }
  }
  return here;
}
