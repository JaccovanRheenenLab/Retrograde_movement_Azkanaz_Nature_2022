const int N1=20; // number of cells in the crypt-villus axis (vertical)
const int N2=5;  // number of cells in the horizontal axis
const int N=100; // total number of cells (N1*N2)
double clones_num=200.; //number of clones being simulated per repetition (should be matched to actual number being considered in experiments)
const int rep=1000;    // number of repetitions
int T=400;      // total timesteps of simulations


// key parameters of the simulations
double divrate=0.1; // rate of cell division
double inter=0.025; // rate of random cell intercalation. Taken for values of cecum

// indices used during the simulatoins
int counts=0;
int counts2=0;
int counts3=0;
int ch1=0;
int ch2=0;
int ch1_new=0;
int ch2_new=0;
int t;      // time variable

// arrys
double surv[4][400];   // array to store the survival probability as a function of starting position (row number from 0 to 3) and time (0 to 400)
double monoclonal[400];     // array to store the fraction of monoclonal crypts
double av_size[400][2];     // array to store the average size of clones in the crypts as a function of time
double count_pers[4];       // array to store the number of time a clone is induced from position/row i in a crypt (from 0 to 3)
int type[N1][N2];   // array representing the crypt. 0 indicates unlabelled cell, 1 indicates labelled cells
int type_av[N1][N2];    // same array serving as a memory between division events

