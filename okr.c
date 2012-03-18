/*
  Title: okr.c: The program used in the following article:
                Yamazaki T, Nagao S (2012)
                A computational mechanism for unified gain and timing control
                in the cerebellum
                PLoS ONE 7(3): e33319. doi:10.1371/ journal.pone.0033319
                http://dx.plos.org/10.1371/journal.pone.0033319
  Author: YAMAZAKI, Tadashi (Neuralgorithm)
  Date: March 14, 2012 
  License: GPLv2
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_errno.h>
#include<gsl/gsl_odeiv.h>
#include<time.h>
#include<sys/times.h>
#include<unistd.h>

#define TH_GR		-35.0
#define C_GR		3.1
#define R_GR		2.3
#define E_leak_GR	-58.0
#define g_leak_GR	0.43
#define E_ex_GR		0
#define g_ex_GR		0.18
#define E_inh_GR	-82.0
#define g_inh_GR	0.028
#define E_ahp_GR	-82.0
#define g_ahp_GR	1.0
#define TH_GO		-52.0
#define C_GO		28.0
#define R_GO		0.428
#define E_leak_GO	-55.0
#define g_leak_GO	(1.0/(R_GO))
#define E_ex_GO		0.0
#define g_ex_GO		45.5
#define E_ahp_GO	-72.7
#define g_ahp_GO	20.0
#define I		1.0

#define tau_ampa_gogr	1.5
#define tau1_nmda_gogr	31.0
#define tau2_nmda_gogr	170.0
#define r_ampa_gogr	0.8
#define r1_nmda_gogr	(0.2*0.33)
#define r2_nmda_gogr	(0.2*0.67)
#define tau1_gaba_grgo	7.0
#define tau2_gaba_grgo	59.0
#define r1_gaba_grgo	0.43
#define r2_gaba_grgo	0.57

#define R_N2	100
#define X_GR	(10*(X_GO))
#define Y_GR    (10*(Y_GO))
#define N_GR	((X_GR)*(Y_GR))
#define X_GO	32
#define Y_GO	32
#define N_GO	((X_GO)*(Y_GO))
#define RX	((X_GR)/(X_GO))
#define RY	((Y_GR)/(Y_GO))
#define DT	1.0

#define r_ahp_gr	1.0
#define r_ahp_go	1.0
#define tau_ahp_gr	5.0
#define tau_ahp_go	5.0

#define X_PKJ	16
#define Y_PKJ	1
#define N_PKJ	((X_PKJ)*(Y_PKJ))

#define TH_PKJ		-55.0
#define C_PKJ		106.0
#define R_PKJ		0.431
#define g_leak_PKJ	(1.0/(R_PKJ))
#define E_leak_PKJ	-68.0
#define g_ex_PKJ	0.7
#define E_ex_PKJ	0
#define g_ahp_PKJ	100.0
#define E_ahp_PKJ	-70.0
#define g_inh_PKJ       1.0
#define E_inh_PKJ       -75.0

#define N_BS            N_PKJ
#define TH_BS		-55.0
#define C_BS		106.0
#define R_BS		0.431
#define g_leak_BS	(1.0/(R_PKJ))
#define E_leak_BS	-68.0
#define g_ex_BS	        0.7
#define E_ex_BS	        0
#define g_ahp_BS	100.0
#define E_ahp_BS	-70.0

#define N_VN	1
#define TH_VN		-38.8
#define C_VN		122.3
#define R_VN		0.61
#define g_leak_VN	(1.0/(R_VN))
#define E_leak_VN	-56.0
#define g_ex_VN		50.0
#define E_ex_VN		0
#define g_inh_VN	30.0
#define E_inh_VN	-88.0
#define g_ahp_VN	50.0
#define E_ahp_VN	-70.0

#define N_MF_PER_VN	100

#define N_IO		1
#define TH_IO		-50.0
#define C_IO		1.0
#define g_leak_IO	0.015
#define R_IO		(1.0/(g_leak_IO))
#define E_leak_IO	-60.0
#define g_ex_IO		0.1
#define E_ex_IO		0.0
#define g_inh_IO	0.018
#define E_inh_IO	-75.0
#define g_ahp_IO	1.0
#define E_ahp_IO	-70.0

#define tau_ampa_pkjpf	8.3
#define r_ampa_pkjpf	1.0
#define tau_ahp_pkj	2.5
#define r_ahp_pkj	1.0
#define tau_gaba_pkjbs  10.0
#define r_gaba_pkjbs    1.0

#define tau_ampa_bspf	8.3
#define r_ampa_bspf	1.0
#define tau_ahp_bs	2.5
#define r_ahp_bs	1.0

#define tau_gaba_vnpkj	42.3
#define r_gaba_vnpkj	1.0
#define tau_ampa_vnmf	9.9
#define r_ampa_vnmf	0.66
#define tau_nmda_vnmf	30.5
#define r_nmda_vnmf	(1.0-(r_ampa_vnmf))
#define tau_ahp_vn	5.0
#define r_ahp_vn	1.0

#define r_ampa_grmf	0.88
#define tau_ampa_grmf	1.2
#define r_nmda_grmf	0.12
#define tau_nmda_grmf	52.0
#define N_MF_PER_GR 4
#define lambda 4.0

#define tau_ampa_preio	10.0
#define tau_gaba_iovn	10.0
#define r_ampa_preio	1.0
#define r_gaba_iovn	1.0
#define kappa_preio	1.0
#define tau_ahp_io	5.0
#define r_ahp_io	1.0

#define tau_ampa_vnio	9.9
#define r_ampa_vnio	1.0
#define tau_nmda_vnio	30.6
#define r_nmda_vnio	(1.0-(r_ampa_vnio))

#define GR(i)	((i))
#define GO(i)	((N_GR)+(i))
#define PKJ(x)	((N_GR)+(N_GO)+(x))
#define VN(x)	((N_GR)+(N_GO)+(N_PKJ)+(x))
#define IO(x)	((N_GR)+(N_GO)+(N_PKJ)+(N_VN)+(x))
#define BS(x)	((N_GR)+(N_GO)+(N_PKJ)+(N_VN)+(N_IO)+(x))
#define N_ALL	((N_GR)+(N_GO)+(N_PKJ)+(N_VN)+(N_IO)+(N_BS))

#define N_MOL	((N_PKJ)+(N_VN)+(N_IO)+(N_BS))
#define MPKJ(x)	((x))
#define MVN(x)	((N_PKJ)+(x))
#define MIO(x)	((N_PKJ)+(N_VN)+(x))
#define MBS(x)	((N_PKJ)+(N_VN)+(N_IO)+(x))

#define conn_prob_gogr 0.05 // go <- gr
#define conn_prob_grgo 0.025 // gr <- go
#define conn_weight_grgo 10.0 // gr <- go
#define conn_weight_gogr (0.2/(49.0*(R_N2))) // go <- gr
#define conn_weight_vnpkj (0.12/(N_PKJ)) // vn <- pkj
#define conn_weight_pkjpf 0.003 // pkj <- gr
#define conn_weight_bspf 300.0 // bs <- pf
#define conn_weight_vnmf (0.2/(N_MF_PER_VN)) // vn <- mf
#define conn_weight_iovn 4.0 // io <- vn
#define conn_weight_vnio 1.0 // vn <- io

#define n_trials 300

// stimulus
#define oscillation_period 2000 // ms
#define MF_maxfreq 0.03 // spikes/ms == 30 spikes/s
#define CF_maxfreq 0.003 // spikes/ms == 3 spikes/s

// synaptic plasticity
#define c_ltd 0.995
#define c_ltp 0.0001
#define weight_file_prefix "w"

double **psp_ampa_grmf, **psp_nmda_grmf, *g_grmf;

time_t time0, time1;
double ahp_go[N_GO], *ahp_gr; //, adp_go[N_GO];
double *g_grgo, g_gogr[N_GO];

int **list_gogr, **list_grgo, nlist_gogr[N_GO], *nlist_grgo;

int *grs_for_ltd, **grs_in_window;
const int ltd_window_size = 50;

double **w_pkjpf, **w_bspf;

double ex[N_MOL], inh[N_MOL], ahp[N_MOL];

// 疑似乱数
extern void init_genrand(unsigned long);
extern double genrand_real2(void);

void init_weight_pkjpf(void)
{
  int i;
  w_pkjpf = (double **)malloc(N_PKJ*sizeof(double *));

  for(i = 0; i < N_PKJ; i++){
    w_pkjpf[i] = (double *)malloc(N_GR*sizeof(double));
  }
}
void free_weight_pkjpf(void)
{
  int i;

  for(i = 0; i < N_PKJ; i++){
    free(w_pkjpf[i]);
  }
  free(w_pkjpf);
}
void init_weight_bspf(void)
{
  int i;
  w_bspf = (double **)malloc(N_BS*sizeof(double *));

  for(i = 0; i < N_BS; i++){
    w_bspf[i] = (double *)malloc(N_GR*sizeof(double));
  }
}
void free_weight_bspf(void)
{
  int i;

  for(i = 0; i < N_BS; i++){
    free(w_bspf[i]);
  }
  free(w_bspf);
}
int read_weight_pkjpf(const int nt)
{
  FILE *file;
  int i, j;
  char buf[1024];

  sprintf(buf, "%s.%d", weight_file_prefix, nt);
  file = fopen(buf, "r");
  if (!file){
    fprintf(stderr, "cannot open %s\n", buf);
    return 1;
  }
  for(i = 0; i < N_PKJ; i++){
    for(j = 0; j < N_GR; j++){
      fgets(buf, 1024, file);
      w_pkjpf[i][j] = atof(buf);
    }
  }
  fclose(file);
  return 0;
}
int read_weight_bspf(const int nt)
{
  FILE *file;
  int i, j;
  char buf[1024];

  sprintf(buf, "%s.%d", weight_file_prefix, nt);
  file = fopen(buf, "r");
  if (!file){
    fprintf(stderr, "cannot open %s\n", buf);
    return 1;
  }
  for(i = 0; i < N_BS; i++){
    for(j = 0; j < N_GR; j++){
      fgets(buf, 1024, file);
      w_bspf[i][j] = atof(buf);
    }
  }
  fclose(file);
  return 0;
}
void write_weight_pkjpf(const int nt)
{
  FILE *file;
  int i, j;
  char buf[1024];

  sprintf(buf, "%s.%d", weight_file_prefix, nt+1);
  file = fopen(buf, "w");
  if (!file){
    fprintf(stderr, "cannot open %s\n", buf);
    return;
  }
  for(i = 0; i < N_PKJ; i++){
    for(j = 0; j < N_GR; j++){
      fprintf(file, "%f\n", w_pkjpf[i][j]);
    }
  }
  fclose(file);
}
int GR_read(void)
{
  int i, t;

  grs_for_ltd = (int *)malloc(1000*10240*sizeof(int));
  grs_in_window = (int **)malloc(ltd_window_size*sizeof(int *));
  
  for(i = 0; i < 1000*10240; i++){
    grs_for_ltd[i] = -1;
  }
  for(t = 0; t < ltd_window_size; t++){
    grs_in_window[t] = (int *)malloc(10240*sizeof(int));
  }
  return 0;
}
void GR_free(void)
{
  int t;

  free(grs_for_ltd);

  for(t = 0; t < ltd_window_size; t++){
    free(grs_in_window[t]);
  }
  free(grs_in_window);

  //free(ngrspikes);
}
void LTD(void)
{
  int i, j;

  for(j = 0; grs_for_ltd[j] != -1; j++){
    for(i = 0; i < X_PKJ; i++){
      if (w_pkjpf[i][grs_for_ltd[j]] > 0){
	w_pkjpf[i][grs_for_ltd[j]] *= c_ltd;
      }
    }
  }
}
void init_wlist(void)
{
  int i;
  list_gogr = (int **)malloc(N_GO*sizeof(int *));
  for(i = 0; i < N_GO; i++){
    list_gogr[i] = (int *)malloc(N_GR*sizeof(int));
    nlist_gogr[i] = 0;
  }
  list_grgo = (int **)malloc(N_GR*sizeof(int *));
  for(i = 0; i < N_GR; i++){
    list_grgo[i] = (int *)malloc(N_GO*sizeof(int));
    nlist_grgo[i] = 0;
  }
}
void free_wlist(void)
{
  int i;
  for(i = 0; i < N_GO; i++){
    free(list_gogr[i]);
  }
  free(list_gogr);
  for(i = 0; i < N_GR; i++){
    free(list_grgo[i]);
  }
  free(list_grgo);
}
void init_w_gogr(void)
{
  int gox, goy;
  int goax, goay, godx, gody;
  int gon, grn, goan;
  int i;

  for(gox = 0; gox < X_GO; gox++){
    for(goy = 0; goy < Y_GO; goy++){
      gon = goy + Y_GO*gox;
      for(godx = -3; godx <= 3; godx++){
	goax = gox + godx;
	if (goax >= X_GO){ goax -= X_GO; }
	if (goax < 0){ goax += X_GO; }
	for(gody = -3; gody <= 3; gody++){
	  goay = goy + gody;
	  if (goay >= Y_GO){ goay -= Y_GO; }
	  if (goay < 0){ goay += Y_GO; }
	  goan = goay + Y_GO*goax;
	  if (genrand_real2() < conn_prob_gogr){
	    for(i = 0; i < R_N2; i++){
	      grn = i+R_N2*goan;
	      list_gogr[gon][nlist_gogr[gon]] = grn;
	      nlist_gogr[gon]++;
	    }
	  }
	}
      }
    }
  }
}

void init_w_grgo(void)
{
  int gox, goy;
  int goax, goay, godx, gody;
  int gon, grn, goan;
  int i;

  int **arr, naxons[N_GO];

  arr = (int **)malloc(sizeof(int *)*N_GO);
  for(i = 0; i < N_GO; i++){
    arr[i] = (int *)malloc(sizeof(int)*N_GO);
  }
  for(i = 0; i < N_GO; i++){
    naxons[i] = 0;
  }
  for(gox = 0; gox < X_GO; gox++){
    for(goy = 0; goy < Y_GO; goy++){
      gon = goy + Y_GO*gox;
      for(godx = -4; godx <= 4; godx++){
	goax = gox + godx;
	if (goax >= X_GO){ goax -= X_GO; }
	if (goax < 0){ goax += X_GO; }
	for(gody = -4; gody <= 4; gody++){
	  goay = goy + gody;
	  if (goay >= Y_GO){ goay -= Y_GO; }
	  if (goay < 0){ goay += Y_GO; }
	  goan = goay + Y_GO*goax;
	  if (genrand_real2() < conn_prob_grgo){
	    arr[goan][naxons[goan]] = gon;
	    naxons[goan]++;
	  }
	}
      }
    }
  }
  for(gox = 0; gox < X_GO; gox++){
    for(goy = 0; goy < Y_GO; goy++){
      gon = goy + Y_GO*gox;
      grn = gon;
      for(godx = 0; godx <= 1; godx++){
	goax = gox + godx;
	if (goax >= X_GO){ goax -= X_GO; }
	for(gody = 0; gody <= 1; gody++){
	  goay = goy + gody;
	  if (goay >= Y_GO){ goay -= Y_GO; }
	  goan = goay + Y_GO*goax;
	  for(i = 0; i < naxons[goan]; i++){
	    list_grgo[grn][nlist_grgo[grn]] = arr[goan][i];
	    nlist_grgo[grn]++;
	  }
	}
      }
    }
  }
  for(i = 0; i < N_GO; i++){
    free(arr[i]);
  }
  free(arr);
}
void initialize(unsigned long seed)
{
  init_genrand(seed);
  init_wlist();
  init_w_grgo();
  init_w_gogr();
  return;
}
int func(double t, const double u[], double du[], void *params)
{
  int i;

  for(i = 0; i < N_GR; i++){
    du[GR(i)] = (1.0/C_GR)*(-g_leak_GR*(u[GR(i)]-E_leak_GR)
			    -g_inh_GR*g_grgo[GR(i)]*(u[GR(i)]-E_inh_GR)
			    -g_ahp_GR*ahp_gr[GR(i)]*(u[GR(i)]-E_ahp_GR)
			    -g_ex_GR*g_grmf[GR(i)]*(u[GR(i)]-E_ex_GR));
  }
  for(i = 0; i < N_GO; i++){
    du[GO(i)] = (1.0/C_GO)*(-g_leak_GO*(u[GO(i)]-E_leak_GO)
			    -g_ex_GO*g_gogr[i]*(u[GO(i)]-E_ex_GO)
			    -g_ahp_GO*ahp_go[i]*(u[GO(i)]-E_ahp_GO));;
  }
  for(i = 0; i < N_PKJ; i++){
    du[PKJ(i)] = (1.0/C_PKJ)*(-g_leak_PKJ*(u[PKJ(i)] - E_leak_PKJ)
			      -g_ex_PKJ*ex[MPKJ(i)]*(u[PKJ(i)] - E_ex_PKJ)
			      -g_inh_PKJ*inh[MPKJ(i)]*(u[PKJ(i)] - E_inh_PKJ)
			      -g_ahp_PKJ*r_ahp_pkj*
			      ahp[MPKJ(i)]*(u[PKJ(i)] - E_ahp_PKJ)
			      +250.0
			      );
  }
  for(i = 0; i < N_VN; i++){
    du[VN(i)] = (1.0/C_VN)*(-g_leak_VN*(u[VN(i)] - E_leak_VN)
			    -g_ex_VN*ex[MVN(i)]*(u[VN(i)] - E_ex_VN)
			    -g_inh_VN*inh[MVN(i)]*(u[VN(i)] - E_inh_VN)
			    -g_ahp_VN*r_ahp_vn*
			    ahp[MVN(i)]*(u[VN(i)] - E_ahp_VN)
			    +700.0
			    );
  }
  for(i = 0; i < N_IO; i++){
    du[IO(i)] = (1.0/C_IO)*(-g_leak_IO*(u[IO(i)] - E_leak_IO)
			    -g_ex_IO*ex[MIO(i)]*(u[IO(i)] - E_ex_IO)
			    -g_inh_IO*inh[MIO(i)]*(u[IO(i)] - E_inh_IO)
			    -g_ahp_IO*r_ahp_io*
			    ahp[MIO(i)]*(u[IO(i)] - E_ahp_IO)
			    );
  }
  for(i = 0; i < N_BS; i++){
    du[BS(i)] = (1.0/C_BS)*(-g_leak_BS*(u[BS(i)] - E_leak_BS)
			    -g_ex_BS*ex[MBS(i)]*(u[BS(i)] - E_ex_BS)
			    -g_ahp_BS*r_ahp_bs*
			    ahp[MBS(i)]*(u[BS(i)] - E_ahp_BS)
			    );
  }

  return GSL_SUCCESS;
}
int main(int argc, char *argv[])
{
  double t, r;
  int i, j, k;
  //FILE *fgr_log, *fgo_log, *fpkjvnio_log;
  FILE *fgr_spk, *fgo_spk, *fgr_spk_p, *fpkjvnio_spk;
  int nspikegr, nspikego;
  char filename[1024];
  double mfseed;

  double *u;
  double *dudt_in, *dudt_out, *u_err;

  double psp_go_gaba1[N_GO], psp_go_gaba2[N_GO], psp_go_gaba3[N_GO];
  double *psp_gr_ampa, *psp_gr_nmda1, *psp_gr_nmda2;
  double psp_go[N_GO], *psp_gr;
  double ahp_go1[N_GO], *ahp_gr1;
  int *spikep;
  double firing_rate;
  double *psp_pf;
  double psp_ampa_vnmf[N_VN][N_MF_PER_VN], psp_nmda_vnmf[N_VN][N_MF_PER_VN];
  double *psp_bspf;
  double psp_pkj[N_PKJ], psp_vn[N_VN], psp_preio[N_IO], psp_bs[N_BS];
  int nt, t_each;

  const double decay_ahp_go = exp(-(double)DT/tau_ahp_go);
  const double decay_ahp_gr = exp(-(double)DT/tau_ahp_gr);
  const double decay1_gaba_grgo = exp(-(double)DT/tau1_gaba_grgo);
  const double decay2_gaba_grgo = exp(-(double)DT/tau2_gaba_grgo);
  const double decay_ampa_gogr = exp(-(double)DT/tau_ampa_gogr);
  const double decay1_nmda_gogr = exp(-(double)DT/tau1_nmda_gogr);
  const double decay2_nmda_gogr = exp(-(double)DT/tau2_nmda_gogr);
  const double decay_ampa_grmf = exp(-(double)DT/tau_ampa_grmf);
  const double decay_nmda_grmf = exp(-(double)DT/tau_nmda_grmf);
  const double decay_ampa_pkjpf = exp(-(double)DT/tau_ampa_pkjpf);
  const double decay_ahp_pkj = exp(-(double)DT/tau_ahp_pkj);
  const double decay_ampa_vnmf = exp(-DT/tau_ampa_vnmf);
  const double decay_nmda_vnmf = exp(-DT/tau_nmda_vnmf);
  const double decay_gaba_vnpkj = exp(-DT/tau_gaba_vnpkj);
  const double decay_gaba_iovn = exp(-DT/tau_gaba_iovn);
  const double decay_ampa_preio = exp(-DT/tau_ampa_preio);
  const double decay_gaba_pkjbs = exp(-DT/tau_gaba_pkjbs);
  const double decay_ampa_bspf = exp(-DT/tau_ampa_bspf);
  const double decay_ahp_bs = exp(-(double)DT/tau_ahp_bs);
  const double decay_ahp_vn = exp(-DT/tau_ahp_vn);
  const double decay_ahp_io = exp(-DT/tau_ahp_io);

  int grs_in_window_tail = 0, grs_for_ltd_tail = 0;

  gsl_odeiv_step *odestep = gsl_odeiv_step_alloc(gsl_odeiv_step_rk4, N_ALL);

  mfseed = 1000;

  GR_read();
  init_weight_pkjpf();
  if (read_weight_pkjpf(0)){
    free_weight_pkjpf();
    exit(1);
  }
  init_weight_bspf();
  if (read_weight_bspf(0)){
    free_weight_bspf();
    exit(1);
  }

  gsl_odeiv_system sys = {func, NULL, N_ALL, NULL};

  spikep = (int *)malloc(N_ALL*sizeof(int));
  dudt_in = (double *)malloc(N_ALL*sizeof(double));
  dudt_out = (double *)malloc(N_ALL*sizeof(double));
  u_err = (double *)malloc(N_ALL*sizeof(double));
  psp_gr_ampa = (double *)malloc(N_GR*sizeof(double));
  psp_gr_nmda1 = (double *)malloc(N_GR*sizeof(double));
  psp_gr_nmda2 = (double *)malloc(N_GR*sizeof(double));
  psp_gr = (double *)malloc(N_GR*sizeof(double));
  ahp_gr1 = (double *)malloc(N_GR*sizeof(double));
  ahp_gr = (double *)malloc(N_GR*sizeof(double));
  g_grgo = (double *)malloc(N_GR*sizeof(double));
  nlist_grgo = (int *)malloc(N_GR*sizeof(int));
  //
  u = (double *)malloc(sizeof(double)*N_ALL);

  initialize(23L);
  
  for(i = 0; i < N_ALL; i++){
    u[i] = 0;
    dudt_in[i] = 0;
    dudt_out[i] = 0;
    u_err[i] = 0;
    spikep[i] = 0;
  }
  for(i = 0; i < N_GR; i++){
    u[GR(i)] = E_leak_GR;
  }
  for(i = 0; i < N_GO; i++){
    u[GO(i)] = E_leak_GO;
  }
  psp_ampa_grmf = (double **)malloc(N_GR*sizeof(double *));
  psp_nmda_grmf = (double **)malloc(N_GR*sizeof(double *));
  g_grmf = (double *)malloc(N_GR*sizeof(double));
  for(i = 0; i < N_GR; i++){
    psp_ampa_grmf[i] = (double *)malloc(N_MF_PER_GR*sizeof(double));
    psp_nmda_grmf[i] = (double *)malloc(N_MF_PER_GR*sizeof(double));
  }
  for(i = 0; i < N_GR; i++){
    for(j = 0; j < N_MF_PER_GR; j++){
      psp_ampa_grmf[i][j] = 0;
      psp_nmda_grmf[i][j] = 0;
    }
    g_grmf[i] = 0;
  }
  for(i = 0; i < N_GO; i++){
    psp_go_gaba1[i] = 0;
    psp_go_gaba2[i] = 0;
    psp_go_gaba3[i] = 0;
    psp_go[i] = 0;
    ahp_go1[i] = 0;
    ahp_go[i] = 0;
  }
  for(i = 0; i < N_GR; i++){
    psp_gr_ampa[i] = 0;
    psp_gr_nmda1[i] = 0;
    psp_gr_nmda2[i] = 0;
    psp_gr[i] = 0;
    ahp_gr1[i] = 0;
    ahp_gr[i] = 0;
  }

  psp_pf = (double *)malloc(N_GR*sizeof(double));
  for(i = 0; i < N_GR; i++){
    psp_pf[i] = 0;
  }
  psp_bspf = (double *)malloc(N_GR*sizeof(double));
  for(i = 0; i < N_GR; i++){
    psp_bspf[i] = 0;
  }

  for(i = 0; i < N_VN; i++){
    for(j = 0; j < N_MF_PER_VN; j++){
      psp_ampa_vnmf[i][j]  = 0;
      psp_nmda_vnmf[i][j]  = 0;
    }
  }
  for(i = 0; i < N_PKJ; i++){
    psp_pkj[i]  = 0;
  }
  for(i = 0; i < N_VN; i++){
    psp_vn[i]  = 0;
  }
  for(i = 0; i < N_IO; i++){
    psp_preio[i]  = 0;
  }
  for(i = 0; i < N_BS; i++){
    psp_bs[i]  = 0;
  }
  for(i = 0; i < N_MOL; i++){
    ex[i] = 0;
    inh[i] = 0;
    ahp[i] = 0;
  }
  for(i = 0; i < N_PKJ; i++){
    u[PKJ(i)] = E_leak_PKJ;
  }
  for(i = 0; i < N_VN; i++){
    u[VN(i)] = E_leak_VN;
  }
  for(i = 0; i < N_IO; i++){
    u[IO(i)] = E_leak_IO;
  }
  for(i = 0; i < N_BS; i++){
    u[BS(i)] = E_leak_BS;
  }

  t = 0;

  GSL_ODEIV_FN_EVAL(&sys, t, u, dudt_in);

  init_genrand(mfseed);

  for(nt = 0; nt <= n_trials; nt++){
    time0 = time(NULL);

    sprintf(filename, "gr.spk.%d", nt);
    fgr_spk = fopen(filename, "w");
    sprintf(filename, "go.spk.%d", nt);
    fgo_spk = fopen(filename, "w");
    sprintf(filename, "gr.spk.a%d", nt);
    fgr_spk_p = fopen(filename, "w");

    sprintf(filename, "pkjvnio.spk.%d", nt);
    fpkjvnio_spk = fopen(filename, "w");

    printf("==%2d==\n", nt);
    // tyam
    t_each = 0;
    grs_for_ltd_tail = 0;
    grs_for_ltd[0] = -1;
    grs_in_window_tail = 0;
    grs_in_window[grs_in_window_tail][0] = -1;

    while(t_each < oscillation_period){
      for(i = 0; i < N_GO; i++){
	psp_go_gaba1[i] = decay1_gaba_grgo*psp_go_gaba1[i] + spikep[GO(i)];
	psp_go_gaba2[i] = decay2_gaba_grgo*psp_go_gaba2[i] + spikep[GO(i)];
	psp_go[i] = (r1_gaba_grgo*psp_go_gaba1[i]+
		     r2_gaba_grgo*psp_go_gaba2[i]);
	ahp_go1[i] = decay_ahp_go*ahp_go1[i];
	if (spikep[GO(i)]){
	  ahp_go1[i] = 1;
	}
	ahp_go[i] = r_ahp_go*ahp_go1[i];

      }
      for(i = 0; i < N_GR; i++){
	psp_gr_ampa[i] = decay_ampa_gogr*psp_gr_ampa[i] + spikep[GR(i)];
	psp_gr_nmda1[i] = decay1_nmda_gogr*psp_gr_nmda1[i] + spikep[GR(i)];
	psp_gr_nmda2[i] = decay2_nmda_gogr*psp_gr_nmda2[i] + spikep[GR(i)];
	psp_gr[i] = (r_ampa_gogr*psp_gr_ampa[i]+
		     r1_nmda_gogr*psp_gr_nmda1[i]+
		     r2_nmda_gogr*psp_gr_nmda2[i]);
	ahp_gr1[i] = decay_ahp_gr*ahp_gr1[i];
	if (spikep[GR(i)]){
	  ahp_gr1[i] = 1;
	}
	ahp_gr[i] = r_ahp_gr*ahp_gr1[i];
      }
      for(i = 0; i < N_GO; i++){
	g_gogr[i] = 0;
	for(j = 0; j < nlist_gogr[i]; j++){
	  g_gogr[i] += conn_weight_gogr*psp_gr[list_gogr[i][j]];
	}
      }
      for(i = 0; i < N_GO; i++){
	r = 0;
	for(j = 0; j < nlist_grgo[i]; j++){
	  r += 2*conn_weight_grgo*psp_go[list_grgo[i][j]];
	}
	for(k = 0; k < R_N2; k++){
	  j = k+R_N2*i;
	  g_grgo[j] = r;
	}
      }

      firing_rate = MF_maxfreq*(1 - cos((2*M_PI/(double)oscillation_period)*t_each))*0.5;
      for(i = 0; i < N_GR; i++){
	for(j = 0; j < N_MF_PER_GR; j++){
	  psp_ampa_grmf[i][j] *= decay_ampa_grmf;
	  psp_nmda_grmf[i][j] *= decay_nmda_grmf;
	  if (j < 1){ // only 1 MF is stimulated. The rest are spontaneous.
	    if (genrand_real2() < firing_rate){
	      psp_ampa_grmf[i][j] += 1;
	      psp_nmda_grmf[i][j] += 1;
	    }
	  }else{
	    if (genrand_real2() < 0.005){
	      psp_ampa_grmf[i][j] += 1;
	      psp_nmda_grmf[i][j] += 1;
	    }
	  }
	}
	g_grmf[i] = 0;
	for(j = 0; j < N_MF_PER_GR; j++){
	  g_grmf[i] += lambda*(r_ampa_grmf*psp_ampa_grmf[i][j]
			       +r_nmda_grmf*psp_nmda_grmf[i][j]);
	}
      }

      for(i = 0; i < N_GR; i++){
	psp_pf[i] = decay_ampa_pkjpf*psp_pf[i] + spikep[GR(i)];
      }
      for(i = 0; i < N_BS; i++){
	psp_bs[i] = decay_gaba_pkjbs*psp_bs[i] + spikep[BS(i)];
      }
      for(i = 0; i < N_PKJ; i++){
	ex[MPKJ(i)] = 0;
	for(j = 0; j < N_GR; j++){
	  ex[MPKJ(i)] += conn_weight_pkjpf*w_pkjpf[i][j]*r_ampa_pkjpf*psp_pf[j];
	}
	inh[MPKJ(i)] = 0;

	for(j = i; j < i+3; j++){
	  if (j < N_BS){
	    k = j;
	  }else{
	    k = j - N_BS;
	  }
	  inh[MPKJ(i)] += N_BS*r_gaba_pkjbs*psp_bs[k]/3.0;
	}
	if (spikep[PKJ(i)]){
	  ahp[MPKJ(i)] = 1.0;
	}else{
	  ahp[MPKJ(i)] *= decay_ahp_pkj;
	}
      }
      for(i = 0; i < N_GR; i++){
	psp_bspf[i] = decay_ampa_bspf*psp_bspf[i] + spikep[GR(i)];
      }
      for(i = 0; i < N_BS; i++){
	ex[MBS(i)] = 0;
	for(j = 0; j < N_GR; j++){
	  ex[MBS(i)] += conn_weight_bspf*w_bspf[i][j]*r_ampa_bspf*psp_bspf[j]/N_GR;
	}
	if (spikep[BS(i)]){
	  ahp[MBS(i)] = 1.0;
	}else{
	  ahp[MBS(i)] *= decay_ahp_bs;
	}
      }
      firing_rate = MF_maxfreq*(1 - cos((2*M_PI/(double)oscillation_period)*t_each))*0.5;
      for(i = 0; i < N_VN; i++){
	for(j = 0; j < N_MF_PER_VN; j++){
	  psp_ampa_vnmf[i][j] *= decay_ampa_vnmf;
	  psp_nmda_vnmf[i][j] *= decay_nmda_vnmf;
	  if (genrand_real2() < firing_rate){
	    psp_ampa_vnmf[i][j] += 1;
	    psp_nmda_vnmf[i][j] += 1;
	  }
	}
      }
      for(i = 0; i < N_PKJ; i++){
	psp_pkj[i] = spikep[PKJ(i)]+decay_gaba_vnpkj*psp_pkj[i];
      }
      for(i = 0; i < N_VN; i++){
	ex[MVN(i)] = 0;
	for(j = 0; j < N_MF_PER_VN; j++){
	
	  ex[MVN(i)] += conn_weight_vnmf*
	    (r_ampa_vnmf*psp_ampa_vnmf[i][j]+
	     r_nmda_vnmf*psp_nmda_vnmf[i][j]);
	}
	inh[MVN(i)] = 0;
	for(j = 0; j < N_PKJ; j++){
	  inh[MVN(i)] += conn_weight_vnpkj*r_gaba_vnpkj*psp_pkj[j];
	}
	if (spikep[VN(i)]){
	  ahp[MVN(i)] = 1;
	}else{
	  ahp[MVN(i)] *= decay_ahp_vn;
	}
      }
      for(i = 0; i < N_VN; i++){
	psp_vn[i] = spikep[VN(i)]+decay_gaba_iovn*psp_vn[i];
      }
      firing_rate = CF_maxfreq*(1 - cos((2*M_PI/(double)oscillation_period)*t_each))*0.5;
      for(i = 0; i < N_IO; i++){
	psp_preio[i] *= decay_ampa_preio;
	if (genrand_real2() < firing_rate){
	  psp_preio[i] += 1;
	}
	ex[MIO(i)] = kappa_preio*r_ampa_preio*psp_preio[i];
	inh[MIO(i)] = 0;
	for(j = 0; j < N_VN; j++){
	  inh[MIO(i)] += conn_weight_iovn*r_gaba_iovn*psp_vn[j];
	}
	if (spikep[IO(i)]){
	  ahp[MIO(i)] = 1;
	}else{
	  ahp[MIO(i)] *= decay_ahp_io;
	}
      }

      int status = gsl_odeiv_step_apply(odestep, (double)t, DT, u, u_err,
					dudt_in, dudt_out, &sys);
      if (status != GSL_SUCCESS){
	break;
      }
      for(i = 0; i < N_ALL; i++){
	dudt_in[i] = dudt_out[i];
      }

      for(i = 0; i < N_IO; i++){
	if (u[IO(i)] > TH_IO){
	  fprintf(fpkjvnio_spk, "%d %d\n", t_each, MIO(i));
	  spikep[IO(i)] = 1;
	  //
	  for(j = 0; j < ltd_window_size; j++){
	    for(k = 0; grs_in_window[j][k] != -1; k++){
	      grs_for_ltd[grs_for_ltd_tail] = grs_in_window[j][k];
	      grs_for_ltd_tail++;
	      grs_for_ltd[grs_for_ltd_tail] = -1;
	    }
	  }
	}else{
	  spikep[IO(i)] = 0;
	}
      }

      nspikegr = 0;
      for(i = 0; i < N_GR; i++){
	if (u[GR(i)] > TH_GR){
	  if (GR(i) % (R_N2) == 0){
	    fprintf(fgr_spk, "%d %d\n", t_each, i/(R_N2));
	  }
	  fprintf(fgr_spk_p, "%d %d\n", t_each, i);
	  if (u[GR(i)] > 0){
	    fprintf(stderr, "ERR! (%d)\n", i);
	    exit(1);
	  }
	  u[GR(i)] = E_leak_GR;
	  spikep[GR(i)] = 1;
	  grs_in_window[grs_in_window_tail][nspikegr] = i;
	  nspikegr++;
	}else{
	  spikep[GR(i)] = 0;
	}
      }
      grs_in_window[grs_in_window_tail][nspikegr] = -1;
      grs_in_window_tail += 1;
      if (ltd_window_size <= grs_in_window_tail){
	grs_in_window_tail = 0;
      }

      nspikego = 0;
      for(i = 0; i < N_GO; i++){
	if (u[GO(i)] > TH_GO){
	  fprintf(fgo_spk, "%d %d\n", t_each, i);
	  nspikego++;
	  u[GO(i)] = E_leak_GO;
	  spikep[GO(i)] = 1;
	}else{
	  spikep[GO(i)] = 0;
	}
      }
      for(i = 0; i < N_PKJ; i++){
	if (u[PKJ(i)] > TH_PKJ){
	  fprintf(fpkjvnio_spk, "%d %d\n", t_each, MPKJ(i));
	  spikep[PKJ(i)] = 1;
	  //u[PKJ(i)] = E_leak_PKJ;
	}else{
	  spikep[PKJ(i)] = 0;
	}
      }
      for(i = 0; i < N_BS; i++){
	if (u[BS(i)] > TH_BS){
	  fprintf(fpkjvnio_spk, "%d %d\n", t_each, MBS(i));
	  spikep[BS(i)] = 1;
	}else{
	  spikep[BS(i)] = 0;
	}
      }
      for(i = 0; i < N_VN; i++){
	if (u[VN(i)] > TH_VN){
	  fprintf(fpkjvnio_spk, "%d %d\n", t_each, MVN(i));
	  spikep[VN(i)] = 1;
	  //u[VN(i)] = E_leak_VN;
	}else{
	  spikep[VN(i)] = 0;
	}
      }

      fflush(fpkjvnio_spk);

      //LTP
      for(j = 0; j < N_GR; j++){
	if (spikep[GR(j)]){
	  for(i = 0; i < N_PKJ; i++){
	    if (w_pkjpf[i][j] > 0){
	      w_pkjpf[i][j] = w_pkjpf[i][j] + c_ltp*(1-w_pkjpf[i][j]);
	    }
	  }
	}
      }
      fflush(fgr_spk);
      fflush(fgr_spk_p);
      fflush(fgo_spk);
      t += DT;
      t_each++;
    }
    time1 = time(NULL);
    fprintf(stderr, "time: %f\n", difftime(time1, time0));

    LTD();
    if (nt == n_trials-1){
      write_weight_pkjpf(nt);
    }

    fclose(fgr_spk);
    fclose(fgr_spk_p);
    fclose(fgo_spk);
    fclose(fpkjvnio_spk);
  }

  gsl_odeiv_step_free(odestep);

  free(u);

  free_wlist();
  for(i = 0; i < N_GR; i++){
    free(psp_ampa_grmf[i]);
    free(psp_nmda_grmf[i]);
  }
  free(psp_ampa_grmf);
  free(psp_nmda_grmf);
  free(g_grmf);	
  free(spikep);
  free(dudt_in);
  free(dudt_out);
  free(u_err);
  free(psp_gr_ampa);
  free(psp_gr_nmda1);
  free(psp_gr_nmda2);
  free(psp_gr);
  free(ahp_gr1);
  free(ahp_gr);
  free(g_grgo);
  free(nlist_grgo);

  GR_free();
  free_weight_pkjpf();
  free_weight_bspf();
  free(psp_pf);
  free(psp_bspf);

  return 0;
}
