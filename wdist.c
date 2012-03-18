#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define R_N2	(100)
#define N	((32)*(32)) // number of granule-cell clusters
#define T      	2000
#define DT	1.0
#define T_ALL	((int)((T)/(DT)))
#define X_PKJ	16
#define Y_PKJ	1
#define N_PKJ	((X_PKJ)*(Y_PKJ))
#define X_GR	(10*32)
#define Y_GR    (10*32)
#define N_GR	((X_GR)*(Y_GR))

double **w_pkjpf, *ex;

void init_w(void)
{
  int i;
  w_pkjpf = (double **)malloc(N_PKJ*sizeof(double *));

  for(i = 0; i < N_PKJ; i++){
    w_pkjpf[i] = (double *)malloc(N_GR*sizeof(double));
  }
}
void free_w(void)
{
  int i;

  for(i = 0; i < N_PKJ; i++){
    free(w_pkjpf[i]);
  }
  free(w_pkjpf);
}
int read_w(char *filename)
{
  FILE *file;
  int i, j;
  char buf[1024];

  file = fopen(filename, "r");
  if (!file){
    fprintf(stderr, "cannot open %s\n", filename);
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

void print_w(char *in, char *out)
{
  FILE *fin, *fout;
  char buf[1024], s[1024];
  int i, t, j;

  fin = fopen(in, "r");
  if (!fin){
    fprintf(stderr, "cannot open %s\n", in);
  }
  fout = fopen(out, "w");
  while(fgets(buf, 1024, fin)){
    sscanf(buf, "%s %d", s, &j);
    t = (int)(atof(s) / DT);
    if (t >= T){
      break;
    }
    if (j%100 == 0){
      for(i = 0; i < N_PKJ; i++){
	if (w_pkjpf[i][j] != 0){
	  fprintf(fout, "%d %f\n", t, w_pkjpf[i][j]);
	}
      }
    }
  }
  fclose(fin);
  fclose(fout);
}

int main(int argc, char *argv[])
{
  if (argc < 4){
    fprintf(stderr, "usage: %s <i:w.?> <i:gr.spk.a?> <o:w.dat>\n", argv[0]);
    exit(1);
  }

  init_w();
  read_w(argv[1]);
  print_w(argv[2], argv[3]);
  free_w();

  return 0;
}
