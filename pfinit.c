#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define X_GO	32
#define Y_GO	32
#define N_GO	((X_GO)*(Y_GO))
#define R_N	10
#define R_N2	((R_N)*(R_N))
#define X_GR	((R_N)*X_GO)
#define Y_GR	((R_N)*Y_GO)
#define N_GR	((X_GR)*(Y_GR))

#define X_PKJ	16
#define Y_PKJ	1
#define N_PKJ	((X_PKJ)*(Y_PKJ))

double **w_pkjpf;

void init_w(void)
{
  int i;
  w_pkjpf = (double **)malloc(N_PKJ*sizeof(double *));

  for(i = 0; i < N_PKJ; i++){
    w_pkjpf[i] = (double *)malloc(N_GR*sizeof(double));
  }
}
void set_w(void)
{
  int ngr, xpkj, ypkj, npkj, axpkj, dx, xgo, ygo, ngo, n;
  int i, j;

  for(i = 0; i < N_PKJ; i++){
    for(j = 0; j < N_GR; j++){
      w_pkjpf[i][j] = 0;
    }
  }
  for(ypkj = 0; ypkj < Y_PKJ; ypkj++){
    for(xpkj = 0; xpkj < X_PKJ; xpkj++){
      npkj = ypkj + Y_PKJ*xpkj;
      axpkj = ((X_GO)/(X_PKJ))*xpkj;
      for(dx = -4; dx <= 4; dx++){
	xgo = axpkj+dx;
	if (xgo >= X_GO){
	  xgo -= X_GO;
	}
	if (xgo < 0){
	  xgo += X_GO;
	}
	for(ygo = 0; ygo < Y_GO; ygo++){
	  ngo = ygo + Y_GO*xgo;
	  for(n = 0; n < R_N2; n++){
	    ngr = n + R_N2*ngo;
	    w_pkjpf[npkj][ngr] = 1.0;
	  }
	}
      }
    }
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
void write_w(char *filename)
{
  FILE *file;
  int i, j;

  file = fopen(filename, "w");
  if (!file){
    fprintf(stderr, "cannot open %s\n", filename);
    return;
  }
  for(i = 0; i < N_PKJ; i++){
    for(j = 0; j < N_GR; j++){
      fprintf(file, "%f\n", w_pkjpf[i][j]);
    }
  }
  fclose(file);
}

int main(int argc, char *argv[])
{
  if (argc < 2){
    fprintf(stderr, "usage: %s <output>\n", argv[0]);
    exit(1);
  }
  init_w();
  set_w();
  write_w(argv[1]);
  free_w();
  exit(0);
}
