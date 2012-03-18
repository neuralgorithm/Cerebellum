#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gd.h>
#include<assert.h>
#include<string.h>

#define sigma	8.3
#define R_N2	(100)
#define N	(32*32)
#define T	2000

#define nGrayTcale	256			// grayscale

int main(int argc, char *argv[])
{
  FILE *file;
  char buf[1024];
  int t, s, t1, t2, i, n;
  double **s1, **s2;
  double **spike1, **spike2;
  double r, norm1[T], norm2[T], **c;
  double texp[T];

  gdImagePtr im;
  int gray[nGrayTcale];

  if (argc < 4){
    fprintf(stderr, "usage %s <input> <input2> <output>\n", argv[0]);
    exit(1);
  }

  s1 = (double **)malloc(T*sizeof(double *));
  s2 = (double **)malloc(T*sizeof(double *));
  spike1 = (double **)malloc(T*sizeof(double *));
  spike2 = (double **)malloc(T*sizeof(double *));
  for(t = 0; t < T; t++){
    s1[t] = (double *)malloc(N*sizeof(double));
    s2[t] = (double *)malloc(N*sizeof(double));
    spike1[t] = (double *)malloc(N*sizeof(double));
    spike2[t] = (double *)malloc(N*sizeof(double));
    for(i = 0; i < N; i++){
      s1[t][i] = 0;
      s2[t][i] = 0;
      spike1[t][i] = 0;
      spike2[t][i] = 0;
    }
  }

  file = fopen(argv[1], "r");
  assert(file != NULL);
  while((fgets(buf, 1024, file))){
    sscanf(buf, "%d %d", &t, &i);
    n = i / R_N2;
    s1[t][n] += 1.0/R_N2;
  }
  fclose(file);

  file = fopen(argv[2], "r");
  assert(file != NULL);
  while((fgets(buf, 1024, file))){
    sscanf(buf, "%d %d", &t, &i);
    n = i / R_N2;
    s2[t][n] += 1.0/R_N2;
  }
  fclose(file);

  for(t = 0; t < T; t++){
    texp[t] = exp(-t/sigma);
  }

  for(t = 0; t < T; t++){
    for(i = 0; i < N; i++){
      r = 0;
      for(s = 0; s < t; s++){
	r += texp[t-s]*s1[s][i];
      }
      spike1[t][i] = r;
    }
  }
  for(t = 0; t < T; t++){
    for(i = 0; i < N; i++){
      r = 0;
      for(s = 0; s < t; s++){
	r += texp[t-s]*s2[s][i];
      }
      spike2[t][i] = r;
    }
  }

  for(t = 0; t < T; t++){
    r = 0;
    for(i = 0; i < N; i++){
      r += spike1[t][i]*spike1[t][i];
    }
    norm1[t] = sqrt(r);
  }
  for(t = 0; t < T; t++){
    r = 0;
    for(i = 0; i < N; i++){
      r += spike2[t][i]*spike2[t][i];
    }
    norm2[t] = sqrt(r);
  }
  
  c = (double **)malloc(T*sizeof(double *));
  for(t1 = 0; t1 < T; t1++){
    c[t1] = (double *)malloc(T*sizeof(double));
  }
  for(t1 = 0; t1 < T; t1++){
    for(t2 = 0; t2 < T; t2++){
      r = 0;
      for(i = 0; i < N; i++){
	r += spike1[t1][i]*spike2[t2][i];
      }
      if (norm1[t1] > 0 && norm2[t2] > 0){
	c[t1][t2] = r/(norm1[t1]*norm2[t2]);
      }else{
	c[t1][t2] = 0;
      }
    }
  }
  im = gdImageCreate(T, T);
  for(i = 0; i < nGrayTcale; i++){
    gray[i] = gdImageColorAllocate(im, i, i, i);
  }
  for(t1 = 0; t1 < T; t1++){
    for(t2 = 0; t2 < T; t2++){
      i = floor((nGrayTcale-1)*c[t1][t2]);
      gdImageSetPixel(im, t1, t2, gray[i]);
    }
  }
  sprintf(buf, "%s.png", argv[3]);
  file = fopen(buf, "wb");
  gdImagePng(im, file);
  fclose(file);
  gdImageDestroy(im);

  if (strcmp(argv[1], argv[2]) != 0){
    sprintf(buf, "%s.d", argv[3]);
    file = fopen(buf, "w");
    for(t1 = 0; t1 < T; t1++){
      fprintf(file, "%d %f\n", t1, c[t1][t1]);
    }
    fclose(file);
  }

  for(t = 0; t < T; t++){
    free(s1[t]);
    free(s2[t]);
    free(spike1[t]);
    free(spike2[t]);
  }
  free(s1);
  free(s2);
  free(spike1);
  free(spike2);
  for(t = 0; t < T; t++){
    free(c[t]);
  }
  free(c);
  return 0;
}
