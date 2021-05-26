#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <cuda.h>

#define PI (3.14159265358979323)

void blur_kernel_old(float *H, uint8_t *vals, uint8_t *blur, int N, int k, int nrows, int ncols);
void blur_kernel(float *H, uint8_t *vals, uint8_t *blur, int N, int k, int nrows, int ncols, int ii, int jj)

int main(int argc, char *argv[])
{
    float sigma, ss, *kernel_matrix;
    int N, k, nrows, ncols, max_pixel;
    uint8_t *vals, *blur;

    if (argc != 4) {
        fprintf(stderr, "Usage: ./gaussian_blur_serial <input_pgm> output_pgm> <sigma>\n");
        return 1;
    }

    sigma = atof(argv[3]);
    ss = sigma * sigma;
    N = ceil(6*sigma);
    if (N % 2 == 0) N++;
    k = N / 2;

    char header[128], n1[32], n2[32];
    FILE *fp1 = fopen(argv[1], "rb");
    fgets(header, 128, fp1);
    char *eol = strchr(header, '\n');
    *eol = '\0';

    if (strcmp(header, "P5")) {
        fprintf(stderr, "'%s' must be in P5 format\n", argv[1]);
        fclose(fp1);
        return 1;
    }

    fgets(header, 128, fp1);
    eol = strchr(header, '\n');
    *eol = '\0';
    char *hp = header;
    for (; *hp != ' '; ++hp);
    *hp = '\0';

    strcpy(n1, header);
    strcpy(n2, hp + 1);

    ncols = atoi(n1);
    nrows = atoi(n2);

    fgets(header, 128, fp1);
    eol = strchr(header, '\n');
    *eol = '\0';
    max_pixel = atoi(header);
    vals = malloc(nrows * ncols);
    fread(vals, 1, nrows * ncols, fp1);
    fclose(fp1);

    float H[N*N];
    for (int i = 0, x = -k; i < N; ++i, ++x)
        for (int j = 0, y = -k; j < N; ++j, ++y)
            H[i*N + j] = exp(-(x*x + y*y)/(2.*ss))/(2.*PI*ss);

    blur = malloc(nrows * ncols);

    for (int ii = 0; ii < nrows; ++ii)
        for (int jj =0; jj < ncols; ++jj)
            blur_kernel(H, vals, blur, N, k, nrows, ncols, ii, jj);

    FILE *fp2 = fopen(argv[2], "wb");
    fprintf(fp2, "P5\n%d %d\n%d\n", ncols, nrows, max_pixel);
    fwrite(blur, 1, nrows * ncols, fp2);
    fclose(fp2);
    return 0;
}

#define vidx(rr, r, R) (((rr) + (r) < 0) ? 0 : (((rr) + (r) >= (R)) ? (R) - 1 : (rr) + (r)))

void blur_kernel(float *H, uint8_t *vals, uint8_t *blur, int N, int k, int nrows, int ncols, int ii, int jj)
{
    int xv, yv, i, j;
    float val = 0;
    for (i = -k; i <= k; ++i) {
        xv = vidx(ii, i, nrows);
        for (j = -k; j <= k; ++j) {
            yv = vidx(jj, j, ncols);
            val += H[(i+k)*N + (j+k)]k*vals[xv*ncols + yv];
        }
    }

    blur[ii*ncols + jj] = val;
}
