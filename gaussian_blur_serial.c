#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#define PI (3.14159265358979323)

static float sigma;
static float ss;
static int N;
static int k;
static int nrows;
static int ncols;
static int max_pixel;
static uint8_t *vals;
static uint8_t *blur;
static float *kernel_matrix;

int main(int argc, char *argv[])
{
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
    FILE *fp = fopen(argv[1], "rb");
    fgets(header, 128, fp);
    char *eol = strchr(header, '\n');
    *eol = '\0';

    if (strcmp(header, "P5")) {
        fprintf(stderr, "'%s' must be in P5 format\n", argv[1]);
        fclose(fp);
        return 1;
    }

    fgets(header, 128, fp);
    eol = strchr(header, '\n');
    *eol = '\0';
    char *hp = header;
    for (; *hp != ' '; ++hp);
    *hp = '\0';

    strcpy(n1, header);
    strcpy(n2, hp + 1);

    ncols = atoi(n1);
    nrows = atoi(n2);

    fgets(header, 128, fp);
    eol = strchr(header, '\n');
    *eol = '\0';
    max_pixel = atoi(header);
    vals = malloc(nrows * ncols);
    fread(vals, 1, nrows * ncols, fp);
    fclose(fp);


    float H[N][N];
    for (int i = 0, x = -k; i < N; ++i, ++x)
        for (int j = 0, y = -k; j < N; ++j, ++y)
            H[i][j] = exp(-(x*x + y*y)/(2.*sigma*sigma))/(2.*PI*sigma*sigma);

    /* start convolution */
    blur = malloc(nrows * ncols);

    #define vidx(rr, r, R) (((rr) + (r) < 0) ? 0 : (((rr) + (r) >= (R)) ? (R) - 1 : (rr) + (r)))

    int xv, yv;
    for (int ii = 0; ii < nrows; ++ii) {
        for (int jj = 0; jj < ncols; ++jj) {
            float val = 0;
            for (int i = -k; i <= k; ++i) {
                xv = vidx(ii, i, nrows);
                for (int j = -k; j <= k; ++j) {
                    yv = vidx(jj, j, ncols);
                    val += H[i+k][j+k]*vals[xv*ncols + yv];
                }
            }
            blur[ii*ncols + jj] = floor(val);
        }
    }

    FILE *fp2 = fopen(argv[2], "wb");
    fprintf(fp2, "P5\n%d %d\n%d\n", ncols, nrows, max_pixel);
    fwrite(blur, 1, nrows * ncols, fp2);
    fclose(fp2);
    return 0;
}

