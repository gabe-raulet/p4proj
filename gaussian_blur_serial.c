#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#define GAUSSIAN_NORMALIZER (0.15915494309189535)

#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))

static float sigma;
static float ss;
static int kernel_order;
static int kernel_offset;
static int nrows;
static int ncols;
static int max_pixel;
static uint8_t *vals;
static uint8_t *blur;
static float *kernel_matrix;

float gauss(float x, float y);

int main(int argc, char *argv[])
{
    if (argc != 4) {
        fprintf(stderr, "Usage: ./gaussian_blur_serial <input_pgm> output_pgm> <sigma>\n");
        return 1;
    }

    sigma = atof(argv[3]);
    ss = sigma * sigma;
    kernel_order = ceil(6*sigma);
    if (kernel_order % 2 == 0) kernel_order++;
    kernel_offset = kernel_order / 2;

    char header[128], n1[32], n2[32];
    FILE *fp = fopen(argv[1], "rb");
    fgets(header, 128, fp);
    char *eol = strchr(header, '\n');
    *eol = '\0';

    if (strcmp(header, "P5")) {
        fprintf(stderr, "'%s' must be in P5 format\n", argv[1]);
        fclose(fp);
        exit(1);
    }

    fgets(header, 128, fp);
    eol = strchr(header, '\n');
    *eol = '\0';
    char *hp = header;
    for (; *hp != ' '; ++hp);
    *hp = '\0';

    strcpy(n1, header);
    strcpy(n2, hp + 1);

    nrows = atoi(n1);
    ncols = atoi(n2);

    fgets(header, 128, fp);
    eol = strchr(header, '\n');
    *eol = '\0';
    max_pixel = atoi(header);
    vals = malloc(nrows * ncols);
    fread(vals, 1, nrows * ncols, fp);
    fclose(fp);

    kernel_matrix = malloc(sizeof(float)*kernel_order*kernel_order);

    for (int i = 0; i < kernel_order; ++i)
        for (int j = 0; j < kernel_order; ++j)
            kernel_matrix[i*kernel_order + j] = gauss(i-kernel_offset, j-kernel_offset);

    /* start convolution */
    blur = malloc(nrows * ncols);

    int xv, yv;
    for (int i = 0; i < nrows; ++i) {
        for (int j = 0; j < ncols; ++j) {
            double val = 0;
            for (int x = 0; x < kernel_order; ++x) {
                for (int y = 0; y < kernel_order; ++y) {
                    if (i+x-kernel_offset < 0) xv = 0;
                    else if (i+x-kernel_offset >= ncols) xv = ncols - 1;
                    else xv = i+x-kernel_offset;
                    if (j+y-kernel_offset < 0) yv = 0;
                    else if (j+y-kernel_offset >= nrows) yv = nrows - 1;
                    else yv = j+y-kernel_offset;
                    val += (vals[xv*ncols + yv]+0.) * kernel_matrix[x*kernel_order + y];
                }
            }
            blur[i*ncols + j] = roundf(val);
        }
    }
    /* end convolution */

    FILE *fp2 = fopen(argv[2], "wb");
    fprintf(fp2, "P5\n%d %d\n%d\n", nrows, ncols, max_pixel);
    fwrite(blur, 1, nrows * ncols, fp2);
    fclose(fp2);
    return 0;
}

float gauss(float x, float y)
{
    return GAUSSIAN_NORMALIZER*exp(-(x*x + y*y)/(2*ss))/ss;
}
