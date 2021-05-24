#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#define GAUSSIAN_NORMALIZER (0.15915494309189535)

static float sigma;
static int order;
static int nrows;
static int ncols;
static int max_pixel;
static uint8_t *vals;
static uint8_t *blur;
static float *kernel_matrix;

float gauss(float x, float y, float sigma)
{
    float ss = sigma*sigma;
    return (GAUSSIAN_NORMALIZER/(ss))*exp(-(x*x + y*y)/(2*ss));
}

void pgm_image_read(const char *filename);
void pgm_image_write(const char *filename);

int main(int argc, char *argv[])
{
    if (argc != 4) {
        fprintf(stderr, "Usage: ./gaussian_blur_serial <input_pgm> output_pgm> <sigma>\n");
        return 1;
    }

    pgm_image_read(argv[1]);
    sigma = atof(argv[3]);

    order = ceil(6*sigma);
    if (order % 2 == 0) order++;

    kernel_matrix = malloc(order * order);

    int offset = order / 2;
    for (int i = 0; i < order; ++i)
        for (int j = 0; j < order; ++j)
            kernel_matrix[i*order + j] = gauss(i-offset, j-offset, sigma);

    blur = calloc(1, nrows * ncols);

    for (int i = 0; i < nrows; ++i) {
        for (int j = 0; j < ncols; ++j) {
            float val = 0;
            for (int x = 0; x < order; ++x) {
                for (int y = 0; y < order; ++y) {
                    val += (kernel_matrix[x*order + y]*vals[(i+x-offset)*ncols + (j+y-offset)]);
                }
                blur[i*ncols + j] = ((uint8_t)val);
            }
        }
    }

    pgm_image_write(argv[2]);
    return 0;
}

void pgm_image_read(const char *filename)
{
    char header[128], n1[32], n2[32];
    FILE *fp = fopen(filename, "rb");
    fgets(header, 128, fp);
    char *eol = strchr(header, '\n');
    *eol = '\0';

    if (strcmp(header, "P5")) {
        fprintf(stderr, "'%s' must be in P5 format\n", filename);
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
}

void pgm_image_write(const char *filename)
{
    FILE *fp = fopen(filename, "wb");
    fprintf(fp, "P5\n%d %d\n%d\n", nrows, ncols, max_pixel);
    fwrite(blur, 1, nrows * ncols, fp);
    fclose(fp);

}
