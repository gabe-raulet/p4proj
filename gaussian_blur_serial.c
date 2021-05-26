#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#define PI (3.14159265358979323)

#define vidx(rr, r, R) (((rr) + (r) < 0) ? 0 : (((rr) + (r) >= (R)) ? (R) - 1 : (rr) + (r)))

int main(int argc, char *argv[])
{
    float sigma, *K, *H;
    uint8_t *I, *J;
    int N, k, nrows, ncols, max_pixel;
    FILE *ifp, *ofp;

    if (argc != 4) {
        fprintf(stderr, "Usage: ./gaussian_blur_serial <input_pgm> output_pgm> <sigma>\n");
        return 1;
    }

    const char *iname = argv[1];
    const char *oname = argv[2];
    sigma = atof(argv[3]);

    N = ceil(6*sigma);
    if (!(N % 2)) N++;
    k = N / 2;

    char header[128], n1[32], n2[32];

    ifp = fopen(iname, "r");
    fscanf(ifp, "P5\n%d %d\n%d\n", &ncols, &nrows, &max_pixel);

    I = (uint8_t *)malloc(nrows*ncols);
    J = (uint8_t *)malloc(nrows*ncols);
    H = (float *)malloc(sizeof(float)*nrows*ncols);
    K = (float *)malloc(sizeof(float)*N);

    fread(I, 1, nrows*ncols, ifp);
    fclose(ifp);

    for (int i = 0, x = -k; i < N; ++i, ++x)
        K[i] = exp(-(x*x)/(2.*sigma*sigma))/(sqrt(2.*PI)*sigma);

    for (int m = 0; m < nrows; ++m) {
        for (int n = 0; n < ncols; ++n) {
            float val = 0;
            for (int j = -k; j <= k; ++j) {
                val += K[j+k]*I[m*ncols + vidx(n, j, ncols)];
            }
            H[m*ncols + n] = val;
        }
    }

    for (int m = 0; m < nrows; ++m) {
        for (int n = 0; n < ncols; ++n) {
            float val = 0;
            for (int i = -k; i <= k; ++i) {
                val += K[i+k]*H[vidx(m, i, nrows)*ncols + n];
            }
            J[m*ncols + n] = val;
        }
    }

    ofp = fopen(oname, "w");
    fprintf(ofp, "P5\n%d %d\n%d\n", ncols, nrows, max_pixel);
    fwrite(J, 1, nrows*ncols, ofp);
    fclose(ofp);

    free(H);
    free(I);
    free(J);
    free(K);
    return 0;
}

