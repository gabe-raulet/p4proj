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
#define BLOCK_DIM (32)

#define vidx(rr, r, R) (((rr) + (r) < 0) ? 0 : (((rr) + (r) >= (R)) ? (R) - 1 : (rr) + (r)))

//void blur_kernel(float *H, uint8_t *vals, uint8_t *blur, int N, int k, int nrows, int ncols, int ii, int jj)
//{
//    int xv, yv, i, j;
//    float val = 0;
//    for (i = -k; i <= k; ++i) {
//        xv = vidx(ii, i, nrows);
//        for (j = -k; j <= k; ++j) {
//            yv = vidx(jj, j, ncols);
//            val += H[(i+k)*N + (j+k)]*vals[xv*ncols + yv];
//        }
//    }
//
//    blur[ii*ncols + jj] = val;
//}

__global__ void blur_kernel_cuda(float *H, uint8_t *vals, uint8_t *blur, int N, int k, int nrows, int ncols)
{
    int ii = blockIdx.x * blockDim.x + threadIdx.x;
    int jj = blockIdx.y * blockDim.y + threadIdx.y;

    if (ii >= nrows || jj >= ncols) return;

    int xv, yv;
    float val = 0;
    for (int i = -k; i <= k; ++i) {
        xv = vidx(ii, i, nrows);
        for (int j = -k; j <= k; ++j) {
            yv = vidx(jj, j, ncols);
            val += H[(i+k)*N + (j+k)]*vals[xv*ncols + yv];
        }
    }
    blur[ii*ncols + jj] = val;
}

int main(int argc, char *argv[])
{
    float sigma, ss;
    int N, k, nrows, ncols, max_pixel;
    uint8_t *vals_d, *vals_h, *blur_d, *blur_h;

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
    vals_h = (uint8_t *) malloc(nrows * ncols);
    fread(vals_h, 1, nrows * ncols, fp1);
    fclose(fp1);

    float H[N*N];
    float *H_d;
    for (int i = 0, x = -k; i < N; ++i, ++x)
        for (int j = 0, y = -k; j < N; ++j, ++y)
            H[i*N + j] = exp(-(x*x + y*y)/(2.*ss))/(2.*PI*ss);

    blur_h = (uint8_t *) malloc(nrows * ncols);
    cudaMalloc((void **)&vals_d, nrows * ncols);
    cudaMalloc((void **)&blur_d, nrows * ncols);
    cudaMalloc((void **)&H_d, N * N * sizeof(float));

    cudaMemcpy(vals_d, vals_h, nrows * ncols, cudaMemcpyHostToDevice);
    cudaMemcpy(H_d, H, N * N * sizeof(float), cudaMemcpyHostToDevice);

    dim3 block(BLOCK_DIM, BLOCK_DIM);
    dim3 grid(ceil(nrows/(BLOCK_DIM + 0.0)), ceil(ncols/(BLOCK_DIM + 0.0)));

    blur_kernel_cuda<<<grid, block>>>(H_d, vals_d, blur_d, N, k, nrows, ncols);

    //for (int ii = 0; ii < nrows; ++ii)
    //    for (int jj = 0; jj < ncols; ++jj)
    //        blur_kernel(H, vals_h, blur_h, N, k, nrows, ncols, ii, jj);

    cudaMemcpy(blur_h, blur_d, nrows * ncols, cudaMemcpyDeviceToHost);

    cudaFree(blur_d);
    cudaFree(vals_d);
    cudaFree(H_d);

    FILE *fp2 = fopen(argv[2], "wb");
    fprintf(fp2, "P5\n%d %d\n%d\n", ncols, nrows, max_pixel);
    fwrite(blur_h, 1, nrows * ncols, fp2);
    free(blur_h);
    free(vals_h);
    fclose(fp2);
    return 0;
}



