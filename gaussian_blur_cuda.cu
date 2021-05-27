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

__global__ void blur_kernel_cuda_pass_one(float *H, float *K, uint8_t *I, int k, int nrows, int ncols)
{
    int m = blockIdx.x * blockDim.x + threadIdx.x;
    int n = blockIdx.y * blockDim.y + threadIdx.y;

    if (m >= nrows || n >= ncols) return;

    float val = 0;
    for (int j = -k; j <= k; ++j)
        val += K[j+k]*I[m*ncols + vidx(n, j, ncols)];

    H[m*ncols + n] = val;
}

__global__ void blur_kernel_cuda_pass_two(float *H, float *K, uint8_t *J, int k, int nrows, int ncols)
{
    int m = blockIdx.x * blockDim.x + threadIdx.x;
    int n = blockIdx.y * blockDim.y + threadIdx.y;

    if (m >= nrows || n >= ncols) return;

    float val = 0;
    for (int i = -k; i <= k; ++i)
        val += K[i+k]*H[vidx(m, i, nrows)*ncols + n];

    J[m*ncols + n] = val;
}

int main(int argc, char *argv[])
{
    float sigma, *K_h, *K_d, *H_h, *H_d;
    uint8_t *I_h, *I_d, *J_h, *J_d;
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

    I_h = (uint8_t *)malloc(nrows*ncols);
    J_h = (uint8_t *)malloc(nrows*ncols);
    H_h = (float *)malloc(sizeof(float)*nrows*ncols);
    K_h = (float *)malloc(sizeof(float)*N);

    fread(I_h, 1, nrows*ncols, ifp);
    fclose(ifp);

    for (int i = 0, x = -k; i < N; ++i, ++x)
        K_h[i] = exp(-(x*x)/(2.*sigma*sigma))/(sqrt(2.*PI)*sigma);

    cudaMalloc((void **)&I_d, nrows*ncols);
    cudaMalloc((void **)&J_d, nrows*ncols);
    cudaMalloc((void **)&H_d, sizeof(float)*nrows*ncols);
    cudaMalloc((void **)&K_d, sizeof(float)*N);

    cudaMemcpy(I_d, I_h, nrows*ncols, cudaMemcpyHostToDevice);
    cudaMemcpy(K_d, K_h, sizeof(float), cudaMemcpyHostToDevice);

    dim3 block(BLOCK_DIM, BLOCK_DIM);
    dim3 grid(ceil(nrows/(BLOCK_DIM + 0.0)), ceil(ncols/(BLOCK_DIM + 0.0)));

    blur_kernel_cuda_pass_one<<<grid, block>>>(H_d, K_d, I_d, k_d, nrows, ncols);
    blur_kernel_cuda_pass_two<<<grid, block>>>(H_d, K_d, J_d, k_d, nrows, ncols);

    cudaMemcpy(J_h, J_d, nrows*ncols, cudaMemcpyDeviceToHost);

    cudaFree(I_d);
    cudaFree(J_d);
    cudaFree(H_d);
    cudaFree(K_d);

    ofp = fopen(oname, "w");
    fprintf(ofp, "P5\n%d %d\n%d\n", ncols, nrows, max_pixel);
    fwrite(J_h, 1, nrows*ncols, ofp);
    fclose(ofp);

    free(H_h);
    free(I_h);
    free(J_h);
    free(K_h);
    return 0;
}

