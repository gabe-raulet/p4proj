all: gaussian_blur_serial gaussian_blur_cuda

gaussian_blur_serial: gaussian_blur_serial.c
	gcc -O2 -o $@ $^ -lm

gaussian_blur_cuda: gaussian_blur_cuda.cu
	nvcc -O2 -o $@ $^ -lm
