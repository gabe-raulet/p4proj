all: gaussian_blur_serial gaussian_blur_cuda

gaussian_blur_serial: gaussian_blur_serial.c
	gcc -O2 -g -o $@ $^ -lm

gaussian_blur_cuda: gaussian_blur_cuda.cu
	nvcc -O2 -g -o $@ $^ -lm
