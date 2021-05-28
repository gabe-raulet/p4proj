all: gaussian_blur_serial gaussian_blur_cuda

gaussian_blur_serial: gaussian_blur_serial.c
	gcc -Wall -Werror -O2 -o $@ $^ -lm

gaussian_blur_cuda: gaussian_blur_cuda.cu
	nvcc -Xcompiler -Wall -Xcompiler -Werror -Xcompiler -O2 -o $@ $^ -lm

clean:
	rm gaussian_blur_serial
	rm gaussian_blur_cuda

