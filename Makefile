all: gaussian_blur_serial

gaussian_blur_serial: gaussian_blur_serial.c
	gcc -O2 -g -o $@ $^ -lm

