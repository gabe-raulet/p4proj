all: gaussian_blur_serial

gaussian_blur_serial: gaussian_blur_serial.c
	clang -o $@ $^

