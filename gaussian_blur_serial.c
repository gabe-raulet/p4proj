#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

struct pgm_image {
    int nrows;
    int ncols;
    int max_pixel;
    uint8_t *vals;
};

void pgm_image_read(struct pgm_image *pgm, const char *filename);
void pgm_image_write(struct pgm_image *pgm, const char *filename);

int main(int argc, char *argv[])
{
    if (argc != 4) {
        fprintf(stderr, "Usage: ./gaussian_blur_serial <input_pgm> output_pgm> <sigma>\n");
        return 1;
    }

    struct pgm_image pgm;
    pgm_image_read(&pgm, argv[1]);
    pgm_image_write(&pgm, argv[2]);

    return 0;
}

void pgm_image_read(struct pgm_image *pgm, const char *filename)
{
    char header[128], n1[32], n2[32];
    FILE *fp = fopen(filename, "r");
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

    pgm->nrows = atoi(n1);
    pgm->ncols = atoi(n2);

    fgets(header, 128, fp);
    eol = strchr(header, '\n');
    *eol = '\0';
    pgm->max_pixel = atoi(header);

    pgm->vals = malloc(pgm->nrows * pgm->ncols);
    fread(pgm->vals, 1, pgm->nrows * pgm->ncols, fp);
    fclose(fp);
}

void pgm_image_write(struct pgm_image *pgm, const char *filename)
{
    FILE *fp = fopen(filename, "w");
    fprintf(fp, "P5\n%d %d\n%d\n", pgm->nrows, pgm->ncols, pgm->max_pixel);
    fwrite(pgm->vals, 1, pgm->nrows * pgm->ncols, fp);
    fclose(fp);
}

