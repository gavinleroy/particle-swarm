#include "logging.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <stdlib.h>
#include <stdio.h>

static int logging_dir_fd = -1;

void set_logging_directory(char const * logdir)
{
    if( (logging_dir_fd = open(logdir, O_RDONLY |O_DIRECTORY)) < 0)
    {
        perror("Error opening logdir");
        exit(1);
    }
}


void log_surrogate(char const * fname,
    double * lambda,
    double * p,
    double * x,
    size_t t,
    size_t dim,
    size_t popsize
)
{
    if (logging_dir_fd < 0)
    {
        return;
    }
    
    int fd = -1;
    FILE * fp = NULL;

    if ( (fd = openat(logging_dir_fd, fname, O_WRONLY | O_CREAT | O_TRUNC, 0644)) < 0)
    {
        perror("Error opening file descriptor");
        exit(1);
    }

    if ( (fp = fdopen(fd, "w")) == NULL)
    {
        perror("Error opening file pointer");
        exit(1);
    }

    fwrite(&t, 1, sizeof(size_t), fp);
    fwrite(&dim, 1, sizeof(size_t), fp);
    fwrite(&popsize, 1, sizeof(size_t), fp);

    size_t lambda_s = (t+1) * popsize;
    size_t p_s = dim + 1;
    size_t x_s = (t + 1) * popsize * dim;

    fwrite(lambda, lambda_s, sizeof(double), fp);
    fwrite(p, p_s, sizeof(double), fp);
    fwrite(x, x_s, sizeof(double), fp);

    fclose(fp);
}