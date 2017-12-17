/*
 * Compile with:    gcc -O3 oeis.c -lm -lz
 */

/*
   Change the following settings according to your needs:
     NBR is the number of sequences to match through PSLQ
     ATLEAST is the expected number of related sequences
       (can be used for discarding results with only two
       related sequences which is unlikely to be very relevant).
     MAXNORM is the maximal value of the norm for the resulting
       vector in the PSLQ algorithm.
*/

#define NBR 4
#define ATLEAST 3
#define MAXNORM 50000

#include <stdlib.h>
#include <stdio.h>
#include <regex.h>
#include <zlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "pslq.h"

#define BUFLENGTH 8192
#define XVAL 0.3183098861837906715377675268
#define CONVERGENCE 8

int cmpfunc(const void * a, const void * b) {
    return ( *(int*)a - *(int*)b );
}

int main(int argc, char **argv) {
    regex_t regex;
    gzFile file;
    unsigned char buf[BUFLENGTH];
    unsigned long int *id;
    unsigned char **names;
    unsigned char *names_buf;
    double *coeffs;
    double *values;
    int i, j, k, m, n, o, t, u;
    void *v;
    int nbr;
    unsigned long int a[NBR];
    double x[NBR];
    double y[NBR];
    double norm;
    double r;

    if(argc != 2) {
        fprintf(stderr,"Need exactly one argument being a regex\n");
        exit(EXIT_FAILURE);
    }

    srand(time(NULL));

    fprintf(stdout,"Compiling regex \"%s\" (ignore case)\n", argv[1]);
    if (regcomp(&regex, argv[1], REG_ICASE)) {
        fprintf(stderr,"Could not compile regex\n");
        exit(EXIT_FAILURE);
    }


    /* Read "names.gz " */
    if (!(file = gzopen("names.gz", "r"))) {
        fprintf(stderr,"Could not open names.gz\n");
        exit(EXIT_FAILURE);
    }
    m = 1;
    j = BUFLENGTH - 1; /* fake next location */
    k = BUFLENGTH;
    n = 8192; id = malloc(n*sizeof(unsigned long int)); nbr = 0;
    names = malloc(n*sizeof(unsigned char *));
    o = 131072; names_buf = malloc(o*sizeof(unsigned char)); t = 0;
    while (1) {
        i = j+1;
        if ((i==k)&&(m=0)&&(k < BUFLENGTH)) break;
        for(j=i; buf[j]!='\n'; j++) {
            if(j==BUFLENGTH) {
                memmove(buf, buf+i, BUFLENGTH-i);
                j -= i; i = 0;
                k = j + gzread(file, buf+j, BUFLENGTH-j);
                if (k < BUFLENGTH) {
                    if (gzeof (file)) {
                        m = 0;
                        if(k==j) {
                            i = 0; j = 0; buf[0] = '\n'; break;
                        }
                    } else {
                        fprintf(stderr,"Error while reading names.gz\n");
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }
        if((m==0)&&(i==0)&&(j==0)) break;
        buf[j] = '\0';
        if(!regexec(&regex, buf+i+8, 0, NULL, 0)) {
            if(t + (j-i-8+1) > o) {
                o += 65536; v = names_buf;
                names_buf = realloc(names_buf,o*sizeof(unsigned char));
                for(u=0;u<nbr;u++) names[u] += names_buf - (unsigned char *)v;
            }
            memmove(names_buf + t, buf+i+8, j-i-8+1);
            names[nbr] = names_buf+t; t += j-i-8+1;
            buf[i+7] = '\0';
            sscanf(buf+i, "A%lu", id+nbr);
            nbr++;
            if(nbr==n) {
                n += 8192;
                id = realloc(id, n*sizeof(unsigned long int));
                names = realloc(names, n*sizeof(unsigned char *));
            }
        }
    }
    names_buf = realloc(names_buf, t*sizeof(unsigned char));
    gzclose(file);
    fprintf(stdout,"Found %d sequences matching the regex \"%s\"\n",
            nbr, argv[1]);

    // TODO: bubble sort


    /* Read "stripped.gz " */
    if (!(file = gzopen("stripped.gz", "r"))) {
        fprintf(stderr,"Could not open stripped.gz\n");
        exit(EXIT_FAILURE);
    }
    m = 1;
    j = BUFLENGTH - 1; /* fake next location */
    k = BUFLENGTH;
    coeffs = malloc(nbr*32*sizeof(double));
    for(i=0;i<nbr;i++) coeffs[32*i] = NAN;
    while (1) {
        i = j+1;
        if ((i==k)&&(m=0)&&(k < BUFLENGTH)) break;
        for(j=i; buf[j]!='\n'; j++) {
            if(j==BUFLENGTH) {
                memmove(buf, buf+i, BUFLENGTH-i);
                j -= i; i = 0;
                k = j + gzread(file, buf+j, BUFLENGTH-j);
                if (k < BUFLENGTH) {
                    if (gzeof (file)) {
                        m = 0;
                        if(k==j) {
                            i = 0; j = 0; buf[0] = '\n'; break;
                        }
                    } else {
                        fprintf(stderr,"Error while reading stripped.gz\n");
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }
        if((m==0)&&(i==0)&&(j==0)) break;
        buf[j] = '\0';

        buf[i+7] = '\0';
        sscanf(buf+i, "A%d", &n);
        v = bsearch(&n, id, nbr, sizeof(unsigned long int), cmpfunc);
        if(v) {
            n =  ((unsigned long int *) v) - id;
            t = i+8+1; v = (void *) buf +t;
            for(u=0;u<32;u++) {
                while(buf[t]!=',') t++;
                buf[t] = '\0';
                sscanf((char *)v, "%lf", coeffs + 32*n + u);
                if(buf[t+1]=='\0') {
                    if (u<31) coeffs[32*n] = NAN;
                    break;
                }
                t++; v = (void *) (buf + t);
            }
        }
    }
    gzclose(file);

    /* check for remaining NaNs */
    for(i=0;i<nbr;i++) {
        if(isnan(coeffs[32*i])) {
            for(j=i+1;j<nbr;j++) {
                id[j-1] = id[j];
                names[j-1] = names[j];
                for(k=0;k<32;k++) coeffs[(32*(j-1))+k] = coeffs[(32*j)+k];
            }
            nbr--; i--;
        }
    }
    fprintf(stdout,"Found %d sequences in stripped.gz", nbr);
    fprintf(stdout," containing at least 32 values\n");

    /* computing the g.f. */
    values = malloc(nbr*sizeof(double));
    for(i=0;i<nbr;i++) {
        values[i] = 0.0;
        for(j=0;j<32;j++) values[i] = values[i] * XVAL + coeffs[32*i+31-j];
    }

    /* check for convergence */
    for(i=0;i<nbr;i++) {
        if (fabs(values[i]) > CONVERGENCE) {
            for(j=i+1;j<nbr;j++) {
                id[j-1] = id[j];
                names[j-1] = names[j];
                values[j-1] = values[j];
                for(k=0;k<32;k++) coeffs[(32*(j-1))+k] = coeffs[(32*j)+k];
            }
            nbr--; i--;
        }
    }
    fprintf(stdout,"Found %d converging generating functions\n", nbr);

    id = realloc(id, nbr*sizeof(unsigned long int));
    names = realloc(names, nbr*sizeof(char *));
    coeffs = realloc(coeffs, nbr*32*sizeof(double));
    values = realloc(values, nbr*sizeof(double));

    while(1) {
        for(i=0;i<NBR;i++) {
            m = 1;
            while(m) {
                a[i] = rand() % nbr;
                m = 0;
                for(j=0;j<i;j++) { if(a[j]==a[i]) { m = 1; break; } }
            }
        }
        qsort(a, NBR, sizeof(unsigned long int), cmpfunc);
        for(i=0;i<NBR;i++) x[i] = values[a[i]];
        pslq(x, NBR, y);
        // compute the norm into 'norm'
        norm = 0.0;
        for(i=0;i<NBR;i++) norm += y[i]*y[i];
        if(norm < MAXNORM) {
            m = 0;
            for(i=0;i<NBR;i++) { if(y[i]!=0.0) m++; }
            if(m<ATLEAST) continue;
            m = 1;
            for(i=0;i<32;i++) {
                r = 0.0;
                for(j=0;j<NBR;j++) r += y[j] * coeffs[32*a[j]+i];
                if (r != 0.0) { m = 0; break; }
            }
            if(m) {
                fprintf(stdout,"\n");
                for(i=0;i<NBR;i++) {
                    if(y[i] != 0.0) {
                        fprintf(stdout,"A%06lu ", id[a[i]]);
                    }
                }
                fprintf(stdout,"   -->    ");
                for(i=0;i<NBR;i++) {
                    if(y[i] != 0.0) {
                        fprintf(stdout,"%ld ", (long int) y[i]);
                    }
                }
                fprintf(stdout,"   (%lu)\n", (unsigned long int) norm);
                for(i=0;i<NBR;i++) {
                    if(y[i] != 0.0) {
                        fprintf(stdout,"A%06lu %s\n", id[a[i]], names[a[i]]);
                    }
                }
            }
        }
    }
}
