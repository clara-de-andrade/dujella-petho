#include <locale.h>

#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>
#include <mpfr.h>

#include "dptest.h"

/* ************************************************************************* */

#define l_min  2
#define k_min  2
#define q_min  3

#define a_min  1
#define a_max  9

/* ************************************************************************* */

void progbar(size_t i, size_t i_min, size_t i_max)
{
    int progression = 100 * (i - i_min + 1) / (i_max - i_min + 1);
    int completed = progression / 2;
    int remaining = (100 - progression + 1) / 2;

    putchar('[');
    for (int i = 0; i < completed; i++)
        putchar('#');
    for (int i = 0; i < remaining; i++)
        putchar(' ');
    putchar(']');

    printf(" %i%%   \r", progression);
}

/* ************************************************************************* */

#define prec   (1 << 12)
#define q_max  10000

int main(void)
{
    #ifdef _WIN32
        setlocale(LC_ALL, ".UTF8");
    #else
        setlocale(LC_ALL, "");
    #endif

    mpfr_set_default_prec(prec);

    puts("Calculando cotas iniciais...");

    mpz_t* l_max = calloc(q_max + 1, sizeof(mpz_t));
    {
        mpfr_t c, l;
        mpfr_inits(c, l, NULL);
        mpfr_set_str(c, "1.2e95", 10, MPFR_RNDN);
    
        for (size_t q = q_min; q <= q_max; q++)
        {
            mpfr_log_ui(l, q + 1, MPFR_RNDN);
            mpfr_pow_ui(l, l, 21, MPFR_RNDN);
            mpfr_mul(l, l, c, MPFR_RNDN);
            mpfr_get_z(l_max[q], l, MPFR_RNDD);

            progbar(q, q_min, q_max);
        }
        puts("\nPRONTO!\n");
    
        mpfr_clears(c, l, NULL);
    }

    puts("Calculando alpha e Gamma...");

    mpfr_t* alpha = calloc(q_max + 1, sizeof(mpfr_t));
    mpfr_t* Gamma = calloc(q_max + 1, sizeof(mpfr_t));
    {
        mpfr_t D;
        mpfr_init(D);
        
        for (size_t q = q_min; q <= q_max; q++)
        {
            mpfr_sqrt_ui(D, (q - 1) * (q - 1) + 4, MPFR_RNDN);
            
            mpfr_init(alpha[q]);
            mpfr_add_ui(alpha[q], D, q + 1, MPFR_RNDN);
            mpfr_div_ui(alpha[q], alpha[q], 2, MPFR_RNDN);

            mpfr_init_set_ui(Gamma[q], q - 1, MPFR_RNDN);
            mpfr_div(Gamma[q], Gamma[q], D, MPFR_RNDN);
            mpfr_add_ui(Gamma[q], Gamma[q], 1, MPFR_RNDN);

            progbar(q, q_min, q_max);
        }
        puts("\nPRONTO!\n");

        mpfr_clear(D);
    }

    puts("Calculando epsilon, q e N...");

    dp_t** dp = calloc(q_max + 1, sizeof(dp_t*));
    { 
        mpfr_t log_10;
        mpfr_init(log_10);
        mpfr_log_ui(log_10, 10, MPFR_RNDN);

        mpfr_t A, B, gamma, log_alpha;
        mpfr_inits(A, B, gamma, log_alpha, NULL);

        for (size_t q = q_min; q <= q_max; q++)
        {
            dp[q] = calloc(a_max + 1, sizeof(dp_t));

            mpfr_inits(A, B, NULL);
            mpfr_set_ui(A, 2 * q, MPFR_RNDN);
            mpfr_sqrt_ui(B, q, MPFR_RNDN);

            mpfr_log(log_alpha, alpha[q], MPFR_RNDN);
            mpfr_div(gamma, log_10, log_alpha, MPFR_RNDN);

            for (size_t a = a_min; a <= a_max; a++)
            {
                dp_init(dp[q][a], prec);
                mpfr_set(dp[q][a]->A, A, MPFR_RNDN);
                mpfr_set(dp[q][a]->B, B, MPFR_RNDN);
                mpfr_set(dp[q][a]->gamma, gamma, MPFR_RNDN);
                mpz_set(dp[q][a]->M, l_max[q]);

                mpfr_mul_ui(dp[q][a]->mu, Gamma[q], 9, MPFR_RNDN);
                mpfr_ui_div(dp[q][a]->mu, a, dp[q][a]->mu, MPFR_RNDN);
                mpfr_log(dp[q][a]->mu, dp[q][a]->mu, MPFR_RNDN);

                dptest(dp[q][a]);
            }

            progbar(q, q_min, q_max);
        }
        puts("\nPRONTO!\n");

        mpfr_clears(A, B, gamma, log_10, log_alpha, NULL);
    }
        
    puts("Escrevendo dados em test_2.csv...");
    {
        FILE* fp = fopen("test_2.csv", "w");
        if (fp == NULL)
            exit(EXIT_FAILURE);

        mpfr_t M_fr, q_fr;
        mpfr_inits2(256, M_fr, q_fr, (mpfr_ptr) NULL);

        fprintf(fp, "q,a,M,A,B,gamma,mu,epsilon,N,q_conv\n");

        for (size_t q = q_min; q <= q_max; q++)
        {
            for (size_t a = a_min; a <= a_max; a++)
            {
                mpfr_set_z(M_fr, dp[q][a]->M, MPFR_RNDN);
                mpfr_set_z(q_fr, dp[q][a]->q, MPFR_RNDN);

                fprintf(fp, "%zu,%zu,", q, a);
                mpfr_fprintf(fp, "%.49Re,%.50Rg,%.50Rg,%.50Rg,%.50Rg,%.49Re,",
                             M_fr,
                             dp[q][a]->A,
                             dp[q][a]->B,
                             dp[q][a]->gamma,
                             dp[q][a]->mu,
                             dp[q][a]->epsilon);
                gmp_fprintf(fp, "%Zd,", dp[q][a]->N);
                mpfr_fprintf(fp, "%.49Re\n", q_fr);
            }

            progbar(q, q_min, q_max);
        }
        puts("\nPRONTO!\n");

        mpfr_clears(M_fr, q_fr, NULL);
        fclose(fp);
    }

    for (size_t q = q_min; q <= q_max; q++)
    {
        mpz_clear(l_max[q]);
        mpfr_clears(alpha[q], Gamma[q], NULL);
        for (size_t a = a_min; a <= a_max; a++)
            dp_clear(dp[q][a]);
        free(dp[q]);
    }
    free(l_max);
    free(alpha);
    free(Gamma);
    free(dp);

    return 0;
}

/* ************************************************************************* */