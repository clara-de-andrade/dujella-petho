#include <locale.h>

#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>
#include <mpfr.h>

#include "polynom.h"
#include "findroot.h"

#include "dptest.h"

/* ************************************************************************* */

#define m_min  3
#define k_min  3

/* ************************************************************************* */

#define min(a, b) ((a) <= (b) ? (a) : (b))

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
#define m_max  1394

int main(void)
{
    #ifdef _WIN32
        setlocale(LC_ALL, ".UTF8");
    #else
        setlocale(LC_ALL, "");
    #endif

    mpfr_set_default_prec(prec);

    puts("Calculando cotas iniciais...");

    mpz_t n_max;
    size_t k_max;
    {
        mpfr_t s_max;
        mpfr_init(s_max);

        mpfr_t temp;
        mpfr_init(temp);
        
        mpfr_set_ui(s_max, m_max, MPFR_RNDN);
        mpfr_log(s_max, s_max, MPFR_RNDN);
        mpfr_pow_ui(s_max, s_max, 4, MPFR_RNDN);
        mpfr_init_set_ui(temp, m_max, MPFR_RNDN);
        mpfr_pow_ui(temp, temp, 5, MPFR_RNDN);
        mpfr_mul(s_max, temp, s_max, MPFR_RNDN);
        mpfr_set_str(temp, "4.7e14", 10, MPFR_RNDN);
        mpfr_mul(s_max, temp, s_max, MPFR_RNDN);

        mpfr_clear(temp);

        mpz_init(n_max);
        mpfr_get_z(n_max, s_max, MPFR_RNDD);
        mpz_mul_ui(n_max, n_max, m_max);
        mpz_add_ui(n_max, n_max, 2);

        mpfr_log(s_max, s_max, MPFR_RNDN);
        k_max = mpfr_get_ui(s_max, MPFR_RNDD);

        mpfr_clear(s_max);
    }
    
    puts("PRONTO!\n");

    puts("Calculando F...");

    mpz_t** F = calloc(k_max + 1, sizeof(mpz_t*));
    for (size_t k = k_min; k <= k_max; k++)
    {
        F[k] = calloc(m_max + 2, sizeof(mpz_t));

        mpz_init_set_ui(F[k][0], 0);
        mpz_init_set_ui(F[k][1], 1);

        mpz_t prev[k];
        for (size_t i = 0; i < k; i++)
            mpz_init_set_ui(prev[i], (i == (k - 1)) ? 1 : 0);
        
        for (size_t m = 2; m <= m_max + 1; m++)
        {
            mpz_init_set_ui(F[k][m], 0);
            for (size_t i = 0; i < k; i++)
                mpz_add(F[k][m], F[k][m], prev[i]);

            for (size_t i = 0; i < k; i++)
                mpz_set(prev[i], (i == (k - 1)) ? F[k][m] : prev[i + 1]);
        }

        progbar(k, k_min, k_max);
    }
    puts("\nPRONTO!\n");

    puts("Calculando alpha e g...");

    mpfr_t* alpha = calloc(k_max + 1, sizeof(mpfr_t));
    mpfr_t* g = calloc(k_max + 1, sizeof(mpfr_t));
    {
        mpfr_t two;
        mpfr_init_set_ui(two, 2, MPFR_RNDN);

        poly_t P;
        poly_init(P);

        mpfr_t temp;
        mpfr_init(temp);
        
        for (size_t k = k_min; k <= k_max; k++)
        {
            poly_set_deg(P, k);

            for (size_t i = 0; i <= k; i++)
                poly_set_coeff_si(P, (i == k) ? 1 : -1, i);

            mpfr_init(alpha[k]);
            findroot(alpha[k], P, two);

            mpfr_init(g[k]);
            mpfr_sub_ui(g[k], alpha[k], 2, MPFR_RNDN);
            mpfr_mul_ui(g[k], g[k], k + 1, MPFR_RNDN);
            mpfr_add_ui(g[k], g[k], 2, MPFR_RNDN);
            mpfr_sub_ui(temp, alpha[k], 1, MPFR_RNDN);
            mpfr_div(g[k], temp, g[k], MPFR_RNDN);

            progbar(k, k_min, k_max);
        }
        puts("\nPRONTO!\n");

        mpfr_clears(two, temp, NULL);
        poly_clear(P);
    }
    

    puts("Calculando epsilon, q e N...");

    dp_t** dp = calloc(m_max + 1, sizeof(dp_t*));
    { 
        mpfr_t log_F;
        mpfr_init(log_F);

        mpfr_t A, B;
        mpfr_init_set_str(A, "2.01", 10, MPFR_RNDN);
        mpfr_init_set_str(B, "1.65", 10, MPFR_RNDN);

        for (size_t m = m_min; m <= m_max; m++)
        {
            dp[m] = calloc(min(m, k_max) + 1, sizeof(dp_t));
            for (size_t k = k_min; k <= min(m, k_max); k++)
            {
                dp_init(dp[m][k], prec);

                mpfr_set(dp[m][k]->A, A, MPFR_RNDN);
                mpfr_set(dp[m][k]->B, B, MPFR_RNDN);
                mpz_set(dp[m][k]->M, n_max);
                
                mpfr_set_z(log_F, F[k][m + 1], MPFR_RNDN);
                mpfr_log(log_F, log_F, MPFR_RNDN);

                mpfr_log(dp[m][k]->gamma, alpha[k], MPFR_RNDN);
                mpfr_div(dp[m][k]->gamma, dp[m][k]->gamma, log_F, MPFR_RNDN);

                mpfr_div(dp[m][k]->mu, alpha[k], g[k], MPFR_RNDN);
                mpfr_log(dp[m][k]->mu, dp[m][k]->mu, MPFR_RNDN);
                mpfr_div(dp[m][k]->mu, dp[m][k]->mu, log_F, MPFR_RNDN);

                dptest(dp[m][k]);
                mpz_mul_ui(dp[m][k]->N, dp[m][k]->N, m_max);
            }
            progbar(m, m_min, m_max);
        }
        puts("\nPRONTO!\n");

        mpfr_clears(log_F, A, B, NULL);
    }
        
    puts("Escrevendo dados em test_1.csv...");
    {
        FILE* fp = fopen("test_1.csv", "w");
        if (fp == NULL)
            exit(EXIT_FAILURE);

        mpfr_t M_fr, q_fr;
        mpfr_inits2(256, M_fr, q_fr, (mpfr_ptr) NULL);

        fprintf(fp, "m,k,M,A,B,gamma,mu,epsilon,N,q_conv\n");

        for (size_t m = m_min; m <= m_max; m++)
        {
            for (size_t k = k_min; k <= min(m, k_max); k++)
            {
                mpfr_set_z(M_fr, dp[m][k]->M, MPFR_RNDN);
                mpfr_set_z(q_fr, dp[m][k]->q, MPFR_RNDN);

                fprintf(fp, "%zu,%zu,", m, k);
                mpfr_fprintf(fp, "%.49Re,%.50Rg,%.50Rg,%.50Rg,%.50Rg,%.49Re,",
                             M_fr,
                             dp[m][k]->A,
                             dp[m][k]->B,
                             dp[m][k]->gamma,
                             dp[m][k]->mu,
                             dp[m][k]->epsilon);
                gmp_fprintf(fp, "%Zd,", dp[m][k]->N);
                mpfr_fprintf(fp, "%.49Re\n", q_fr);
            }

            progbar(m, m_min, m_max);
        }
        puts("\nPRONTO!\n");

        mpfr_clears(M_fr, q_fr, NULL);
        fclose(fp);
    }

    for (size_t k = k_min; k <= k_max; k++)
    {
        for (size_t m = 0; m <= m_max + 1; m++)
            mpz_clear(F[k][m]);
        free(F[k]);
        mpfr_clears(alpha[k], g[k], NULL);
    }
    free(F);
    free(alpha);
    free(g);

    for (size_t m = m_min; m <= m_max; m++)
    {
        for (size_t k = k_min; k <= min(m, k_max); k++)
            dp_clear(dp[m][k]);
        free(dp[m]);
    }
    free(dp);

    return 0;
}

/* ************************************************************************* */