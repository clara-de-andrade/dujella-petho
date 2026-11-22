#include <locale.h>

#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>
#include <mpfr.h>

#include "polynomials.h"
#include "dujella_petho.h"
#include "matveev.h"

/* ************************************************************************* */

#define m_min 0
#define k_min 3
#define n_min (k + 2)

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
    fflush(stdout);
}

/* ************************************************************************* */

#define prec (1 << 12)
#define k_max 500

int main(void)
{
#ifdef _WIN32
    setlocale(LC_ALL, ".UTF8");
#else
    setlocale(LC_ALL, "");
#endif

    mpfr_set_default_prec(prec);

    puts("Calculando alpha e g...");

    mpfr_t *alpha = calloc(k_max + 1, sizeof(mpfr_t));
    mpfr_t *g = calloc(k_max + 1, sizeof(mpfr_t));
    {
        mpfr_t two;
        mpfr_init_set_ui(two, 2, MPFR_RNDN);

        poly_t P;
        poly_init(P);

        mpfr_t temp;
        mpfr_init(temp);

        for (size_t k = k_min; k <= k_max; k++)
        {
            for (size_t i = 0; i <= k; i++)
                poly_set_coeff_si(P, (i == k) ? 1 : -1, i);

            mpfr_init(alpha[k]);
            poly_root(alpha[k], P, two);

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

    puts("Calculando M0 e N0...");

    mat_t *mat = calloc(k_max + 1, sizeof(mat_t));
    mpz_t *M0 = calloc(k_max + 1, sizeof(mpz_t));
    mpz_t *N0 = calloc(k_max + 1, sizeof(mpz_t));
    {
        mpfr_t temp;
        mpfr_init(temp);

        for (size_t k = k_min; k <= k_max; k++)
        {
            mat_init(mat[k], prec, 3);
            mpfr_set_ui(mat[k]->alpha[0], 2, MPFR_RNDU);
            mpfr_set(mat[k]->alpha[1], g[k], MPFR_RNDU);
            mpfr_set(mat[k]->alpha[2], alpha[k], MPFR_RNDU);

            mat[k]->degree = k;
            mpfr_log_ui(temp, 2, MPFR_RNDN);
            mpfr_mul_ui(mat[k]->A[0], temp, k, MPFR_RNDU);
            mpfr_log_ui(mat[k]->A[1], 2 * k, MPFR_RNDU);
            mpfr_mul_ui(mat[k]->A[1], mat[k]->A[1], k, MPFR_RNDU);
            mpfr_set(mat[k]->A[2], temp, MPFR_RNDU);

            matveev_C(mat[k]);

            mpz_init(M0[k]);
            mpfr_div(temp, mat[k]->C, temp, MPFR_RNDU);
            mpfr_get_z(M0[k], temp, MPFR_RNDU);
            mpz_pow_ui(M0[k], M0[k], 2);

            mpz_init(N0[k]);
            mpz_mul_ui(N0[k], M0[k], 2);
            mpz_add_ui(N0[k], N0[k], 2);

            progbar(k, k_min, k_max);
        }
        puts("\nPRONTO!\n");

        mpfr_clear(temp);
    }

    puts("Calculando eps e q para N0...");

    dp_t *dp1 = calloc(k_max + 1, sizeof(dp_t));
    mpz_t *M1 = calloc(k_max + 1, sizeof(mpz_t));
    mpz_t *N1 = calloc(k_max + 1, sizeof(mpz_t));
    {
        mpfr_t temp, log_2;
        mpfr_inits(temp, log_2, NULL);
        mpfr_log_ui(log_2, 2, MPFR_RNDN);

        for (size_t k = k_min; k <= k_max; k++)
        {
            dp_init(dp1[k], prec);

            mpfr_ui_div(dp1[k]->A, 1, log_2, MPFR_RNDN);
            mpfr_set_ui(dp1[k]->B, 2, MPFR_RNDN);
            mpz_set(dp1[k]->M, N0[k]);

            mpfr_log(dp1[k]->gamma, alpha[k], MPFR_RNDN);
            mpfr_div(dp1[k]->gamma, dp1[k]->gamma, log_2, MPFR_RNDN);

            mpfr_log(dp1[k]->mu, g[k], MPFR_RNDN);
            mpfr_div(dp1[k]->mu, dp1[k]->mu, log_2, MPFR_RNDN);

            dujella_petho(dp1[k]);

            mpz_init_set(M1[k], dp1[k]->N);
            mpz_init(N1[k]);
            mpz_mul_ui(N1[k], M1[k], 2);
            mpz_add_ui(N1[k], N1[k], 2);

            progbar(k, k_min, k_max);
        }
        puts("\nPRONTO!\n");

        mpfr_clears(temp, log_2, NULL);
    }

    puts("Calculando eps e q para N1...");

    dp_t *dp2 = calloc(k_max + 1, sizeof(dp_t));
    mpz_t *M2 = calloc(k_max + 1, sizeof(mpz_t));
    mpz_t *N2 = calloc(k_max + 1, sizeof(mpz_t));
    for (size_t k = k_min; k <= k_max; k++)
    {
        dp_init(dp2[k], prec);

        mpfr_set(dp2[k]->A, dp1[k]->A, MPFR_RNDN);
        mpfr_set(dp2[k]->B, dp1[k]->B, MPFR_RNDN);
        mpfr_set(dp2[k]->gamma, dp1[k]->gamma, MPFR_RNDN);
        mpfr_set(dp2[k]->mu, dp1[k]->mu, MPFR_RNDN);

        mpz_set(dp2[k]->M, N1[k]);

        dujella_petho(dp2[k]);

        mpz_init_set(M2[k], dp2[k]->N);
        mpz_init(N2[k]);
        mpz_mul_ui(N2[k], M2[k], 2);
        mpz_add_ui(N2[k], N2[k], 2);

        progbar(k, k_min, k_max);
    }
    puts("\nPRONTO!\n");

    puts("Testando F(k)_n = 2**m para m < M2, n < N2, ...");
    for (size_t k = k_min; k <= k_max; k++)
    {
        int flag = 0;

        mpz_t *F = calloc(k, sizeof(mpz_t));
        for (size_t i = 0; i < k; i++)
        {
            mpz_init_set_ui(F[i], 1);
            mpz_mul_2exp(F[i], F[i], i);
        }

        mpz_t temp;
        mpz_init(temp);

        size_t m_max = (size_t)mpz_get_ui(M2[k]);
        size_t n_max = (size_t)mpz_get_ui(N2[k]);

        for (size_t n = n_min; n <= n_max; n++)
        {
            mpz_set_ui(temp, 0);
            for (size_t i = 0; i < k; i++)
                mpz_add(temp, temp, F[i]);
            for (size_t i = 0; i < k; i++)
                mpz_set(F[i], (i < k - 1) ? F[i + 1] : temp);

            for (size_t m = m_min; m <= m_max; m++)
            {
                mpz_ui_pow_ui(temp, 2, m);
                if (mpz_cmp(F[k - 1], temp) == 0)
                {
                    printf(
                        "\nSolução encontrada: (k, m, n) = (%lu, %lu, %lu)\n",
                        k,
                        m,
                        n);
                    flag = 1;
                    break;
                }
            }

            if (flag)
                break;
        }

        mpz_clear(temp);
        for (size_t i = 0; i < k; i++)
            mpz_clear(F[i]);
        free(F);

        if (flag)
            break;

        progbar(k, k_min, k_max);
    }
    puts("\nPRONTO!\n");

    puts("Escrevendo dados em test_3.csv...");
    {
        FILE *fp = fopen("test_3.csv", "w");
        if (fp == NULL)
            exit(EXIT_FAILURE);

        mpfr_t M0_fr, N0_fr, q_fr, M1_fr, N1_fr, M2_fr, N2_fr;
        mpfr_inits2(prec, M0_fr, N0_fr, q_fr, M1_fr, N1_fr, M2_fr, N2_fr, (mpfr_ptr)NULL);

        fprintf(fp, "k,alpha_k,g_k,C_k,M0_k,N0_k,eps_k,q_k,M1_k,N1_k,M2_k,N2_k\n");

        for (size_t k = k_min; k <= k_max; k++)
        {
            mpfr_set_z(M0_fr, M0[k], MPFR_RNDN);
            mpfr_set_z(N0_fr, N0[k], MPFR_RNDN);
            mpfr_set_z(q_fr, dp1[k]->q, MPFR_RNDN);
            mpfr_set_z(M1_fr, M1[k], MPFR_RNDN);
            mpfr_set_z(N1_fr, N1[k], MPFR_RNDN);
            mpfr_set_z(M2_fr, M2[k], MPFR_RNDN);
            mpfr_set_z(N2_fr, N2[k], MPFR_RNDN);

            fprintf(fp, "%zu,", k);
            mpfr_fprintf(fp, "%.50Rg,%.50Rg,%.49Re,%.49Re,%.49Re,%.49Re,%.49Re,%.49Re,%.49Re,%.49Re,%.49Re\n",
                         alpha[k],
                         g[k],
                         mat[k]->C,
                         M0_fr,
                         N0_fr,
                         dp1[k]->eps,
                         q_fr,
                         M1_fr,
                         N1_fr,
                         M2_fr,
                         N2_fr);

            progbar(k, k_min, k_max);
        }
        puts("\nPRONTO!\n");

        mpfr_clears(M0_fr, N0_fr, q_fr, M1_fr, N1_fr, M2_fr, N2_fr, (mpfr_ptr)NULL);

        fclose(fp);
    }

    for (size_t k = k_min; k <= k_max; k++)
    {
        mpfr_clears(alpha[k], g[k], NULL);
        mpz_clears(M0[k], M1[k], M2[k], N0[k], N1[k], N2[k], NULL);
        mat_clear(mat[k]);
        dp_clear(dp1[k]);
        dp_clear(dp2[k]);
    }
    free(alpha);
    free(g);
    free(mat);
    free(M0);
    free(N0);
    free(dp1);
    free(M1);
    free(N1);
    free(dp2);
    free(M2);
    free(N2);

    return 0;
}

/* ************************************************************************* */