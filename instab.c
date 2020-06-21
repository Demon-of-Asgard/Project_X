#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include "global.h"
#include "profile_generator.h"
#include "subroutines.h"
#include "DR_function.h"

int main(void)
{
  int xgrids, zgrids;
  int countx = 0, countz = 0, count = 0;
  double Rd1, Rd2, Hd1, Hd2;
  double umax1, umax2, pmax1, pmax2, umin1, umin2, phase1 = 0.0, phase2 = 0.0;
  long double lambda_xz = 0.0, mu_xz = 0.0;
  long double n_nu = 0.0;
  long double rho_xz = 0.0, ye_xz = 0.0;
  double x[ND1], z[ND1], tmp = 0.0;
  double xSurf = 0.0, zSurf = 0.0, xx_max = 0.0;
  struct rparams p = {1.0, 0.0};
  const size_t n = 2;
  long double xinit_tmp = 0.0;
  double x_init[2] = {0.0};

  int status, i, j, k, ii;
  FILE *fp1, *fp2, *fparamsi, *fpmu;
  fp1 = fopen(surf_nu, "r");
  for (i = 0; i < ND1; i++)
  {
    fscanf(fp1, "%lf %lf %lf", &x[i], &z[i], &tmp);
    x[i] = (x[i] - 145) * 0.739 / 100.0;
    z[i] = (z[i] - 145) * 0.739 / 100.0;
    printf("%d.\t %lf\t %lf \n", i, x[i], z[i]);
  }
  fclose(fp1);

  xgrids = 55;
  zgrids = 100;
  dz = 0.739 / zgrids;
  long double IOmega[xgrids + 1][zgrids + 1];
  long double mu[xgrids + 1][zgrids + 1];

  for (i = 0; i <= xgrids; i++)
  {
    for (j = 0; j <= zgrids; j++)
    {
      IOmega[i][j] = 0.0;
      mu[i][j] = 0.0;
    }
  }

  fp1 = fopen(out_fname, "w");
  if (fp1 == NULL)
  {
    printf("fp1 == NULL \n");
    exit(0);
  }
  for (countx = 0; countx <= xgrids; countx++)
  {
    if (countx <= (ND1 - 1))
    {
      xx = x[countx + 1];
      zSurf = z[countx + 1];
    }
    else
    {
      xx = x[ND1 - 1] + (0.739 / 100.0);
      zSurf = 0.0;
    }
    x_init[0] = ORE;
    x_init[1] = OIM;
    printf("\n****************\txx = %lf\tZSurf = %lf\t ORE = %e\t OIM = %e\t**********************\n\n", xx, zSurf, ORE, OIM);
    for (countz = 0; countz <= zgrids; countz++)
    {
      zz = 0.739 - (countz * dz);
      if (zz < zSurf)
      {
        fprintf(fp1, "%d\t%d\t%4.3e\t%4.3Le\t%4.3Le\n", countx, countz, KZ, xinit_tmp, xinit_tmp);
        printf("%d\t%d\t%4.3e\t%4.3Le\t%4.3Le\n", countx, countz, KZ, xinit_tmp, xinit_tmp);
      }
      else
      {

        const gsl_multiroot_fsolver_type *T;
        gsl_multiroot_fsolver *s;
        size_t iter = 0;

        lambda_xz = 0.0;
        mu_xz = 0.0;
        n_nu = 0.0;
        rho_xz = 0.0;
        ye_xz = 0.0;
        eps = 0.0;
        epsx = 0.0;
        epsy = 0.0;
        epsz = 0.0;
        phase1 = 0.0;
        phase2 = 0.0;

        C_Hnu = sqrt(2.0) * GFIHBARC * HBARC * HBARC * 1e-26;
        C_Hm = sqrt(2.0) * GFIHBARC * HBARC * HBARC / AMU * 1e-26;

        dw = 2.0;
        du = 2.0 / ((double)NU);
        dp = 2.0 * M_PI / ((double)NP);
        profile_generator(xx, zz);
        n_nu = read_nu_anu_profile();
        // ye_xz = read_ye(xx, zz);
        //rho_xz = read_rho(xx, zz);

        lambda_xz = rho_xz * ye_xz * C_Hm;
        mu_xz = n_nu * C_Hnu;
        p.lamr = 0.0;
        p.mur = 1.0;

        for (i = 0; i < NW; i++)
        {
          for (k = 0; k < NP; k++)
          {
            for (j = 0; j < NU; j++)
            {
              eps = eps + dw * du * dp / (4.0 * M_PI) * (gp[i][j][k] + gn[i][j][k]);                                            //[MeV]
              epsx = epsx + dw * du * dp / (4.0 * M_PI) * (gp[i][j][k] + gn[i][j][k]) * cos(pg[k]) * sqrt(1.0 - ug[j] * ug[j]); //[MeV]
              epsy = epsy + dw * du * dp / (4.0 * M_PI) * (gp[i][j][k] + gn[i][j][k]) * sin(pg[k]) * sqrt(1.0 - ug[j] * ug[j]); //[MeV]
              epsz = epsz + dw * du * dp / (4.0 * M_PI) * (gp[i][j][k] + gn[i][j][k]) * ug[j];                                  //[MeV]
              phase1 = phase1 + dw * du * dp / (4.0 * M_PI) * gp[i][j][k];                                                      //[MeV]
            }
          }
        }

        KZ = 0.0;
        gsl_multiroot_function f = {&rosenbrock_f, n, &p};
        gsl_vector *x = gsl_vector_alloc(n);
        T = gsl_multiroot_fsolver_hybrids;
        s = gsl_multiroot_fsolver_alloc(T, 2);
        gsl_vector_set(x, 0, x_init[0]);
        gsl_vector_set(x, 1, x_init[1]);
        gsl_multiroot_fsolver_set(s, &f, x);

        do
        {
          // finding the root from an initial guess
          iter++;
          status = gsl_multiroot_fsolver_iterate(s);
          if (status)
          {
            break;
          }
          status = gsl_multiroot_test_residual(s->f, 1e-12);
        } while (status == GSL_CONTINUE && iter < 200);
        print_state(iter, s);
        printf("status = %s\n", gsl_strerror(status));

        if (!status)
        {
          if (fabs(gsl_vector_get(s->x, 1)) > 8.0E-4)
          {
            x_init[0] = gsl_vector_get(s->x, 0);
            x_init[1] = fabs(gsl_vector_get(s->x, 1));
            IOmega[countx][countz] = fabs(gsl_vector_get(s->x, 1));
            mu[countx][countz] = mu_xz;
            fprintf(fp1, "%d\t%d\t%4.3e\t%4.3Le\t%4.3Le\n", countx, countz, KZ, (mu_xz * x_init[0] - lambda_xz), mu_xz * fabs(x_init[1]));
            printf("%d\t%d\t%4.3e\t%4.3Le\t%4.3Le\n", countx, countz, KZ, (mu_xz * x_init[0] - lambda_xz), mu_xz * fabs(x_init[1]));
          }
          else
          {
            fprintf(fp1, "%d\t%d\t%4.3e\t%4.3Le\t%4.3Le\n", countx, countz, KZ, (gsl_vector_get(s->x, 0) - lambda_xz), mu_xz * fabs(gsl_vector_get(s->x, 1)));
            printf("%d\t%d\t%4.3e\t%4.3Le\t%4.3Le\n", countx, countz, KZ, (mu_xz * gsl_vector_get(s->x, 0) - lambda_xz), mu_xz * fabs(gsl_vector_get(s->x, 1)));
            IOmega[countx][countz] = fabs(gsl_vector_get(s->x, 1));
            mu[countx][countz] = mu_xz;
          }
          if (count == 0)
          {
            ORE = x_init[0];
            OIM = x_init[1];
          }
        }
        else
        {
          long double solu = 0.0;
          fprintf(fp1, "%d\t%d\t%4.3e\t%4.3Le\t%4.3Le\n", countx, countz, KZ, solu, solu);
          printf("%d\t%d\t%4.3e\t%4.3Le\t%4.3Le\n", countx, countz, KZ, solu, solu);
          IOmega[countx][countz] = solu;
          mu[countx][countz] = mu_xz;
        }

        gsl_multiroot_fsolver_free(s);
        gsl_vector_free(x);
      }
      printf("\n");
      count++;
    }
    fprintf(fp1, "\n");
    count = 0;
  }
  fclose(fp1);
  for (countx = 0; countx <= xgrids; countx++)
  {
    for (countz = 1; countz <= (zgrids - 1); countz++)
    {
      if (IOmega[countx][countz] == 0.0)
      {
        if ((IOmega[countx][countz - 1] > 0.0) && (IOmega[countx][countz + 1]) > 0.0)
        {
          IOmega[countx][countz] = (IOmega[countx][countz - 1] + IOmega[countx][countz + 1]) / 2.0;
        }
      }
    }
  }
  fp2 = fopen(interpole_fname, "w");
  fpmu = fopen(mufname, "w");
  for (countx = 0; countx <= xgrids; countx++)
  {
    for (countz = 0; countz <= zgrids; countz++)
    {
      fprintf(fp2, "%d\t%d\t%LE\n", countx, countz, IOmega[countx][countz]);
      fprintf(fp2, "%d\t%d\t%LE\n", countx, countz, mu[countx][countz]);
    }
    fprintf(fp2, "\n");
  }
  fclose(fp2);
  fclose(fpmu);
  return (0);
}
