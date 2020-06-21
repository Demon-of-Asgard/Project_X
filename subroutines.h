// User defined function
long double read_nu_anu_profile()
{
  int i, j, k, status, iter = 0;
  int count = 0, file_length = 0;
  long double tmp_data;
  long double tmp_ctheta, tmp_phi, tmp_nnu, tmp_nanu;
  long double n_nu = 0.0, n_anu = 0.0;
  char tmp_char, last_char, *file_read_status;
  FILE *fpr, *fpw;
  //
  //
  fpr = fopen(profname, "r");
  if (fpr == NULL)
  {
    printf("Unable to open profile. Exiting.\n");
    fclose(fpw);
    exit(0);
  }
  else
  {
    while (fscanf(fpr, "%LE", &tmp_data) != EOF)
    {
      // printf("%LE\n", tmp_data);
      count++;
    }
  }
  fclose(fpr);
  file_length = count;
  //
  count = 0;
  fpr = fopen(profname, "r");
  if (fpr == NULL)
  {
    printf("Unable to open profile. Exiting.\n");
    fclose(fpw);
    exit(0);
  }
  //
  while ((tmp_char = fgetc(fpr)) != EOF)
  {

    //printf("%c\n", tmp_char);

    if (tmp_char == ' ')
      count++;
    if (tmp_char == '\n')
    {
      count++;
      break;
    }
    last_char = tmp_char;
  }
  fclose(fpr);
  //
  int columns = count;
  int rows = file_length / columns;

  fpr = fopen(profname, "r");
  if (fpr == NULL)
  {
    printf("Unable to open profile. Exiting.\n");
    fclose(fpw);
    exit(0);
  }
  iter = 0;
  for (i = 0; i < NW; i++)
  {
    for (k = 0; k < NP; k++)
    {
      for (j = 0; j < NU; j++)
      {
        iter++;
        fscanf(fpr, "%Le %Le %Le %Le", &tmp_ctheta, &tmp_phi, &tmp_nnu, &tmp_nanu);
        if (k == 0)
        {
          ug[j] = tmp_ctheta;
        }
        gp[i][j][k] = tmp_nnu;
        gn[i][j][k] = tmp_nanu;
      }
      pg[k] = tmp_phi;
    }
  }
  fclose(fpr);
  //
  if (iter != NU * NP * NW)
  {
    printf("ERROR!!!... File length does not match with number of iterations. Quitting the program. :( \n\n");
    exit(0);
  }

  for (i = 0; i < NW; i++) //Calculating neutrino number densities.
  {
    for (k = 0; k < NP; k++)
    {
      for (j = 0; j < NU; j++)
      {
        n_nu += gp[i][j][k] * du * dp;  // / (2.0 * M_PI);
        n_anu += gn[i][j][k] * du * dp; // / (2.0 * M_PI);
      }
    }
  }

  for (i = 0; i < NW; i++) //Normalizing profiles.
  {
    for (k = 0; k < NP; k++)
    {
      for (j = 0; j < NU; j++)
      {
        gp[i][j][k] = gp[i][j][k] / n_nu;
        gn[i][j][k] = -1 * gn[i][j][k] / n_nu;
      }
    }
  }

  double int_gp = 0.0;
  double int_gn = 0.0;
  for (i = 0; i < NW; i++) //Testing normalization
  {
    for (k = 0; k < NP; k++)
    {
      for (j = 0; j < NU; j++)
      {
        int_gp += gp[i][j][k] * du * dp;
        int_gn += gn[i][j][k] * du * dp;
      }
    }
  }
  printf("int_gp = %lf   int_gn = %lf, alpha = %lf \n", int_gp, -int_gn, -int_gn / int_gp);
  //
  alpha = -int_gn / int_gp;
  return (n_nu);
}
//
//                                                      BI_LINEAR SUBROUTINE
//
long double bi_linear(double x1, double x2, double z1, double z2, long double f11, long double f12, long double f21, long double f22)
{
  long double f;
  f = (f11 * (x2 - xx) * (z2 - zz) + f21 * (xx - x1) * (z2 - zz) + f12 * (x2 - xx) * (zz - z1) + f22 * (xx - x1) * (zz - z1)) / ((x2 - x1) * (z2 - z1));
  return (f);
}
//
//
long double read_rho(double xx, double zz) //                                             READ  RHO
{
  FILE *fp1;
  int i, j;
  double tmp;
  double xg[NR], zg[NT];
  long double rho[NR][NT];
  char *fname;

  long double rhoxz = 0.0;

  // rewrite this part to read the rho and ye on your cartesian grids

  // fp1 = fopen("/home/manu/MEGAsync/Work_IOP/NSNS_MERGER_SNAPSHOT/INPUT_DATA_FILES/dd21357p5/dd21357p5rho.dat", "r");
  fp1 = fopen("sfho1355rho.dat", "r");

  for (i = 0; i < NR; i++)
  {
    for (j = 0; j < NT; j++)
    {
      fscanf(fp1, "%lf %lf %LE\n", &xg[i], &zg[j], &rho[i][j]); //x[i]----->rg[i],z[j]----->tg[j],rho--->g/cm^-3  (edited---MG)
      xg[i] = (xg[i] - 145.0) * 0.793 / 100.0;                  //km
      zg[j] = (zg[j] - 145.0) * 0.793 / 100.0;                  //km
      if (i > 0 && j > 0)
      {
        if (xg[i - 1] < xx && xg[i] >= xx)
        {
          if (zg[j - 1] < zz && zg[j] >= zz)
          {
            rhoxz = bi_linear(xg[i - 1], xg[i], zg[j - 1], zg[j], rho[i - 1][j - 1], rho[i - 1][j], rho[i][j - 1], rho[i][j]); //[gm per cm^-3]
          }
        }
      }
    }
  }
  fclose(fp1);
  return (rhoxz);
}
//
long double read_ye(double xx, double zz) //                                        READ YE
{
  FILE *fp1;
  int i, j;
  double tmp;
  double xg[NR], zg[NT];
  long double ye[NR][NT];
  char *fname;

  long double yexz = 0.0;

  // rewrite this part to read the rho and ye on your cartesian grids

  //fp1 = fopen("/home/manu/MEGAsync/Work_IOP/NSNS_MERGER_SNAPSHOT/INPUT_DATA_FILES/dd21357p5/dd21357p5Ye.dat", "r");
  fp1 = fopen("sfho1355Ye.dat", "r");

  for (i = 0; i < NR; i++)
  {
    for (j = 0; j < NT; j++)
    {
      fscanf(fp1, "%lf %lf %LE\n", &xg[i], &zg[j], &ye[i][j]); //x[i]----->rg[i],z[j]----->tg[j],rho--->g/cm^-3  (edited---MG)
      xg[i] = (xg[i] - 145.0) * 0.793 / 100.0;                 //km
      zg[j] = (zg[j] - 145.0) * 0.793 / 100.0;                 //km
      if (i > 0 && j > 0)
      {
        if (xg[i - 1] < xx && xg[i] >= xx)
        {
          if (zg[j - 1] < zz && zg[j] >= zz)
          {
            yexz = bi_linear(xg[i - 1], xg[i], zg[j - 1], zg[j], ye[i - 1][j - 1], ye[i - 1][j], ye[i][j - 1], ye[i][j]); //DIMENSION LESS
          }
        }
      }
    }
  }
  fclose(fp1);
  return (yexz);
}
