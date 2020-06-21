int rosenbrock_f(const gsl_vector *x, void *params, gsl_vector *f)
{
  double mur = ((struct rparams *)params)->mur;
  double lamr = ((struct rparams *)params)->lamr;
  //  double eps = ((struct rparams *) params)->eps;
  const double x0 = gsl_vector_get(x, 0);
  const double x1 = gsl_vector_get(x, 1);

  double y0 = 0.0;
  double y1 = 0.0;
  int i, j, k;
  double lambar, vx, vy, vz, v0;
  double complex omega;
  double complex R00, R01, R02, R03, R11, R12, R13, R22, R23, R33, cfactor, yc;

  R00 = 0.0 + 0.0 * I;
  R01 = 0.0 + 0.0 * I;
  R02 = 0.0 + 0.0 * I;
  R03 = 0.0 + 0.0 * I;
  R11 = 0.0 + 0.0 * I;
  R12 = 0.0 + 0.0 * I;
  R13 = 0.0 + 0.0 * I;
  R22 = 0.0 + 0.0 * I;
  R23 = 0.0 + 0.0 * I;
  R33 = 0.0 + 0.0 * I;
  omega = x0 + I * x1;
  lambar = lamr + eps * mur;

  for (i = 0; i < NW; i++)
  {
    v0 = wg[i];
    for (k = 0; k < NP; k++)
    {
      //vz = ug[j];
      for (j = 0; j < NU; j++)
      {
        vz = ug[j];
        vx = cos(pg[k]) * sqrt(1.0 - vz * vz);
        vy = sin(pg[k]) * sqrt(1.0 - vz * vz);
        //        cfactor=v0-W0+lambar-vx*mur*epsx-vy*mur*epsy+vz*(omega-mur*epsz);
        cfactor = v0 + lambar - omega - vx * mur * epsx - vy * mur * epsy - vz * (-KZ + mur * epsz); // check this is correct
        R00 = R00 + dw * du * dp / (4.0 * M_PI) * gp[i][j][k] / cfactor;
        R01 = R01 + dw * du * dp / (4.0 * M_PI) * gp[i][j][k] * vx / cfactor;
        R02 = R02 + dw * du * dp / (4.0 * M_PI) * gp[i][j][k] * vy / cfactor;
        R03 = R03 + dw * du * dp / (4.0 * M_PI) * gp[i][j][k] * vz / cfactor;
        R11 = R11 + dw * du * dp / (4.0 * M_PI) * gp[i][j][k] * vx * vx / cfactor;
        R12 = R12 + dw * du * dp / (4.0 * M_PI) * gp[i][j][k] * vx * vy / cfactor;
        R13 = R13 + dw * du * dp / (4.0 * M_PI) * gp[i][j][k] * vx * vz / cfactor;
        R22 = R22 + dw * du * dp / (4.0 * M_PI) * gp[i][j][k] * vy * vy / cfactor;
        R23 = R23 + dw * du * dp / (4.0 * M_PI) * gp[i][j][k] * vy * vz / cfactor;
        R33 = R33 + dw * du * dp / (4.0 * M_PI) * gp[i][j][k] * vz * vz / cfactor;
        cfactor = cfactor - 2.0 * v0;
        R00 = R00 + dw * du * dp / (4.0 * M_PI) * gn[i][j][k] / cfactor;
        R01 = R01 + dw * du * dp / (4.0 * M_PI) * gn[i][j][k] * vx / cfactor;
        R02 = R02 + dw * du * dp / (4.0 * M_PI) * gn[i][j][k] * vy / cfactor;
        R03 = R03 + dw * du * dp / (4.0 * M_PI) * gn[i][j][k] * vz / cfactor;
        R11 = R11 + dw * du * dp / (4.0 * M_PI) * gn[i][j][k] * vx * vx / cfactor;
        R12 = R12 + dw * du * dp / (4.0 * M_PI) * gn[i][j][k] * vx * vy / cfactor;
        R13 = R13 + dw * du * dp / (4.0 * M_PI) * gn[i][j][k] * vx * vz / cfactor;
        R22 = R22 + dw * du * dp / (4.0 * M_PI) * gn[i][j][k] * vy * vy / cfactor;
        R23 = R23 + dw * du * dp / (4.0 * M_PI) * gn[i][j][k] * vy * vz / cfactor;
        R33 = R33 + dw * du * dp / (4.0 * M_PI) * gn[i][j][k] * vz * vz / cfactor;
      }
    }
  }
  R00 = 1.0 - mur * R00;
  R11 = -1.0 - mur * R11;
  R22 = -1.0 - mur * R22;
  R33 = -1.0 - mur * R33;
  R01 = mur * R01;
  R02 = mur * R02;
  R03 = mur * R03;
  R12 = -mur * R12;
  R13 = -mur * R13;
  R23 = -mur * R23;

  //printf("%8.7e %8.7e\n",mur,lambar);
  //printf("%8.7e %8.7e\n",creal(R22),cimag(R22));
  //getchar();

  yc = R22; // axial-symmtry breaking solution
            //  yc=-R03*R03+R00*R33;   // axial-symmtry conserving solution

  // below is the full solution
  //  yc=R02*R02*R13*R13-R00*R13*R13*R22+R03*R03*(R12*R12-R11*R22)
  //     -2.0*R01*R02*R13*R23+2.0*R00*R12*R13*R23+ R01*R01*R23*R23
  //     -R00*R11*R23*R23+2.0*R03*(-R02*R12*R13+R01*R13*R22+R02*R11*R23-R01*R12*R23)
  //     -(R02*R02*R11-2.0*R01*R02*R12+R00*R12*R12+R01*R01*R22-R00*R11*R22)*R33;

  y0 = creal(yc);
  y1 = cimag(yc);

  gsl_vector_set(f, 0, y0);
  gsl_vector_set(f, 1, y1);
  return GSL_SUCCESS;
}

int print_state(size_t iter, gsl_multiroot_fsolver *s)
{
  printf("\n iter = %3lu x = %5.4e %5.4e f(x) = % .3e % .3e\n", iter, gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1), gsl_vector_get(s->f, 0), gsl_vector_get(s->f, 1));
  return (0);
}
