//CREATE PROFILE

int profile_generator(float x_point, float z_point)
{
  int i, j, k, id1m[8], id2m[8], ii, jj, nar;
  double p, u, du, dp, nnu1, nnu2, uf1, uf2, fnu1 = 0.0, fnu2 = 0.0, dx, dz, Z0[NX], X0[NX];
  double rd1[ND1], zd1[ND1], ud1[ND1][2], ufd1[ND1][2];
  double nnud1[ND1];
  double cad1[ND1], sad1[ND1], cbd1[ND1][2], sbd1[ND1][2], sgnp1[ND1];
  double rd2[ND2], zd2[ND2], ud2[ND2][2], ufd2[ND2][2];
  double nnud2[ND2];
  double cad2[ND2], sad2[ND2], cbd2[ND2][2], sbd2[ND2][2], sgnp2[ND2];
  double rd1o[ND1], zd1o[ND1], nnud1o[ND1];
  double cad1o[ND1], sad1o[ND1], sgnp1o[ND1];
  double rd2o[ND2], zd2o[ND2], nnud2o[ND2];
  double cad2o[ND2], sad2o[ND2], sgnp2o[ND2];
  double dsint, fac, tmp, att, relang, ud1m[8], ud2m[8];
  double *xar, *yar1, *yar2, *yar12, *yar22;
  FILE *fp1, *fp2, *fp3, *fp4;

  xx = x_point;
  zz = z_point;
  fp1 = fopen(surf_nu, "r");
  fp4 = fopen(surf_anu, "r");

  fp2 = fopen(profname, "w");
  fp3 = fopen("deb.dat", "w");
  for (i = 0; i < ND1; i++)
  { // read in the grid data of neutrino surface and compute the normal direction associated with the grids, nue
    fscanf(fp1, "%lf %lf %lf", &rd1o[i], &zd1o[i], &nnud1o[i]);
    rd1o[i] = (rd1o[i] - 145.0 + 0.001) * 0.739 / 100.0;
    zd1o[i] = (zd1o[i] - 145.0) * 0.739 / 100.0;
  }
  tmp = sqrt((zd1o[1] - zd1o[0]) * (zd1o[1] - zd1o[0]) + (rd1o[1] - rd1o[0]) * (rd1o[1] - rd1o[0])); // first point
  sad1o[0] = fabs(zd1o[1] - zd1o[0]) / tmp;
  cad1o[0] = fabs(rd1o[1] - rd1o[0]) / tmp;
  sgnp1o[0] = (zd1o[1] - zd1o[0]) / fabs(zd1o[1] - zd1o[0]);
  tmp = sqrt((zd1o[ND1 - 1] - zd1o[ND1 - 2]) * (zd1o[ND1 - 1] - zd1o[ND1 - 2]) + (rd1o[ND1 - 1] - rd1o[ND1 - 2]) * (rd1o[ND1 - 1] - rd1o[ND1 - 2])); // last point
  sad1o[ND1 - 1] = fabs(zd1o[ND1 - 1] - zd1o[ND1 - 2]) / tmp;
  cad1o[ND1 - 1] = fabs(rd1o[ND1 - 1] - rd1o[ND1 - 2]) / tmp;
  sgnp1o[ND1 - 1] = (zd1o[ND1 - 1] - zd1o[ND1 - 2]) / fabs(zd1o[ND1 - 1] - zd1o[ND1 - 2]);
  for (i = 1; i < ND1 - 1; i++)
  { // all the intermediate points
    tmp = sqrt((zd1o[i + 1] - zd1o[i - 1]) * (zd1o[i + 1] - zd1o[i - 1]) + (rd1o[i + 1] - rd1o[i - 1]) * (rd1o[i + 1] - rd1o[i - 1]));
    sad1o[i] = fabs(zd1o[i + 1] - zd1o[i - 1]) / tmp;
    cad1o[i] = fabs(rd1o[i + 1] - rd1o[i - 1]) / tmp;
    sgnp1o[i] = (zd1o[i + 1] - zd1o[i - 1]) / fabs(zd1o[i + 1] - zd1o[i - 1]);
  }
  for (i = 0; i < ND2; i++)
  { // same thing for nuebar
    fscanf(fp4, "%lf %lf %lf", &rd2o[i], &zd2o[i], &nnud2o[i]);
    rd2o[i] = (rd2o[i] - 145.0 + 0.001) * 0.739 / 100.0;
    zd2o[i] = (zd2o[i] - 145.0) * 0.739 / 100.0;
  }
  tmp = sqrt((zd2o[1] - zd2o[0]) * (zd2o[1] - zd2o[0]) + (rd2o[1] - rd2o[0]) * (rd2o[1] - rd2o[0]));
  sad2o[0] = fabs(zd2o[1] - zd2o[0]) / tmp;
  cad2o[0] = fabs(rd2o[1] - rd2o[0]) / tmp;
  sgnp2o[0] = (zd2o[1] - zd2o[0]) / fabs(zd2o[1] - zd2o[0]);
  tmp = sqrt((zd2o[ND2 - 1] - zd2o[ND2 - 2]) * (zd2o[ND2 - 1] - zd2o[ND2 - 2]) + (rd2o[ND2 - 1] - rd2o[ND2 - 2]) * (rd2o[ND2 - 1] - rd2o[ND2 - 2]));
  sad2o[ND2 - 1] = fabs(zd2o[ND2 - 1] - zd2o[ND2 - 2]) / tmp;
  cad2o[ND2 - 1] = fabs(rd2o[ND2 - 1] - rd2o[ND2 - 2]) / tmp;
  sgnp2o[ND2 - 1] = (zd2o[ND2 - 1] - zd2o[ND2 - 2]) / fabs(zd2o[ND2 - 1] - zd2o[ND2 - 2]);
  for (i = 1; i < ND2 - 1; i++)
  {
    tmp = sqrt((zd2o[i + 1] - zd2o[i - 1]) * (zd2o[i + 1] - zd2o[i - 1]) + (rd2o[i + 1] - rd2o[i - 1]) * (rd2o[i + 1] - rd2o[i - 1]));
    sad2o[i] = fabs(zd2o[i + 1] - zd2o[i - 1]) / tmp;
    cad2o[i] = fabs(rd2o[i + 1] - rd2o[i - 1]) / tmp;
    sgnp2o[i] = (zd2o[i + 1] - zd2o[i - 1]) / fabs(zd2o[i + 1] - zd2o[i - 1]);
  }
  du = 2.0 / ((double)NU);        // local angular grid size in cos\theta
  dp = 2.0 * M_PI / ((double)NP); // local angular grid size in phi
  fnu1 = 0.0;
  fnu2 = 0.0;
  for (k = 0; k < NP; k++)
  {                                         // go through the phi direction loop first
                                            //for(k=174;k<175;k++){
    p = -M_PI + dp * (2.0 * k + 1.0) / 2.0; // phi value
                                            //p=0.0;
    for (i = 0; i < ND1; i++)
    {
      rd1[i] = rd1o[i];
      zd1[i] = zd1o[i];
      nnud1[i] = nnud1o[i];
      cad1[i] = cad1o[i];
      sad1[i] = sad1o[i];
      sgnp1[i] = sgnp1o[i];
    }
    if (xx > rd1o[0] && xx < rd1o[ND1 - 1])
    { // add one extra data point with same x-coordinate as xx
      for (i = 0; i < ND1 - 1; i++)
      {
        if (xx > rd1o[i] && xx <= rd1o[i + 1])
        {
          rd1[i + 1] = xx;
          tmp = (xx - rd1o[i]) / (rd1o[i + 1] - rd1o[i]);
          zd1[i + 1] = zd1o[i] + (zd1o[i + 1] - zd1o[i]) * tmp;
          nnud1[i + 1] = nnud1o[i] + (nnud1o[i + 1] - nnud1o[i]) * tmp;
          cad1[i + 1] = cad1o[i] + (cad1o[i + 1] - cad1o[i]) * tmp;
          sad1[i + 1] = sqrt(1.0 - cad1[i + 1] * cad1[i + 1]);
          break;
        }
      }
    }
    if (fabs(xx * sin(p)) > rd1o[0] && fabs(xx * sin(p)) < rd1o[ND1 - 1])
    { // MRW: ?
      for (i = 0; i < ND1 - 1; i++)
      {
        if (fabs(xx * sin(p)) > rd1o[i] && fabs(xx * sin(p)) <= rd1o[i + 1])
        {
          rd1[i] = fabs(xx * sin(p));
          tmp = (fabs(xx * sin(p)) - rd1o[i]) / (rd1o[i + 1] - rd1o[i]);
          zd1[i] = zd1o[i] + (zd1o[i + 1] - zd1o[i]) * tmp;
          nnud1[i] = nnud1o[i] + (nnud1o[i + 1] - nnud1o[i]) * tmp;
          cad1[i] = cad1o[i] + (cad1o[i + 1] - cad1o[i]) * tmp;
          sad1[i] = sqrt(1.0 - cad1[i + 1] * cad1[i + 1]);
          break;
        }
      }
    }
    for (i = 0; i < ND1; i++)
    {
      fac = xx * xx * cos(p) * cos(p) - xx * xx + rd1[i] * rd1[i];
      if (fabs(fac) < 1.0e-12)
      {
        fac = 0.0;
      }
      if (fac >= 0.0)
      {                                  // take only the points that can reach (xx,zz) with a given phi
        dsint = xx * cos(p) + sqrt(fac); // first solution for a given (xd, zd)
        if (fabs(dsint) < 1.0e-12)
        {
          dsint = 0.0;
        }
        if (dsint >= 0.0)
        { // if first solution exists
          att = atan(dsint / (zz - zd1[i]));
          if (att < 0)
          {
            att = att + M_PI;
          }
          tmp = sqrt((xx - dsint * cos(p)) * (xx - dsint * cos(p)) + dsint * sin(p) * dsint * sin(p));
          cbd1[i][0] = -(xx - dsint * cos(p)) / tmp * sgnp1[i];
          sbd1[i][0] = dsint * sin(p) / tmp * sgnp1[i];
          relang = sad1[i] * cbd1[i][0] * sin(att) * cos(p) + sad1[i] * sbd1[i][0] * sin(att) * sin(p) + cad1[i] * cos(att);
          if (relang >= 0.0)
          {                       // take the neutrinos that do not cross the disk itself
            ud1[i][0] = cos(att); // cos\theta for that (xd, zd)
            ufd1[i][0] = relang;  // cos(angle) relative to the normal of (xd, zd)
          }
          else
          {
            ud1[i][0] = NAN;
            ufd1[i][0] = NAN;
          }
        }
        else
        {
          ud1[i][0] = NAN;
          ufd1[i][0] = NAN;
        }
        dsint = xx * cos(p) - sqrt(fac); // second solution for a given (xd, zd)
        if (fabs(dsint) < 1.0e-12)
        {
          dsint = 0.0;
        }
        if (dsint >= 0.0)
        { // if second solution exists
          att = atan(dsint / (zz - zd1[i]));
          if (att < 0)
          {
            att = att + M_PI;
          }
          tmp = sqrt((xx - dsint * cos(p)) * (xx - dsint * cos(p)) + dsint * sin(p) * dsint * sin(p));
          cbd1[i][1] = -(xx - dsint * cos(p)) / tmp * sgnp1[i];
          sbd1[i][1] = dsint * sin(p) / tmp * sgnp1[i];
          relang = sad1[i] * cbd1[i][1] * sin(att) * cos(p) + sad1[i] * sbd1[i][1] * sin(att) * sin(p) + cad1[i] * cos(att);
          if (relang >= 0.0)
          {
            ud1[i][1] = cos(att);
            ufd1[i][1] = relang;
          }
          else
          {
            ud1[i][1] = NAN;
            ufd1[i][1] = NAN;
          }
        }
        else
        {
          ud1[i][1] = NAN;
          ufd1[i][1] = NAN;
        }
      }
      else
      {
        ud1[i][0] = NAN;
        ud1[i][1] = NAN;
        ufd1[i][0] = NAN;
        ufd1[i][1] = NAN;
      }
    }
    ud1m[0] = 0.0;
    ud1m[1] = 0.0;
    ud1m[2] = 0.0;
    ud1m[3] = 0.0;
    ud1m[4] = 0.0;
    ud1m[5] = 0.0;
    ud1m[6] = 0.0;
    ud1m[7] = 0.0;
    id1m[0] = 0;
    id1m[1] = 0;
    id1m[2] = 0;
    id1m[3] = 0;
    id1m[4] = 0;
    id1m[5] = 0;
    id1m[6] = 0;
    id1m[7] = 0;
    if (isnan(ud1[1][0]) == 1)
    {
      ud1[0][0] = NAN;
      ufd1[0][0] = NAN;
    }
    if (isnan(ud1[ND1 - 2][0]) == 1)
    {
      ud1[ND1 - 1][0] = NAN;
      ufd1[ND1 - 1][0] = NAN;
    }
    for (i = 1; i < ND1 - 1; i++)
    {
      if (isnan(ud1[i + 1][0]) == 1 && isnan(ud1[i - 1][0]) == 1)
      {
        ud1[i][0] = NAN;
        ufd1[i][0] = NAN;
      }
      if (isnan(ud1[i + 1][1]) == 1 && isnan(ud1[i - 1][1]) == 1)
      {
        ud1[i][1] = NAN;
        ufd1[i][1] = NAN;
      }
    }
    for (i = 0; i < ND1; i++)
    {
      fprintf(fp3, "%lf %lf %lf %lf %lf %lf %d\n", rd1[i], zd1[i], ud1[i][0], ufd1[i][0], ud1[i][1], ufd1[i][1], i);
    }
    for (i = ND1 - 1; i >= 0; i--)
    {
      if (ud1[i][1] < 50.0)
      {
        ud1m[7] = ud1[i][1];
        id1m[7] = i;
        break;
      }
    }
    for (i = id1m[7]; i > 0; i--)
    {
      if (isnan(ud1[i - 1][1]) == 1)
      {
        ud1m[6] = ud1[i][1];
        id1m[6] = i;
        break;
      }
    }
    if (i == 0)
    {
      ud1m[6] = ud1[0][1];
      id1m[6] = 0;
    }
    for (i = 0; i < ND1; i++)
    {
      if (ud1[i][1] < 50.0)
      {
        ud1m[4] = ud1[i][1];
        id1m[4] = i;
        break;
      }
    }
    for (i = id1m[4]; i < ND1 - 1; i++)
    {
      if (isnan(ud1[i + 1][1]) == 1)
      {
        ud1m[5] = ud1[i][1];
        id1m[5] = i;
        break;
      }
    }
    if (i == ND1 - 1)
    {
      ud1m[5] = ud1[ND1 - 1][1];
      id1m[5] = ND1 - 1;
    } // MRW: ?
    if (((ud1m[4] - ud1m[6]) * (ud1m[4] - ud1m[7]) < -1.0e-12) && ((ud1m[5] - ud1m[6]) * (ud1m[5] - ud1m[7]) < -1.0e-12))
    {
      id1m[4] = 0;
      id1m[5] = 0;
      ud1m[4] = 0.0;
      ud1m[5] = 0.0;
    }
    if (((ud1m[6] - ud1m[4]) * (ud1m[6] - ud1m[5]) < -1.0e-12) && ((ud1m[7] - ud1m[4]) * (ud1m[7] - ud1m[5]) < -1.0e-12))
    {
      id1m[6] = 0;
      id1m[7] = 0;
      ud1m[6] = 0.0;
      ud1m[7] = 0.0;
    }
    /*        for(i=0;i<ND1;i++){
          if(ud1[i][1]<50.0){ud1m[2]=ud1[i][1];id1m[2]=i;break;}
        }
        for(i=ND1-1;i>=0;i--){
          if(ud1[i][1]<50.0){ud1m[3]=ud1[i][1];id1m[3]=i;break;}
        }*/
    for (i = 0; i < ND1; i++)
    {
      if ((ud1[i][0] - ud1m[4]) * (ud1[i][0] - ud1m[5]) < -1.0e-12)
      {
        ud1[i][0] = NAN;
        ufd1[i][0] = NAN;
      }
      if ((ud1[i][0] - ud1m[6]) * (ud1[i][0] - ud1m[7]) < -1.0e-12)
      {
        ud1[i][0] = NAN;
        ufd1[i][0] = NAN;
      }
    }
    /*        for(i=0;i<ND1;i++){
          if((ud1[i][0]-ud1m[2])*(ud1[i][0]-ud1m[3])<-1.0e-12){ud1[i][0]=NAN;ufd1[i][0]=NAN;}
//fprintf(fp2,"%lf %lf %lf %lf %lf %lf\n",rd1[i],zd1[i],ud1[i][0],ufd1[i][0],ud1[i][1],ufd1[i][1]);
        }*/
    for (i = 0; i < ND1; i++)
    {
      if (ud1[i][0] < 50.0)
      {
        ud1m[0] = ud1[i][0];
        id1m[0] = i;
        break;
      }
    }
    for (i = ND1 - 1; i >= 0; i--)
    {
      if (ud1[i][0] < 50.0)
      {
        ud1m[1] = ud1[i][0];
        id1m[1] = i;
        break;
      }
    }
    //printf("%lf %lf %lf %lf %lf %lf %lf %lf\n",ud1m[0],ud1m[1],ud1m[2],ud1m[3],ud1m[4],ud1m[5],ud1m[6],ud1m[7]);
    for (i = 0; i < ND2; i++)
    {
      rd2[i] = rd2o[i];
      zd2[i] = zd2o[i];
      nnud2[i] = nnud2o[i];
      cad2[i] = cad2o[i];
      sad2[i] = sad2o[i];
      sgnp2[i] = sgnp2o[i];
    }
    if (xx > rd2o[0] && xx < rd2o[ND2 - 1])
    {
      for (i = 0; i < ND2 - 1; i++)
      {
        if (xx > rd2o[i] && xx <= rd2o[i + 1])
        {
          rd2[i + 1] = xx;
          tmp = (xx - rd2o[i]) / (rd2o[i + 1] - rd2o[i]);
          zd2[i + 1] = zd2o[i] + (zd2o[i + 1] - zd2o[i]) * tmp;
          nnud2[i + 1] = nnud2o[i] + (nnud2o[i + 1] - nnud2o[i]) * tmp;
          cad2[i + 1] = cad2o[i] + (cad2o[i + 1] - cad2o[i]) * tmp;
          sad2[i + 1] = sqrt(1.0 - cad2[i + 1] * cad2[i + 1]);
          break;
        }
      }
    }
    if (fabs(xx * sin(p)) > rd2o[0] && fabs(xx * sin(p)) < rd2o[ND2 - 1])
    {
      for (i = 0; i < ND2 - 1; i++)
      {
        if (fabs(xx * sin(p)) > rd2o[i] && fabs(xx * sin(p)) <= rd2o[i + 1])
        {
          rd2[i] = fabs(xx * sin(p));
          tmp = (fabs(xx * sin(p)) - rd2o[i]) / (rd2o[i + 1] - rd2o[i]);
          zd2[i] = zd2o[i] + (zd2o[i + 1] - zd2o[i]) * tmp;
          nnud2[i] = nnud2o[i] + (nnud2o[i + 1] - nnud2o[i]) * tmp;
          cad2[i] = cad2o[i] + (cad2o[i + 1] - cad2o[i]) * tmp;
          sad2[i] = sqrt(1.0 - cad2[i + 1] * cad2[i + 1]);
          break;
        }
      }
    }
    for (i = 0; i < ND2; i++)
    {
      fac = xx * xx * cos(p) * cos(p) - xx * xx + rd2[i] * rd2[i];
      if (fabs(fac) < 1.0e-12)
      {
        fac = 0.0;
      }
      if (fac >= 0.0)
      {
        dsint = xx * cos(p) + sqrt(fac);
        if (fabs(dsint) < 1.0e-12)
        {
          dsint = 0.0;
        }
        if (dsint >= 0.0)
        {
          att = atan(dsint / (zz - zd2[i]));
          if (att < 0)
          {
            att = att + M_PI;
          }
          tmp = sqrt((xx - dsint * cos(p)) * (xx - dsint * cos(p)) + dsint * sin(p) * dsint * sin(p));
          cbd2[i][0] = -(xx - dsint * cos(p)) / tmp * sgnp2[i];
          sbd2[i][0] = dsint * sin(p) / tmp * sgnp2[i];
          relang = sad2[i] * cbd2[i][0] * sin(att) * cos(p) + sad2[i] * sbd2[i][0] * sin(att) * sin(p) + cad2[i] * cos(att);
          if (relang >= 0.0)
          {
            ud2[i][0] = cos(att);
            ufd2[i][0] = relang;
          }
          else
          {
            ud2[i][0] = NAN;
            ufd2[i][0] = NAN;
          }
        }
        else
        {
          ud2[i][0] = NAN;
          ufd2[i][0] = NAN;
        }
        dsint = xx * cos(p) - sqrt(fac);
        if (fabs(dsint) < 1.0e-12)
        {
          dsint = 0.0;
        }
        if (dsint >= 0.0)
        {
          att = atan(dsint / (zz - zd2[i]));
          if (att < 0)
          {
            att = att + M_PI;
          }
          tmp = sqrt((xx - dsint * cos(p)) * (xx - dsint * cos(p)) + dsint * sin(p) * dsint * sin(p));
          cbd2[i][1] = -(xx - dsint * cos(p)) / tmp * sgnp2[i];
          sbd2[i][1] = dsint * sin(p) / tmp * sgnp2[i];
          relang = sad2[i] * cbd2[i][1] * sin(att) * cos(p) + sad2[i] * sbd2[i][1] * sin(att) * sin(p) + cad2[i] * cos(att);
          if (relang >= 0.0)
          {
            ud2[i][1] = cos(att);
            ufd2[i][1] = relang;
          }
          else
          {
            ud2[i][1] = NAN;
            ufd2[i][1] = NAN;
          }
        }
        else
        {
          ud2[i][1] = NAN;
          ufd2[i][1] = NAN;
        }
      }
      else
      {
        ud2[i][0] = NAN;
        ud2[i][1] = NAN;
        ufd2[i][0] = NAN;
        ufd2[i][1] = NAN;
      }
    }
    ud2m[0] = 0.0;
    ud2m[1] = 0.0;
    ud2m[2] = 0.0;
    ud2m[3] = 0.0;
    ud2m[4] = 0.0;
    ud2m[5] = 0.0;
    ud2m[6] = 0.0;
    ud2m[7] = 0.0;
    id2m[0] = 0;
    id2m[1] = 0;
    id2m[2] = 0;
    id2m[3] = 0;
    id2m[4] = 0;
    id2m[5] = 0;
    id2m[6] = 0;
    id2m[7] = 0;
    if (isnan(ud2[1][0]) == 1)
    {
      ud2[0][0] = NAN;
      ufd2[0][0] = NAN;
    }
    if (isnan(ud2[ND2 - 2][0]) == 1)
    {
      ud2[ND2 - 1][0] = NAN;
      ufd2[ND2 - 1][0] = NAN;
    }
    for (i = 1; i < ND2 - 1; i++)
    {
      if (isnan(ud2[i + 1][0]) == 1 && isnan(ud2[i - 1][0]) == 1)
      {
        ud2[i][0] = NAN;
        ufd2[i][0] = NAN;
      }
      if (isnan(ud2[i + 1][1]) == 1 && isnan(ud2[i - 1][1]) == 1)
      {
        ud2[i][1] = NAN;
        ufd2[i][1] = NAN;
      }
    }
    for (i = 0; i < ND2; i++)
    {
      fprintf(fp3, "%lf %lf %lf %lf %lf %lf %d\n", rd2[i], zd2[i], ud2[i][0], ufd2[i][0], ud2[i][1], ufd2[i][1], i);
    }
    for (i = ND2 - 1; i >= 0; i--)
    {
      if (ud2[i][1] < 50.0)
      {
        ud2m[7] = ud2[i][1];
        id2m[7] = i;
        break;
      }
    }
    for (i = id2m[7]; i > 0; i--)
    {
      if (isnan(ud2[i - 1][1]) == 1)
      {
        ud2m[6] = ud2[i][1];
        id2m[6] = i;
        break;
      }
    }
    if (i == 0)
    {
      ud2m[6] = ud2[0][1];
      id2m[6] = 0;
    }
    for (i = 0; i < ND2; i++)
    {
      if (ud2[i][1] < 50.0)
      {
        ud2m[4] = ud2[i][1];
        id2m[4] = i;
        break;
      }
    }
    for (i = id2m[4]; i < ND2 - 1; i++)
    {
      if (isnan(ud2[i + 1][1]) == 1)
      {
        ud2m[5] = ud2[i][1];
        id2m[5] = i;
        break;
      }
    }
    if (i == ND2 - 1)
    {
      ud2m[5] = ud2[ND2 - 1][1];
      id2m[5] = ND2 - 1;
    }
    if (((ud2m[4] - ud2m[6]) * (ud2m[4] - ud2m[7]) < -1.0e-12) && ((ud2m[5] - ud2m[6]) * (ud2m[5] - ud2m[7]) < -1.0e-12))
    {
      id2m[4] = 0;
      id2m[5] = 0;
      ud2m[4] = 0.0;
      ud2m[5] = 0.0;
    }
    if (((ud2m[6] - ud2m[4]) * (ud2m[6] - ud2m[5]) < -1.0e-12) && ((ud2m[7] - ud2m[4]) * (ud2m[7] - ud2m[5]) < -1.0e-12))
    {
      id2m[6] = 0;
      id2m[7] = 0;
      ud2m[6] = 0.0;
      ud2m[7] = 0.0;
    }
    /*        for(i=0;i<ND2;i++){
          if(ud2[i][1]<50.0){ud2m[2]=ud2[i][1];id2m[2]=i;break;}
        }
        for(i=ND2-1;i>=0;i--){
          if(ud2[i][1]<50.0){ud2m[3]=ud2[i][1];id2m[3]=i;break;}
        }*/
    for (i = 0; i < ND2; i++)
    {
      if ((ud2[i][0] - ud2m[4]) * (ud2[i][0] - ud2m[5]) < -1.0e-12)
      {
        ud2[i][0] = NAN;
        ufd2[i][0] = NAN;
      }
      if ((ud2[i][0] - ud2m[6]) * (ud2[i][0] - ud2m[7]) < -1.0e-12)
      {
        ud2[i][0] = NAN;
        ufd2[i][0] = NAN;
      }
    }
    /*        for(i=0;i<ND2;i++){
          if((ud2[i][0]-ud2m[2])*(ud2[i][0]-ud2m[3])<-1.0e-12){ud2[i][0]=NAN;ufd2[i][0]=NAN;}
//fprintf(fp2,"%lf %lf %lf %lf %lf %lf\n",rd2[i],zd2[i],ud2[i][0],ufd2[i][0],ud2[i][1],ufd2[i][1]);
        }*/
    for (i = 0; i < ND2; i++)
    {
      if (ud2[i][0] < 50.0)
      {
        ud2m[0] = ud2[i][0];
        id2m[0] = i;
        break;
      }
    }
    for (i = ND2 - 1; i >= 0; i--)
    {
      if (ud2[i][0] < 50.0)
      {
        ud2m[1] = ud2[i][0];
        id2m[1] = i;
        break;
      }
    }
    //printf("%lf %lf %lf %lf\n",ud2m[0],ud2m[1],ud2m[2],ud2m[3]);

    for (j = 0; j < NU; j++)
    { // go through cos\theta loop
      u = -1.0 + du * (2.0 * j + 1.0) / 2.0;
      //u=-0.5;

      if ((u - ud1m[0]) * (u - ud1m[1]) < -1.0e-12)
      {
        for (i = id1m[0]; i < id1m[1]; i++)
        {
          if ((u - ud1[i][0]) * (u - ud1[i + 1][0]) < -1.0e-12)
          {
            nnu1 = nnud1[i] + (nnud1[i + 1] - nnud1[i]) / (ud1[i + 1][0] - ud1[i][0]) * (u - ud1[i][0]);
            uf1 = ufd1[i][0] + (ufd1[i + 1][0] - ufd1[i][0]) / (ud1[i + 1][0] - ud1[i][0]) * (u - ud1[i][0]);
            fnu1 = fnu1 + nnu1 * (uf1 + 1.0) / 2.0;
            fprintf(fp2, "%6.5e %6.5e %6.5e ", u, p, nnu1 * (uf1 + 1.0) / 2.0);
            //fprintf(fp2,"%6.5e %6.5e %6.5e %6.5e %6.5e %d %lf %lf %d\n",u,p,nnu1*(uf1+1.0)/2.0,nnu1,uf1,i,rd1[i],zd1[i],1);
            break;
          }
        }
      }
      //          if((u-ud1m[2])*(u-ud1m[3])<-1.0e-12){
      else if ((u - ud1m[2]) * (u - ud1m[3]) < -1.0e-12)
      {
        for (i = id1m[2]; i < id1m[3]; i++)
        {
          if ((u - ud1[i][1]) * (u - ud1[i + 1][1]) < -1.0e-12)
          {
            nnu1 = nnud1[i] + (nnud1[i + 1] - nnud1[i]) / (ud1[i + 1][1] - ud1[i][1]) * (u - ud1[i][1]);
            uf1 = ufd1[i][1] + (ufd1[i + 1][1] - ufd1[i][1]) / (ud1[i + 1][1] - ud1[i][1]) * (u - ud1[i][1]);
            fnu1 = fnu1 + nnu1 * (uf1 + 1.0) / 2.0;
            fprintf(fp2, "%6.5e %6.5e %6.5e ", u, p, nnu1 * (uf1 + 1.0) / 2.0 / (2.0 * M_PI));
            //fprintf(fp2,"%6.5e %6.5e %6.5e %6.5e %6.5e %d %lf %lf %d\n",u,p,nnu1*(uf1+1.0)/2.0,nnu1,uf1,i,rd1[i],zd1[i],2);
            break;
          }
        }
      }
      else if ((u - ud1m[4]) * (u - ud1m[5]) < -1.0e-12)
      {
        for (i = id1m[4]; i < id1m[5]; i++)
        {
          if ((u - ud1[i][1]) * (u - ud1[i + 1][1]) < -1.0e-12)
          {
            nnu1 = nnud1[i] + (nnud1[i + 1] - nnud1[i]) / (ud1[i + 1][1] - ud1[i][1]) * (u - ud1[i][1]);
            uf1 = ufd1[i][1] + (ufd1[i + 1][1] - ufd1[i][1]) / (ud1[i + 1][1] - ud1[i][1]) * (u - ud1[i][1]);
            fnu1 = fnu1 + nnu1 * (uf1 + 1.0) / 2.0;
            fprintf(fp2, "%6.5e %6.5e %6.5e ", u, p, nnu1 * (uf1 + 1.0) / 2.0);
            //fprintf(fp2,"%6.5e %6.5e %6.5e %6.5e %6.5e %d %lf %lf %d\n",u,p,nnu1*(uf1+1.0)/2.0,nnu1,uf1,i,rd1[i],zd1[i],2);
            break;
          }
        }
      }
      else if ((u - ud1m[6]) * (u - ud1m[7]) < -1.0e-12)
      {
        for (i = id1m[6]; i < id1m[7]; i++)
        {
          if ((u - ud1[i][1]) * (u - ud1[i + 1][1]) < -1.0e-12)
          {
            nnu1 = nnud1[i] + (nnud1[i + 1] - nnud1[i]) / (ud1[i + 1][1] - ud1[i][1]) * (u - ud1[i][1]);
            uf1 = ufd1[i][1] + (ufd1[i + 1][1] - ufd1[i][1]) / (ud1[i + 1][1] - ud1[i][1]) * (u - ud1[i][1]);
            fnu1 = fnu1 + nnu1 * (uf1 + 1.0) / 2.0;
            fprintf(fp2, "%6.5e %6.5e %6.5e ", u, p, nnu1 * (uf1 + 1.0) / 2.0);
            //fprintf(fp2,"%6.5e %6.5e %6.5e %6.5e %6.5e %d %lf %lf %d\n",u,p,nnu1*(uf1+1.0)/2.0,nnu1,uf1,i,rd1[i],zd1[i],2);
            break;
          }
        }
      }
      else
      {
        fprintf(fp2, "%6.5e %6.5e %6.5e ", u, p, 0.0);
      }
      if ((u - ud2m[0]) * (u - ud2m[1]) < -1.0e-12)
      {
        for (i = id2m[0]; i < id2m[1]; i++)
        {
          if ((u - ud2[i][0]) * (u - ud2[i + 1][0]) < -1.0e-12)
          {
            nnu2 = nnud2[i] + (nnud2[i + 1] - nnud2[i]) / (ud2[i + 1][0] - ud2[i][0]) * (u - ud2[i][0]);
            uf2 = ufd2[i][0] + (ufd2[i + 1][0] - ufd2[i][0]) / (ud2[i + 1][0] - ud2[i][0]) * (u - ud2[i][0]);
            fnu2 = fnu2 + nnu2 * (uf2 + 1.0) / 2.0;
            fprintf(fp2, "%6.5e\n", nnu2 * (uf2 + 1.0) / 2.0);
            //fprintf(fp2,"%6.5e %6.5e %6.5e %6.5e %6.5e\n",u,p,nnu2*(uf2+1.0)/2.0,nnu2,uf2);
            break;
          }
        }
      }
      else if ((u - ud2m[2]) * (u - ud2m[3]) < -1.0e-12)
      {
        for (i = id2m[2]; i < id2m[3]; i++)
        {
          if ((u - ud2[i][1]) * (u - ud2[i + 1][1]) < -1.0e-12)
          {
            nnu2 = nnud2[i] + (nnud2[i + 1] - nnud2[i]) / (ud2[i + 1][1] - ud2[i][1]) * (u - ud2[i][1]);
            uf2 = ufd2[i][1] + (ufd2[i + 1][1] - ufd2[i][1]) / (ud2[i + 1][1] - ud2[i][1]) * (u - ud2[i][1]);
            fnu2 = fnu2 + nnu2 * (uf2 + 1.0) / 2.0;
            fprintf(fp2, "%6.5e\n", nnu2 * (uf2 + 1.0) / 2.0);
            //fprintf(fp2,"%6.5e %6.5e %6.5e %6.5e %6.5e\n",u,p,nnu2*(uf2+1.0)/2.0,nnu2,uf2);
            break;
          }
        }
      }
      else if ((u - ud2m[4]) * (u - ud2m[5]) < -1.0e-12)
      {
        for (i = id2m[4]; i < id2m[5]; i++)
        {
          if ((u - ud2[i][1]) * (u - ud2[i + 1][1]) < -1.0e-12)
          {
            nnu2 = nnud2[i] + (nnud2[i + 1] - nnud2[i]) / (ud2[i + 1][1] - ud2[i][1]) * (u - ud2[i][1]);
            uf2 = ufd2[i][1] + (ufd2[i + 1][1] - ufd2[i][1]) / (ud2[i + 1][1] - ud2[i][1]) * (u - ud2[i][1]);
            fnu2 = fnu2 + nnu2 * (uf2 + 1.0) / 2.0;
            fprintf(fp2, "%6.5e\n", nnu2 * (uf2 + 1.0) / 2.0);
            //fprintf(fp2,"%6.5e %6.5e %6.5e %6.5e %6.5e\n",u,p,nnu2*(uf2+1.0)/2.0,nnu2,uf2);
            break;
          }
        }
      }
      else if ((u - ud2m[6]) * (u - ud2m[7]) < -1.0e-12)
      {
        for (i = id2m[6]; i < id2m[7]; i++)
        {
          if ((u - ud2[i][1]) * (u - ud2[i + 1][1]) < -1.0e-12)
          {
            nnu2 = nnud2[i] + (nnud2[i + 1] - nnud2[i]) / (ud2[i + 1][1] - ud2[i][1]) * (u - ud2[i][1]);
            uf2 = ufd2[i][1] + (ufd2[i + 1][1] - ufd2[i][1]) / (ud2[i + 1][1] - ud2[i][1]) * (u - ud2[i][1]);
            fnu2 = fnu2 + nnu2 * (uf2 + 1.0) / 2.0;
            fprintf(fp2, "%6.5e\n", nnu2 * (uf2 + 1.0) / 2.0);
            //fprintf(fp2,"%6.5e %6.5e %6.5e %6.5e %6.5e\n",u,p,nnu2*(uf2+1.0)/2.0,nnu2,uf2);
            break;
          }
        }
      }
      else
      {
        fprintf(fp2, "%6.5e\n", 0.0);
      }
    }
    fprintf(fp2, "\n");
  }
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  fclose(fp4);
  printf("ELN profile of the point (xx = %f, zz = %f) is created.\n", xx, zz);
  return (0);
}
