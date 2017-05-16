double adjust3d(PointPtr scan, int num, PosturePtr initial, int target)
{
  // double gsum[6], Hsum[6][6],Hsumh[6][6],Hinv[6][6],g[6],gd[6],ge[6][6],H[6][6],hH[6][6];
  double gsum[6], Hsum[6][6], Hsumh[6][6], Hinv[6][6], g[6], H[6][6], hH[6][6];
  // double sc[3][3],sc_d[3][3][3],sc_dd[3][3][3][3],sce[3][3][3];
  double sc[3][3], sc_d[3][3][3], sc_dd[3][3][3][3];
  // double *work,*work2,*work3;
  double *work;
  double esum = 0, gnum = 0;
  NDPtr nd[8];
  NDMapPtr nd_map;
  int i, j, n, m, k, layer;
  double x, y, z;  //,sa,ca,sb,cb,sg,cg;
  PosturePtr pose;
  // Point p,pe[6],pd;
  Point p;
  PointPtr scanptr;
  // int inc,count;
  int inc;
  int ndmode;
  double dist, weight_total, weight_sum, weight_next;

  /*initialize*/
  gsum[0] = 0;
  gsum[1] = 0;
  gsum[2] = 0;
  gsum[3] = 0;
  gsum[4] = 0;
  gsum[5] = 0;
  j = 0;
  zero_matrix6d(Hsum);
  zero_matrix6d(Hsumh);
  pose = initial;

  /*�Ѵ������1����ʬʬ��ޤ�ˤβ�žʬ��׻�*/
  set_sincos(pose->theta, pose->theta2, pose->theta3, sc_d);
  set_sincos(pose->theta + E_THETA, pose->theta2, pose->theta3, sc_dd[0]);
  set_sincos(pose->theta, pose->theta2 + E_THETA, pose->theta3, sc_dd[1]);
  set_sincos(pose->theta, pose->theta2, pose->theta3 + E_THETA, sc_dd[2]);

  /*��ɸ�Ѵ���*/
  set_sincos2(pose->theta, pose->theta2, pose->theta3, sc);

  /*�켡��ʬ������Ѳ����ʤ���ʬ�η׻�*/
  qd3[0][0] = 1;
  qd3[0][1] = 0;
  qd3[0][2] = 0;

  qd3[1][0] = 0;
  qd3[1][1] = 1;
  qd3[1][2] = 0;

  qd3[2][0] = 0;
  qd3[2][1] = 0;
  qd3[2][2] = 1;
  for (n = 0; n < 6; n++)
  {
    for (m = 0; m < 6; m++)
    {
      for (k = 0; k < 3; k++)
      {
        qdd3[n][m][k] = 0;
      }
    }
  }

  //#if WEIGHTED_SELECT
  if (_downsampler_num == 0)
  {
    /*�ǡ��������Ф�����1=��ĤŤġ�*/
    switch (target)
    {
      case 3:
        inc = 1;
        ndmode = 0;
        break;
      case 2:
        inc = 500;
        ndmode = 1;
        break;
      case 1:
        inc = 5000;
        ndmode = 0;
        break;
      default:
        inc = 5000;
        ndmode = 0;
        break;
    }
  }
  //#else
  if (_downsampler_num == 1)
  {
    /*�ǡ��������Ф�����1=��ĤŤġ�*/
    switch (target)
    {
      case 3:
        inc = 1;
        ndmode = 0;
        break;
      case 2:
        inc = 1;
        ndmode = 1;
        break;
      case 1:
        inc = 1;
        ndmode = 0;
        break;
      default:
        inc = 1;
        ndmode = 0;
        break;
    }
  }
  //#endif

  scanptr = scan;

  /*����ˤĤ��Ʒ����֤��׻�*/

  //#if WEIGHTED_SELECT
  if (_downsampler_num == 0)
  {
    weight_total = scan_points_totalweight;
    ;
    weight_next = 0;
    weight_sum = 0;

    //  FILE *point_fp;
    // point_fp=fopen("/tmp/range","w");
    for (i = 0; i < num; i++)
    {
      weight_sum += scan_points_weight[i];
      if (weight_sum < weight_next)
      {
        scanptr++;
        continue;
      }

      /*���κ�ɸ�Ѵ��׻�*/
      x = scanptr->x;
      y = scanptr->y;
      z = scanptr->z;
      //    fprintf(point_fp,"%f %f %f \n",x,y,z);

      scanptr++;
      weight_next += weight_total / (double)inc;  // 1000;
      dist = 1;

      p.x = x * sc[0][0] + y * sc[0][1] + z * sc[0][2] + pose->x;
      p.y = x * sc[1][0] + y * sc[1][1] + z * sc[1][2] + pose->y;
      p.z = x * sc[2][0] + y * sc[2][1] + z * sc[2][2] + pose->z;

      /*���ϥ������ˤ�����������*/
      if (ndmode == 1)
        layer = 1;  // layer_select;
      if (ndmode == 0)
        layer = 0;  // layer_select;
      nd_map = NDmap;

      while (layer > 0)
      {
        if (nd_map->next)
          nd_map = nd_map->next;
        layer--;
      }

      /*�����б�����ND�ܥ�����������Ʊ���˼�������ND�ܥ�����򹹿���
        �٤����Τ�Ĥ������Ӥ���Ĥ�Ĥ�����*/

      if (!get_ND(nd_map, &p, nd, target))
        continue;

      /*q�ΰ켡��ʬ(�Ѳ�������Τ�)*/
      work = (double *)sc_d;
      for (m = 0; m < 3; m++)
      {
        for (k = 0; k < 3; k++)
        {
          // qd3[txtytzabg][xyz]
          qd3[m + 3][k] = x * (*work) + y * (*(work + 1)) + z * (*(work + 2));
          // x*sc_d[m][k][0] + y*sc_d[m][k][1] + z*sc_d[m][k][2];
          work += 3;
        }
      }

      /*q������ʬ���Ѳ�������Τߡ�*/
      work = (double *)sc_dd;
      for (n = 0; n < 3; n++)
      {
        for (m = 0; m < 3; m++)
        {
          for (k = 0; k < 3; k++)
          {
            qdd3[n + 3][m + 3][k] = (*work * x + *(work + 1) * y + *(work + 2) * z - qd3[m + 3][k]) / E_THETA;
            work += 3;
          }
        }
      }

      /*�����̷׻�*/
      if (nd[j])
      {
        if (nd[j]->num > 10 && nd[j]->sign == 1)
        {
          //	double e;
          esum += calc_summand3d(&p, nd[j], pose, g, hH, qd3, dist);
          add_matrix6d(Hsumh, hH, Hsumh);

          //	  dist =1;
          gsum[0] += g[0];                //*nd[j]->w;
          gsum[1] += g[1];                //*nd[j]->w;
          gsum[2] += g[2] + pose->z * 0;  //*nd[j]->w;
          gsum[3] += g[3];                //*nd[j]->w;
          gsum[4] += g[4];                //+(pose->theta2-(0.0))*1;//*nd[j]->w;
          gsum[5] += g[5];                //*nd[j]->w;
          gnum += 1;  // nd[j]->w;
        }
      }
    }
  }

  //#else
  if (_downsampler_num == 1)
  {
    for (i = 0; i < num; i += inc)
    {
      //    dist = (x*x+y*y+z*z);
      // dist *= (1.2-exp(-1*(-1 - z)*(-1 - z)/4.0));
      //    if(dist>2500)dist=2500;
      /*���κ�ɸ�Ѵ��׻�*/
      x = scanptr->x;
      y = scanptr->y;
      z = scanptr->z;
      dist = 1;
      scanptr += inc;

      p.x = x * sc[0][0] + y * sc[0][1] + z * sc[0][2] + pose->x;
      p.y = x * sc[1][0] + y * sc[1][1] + z * sc[1][2] + pose->y;
      p.z = x * sc[2][0] + y * sc[2][1] + z * sc[2][2] + pose->z;

      /*���ϥ������ˤ�����������*/
      if (ndmode == 1)
        layer = 1;  // layer_select;
      if (ndmode == 0)
        layer = 0;  // layer_select;
      nd_map = NDmap;

      while (layer > 0)
      {
        if (nd_map->next)
          nd_map = nd_map->next;
        layer--;
      }

      /*�����б�����ND�ܥ�����������Ʊ���˼�������ND�ܥ�����򹹿���
        �٤����Τ�Ĥ������Ӥ���Ĥ�Ĥ�����*/

      if (!get_ND(nd_map, &p, nd, target))
        continue;

      /*q�ΰ켡��ʬ(�Ѳ�������Τ�)*/
      work = (double *)sc_d;
      for (m = 0; m < 3; m++)
      {
        for (k = 0; k < 3; k++)
        {
          // qd3[txtytzabg][xyz]
          qd3[m + 3][k] = x * (*work) + y * (*(work + 1)) + z * (*(work + 2));
          // x*sc_d[m][k][0] + y*sc_d[m][k][1] + z*sc_d[m][k][2];
          work += 3;
        }
      }

      /*q������ʬ���Ѳ�������Τߡ�*/
      work = (double *)sc_dd;
      for (n = 0; n < 3; n++)
      {
        for (m = 0; m < 3; m++)
        {
          for (k = 0; k < 3; k++)
          {
            qdd3[n + 3][m + 3][k] = (*work * x + *(work + 1) * y + *(work + 2) * z - qd3[m + 3][k]) / E_THETA;
            work += 3;
          }
        }
      }

      /*�����̷׻�*/
      if (nd[j])
      {
        if (nd[j]->num > 10 && nd[j]->sign == 1)
        {
          //	double e;
          esum += calc_summand3d(&p, nd[j], pose, g, hH, qd3, dist);
          add_matrix6d(Hsumh, hH, Hsumh);

          //	  dist =1;
          gsum[0] += g[0];                //*nd[j]->w;
          gsum[1] += g[1];                //*nd[j]->w;
          gsum[2] += g[2] + pose->z * 0;  //*nd[j]->w;
          gsum[3] += g[3];                //*nd[j]->w;
          gsum[4] += g[4];                //+(pose->theta2-(0.0))*1;//*nd[j]->w;
          gsum[5] += g[5];                //*nd[j]->w;
          gnum += 1;  // nd[j]->w;
        }
      }
    }
  }
  //#endif

  if (gnum > 1)
  {
    //  printf("gnum=%lf\n",gnum);
    //    fclose(point_fp);
    identity_matrix6d(H);
    H[0][0] = H[0][0] / (gnum * gnum * 1000.001);
    H[1][1] = H[1][1] / (gnum * gnum * 1000.001);
    H[2][2] = H[2][2] / (gnum * gnum * 1000.001);
    H[3][3] = H[3][3] / (gnum * gnum * 0.001);
    H[4][4] = H[4][4] / (gnum * gnum * 0.001);
    H[5][5] = H[5][5] / (gnum * gnum * 0.001);

    add_matrix6d(Hsumh, H, Hsumh);

    ginverse_matrix6d(Hsumh, Hinv);

    /*----------------����------------------------*/
    pose->x -= (Hinv[0][0] * gsum[0] + Hinv[0][1] * gsum[1] + Hinv[0][2] * gsum[2] + Hinv[0][3] * gsum[3] +
                Hinv[0][4] * gsum[4] + Hinv[0][5] * gsum[5]);
    pose->y -= (Hinv[1][0] * gsum[0] + Hinv[1][1] * gsum[1] + Hinv[1][2] * gsum[2] + Hinv[1][3] * gsum[3] +
                Hinv[1][4] * gsum[4] + Hinv[1][5] * gsum[5]);
    pose->z -= (Hinv[2][0] * gsum[0] + Hinv[2][1] * gsum[1] + Hinv[2][2] * gsum[2] + Hinv[2][3] * gsum[3] +
                Hinv[2][4] * gsum[4] + Hinv[2][5] * gsum[5]);
    pose->theta -= (Hinv[3][0] * gsum[0] + Hinv[3][1] * gsum[1] + Hinv[3][2] * gsum[2] + Hinv[3][3] * gsum[3] +
                    Hinv[3][4] * gsum[4] + Hinv[3][5] * gsum[5]);
    pose->theta2 -= (Hinv[4][0] * gsum[0] + Hinv[4][1] * gsum[1] + Hinv[4][2] * gsum[2] + Hinv[4][3] * gsum[3] +
                     Hinv[4][4] * gsum[4] + Hinv[4][5] * gsum[5]);
    pose->theta3 -= (Hinv[5][0] * gsum[0] + Hinv[5][1] * gsum[1] + Hinv[5][2] * gsum[2] + Hinv[5][3] * gsum[3] +
                     Hinv[5][4] * gsum[4] + Hinv[5][5] * gsum[5]);
  }
  return esum;
}
