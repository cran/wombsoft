#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "R.h"
#include "womb_mv.h"
#include "wombling.h"

double * calc_W_rg(double *rgs_coord, int *nb_indiv,int i, int j,double *h, int nbx)
	{
	double dist;
	double *W;
	int z;
	W=(double *)malloc(*nb_indiv*sizeof(double));
	if (W!=NULL)
	 {
  for (z=0;z<*nb_indiv;z++)
		{
		*(W+z)=exp(-*(rgs_coord+z+i*(*nb_indiv)+j*(*nb_indiv*nbx))/(*h));
		}
    return(W);
   }
	if (W==NULL)
	 {
    double * erreur;
    *erreur=-1;
    return(erreur);
	 }
  }


void wombling_aflp_rg(double * coord_rg, double * data, double * coords, int * nb_indiv, int * nb_data, double * grid_x, int *nb_x, double * grid_y, int * nb_y, double * cvxhull, double *h, double * syst, double * dir)
  {
  double * XX;
  double * W;
  double * tXW;
  double * tXWX;
  double * BB;
  double * B;
  double * syst_tmp;
  double x;
  double y;
  int i,j,k;
  int l=3;
  XX=calc_X_ps(coords,nb_indiv);
  for (i=0;i<*nb_x;i++)
    {
    for (j=0;j<*nb_y;j++)
       {
      if (*(cvxhull+i+j*(*nb_x))==1)
          {
          x=*(grid_x+i);
          y=*(grid_y+j);
        //   W=calc_Q_rg(coord_rg, coords, nb_indiv,i, x,j, y,h,*nb_x);
          W=calc_W_rg(coord_rg,nb_indiv,i,j,h,*nb_x);
          if ((*XX==-1)||(*W==-1))
            {*(syst+i+j*(*nb_x))=-1;}
          else
            {
            tXW=produit_diag_transp(XX,nb_indiv,&l,W);           // supprimer taille de diag, de 2ème matrice, là et dans fct "produit_..."
            if (*tXW==-1)
              {*(syst+i+j*(*nb_x))=-2;}
            else
              {
              tXWX=produit(tXW,&l,nb_indiv,XX,nb_indiv,&l);
              if (*tXWX==-1)
                {*(syst+i+j*(*nb_x))=-3;}
              else
                {
                BB=inverse(tXWX,&l,&l);
                if (*BB==-1)
                  {*(syst+i+j*(*nb_x))=-4;}
                else
                  {
                  B=produit_restr(BB,&l,&l,tXW,&l,nb_indiv);
                  if (*B==-1)
                    {*(syst+i+j*(*nb_x))=-5;}
                  else
                    {
                    syst_tmp=systemic_aflp(data,nb_indiv,nb_data,B);
                    if (*syst_tmp==-1)
                      {*(syst+i+j*(*nb_x))=-6;}
                    else
                      {
                      *(syst+i+j*(*nb_x))=*syst_tmp;
                      *(dir+i+j*(*nb_x))=*(syst_tmp+1);
                      free(W);
                      free(tXW);
                      free(tXWX);
                      free(BB);
                      free(B);
                      free(syst_tmp);
                      }
                    }
                  }
                }
              }
            }
          }
       }
    }
   }
//---------------------------------------------------------------------------------------------------------------------------------------------


void wombling_codo_rg(double * coord_rg,double * data, double * coords, int * nr_ind, int * nb_data, int * nb_alleles_col, double * grid_x, int *nb_x, double * grid_y, int * nb_y, double * cvxhull, double *h, double * syst, double * dir)
  {
  double * XX;
  double * W;
  double * tXW;
  double * tXWX;
  double * BB;
  double * B;
  double * syst_tmp;
  double x;
  double y;
  int i,j,k,nb_alleles=0;
  int l=3;
  XX=calc_X_ps(coords,nr_ind);
  for (k=0;k<(*nb_data)/2;k++)
    { nb_alleles+=*(nb_alleles_col+k);}
  for (i=0;i<*nb_x;i++)
    {
    for (j=0;j<*nb_y;j++)
       {
        if (*(cvxhull+i+(*nb_x)*j)==1)
            {
            x=*(grid_x+i);
            y=*(grid_y+j);
            *(syst+i+j*(*nb_x))=-1;
            //W=calc_Q_rg(coord_rg, coords, nr_ind,i, x,j, y,h,*nb_x);
            W=calc_W_rg(coord_rg,nr_ind,i,j,h,*nb_x);
            if ((*XX==-1)||(*W==-1))
              {*(syst+i+j*(*nb_x))=-1;}
            else
              {
              tXW=produit_diag_transp(XX,nr_ind,&l,W);           // supprimer taille de diag, de 2ème matrice, là et dans fct "produit_..."
              if (*tXW==-1)
                {*(syst+i+j*(*nb_x))=-2;}
              else
                {
                tXWX=produit(tXW,&l,nr_ind,XX,nr_ind,&l);
                if (*tXWX==-1)
                  {*(syst+i+j*(*nb_x))=-3;}
                else
                  {
                  BB=inverse(tXWX,&l,&l);
                  if (*BB==-1)
                    {*(syst+i+j*(*nb_x))=-4;}
                  else
                    {
                    B=produit_restr(BB,&l,&l,tXW,&l,nr_ind);
                    if (*B==-1)
                      {*(syst+i+j*(*nb_x))=-5;}
                    else
                      {
                      syst_tmp=systemic_codo(data,nr_ind,nb_data,nb_alleles,nb_alleles_col,B);
                      if (*syst_tmp==-1)
                        {*(syst+i+j*(*nb_x))=-6;}
                      else
                        {
                        *(syst+i+j*(*nb_x))=*syst_tmp;
                        *(dir+i+j*(*nb_x))=*(syst_tmp+1);
                        free(W);
                        free(tXW);
                        free(tXWX);
                        free(BB);
                        free(B);
                        free(syst_tmp);
                        }
                     }
                  }
                }
              }
            }
          }
        }
     }
  }

