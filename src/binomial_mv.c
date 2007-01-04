#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "R.h"  // ATTENTION pour que ça marche, n'inclue que RS.h pour les error et les free !!!!!! le reste : abandon !!!
#include "womb_mv.h"
#include "binomial_mv.h"
#include "wombling.h"

// une fonction pour calculer la fonction systémique d'un allèle, et une autre qui fait la boucle sur les allèles


void systemic_allelic_codom_rg(double * syst, double * dir, double *coords, double * rg, double *data,double *h,int *allele_nb,int *allele_col,double *grid_x,double *grid_y,int *nb_x,int *nb_y,int *nr_ind, double *cvxhull)
  {
  int i,j,k;
  int l=3;
  double * XX;
  double * W;
  double * tXW;
  double * tXWX;
  double * BB;
  double * B;
  double beta_x,beta_y;
  double x;
  double y;
  int v=2*(*nr_ind);
//  XX=calc_X_ps(coords,nr_ind);
  for (i=0; i<*nb_x;i++)
    {
    for (j=0;j<*nb_y;j++)
      {
      beta_x=0;
      beta_y=0;
      *(syst+i+j*(*nb_x))=0;
      *(dir+i+j*(*nb_x))=0;
      if (*(cvxhull+i+j*(*nb_x))==1)
        {
        x=*(grid_x+i);
        y=*(grid_y+j);
        W=calc_W_rg(rg,nr_ind,i,j,h,*nb_x);
	XX=calc_X(coords,nr_ind,&x,&y);
         if ((*XX==-1)||(*W==-1))
          {*(syst+i+j*(*nb_x))=-1;}
        else
          {
          tXW=produit_diag_transp(XX,nr_ind,&l,W);           // supprimer taille de diag, de 2ème matrice, là et dans fct "produit_..."
          if (*tXW==-1)
            {*(syst+i+j*(*nb_x))=-1;}
          else
            {
            tXWX=produit(tXW,&l,nr_ind,XX,nr_ind,&l);
            if (*tXWX==-1)
              {*(syst+i+j*(*nb_x))=-1;}
            else
              {
              BB=inverse(tXWX,&l,&l);
              if (*BB==-1)
                {*(syst+i+j*(*nb_x))=-1;}
              else
                {
                B=produit_restr(BB,&l,&l,tXW,&l,nr_ind);
                if (*B==-1)
                  {*(syst+i+j*(*nb_x))=-1;}
                else
                  {
                  for (k=0;k<*nr_ind;k++)
                   {
                   if (*(data+(*allele_col-1)*v+k)==*allele_nb)
                     {
                     beta_x+=*(B+2*k);
                     beta_y+=*(B+2*k+1);
                     }
                   }
                  for (k=*nr_ind;k<v;k++)
                   {
                   if (*(data+(*allele_col-1)*v+k)==*allele_nb)
                     {
                     beta_x+=*(B+2*(k-*nr_ind));
                     beta_y+=*(B+2*(k-*nr_ind)+1);
                     }
                   }
                   free(XX);
                   free(W);
                   free(tXW);
                   free(tXWX);
                   free(BB);
                   free(B);
                  }
                }
              }
            }
          }
        *(syst+i+j*(*nb_x))=sqrt((beta_x*beta_x)+(beta_y*beta_y));
        *(dir+i+j*(*nb_x))=direction(beta_x,beta_y);
        }
      }
    }
  }



  
void systemic_allelic_dom_rg(double * syst_p, double * dir_p, double * syst_a, double * dir_a, double *coords, double * rg, double *data, double *h, int *allele_col, double *grid_x, double *grid_y, int *nb_x, int *nb_y, int *nr_ind, double *cvxhull)
  {                     // on va gérer colonne par colonne, les deux allèles à chaque fois !  la colonne concernée est "allele_col".
  int i,j,k;            // syst, dir _p : allèle présence. syst, dir _a : allèle absence.
  int l=3;
  double * XX;
  double * W;
  double * tXW;
  double * tXWX;
  double * BB;
  double * B;
  double beta_xa,beta_ya,beta_xp,beta_yp;
  double x;
  double y;
//  XX=calc_X_ps(coords,nr_ind);
  for (i=0; i<*nb_x;i++)
    {
    for (j=0;j<*nb_y;j++)
      {
      beta_xa=0;
      beta_ya=0;
      beta_xp=0;
      beta_yp=0;
      *(syst_a+i+j*(*nb_x))=0;
      *(dir_a+i+j*(*nb_x))=0;
      *(syst_p+i+j*(*nb_x))=0;
      *(dir_p+i+j*(*nb_x))=0;
      if (*(cvxhull+i+j*(*nb_x))==1)
        {
        x=*(grid_x+i);
        y=*(grid_y+j);
        XX=calc_X(coords,nr_ind,&x,&y);
        W=calc_W_rg(rg,nr_ind,i,j,h,*nb_x);
        if ((*XX==-1)||(*W==-1))
          {*(syst_a+i+j*(*nb_x))=-1;}
        else
          {
          tXW=produit_diag_transp(XX,nr_ind,&l,W);           // supprimer taille de diag, de 2ème matrice, là et dans fct "produit_..."
          if (*tXW==-1)
            {*(syst_a+i+j*(*nb_x))=-1;}
          else
            {
            tXWX=produit(tXW,&l,nr_ind,XX,nr_ind,&l);
            if (*tXWX==-1)
              {*(syst_a+i+j*(*nb_x))=-1;}
            else
              {
              BB=inverse(tXWX,&l,&l);
              if (*BB==-1)
                {*(syst_a+i+j*(*nb_x))=-1;}
              else
                {
                B=produit_restr(BB,&l,&l,tXW,&l,nr_ind);
                if (*B==-1)
                  {*(syst_a+i+j*(*nb_x))=-1;}
                else
                  {
                  for (k=0;k<*nr_ind;k++)
                   {
                   if (*(data+(*allele_col-1)*(*nr_ind)+k)==0)
                     {
                     beta_xa+=*(B+2*k);
                     beta_ya+=*(B+2*k+1);
                     }
                   if (*(data+(*allele_col-1)*(*nr_ind)+k)==1)
                     {
                     beta_xp+=*(B+2*k);
                     beta_yp+=*(B+2*k+1);
                     }
                    }
                   free(XX);
                   free(W);
                   free(tXW);
                   free(tXWX);
                   free(BB);
                   free(B);
                   }
                }
              }
            }
          }
        *(syst_a+i+j*(*nb_x))=sqrt((beta_xa*beta_xa)+(beta_ya*beta_ya));
        *(dir_a+i+j*(*nb_x))=direction(beta_xa,beta_ya);
        *(syst_p+i+j*(*nb_x))=sqrt((beta_xp*beta_xp)+(beta_yp*beta_yp));
        *(dir_p+i+j*(*nb_x))=direction(beta_xp,beta_yp);
        }
      }
    }
  }
  
