#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "R.h"  // ATTENTION pour que ça marche, n'inclue que RS.h pour les error et les free !!!!!! le reste : abandon !!!
#include "quantitative.h"
#include "wombling.h"

double * systemic_quantit(double * data, int * nb_indiv, int * nb_data, double * b)	// AFLP !!!!  nb_indiv=nb_col_b !!!!!
	{
	int i,j,k;
	double * beta;
	double *syst;
	double th1,th2;
	syst=(double *)malloc(2*sizeof(double));
  beta=(double *)malloc(2*sizeof(double));
  if ((syst!=NULL)&&(beta!=NULL))
  {
	*syst=0;
	*(syst+1)=0;
  for (j=0;j<*nb_data;j++)			// Parcours (i,j) de la matrice des données (i : indiv, j : allèles) !!!
		{
    for (k=0;k<2;k++)
      {*(beta+k)=0;}
		for (i=0;i<*nb_indiv;i++)
			{
			if (*(data+j*(*nb_indiv)+i)!=-1)
			{
		  *beta+=*(b+2*i)*(*(data+j*(*nb_indiv)+i));
		  *(beta+1)+=*(b+2*i+1)*(*(data+j*(*nb_indiv)+i));
      }
		  }
		*syst+=sqrt((*beta)*(*beta)+*(beta+1)*(*(beta+1)));
		*(syst+1)+=direction(*beta,*(beta+1)); 
 		}   
  *(syst+1)/=2*(*nb_data);
  free(beta);
	return(syst);
	}
	else
	 {double *erruer;
    *erruer=-1;
    return(erruer);
    }
	}
	
//---------------------------------------------------------------------------------------------------------------------------------------------

void wombling_quantit(double * data, double * coords, int * nb_indiv, int * nb_data, double * grid_x, int *nb_x, double * grid_y, int * nb_y, double * cvxhull, double *h, double * syst, double * dir)
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
  for (i=0;i<*nb_x;i++)
    { 
    for (j=0;j<*nb_y;j++)        
       { 
      if (*(cvxhull+i+j*(*nb_x))==1)
          {
          x=*(grid_x+i);
          y=*(grid_y+j);
          XX=calc_X(coords,nb_indiv,&x,&y);
          W=calc_W(coords,nb_indiv,&x,&y,h);
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
                    syst_tmp=systemic_quantit(data,nb_indiv,nb_data,B);
                    if (*syst_tmp==-1)
                      {*(syst+i+j*(*nb_x))=-6;}
                    else
                      {
                      *(syst+i+j*(*nb_x))=*syst_tmp;
                      *(dir+i+j*(*nb_x))=*(syst_tmp+1);
                      free(XX);
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

//-------------------------------------------------------------------------------
 void systemic_allelic_quantit(double * syst, double * dir, double *coords,double *data,double *h,int *data_col,double *grid_x,double *grid_y,int *nb_x,int *nb_y,int *nr_ind, double *cvxhull)
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
        XX=calc_X(coords,nr_ind,&x,&y);
        W=calc_W(coords,nr_ind,&x,&y,h);
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
                   if (*(data+(*data_col-1)*(*nr_ind)+k)!=-1)
                   {
                   beta_x+=*(B+2*k)*(*(data+(*data_col-1)*(*nr_ind)+k));
                   beta_y+=*(B+2*k+1)*(*(data+(*data_col-1)*(*nr_ind)+k));
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

