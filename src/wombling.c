 #include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "R.h"  
#include "wombling.h"

double * calc_W(double *coord, int *nb_indiv,double *x,double *y,double *h)
	{
	double dist;
	double *W;
	int i;
	W=(double *)malloc(*nb_indiv*sizeof(double));
	if (W!=NULL)
	 {
  for (i=0;i<*nb_indiv;i++)
		{
		dist=(*(coord+i)-*x)*(*(coord+i)-*x)+(*(coord+*nb_indiv+i)-*y)*(*(coord+*nb_indiv+i)-*y);
		*(W+i)=exp(-dist/(2*(*h)*(*h)));
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

//---------------------------------------------------------------------------------------------------------------------------------------------

/*void essai(double *x, int * nrow_x, int * ncol_x, double *res)
	{
	int i,j;
	double * p;
	p=calc_W(x,nrow_x,ncol_x);
	for (i=0;i<*ncol_x;i++)
		{
		for (j=0; j<*nrow_x;j++)
			{
			*(res+j*(*ncol_x)+i)=*(p+j*(*ncol_x)+i);
			}
		}
	}*/

//---------------------------------------------------------------------------------------------------------------------------------------------


double * calc_X(double *coord, int *nb_indiv, double *x, double *y)
	{
	int i;
	double *X;
	X=(double *)malloc(*nb_indiv*3*sizeof(double));
	if (X!=NULL)
	 {
    for (i=0;i<*nb_indiv; i++)
		  {
  		*(X+i)=1;
	   	*(X+i+*nb_indiv)=*(coord+i)-*x;
	 	  *(X+i+2*(*nb_indiv))=*(coord+i+*nb_indiv)-*y;
		  }
  return(X);
   }
  if (X==NULL)
    {
    double * erreur;
    *erreur=-1;
    return(erreur);
    }
	}

//---------------------------------------------------------------------------------------------------------------------------------------------

double * produit(double *a, int *nrow_a, int *ncol_a, double *b, int *nrow_b, int *ncol_b)
	{
	int i,j,k;
	double *c;
	double g,r;
	c=(double *)malloc(*nrow_a*(*ncol_b)*sizeof(double));
	if (c!=NULL)
	 {
    for (i=0;i<*nrow_a;i++)
	   {
	     for (j=0;j<*ncol_b;j++)
	       {
          *(c+i+j*(*nrow_a))=0;
         }
     }
    if (*ncol_a==*nrow_b)
		 {
		  for (i=0; i<*nrow_a;i++)
			 {
			 for (j=0; j<*ncol_b;j++)
		    {
				for (k=0;k<*nrow_b;k++)
					{*(c+i+j*(*nrow_a))+=*(a+i+k*(*nrow_a))*(*(b+k+j*(*nrow_b)));}
				}
			 }
		  }
    return(c);
   }
  if (c==NULL)
    {
      double * erruer;
      *erruer=-1;
      return(erruer);
    }
	}
//---------------------------------------------------------------------------------------------------------------------------------------------

double * produit_restr(double * a, int * nrow_a, int * ncol_a, double * b, int * nrow_b, int * ncol_b)
	{
	int i,j,k;
	int x,y,z,t;
	double *c;
	x=*nrow_a;
	y=*nrow_b;
	z=*ncol_a;
	t=*ncol_b;
	double g,r;
	c=(double *)malloc(2*t*sizeof(double));
	if (c!=NULL)
	 {
    if (z==y)
		  {
		  for (i=1; i<3;i++)
			 {
			 for (j=0; j<t;j++)
			 	{
			 	*(c+(i-1)+j*2)=0;
				for (k=0;k<y;k++)
					{g=*(a+i+k*x);
					r=*(b+k+j*y);
					*(c+(i-1)+j*2)+=g*r;}
				}
			 }
		  }
   	return(c);
    }
  if (c==NULL)
    {
    double * erruer;
    *erruer=-1;
    return(erruer);
    }
	}

//---------------------------------------------------------------------------------------------------------------------------------------------

/*void produit_diag(double *a, int *nrow_a, int *ncol_a, double *diag, int *taille_diag)      //  ATTENTION : transforme la matrice d'entrée !!!
	{
	int i,j;
	if (*ncol_a==*taille_diag)
		{
		for (i=0;i<*nrow_a;i++)
			for (j=0;j<*ncol_a;j++)
				{
				*(a+i+j*(*nrow_a))*=*(diag+j);
				}
		}
	} */

double * produit_diag_transp(double *x, int * nrow_x, int * ncol_x, double * diag)
	{
	int i,j;
	double *c;
	c=(double *)malloc(*nrow_x*(*ncol_x)*sizeof(double));
	if (c!=NULL)
	 {
   for (i=0;i<*nrow_x;i++)
	   {
      for (j=0;j<*ncol_x;j++)
       {
        *(c+i+j*(*nrow_x))=0;
       }
     }
	 for (i=0;i<*nrow_x; i++)
		{
		for (j=0;j<*ncol_x;j++)
			{
			*(c+i*(*ncol_x)+j)=*(diag+i)*(*(x+j*(*nrow_x)+i));
			}
		}
	 return(c);
   }
  if (c==NULL)
    {
     double *erruer;
     *erruer=-1;
     return(erruer);
    }
	}

//---------------------------------------------------------------------------------------------------------------------------------------------

double * inverse(double * a, int *nrow, int * ncol) // SI LA MATRICE N'EST PAS INVERSIBLE, RESSORT L'IDENTITE !!!!!!!!!!!!!!!!!!!!!!!
	{
	double * id;
	double tmp,p;
	int i,j,k;
	if (*nrow==*ncol)
		{
		id=(double *)malloc(*nrow*(*nrow)*sizeof(double));
		if (id != NULL)
		{
		for (i=0;i<*nrow;i++)				    // création de la matrice identité
			{
			*(id+i*(*nrow)+i)=1;
			for (j=0;j<*nrow;j++)
				{
				if (j!=i)
					{*(id+j*(*nrow)+i)=0;}
				}
			}				                     	// Id est prête !!!!!!!!!
		for (j=0;j<*nrow-1;j++)         // pour toutes les lignes
			{
			i=j;
			for (i=j;(i<*nrow)&&(*(a+j*(*nrow)+i)==0);i++)   // recherche du pivot : ligne i
				{
				}
			if (i==*nrow)                  // <=> matrice non inversible
				{                            //
				for (i=0;i<*nrow;i++)				 // création de la matrice identité
					{                          //
					*(id+i*(*nrow)+i)=1;       //
					for (j=0;j<*nrow;j++)      //
						{                        //
						if (j!=i)                //
							{*(id+j*(*nrow)+i)=0;} //
						}                        //
					}				                   //
				return(id);                  // on retourne la matrice Id
				}
			else
				{
				p=*(a+i+j*(*nrow));          // pivot
				for (k=0;k<*nrow;k++)                  //
					{                                    //
					tmp=*(a+k*(*nrow)+i);                //
					*(a+k*(*nrow)+i)=*(a+k*(*nrow)+j);   //    échange des lignes i et j
					*(a+k*(*nrow)+j)=tmp/p;              //    et division par le pivot
					tmp=*(id+k*(*nrow)+i);               //
					*(id+k*(*nrow)+i)=*(id+k*(*nrow)+j); //
					*(id+k*(*nrow)+j)=tmp/p;             //
					}                                    //
				}
			for (i=j+1;i<*nrow;i++)
				{                                            //
				p=*(a+i+j*(*nrow));                          //
				for (k=0; k<*nrow;k++)                       //
					{                                          //   on met des 0 sur le bas
					*(a+i+k*(*nrow))-=*(a+j+k*(*nrow))*p;      //   de la colonne
					*(id+i+k*(*nrow))-=*(id+j+k*(*nrow))*p;    //
					}                                          //
				}                                            //
			}                                              //
		p=*(a+(*nrow-1)*(*nrow+1));
		if (p!=0)
		  {
  		*(a+(*nrow-1)*(*nrow+1))=1;                      // cas particulier de la dernière colonne
	    for (k=0;k<*nrow;k++)                            //
			 {                                              //
			 *(id+k*(*nrow)+*nrow-1)/=p;                    //
			 }				//OUF, y'a des zéros en bas !!! (confondu lignes et colonnes bien sûr !!!)
			}
    if(p==0)                         // <=> la matrice n'est pas inversible
      {for (i=0;i<*nrow;i++)				 //
					{                          //
					*(id+i*(*nrow)+i)=1;       //
					for (j=0;j<*nrow;j++)      //   création de la matrice identité
						{                        //
						if (j!=i)                //
							{*(id+j*(*nrow)+i)=0;} //
						}                        //
					}				                   //
				return(id);                  //    on retourne la matrice Id
				}
// PAUSE !!!!!!!!!!!!!!!!!!!!!!!!
		for (j=*nrow-1;j>0;j--)
			{
			for (i=0;i<j;i++)
				{
				p=*(a+j*(*nrow)+i);
				for (k=0;k<*nrow;k++)                         // 0 au dessus de la diagonale
					{
					*(a+k*(*nrow)+i)-=*(a+k*(*nrow)+j)*p;
					*(id+k*(*nrow)+i)-=*(id+k*(*nrow)+j)*p;
					}
				}
			}
	return(id);
		}
	if (id==NULL)
	 {
    double * erreur;
    *erreur=-1;
    return(erreur);
   }
  }
	else
	  {
    double * erreur;
    *erreur=-1;
    return(erreur);
    }
  }

//---------------------------------------------------------------------------------------------------------------------------------------------

double direction(double x, double y)
	{
	double pi,tmp;
  pi=3.1415926;
	if (x==0)
		{
		if (y>0)
			{return(pi);}
		else
			{return(-pi);}
		}
	if (x>0)
		{
    tmp=atan(y/x);
    if (2*tmp>pi)
      {
      return(2*tmp-2*pi);
      }
    else
      {
      return(2*tmp);
      }
    }
	if (x<0)
		{
    tmp=atan(y/x)+pi;
    if (2*tmp>pi)
      {
      return(2*tmp-2*pi);
      }
    else
      {
      return(2*tmp);
      }
    }
	}

//---------------------------------------------------------------------------------------------------------------------------------------------

double * systemic_aflp(double * data, int * nb_indiv, int * nb_data, double * b)	// AFLP !!!!  nb_indiv=nb_col_b !!!!!
	{
	int i,j,k;
	double * beta;
	double *syst;
	double th1,th2;
	syst=(double *)malloc(2*sizeof(double));
  beta=(double *)malloc(4*sizeof(double));
  if ((syst!=NULL)&&(beta!=NULL))
  {
	*syst=0;
	*(syst+1)=0;
  for (j=0;j<*nb_data;j++)			// Parcours (i,j) de la matrice des données (i : indiv, j : allèles) !!!
		{
    for (k=0;k<4;k++)
      {*(beta+k)=0;}
		for (i=0;i<*nb_indiv;i++)
			{
			if (*(data+i+j*(*nb_indiv))==0)
				{
				(*beta)+=*(b+2*i);
		    (*(beta+1))+=*(b+2*i+1);
				}
			if (*(data+i+j*(*nb_indiv))==1)
				{
				(*(beta+2))+=*(b+2*i);
				(*(beta+3))+=*(b+2*i+1);
				}
			}
			*syst+=sqrt((*beta)*(*beta)+*(beta+1)*(*(beta+1)))+sqrt(*(beta+2)*(*(beta+2))+(*(beta+3))*(*(beta+3)));
			*(syst+1)+=direction(*beta,*(beta+1))+direction(*(beta+2),*(beta+3));
		}
  *(syst+1)/=4*(*nb_data);
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

void wombling_aflp(double * data, double * coords, int * nb_indiv, int * nb_data, double * grid_x, int *nb_x, double * grid_y, int * nb_y, double * cvxhull, double *h, double * syst, double * dir)
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
                    syst_tmp=systemic_aflp(data,nb_indiv,nb_data,B);
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

//---------------------------------------------------------------------------------------------------------------------------------------------

void wombling_codo(double * data, double * coords, int * nr_ind, int * nb_data, int * nb_alleles_col, double * grid_x, int *nb_x, double * grid_y, int * nb_y, double * cvxhull, double *h, double * syst, double * dir)
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
            *(syst+i+j*(*nb_x))=0;
            XX=calc_X(coords,nr_ind,&x,&y);
            W=calc_W(coords,nr_ind,&x,&y,h);
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

//---------------------------------------------------------------------------------------------------------------------------------------------

double * systemic_codo(double * data, int * nr_ind, int * nb_data, int nb_alleles, int * nb_alleles_col, double * b)     // nb indiv = le faux : nb_lignes !!!
  {                                                                                                     // nb_data = le vrai (=2*nb_loci=2*nb_colonnes)
   double *beta;                                                                                        // nb_alleles = nb total
   double *syst;
   int v;
   v=2*(*nr_ind);                                                                                      // nb_alleles_col = vecteur des nbs d'allèles différents dans chaque colonne.
   syst=(double *)malloc(2*sizeof(double));
   if (syst!=NULL)
    {
    *syst=0;
    *(syst+1)=0;
    int l,i,j,k,m=1;       // l = le numéro de l'allèle considéré !    //m le premier de la colonne, pour faire -m et numéroter de 0 à nj la colonne j (logique)
    for (j=0;j<((*nb_data)/2);j++)
      {
      beta=(double *)malloc(*(nb_alleles_col+j)*2*sizeof(double));
      if (beta!=NULL)
        {
        for(i=0;i<(*(nb_alleles_col+j)*2);i++)
          {*(beta+i)=0;}
        for(i=0;i<(*nr_ind);i++)
          {
          if ((*(data+i+j*v))!=0)
            {
            l=*(data+i+j*v)-m;
            *(beta+2*l)+=*(b+2*i);
            *(beta+2*l+1)+=*(b+2*i+1);
            }   
          }         
        for(i=(*nr_ind);i<v;i++)
          {
          if ((*(data+i+j*v))!=0)
            {
            l=*(data+i+j*v)-m;
            *(beta+2*l)+=*(b+2*(i-v/2));
            *(beta+2*l+1)+=*(b+2*(i-v/2)+1);
            }     
          } 
         for(k=0;k<(*(nb_alleles_col+j));k++)
          {
          *syst+=sqrt(*(beta+2*k)*(*(beta+2*k))+*(beta+2*k+1)*(*(beta+2*k+1)));
          *(syst+1)+=direction(*(beta+2*k),*(beta+2*k+1));
          }
        free(beta);
        m+=*(nb_alleles_col+j);
        }
      else
        {
        double *erruer;
        *erruer=-1;
        return(erruer);
        }
      }
    *(syst+1)/=2*nb_alleles;
    return(syst);
    }
  else
    {
    double *erruer;
    *erruer=-1;
    return(erruer);
    }
  }

//---------------------------------------------------------------------------------------------------------------------------------------------


/*void essai_syst(double *x, int * nrow_x, int * ncol_x, double * y, double *res)
	{
	int i,j;
	double * p;
	p=systemic(x,nrow_x,ncol_x,y);
/*	for (i=0;i<2;i++)
		{
		for (j=0; j<*ncol_y;j++)
			{
			*(res+j*2+i)=*(p+j*2+i);
			}
		}*/
 /* *res=*p;
  *(res+1)=*(p+1);
	}
     */

//---------------------------------------------------------------------------------------------------------------------------------------------
/* void main()
	{
	int i,j;
	double r;
	double *essai2, *essai;
	essai=(double *)malloc(9*sizeof(double));
	essai2=(double *)malloc(3*sizeof(double));
	for (i=0;i<3;i++)
		{
		for (j=0; j<3; j++)
			{*(essai+3*i+j)=i+j+1;
			printf("%f",*(essai+3*i+j));}
		*(essai2+i)=i;
		}
	*(essai+8)=0;
	produit_diag(essai,3,3,essai2,3);
	for (i=0;i<3;i++)
		for (j=0; j<3; j++)
			printf("%f",*(essai+3*i+j));
	scanf("%f",r);
	free(essai);
	free(essai2);
	}
*/
//---------------------------------------------------------------------------------------------------------------------------------------------
