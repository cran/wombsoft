/* # include <stdio.h> */
/* # include <stdlib.h> */

double * calc_W(double *coord, int *nb_indiv,double *x,double *y,double *h);
double * calc_X(double *coord, int *nb_indiv, double *x, double *y);
double * produit(double *a, int *nrow_a, int *ncol_a, double *b, int *nrow_b, int *ncol_b);
double * produit_restr(double * a, int * nrow_a, int * ncol_a, double * b, int * nrow_b, int * ncol_b);
double * produit_diag_transp(double *x, int * nrow_x, int * ncol_x, double * diag) ;
double * inverse(double * a, int *nrow, int * ncol);
double direction(double x, double y);
double * systemic_aflp(double * data, int * nb_indiv, int * nb_data, double * b);
void wombling_aflp(double * data, double * coords, int * nb_indiv, int * nb_data, double * grid_x, int *nb_x, double * grid_y, int * nb_y, double * cvxhull, double *h, double * syst, double * dir);
void wombling_codo(double * data, double * coords, int * nb_indiv, int * nb_data, int * nb_alleles_col, double * grid_x, int *nb_x, double * grid_y, int * nb_y, double * cvxhull, double *h, double * syst, double * dir);
double * systemic_codo(double * data, int * nb_indiv, int * nb_data, int nb_alleles, int * nb_alleles_col, double * b);

