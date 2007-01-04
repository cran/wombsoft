/*# include <stdio.h> 
# include <stdlib.h>*/

double * calc_W_rg(double *rgs_coord, int *nb_indiv,int i, int j,double *h, int nbx);
void wombling_aflp_rg(double * coord_rg, double * data, double * coords, int * nb_indiv, int * nb_data, double * grid_x, int *nb_x, double * grid_y, int * nb_y, double * cvxhull, double *h, double * syst, double * dir);
void wombling_codo_rg(double * coord_rg,double * data, double * coords, int * nr_ind, int * nb_data, int * nb_alleles_col, double * grid_x, int *nb_x, double * grid_y, int * nb_y, double * cvxhull, double *h, double * syst, double * dir);
