/* # include <stdio.h> */
/* # include <stdlib.h> */
/* # include "wombling.h" */
/* # include "wombling.c" */

double * systemic_quantit(double * data, int * nb_indiv, int * nb_data, double * b);
void wombling_quantit(double * data, double * coords, int * nb_indiv, int * nb_data, double * grid_x, int *nb_x, double * grid_y, int * nb_y, double * cvxhull, double *h, double * syst, double * dir);
void systemic_allelic_quantit(double * syst, double * dir, double *coords,double *data,double *h,int *data_col,double *grid_x,double *grid_y,int *nb_x,int *nb_y,int *nr_ind, double *cvxhull);

