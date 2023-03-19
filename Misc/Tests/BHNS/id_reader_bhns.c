#include <stdio.h>
#include <assert.h>
#include "elliptica_id_reader_lib.h"


int main(int argn, char **argv)
{
    assert(argn == 2 && "$./me /path/to/checkpoint/file");
    
    const char *f = argv[1]; /* /path/to/checkpoint/file */
    
    Elliptica_ID_Reader_T *idr = elliptica_id_reader_init(f,"generic");
    idr->ifields = "alpha,adm_gxx,adm_gxy";
    idr->npoints = 2;
    double x[2]  = {0.1 ,  5};
    double y[2]  = {-0.6,  3};
    double z[2]  = {-2  , -1};
    idr->x_coords = x;
    idr->y_coords = y;
    idr->z_coords = z;
    idr->set_param("BHNS_filler_method","ChebTn_Ylm_perfect_s2",idr);
    idr->set_param("ADM_B1I_form","zero",idr);
    
    elliptica_id_reader_interpolate(idr);
    
    printf("\nsome tests:\n");
    
    printf("alpha(%g,%g,%g) = %g\n",
            x[0],y[0],z[0],idr->field[idr->indx("alpha")][0]);

    printf("alpha(%g,%g,%g) = %g\n",
            x[1],y[1],z[1],idr->field[idr->indx("alpha")][1]);

    printf("adm_gxx(%g,%g,%g) = %g\n",
            x[0],y[0],z[0],idr->field[idr->indx("adm_gxx")][0]);

    printf("adm_gxy(%g,%g,%g) = %g\n",
            x[1],y[1],z[1],idr->field[idr->indx("adm_gxy")][1]);

    printf ("BHNS_angular_velocity = %g\n",
            idr->get_param_dbl("BHNS_angular_velocity",idr));

    elliptica_id_reader_free(idr);
}
