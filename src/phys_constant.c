#include "phys_constant.h"
#include <stdio.h>

double G_grav = G_DEF;
double H_0 = H_0_DEF;
double c = C_DEF;

void export_phys_const(char *file_name)
{
    FILE *f;
    f = fopen(file_name,"w+");
    fprintf(f,"G=%g\n",G_grav);
    fprintf(f,"H0=%g\n",H_0);
    fprintf(f,"Gyr_per_Gpc=%g\n",Gyr_per_Gpc);
    fprintf(f,"kappa=%g\n",8*pi*G_grav/3.0);
    fprintf(f,"c=%g\n",c);
    fclose(f);
}

void set_phys_const(double G_in, double H_0_in, double c_in)
{
    G_grav = G_in;
    H_0 = H_0_in;
    c = c_in;
}
