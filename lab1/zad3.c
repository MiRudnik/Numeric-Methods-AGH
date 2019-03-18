#include <stdio.h>
#include <gsl/gsl_ieee_utils.h>

// link with -lm -lgsl -lgslcblas

int main(int argc, char **argv){

    float f = 1.0;
    double d = 1.0;

    for(int i=0;i<100;i++){
        
        printf("Przejscie %d:\nf = ",i+1);
        gsl_ieee_printf_float(&f);
        printf("\nd = ");
        gsl_ieee_printf_double(&d);
        printf("\n");
        f = f / 10.0;
        d = d / 10.0;
    }
    return 0;
}