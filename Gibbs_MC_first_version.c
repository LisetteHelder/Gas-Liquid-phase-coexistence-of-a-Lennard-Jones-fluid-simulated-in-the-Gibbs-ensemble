#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3
#define N 513

/* Initialization variables */
const int    mc_steps      = 100000;
const int    output_steps  = 1000;
const double overall_density= 0.2; // you may actually not have to use this one 
const double delta         = 0.1; // to be tuned 
const double r_cut         = 2.5;
const double Temperature=2;
const double beta=1.0/Temperature;
const char*  file_g = "gas.dat";
const char*  file_l = "liquid.dat";


/* Simulation variables */
int n_gas;
int n_liq;
int n_tot;
double e_cut;
double r_gas[N][NDIM];
double r_liq[N][NDIM];
double box_g[NDIM], box_l[NDIM];

double energy = 0.0;
double virial = 0.0;



void read_data(void){
    FILE* fp_1 = fopen(file_g, "r");
    FILE* fp_2 = fopen(file_l, "r");

    int n, d;
    double dmin,dmax;

    fscanf(fp_1, "%d\n", &n_gas);
    for(d = 0; d < NDIM; ++d){
        fscanf(fp_1, "%lf %lf\n", &dmin, &dmax);
        box_g[d] = fabs(dmax-dmin);
    }
    for(n = 0; n < n_gas; ++n){
        for(d = 0; d < NDIM; ++d) fscanf(fp_1, "%lf\t", &r_gas[n][d]);
        double diameter;
        fscanf(fp_1, "%lf\n", &diameter);
    }

    fscanf(fp_2, "%d\n", &n_liq);
    for(d = 0; d < NDIM; ++d){
        fscanf(fp_2, "%lf %lf\n", &dmin, &dmax);
        box_l[d] = fabs(dmax-dmin);
    }
    for(n = 0; n < n_liq; ++n){
        for(d = 0; d < NDIM; ++d) fscanf(fp_2, "%lf\t", &r_liq[n][d]);
        double diameter;
        fscanf(fp_2, "%lf\n", &diameter);
    }

    n_tot=n_gas+n_liq;

    fclose(fp_1);
    fclose(fp_2);
}

int displacement(void){
    // our code goes here 
    return 0;
}

int displacement(void){
    // our code goes here 
    return 0;
}

int displacement(void){
    // our code goes here 
    return 0;
}

void write_data(int step){
    char buffer_1[128];
    sprintf(buffer_1, "coords_gas_step%07d.dat", step);
    FILE* fp_1 = fopen(buffer_1, "w");
    int d, n;
    fprintf(fp_1, "%d\n", n_gas);
    for(d = 0; d < NDIM; ++d){
        fprintf(fp_1, "%lf %lf\n",0.0,box_g[d]);
    }
    for(n = 0; n < n_gas; ++n){
        for(d = 0; d < NDIM; ++d) fprintf(fp_1, "%f\t", r_gas[n][d]);
        fprintf(fp_1, "%lf\n", 1.0);
    }
    fclose(fp_1);

     char buffer_2[128];
    sprintf(buffer_2, "coords_liq_step%07d.dat", step);
    FILE* fp_2 = fopen(buffer_2, "w");

    fprintf(fp_2, "%d\n", n_liq);
    for(d = 0; d < NDIM; ++d){
        fprintf(fp_2, "%lf %lf\n",0.0,box_l[d]);
    }
    for(n = 0; n < n_liq; ++n){
        for(d = 0; d < NDIM; ++d) fprintf(fp_2, "%f\t", r_liq[n][d]);
        fprintf(fp_2, "%lf\n", 1.0);
    }
    fclose(fp_2);

}

void set_overall_density(void){ // function to be used after read data() if necessary

    double volume_g = 1.0;
    int d, n;
    for(d = 0; d < NDIM; ++d) volume_g *= box_g[d];
    double volume_l = 1.0;
    for(d = 0; d < NDIM; ++d) volume_l *= box_l[d];
    double volume= volume_g + volume_l;
    

    double target_volume = (n_tot) / overall_density;
    double scale_factor = pow(target_volume / volume, 1.0 / NDIM);

    for(n = 0; n < n_gas; ++n){
        for(d = 0; d < NDIM; ++d) r_gas[n][d] *= scale_factor;
    }
    for(d = 0; d < NDIM; ++d) box_g[d] *= scale_factor;
    
    for(n = 0; n < n_liq; ++n){
        for(d = 0; d < NDIM; ++d) r_liq[n][d] *= scale_factor;
    }
    for(d = 0; d < NDIM; ++d) box_l[d] *= scale_factor;
}

int main(int argc, char* argv[]){

    assert(delta > 0.0);

    e_cut = 4.0 * (pow(1.0 / r_cut, 12.0) - pow(1.0 / r_cut, 6.0));

    read_data();

    if(n_gas == 0){
        printf("Error: Number of particles, n_gas = 0.\n");
        return 0;
    }
        if(n_liq == 0){
        printf("Error: Number of particles, n_liq = 0.\n");
        return 0;
    }

    int d;
    for(d = 0; d < NDIM; ++d) assert(r_cut <= 0.5 * box_g[d]);
    for(d = 0; d < NDIM; ++d) assert(r_cut <= 0.5 * box_l[d]);
    int step, n;

    size_t seed = time(NULL);
    dsfmt_seed(seed);

    double volume_g = 1.0;
    int d, n;
    for(d = 0; d < NDIM; ++d) volume_g *= box_g[d];
    double volume_l = 1.0;
    for(d = 0; d < NDIM; ++d) volume_l *= box_l[d];
    double volume= volume_g + volume_l;

    printf("Starting volume: %f\n", volume);
    printf("Starting energy: %f\n", energy); // to be completed 
    printf("Starting virial: %f\n", virial); // to be completed 
    printf("Starting seed: %lu\n", seed);

    FILE* fp = fopen("measurements.dat", "w");

    // make a cycle where the 3 types of MC steps are performed after being extracted with a certain probability 
// compute the relevant physical quantities and print them on measurements.dat file 
// maybe print the files with the boxes (you have the functions)

/* possibly rememebering 
        if(step % output_steps == 0){
            printf("Step %d. Move acceptance: %f.\n",
                step, (double)accepted / (n_particles * output_steps)
            );
            accepted = 0;
            write_data(step);
        }
    }
*/



    fclose(fp);

    return 0;
}
