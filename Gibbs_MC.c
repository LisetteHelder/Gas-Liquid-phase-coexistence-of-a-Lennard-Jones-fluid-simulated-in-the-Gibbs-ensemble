#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3
#define N 513

const int    mc_steps        = 100000;
const int    output_steps    = 1000;
const double overall_density = 0.2;
const double delta           = 0.1;
const double r_cut           = 2.5;
const double Temperature     = 2.0;
const double beta            = 1.0 / Temperature;

typedef struct {
    int n;
    double energy; 
    double r[N][NDIM];
    double box[NDIM];
    const char* input_file;
    const char* label;
} Box;

Box gas = { .n = 0, .energy=0, .input_file = "gas.dat",    .label = "gas" };
Box liq = { .n = 0, .energy=0, .input_file = "liquid.dat", .label = "liq" };


double E_tot=0 ;
double n_tot=0; 
double e_cut;
double dV = 0.0;

double box_volume(const Box* b)
{
    double volume = 1.0;

    for(int d = 0; d < NDIM; ++d){
        volume *= b->box[d];
    }

    return volume;
}

void read_box(Box* b)
{
    FILE* fp = fopen(b->input_file, "r");
    if(fp == NULL){
        fprintf(stderr, "Error: could not open input file %s\n", b->input_file);
        exit(EXIT_FAILURE);
    }

    if(fscanf(fp, "%d\n", &b->n) != 1){
        fprintf(stderr, "Error: could not read number of particles from %s\n", b->input_file);
        exit(EXIT_FAILURE);
    }

    if(b->n < 0 || b->n > N){
        fprintf(stderr, "Error: file %s contains n = %d, but maximum allowed is N = %d and n must be > 0\n",
                b->input_file, b->n, N);
        exit(EXIT_FAILURE);
    }

    for(int d = 0; d < NDIM; ++d){
        double dmin, dmax;

        if(fscanf(fp, "%lf %lf\n", &dmin, &dmax) != 2){
            fprintf(stderr, "Error: could not read box size from %s\n", b->input_file);
            exit(EXIT_FAILURE);
        }

        b->box[d] = fabs(dmax - dmin);
    }

    for(int n = 0; n < b->n; ++n){
        for(int d = 0; d < NDIM; ++d){
            if(fscanf(fp, "%lf", &b->r[n][d]) != 1){
                fprintf(stderr, "Error: could not read particle coordinates from %s\n", b->input_file);
                exit(EXIT_FAILURE);
            }
        }

        double diameter;

        if(fscanf(fp, "%lf\n", &diameter) != 1){
            fprintf(stderr, "Error: could not read particle diameter from %s\n", b->input_file);
            exit(EXIT_FAILURE);
        }
    }

    fclose(fp);
}

void read_data(void)
{
    read_box(&gas);
    read_box(&liq);

    n_tot = gas.n + liq.n;
}

void write_box(const Box* b, int step)
{
    char filename[128];
    snprintf(filename, sizeof(filename), "coords_%s_step%07d.dat", b->label, step);

    FILE* fp = fopen(filename, "w");
    if(fp == NULL){
        fprintf(stderr, "Error: could not write output file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    fprintf(fp, "%d\n", b->n);

    for(int d = 0; d < NDIM; ++d){
        fprintf(fp, "%lf %lf\n", 0.0, b->box[d]);
    }

    for(int n = 0; n < b->n; ++n){
        for(int d = 0; d < NDIM; ++d){
            fprintf(fp, "%f\t", b->r[n][d]);
        }

        fprintf(fp, "%lf\n", 1.0);
    }

    fclose(fp);
}

void write_data(int step)
{
    write_box(&gas, step);
    write_box(&liq, step);
}

void rescale_box(Box* b, double scale_factor)
{
    for(int n = 0; n < b->n; ++n){
        for(int d = 0; d < NDIM; ++d){
            b->r[n][d] *= scale_factor;
        }
    }

    for(int d = 0; d < NDIM; ++d){
        b->box[d] *= scale_factor;
    }
}

void set_overall_density(void)
{
    double current_volume = box_volume(&gas) + box_volume(&liq);
    double target_volume  = ((double)n_tot) / overall_density;
    double scale_factor   = pow(target_volume / current_volume, 1.0 / NDIM);

    rescale_box(&gas, scale_factor);
    rescale_box(&liq, scale_factor);
}

double distance_squared_pbc(const Box* b, const double a[NDIM], const double c[NDIM])
{
    double r2 = 0.0;

    for(int d = 0; d < NDIM; ++d){
        double dr = a[d] - c[d];

        dr -= b->box[d] * round(dr / b->box[d]);
        r2 += dr * dr;
    }

    return r2;
}

double particle_energy_at_position(const Box* b, const double pos[NDIM], int skip_index)
{
    double en = 0.0;
    double r_cut2 = r_cut * r_cut;

    for(int j = 0; j < b->n; ++j){
        if(j == skip_index){
            continue;
        }

        double r2 = distance_squared_pbc(b, pos, b->r[j]);

        if(r2 < r_cut2){
            double inv_r2  = 1.0 / r2;
            double inv_r6  = inv_r2 * inv_r2 * inv_r2;
            double inv_r12 = inv_r6 * inv_r6;

            en += 4.0 * (inv_r12 - inv_r6) - e_cut;
        }
    }

    return en;
}

double particle_energy(const Box* b, int i)
{
    assert(i >= 0 && i < b->n);

    return particle_energy_at_position(b, b->r[i], i);
}

double box_energy(const Box* b)
{
    double en = 0.0;

    for(int i = 0; i < b->n; ++i){
        en += particle_energy(b, i);
    }

    return 0.5 * en;
}

void check_box(const Box* b)
{
    if(b->n == 0){
        fprintf(stderr, "Error: number of particles in %s box is zero.\n", b->label);
        exit(EXIT_FAILURE);
    }

    for(int d = 0; d < NDIM; ++d){
        assert(r_cut <= 0.5 * b->box[d]);
    }
}

int displacement(void)
{
    return 0;
}

int change_volume(void)
{
    double tot_volume_i = box_volume(&gas) + box_volume(&liq);
    double tot_volume_f = exp(log(tot_volume_i)+log(dV));
    



    return 0;
}

int particle_transfer(void)
{
    return 0;
}

int main(int argc, char* argv[])
{
    assert(delta > 0.0);

    e_cut = 4.0 * (pow(1.0 / r_cut, 12.0) - pow(1.0 / r_cut, 6.0));

    read_data();
    check_box(&gas);
    check_box(&liq);

    size_t seed = time(NULL);
    dsfmt_seed(seed);

    dV = 0.05 * (box_volume(&gas) + box_volume(&liq));
    gas.energy = box_energy(&gas);
    liq.energy = box_energy(&liq);
    E_tot=gas.energy+liq.energy;

    printf("Starting gas volume:    %f\n", box_volume(&gas));
    printf("Starting liquid volume: %f\n", box_volume(&liq));
    printf("Starting total volume:  %f\n", box_volume(&gas) + box_volume(&liq));
    printf("Starting total energy:        %f\n", E_tot);
    printf("Starting seed:          %lu\n", seed);

    FILE* fp = fopen("measurements.dat", "w");
    if(fp == NULL){
        fprintf(stderr, "Error: could not write measurements.dat\n");
        exit(EXIT_FAILURE);
    }

    fclose(fp);

    return 0;
}
