#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <random>
#include "classes.h"
#include <omp.h>

ATOMS::ATOMS() = default;

double POTENTIAL::calcEnergy(double r){
    if(r >= r_cut) return 0.;
    else{
        double term1 = C1/pow(r, n1) + C2/pow(r,n2)*cos(k*r+phi);
        double term2 = n1*C1/pow(r_cut, n1)*(r/r_cut - 1./n1 - 1);
        double term3 = C2/pow(r_cut,n2)*cos(k*r_cut + phi)*(r*n2/r_cut - n2 - 1);
        double term4 = C2*k/pow(r_cut, n2)*sin(k*r_cut + phi)*(r-r_cut);
        return term1 + term2 + term3 + term4;
    } 
}

double POTENTIAL::calcForce(double r){
    if(r >= r_cut) return 0.;
    else{
        double term1 = n1*C1*(-pow(r, -(n1+1)) + pow(r_cut, -(n1+1)));
        double term2 = n2*C2*(-cos(k*r+phi)/pow(r, n2+1) + cos(k*r_cut + phi)/pow(r_cut, n2+1));
        double term3 = k*C2*(-sin(k*r+phi)/pow(r, n2) + sin(k*r_cut + phi)/pow(r_cut, n2));
        return (term1 + term2 + term3)/r;
    } 
}

void POTENTIAL::buildTable() {
    force_table.resize(tab_size);
    energy_table.resize(tab_size);
    dr = (r_cut - r_min) / (tab_size - 1);
    inv_dr = 1.0 / dr; 
    for (int i = 0; i < tab_size; i++) {
        double r = r_min + i * dr;
        force_table[i]  = calcForce(r);
        energy_table[i] = calcEnergy(r);
    }
}
        
inline double POTENTIAL::getForceFast(double r) const {
    if (r >= r_cut) return 0.0;
    if (r < r_min) return force_table[0];
    double t = (r - r_min) * inv_dr; 
    int idx = (int)t;
    if (idx >= tab_size - 1) return force_table[tab_size - 1];
    double frac = t - idx;
    return force_table[idx] + frac * (force_table[idx + 1] - force_table[idx]);
}

inline double POTENTIAL::getEnergyFast(double r) const {
    if (r >= r_cut) return 0.0;
    if (r < r_min) return energy_table[0];
    double t = (r - r_min) * inv_dr; 
    int idx = (int)t;
    if (idx >= tab_size - 1) return energy_table[tab_size - 1];
    double frac = t - idx;
    return energy_table[idx] + frac * (energy_table[idx + 1] - energy_table[idx]);
}

CELL::CELL(const PARAMETRY& params) {
    ncx = floor(params.Lx / params.rCut);
    ncy = floor(params.Ly / params.rCut);
    ncz = floor(params.Lz / params.rCut);
    L_cx = params.Lx / ncx;
    L_cy = params.Ly / ncy;
    L_cz = params.Lz / ncz;
    head.assign(ncx * ncy * ncz, -1);
    list.assign(params.N, -1);
}

int CELL::get_index(const ATOMS& part) {
    int ix = int(floor(part.x / L_cx)) % ncx;
    int iy = int(floor(part.y / L_cy)) % ncy;
    int iz = int(floor(part.z / L_cz)) % ncz;
    return ix + iy * ncx + iz * (ncx * ncy);
}

double virial_xx, virial_yy, virial_zz;

void calcPotAndForce(const PARAMETRY& params, std::vector<ATOMS>& particles,
                     const std::vector<POTENTIAL>& pots, CELL& lista) { 
    
    virial_xx = 0.; virial_yy = 0.; virial_zz = 0.;

    #pragma omp parallel for
    for (int i = 0; i < params.N; i++) {
        particles[i].fx = particles[i].fy = particles[i].fz = particles[i].E_pot = 0.0;
    }

    double half_Lx = params.Lx * 0.5;
    double half_Ly = params.Ly * 0.5;
    double half_Lz = params.Lz * 0.5;
    double rCut2 = params.rCut * params.rCut;

    std::vector<int> atom_type(params.N);
    for(int i=0; i<params.N; i++) atom_type[i] = (particles[i].mass < 2.0 ? 0 : 1);

    auto apply_pair = [&](int i, int j, bool update_j) {
        if (i == j) return;
        double dx = particles[j].x - particles[i].x;
        double dy = particles[j].y - particles[i].y;
        double dz = particles[j].z - particles[i].z;
        
        dx += (dx < -half_Lx) ? params.Lx : 0.0;
        dx -= (dx >  half_Lx) ? params.Lx : 0.0;
        dy += (dy < -half_Ly) ? params.Ly : 0.0;
        dy -= (dy >  half_Ly) ? params.Ly : 0.0;
        dz += (dz < -half_Lz) ? params.Lz : 0.0;
        dz -= (dz >  half_Lz) ? params.Lz : 0.0;
        
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 >= rCut2) return;
        
        double r = sqrt(r2);
        int pot_idx = atom_type[i] + atom_type[j];
        
        double f = pots[pot_idx].getForceFast(r);
        double e = pots[pot_idx].getEnergyFast(r);
        double fx = f * dx, fy = f * dy, fz = f * dz;
        
        #pragma omp atomic
        particles[i].fx += fx;
        #pragma omp atomic
        particles[i].fy += fy;
        #pragma omp atomic
        particles[i].fz += fz;
        #pragma omp atomic
        particles[i].E_pot += 0.5 * e;
        
        double virial_factor = update_j ? 1.0 : 0.5;
        
        // 2. ADD ATOMICS FOR VIRIALS (solves the lambda capture bug)
        #pragma omp atomic
        virial_xx -= virial_factor * fx * dx;
        #pragma omp atomic
        virial_yy -= virial_factor * fy * dy;
        #pragma omp atomic
        virial_zz -= virial_factor * fz * dz;
        
        if (update_j) {
            // 3. ADD ATOMICS FOR PARTICLE J
            #pragma omp atomic
            particles[j].fx -= fx;
            #pragma omp atomic
            particles[j].fy -= fy;
            #pragma omp atomic
            particles[j].fz -= fz;
            #pragma omp atomic
            particles[j].E_pot += 0.5 * e;
        }
    };
    
    #pragma omp parallel for collapse(3)
    for (int cx = 0; cx < lista.ncx; cx++) {
        for (int cy = 0; cy < lista.ncy; cy++) {
            for (int cz = 0; cz < lista.ncz; cz++) {
                int cell_index = cx + cy * lista.ncx + cz * lista.ncx * lista.ncy;
                int i = lista.head[cell_index];
                while (i != -1) {
                    int j = lista.list[i];      
                    while (j != -1) {
                        apply_pair(i, j, true); 
                        j = lista.list[j];
                    }
                    i = lista.list[i];
                }
                for (int ddx = -1; ddx <= 1; ddx++) {
                    for (int ddy = -1; ddy <= 1; ddy++) {
                        for (int ddz = -1; ddz <= 1; ddz++) {
                            if (ddx == 0 && ddy == 0 && ddz == 0) continue;
                            if (lista.ncx <= 2 && ddx == 1) continue;
                            if (lista.ncy <= 2 && ddy == 1) continue;
                            if (lista.ncz <= 2 && ddz == 1) continue;
                            int nx = (cx + ddx + lista.ncx) % lista.ncx;
                            int ny = (cy + ddy + lista.ncy) % lista.ncy;
                            int nz = (cz + ddz + lista.ncz) % lista.ncz;
                            int neigh = nx + ny * lista.ncx + nz * lista.ncx * lista.ncy;
                            i = lista.head[cell_index];
                            while (i != -1) {
                                int j = lista.head[neigh];
                                while (j != -1) {
                                    apply_pair(i, j, false); 
                                    j = lista.list[j];
                                }
                                i = lista.list[i];
                            }
                        }
                    }
                }
            }
        }
    }
}

void Verlet(PARAMETRY& params, std::vector<ATOMS>& particles,
            std::vector<POTENTIAL>& pots, CELL& lista) {
    int wypisanie = 1;
    double Temp = params.T, Ekin = 0., Epot = 0.;
    double dt = params.dt;
    std::ofstream dane("wyniki/dane.dat");         

    // --- NPT variables (Nosé-Hoover-Parrinello-Rahman) ---
    double zeta = 0.0;            // Thermostat friction
    double v_eps = 0.0;           // Barostat strain rate
    double P_target = 0.0;        // Target pressure (0.0 for ambient/vacuum in typical MD)
    
    double tau_T2 = params.tau_T * params.tau_T;
    double tau_P2 = params.tau_P * params.tau_P;
    
    // Barostat fictitious mass (W)
    double W = 3.0 * params.N * params.k_B * params.T * tau_P2;
    
    // Initial Volume and Pressure calculation
    double V = params.Lx * params.Ly * params.Lz;
    double P_inst = 0.0; 

    // MSD
    double MSD = 0.;
    std::vector<POSITION> position_old(params.N);

    // cooling scheme 
    int beginning_eq = 1000;
    int dT = 1e-5/params.T_ref;
    

    for (int it = 1; it <= params.Nit; it++) {
        Ekin = 0.; Epot = 0.;
        
        // Calculate instantaneous pressure from the PREVIOUS step's virials and Ekin
        // (For step 1, virials are 0, which is fine)
        P_inst = (2.0 * Ekin + virial_xx + virial_yy + virial_zz) / (3.0 * V);

        // 1. Advance thermostat and barostat variables by a half-step
        double G_zeta = (Temp / params.T - 1.0) / tau_T2;
        double zeta_half = zeta + 0.5 * dt * G_zeta;

        double G_eps = 3.0 * V * (P_inst - P_target) / W - zeta_half * v_eps;
        double v_eps_half = v_eps + 0.5 * dt * G_eps;

        // 2. Scale the simulation box volume by the strain rate
        double scale_factor = exp(dt * v_eps_half);
        params.Lx *= scale_factor;
        params.Ly *= scale_factor;
        params.Lz *= scale_factor;
        V = params.Lx * params.Ly * params.Lz;

        // Update CELL sizes based on new box dimensions
        lista.L_cx = params.Lx / lista.ncx;
        lista.L_cy = params.Ly / lista.ncy;
        lista.L_cz = params.Lz / lista.ncz;

        MSD = 0.;

        // 3. First step: Update velocities (half) and positions (full)
        #pragma omp parallel for reduction(+:MSD)
        for (int i = 0; i < params.N; i++) {
            // MSD with scaling box
            position_old[i].x = particles[i].x*scale_factor;
            position_old[i].y = particles[i].y*scale_factor;
            position_old[i].z = particles[i].z*scale_factor;

            // Apply NPT friction to half-step velocities
            particles[i].vx += 0.5 * dt * (particles[i].ax - (zeta_half + v_eps_half) * particles[i].vx);
            particles[i].vy += 0.5 * dt * (particles[i].ay - (zeta_half + v_eps_half) * particles[i].vy);
            particles[i].vz += 0.5 * dt * (particles[i].az - (zeta_half + v_eps_half) * particles[i].vz);

            // Update positions (Scale coordinates with the box expansion, then add velocity step)
            particles[i].x = particles[i].x*scale_factor + dt * particles[i].vx;
            particles[i].y = particles[i].y*scale_factor + dt * particles[i].vy;
            particles[i].z = particles[i].z*scale_factor + dt * particles[i].vz;

            // Apply Periodic Boundary Conditions (using the NEW box sizes)
            if (particles[i].x < 0) particles[i].x += params.Lx;
            else if (particles[i].x >= params.Lx) particles[i].x -= params.Lx;
            if (particles[i].y < 0) particles[i].y += params.Ly;
            else if (particles[i].y >= params.Ly) particles[i].y -= params.Ly;
            if (particles[i].z < 0) particles[i].z += params.Lz;
            else if (particles[i].z >= params.Lz) particles[i].z -= params.Lz;

            // MSD calculation
            MSD += ((particles[i].x - position_old[i].x)*(particles[i].x - position_old[i].x) + (particles[i].y - position_old[i].y)*(particles[i].y - position_old[i].y) + (particles[i].z - position_old[i].z)*(particles[i].z - position_old[i].z))*((particles[i].x - position_old[i].x)*(particles[i].x - position_old[i].x) + (particles[i].y - position_old[i].y)*(particles[i].y - position_old[i].y) + (particles[i].z - position_old[i].z)*(particles[i].z - position_old[i].z));
        }

        MSD /= params.N;

        // aktualizacja listy
        std::fill(lista.head.begin(), lista.head.end(), -1);
        std::fill(lista.list.begin(), lista.list.end(), -1);
        for (int i = 0; i < params.N; i++) {
            int cell_index = lista.get_index(particles[i]);
            lista.list[i] = lista.head[cell_index];
            lista.head[cell_index] = i;
        }

        // aktualizacja sił (will calculate new virials and Epot)
        calcPotAndForce(params, particles, pots, lista);
        
        // 4. Second step: Update velocities to full step and calculate new Ekin
        double denom = 1.0 + 0.5 * dt * (zeta_half + v_eps_half);
        
        #pragma omp parallel for reduction(+:Ekin, Epot)
        for (int i = 0; i < params.N; i++) {
            double ax = particles[i].fx / particles[i].mass;
            double ay = particles[i].fy / particles[i].mass;
            double az = particles[i].fz / particles[i].mass;
            
            particles[i].ax = ax; particles[i].ay = ay; particles[i].az = az;
            
            // Advance velocities to full step algebraically
            particles[i].vx = (particles[i].vx + 0.5 * dt * ax) / denom;
            particles[i].vy = (particles[i].vy + 0.5 * dt * ay) / denom;
            particles[i].vz = (particles[i].vz + 0.5 * dt * az) / denom;
            
            Ekin += 0.5 * particles[i].mass * (particles[i].vx*particles[i].vx + 
                                               particles[i].vy*particles[i].vy + 
                                               particles[i].vz*particles[i].vz);
            Epot += particles[i].E_pot;
        }

        // Calculate new Temperature and Pressure
        Temp = 2.0 * Ekin / (3.0 * params.N * params.k_B);
        P_inst = (2.0 * Ekin + virial_xx + virial_yy + virial_zz) / (3.0 * V);

        // 5. Advance v_eps and zeta by the remaining half-step based on NEW T and P
        G_eps = 3.0 * V * (P_inst - P_target) / W - zeta_half * v_eps_half;
        v_eps = v_eps_half + 0.5 * dt * G_eps;

        G_zeta = (Temp / params.T - 1.0) / tau_T2;
        zeta = zeta_half + 0.5 * dt * G_zeta;

        // update temperature, set cooling rate
        if(it < 1000){
            continue;
        }
        else if(Temp / 100. == 0)

        //  --- wypisanie wyników ---
        if (wypisanie == 1000) {
            wypisanie = 0;
            double Etot = Ekin + Epot;
            // P_avg requires an external rolling average implementation if you use it in your code.
            // Here, we print P_inst to visualize instantaneous box pressure.
            dane << std::fixed << std::setprecision(3) << std::setw(5) << it
                 << " " << std::fixed << std::setprecision(5) << std::setw(6)
                 << Ekin << " " << Epot << " " << Etot << " " << Temp*params.T_ref 
                 << " " << P_inst << " " << V << std::endl;
        }
        wypisanie++;
    }
}