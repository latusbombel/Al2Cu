#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <random>
#include "classes.h"
#include <omp.h>

void minimalizacja(PARAMETRY& params, std::vector<ATOMS>& particles, std::vector<POTENTIAL>& pots, CELL& lista){
    int wypisanie = 0;
    std::ofstream minimalizacja("wyniki/minimalizacja.dat");
    double Epot_old = 0.;
    double Ekin = 0.;
    double Epot = 0.;

    #pragma omp parallel for
    for (int it = 1; it <= 100; it++) {
        double dt = params.dt;
        Epot_old = Epot;
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
        
        #pragma omp parallel for reduction(+:Ekin, Epot)
        for (int i = 0; i < params.N; i++) {
            double ax = particles[i].fx / particles[i].mass;
            double ay = particles[i].fy / particles[i].mass;
            double az = particles[i].fz / particles[i].mass;
            
            particles[i].ax = ax; particles[i].ay = ay; particles[i].az = az;
            
            // Advance velocities to full step algebraically
            particles[i].vx = (particles[i].vx + 0.5 * dt * ax);
            particles[i].vy = (particles[i].vy + 0.5 * dt * ay);
            particles[i].vz = (particles[i].vz + 0.5 * dt * az);
            
            Ekin += 0.5 * particles[i].mass * (particles[i].vx*particles[i].vx + 
                                               particles[i].vy*particles[i].vy + 
                                               particles[i].vz*particles[i].vz);
            Epot += particles[i].E_pot;
        }

        //  --- wypisanie wyników ---
        if (wypisanie == 10) {
            wypisanie = 0;
            double Etot = Ekin + Epot;
            
            minimalizacja << std::fixed << std::setprecision(3) << std::setw(5) << it
                 << " " << std::fixed << std::setprecision(8) << std::setw(8)
                 << Ekin << " " << Epot << " " << Etot << std::endl;
        }
        wypisanie++;

        if(abs(Epot - Epot_old) < 1e-6) break;
    }
}