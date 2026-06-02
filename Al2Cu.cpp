#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <random>
#include "classes.h"
#include <chrono>
#include <algorithm>


int main(){
    // inicjalizacja parametrów symulacji
    PARAMETRY params;
    params.N = 12000;
    params.r_ref = 2.5; params.rCut = 12./params.r_ref; // w Angstremach
    params.Lx = 58.49/params.r_ref;
    params.Ly = 58.49/params.r_ref;
    params.Lz = 58.49/params.r_ref;
    params.k_B = 1.;
    params.Nit = 10000;
    params.dt = 0.01; 
    params.tau_P = 1000*params.dt;
    params.tau_T = 100*params.dt;

    // inicjalizacja potencjału
    // Al-Al
    POTENTIAL potAl;
    // potAl.C1 = 1.;
    potAl.n1 = 8.3309; potAl.n2 = 4.3472;
    potAl.phi = -10.6377;
    double C1 = 771.62;
    double C2 = -2.4184;
    double k = 3.8752;
    double E_ref = C1/pow(params.r_ref, potAl.n1);
    potAl.C1 = 771.62/(E_ref*pow(params.r_ref, potAl.n1));
    potAl.C2 = C2/(E_ref*pow(params.r_ref, potAl.n2));
    potAl.k = k*params.r_ref;
    potAl.r_cut = params.rCut;
    potAl.r_min = 1./params.r_ref;

    // Cu-Cu
    POTENTIAL potCu;
    potCu.n1 = 7.8303; potCu.n2 = 2.7515;
    potCu.C1 = 304.2326/(E_ref*pow(params.r_ref, potCu.n1));
    potCu.C2 = -0.9824/(E_ref*pow(params.r_ref, potCu.n2));
    potCu.phi = -6.9253;
    potCu.k = 3.1140*params.r_ref;
    potCu.r_cut = params.rCut;
    potCu.r_min = 1./params.r_ref;

    // Al-Cu
    POTENTIAL potAlCu;
    potAlCu.n1 = 7.1380; potAlCu.n2 = 3.5251;
    potAlCu.C1 = 164.5061/(E_ref*pow(params.r_ref, potAlCu.n1));
    potAlCu.C2 = -2.6665/(E_ref*pow(params.r_ref, potAlCu.n2));
    potAlCu.phi = -7.8115;
    potAlCu.k = 3.1879*params.r_ref;
    potAlCu.r_cut = params.rCut;
    potAlCu.r_min = 1./params.r_ref;

    params.T_ref = E_ref/(8.617*1e-5);
    params.T = 1000./params.T_ref;
    params.T_max = params.T;
    params.T_min = 100./params.T_ref;
    params.cryst_temp = 500/params.T_ref;

    // Inicjalizacja atomów
    ATOMS Al;
    Al.mass = 1.;
    
    ATOMS Cu;
    Cu.mass = 2.3549;

    std::vector<ATOMS> particles(params.N);
    std::default_random_engine generator;
    std::normal_distribution normal(0.0, 1.0);

    // Wczytywanie idealnej sieci krystalicznej z pliku
    // std::ifstream wczytaj("start_pos.dat");
    // if (!wczytaj.is_open()) {
    //     std::cerr << "Blad: Nie mozna otworzyc pliku start_pos.dat!" << std::endl;
    //     return 1;
    // }

    // for(int i = 0; i<params.N; i++){
    //     // Ponieważ skrypt w Pythonie zapisał najpierw wszystkie Al, a potem Cu:
    //     if(i < (params.N * 2 / 3)) particles[i] = Al; // 2/3 to aluminium
    //     else particles[i] = Cu;                       // 1/3 to miedź
        
    //     // Wczytanie dokładnych współrzędnych XYZ
    //     wczytaj >> particles[i].x >> particles[i].y >> particles[i].z;
    //     particles[i].x /= params.r_ref;
    //     particles[i].y /= params.r_ref;
    //     particles[i].z /= params.r_ref;
    // }
    // wczytaj.close();

    //losowe zamienienie atomów

    // int n_Al = params.N * 2 / 3;
    // for (int i = 0; i < params.N; i++) {
    //     if (type_indices[i] < n_Al) particles[i] = Al;
    //     else particles[i] = Cu;
    //     // positions are already loaded, copy them back:
    //     wczytaj >> particles[i].x >> particles[i].y >> particles[i].z;
    //     particles[i].x /= params.r_ref;
    //     particles[i].y /= params.r_ref;
    //     particles[i].z /= params.r_ref;
    // }

    std::vector<int> type_indices(params.N);
    for (int i = 0; i < params.N; i++) type_indices[i] = i;
    std::shuffle(type_indices.begin(), type_indices.end(), generator);
    int n_Al = params.N * 2 / 3;
    double delta_xy = params.Lx/ceil(cbrt(params.N));
    int nodes_xy = floor(params.Lx/delta_xy);
    double delta_z = params.Lz/ceil(cbrt(params.N));
    int nodes_z = floor(params.Lz/delta_z);
    double offset_xy = delta_xy/7.;
    double offset_z = delta_z/7.;
    for(int i = 0; i<params.N; i++){
        if (type_indices[i] < n_Al) particles[i] = Al;
        else particles[i] = Cu;
        double rand_x = normal(generator);
        double rand_y = normal(generator);
        double rand_z = normal(generator);
        particles[i].x = delta_xy/2. + (i%nodes_xy)*delta_xy + offset_xy*rand_x;
        particles[i].y = delta_xy/2. + (int(floor(i/nodes_xy))%nodes_xy)*delta_xy + offset_xy*rand_y;
        particles[i].z = delta_z/2. + floor(i/(nodes_z*nodes_z))*delta_z + offset_z*rand_z;
    }

    // Nadanie prędkości początkowych

    for(int i = 0; i<params.N; i++){
        double sigma = sqrt(params.k_B*params.T/particles[i].mass);
        double sx = normal(generator);
        double sy = normal(generator);
        double sz = normal(generator);

        particles[i].vx = sigma * sx;
        particles[i].vy = sigma * sy;
        particles[i].vz = sigma * sz;
    }

    // usunięcie dryfu prędkości środka masy
    double vx_cm = 0., vy_cm = 0., vz_cm = 0.;

    double M = 0;
    for(int i = 0; i<params.N; i++){
        vx_cm += particles[i].vx*particles[i].mass;
        vy_cm += particles[i].vy*particles[i].mass;
        vz_cm += particles[i].vz*particles[i].mass;
        M += particles[i].mass;
    }

    vx_cm /= M;
    vy_cm /= M;
    vz_cm /= M;


    for(int i = 0; i<params.N; i++){
        particles[i].vx -= vx_cm;
        particles[i].vy -= vy_cm;
        particles[i].vz -= vz_cm;
    }

    std::ofstream polozenie_pocz("wyniki/polozenie_pocz.dat");
    int rodzaj;
    for(int i = 0; i<params.N; i++){
        if(particles[i].mass > 2) rodzaj = 1;
        else rodzaj = 0;
        polozenie_pocz << i << " " << rodzaj << " " << 
        std::fixed << std::setprecision(5) << std::setw(6) << 
        particles[i].x << " " << particles[i].y << " " << particles[i].z << " " <<
        particles[i].vx << " " << particles[i].vy << " " << particles[i].vz << std::endl;
    }

    //inicjalizacja połączonej listy

    CELL lista(params);

    for(int i = 0; i<params.N; i++){
        int cell_index = lista.get_index(particles[i]);
        lista.list[i] = lista.head[cell_index];
        lista.head[cell_index] = i;
    }

    std::vector<POTENTIAL> pots = {potAl, potAlCu, potCu};

    auto start = std::chrono::high_resolution_clock::now();

    pots[0].buildTable();
    pots[1].buildTable();
    pots[2].buildTable();

    calcPotAndForce(params, particles, pots, lista);

    for(int i = 0; i < params.N; i++){
    particles[i].ax = particles[i].fx / particles[i].mass;
    particles[i].ay = particles[i].fy / particles[i].mass;
    particles[i].az = particles[i].fz / particles[i].mass;
    }

    Verlet(params, particles, pots, lista);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Execution time: " << elapsed.count() << " s" << std::endl;

    std::ofstream polozenie_end("wyniki/polozenie_end.dat");
    for(int i = 0; i<params.N; i++){
        if(particles[i].mass > 2) rodzaj = 1;
        else rodzaj = 0;
        polozenie_end << i << " " << rodzaj << " " << 
        std::fixed << std::setprecision(5) << std::setw(6) << 
        particles[i].x << " " << particles[i].y << " " << particles[i].z << " " <<
        particles[i].vx << " " << particles[i].vy << " " << particles[i].vz << std::endl;
    }

    return 0;
}