#pragma once


class ATOMS{
    public:
        ATOMS();
        double mass = 1.;
        double x, y, z;
        double vx, vy, vz;
        double vxp, vyp, vzp;
        double ax, ay, az;
        double fx, fy, fz;
        double E_pot, E_kin;
};

class PARAMETRY{
    public:
        int N, Nit;
        double r_ref, rCut;
        double Lx, Ly, Lz;
        double T, T_max, T_min, T_ref, k_B, cryst_temp;
        double m_ref;
        double dt;
        double tau_P, tau_T;
};

class POTENTIAL{
    public:
        double C1, C2;
        double n1, n2;
        double k, phi;
        double r_cut;
        std::vector<double> force_table;
        std::vector<double> energy_table;
        double r_min = 0.5/2.5; 
        double dr;
        double inv_dr; // <--- DODANO ODWROTNOŚĆ
        int tab_size = 10000;

        double calcEnergy(double r);
        double calcForce(double r);
        inline double getForceFast(double r) const; // <--- Dodano const
        inline double getEnergyFast(double r) const; // <--- Dodano const
        void buildTable();
};

class CELL{
    public:
        int ncx, ncy, ncz;
        double L_cx, L_cy, L_cz; 
        std::vector<int> head;
        std::vector<int> list;

        CELL(const PARAMETRY& params);

        int get_index(const ATOMS& part);
};

void calcPotAndForce(const PARAMETRY& params, std::vector<ATOMS>& particles, const std::vector<POTENTIAL>& pots, CELL& lista);

void Verlet(PARAMETRY& params, std::vector<ATOMS>& particles, std::vector<POTENTIAL>& pots, CELL& lista);

struct POSITION{
    double x; double y; double z;
}