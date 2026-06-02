#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <omp.h>

struct Atom{
    double x,y,z;
    int type; // 0 = Al, 1 = Cu
};

double formFactor(double q,
                  const double a[4],
                  const double b[4],
                  double c)
{
    double s = q/(4.0*M_PI);
    double s2 = s*s;

    double f = c;

    for(int i=0;i<4;i++)
        f += a[i]*exp(-b[i]*s2);

    return f;
}

int main(){

    std::ifstream file("wyniki_male/polozenie_pocz.dat");

    std::vector<Atom> atoms;

    int id;
    double x,y,z,vx,vy,vz;

    while(file >> id >> x >> y >> z >> vx >> vy >> vz){

        Atom a;
        a.x = x;
        a.y = y;
        a.z = z;

        // same logic as MD code
        if(id < 1728) a.type = 0;  // Al
        else a.type = 1;           // Cu

        atoms.push_back(a);
    }

    int N = atoms.size();

    std::cout<<"Atoms loaded: "<<N<<"\n";

    double r_ref = 2.5;
    // double L = 36.402/r_ref;
    double Lx = 36.402 / (r_ref*3);
    double Ly = 36.402 / (r_ref*3);
    double Lz = 29.244 / (r_ref*2);

    // ---------- pair distances ----------

    struct Pair{
        double r;
        int t1;
        int t2;
    };

    std::vector<Pair> pairs;

    pairs.reserve((size_t)N*(N-1)/2);

#pragma omp parallel
{
    std::vector<Pair> local_pairs;

#pragma omp for schedule(dynamic)
    for(int i=0;i<N;i++){
        for(int j=i+1;j<N;j++){

            double dx = atoms[j].x - atoms[i].x;
            double dy = atoms[j].y - atoms[i].y;
            double dz = atoms[j].z - atoms[i].z;

            // dx -= Lx*round(dx/Lx);
            // dy -= Ly*round(dy/Ly);
            // dz -= Lz*round(dz/Lz);

            double r = sqrt(dx*dx + dy*dy + dz*dz);

            local_pairs.push_back({r,atoms[i].type,atoms[j].type});
        }
    }

#pragma omp critical
    pairs.insert(pairs.end(), local_pairs.begin(), local_pairs.end());
}

    std::cout<<"Pairs: "<<pairs.size()<<"\n";

    // ---------- form factor parameters ----------

    double aAl[4] = {6.4202,1.9002,1.5936,1.9646};
    double bAl[4] = {3.0387,0.7426,31.5472,85.0886};
    double cAl = 1.1151;

    double aCu[4] = {13.3380,7.1676,5.6158,1.6735};
    double bCu[4] = {3.5828,0.2470,11.3966,64.8126};
    double cCu = 1.1910;

    // ---------- XRD parameters ----------

    double wavelength = 1.540598 ;

    int Nq = 600;

    double theta_min = 8.0;
    double theta_max = 100.0;

    std::vector<double> q(Nq);
    std::vector<double> intensity(Nq);

    for(int i=0;i<Nq;i++){

        double two_theta = theta_min + (theta_max-theta_min)*i/(Nq-1);

        double theta = two_theta*M_PI/360.0;

        q[i] = 4.0*M_PI*sin(theta)/wavelength;
    }

    // ---------- Debye sum ----------

#pragma omp parallel for
    for(int iq=0;iq<Nq;iq++){

        double qv = q[iq];

        double fAl = formFactor(qv,aAl,bAl,cAl);
        double fCu = formFactor(qv,aCu,bCu,cCu);

        double sum = 0.0;

        // self scattering
        int NAl = 2/3*N;
        int NCu = N - NAl;

        sum += NAl*fAl*fAl + NCu*fCu*fCu;

        for(size_t p=0;p<pairs.size();p++){

            double qr = qv*pairs[p].r * r_ref;

            double sinc = (qr>1e-10) ? sin(qr)/qr : 1.0;

            double fi = (pairs[p].t1==0) ? fAl : fCu;
            double fj = (pairs[p].t2==0) ? fAl : fCu;

            sum += 2.0*fi*fj*sinc;
        }

        intensity[iq] = sum;
    }

    // ---------- output ----------

    std::ofstream out("wyniki_male/xrd.dat");

    for(int i=0;i<Nq;i++){

        double two_theta = theta_min + (theta_max-theta_min)*i/(Nq-1);

        out<<two_theta<<" "<<intensity[i]<<"\n";
    }

    double fAl0 = formFactor(0.0, aAl, bAl, cAl);
    double fCu0 = formFactor(0.0, aCu, bCu, cCu);
    std::cout << "fAl(0) = " << fAl0 << ", fCu(0) = " << fCu0 << std::endl;

    std::cout<<"XRD calculation finished\n";
}

// #include <iostream>
// #include <vector>
// #include <fstream>
// #include <cmath>
// #include <omp.h>

// struct Atom {
//     double x, y, z;
//     int type;
// };

// double formFactor(double q, const double a[4], const double b[4], double c) {
//     double s = q / (4.0 * M_PI);
//     double s2 = s * s;
//     double f = c;
//     for (int i = 0; i < 4; i++) f += a[i] * exp(-b[i] * s2);
//     return f;
// }

// int main() {
//     std::ifstream file("wyniki/polozenie_end.dat");
//     std::vector<Atom> atoms;
//     int id; double x, y, z, vx, vy, vz;
    
//     while (file >> id >> x >> y >> z >> vx >> vy >> vz) {
//         atoms.push_back({x, y, z, (id < 1728) ? 0 : 1});
//     }
//     std::cout << "Atoms loaded: " << atoms.size() << "\n";

//     double r_ref = 2.5;
//     double Lx = 36.402/r_ref; // w angstremach!
//     double Ly = 36.402/r_ref;
//     double Lz = 29.244/r_ref;

//     double aAl[4] = {6.4202, 1.9002, 1.5936, 1.9646};
//     double bAl[4] = {3.0387, 0.7426, 31.5472, 85.0886};
//     double cAl = 1.1151;

//     double aCu[4] = {13.3380, 7.1676, 5.6158, 1.6735};
//     double bCu[4] = {3.5828, 0.2470, 11.3966, 64.8126};
//     double cCu = 1.1910;

//     double wavelength = 1.540598;
//     int Nq = 2000;
//     double theta_min = 8.0, theta_max = 100.0;
//     std::vector<double> intensity(Nq, 0.0);

//     int h_max = 10, k_max = 10, l_max = 10;
//     std::cout << "Calculating 3D Structure Factor...\n";

//     #pragma omp parallel for collapse(3)
//     for (int h = -h_max; h <= h_max; h++) {
//         for (int k = -k_max; k <= k_max; k++) {
//             for (int l = -l_max; l <= l_max; l++) {
//                 if (h == 0 && k == 0 && l == 0) continue;

//                 double qx = 2.0 * M_PI * h / Lx;
//                 double qy = 2.0 * M_PI * k / Ly;
//                 double qz = 2.0 * M_PI * l / Lz;
//                 double q = sqrt(qx*qx + qy*qy + qz*qz);

//                 double sin_theta = q * wavelength / (4.0 * M_PI);
//                 if (sin_theta >= 1.0) continue;

//                 double two_theta = 2.0 * asin(sin_theta) * 180.0 / M_PI;
//                 if (two_theta < theta_min || two_theta > theta_max) continue;

//                 double sum_re = 0.0, sum_im = 0.0;
//                 double fAl = formFactor(q, aAl, bAl, cAl);
//                 double fCu = formFactor(q, aCu, bCu, cCu);

//                 for (const auto& a : atoms) {
//                     double phase = qx*(a.x * r_ref) + qy*(a.y * r_ref) + qz*(a.z * r_ref);
//                     double f = (a.type == 0) ? fAl : fCu;
//                     sum_re += f * cos(phase);
//                     sum_im += f * sin(phase);
//                 }
                
//                 double I = sum_re*sum_re + sum_im*sum_im;

//                 // Sztuczne nałożenie drobnego rozmycia Gaussa, żeby wąskie piki były widoczne na wykresie
//                 double sigma = 0.15; 
//                 for (int i = 0; i < Nq; i++) {
//                     double th = theta_min + (theta_max - theta_min) * i / (Nq - 1);
//                     double diff = th - two_theta;
//                     if (std::abs(diff) < 2.0) {
//                         #pragma omp atomic
//                         intensity[i] += I * exp(-(diff*diff)/(2.0*sigma*sigma));
//                     }
//                 }
//             }
//         }
//     }

//     std::ofstream out("wyniki/xrd.dat");
//     for (int i = 0; i < Nq; i++) {
//         double two_theta = theta_min + (theta_max - theta_min) * i / (Nq - 1);
//         out << two_theta << " " << intensity[i] << "\n";
//     }
//     std::cout << "Done!\n";
// }