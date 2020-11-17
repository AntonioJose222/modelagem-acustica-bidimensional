#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "Matriz3D.cpp"

#define M_PI 3.14159265358979323846 
#define M_e 2.71828182845904523536

using namespace std;

int main() {

    cout.precision(17);

    float X = 1100; //Largura do dominio
    float Z = 600;  //Altura do dominio
    float T = 2000;  //Tempo total
    float c = 1.4;  //Celeridade da onda em km/s
    float fp = 80;  //frequencia de pico = 80Hz

    float alpha = 5;
    float beta = 5;

    float dx = c/fp*alpha;  //Passos em x
    float dz = dx;          //Passos em z
    float dt = dx/c*beta;   //Passos no tempo

    int Nx = X/dx + 1; //Iteracoes em x
    int Nz = Z/dz + 1; //Iteracoes em z
    int Nt = T/dt + 1; //Iteracoes no tempo

    int xs = 550; //posicao da fonte em x
    int zs = 300; //posicao da fonte em z
    float cou = c*dt/dx;

    cout << "Nx = " << Nx << endl;
    cout << "Nz = " << Nz << endl;
    cout << "Nt = " << Nt << endl;
    cout << "Numero de Courant = " << cou << endl;    

    float val;

    float fonte(int x, int z, int t, float fp, float xs, float zs);
    float cerjan(int x, int z, int Nx, int Nz, float P0);//funcao usada para impedir reflexoes nas bordas

    Matriz3D u(Nx, Nz, Nt);

    ofstream myfile;

    u.set(xs/dx, zs/dz, 0, fonte(xs/dx, zs/dz, 0, fp, xs/dx, zs/dz));
    
    for (int k = 0; k < Nt - 1; k++){

        for (int j = 0; j < Nz; j++){

            for (int i = 0; i < Nx; i++){

                val = 
                (1.0/12.0)*pow((cou), 2) * 
                (
                    -1*(u.get(i - 2, j, k) + u.get(i, j - 2, k)) + 
                    16*(u.get(i - 1, j, k) + u.get(i, j - 1, k)) - 60 * u.get(i, j, k) +
                    16*(u.get(i + 1, j, k) + u.get(i, j + 1, k)) -
                    (u.get(i + 2, j, k) + u.get(i, j + 2, k)) 
                )
                + 2*u.get(i, j, k) - u.get(i, j, k - 1) - pow(c*dt, 2) * fonte(i, j, k, fp, xs, zs);
                val = cerjan(i, j, Nx, Nz, val);
                u.set(i, j, k + 1, val);
            
            }

        }

    }

    myfile.open("data.dat");
    
    for (int j = 0; j < Nz; j++){
        for (int i = 0; i < Nx; i++){
            myfile << i*dx << " " << j*dz << " " << std::setprecision(17) << u.get(i, j, 10) << "\n";
        }
        myfile << "\n\n";
    }

    myfile.close();
    
}

float fonte(int x, int z, int t, float fp, float xs, float zs){
    if (x != xs || z != zs){
        return 0;
    }
    return (1.0 - 2*pow(M_PI*fp*t, 2))*pow(M_e, pow((-1*M_PI*fp*t), 2));
}

float cerjan(int x, int z, int Nx, int Nz, float P0){
     
    if (x <= 40){
        return P0 * pow(M_e, -1*pow(0.015*(40 - x), 2));
    } else if (x >= Nx - 20){
        return P0 * pow(M_e, -1*pow(0.015*(20 - (Nx - x)), 2));
    } else if (z <= 20){
        return P0 * pow(M_e, -1*pow(0.015*(20 - z), 2));
    } else if (z >= Nz - 20){
        return P0 * pow(M_e, -1*pow(0.015*(20 - (Nz - z)), 2));
    } else {
        return P0;
    }

}