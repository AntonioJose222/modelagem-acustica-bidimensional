#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>
#include "Matriz3D.cpp"

#define M_PI 3.14159265358979323846 
#define M_e 2.71828182845904523536

using namespace std;

int main() {

    float X = 3000; //Largura do dominio em m
    float Z = 3000; //Altura do dominio em m
    float T = 1;    //Tempo total em s
    float c = 2200;  //Celeridade da onda em m/s
    float fp = 40;  //frequencia de pico = 80Hz

    float alpha = 5;
    float beta = 5;

    float dx = 10;  //Passos em x
    float dz = dx;  //Passos em z
    float dt = 0.00025;   //Passos no tempo

    int Nx = 300; //Iteracoes em x
    int Nz = 300; //Iteracoes em z
    int Nt = T/dt+1; //Iteracoes no tempo

    int xs = 150; //posicao da fonte em x
    int zs = 150; //posicao da fonte em z
    float cou = c*dt/dx; //numero de courant

    cout << "Nx = " << Nx << endl;
    cout << "Nz = " << Nz << endl;
    cout << "Nt = " << Nt << endl;
    cout << "Numero de Courant = " << cou << endl;    

    float val;

    float fonte(int x, int z, float t, float fp, float xs, float zs);
    float cerjan(int x, int z, int Nx, int Nz, float P0);//funcao usada para impedir reflexoes nas bordas

    bool yet = false;

    Matriz3D u(Nx, Nz, Nt);

    ofstream myfile;
    
    u.set(xs, zs, 0, fonte(xs, zs, 0, fp, xs, zs));
    
    for (int k = 1; k < Nt - 1; k++){

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
                + 2*u.get(i, j, k) - u.get(i, j, k - 1) - 
                pow(c*dt, 2) * fonte(i, j, (float)k*dt, fp, xs, zs);

                val = cerjan(i, j, Nx, Nz, val);

                u.set(i, j, k + 1, val);
            
            }

        }
        if(u.get(150, 10, k) != 0 && !yet){
            yet = true;
            cout << "Onda atinge o contorno em k = " << k << " valendo: "<< u.get(150, 0, k) << endl;
        }
    }
    
    string base(".dat");
    for(int k = 0; k < Nt; k += 100){
        myfile.open("/home/antonio/Desktop/Modelagem Acustica Bidimensional/data/data" + to_string(k/100) + base);
        for (int j = 0; j < Nz; j++){
            for (int i = 0; i < Nx; i++){
                myfile << i << " " << j << " " 
                << std::setprecision(17) << u.get(i, j, k) << "\n";
            }
            myfile << "\n\n";
        }
        myfile.close();
    }

}

float fonte(int x, int z, float t, float fp, float xs, float zs){

    if (x != xs || z != zs){
        return 0;
    }    
    return (1.0 - 2*pow(M_PI*fp*t, 2))/pow(M_e, pow((M_PI*fp*t), 2));
}

float cerjan(int x, int z, int Nx, int Nz, float P0){
    int n = 20;
    if (x < n){
        return P0 / pow(M_e, pow(0.98*(n - x), 2));
    } else if (x > Nx - n){
        return P0 / pow(M_e, pow(0.98*(n - (Nx - x)), 2));
    } else if (z < n){
        return P0 / pow(M_e, pow(0.98*(n - z), 2));
    } else if (z > Nz - n){
        return P0 / pow(M_e, pow(0.98*(n - (Nz - z)), 2));
    } else {
        return P0;
    }

}