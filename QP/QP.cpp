// ConsoleApplication1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//


#include <clocale>

#include <iomanip>
#include <iostream>
#include <vector>
#include <fstream>
#include <time.h>
#include <algorithm>

#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>


using namespace std;

struct pipe
{
    double L, // протяженность участка
        V, //скорость нефти
        D, //диаметр внутр.
        D_vnesh, //диаметр внешний
        b, //толщина стенки
        abc, //абсолютная шероховатость в м
        z_0, z_L, //высотные отметки начала,конца 
        ro, //плотность
        u, // кинематическая вязкость Стоксы в м2/с
        Q, //расход м3/ч
        p_0,// давление в начале участка
        lambda,//коэфф.гидравл.сопротивления
        n,//кол-во шагов стеки
        t_w,//касательное напряжение трения
        h; //шаг по координате расчетной сетки, в метрах
      
};


double pressure(struct pipe myPipe) {
    ofstream outFile("pressure2.csv");
    double p_0 = myPipe.p_0;
    double p_L;
    outFile << p_0 << "\n";
    for (int i = 0; i < myPipe.n; ++i) {
        p_L = p_0 + myPipe.h * (-4 / myPipe.D * myPipe.t_w - myPipe.ro * M_G * (myPipe.z_L - myPipe.z_0) / ((myPipe.n - 1) * myPipe.h));
        outFile << p_L << "\n";
        p_0 = p_L;
    }
    outFile.close();
    return p_L;
}



int main()
{
    double  Re;

    setlocale(LC_ALL, "Russian");
    //Данные по трубопроводу в СИ
    pipe myPipe;
    myPipe.p_0 = 5.65e6;
    myPipe.L = 80e3;
    myPipe.D_vnesh = 720e-3;
    myPipe.b = 10e-3;
    myPipe.z_0 = 100;
    myPipe.z_L = 50;
    myPipe.ro = 870;
    myPipe.u = 15e-6;
    myPipe.Q = 3500;
    myPipe.abc = 15e-6;
    myPipe.n = 100;
    myPipe.t_w = 0;
    myPipe.h = myPipe.L / myPipe.n;
    myPipe.D = myPipe.D_vnesh - 2 * myPipe.b;
    double Q = myPipe.Q / 3600;                     //перевод в м3/с
    myPipe.V = 4 * Q / (3.14 * pow(myPipe.D, 2));      //вычисляем скорость нефти 
    Re = myPipe.V * myPipe.D / myPipe.u;
    double relative_roughness = myPipe.abc / myPipe.D;    //вычисляем относ.шероховатость
    myPipe.lambda = hydraulic_resistance_isaev(Re, relative_roughness);
    myPipe.t_w = myPipe.lambda / 8 * myPipe.ro * pow(myPipe.V,2);

    double p_L = pressure(myPipe);   //присваиваем значения из функции pressure

}