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
#include <fixed/fixed_nonlinear_solver.h>
#include <pde_solvers/pde_solvers.h>


using namespace std;


//по Ньютону

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
        p_L,// давление в начале участка
        lambda,//коэфф.гидравл.сопротивления
        n,//кол-во шагов стеки
        t_w,//касательное напряжение трения
        h; //шаг по координате расчетной сетки, в метрах
    double relative_roughness, Re;
};
/// инициализация функций для вызова в main 
 
void u_Newton(pipe& myPipe);
void u_Newton_Euler (pipe& myPipe);



int main()
{
    
    setlocale(LC_ALL, "Russian");
    //Данные по трубопроводу в СИ
    pipe myPipe;
    myPipe.p_0 = 5e6;
    myPipe.p_L = 0.8e6;
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
    myPipe.Re = myPipe.V * myPipe.D / myPipe.u;
    myPipe.relative_roughness = myPipe.abc / myPipe.D;    //вычисляем относ.шероховатость
    myPipe.lambda = hydraulic_resistance_isaev(myPipe.Re, myPipe.relative_roughness);
    myPipe.t_w = myPipe.lambda / 8 * myPipe.ro * pow(myPipe.V, 2);

    //double p_L = pressure_EULER(myPipe);   //присваиваем значения из функции pressure
   
    u_Newton(myPipe);
    u_Newton_Euler(myPipe);
    
};





class QP_Newton : public fixed_system_t<1>
{
    /// @brief Ссылка на структуру с параметрами трубы 
    const  pipe& pipe_dannye; //pipe моя структура, pipe_dannye -ссылка на неё

using fixed_system_t<1>::var_type; public:
    /// @brief Констуктор
    /// @param pipe Ссылка на сущность трубы 
    /// @param problem Ссылка на сущность с условием задачи
    //  конструктор QP_Newton принимает  аргумент - ссылку на объект класса pipe. 
    QP_Newton(const pipe& pipe_dannye) : pipe_dannye{ pipe_dannye } {};
    /// @brief Функция невязок - все члены уравнения Бернулли в правой части 
    /// @param v - скорость течения нефти, [м/с] 
    /// @return Значение функции невязок при заданной скорости (из библиотеки fixed solvers )
    var_type residuals(const var_type& v)  
    { //v - искомая скорость
        double Re = v * pipe_dannye.D / pipe_dannye.u;

        double lambda = hydraulic_resistance_isaev(pipe_dannye.Re, pipe_dannye.relative_roughness);
        double p_L = (pipe_dannye.p_0 / (pipe_dannye.ro * M_G) - pipe_dannye.z_0 + pipe_dannye.z_L - pipe_dannye.lambda * (pipe_dannye.L / pipe_dannye.D) * pow(v, 2) / (2 * M_G)) * (pipe_dannye.ro * M_G);
        double delta_p_L;


        return //возвращает разницу между заданным давлением в конце и рассчитанным
        {

           delta_p_L = pipe_dannye.p_L - p_L

        };
    };
};

void u_Newton(struct pipe& myPipe) {

    // Подключение бибилиотеки
// Класс, для системы размерности <2> - Векторный случай
// <2> - Размерность системы уравнений


    // Создание экземпляра класса, который и будет решаемой системой
    QP_Newton test(myPipe);
    // Задание настроек решателя по умолчанию
    fixed_solver_parameters_t<1, 0> parameters;
    // Создание структуры для записи результатов расчета
    fixed_solver_result_t<1> result;
    // Решение системы нелинейныйх уравнений <1> с помощью решателя Ньютона - Рафсона
    // { 0, 0 } - Начальное приближение
    fixed_newton_raphson<1>::solve_dense(test, { 1 }, parameters, &result);
    cout << '\n' << "Скорость изначальная " << '\n' << "u = " << myPipe.V << '\n' << "Классическая задача PP поверх метода Эйлера " << '\n' << "u = " << result.argument;
}





class QP_Newton_Euler : public fixed_system_t<1>
{
    /// @brief Ссылка на структуру с параметрами трубы 
    const  pipe& pipe_dannye;
    
using fixed_system_t<1>::var_type; public:
    /// @brief Констуктор
    /// @param pipe Ссылка на сущность трубы 
    /// @param problem Ссылка на сущность с условием задачи
    QP_Newton_Euler(const pipe& pipe_dannye) : pipe_dannye{ pipe_dannye }{};

    /// @brief Функция невязок - все члены уравнения Бернулли в правой части 
    /// @param v - скорость течения нефти, [м/с] 
    /// @return Значение функции невязок при заданной скорости 
    var_type residuals(const var_type& v)
    { //v - искомая скорость
        double Re = v * pipe_dannye.D / pipe_dannye.u;
        double lambda = hydraulic_resistance_isaev(pipe_dannye.Re, pipe_dannye.relative_roughness);
        double t_w = lambda / 8 * pipe_dannye.ro * pow(v, 2);
        double p_0_rachet;
        double p_L= pipe_dannye.p_L;
        //последовательный расчет p_0 по методу Эйлера 
        for (int i = 0; i < pipe_dannye.n; ++i) {
            p_0_rachet = p_L - pipe_dannye.h * (-4 / pipe_dannye.D * t_w - pipe_dannye.ro * M_G * (pipe_dannye.z_L - pipe_dannye.z_0) / ((pipe_dannye.n - 1) * pipe_dannye.h));
            p_L = p_0_rachet;
        }
        
                
         double delta_p_0;
               
            return
        {

           delta_p_0 =  pipe_dannye.p_0- p_0_rachet

        };

    };
};




/// @brief 
/// @param myPipe 
void u_Newton_Euler(struct pipe& myPipe) {

    
    //В качестве аргумента для конструктора передается myPipe
    // Создание экземпляра класса, который и будет решаемой системой
    QP_Newton_Euler test(myPipe);
    // Задание настроек решателя по умолчанию
    fixed_solver_parameters_t<1, 0> parameters;
    // Создание структуры для записи результатов расчета
    fixed_solver_result_t<1> result;
    // Решение системы нелинейныйх уравнений <1> с помощью решателя Ньютона - Рафсона
    // { 1} - Начальное приближение скорости 1 (не 0 - т.к. иначе получается деление на 0)
    fixed_newton_raphson<1>::solve_dense(test, {1}, parameters, &result);
    cout << '\n' << "Классическая задача PP поверх Эйлера на методе Ньютона " << '\n' << "u = " << result.argument;
}

















