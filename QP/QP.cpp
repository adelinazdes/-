
#include <iostream>
#include <clocale>
#include <iomanip>
#include <iostream>
#include <vector>
#include <fstream>
#include <time.h>
#include <algorithm>
#include <cmath>
#include <locale.h>
#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>
#include <fixed/fixed_nonlinear_solver.h>

#include <gtest/gtest.h>  


using namespace std;
using namespace pde_solvers;


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
		Q, //расход м3/c
		p_0,// давление в начале участка
		p_L,// давление в начале участка
		lambda,//коэфф.гидравл.сопротивления
		n,//кол-во шагов стеки
		t_w,//касательное напряжение трения
		h; //шаг по координате расчетной сетки, в метрах
	double relative_roughness, Re;
};







///Освоение классических гидравлических задач QP, PP

//задача QP Лурье 1
double pressure() {
	setlocale(LC_ALL, "Russian");
	//Данные по трубопроводу в СИ
	pipe myPipe;
	myPipe.p_L = 0.6e6;
	myPipe.L = 80e3;
	myPipe.D_vnesh = 720e-3;
	myPipe.b = 10e-3;
	myPipe.z_0 = 50;
	myPipe.z_L = 100;
	myPipe.ro = 870;
	myPipe.u = 15e-6;
	double Q = 3500;
	myPipe.Q = Q / 3600;
	myPipe.abc = 15e-6;//шероховатость

	myPipe.D = myPipe.D_vnesh - 2 * myPipe.b;
	
	myPipe.V = 4 * myPipe.Q / (3.1415 * pow(myPipe.D, 2));      //вычисляем скорость нефти 

	double Re = myPipe.V * myPipe.D / myPipe.u;
	double relative_roughness = myPipe.abc / myPipe.D;    //вычисляем относ.шероховатость
	myPipe.lambda = hydraulic_resistance_isaev(Re, relative_roughness);

	double  p_0;
	p_0 = (myPipe.p_L / (myPipe.ro * M_G) + myPipe.z_0 - myPipe.z_L + myPipe.lambda * (myPipe.L / myPipe.D) * pow(myPipe.V, 2) / (2 * M_G)) * (myPipe.ro * M_G);

	cout << "Результат: p_0 = " << p_0 << " Па" << "\n";
	return p_0;

};


//Задача Лурье 2
//ЗАДАЧА PP КЛАССИЧЕСКАЯ
//итерация по Лурье

double Q_iteraziia() {
	pipe myPipe;
	myPipe.L = 80e3;
	myPipe.D_vnesh = 720e-3;
	myPipe.b = 10e-3;
	myPipe.z_0 = 50;
	myPipe.z_L = 100;
	myPipe.ro = 870;
	myPipe.u = 15e-6;
	
	myPipe.abc = 15e-6;
	
	myPipe.D = myPipe.D_vnesh - 2 * myPipe.b;
	
	double relative_roughness = myPipe.abc / myPipe.D;    //вычисляем относ.шероховатость
	
	ofstream outFile("Q.csv");
	double p_0 = 5.65e6;
	double p_L = 0.6e6;
	double lambda_iteraziia = 1;
	myPipe.lambda = 0.02;

	while ((abs(lambda_iteraziia - myPipe.lambda) > 0.00001))
	{
		
		lambda_iteraziia = myPipe.lambda;
		myPipe.V = sqrt(((p_0 / (myPipe.ro * M_G) - (p_L / (myPipe.ro * M_G) + myPipe.z_0 - myPipe.z_L)) / (lambda_iteraziia * (myPipe.L / myPipe.D)/ (2 * M_G))));
		cout << "Скорость:   " << myPipe.V << "     " << "\n";
		double Re = myPipe.V * myPipe.D / myPipe.u;
		myPipe.lambda = hydraulic_resistance_isaev(Re, relative_roughness);
		outFile << "коэфф.гидравл.сопротивления" << lambda_iteraziia << "\n"; //начальное значение давления в трубе

	}

	myPipe.V = sqrt(((p_0 / (myPipe.ro * M_G) - (p_L / (myPipe.ro * M_G) + myPipe.z_0 - myPipe.z_L)) / (lambda_iteraziia * (myPipe.L / myPipe.D)) * (2 * M_G)));
	outFile << "скорость, м2/с:   " << myPipe.V / 4 << "\n";
	outFile << "расход, м3/с:   " << myPipe.V * (3.14 * pow(myPipe.D, 2)) / 4 << "\n";
	outFile << "расход, м3/ч:   " << myPipe.V * (3.14 * pow(myPipe.D, 2)) / 4 * 3600 << "\n";
	outFile.close();
	cout << "Расход:   " << myPipe.V * (3.14 * pow(myPipe.D, 2)) / 4 * 3600 << "  м3/ч   " << "\n";
	
	return myPipe.V;
	
};




//ЗАДАЧА 3
//Решение классической задачи QP за счет численного интегрирования дифференциального уравнения методом Эйлера

double pressure_EULER() {
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
	double Q = 3500;
	myPipe.Q = Q/3600;
	myPipe.abc = 15e-6;
	myPipe.n = 100;
	myPipe.h = myPipe.L / myPipe.n;
	myPipe.D = myPipe.D_vnesh - 2 * myPipe.b;
	                
	myPipe.V = 4 * myPipe.Q/ (3.1415 * pow(myPipe.D, 2));      //вычисляем скорость нефти 
	double Re = myPipe.V * myPipe.D / myPipe.u;
	double relative_roughness = myPipe.abc / myPipe.D;    //вычисляем относ.шероховатость
	myPipe.lambda = hydraulic_resistance_isaev(Re, relative_roughness);
	myPipe.t_w = myPipe.lambda / 8 * myPipe.ro * pow(myPipe.V, 2);

		ofstream outFile("pressure.csv");
	double p_0 = myPipe.p_0;
	double p_L;
	outFile << p_0 << "\n"; //начальное значение давления в трубе
	for (int i = 0; i < myPipe.n; ++i) {
		p_L = p_0 + myPipe.h * (-4 / myPipe.D * myPipe.t_w - myPipe.ro * M_G * (myPipe.z_L - myPipe.z_0) / ((myPipe.n - 1) * myPipe.h)); //след.значение давления в трубе
		outFile << p_L << "\n";
		p_0 = p_L;
	}
	outFile.close();
	cout << "Давление:   " << p_L << " Па" << "\n";
	return p_L;
}





//ЗАДАЧА 4
//Классическая задача PP, но вместо простой итерации задействуем метод Ньютона

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









//ЗАДАЧА 5
//Классическая задача PP поверх Эйлера на методе Ньютона

class QP_Newton_Euler : public fixed_system_t<1>
{
	/// @brief Ссылка на структуру с параметрами трубы 
	const  pipe& pipe_dannye;

using fixed_system_t<1>::var_type; public:
	/// @brief Констуктор
	/// @param pipe Ссылка на сущность трубы 
	/// @param problem Ссылка на сущность с условием задачи
	QP_Newton_Euler(const pipe& pipe_dannye) : pipe_dannye{ pipe_dannye } {};

	/// @brief Функция невязок - все члены уравнения Бернулли в правой части 
	/// @param v - скорость течения нефти, [м/с] 
	/// @return Значение функции невязок при заданной скорости 
	var_type residuals(const var_type& v)
	{ //v - искомая скорость
		double Re = v * pipe_dannye.D / pipe_dannye.u;
		double lambda = hydraulic_resistance_isaev(pipe_dannye.Re, pipe_dannye.relative_roughness);
		double t_w = lambda / 8 * pipe_dannye.ro * pow(v, 2);
		double p_0_rachet;
		double p_L = pipe_dannye.p_L;
		//последовательный расчет p_0 по методу Эйлера 
		for (int i = 0; i < pipe_dannye.n; ++i) {
			p_0_rachet = p_L - pipe_dannye.h * (-4 / pipe_dannye.D * t_w - pipe_dannye.ro * M_G * (pipe_dannye.z_L - pipe_dannye.z_0) / ((pipe_dannye.n - 1) * pipe_dannye.h));
			p_L = p_0_rachet;
		}


		double delta_p_0;

		return
		{

		   delta_p_0 = pipe_dannye.p_0 - p_0_rachet

		};

	};
};








TEST(Block_2, Task_1) {

	double pressure_p0 = pressure();
	double abs_error = 0.01e6;
	EXPECT_NEAR(5.65e6, pressure_p0, abs_error);
}


TEST(Block_2, Task_2) {

	double Q = Q_iteraziia();
	double abs_error = 4;
	EXPECT_NEAR(3500, Q, abs_error);
}




TEST(Block_2, Task_3) {

	double  p_L = pressure_EULER();
	double abs_error = 0.01e6;
	EXPECT_NEAR(0.6e6, p_L, abs_error);
}


double u_Newton() {
	pipe myPipe;
	myPipe.p_0 = 5.65e6;
	myPipe.p_L = 0.6e6;
	myPipe.L = 80e3;
	myPipe.D_vnesh = 720e-3;
	myPipe.b = 10e-3;
	myPipe.z_0 = 100;
	myPipe.z_L = 50;
	myPipe.ro = 870;
	myPipe.u = 15e-6;
	double Q = 3500;
	myPipe.Q = Q / 3600;
	myPipe.abc = 15e-6;
	myPipe.n = 100;
	myPipe.t_w = 0;
	myPipe.h = myPipe.L / myPipe.n;
	myPipe.D = myPipe.D_vnesh - 2 * myPipe.b;
	                
	myPipe.V = 4 * myPipe.Q / (3.14 * pow(myPipe.D, 2));      //вычисляем скорость нефти 
	myPipe.Re = myPipe.V * myPipe.D / myPipe.u;
	myPipe.relative_roughness = myPipe.abc / myPipe.D;    //вычисляем относ.шероховатость
	myPipe.lambda = hydraulic_resistance_isaev(myPipe.Re, myPipe.relative_roughness);
	myPipe.t_w = myPipe.lambda / 8 * myPipe.ro * pow(myPipe.V, 2);
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
	cout <<  "Скорость изначальная " << " u = " << myPipe.V << '\n' << "Классическая задача PP поверх метода Эйлера "  << " u = " << result.argument << '\n' << "Расход " << "Q = " << result.argument* (3.14 * pow(myPipe.D, 2))/4*3600 << " м3/ч " << "\n";
	return result.argument;
}



TEST(Block_2, Task_4) {
	
	//Классическая задача PP, но вместо простой итерации задействуем метод Ньютона

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

	double u = u_Newton();
	double abs_error = 10e6;
	EXPECT_NEAR(6.03e6, u, abs_error);

}




double u_Newton_Euler() {
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
	double Q = 3500;
	myPipe.Q = Q / 3600;
	myPipe.abc = 15e-6;
	myPipe.n = 100;
	myPipe.t_w = 0;
	myPipe.h = myPipe.L / myPipe.n;
	myPipe.D = myPipe.D_vnesh - 2 * myPipe.b;
	            
	myPipe.V = 4 * myPipe.Q/ (3.14 * pow(myPipe.D, 2));      //вычисляем скорость нефти 
	myPipe.Re = myPipe.V * myPipe.D / myPipe.u;
	myPipe.relative_roughness = myPipe.abc / myPipe.D;    //вычисляем относ.шероховатость
	myPipe.lambda = hydraulic_resistance_isaev(myPipe.Re, myPipe.relative_roughness);
	myPipe.t_w = myPipe.lambda / 8 * myPipe.ro * pow(myPipe.V, 2);
	//Классическая задача PP, но вместо простой итерации задействуем метод Ньютона

	//В качестве аргумента для конструктора передается myPipe
	// Создание экземпляра класса, который и будет решаемой системой
	QP_Newton_Euler test(myPipe);
	// Задание настроек решателя по умолчанию
	fixed_solver_parameters_t<1, 0> parameters;
	// Создание структуры для записи результатов расчета
	fixed_solver_result_t<1> result;
	// Решение системы нелинейныйх уравнений <1> с помощью решателя Ньютона - Рафсона
	// { 1} - Начальное приближение скорости 1 (не 0 - т.к. иначе получается деление на 0)
	fixed_newton_raphson<1>::solve_dense(test, { 1 }, parameters, &result);
	cout << "Классическая задача PP поверх Эйлера на методе Ньютона " << " u = " << result.argument << "\n";
	return result.argument;
}


TEST(Block_2, Task_5) {
	class QP_Newton_Euler : public fixed_system_t<1>
	{
		/// @brief Ссылка на структуру с параметрами трубы 
		const  pipe& pipe_dannye;
		vector<double>& massiv_dannye;

	using fixed_system_t<1>::var_type; public:
		/// @brief Констуктор
		/// @param pipe Ссылка на сущность трубы 
		/// @param problem Ссылка на сущность с условием задачи
		QP_Newton_Euler(const pipe& pipe_dannye, vector<double>& massiv_dannye) : pipe_dannye{ pipe_dannye }, massiv_dannye{ massiv_dannye } {};

		/// @brief Функция невязок - все члены уравнения Бернулли в правой части 
		/// @param v - скорость течения нефти, [м/с] 
		/// @return Значение функции невязок при заданной скорости 
		var_type residuals(const var_type& v)
		{ //v - искомая скорость
			double Re = v * pipe_dannye.D / pipe_dannye.u;
			double lambda = hydraulic_resistance_isaev(pipe_dannye.Re, pipe_dannye.relative_roughness);
			double t_w = lambda / 8 * pipe_dannye.ro * pow(v, 2);
			double p_0_rachet;
			double p_L = pipe_dannye.p_L;
			//последовательный расчет p_0 по методу Эйлера 
			for (int n = 0; n < pipe_dannye.n; ++n) {

				for (int i = 0; i < pipe_dannye.n; ++i) {
					p_0_rachet = p_L - pipe_dannye.h * (-4 / pipe_dannye.D * t_w - pipe_dannye.ro * M_G * (pipe_dannye.z_L - pipe_dannye.z_0) / ((pipe_dannye.n - 1) * pipe_dannye.h));
					p_L = p_0_rachet;
					massiv_dannye[n] = p_0_rachet;

				}

			}

			double delta_p_0;

			return
			{

			   delta_p_0 = pipe_dannye.p_0 - p_0_rachet

			};

		};
	};

}
