// Vch_pract_zadanie_1.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <cstdlib>
#include <locale>
#include <cmath>
#include <vector>
#include <iomanip>
#define M_PI 3.14159265358979323846
using namespace std;
double f(double x) { return sqrt(x) - 2 * cos(M_PI * x / 2); }
double f_p (double x) { return 0.5/sqrt(x) + M_PI * sin(M_PI * x / 2); }
double f_2_p(double x) { return -0.25*pow(x, -1.5) + M_PI*M_PI*cos(M_PI * x / 2)/2; }

void bisec(double a, double b, double e) {
	//std::cout << std::setprecision(20) << std::fixed;
	cout << "Метод бисекции" << endl;
	int k = 0;
	double x_0 = (a + b) / 2, x, delta, c;
	cout << "Начальное приближение к корню: " << x_0 << endl;
	while ((b - a) > 2 * e) 
	{
		c = (a + b) / 2;
		k++;
		if (f(a) * f(c) <= 0) 
			b = c;
		else
			a = c;
	}
	x = (a + b) / 2;
	delta = (b - a) / 2;
	cout << "Количество шагов: " << k << endl;
	cout << "Приближенное решение уравнения на данном отрезке: " << x << endl;
	cout << "Разница между последними приближениями: " << 2 * delta << endl;
	cout << "Абсолютная величина невязки для приближенного решения: " << fabs(f(x)) << endl;
	return;

}

void newton(double a, double b, double e) {
	cout << "Метод Ньютона" << endl;
	double x_0 = (a + b) / 2, x; //начальное приближение
	int k = 0;
	if (f(x_0) * f_2_p(x_0) <= 0) {
		if (f(a) * f_2_p(a) > 0)
			x_0 = a;
		else
			x_0 = b;
	}
	cout << "Начальное приближение к корню: " << x_0 << endl;
	
	x = x_0 - f(x_0) / f_p(x_0);
	while (fabs(x - x_0) >= e) {
		k++;
		x_0 = x;
		x = x_0 - f(x_0) / f_p(x_0);
	}
	cout << "Количество шагов: " << k << endl;
	cout << "Приближенное решение уравнения на данном отрезке: " << x << endl;
	cout << "Разница между последними приближениями |x_" <<k<<" - x_"<<k-1<<"|= "<<fabs(x-x_0) << endl;
	cout << "Абсолютная величина невязки для приближенного решения: " << fabs(f(x)) << endl;

}

void mod_newton(double a, double b, double e) {
	cout << "Модифицированный метод Ньютона" << endl;
	double x_0 = (a + b) / 2, x, f_x0; //начальное приближение
	int k = 0;
	if (f(x_0) * f_2_p(x_0) <= 0) {
		if (f(a) * f_2_p(a) > 0)
			x_0 = a;
		else
			x_0 = b;
	}
	cout << "Начальное приближение к корню: " << x_0 << endl;
	f_x0 = f_p(x_0);
	x = x_0 - f(x_0) / f_p(x_0);
	while (fabs(x - x_0) >= e) {
		k++;
		x_0 = x;
		x = x_0 - f(x_0) / f_x0;
	}
	cout << "Количество шагов: " << k << endl;
	cout << "Приближенное решение уравнения на данном отрезке: " << x << endl;
	cout << "Разница между последними приближениями |x_" << k << " - x_" << k - 1 << "|= " << fabs(x - x_0) << endl;
	cout << "Абсолютная величина невязки для приближенного решения: " << fabs(f(x)) << endl;

}

void metod_sec (double a, double b, double e) {
	cout << "Метод секущих" << endl;
	double x_0 = a, x_1=b, x; 
	int k = 0;
	cout << "Начальное приближение к корню: " << x_0 << endl;
	x = x_1 - f(x_1) * (x_1 - x_0) / (f(x_1) - f(x_0));
	while (fabs(x - x_0) >= e) {
		k++;
		x_0 = x_1;
		x_1 = x;
		x = x_1 - f(x_1) * (x_1 - x_0) / (f(x_1) - f(x_0));
	}
	cout << "Количество шагов: " << k << endl;
	cout << "Приближенное решение уравнения на данном отрезке: " << x << endl;
	cout << "Разница между последними приближениями |x_" << k << " - x_" << k - 1 << "|= " << fabs(x - x_0) << endl;
	cout << "Абсолютная величина невязки для приближенного решения: " << fabs(f(x)) << endl;

}



int main()
{
	setlocale(LC_ALL, "Rus");
	int N,k=0;
	double A = 0, B = 4.5, h, x_1, x_2, y_1, y_2, e = pow(10, -8);
	cout << "Задана функция f(x)=sqrt(x)-2*cos(P*x/2)" << endl << "Найти количество корней на отрезке [0; 4.5] с точностью e=10^-8" << endl;
	cout << "Введите на сколько частей нужно разбить отрезок" << endl;
	cin >> N;
	vector <pair<double, double>> tab_root;
	h = (B - A) / N;
	x_1 = A; x_2 = x_1 + h; y_1 = f(x_1);

	while (x_2 <= B) {
		y_2 = f(x_2);
		if (y_1 * y_2 <= 0) {
			k++;
			tab_root.push_back(make_pair(x_1, x_2));
		}
		x_1 = x_2;
		x_2 = x_1 + h;
		y_1 = y_2;
	}
	for (vector<pair<double, double>>::iterator it = tab_root.begin(); it != tab_root.end(); it++) {
		cout << "[ " << it->first << " ; " << it->second << " ] " << endl;
	}
		cout << "Количество отрезков, содержащих корень, равняется " << k << endl;
		cout << endl;

		for (int i = 0; i < k; i++) {
			cout <<"[" << tab_root[i].first << ";" << tab_root[i].second << "]" << endl;
			bisec(tab_root[i].first,tab_root[i].second, e);
			cout << endl;
			newton(tab_root[i].first, tab_root[i].second, e);
			cout << endl;
			mod_newton(tab_root[i].first, tab_root[i].second, e);
			cout << endl;
			metod_sec(tab_root[i].first, tab_root[i].second, e);
			cout << endl;
		}
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
