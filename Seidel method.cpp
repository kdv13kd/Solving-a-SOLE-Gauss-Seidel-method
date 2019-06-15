#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <fstream>
using namespace std;


void printSystem(float** m, float* v, int n)	//вывод системы
{
	cout << "Ax = b" << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout.width(10);		cout << setprecision(3) << m[i][j];
		}
		cout.width(15);		cout << "x" << i + 1;
		if (i == (n / 2))
		{
			cout.width(5);		cout << "=";
		}
		else
		{
			cout.width(5);		cout << " ";
		}
		cout.width(10);		cout << v[i] << endl;
	}
}
float searchMax(float* M, int n) //функция находит максимальное значение массива
{
	float  a = M[0];
	for (int k = 0; k < n; k++)
		if (abs(M[k]) > abs(a))
			a = M[k];
	return a;
}

void main()
{
	setlocale(LC_ALL, "Rus");
	ifstream fin;	
	ofstream fout;		
	fin.open("input.txt");
	if (!fin.is_open())
		cout << "Ошибка открытия файла" << endl;
	else
	{
		int n;							//размерность матрицы
		bool errorCondition = true;		//условие оценки погрешности
		float c, max;
		setlocale(LC_ALL, "Rus");
		fin >> n;
		float** a = new float* [n];		//начальная матрица А линейной системы
		float** tA = new float* [n];	//транспонированная матрица А
		float** mA = new float* [n];	//перемноженная транспонир. матрица на матрицу А
		for (int i = 0; i < n; i++)
		{
			a[i] = new float[n];
			tA[i] = new float[n];
			mA[i] = new float[n];
		}
		float* b = new float[n];		// вектор массив В линейной системы
		float* mB = new float[n];		//новый вектор В, полученный перемножением транспонир. м-цы на изначальн. вектор В
		float* x = new float[n];		//начальные приближения
		for (int i = 0; i < n; i++)
			x[i] = 0;					//заполняем введенными числами

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
				fin >> a[i][j];
			fin >> b[i];
		}
		fin.close();

		cout << "Система в матричной форме" << endl;
		printSystem(a, b, n);


		//для улучшенной работы метода зейделя следует изменить начальную систему СЛАУ, умножив
		//матрицу А и вектор В на трнаспонир. матрицу А
		//благодаря этому система станет нормальной, а матрица А положительно определенной, откуда 
		//следует сходимость метода к правильному решению
		for (int i = 0; i < n; i++)
			for (int j = i; j < n; j++)
			{
				tA[i][j] = a[j][i];			//получение транспонированной матрицы А
				tA[j][i] = a[i][j];
			}

		for (int i = 0; i < n; i++)
		{
			mB[i] = 0;
			for (int k = 0; k < n; k++)
			{
				mB[i] = mB[i] + tA[i][k] * b[k];		//умножение трансп. матрицы на массив В
				mA[i][k] = 0;
				for (int j = 0; j < n; j++)
					mA[i][k] = mA[i][k] + tA[i][j] * a[j][k];		//умножение трансп. матрицы на матрицу А изначальную
			}
		}

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
				if (i != j)
					mA[i][j] = -mA[i][j] / mA[i][i];			//согласно формулам вычисляются новые значения коэффициентов матрицы
			mB[i] = mB[i] / mA[i][i];
			mA[i][i] = 0;
		}

		//пока выполняется условие оценки погрешности цикл метода вычисляет на каждом шаге начальные приближения, изменяя и приближая их к решению системы
		while (errorCondition)
		{
			max = searchMax(x, n);  //ищется максимальное значение массива начальных приближений(решения системы) на данном шаге
			for (int i = 0; i < n; i++)
			{
				c = 0;
				for (int j = 0; j < n; j++)
					c = c + mA[i][j] * x[j];		//вычисление новых значений вектора х согласно формулам метода Зейделя
				x[i] = c + mB[i];
			}
			if (abs((searchMax(x, n)) - (max)) < 1E-10)			//оценка погрешности путем сравнения разности максимальных
				errorCondition = false;								//величин массива х на предыдущем шаге и следующем с некоторой малой величиной
		}														//если разность оказывается меньше некоторой маленькой величины, то это значит, что
		
		cout << endl << "Решение системы" << endl;

																//метод сошелся и вычисления прекращаются.
		for (int i = 0; i < n; i++)
		{
			if (abs(x[i]) < 0.0001)					//из-за погрешности численных методов компьютер часто заместо 0 получает очень маленькую величину близкую к 0
				x[i] = 0;							//которую можно заменить нормальным 0 для вывода решений системы
			cout << "x" << i + 1 << "= " << x[i] << endl;
		}

		for (int i = 0; i < n; i++)
		{
			delete[] a[i];	
			delete[] mA[i];	
			delete[] tA[i];
		}
		delete[] x;		
		delete[] b;			
		delete[] mB;
	}
	system("pause");
}
