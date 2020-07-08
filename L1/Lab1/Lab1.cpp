/*Подключение директив*/
#include <iostream>
#include <malloc.h>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <locale>
using namespace std;

/*define-замены*/
#define M 10
#define e 0.00001

int _SizeMatrix = 0; //размер матрицы
int _count_swap = 0;
float MatrixA[M][M + 1]; //матрица A
float CMatrixA[M][M + 1]; //копия матрицы A
float Answers[M]; //вектор ответов
float Discrepancy[M]; //вектор величччин невязки
float InverseMatrix[M][M]; //обратная матрица
float A[M]; //матрица Cpk

/*Чтение матрицы из файла*/
void Init() {
	string num;
	cout << "Введите номер файла с матрицей: "; cin >> num;

	ifstream ifs;
	ifs.open("file" + num + ".txt"); 
	ifs >> _SizeMatrix;

	//заносим матрицу
	for (int i = 0; i < _SizeMatrix; i++)
		for (int j = 0; j < _SizeMatrix; j++)
			ifs >> MatrixA[i][j];
	ifs.close();
}

/*Чтение вектора ответов из файла*/
void InitAnswer() {
	string num;
	cout << "Введите пожалуйста номер файла с ответами: "; cin >> num;

	ifstream ifs;
	ifs.open("file_ans" + num + ".txt");

	//расширяем матрицу
	for (int j = 0; j < _SizeMatrix; j++)
		ifs >> MatrixA[j][_SizeMatrix];
}

/*Копирование расширенной матрицы системы*/
void CopyMatrix() {
	for (int i = 0; i < _SizeMatrix; i++)
		for (int j = 0; j < _SizeMatrix + 1; j++)
			CMatrixA[i][j] = MatrixA[i][j];
}

/*Проверяем на соответсвие точности*/
bool Equality(float a, float b) { return abs(a - b) < e ? true : false; }

/*Перестановка строк k, l местами*/
void SwapLines(int k, int l) {
	for (int i = 0; i < _SizeMatrix + 1; i++)
		swap(CMatrixA[k][i], CMatrixA[l][i]);
	_count_swap++;
}

/*Поиск главного элемента(большего по модулю в столбце)*/
float SelectingMainElement(int s)	{
	float max = abs(CMatrixA[s][s]); //максимальное число по модулю
	int num_row = s;  //номер рассматриваемой строки

	for (int i = s + 1; i < _SizeMatrix; i++)
		//если нашли новое максимальное число по модулю
		if (abs(CMatrixA[i][s]) > max) {
			num_row = i; //запоминаем новую строку
			max = abs(CMatrixA[i][s]); //запоминаем новый главный элемент
		}

	if (Equality(max, 0)) return 0; 

	if (num_row != s) {
		SwapLines(num_row, s);
		return -CMatrixA[s][s];
	}
	return CMatrixA[s][s];
}

/*Вычитание строки l*С из строки k*/
void SubLines(int k, int l, float С) {
	for (int i = 0; i < _SizeMatrix + 1; i++)
		CMatrixA[k][i] -= CMatrixA[l][i] * С;
}

/*Обратный ход Гауса или вычисиление ответа*/
void ReverseCourse() {
	//создаем вектор и зануляем его
	for (int j = 0; j < M; j++) Answers[j] = 0;

	for (int i = _SizeMatrix - 1; i >= 0; i--) {
		float s = 0;

		for (int j = _SizeMatrix - 1; j > i; j--) {
			s += Answers[j] * CMatrixA[i][j];
		}

		//подсчет Xk, вычисление ответа
		Answers[i] = (CMatrixA[i][_SizeMatrix] - s) / CMatrixA[i][i];
	}
}

/*Решение системы уравнений методом Гаусса*/
float Gauss() {
	float det = 1;

	for (int i = 0; i < _SizeMatrix; i++) {
		det *= SelectingMainElement(i);	//выбор главного элемента

		//проверяем на вырожденность
		if (Equality(det, 0)) {
			cout << "Определитель матрицы равен 0! Матрица вырождена.";
			return 0;
		}

		//считаем Cpk по формуле и вычитаем из p строку k умноженное на Cpk
		//зануление нижк главной диагонали
		for (int j = i + 1; j < _SizeMatrix; j++) {
			A[j] = CMatrixA[j][i] / CMatrixA[i][i];
			SubLines(j, i, A[j]);
		}
	}

	ReverseCourse();
	return det;
}


/*Подсчет величины невязки*/
void DiscrepancyCalculation() {
	for (int i = 0; i < _SizeMatrix; i++) {
		float s = 0;
		for (int j = 0; j < _SizeMatrix; j++)
			s += Answers[j] * MatrixA[i][j];
		Discrepancy[i] = MatrixA[i][_SizeMatrix] - s;
	}
}

/*Вычисление обратной матрицы*/
void InverseMatrixCalculation(){
	//решаем m матриц с одинаковой матрицей A и разными спец.частями
	//где k столбец обратной матрицы - это решение A со спец частью

	for (int i = 0; i < _SizeMatrix; i++) {
		CopyMatrix();

		for (int k = 0; k < _SizeMatrix; k++)
			CMatrixA[k][_SizeMatrix] = 0; //зануляем спец часть

		CMatrixA[i][_SizeMatrix] = 1; //ставим 1 в k позицию спец части

		Gauss();	

		for (int k = 0; k < _SizeMatrix; k++)
			InverseMatrix[k][i] = Answers[k];
	}
}

/*Вывод результата*/
void PrintAnswer() {
	ofstream ofs;
	ofs.open("Answer.txt");

	//Вывод системы уравнений
	cout << "Для системы уравнений вида: " << endl;
	for (int i = 0; i < _SizeMatrix; i++) {
		for (int j = 0; j < _SizeMatrix; j++)
			cout <<  MatrixA[i][j] << "*x" << setw(9) << left << j + 1 << " ";
		cout << setw(9) << left << " = " << MatrixA[i][_SizeMatrix];
		cout << endl;
	}
	cout << endl;

	//Вывод треугольной матрицы
	cout << "Треугольная матрица имеет вид: " << endl;
	ofs << "Треугольная матрица имеет вид: " << endl;
	for (int i = 0; i < _SizeMatrix; i++) {
		for (int j = 0; j < _SizeMatrix; j++) {
			cout << setw(9) << left << CMatrixA[i][j] << " ";
			ofs << CMatrixA[i][j] << " ";
		}
		cout << endl;
		ofs << endl;
	}
	cout << endl;
	ofs << endl;

	//Вывод X-ов(ответов)
	cout << "Решение системы уравнений: " << endl;
	ofs << "Решение системы уравнений: " << endl;
	for (int i = 0; i < _SizeMatrix; i++) {
		cout << "x" << i + 1 << " = " << setw(9) << left << Answers[i] << endl;
		ofs << Answers[i] << endl;
	}
	cout << endl;
	ofs << endl;

	//вывод невязки
	cout << "Величина невязки:" << endl;
	ofs << "Величина невязки:" << endl;
	for(int i = 0; i < _SizeMatrix; i++) {
		cout << "x" << i + 1  << " ~ " << setw(9) << left << Discrepancy[i] << endl;
		ofs << Discrepancy[i] << endl;
	}
	cout << endl;
	ofs << endl;

	//вывод обратной матрицы
	InverseMatrixCalculation();
	cout << "Обратная матрица матрице A: " << endl;
	ofs << "Обратная матрица матрице A: " << endl;
	for (int i = 0; i < _SizeMatrix; i++){
		for (int j = 0; j < _SizeMatrix; j++) {
			cout << setw(9) << left << InverseMatrix[i][j] << " ";
			ofs << InverseMatrix[i][j] << " ";
		}
		cout << endl;
		ofs << endl;
	}
	cout << endl;
	ofs << endl;
	ofs.close();
}

/*Основная функция*/
int main() {
	setlocale(LC_ALL, "rus");

	Init(); // чтение матрицы A
	if (_SizeMatrix == 0) {
		cout << "Матрица A - пустая!Заполните файл." << endl;;
		return 1;
	}
	else
		cout << "Матрица A успешно считана из файла." << endl;

	InitAnswer(); // чтение ответов b
	cout << "Матрица ответов b успешно считана из файла." << endl;

	CopyMatrix(); //копирование расширенной матрицы системы

	float Det = Gauss(); //решение методом Гауса(D -определитель матрицы A)

	if (Equality(Det, 0)) return 0;

	cout << endl << "Результат вычислиней по методу Гауса:" << endl;
	if(_count_swap % 2 != 0)
		cout << "Определитель = " << Det*(-1) << endl << endl;
	else
		cout << "Определитель = " << Det << endl << endl;
	DiscrepancyCalculation(); //вычисление величины невязки
	PrintAnswer();

	system("pause");
	return 0;
}

