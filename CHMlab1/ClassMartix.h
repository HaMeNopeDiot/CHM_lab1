#pragma once
#include "ClassVector.h"
#include <ctime>
#include <iomanip>




class ClassMatrix : ClassVector
{
	int m_size;
protected:
	ClassVector a, b, c;
	ClassVector q, p, x;
	ClassVector f;

	ClassVector fot, xot;
public:
	ClassMatrix(vector<vector<long double>> vectors);
	ClassMatrix(int _size = 5, int _start = 1, int _end = 10);
	void sync_vectors();
	void solve_lab1();
	void solve_lab1_var2();
	void delim_cpqf_vectors(long double delim, int i);
	void minus_bcpqf_vectors(long double mult, int i);
	void print_martix(int space = 12, int precision = 8);
	ClassVector calculate_F();
	ClassVector calculate_X();
	long double calculate_falibility(ClassVector a, ClassVector b);
	long double calculate_falibility_X(ClassVector a, ClassVector b);

};