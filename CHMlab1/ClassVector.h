#pragma once
#include <vector>
#include <iostream>
#include <iomanip>
using namespace std;

class ClassVector
{
	int v_size;
protected:
	vector<long double> current_vector;
public:
	ClassVector(int _start_d, int _end_d, int _size = 3);
	ClassVector(int _size = 3, int _mode = 0);
	ClassVector(vector<long double> external_vector);

	void set_elem_vector(int index, long double digit);
	long double get_elem_vector(int index);
	int get_size() { return current_vector.size(); }

	long double get_norm();
	//long double get_elem_from_vector(int num);
	//void set_elem_to_vector(int num, long double elem);
	bool vector_size_equal(ClassVector other_vector);
	void print_vector();
	void change_vector(vector<long double> new_vector);
	vector<long double> get_vector() { return current_vector; }

	void delim_elem_vector(int index, long double digit);
	void multyply_elem_vector(int index, long double digit);
	void sum_elem_vector(int index, long double digit);
	void minus_elem_vector(int index, long double digit);

	ClassVector operator+(const ClassVector& other_vector)
	{
		for (int i = 0; i < v_size; i++) {
			current_vector[i] = current_vector[i] + other_vector.current_vector[i];
		}
		return ClassVector(current_vector);
	}
	

	ClassVector operator*(const ClassVector& other_vector)
	{
		long double sum = 0;
		for (int i = 0; i < v_size; i++) {
			sum += current_vector[i] * other_vector.current_vector[i];
		}
		return ClassVector(current_vector);
	}
	ClassVector operator-(const ClassVector& other_vector)
	{
		for (int i = 0; i < v_size; i++) {
			current_vector[i] = current_vector[i] - other_vector.current_vector[i];
		}
		return ClassVector(current_vector);
	}

	

};


