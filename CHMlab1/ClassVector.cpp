#include "ClassVector.h"
#include <cmath>

ClassVector::ClassVector(int _size, int _mode)
{
	v_size = _size + 1;
	current_vector.push_back(0);
	for (int i = 0; i < v_size; i++) {
		current_vector.push_back(_mode);
	}
}

ClassVector::ClassVector(int _start_d, int _end_d, int _size)
{
	v_size = _size + 1;
	current_vector.push_back(0);
	for (int i = 0; i < v_size; i++) {
		current_vector.push_back(_start_d + rand() % _end_d);
	}
}

ClassVector::ClassVector(vector<long double> external_vector)
{
	v_size = external_vector.size();
	current_vector.push_back(0);
	for (int i = 0; i < v_size; i++) {
		current_vector.push_back(external_vector[i]);
	}
}

void ClassVector::print_vector()
{
	for (int i = 1; i < v_size; i++)
	{
		cout << current_vector[i] << ' ';
	}
	cout << '\n';
}

bool ClassVector::vector_size_equal(ClassVector other_vector)
{
	return other_vector.current_vector.size() == current_vector.size();
}

//long double ClassVector::get_elem_from_vector(int num)
//{
//	if (num < current_vector.size())
//	{
//		return current_vector[num];
//	}
//	cout << "Get incorrect num!\n";
//	return -1;
//}

//void ClassVector::set_elem_to_vector(int num, long double elem)
//{
//	if (num < current_vector.size())
//	{
//		current_vector[num] = elem;
//	}
//	cout << "Get incorrect num!\n";
//}

long double ClassVector::get_norm()
{
	long double norm = 0;
	for (int i = 1; i < current_vector.size(); i++)
	{
		norm = norm + current_vector[i] * current_vector[i] * 1.0f;
	}
	norm = sqrt(norm);
	return norm;
}

void ClassVector::change_vector(vector<long double> new_vector)
{
	current_vector.clear();
	v_size = new_vector.size();
	current_vector.push_back(0);
	for (int i = 1; i < v_size; i++)
	{
		current_vector.push_back(new_vector[i]);
	}
}

void ClassVector::set_elem_vector(int index, long double digit)
{
	if (index < v_size) {
		current_vector[index] = digit;
	}
}

long double ClassVector::get_elem_vector(int index)
{
	if (index < v_size) {
		return current_vector[index];
	}
}

void ClassVector::delim_elem_vector(int index, long double digit)
{
	if (index < v_size) {
		//long double k = 1.0f / digit;
		//current_vector[index] *= k;
		if (digit < 0) {
			current_vector[index] = ((long double)current_vector[index] * 1.0f) / ((long double)abs(digit) * 1.0f);
			current_vector[index] *= -1;
		} else {
			current_vector[index] = ((long double)current_vector[index] * 1.0f) / ((long double)digit * 1.0f);
		}
		
	}
}
void ClassVector::multyply_elem_vector(int index, long double digit)
{
	if (index < v_size) {
		current_vector[index] = current_vector[index] * 1.0f * digit;
	}
}
void ClassVector::sum_elem_vector(int index, long double digit)
{
	if (index < v_size) {
		current_vector[index] = current_vector[index] + digit;
	}
}
void ClassVector::minus_elem_vector(int index, long double digit)
{
	if (index < v_size) {
		current_vector[index] = current_vector[index] - digit;
	}
}
