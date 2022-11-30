#include "ClassMartix.h"

ClassMatrix::ClassMatrix(vector<vector<long double>> vectors)
{
	m_size = vectors[1].size();
	a.change_vector(vectors[1]); // upper   diagonal line
	b.change_vector(vectors[2]); // central   diagonal line
	c.change_vector(vectors[3]); // downer  diagonal line

	/*
	//sync p
	vectors[4][size - 1] = vectors[1][size - 1];
	vectors[4][size - 2] = vectors[2][size - 2];
	vectors[4][size - 3] = vectors[3][size - 3];
	//sync q
	vectors[5][size - 1] = vectors[2][size - 1];
	vectors[5][size - 2] = vectors[3][size - 2];
	*/

	p.change_vector(vectors[4]); // first  vertical line
	q.change_vector(vectors[5]); // second vertical line
	x.change_vector(vectors[0]); 

	sync_vectors();
}

ClassMatrix::ClassMatrix(int _size, int _start, int _end) 
{
	m_size = _size;
	srand(time(0));
	vector<vector<long double>> vectors(6);

	for (int i = 0; i < 6; i++) {
		vectors[i].push_back(-1);
		if (i == 2) {
			for (int j = 1; j <= m_size; j++) {
				vectors[i].push_back(_start + rand() % _end);
				//std::cout << std::setw(3) << vectors[i][j + 1] << " ";
			}
		} else {
			for (int j = 0; j < m_size; j++) {
				vectors[i].push_back(_start + rand() % _end);
				//std::cout << std::setw(3) << vectors[i][j + 1] << " ";
			}
			
		}
		//std::cout << '\n';
	}

	a.change_vector(vectors[0]);
	c.change_vector(vectors[1]);
	b.change_vector(vectors[2]);
	p.change_vector(vectors[3]);
	q.change_vector(vectors[4]);
	x.change_vector(vectors[5]);

	sync_vectors();
	f = calculate_F();
}

void ClassMatrix::sync_vectors()
{
	vector<long double> _a = a.get_vector();
	vector<long double> _b = b.get_vector();
	vector<long double> _c = c.get_vector();
	vector<long double> _p = p.get_vector();
	vector<long double> _q = q.get_vector();

	_b[m_size] = _q[m_size];
	_a[m_size] = _p[m_size];

	_c[m_size - 1] = _q[m_size - 1];
	_b[m_size - 1] = _p[m_size - 1];

	_c[m_size - 2] = _p[m_size - 2];
	//_a[m_size] = _a[m_size];
	//_b[m_size - 2] = _b[m_size - 2];
	//_c[m_size - 3] = _c[m_size - 3];

	//_q[m_size - 1] = _b[m_size - 1];
	//_q[m_size - 2] = _c[m_size - 2];

	p.change_vector(_p);
	q.change_vector(_q);
	a.change_vector(_a);
	b.change_vector(_b);
	c.change_vector(_c);
}

void ClassMatrix::solve_lab1()
{
}

void ClassMatrix::delim_cpqf_vectors(long double delim, int i)
{
	c.delim_elem_vector(i, delim);
	p.delim_elem_vector(i, delim);
	q.delim_elem_vector(i, delim);
	f.delim_elem_vector(i, delim);
}

void ClassMatrix::minus_bcpqf_vectors(long double mult, int i)
{

	//c.minus_elem_vector(i, c.get_elem_vector(i) - mult * 1.0f * c.get_elem_vector(i - 1));
	p.minus_elem_vector(i, mult * 1.0f * p.get_elem_vector(i - 1));
	q.minus_elem_vector(i, mult * 1.0f * q.get_elem_vector(i - 1));
	f.minus_elem_vector(i, mult * 1.0f * f.get_elem_vector(i - 1));
	if (m_size - 2 == i) {
		c.set_elem_vector(i, p.get_elem_vector(i));
	}
	if (m_size - 1 == i) {
		c.set_elem_vector(i, q.get_elem_vector(i));
	}
	b.minus_elem_vector(i, mult * 1.0f * c.get_elem_vector(i - 1));
}

void ClassMatrix::solve_lab1_var2()
{
	std::cout << "### Situation with lenght: " << m_size << " ###\n";
	//c.print_vector();
	//p.print_vector();
	print_martix();
	//std::cout << " ######### STEP 1 #########\n";
	long double b_del = b.get_elem_vector(1);
	b.set_elem_vector(1, 1.0f);
	delim_cpqf_vectors(b_del, 1);
	//print_martix();

	for (int i = 2; i < m_size; i++) {
		long double a_mult = a.get_elem_vector(i) * 1.0f / b.get_elem_vector(i - 1);
		a.set_elem_vector(i, 0);
		minus_bcpqf_vectors(a_mult, i);
		//std::cout << "minus";
		//print_martix();
		long double b_del_tmp = b.get_elem_vector(i);
		b.set_elem_vector(i, 1.0f);
		delim_cpqf_vectors(b_del_tmp, i);
		//std::cout << "multy";
		//print_martix();
	}
	fot = calculate_F();
	std::cout << "Step 1: Fallibility to F - Fot = " << setprecision(15) << calculate_falibility(f, fot) << '\n';
	//std::cout << "########### STEP 2: P to null ###########\n";
	
	long double tmp_mult = a.get_elem_vector(m_size);// *1.0f / b.get_elem_vector(m_size - 1);
	/*
	c.print_vector();
	p.print_vector();
	q.print_vector();
	a.print_vector();
	b.print_vector();
	sync_vectors();
	*/
	//std::cout << "Fallibility to F - Fot = " << setprecision(15) << calculate_falibility(f, fot) << '\n';

	a.set_elem_vector(m_size, 0);
	b.minus_elem_vector(m_size, c.get_elem_vector(m_size - 1) * tmp_mult);
	f.minus_elem_vector(m_size, f.get_elem_vector(m_size - 1) * tmp_mult);
	

	f.delim_elem_vector(m_size, b.get_elem_vector(m_size));
	b.set_elem_vector(m_size, 1.0f);
	//print_martix();
	for (int i = m_size - 3; i > 0; i--) {
		long double tmp_mult2 = p.get_elem_vector(i);// / p.get_elem_vector(m_size - 1);
		p.set_elem_vector(i, 0);
		f.minus_elem_vector(i, f.get_elem_vector(m_size - 1) * tmp_mult2);
		q.minus_elem_vector(i, q.get_elem_vector(m_size - 1) * tmp_mult2);
		//print_martix();
	}
	fot = calculate_F();
	std::cout << "Step 2: Fallibility to F - Fot = " << setprecision(15) << calculate_falibility(f, fot) << '\n';
	//std::cout << "########### STEP 3: Q to null ###########\n";
	for (int i = m_size - 2; i > 0; i--) {
		long double tmp_mult2 = q.get_elem_vector(i);
		q.set_elem_vector(i, 0);
		f.minus_elem_vector(i, f.get_elem_vector(m_size) * tmp_mult2);
		print_martix();
	}
	fot = calculate_F();
	xot = calculate_X();
	std::cout << "Step 3: Fallibility to F - Fot = " << setprecision(15) << calculate_falibility(f, fot) << '\n';
	std::cout << "Fallibility to X - Xot = " << setprecision(15) << calculate_falibility_X(x, xot) << '\n';
	std::cout << "###################################################\n";
}


void ClassMatrix::print_martix(int space, int precision)
{
	bool status = false;
	if (status) {
		std::cout << std::fixed;
		std::cout << std::setprecision(precision);
		std::cout << " ---------- Matrix! ----------\n";
		std::cout << setw(space) << b.get_elem_vector(1) << setw(space) << c.get_elem_vector(1) << setw(space);
		for (int i = 0; i < m_size - 4; i++) {
			std::cout << "0" << setw(space);
		}
		std::cout << p.get_elem_vector(1) << setw(space);
		std::cout << q.get_elem_vector(1) << setw(space);
		std::cout << setw(space * 2) << x.get_elem_vector(1) << setw(space * 2) << f.get_elem_vector(1) << '\n';
		for (int i = 2; i < m_size; i++) {
			std::cout << setw(space);
			for (int j = 0; j < i - 2; j++) {
				std::cout << "0" << setw(space);
			}
			std::cout << a.get_elem_vector(i) << setw(space) << b.get_elem_vector(i) << setw(space) << c.get_elem_vector(i) << setw(space);
			for (int j = 0; j < m_size - (int)(3 + i); j++) {
				std::cout << "0" << setw(space);
			}
			if (i < m_size - 2)
			{
				std::cout << p.get_elem_vector(i) << setw(space);
			}
			if (i < m_size - 1)
			{
				std::cout << q.get_elem_vector(i) << setw(space);
			}
			std::cout << setw(space * 2) << x.get_elem_vector(i) << setw(space * 2) << f.get_elem_vector(i) << '\n';
		}
		std::cout << setw(space);
		for (int i = 0; i < m_size - 2; i++) {
			std::cout << "0" << setw(space);
		}
		std::cout << a.get_elem_vector(m_size) << setw(space) << b.get_elem_vector(m_size) << setw(space);
		std::cout << setw(space * 2) << x.get_elem_vector(m_size) << setw(space * 2) << f.get_elem_vector(m_size) << '\n' << " ---------- ------- ----------\n";
	}
}

ClassVector ClassMatrix::calculate_F()
{
	vector<long double> _f(m_size + 1);
	for (int i = 1; i <= m_size; i++) {
		if (i == 1) {
			_f[i - 1] = b.get_elem_vector(i) * x.get_elem_vector(1) + 
					c.get_elem_vector(i) * x.get_elem_vector(2) + 
					p.get_elem_vector(i) * x.get_elem_vector(m_size - 1) + 
					q.get_elem_vector(i) * x.get_elem_vector(m_size);
		} else {
			if (i <= m_size - 3) {
				_f[i - 1] = a.get_elem_vector(i) * x.get_elem_vector(i - 1) +
					    b.get_elem_vector(i) * x.get_elem_vector(i) +
						c.get_elem_vector(i) * x.get_elem_vector(i + 1) +
						p.get_elem_vector(i) * x.get_elem_vector(m_size - 1) +
						q.get_elem_vector(i) * x.get_elem_vector(m_size);
			} else {
				if (i == m_size - 2) {
					_f[i - 1] = a.get_elem_vector(i) * x.get_elem_vector(i - 1) +
							b.get_elem_vector(i) * x.get_elem_vector(i) +
							c.get_elem_vector(i) * x.get_elem_vector(i + 1) +
							q.get_elem_vector(i) * x.get_elem_vector(m_size);
				}
				if (i == m_size - 1) {
					_f[i - 1] = a.get_elem_vector(i) * x.get_elem_vector(i - 1) +
							b.get_elem_vector(i) * x.get_elem_vector(i) +
							c.get_elem_vector(i) * x.get_elem_vector(i + 1);
				}
				if (i == m_size) {
					_f[i - 1] = a.get_elem_vector(i) * x.get_elem_vector(i - 1) +
						    b.get_elem_vector(i) * x.get_elem_vector(i);
				}
			}
		}
	}
	ClassVector f_tmp(_f);
	return f_tmp;
}

ClassVector ClassMatrix::calculate_X()
{
	vector<long double> _x(m_size + 1);
	_x[m_size] = f.get_elem_vector(m_size);
	for (int i = m_size - 1; i > 0; i--) {
		_x[i] = f.get_elem_vector(i) - c.get_elem_vector(i) * _x[i + 1];
	}

	/*_x.push_back(0);
	_x.push_back(f.get_elem_vector(m_size));
	long double temp = 0;
	for (int i = 2; i <= m_size; i++)
	{
		temp = f.get_elem_vector(m_size - i + 1) - a.get_elem_vector(m_size - i + 1) * _x[i - 2];
		_x.push_back(temp);
	}*/
	ClassVector x_tmp(_x);
	return x_tmp;
}

long double ClassMatrix::calculate_falibility(ClassVector a, ClassVector b)
{
	ClassVector fallib(m_size);
	//std::cout << "f   size " << a.get_size() << "  :";
	//a.print_vector();
	//std::cout << "fot size " << b.get_size() << "  :";
	//b.print_vector();
	for (int i = 1; i <= m_size; i++) {
		fallib.set_elem_vector(i, abs(a.get_elem_vector(i) - b.get_elem_vector(i)));
	}
	return fallib.get_norm();
}

long double ClassMatrix::calculate_falibility_X(ClassVector a, ClassVector b)
{
	ClassVector fallib(m_size);
	//std::cout << "f   size " << a.get_size() << "  :";
	//a.print_vector();
	//std::cout << "fot size " << b.get_size() << "  :";
	//b.print_vector();
	for (int i = 1; i <= m_size; i++) {
		fallib.set_elem_vector(i, abs(a.get_elem_vector(i) - b.get_elem_vector(i + 1)));
	}
	return fallib.get_norm();
}