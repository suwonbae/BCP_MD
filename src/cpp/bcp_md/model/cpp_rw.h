#ifndef CPP_RW_H
#define CPP_RW_H

#include <cstdio>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <random>
#include <algorithm>
#include <cmath>

void rw_test();

class randomwalk {
public:
	randomwalk();

	virtual ~randomwalk() {}

	void setN( int );
	void printN();

	void readParams( std::string );
	void generate();
	void write();

	void *smalloc(int nbytes){
		void *ptr = malloc(nbytes);

		return ptr;
	}

	void *srealloc(void *ptr, int nbytes) {
		ptr = realloc(ptr, nbytes);

		return ptr;
	}

	template <typename T> T **create(T **&array, int n1, int n2) {

		int nbytes = sizeof(T) * n1 * n2;
		T *data = (T *) smalloc(nbytes);
		nbytes = sizeof(T *) * n1;
		array = (T **) smalloc(nbytes);

		int n = 0;
		for (int i = 0; i < n1; i++) {
			array[i] = &data[n];
			n += n2;
		}

		return array;
	}

	template <typename T> T **grow(T **&array, int n1, int n2) {

		if (array == nullptr) return create(array, n1, n2);

		int nbytes = sizeof(T) * n1 * n2;
		T *data = (T *) srealloc(array[0], nbytes);
		nbytes = sizeof(T *) * n1;
		array = (T **) srealloc(array, nbytes);

		int n = 0;
		for (int i = 0; i < n1; i++) {
			array[i] = &data[n];
			n += n2;
		}

		return array;
	}

protected:
	int N;

	int flag, ind, count;
	int num_components, line;

	double R, cutoff;

	double dist, res_cos, b_new, b_old;
	
	int **chains;
	double **atoms, **bonds, **components;
	double *bond_new, *bond_old, *vec;

	FILE *fout, *fparams;
	char dummy[10];

	int n_total, N_total, num_bonds;
	int bond_1, bond_2, bond_3, bond_type;

	double xlo, xhi, ylo, yhi, zlo, zhi, lx, ly, lz;
	int bond_types, tp;

	std::vector <int> atom_types;

};

#endif
