#include "cpp_rw.h"
#include <pybind11/pybind11.h>

// helper function
long double next_rand_val(std::mt19937& gen, int a, int b) {
    std::uniform_real_distribution<long double> uni(a, b);

    return uni(gen);
}

void fn_rw() {

    //int N, n;
    int flag;
    int ind;
    int count;

    double lx, ly, lz;

    double R = 0.96; // 0.95
    double cutoff = 1.05; //1.121; // 0.8
    double **atoms = nullptr;
    double **bonds = nullptr;
    double *bond_new = nullptr;
    double *bond_old = nullptr;
    double *vec = nullptr;
    double denom, dist, res_cos, b_new, b_old;

    FILE *fout;

    std::random_device rd;
    std::mt19937 gen(rd() + time(0));

	FILE *fparams;
	char dummy[10];
	int n_total, N_total, num_bonds;
	int bond_1, bond_2, bond_3, bond_type;
	float f_A;
	float **components = nullptr;
	int num_components; 
	int line = 0;
	int **chains = nullptr;
	double xlo, xhi;
	double ylo, yhi;
	double zlo, zhi;
	std::vector <int> atom_types;
	int bond_types;
	int tp;

	fparams = fopen("in.parameters", "r");
	fscanf(fparams, "%s %d\n", dummy, &num_components);
	for (int i = 0; i < num_components; i++) {
		components = grow <float> (components, line + 1, 5);
		fscanf(fparams, "%s %s %s %s %s %f %f %f %f %f\n", dummy, dummy, dummy, dummy, dummy, &components[line][0], &components[line][1], &components[line][2], &components[line][3], &components[line][4]);

		std::vector<int>::iterator it;
		for (int j = 0; j < 2; j++) {
			tp = (int)components[line][2 + j];
			it = find(atom_types.begin(), atom_types.end(), tp);
			if (it == atom_types.end())
				atom_types.push_back(tp);
		}

		line++;
	}
	bond_types = (int)((float)atom_types.size()*3/2);

	printf("%d atom types, %d bond types\n", (int)atom_types.size(), bond_types);

	fscanf(fparams, "%s %s %lf %lf\n", dummy, dummy, &xlo, &xhi);
	fscanf(fparams, "%s %s %lf %lf\n", dummy, dummy, &ylo, &yhi);
	fscanf(fparams, "%s %s %lf %lf\n", dummy, dummy, &zlo, &zhi);
	fclose(fparams);

	n_total = 0;
	N_total = 0;
	num_bonds = 0;
	for (int i = 0; i < num_components; i++) {
		n_total += (int)components[i][4];
		N_total += (int)components[i][4]*(int)components[i][0];
		num_bonds += (int)components[i][4]*((int)components[i][0]-1);
	}
	printf("%d chains, %d monomer long, %d bonds\n", n_total, N_total, num_bonds);

	std::vector <int> labels;
	int count_chain = 0;
	int label;

	label = (int)next_rand_val(gen, 0, n_total);
	labels.push_back(label);
	count_chain++;

	while ((int)labels.size() < (int)components[0][4]) {
		std::vector<int>::iterator it;
		label = (int)next_rand_val(gen, 0, n_total);
		it = find(labels.begin(), labels.end(), label);

		if (it == labels.end())
			labels.push_back(label);
			count_chain++;
	}

	chains = grow <int> (chains, n_total, 2);
	for (int i = 0; i < n_total; i++) {
		chains[i][0] = i;
		chains[i][1] = 1;
	}
	for (std::vector<int>::iterator it = labels.begin(); it != labels.end(); ++it)
		chains[*it][1] = 0;

	int sum = 0;
	for (int i = 0; i < n_total; i++) {
		sum += chains[i][1];
		//printf("%d %d %d\n", chains[i][0], chains[i][1], sum);
	}

    //xlo = 0;
    //xhi = 49.99;//57;
    //ylo = 0;
    //yhi = 79.66;//90;
    //zlo = 0;
    //zhi = 14.0;// 18;

    //n = 1474;//1000; //1286;
    //N = 20;//28;

    lx = xhi - xlo;
    ly = yhi - ylo;
    lz = zhi - zlo;

    atoms = grow <double> (atoms, N_total, 6);
    bonds = grow <double> (bonds, num_bonds, 4);
    bond_new = new double[3];
    bond_old = new double[3];
    vec = new double[3];

    atoms[0][0] = 1;

	ind = 0;
    for (int i = 0; i < n_total; i++) {

		// Long running iteration; pybind11
        if (PyErr_CheckSignals() != 0)
            throw pybind11::error_already_set();

		ind = 0;
		for (int j = 0; j < i; j++) {
        	ind += (int)components[chains[j][1]][0];
		}
        printf("%.2f\n", (float)i/(float)n_total);

        // The first atom for the i-th molecule
        flag = 0;
        while (flag == 0) {

            atoms[ind][0] = ind + 1;
            atoms[ind][1] = i + 1;
            atoms[ind][2] = (int)components[chains[i][1]][2];
            atoms[ind][3] = next_rand_val(gen, xlo, xhi);
            atoms[ind][4] = next_rand_val(gen, ylo, yhi);
            atoms[ind][5] = next_rand_val(gen, zlo, zhi); //7

            // ind = 0; the first molecule, not requiring distance check
            if (ind == 0) flag = 1;
            for (int k = 0; k < ind; k++) {
                flag = 1;
                vec[0] = atoms[ind][3] - atoms[k][3];
                vec[0] = vec[0] - lx*round(vec[0]/lx); // pbc
                vec[1] = atoms[ind][4] - atoms[k][4];
                vec[1] = vec[1] - ly*round(vec[1]/ly); // pbc
                vec[2] = atoms[ind][5] - atoms[k][5];

                dist = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

                if (dist < cutoff) {
                    flag = 0;
                    break;
                }
            }
        }

        // From the second to N-th atom for the i-th molecule
        for (int j = 1; j < (int)components[chains[i][1]][0]; j++) {

            flag = 0;
            count = 0;
            while (flag == 0) {
                bond_new[0] = next_rand_val(gen, -1, 1);
                bond_new[1] = next_rand_val(gen, -1, 1);
                bond_new[2] = next_rand_val(gen, -1, 1);

                b_new = sqrt(bond_new[0]*bond_new[0] + bond_new[1]*bond_new[1] + bond_new[2]*bond_new[2]);

                bond_new[0] /= (b_new / R);
                bond_new[1] /= (b_new / R);
                bond_new[2] /= (b_new / R);

                atoms[ind + j][0] = atoms[ind + j - 1][0] + 1;
                atoms[ind + j][1] = i + 1;
				//atoms[ind + j][2] = 1;
                atoms[ind + j][2] = (j < (int) ((int)components[chains[i][1]][0] * components[chains[i][1]][1])) ? (int)components[chains[i][1]][2] : (int)components[chains[i][1]][3];
                atoms[ind + j][3] = atoms[ind + j - 1][3] + bond_new[0];
                atoms[ind + j][4] = atoms[ind + j - 1][4] + bond_new[1];
                atoms[ind + j][5] = atoms[ind + j - 1][5] + bond_new[2];

                if (atoms[ind + j][5] < cutoff) {
                    count++;
                    if (count == 100) flag = 1;
                    continue;
                }
                if (j < 2) {
                    //printf("p2");

                    /*
                    * The second atom for the i-th molecule
                    */

                    // ind = 0; the first molecule, not requiring distance check
                    if (ind == 0) flag = 1;
                    for (int k = 0; k < ind; k++) {
                        flag = 1;
                        vec[0] = atoms[ind + j][3] - atoms[k][3];
                        vec[0] = vec[0] - lx*round(vec[0]/lx); // pbc
                        vec[1] = atoms[ind + j][4] - atoms[k][4];
                        vec[1] = vec[1] - ly*round(vec[1]/ly); // pbc
                        vec[2] = atoms[ind + j][5] - atoms[k][5];

                        dist = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

                        if (dist < cutoff) {
                            flag = 0;
                            break;
                        }
                    }
                } else {
                    //printf("p3");

                    /*
                    * The third and above for the i-th molecule
                    */

                    bond_old[0] = atoms[ind + j - 1][3] - atoms[ind + j - 2][3];
                    bond_old[1] = atoms[ind + j - 1][4] - atoms[ind + j - 2][4];
                    bond_old[2] = atoms[ind + j - 1][5] - atoms[ind + j - 2][5];

                    b_new /= R; // re-using b_new already calculated before
                    b_old = sqrt(bond_old[0]*bond_old[0] + bond_old[1]*bond_old[1] + bond_old[2]*bond_old[2]);

                    // cos of the angle made by bond_new and bond_old vectors
                    res_cos = (bond_new[0]*bond_old[0] + bond_new[1]*bond_old[1] + bond_new[2]*bond_old[2])/b_new/b_old;

                    // criterion
                    if (res_cos > 0.3) {
                        if (ind == 0) flag = 1;
                        for (int k = 0; k < ind; k++) {
                            flag = 1;
                            vec[0] = atoms[ind + j][3] - atoms[k][3];
                            vec[0] = vec[0] - lx*round(vec[0]/lx); // pbc
                            vec[1] = atoms[ind + j][4] - atoms[k][4];
                            vec[1] = vec[1] - ly*round(vec[1]/ly); // pbc
                            vec[2] = atoms[ind + j][5] - atoms[k][5];

                            dist = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

                            if (dist < cutoff) {
                                flag = 0;
                                break;
                            }
                        }                        
                    }                                       
                }

                count++;

                // try 100 times
                if (count == 100) {
                    flag = 1;
                }
            }

            // If tried 100 times, redo from the previous molecule
            if (flag == 1 && count == 100) {
                i -= 1;
                break;
            }
        }
    }

    /*
    * Generate bonds
    */
    for (int i = 0; i < n_total; i++) {

		ind = 0;
		for (int j = 0; j < i; j++) {
			ind += (int)(components[chains[j][1]][0]) - 1;
		}
		//printf("%d %d\n", i, ind);

		if ((int)components[chains[i][1]][2] == 1) {
			bond_1 = 1;
			bond_2 = 2;
			bond_3 = 3;
		}
		if ((int)components[chains[i][1]][2] == 3) {
			bond_1 = 4;
			bond_2 = 5;
			bond_3 = 6;
		}
        for (int j = 0; j < (int)(components[chains[i][1]][0]) - 1; j++) {
			if (j < (int) ((int)components[chains[i][1]][0] * components[chains[i][1]][1]) - 1) {
				bond_type = bond_1;
			} else if (j == (int) ((int)components[chains[i][1]][0] * components[chains[i][1]][1]) - 1) {
				bond_type = bond_3;
			} else {
				bond_type = bond_2;
			}
            bonds[ind + j][0] = ind + j + 1;
            bonds[ind + j][1] = bond_type; 
            bonds[ind + j][2] = ind + i + j + 1;
            bonds[ind + j][3] = ind + i + j + 2;
        }
    }

    fout = fopen("data.txt", "w");
    fprintf(fout, "# data\n\n");
    fprintf(fout, "%d atoms\n", N_total);
    fprintf(fout, "%d atom types\n", (int)atom_types.size());
    fprintf(fout, "%d bonds\n", num_bonds);
    fprintf(fout, "%d bond types\n", bond_types);
    fprintf(fout, "\n");
    fprintf(fout, "%lf %lf xlo xhi\n", xlo, xhi);
    fprintf(fout, "%lf %lf ylo yhi\n", ylo, yhi);
    fprintf(fout, "%lf %lf zlo zhi\n", zlo, zhi + 20);
    fprintf(fout, "\n");

    fprintf(fout, "Atoms\n\n");
    for (int i = 0; i < N_total; i++) {
        if (atoms[i][3] < xlo) atoms[i][3] += lx;
        if (atoms[i][4] < xlo) atoms[i][4] += ly;
	    fprintf(fout, "%d %d %d %lf %lf %lf\n", (int) atoms[i][0], (int) atoms[i][1], (int) atoms[i][2], std::fmod(atoms[i][3], lx), std::fmod(atoms[i][4], ly), atoms[i][5]);
    }
    fprintf(fout, "\n");

    fprintf(fout, "Bonds\n\n");
    for (int i = 0; i < num_bonds; i++) {
        fprintf(fout, "%d %d %d %d\n", (int) bonds[i][0], (int) bonds[i][1], (int) bonds[i][2], (int) bonds[i][3]);
    }
    fclose(fout);

	free(&(components[0][0]));
	free(components);
	free(&(chains[0][0]));
	free(chains);
    free(&(atoms[0][0]));
    free(atoms);
    free(&(bonds[0][0]));
    free(bonds);
    free(bond_new);
    free(bond_old);
    free(vec);

}

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

