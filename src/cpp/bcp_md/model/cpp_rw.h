#ifndef CPP_RW_H
#define CPP_RW_H

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cmath>

#include <vector>
#include <algorithm>
#include <random>

void fn_rw();

void *smalloc(int nbytes);
void *srealloc(void *ptr, int nbytes); 
template <typename T> T **create(T **&array, int n1, int n2);
template <typename T> T **grow(T **&array, int n1, int n2);

#endif
