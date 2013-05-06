#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS

#include "Rng.h"
#include "Lib.h"
#include <iostream>

using namespace std;

int main() {
	InitWELLRNG512a();
	vector<double> tb = readInput("input.txt");
	vector< vector<int> > genotype;
	vector< int > phenotype;

	const int K = 3;
	const double maf = 0.1;
	const int nSNP = 2000;
	const int nCase = 2000;
	const int nCont = 2000;
	
	generateTable( K, maf, nSNP, nCase, nCont, tb, genotype, phenotype );
	writeOutput( "output.txt", nSNP, genotype, phenotype );
}