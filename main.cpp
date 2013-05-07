#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS

#include "Rng.h"
#include "Lib.h"
#include <iostream>
#include <Windows.h>
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
	
	DWORD st = GetTickCount();
	//generateData( K, maf, nSNP, nCase, nCont, tb, genotype, phenotype );
	generateBalancedData( K, maf, nSNP, nCase+nCont, tb, genotype, phenotype );
	DWORD ed = GetTickCount();

	cerr << "Time Elapsed : " << (ed-st) / 1000.0 << "sec" << endl;
	// writeOutput( "output.txt", nSNP, genotype, phenotype );

	
	vector< double > penetrance = checkPenetrance( K, genotype, phenotype );
	/*
	for( size_t i = 0 ; i < penetrance.size() ; ++i ) {
		printf("%f\n", penetrance[i]);
	}*/
}