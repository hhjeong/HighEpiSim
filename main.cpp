#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS

#include "Rng.h"
#include "Lib.h"
#include "Test.h"
#include <iostream>
#include <Windows.h>

using namespace std;

int main() {
	InitWELLRNG512a();
	vector<double> tb = readInput("input.txt");
	vector< vector<int> > genotype;
	vector< int > phenotype;

	const int dim = 3;
	const int nSNPs = 2000;
	const int nTot = 4000;
	const int nReplication = 100;

	
	string modelName[] = { "Ep-1", "Ep-3", "Ep-5", "Ep-6", "Het-1", "Het-3", "S-1" };
	double het[] = { 0.02, 0.03 };
	double maf[] = { 0.1, 0.2, 0.4, 0.5 };

	for( int i = 0 ; i < 7 ; ++i ) {
		for( int j = 0 ; j < 2 ; ++j ) {
			for( int k = 0 ; k < 4 ; ++k ) {
				char inpName[100];
				sprintf( inpName, "Table//%s_%.2f_%.1f.txt", modelName[i].c_str(), het[j], maf[k] );
				vector< double > tb = readInput( inpName );
				vector< vector<int> > genotype;
				vector< int > phenotype;

				for( int rep = 1 ; rep <= nReplication ; ++rep ) {
					char oupName[100];
					sprintf( oupName, "Output//%s//%.2f_%.1f_%03d.txt", modelName[i].c_str(), het[j], maf[k], rep );
					generateBalancedData( dim, maf[k], nSNPs, nTot, tb, genotype, phenotype, rep == 1 );
					writeOutput( oupName, nSNPs, genotype, phenotype );
					fprintf( stderr,"Model[%s] with (%.2f,%.1f) %dth data ", modelName[i].c_str(), het[j], maf[k], rep );
					checkPenetrance( dim , genotype, phenotype );
				}
			}
		}

	}

	/*
	double maf[] = { 0.1, 0.2, 0.5 };

	for( int i = 0 ; i < 3 ; ++i ) {
		char inpName[100];
		sprintf( inpName, "Table2//%.1f.txt", maf[i] );
		vector< double > tb = readInput( inpName );
		vector< vector<int> > genotype;
		vector< int > phenotype;

		for( int rep = 1 ; rep <= nReplication ; ++rep ) {
			char oupName[100];
			sprintf( oupName, "Output2//%.1f//%03d.txt", maf[i], rep );
			generateBalancedData( dim, maf[i], nSNPs, nTot, tb, genotype, phenotype );
			writeOutput( oupName, nSNPs, genotype, phenotype );
			checkPenetrance( dim , genotype, phenotype );
		}
		
	}*/
	/*
	double maf[] = { 0.5 };

	for( int i = 0 ; i < 1 ; ++i ) {
		char inpName[100];
		sprintf( inpName, "Table2//%.1f.txt", maf[i] );
		vector< double > tb = readInput( inpName );
		vector< vector<int> > genotype;
		vector< int > phenotype;
		generateBalancedData( dim, maf[i], nSNPs, nTot, tb, genotype, phenotype );

		
		char oupName[100];
		sprintf( oupName, "%.1f.txt", maf[i] );

		checkPenetrance( dim, genotype, phenotype );
	}*/

	/*
	char oupName[100];
	sprintf( oupName, "test.txt" );
	getDataFromDefaultFormat( "data1.txt", genotype, phenotype );
	FILE *oup = fopen(oupName,"w");
	for( int j = 0 ; j < genotype.front().size() ; ++j ) {
		fprintf( oup, "%d %f\n", j, getSingleMutualInformation( j, genotype, phenotype ) );
	}
	fclose(oup);
	*/
}