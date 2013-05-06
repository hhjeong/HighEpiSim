#include <vector>
#include <map>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <cstdio>
#include <string>
#include <ppl.h>

#include "Rng.h"

using namespace concurrency;
using namespace std;

void exitWithError( const string &msg ) {
	cerr << "Error: " << msg << endl;
	exit(1);
}

vector< double > readInput( char *fileName ) {
	map< int, double > freq;

	int maxi = -1;


	ifstream inp(fileName);
	string sym;
	double val;
	while( inp >> sym >> val ) {
		int c = 0;
		for( int i = 0 ; i < sym.size() ; i += 2 ) {
			c = 3*c + (int)isupper(sym[i]) + (int)isupper(sym[i+1]);
		}
		freq[c] = val;
		maxi = max( maxi, c );
	}
	inp.close();

	vector< double > ret( maxi+1, 0 );

	for( int i = 0 ; i <= maxi ; ++i ) {
		if( freq.count(i) == 0 ) {
			exitWithError("Invalid inputfile!");
		}
		ret[i] = freq[i];
	}

	return ret;
}

int getMarkerGenotype( double maf ) {
	double MM = (1-maf) * (1-maf);
	double Mm = 2 * maf * (1-maf);
	double p = Random();

	if( p <= MM ) return 2;
	else if( p <= MM + Mm ) return 1;
	else return 0;
}

int getNormalGenotype() {
	double p = Random();
	if( p <= 0.25 ) return 2;
	else if( p <= 0.75 ) return 1;
	else return 0;
}

vector< int > generateRow( int dim, int len, double maf ) {
	vector<int> ret = vector< int >( len );

	for( int i = 0 ; i < dim ; ++i ) ret[i] = getMarkerGenotype(maf);
	for( int i = dim ; i < len ; ++i ) ret[i] = getNormalGenotype();
	
	return ret;
}

void generateTable( const int &dim, const double &maf, const int &nSNP,  int nCase, int nCont, const vector<double> &tb, vector< vector<int> > &genotype, vector< int > &phenotype ) {
	genotype.reserve( nCase + nCont );
	phenotype.reserve( nCase + nCont );
	genotype.clear();
	phenotype.clear();

	int iter = 0;
	int num[2] = { nCont, nCase };
	while( num[0] > 0 || num[1] > 0 ) {
		++iter;
		vector< int > row = generateRow( dim, nSNP, maf );

		int wh = 0;
		for( int i = 0 ; i < dim ; ++i ) {
			wh = 3 * wh + row[i];
		}

		int label = Random() <= tb[wh];

		if( num[label] > 0 ) {
			genotype.push_back( row );
			phenotype.push_back( label );
			num[label]--;
		}
	}
	cerr << "Iteration : " << iter << endl;
}

void writeOutput( char *fileName, int nSNP, vector< vector<int> > &genotype, vector< int > &phenotype ) {
	FILE *oup = fopen(fileName,"w");

	for( int i = 0 ; i < nSNP ; ++i ) {
		fprintf( oup, "X%d\t", i );
	}

	fprintf( oup, "Class\n");

	for( int i = 0 ; i < genotype.size() ; ++i ) {
		for( int j = 0 ; j < genotype[i].size() ; ++j ) {
			fprintf(oup,"%d\t",genotype[i][j]);
		}
		fprintf(oup,"%d\n",phenotype[i]);
	}
	fclose(oup);
}

vector<double> checkPenetrance( const int &dim, const vector< vector<int> > &genotype, vector< int > &phenotype ) {
	int maxComb = 1;
	for( int i = 0 ; i < dim ; ++i ) maxComb *= 3;

	vector<int> nCase( maxComb, 0 );
	vector<int> nCont( maxComb, 0 );


	for( int i = 0 ; i < genotype.size() ; ++i ) {
		int c = 0;
		for( int j = 0 ; j < dim ; ++j ) c = 3 * c + genotype[i][j];
		if( phenotype[i] == 1 ) nCase[c]++;
		else nCont[c]++;
	}

	vector< double > ret;
	
	for( int i = 0 ; i < maxComb ; ++i ) {
		int tot = nCase[i] + nCont[i];
		
		cerr << nCase[i] << " " << nCont[i] << endl;
		if( tot == 0 ) ret.push_back(0);
		else {
			ret.push_back( double(nCase[i]) / tot );
		}
	}
	return ret;
}
