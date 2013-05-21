#include <vector>
#include <map>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <cstdio>
#include <string>
#include <cassert>
#include <numeric>
#include <sstream>
#include "Rng.h"

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
		for( size_t i = 0 ; i < sym.size() ; i += 2 ) {
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

int getGenotype( double maf ) {
	double MM = (1-maf) * (1-maf);
	double Mm = 2 * maf * (1-maf);
	double p = Random();

	if( p <= MM ) return 2;
	else if( p <= MM + Mm ) return 1;
	else return 0;
}

vector< int > generateRow( int dim, int len, const vector< double > &maf  ) {
	vector<int> ret = vector< int >( len );

	for( int i = 0 ; i < dim ; ++i ) ret[i] = getGenotype(maf[i]);
	for( int i = dim ; i < len ; ++i ) ret[i] = getGenotype(maf[i]);
	
	return ret;
}

vector< double > getMafs( const int &nSNPs, const int &dim, const double &maf ) {
	vector< double > ret( nSNPs );
	for( int i = 0 ; i < dim ; ++i ) ret[i] = maf;
	for( int i = dim ; i < nSNPs ; ++i ) ret[i] = Random() * 0.49 + 0.01;
	return ret;
}

void generateData( const int &dim, const double &maf, const int &nSNP,  int nCase, int nCont, const vector<double> &tb, vector< vector<int> > &genotype, vector< int > &phenotype ) {
	genotype.reserve( nCase + nCont );
	phenotype.reserve( nCase + nCont );
	genotype.clear();
	phenotype.clear();

	vector< double > mafs = getMafs( nSNP, dim, maf );

	
	int iter = 0;
	int num[2] = { nCont, nCase };
	while( num[0] > 0 || num[1] > 0 ) {
		++iter;
		vector< int > row = generateRow( dim, nSNP, mafs );

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

double getFreq( int x, double maf ) {
	if( x == 0 ) return maf * maf;
	else if( x == 1 ) return 2 * maf * (1-maf);
	else return (1-maf) * (1-maf);
}

void showTxtTable( const vector<double> &tb, const int dim = 3, FILE *oup = stderr ) {

	if( dim == 3 ) {
		const string A[] = { "AA", "Aa", "aa" };
		const string B[] = { "BB", "Bb", "bb" };
		const string C[] = { "CC", "Cc", "cc" };

		fprintf( oup, "    " );
		for( int k = 0 ; k < 3 ; ++k ) {
			
			for( int j = 0 ; j < 3 ; ++j ) fprintf( oup, "%s    ", B[j].c_str() );
			fprintf( oup, "\t" );
		}
		fprintf( oup, "\n" );
		for( int i = 2 ; i >= 0 ; --i ) {
				fprintf( oup, "%s  ", A[i].c_str() );

				for( int k = 2 ; k >= 0 ; --k ) {
					for( int j = 2 ; j >= 0 ; --j ) {
						fprintf( oup, "%.3f ", tb[k*9+i*3+j] );
					}
					fprintf( oup, "\t" );
				}
				fprintf( oup, "\n" );
		}
		
	}

}

void generateBalancedData( const int &dim, const double &maf, const int &nSNPs,  int nTot, const vector<double> &tb, vector< vector<int> > &genotype, vector< int > &phenotype, bool debug = false ) {
	assert( nTot % 2 == 0 );

	int nCase = nTot / 2;
	int nCont = nTot / 2;

	genotype.clear();
	phenotype.clear();
	genotype.resize( nTot, vector< int > ( nSNPs, 0 ) );
	phenotype.resize( nTot, 0 );

	vector< double > mafs = getMafs( nSNPs, dim, maf );
	for( int i = 0 ; i < nCase ; ++i ) {
		for( int j = dim ; j < nSNPs ; ++j ) {
			genotype[i][j] = getGenotype(mafs[j]);
		}
		phenotype[i] = 1;
	}

	for( int i = nCase ; i < nTot ; ++i ) {
		for( int j = dim ; j < nSNPs ; ++j ) {
			genotype[i][j] = getGenotype(mafs[j]);
		}
		phenotype[i] = 0;
	}

	int maxComb = 1;
	for( int i = 0 ; i < dim ; ++i ) maxComb *= 3;
	vector< double > bucket( maxComb, 0 );
	vector< double > genoFreq( maxComb, 0 );
	vector< double > genoDotPene( maxComb, 0 );
	vector< double > relGenoDotPeneCase( maxComb, 0 );
	vector< double > relGenoDotPeneCont( maxComb, 0 );

	for( int i = 0 ; i < maxComb ; ++i ) {
		double freq = 1.0;
		int _i = i;
		for( int j = 0 ; j < dim ; ++j ) {
			freq *= getFreq( _i % 3, maf );
			_i /= 3;
		}
		genoFreq[i] = freq;
		genoDotPene[i] = tb[i] * genoFreq[i];
	}

	double prevalence = accumulate( begin(genoDotPene), end(genoDotPene), 0.0 );
	for( int i = 0 ; i < maxComb ; ++i ) {
		relGenoDotPeneCase[i] = tb[i] * genoFreq[i] / prevalence;
		relGenoDotPeneCont[i] = (1-tb[i]) * genoFreq[i] / (1-prevalence);
	}
	
	bucket[0] = relGenoDotPeneCase[0];
	for( int i = 1 ; i < maxComb ; ++i ) {
		bucket[i] = bucket[i-1] + relGenoDotPeneCase[i];
	}


	for( int i = 0 ; i < nCase ; ++i ) {
		double p = Random();
		int c;
		for( c = 0 ; c < maxComb ; ++c ) {
			if( p <= bucket[c] ) break;
		}
		for( int j = 0 ; j < dim ; ++j ) {
			genotype[i][j] = c%3; 
			c /= 3;
		}
	}

	bucket[0] = relGenoDotPeneCont[0];

	for( int i = 1 ; i < maxComb ; ++i ) {
		bucket[i] = bucket[i-1] + relGenoDotPeneCont[i];
	}

	for( int i = nCase ; i < nTot ; ++i ) {
		double p = Random();
		int c;
		for( c = 0 ; c < maxComb ; ++c ) {
			if( p <= bucket[c] ) break;
		}
		for( int j = 0 ; j < dim ; ++j ) {
			genotype[i][j] = c%3; 
			c /= 3;
		}
	}	

	if( debug ) {
		showTxtTable( relGenoDotPeneCase, dim );
		showTxtTable( relGenoDotPeneCont, dim );

		vector< double > fCase(maxComb,0), fCont(maxComb,0);

		for( int i = 0 ; i < nCase ; ++i ) {
			int c = 0;
			for( int j = 0 ; j < dim ; ++j ) c = 3 * c + genotype[i][j];
			fCase[c]++;
		}

		for( int i = nCase ; i < nTot ; ++i ) {
			int c = 0;
			for( int j = 0 ; j < dim ; ++j ) c = 3 * c + genotype[i][j];
			fCont[c]++;
		}

		for( int i = 0 ; i < maxComb ; ++i ) fCase[i] /= nCase, fCont[i] /= nCont;

		showTxtTable( fCase, dim );
		showTxtTable( fCont, dim );

	}
}




void writeOutput( char *fileName, int nSNP, vector< vector<int> > &genotype, vector< int > &phenotype ) {
	FILE *oup = fopen(fileName,"w");

	for( int i = 0 ; i < nSNP ; ++i ) {
		fprintf( oup, "X%d\t", i );
	}

	fprintf( oup, "Class\n");

	for( size_t i = 0 ; i < genotype.size() ; ++i ) {
		for( size_t j = 0 ; j < genotype[i].size() ; ++j ) {
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


	for( size_t i = 0 ; i < genotype.size() ; ++i ) {
		int c = 0;
		for( int j = 0 ; j < dim ; ++j ) c = 3 * c + genotype[i][j];
		if( phenotype[i] == 1 ) nCase[c]++;
		else nCont[c]++;
	}

	vector< double > ret;
	vector< double > genof;
	
	double prevalence = 0.0;
	for( int i = 0 ; i < maxComb ; ++i ) {
		int tot = nCase[i] + nCont[i];	
		if( tot == 0 ) {
			genof.push_back( tot / (double)genotype.size()  );
			ret.push_back(0);
		}
		else {
			ret.push_back( double(nCase[i]) / tot );
			genof.push_back( tot / (double)genotype.size()  );
		}
		prevalence += ret.back() * genof.back();
	}

	double het = 0.0;
	for( int i = 0 ; i < maxComb ; ++i ) {
		het += ( ret[i] - prevalence ) * ( ret[i] - prevalence ) * genof[i];
	}

	//for( int i = 0 ; i < maxComb ; ++i ) cerr << nCase[i] << " " << nCont[i] << endl;
	het /= prevalence * (1-prevalence);
	cerr << prevalence << " " << het << endl;
	return ret;
}


template <typename T>
vector< T >  split( const string &s ) {
	vector< T > ret;
	istringstream inp(s);
	T tok;
	while(inp>>tok) ret.push_back(tok);
	return ret;
}