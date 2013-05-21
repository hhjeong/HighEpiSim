#include <vector>
#include <cmath>
using namespace std;

inline float log2( float x ) {
	return logf(x) / logf(2.0);
}

inline float entropy( float x ) {
	return x == 0 || x == 1 ? 0.0 : -x * log(x);
}
float getSingleMutualInformation( int locus, const vector< vector< int > > &G, const vector<int> &P ) {
	float ret = 0.0;

	int fCase[3] = {0};
	int fCont[3] = {0};

	int nCase = 0, nCont = 0;
	for( int i = 0 ; i < G.size() ; ++i ) {
		if( P[i] == 1 ) fCase[ G[i][locus] ]++, nCase++;
		else fCont[ G[i][locus] ]++, nCont++;
	}

	float nTot = nCase + nCont;

	ret += entropy( nCase/nTot ) + entropy( nCont/nTot );

	for( int i = 0 ; i < 3 ; ++i ) {
		ret += entropy( (fCase[i]+fCont[i])/nTot );
		ret -= entropy( fCase[i]/nTot );
		ret -= entropy( fCont[i]/nTot );
	}

	fprintf( stderr, "%d", locus );
	for( int i = 0 ; i < 3 ; ++i ) fprintf( stderr," %d", fCase[i] );
	for( int i = 0 ; i < 3 ; ++i ) fprintf( stderr," %d", fCont[i] );
	fprintf( stderr, "\n" );
	return ret;
}

void getDataFromDefaultFormat( char *fileName, vector< vector<int> > &G, vector< int > &P ) {
	G.clear(), P.clear();
	ifstream inp(fileName);

	string line;
	getline( inp, line );

	vector< string > label = split<string>( line );

	int nSNP = (int)label.size();

	while( getline(inp,line) ) {
		vector<int> row = split<int>( line );
		P.push_back( row.back() );
		row.pop_back();
		G.push_back(row);
	}

	inp.close();
}