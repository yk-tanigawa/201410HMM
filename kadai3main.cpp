#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cfloat>
#include "hmm.hpp"

using namespace std;


int main(int argc, char *argv[]){
  if(argc < 2){
    cerr << "usage: $" << argv[0] 
	 << " <parameter file> <FASTA file>" << endl;
    return EXIT_FAILURE;
  }else{
    job *myjob = read_from_input_file(argv[1], argv[2]);

#if 1
    myjob -> forback_prep();
    int BaumWelch_repeatNum = myjob -> BaumWelch();
#if 0
    myjob -> dump(); /* Baum Welchの実行後の確率を表示する。 */
#endif
    cout << "repeat num : " << BaumWelch_repeatNum << endl;
#endif
    myjob -> viterbi();  myjob -> viterbi_dump(); myjob -> viterbi_delete();
  }
}

