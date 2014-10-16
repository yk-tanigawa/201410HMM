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

#if 0
    myjob -> dump();
    myjob -> viterbi();  myjob -> viterbi_dump();
#endif

    myjob -> forback_prep();

#if 0
    myjob -> forward(); myjob -> backward();
    myjob -> lpx_dump(); /* log P(x) が正しく計算できているか確認 */
#endif

    int BaumWelch_repeatNum = myjob -> BaumWelch();
    cout << "repeat num : " << BaumWelch_repeatNum << endl;
    myjob -> viterbi();  myjob -> viterbi_dump();
  }
}

