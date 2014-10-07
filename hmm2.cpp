#include <iostream>
#include <string>
#include <fstream>

using namespace std;


class hmm{
  /* HMM の構造体 */
  long double **trans;  /* 遷移確率 */
  long double **emit;   /* 出力確率 */
  long double **ltrans; /* 遷移確率のlog */
  long double **lemit;  /* 出力確率のlog */
  int state_size;  /* 状態数 */
  int alph_size;   /* アルファベットの数 */
  char *alph;      /* 出力アルファベット記号 */
public:
  hmm(int s_size, int a_size){
    /*
     * 構造体hmmのメモリ領域を確保する。
     * 引数には状態数と、アルファベットの数を与える
     * この関数では、callocを行うだけで、
     * パラメータのセットは別の関数で行う。
     */
    alph   = new char [a_size];
    trans  = new long double * [s_size];
    ltrans = new long double * [s_size];
    emit   = new long double * [s_size];
    lemit  = new long double * [s_size];
    for(int i = 0; i < s_size; i++){
      trans[i]  = new long double [s_size];
      ltrans[i] = new long double [s_size];
      emit[i]   = new long double [a_size];
      lemit[i]  = new long double [a_size];
    }
    state_size = s_size;
    alph_size = a_size;
    return;
  }
  friend hmm *read_params(ifstream ifs);
};


istream &getline_wocomment(istream &is, string &str, char c){
  /* getlineを行ってコメントを削除する。
   * char c以降をコメントとみなす。*/
  getline(is, str);
  int comment_start = str.find(c);
  if(comment_start >= 0)
    str.erase(str.find(c));
  return is;
}


hmm *read_params(ifstream &ifs){
  std::string str;
  char comment_ch = '%';

  getline_wocomment(ifs, str, comment_ch);
  

#if 0
  while(getline_wocomment(ifs, str, comment_ch)){
    std::cout << str << std::endl;
  }
#endif
  return NULL;
}

int viterbi_prepare(char *param_file, char *data_file){
  std::ifstream param_fs(param_file);
  if (param_fs.fail()){
    std::cerr << "cannot open params file" << std::endl;
    return EXIT_FAILURE;
  }
  std::ifstream data_fs(data_file);
  if (data_fs.fail()){
    std::cerr << "cannot open data file" << std::endl;
    return EXIT_FAILURE;
  }

  hmm *model = read_params(param_fs);

  return EXIT_SUCCESS;
}

int main(int argc, char *argv[]){
  return viterbi_prepare((char *)"params.txt", (char *)"sample-RNA.fa");
  //return viterbi_prepare(argv[1], argv[2]);
}

