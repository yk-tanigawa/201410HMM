#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

using namespace std;

class hmm{
  /* HMM の構造体 */
  long double **trans;  /* 遷移確率 */
  long double **emit;   /* 出力確率 */
  long double **ltrans; /* 遷移確率のlog */
  long double **lemit;  /* 出力確率のlog */
  int state_size;   /* 状態数 */
  int alph_size;    /* アルファベットの数 */
  string alph;      /* アルファベット */
public:
  void dump();
  void init(int a_size, int s_size){
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
  friend hmm &read_params(hmm &model, ifstream &ifs);
};

vector<string> split(const string &str, char delim);

template <class T> 
void show_matrix(T **matrix, string format, int n, int m){
  /* 行列の内容をprintfを使って表示する。
   * 各要素のformatも引数として与える。*/
  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++){
      printf(format.c_str(), matrix[i][j]);
    }
    printf("\n");
  }
  return;
}

vector<string> split(const string &str, char delim){
  /* strを受け取って delim でsplitしてvector<string>として返す */
  std::istringstream iss(str);
  string tmp;
  vector<string> res;
  while(getline(iss, tmp, delim)){
    res.push_back(tmp);
  }
  return res;
}

void hmm::dump(){
  /* HMM構造体の内容を表示する */
  std::cout << "number of alphabet :" << alph_size << std::endl;
  std::cout << "alphabet are       :" << alph << std::endl;
  std::cout << "number of states   :" << state_size << std::endl;
  std::cout << "transition probability matrix is as follows: " << std::endl;
  show_matrix(trans,  "%10Lf", state_size, state_size);
  std::cout << "log transition probability matrix is as follows: " << std::endl;
  show_matrix(ltrans, "%10Lf", state_size, state_size);
  std::cout << "emittion probabiliry matrix is as follows: " << std::endl;
  show_matrix(emit,   "%10Lf", state_size, alph_size);
  std::cout << "log emittion probabiliry matrix is as follows: " << std::endl;
  show_matrix(lemit,  "%10Lf", state_size, alph_size);
  return;
}

istream &getline_wocomment(char c, istream &is, string &str){
  /* getlineを行ってコメントを削除する。
   * char c以降をコメントとみなす。*/
  getline(is, str);
  int comment_start = str.find(c);
  if(comment_start >= 0)
    str.erase(str.find(c));
  return is;
}

hmm &read_params(hmm &model, ifstream &ifs){
  std::string str;
  char comment_ch = '%';

  int a_size = -1;
  getline_wocomment(comment_ch, ifs, str);
  a_size = std::stoi(str);

  std::string alphabet;
  getline_wocomment(comment_ch, ifs, alphabet);
  
  int s_size = -1;
  getline_wocomment(comment_ch, ifs, str);
  s_size = std::stoi(str);

  model.init(a_size, s_size);
  model.alph = alphabet;

  /* 状態遷移確率を読み込む */
  for(int i = 0; i < s_size; i++){
    getline_wocomment(comment_ch, ifs, str);
    vector<string> input_str = split(str, ' ');
    long double input;
    for(int j = 0; j < s_size; j++){
      sscanf((input_str[j]).c_str(), "%Lf", &input);
      model.trans[i][j] = input;
      model.ltrans[i][j] = logl(input);
      //      cout << input[j] << endl;
    }
  }

  /* 出力確率を読み込む 
   * s0の出力確率は0であることに注意。s1から読み込む */
  for(int j = 0; j < a_size; j++){
    model.emit[0][j]  = 0;
    model.lemit[0][j] = logl(0);
  }  
  for(int i = 1; i < s_size; i++){
    getline_wocomment(comment_ch, ifs, str);
    vector<string> input_str = split(str, ' ');
    long double input;
    for(int j = 0; j < a_size; j++){
      sscanf((input_str[j]).c_str(), "%Lf", &input);
      model.emit[i][j]  = input;
      model.lemit[i][j] = logl(input);
    }
  }
  return model;
}

int prepare(char *param_file, char *data_file){
  /* FILE stream を開き，パラメータファイルを読み込む 
   * その後，viterbi本体に投げる */

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

  hmm model;
  model = read_params(model, param_fs);
  model.dump();
  param_fs.close();


  return EXIT_SUCCESS;
}

int main(int argc, char *argv[]){
  return prepare((char *)"params.txt", (char *)"sample-RNA.fa");
  //return viterbi_prepare(argv[1], argv[2]);
}

