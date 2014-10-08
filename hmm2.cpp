#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstdlib>

using namespace std;

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

template <class T>
void ary_cpy(T *target, const T *source, int len){
  for(int i = 0; i < len; i++){
    target[i] = source[i];
  }
  return;
}

template <class T>
void ary_dump(T *ary, string format, int len){
  for(int i = 0; i < len; i++){
    printf(format.c_str(), ary[i]);
  }
  printf("\n");
  return;
}

template <class T>
void ary_dump_expl(T *ary, string format, int len){
  /* 配列を表示, 各要素のexplをとってから表示
   * debug用関数 */
  for(int i = 0; i < len; i++){
    printf(format.c_str(), expl(ary[i]));
  }
  printf("\n");
  return;
}

template <class T>
int find_max_index(T *ary, int len){
  T max = ary[0];
  int index = 0;
  for(int i = 1; i < len; i++){
    if(ary[i] > max){ max = ary[i]; index = i; }
  }
  return index;
}

class hmm;
class job;
class viterbi;

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
  int a_size(){ return alph_size; }  
  int s_size(){ return state_size; }  
  long double get_ltrans(const int i, const int j){return ltrans[i][j]; }
  long double get_lemit(const int i, const int c) {return lemit[i][c];  }
  void init(const int a_size, const int s_size){
    /*
     * 構造体hmmのメモリ領域を確保する。
     * 引数には状態数と、アルファベットの数を与える
     * この関数では、callocを行うだけで、
     * パラメータのセットは別の関数で行う。
     */
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
  int alph_to_digit(char c){
    /* hmm の中に含まれるアルファベットの文字を受け取って、
     * 何番目の文字であるかindexを返す*/
    for(int i = 0; i < alph.length(); i++){
      if(alph[i] == c){ return i; }
    }
    /* HMMのアルファベットテーブルの中に見つからない */
    cerr << "data file conteins unknown alphabet " << c << endl;
    exit(EXIT_FAILURE);
  }
  int *data_convert(const string str){
    int *data = new int [str.length()];
    for(int i = 0; i < str.length(); i++){
      data[i] = this -> alph_to_digit(str[i]);
    }
    return data;
  }
  void dump(){
    /* HMM構造体の内容を表示する */
    cout << "number of alphabet :" << alph_size << endl;
    cout << "alphabet are       :" << alph << endl;
    cout << "number of states   :" << state_size << endl;
    cout << "transition probability matrix is as follows: " << endl;
    show_matrix(trans,  "%10Lf", state_size, state_size);
    cout << "log transition probability matrix is as follows: " << endl;
    show_matrix(ltrans, "%10Lf", state_size, state_size);
    cout << "emittion probabiliry matrix is as follows: " << endl;
    show_matrix(emit,   "%10Lf", state_size, alph_size);
    cout << "log emittion probabiliry matrix is as follows: " << endl;
    show_matrix(lemit,  "%10Lf", state_size, alph_size);
    return;
  }
  void destroy(){
    for(int i = 0; i < state_size; i++){
      delete [] trans[i]; delete [] ltrans[i]; 
      delete [] emit[i];  delete [] lemit[i];
    } 
    delete [] trans; delete [] ltrans;
    delete [] emit;  delete [] lemit;
  }
  friend hmm &params_read(hmm &model, ifstream &ifs);
};

class viterbi{
  long double *lv;
  long double *lv_before;
  int **tracebk;
  int *path;
  int s_size;
  int a_size;
  int length;
  hmm *model;
public:
  int init(int l, hmm *m){
    lv        = new long double [m -> s_size()];
    lv_before = new long double [m -> s_size()];
    tracebk   = new int * [l];
    for(int t = 0; t < l; t++){ tracebk[t] = new int [m -> s_size()];}
    path      = new int [l];
    s_size    = m -> s_size();
    a_size    = m -> a_size();
    length    = l;
    model = m;
    return 0;
  }
  ~viterbi(){
    delete [] lv;
    delete [] lv_before;
    for(int t = 0; t < length; t++){delete [] tracebk[t]; }
    delete [] tracebk;
    delete [] path;
  }
  int repeat(int t, int c);
  friend int viterbi_body(job &j);
};

class job{
  hmm   *model;
  int   *data;
  int    data_len;
  string data_head;
  long double *forward;
  long double *backward;
public:
  void init(hmm m, string head, int *ary, int len){
    model = &m;   data_head = head;
    data = ary;  data_len = len;  return;
  }
  ~job(){
    delete [] data;
    //model->destroy();
    return;
  }
  void dump(){ /* job classの内容を表示 */
    model->dump();
    cout << data_head << endl;
    cout << "length : " << data_len << endl;
    for(int i = 0; i < data_len; i++){ cout << data[i];}
    cout <<endl;
  }
  friend int viterbi_body(job &j);
};


vector<string> split(const string &str, char delim);

vector<string> split(const string &str, char delim){
  /* strを受け取って delim でsplitしてvector<string>として返す */
  istringstream iss(str);
  string tmp;
  vector<string> res;
  while(getline(iss, tmp, delim)){
    res.push_back(tmp);
  }
  return res;
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

hmm &params_read(hmm &model, ifstream &ifs){
  string str;
  char comment_ch = '%';

  int a_size = -1;
  getline_wocomment(comment_ch, ifs, str);
  sscanf(str.c_str(), "%d", &a_size);
  //a_size = stoi(str);

  string alphabet;
  getline_wocomment(comment_ch, ifs, alphabet);      
  
  int s_size = -1;
  getline_wocomment(comment_ch, ifs, str);
  sscanf(str.c_str(), "%d", &s_size);
  //s_size = stoi(str);

  model.init(a_size, s_size);
  for(int i = 0; i < alphabet.length(); i++){
    if(alphabet[i] != ' '){ model.alph += alphabet[i]; }
  }

  /* 状態遷移確率を読み込む */
  for(int i = 0; i < s_size; i++){
    getline_wocomment(comment_ch, ifs, str);
    vector<string> input_str = split(str, ' ');
    long double input;
    for(int j = 0; j < s_size; j++){
      sscanf((input_str[j]).c_str(), "%Lf", &input);
      model.trans[i][j] = input;
      model.ltrans[i][j] = logl(input);
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


std::string &data_read(string &data, ifstream &ifs){
  std::string str;
  while(getline(ifs, str)){ data += str; }
  return data;
}

int prepare(job &myjob, char *param_file, char *data_file){
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
  model = params_read(model, param_fs);
  param_fs.close();

  string data_str, data_head;
  getline(data_fs, data_head);
  data_str = data_read(data_str, data_fs);
  int *data =  model.data_convert(data_str);  
  
  /* job構造体をつくる */
  myjob.init(model, data_head, data, data_str.length());

  return EXIT_SUCCESS;
}

int viterbi_body(job &j){
  viterbi vit;
  vit.init(j.data_len, j.model);

  /* viterbi変数(log)の初期可 */
  vit.lv[0] = logl(1);
  for(int s = 1; s < vit.s_size; s++){
    vit.lv[s] = logl(0);
  }

  /* アルゴリズム本体を回す */
  for(int t = 0; t < vit.length; t++){
    vit.repeat(t, j.data[t]);
  }
  
  /* trace back をたどる */
  vit.path[vit.length - 1] = find_max_index(vit.lv, vit.s_size);
  for(int t = vit.length - 1; t > 0; t--){
    vit.path[t - 1] = vit.tracebk[t][vit.path[t]];
  }

  /* 結果を表示 */
  for(int t = 0; t < vit.length; t++){
    cout << vit.path[t];
  }
  cout << endl;
  return 0;
}

inline int viterbi::repeat(int t, int c){
  /* アルゴリズムの繰り返しステップ 
   * 時刻tでのviterbi変数(対数)を計算してlvに格納して，
   * traceback pointerをセットする */
  ary_cpy(lv_before, lv, s_size);  /* log viterbi変数をひとつずらす */
  for(int l = 0; l < s_size; l++){ /* 状態が k => l に遷移した */
    long double max = -1 * INFINITY;
    int max_index = -1;
    for(int k = 0; k < s_size; k++){
      long double ltrans_kl = model -> get_ltrans(k, l);
      long double temp = lv_before[k] + ltrans_kl;
      if(temp > max){ max = temp; max_index = k; }
    }
    lv[l] = model->get_lemit(l, c) + max; /* viterbi変数を計算 */
    tracebk[t][l] = max_index; /* trace back ポインタをセット */
  }
  return 0;
}

int main(int argc, char *argv[]){
  job myjob;
  prepare(myjob, (char *)"params.txt", (char *)"sample-RNA.fa");
  //myjob.dump();
  viterbi_body(myjob);
  //myjob.dump();
  return 0;
  //return viterbi_prepare(argv[1], argv[2]);
}

