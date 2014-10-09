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
void show_matrix_expl(T **matrix, string format, int n, int m){
  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++){
      printf(format.c_str(), expl(matrix[i][j]));
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
  /* 長さlenの配列の要素で，最大のものが格納されているindexを返す */
  T max = ary[0];
  int index = 0;
  for(int i = 1; i < len; i++){
    if(ary[i] > max){ max = ary[i]; index = i; }
  }
  return index;
}

class hmm;
class sequence;
class job;
class viterbi;
class forward_backward;

class hmm{
  /* HMM の構造体 */
  long double **trans;  /* 遷移確率 */
  long double **emit;   /* 出力確率 */
  long double **ltrans; /* 遷移確率のlog */
  long double **lemit;  /* 出力確率のlog */
  int state_size;   /* 状態数 */
  int alph_size;    /* アルファベットの数 */
  string alph;      /* アルファベット */
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
public:
  int a_size(){ return alph_size; }  
  int s_size(){ return state_size; }  
  long double get_ltrans(const int i, const int j){return ltrans[i][j]; }
  long double get_trans(const int i, const int j){return trans[i][j]; }
  long double get_lemit(const int i, const int c) {return lemit[i][c];  }
  long double get_emit(const int i, const int c) {return emit[i][c];  }
  int *data_convert(const string str){ /* str のデータをint *に変換 */
    int *data = new int [str.length()];
    for(int i = 0; i < str.length(); i++){
      data[i] = this -> alph_to_digit(str[i]);
    }
    return data;
  }
  void init(const int a_size, const int s_size){
    /* 構造体hmmのメモリ領域を確保する。
     * 引数には状態数と、アルファベットの数を与える
     * この関数では、callocを行うだけでパラメータのセットは別の関数で行う */
    trans  = new long double * [s_size]; ltrans = new long double * [s_size];
    emit   = new long double * [s_size]; lemit  = new long double * [s_size];
    for(int i = 0; i < s_size; i++){
      trans[i]  = new long double [s_size];
      ltrans[i] = new long double [s_size];
      emit[i]   = new long double [a_size];
      lemit[i]  = new long double [a_size];
    }
    state_size = s_size; alph_size = a_size;
    return;
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
  friend hmm *params_read(ifstream &ifs);
};

class viterbi{
  /* viterbi algorithmのためのクラス */
  long double *lv;
  long double *lv_before;
  int **tracebk;
  int *path; /* 求めたViterbi pathを入れる */
  int s_size;
  int a_size;
  int length; /* データの長さ */
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

class forward_backward{
  long double **tbl; /*logをとった後の値を保存*/
  long double *scale;
  int len;
  int s_size;
  hmm *model;
public:
  forward_backward(const class job);
  int forward(const sequence data, int t_len);
  int backward(const sequence data, int t_len);
  long double calc(int t_len);
  void show(int t_len);
  long double get_tbl(int i, int j){ return tbl[i][j];}
};

class sequence{
  int  len;
  string head;
public:
  int *array;
  sequence(){}
  sequence(string header, int *ary, int length){
    len = length; head = header; array = ary;
  }
  int length(){ return len; }
  string header(){ return head;}
  int *ary(){ return array; }
  int ary(int i){ return array[i]; }
  void dump(){
    cout << head << endl << "length : " << len << endl;
    for(int i = 0; i < len; i++){ cout << array[i];}
    cout <<endl;
  }
  void destroy(){
    delete [] array;
  }
  friend sequence *seq_init(string header, int *ary, int length);
};

class job{
  hmm   *model;
  sequence *data;
  forward_backward *forward;
  forward_backward *backward;
public:
  void init(hmm *m, sequence *seq){ model = m; data = seq;}
  void destroy(){
    data->destroy();
    model->destroy();
    return;
  }
  void dump(){ /* job classの内容を表示 */
    model->dump();
    data->dump();
  }
  sequence get_data(){return *data;}
  friend int viterbi_body(job &j);
  friend forward_backward::forward_backward(const job j);
};

vector<string> split(const string &str, char delim);
istream &getline_wocomment(char c, istream &is, string &str);


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
    str.erase(comment_start);
  return is;
}

hmm *params_read(ifstream &ifs){
  hmm *model = new hmm;
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

  model->init(a_size, s_size);
  for(int i = 0; i < alphabet.length(); i++){
    if(alphabet[i] != ' '){ model->alph += alphabet[i]; }
  }

  /* 状態遷移確率を読み込む */
  for(int i = 0; i < s_size; i++){
    getline_wocomment(comment_ch, ifs, str);
    vector<string> input_str = split(str, ' ');
    long double input;
    for(int j = 0; j < s_size; j++){
      sscanf((input_str[j]).c_str(), "%Lf", &input);
      model->trans[i][j] = input;
      model->ltrans[i][j] = logl(input);
    }
  }

  /* 出力確率を読み込む 
   * s0の出力確率は0であることに注意。s1から読み込む */
  for(int j = 0; j < a_size; j++){
    model->emit[0][j]  = 0;
    model->lemit[0][j] = logl(0);

  }  
  for(int i = 1; i < s_size; i++){
    getline_wocomment(comment_ch, ifs, str);
    vector<string> input_str = split(str, ' ');
    long double input;
    for(int j = 0; j < a_size; j++){
      sscanf((input_str[j]).c_str(), "%Lf", &input);
      model->emit[i][j]  = input;
      model->lemit[i][j] = logl(input);
    }
  }
  return model;
}


std::string &data_read(string &data, ifstream &ifs){
  std::string str;
  while(getline(ifs, str)){ data += str; }
  return data;
}

sequence *seq_init(string header, int *ary, int length){
  sequence *seq = new sequence;
  seq -> len = length;  seq -> head = header;
  seq -> array = new int [length];
  for(int i = 0; i < length; i++){(seq -> array)[i] = ary[i];}
  return seq;
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

  /* hmm を構成 */
  hmm *model = params_read(param_fs);
  param_fs.close();
  
  /* データを読み込む */
  string data_str, data_head;
  getline(data_fs, data_head);
  data_str = data_read(data_str, data_fs);
  int *data =  model -> data_convert(data_str);
  /* 読み込んだ文字列をintの配列に変換 */
  data_fs.close();

  /* sequence構造体をつくる */
  sequence *seq_data = seq_init(data_head, data, data_str.length());

  /* job構造体をつくる */
  myjob.init(model, seq_data);

  return EXIT_SUCCESS;
}

int viterbi_body(job &j){
  viterbi vit;
  vit.init(j.data->length(), j.model);

  /* viterbi変数(log)の初期可 */
  vit.lv[0] = logl(1);
  for(int s = 1; s < vit.s_size; s++){
    vit.lv[s] = logl(0);
  }

  /* アルゴリズム本体を回す */
  for(int t = 0; t < vit.length; t++){
    vit.repeat(t, (j.data->ary())[t]);
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

forward_backward::forward_backward(const job j){
  len = j.data->length();
  s_size = j.model -> s_size();
  tbl = new long double * [len + 2];
  for(int t = 0; t < len + 2; t++){
    tbl[t] = new long double [s_size];
    for(int s = 0; s  < s_size ; s++){ tbl[t][s] = 0; }
  }
  scale = new long double [len + 2];
  for(int t = 0; t < len + 2; t++){ scale[t] = 0; }
  model = j.model;
}

int forward_backward::forward(const sequence data, int t_len){
  /* forward変数(\forall t, \sum_s tbl[t][s] = 1 と規格化)の初期化 */
  tbl[0][0] = 1; scale[0] = logl(1);
  for(int s = 1; s < s_size; s++){ tbl[0][s] = 0; }
  /* アルゴリズム本体を回す */
  for(int t = 1; t <= t_len; t++){
    int c = data.array[t - 1]; // 1文字目はary[0]に入っている
    for(int l = 1; l < s_size; l++){
      volatile long double sum = 0;
      for(int k = 0; k < s_size; k++){
	volatile long double trans_kl = model -> get_trans(k, l);
	//printf(" k = %Lf, trans = %Lf\n", tbl[t -1][k], trans_kl);
	sum +=  tbl[t - 1][k] * trans_kl ;
      }
      //printf("(t = %d, l = %d), sum = %Lf", t, l, sum);
      tbl[t][l] = model -> get_emit(l, c) * sum;
      //printf(", tbl[%d][%d] = %Lf\n", t, l, tbl[t][l]);
    }
    /* 規格化する */
    for(int s = 0; s < s_size; s++){ scale[t] += tbl[t][s]; }
    for(int s = 0; s < s_size; s++){ tbl[t][s] /= scale[t]; }
    /*scaling factor(log)の和を計算しておく*/
    scale[t] = logl(scale[t]);  scale[t] += scale[t - 1]; 
    printf("Forward probability (t = %d) = %Lf\n", t, scale[t]);
  }
  return 0;
}

int forward_backward::backward(const sequence data, int t_len){
  /* forward変数(\forall t, \sum_s tbl[t][s] = 1 と規格化)の初期化 */
  tbl[t_len][0] = 0; scale[t_len] = logl(1);
  for(int s = 1; s < s_size; s++){ tbl[t_len][s] = 1; }

  /* アルゴリズム本体を回す */
  for(int t = t_len - 1; t >= 0; t--){
    int c = data.array[t - 1]; // 1文字目はary[0]に入っている
    for(int k = 0; k < s_size; k++){
      volatile long double sum = 0;
      for(int l = 0; l < s_size; l++){
	volatile long double trans_kl = model -> get_trans(k, l);
	volatile long double emit_lc  = model -> get_emit(l,c);
	sum += trans_kl * emit_lc * tbl[t + 1][l];
	printf("trans_%d%d = %Lf, ", k, l, trans_kl);
	printf("emit_%d%d = %Lf, ", l, c, emit_lc);
	printf(", tbl[%d][%d] = %Lf\n", t + 1, l, tbl[t + 1][l]);
      }
      tbl[t][k] = sum;
      printf(", tbl[%d][%d] = %Lf\n", t, k, tbl[t][k]);
    }
    /* 規格化する */
    for(int s = 0; s < s_size; s++){ scale[t] += tbl[t][s]; }
    for(int s = 0; s < s_size; s++){ tbl[t][s] /= scale[t]; }
    /*scaling factor(log)の和を計算しておく*/
    scale[t] = logl(scale[t]);  scale[t] += scale[t + 1];
    printf("Backward probability (t = %d) = %Lf\n", t, scale[t]);
  }
  return 0;
}

inline long double forward_backward::calc(const int t_len){
  return scale[t_len];
}

void forward_backward::show(int t_len){
  long double lresults = this -> calc(t_len);
  cout << expl(lresults) << ", "<<lresults << endl;
  return;
}

int main(int argc, char *argv[]){
  if(argc < 2){
    cerr << "usage: $" << argv[0] << " <parameter file> <FASTA file>" << endl;
    return EXIT_FAILURE;
  }else{
    job myjob;
    prepare(myjob, argv[1], argv[2]);
    //myjob.dump();
    viterbi_body(myjob);

    forward_backward forward(myjob);
    //myjob.dump();
    forward.forward(myjob . get_data(), myjob. get_data() . length());
    forward.show(myjob . get_data() . length());

    forward_backward backward(myjob);
    backward.backward(myjob . get_data(), myjob. get_data() . length());
    backward.show(0);

    myjob.destroy();
    return 0;
  }
}
