#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <boost/algorithm/string.hpp>

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
  void replace_trans(const int i, const int j, long double ltij){
    ltrans[i][j] = ltij; trans[i][j]  = expl(ltij); return;
  }
  void replace_emit(const int i, const int c, long double leic){
    lemit[i][c] = leic; emit[i][c]  = expl(leic); return;
  }
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
  vector <sequence *> get_data(ifstream &data_fs);
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
  friend int viterbi_body(job &j, int i);
};

class forward_backward{
  long double **tbl; /*logをとった後の値を保存*/
  long double *scale;
  int len; int s_size; hmm *model;
  long double forward_chk_calc(const int t_len){
    /* forward アルゴリズムの動作の確認のため, log( P(x_1..t_len) )を計算 
     * アルゴリズムの"終了処理"に相当する計算を行う */
    return scale[t_len];
  }
  long double backward_chk_calc(const int t_len){
    /* backward アルゴリズムの動作の確認のため, log( P(x_t_len..T) )を計算
     * アルゴリズムの"終了処理"に相当する計算を行う */
    return scale[t_len - 1] + logl(tbl[t_len - 1][0]);
  }
public:
  void init(const class job, int i);
  int forward(const sequence data, int t_len);
  int backward(const sequence data, int t_len);
  void forward_chk(){
    long double lresults = this -> forward_chk_calc(len);
    cout << "P(X) = " << lresults << ", "<< expl(lresults) << " (forward)" << endl;
    return;
  }
  void backward_chk(){
    long double lresults = this -> backward_chk_calc(1);
    cout << "P(X) = " << lresults << ", "<< expl(lresults) << " (backward)" << endl;
    return;
  }
  /* 以下の2つの関数はforward, backward変数のlogを返す */
  long double lf(int t, int k){ return log(tbl[t][k]) + scale[t]; }
  long double lb(int t, int k){ return log(tbl[t][k]) + scale[t]; }
  /* P(x) を返す */
  long double lp_x(){
    return forward_chk_calc(len);
  }
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
  vector <sequence *> data;
  vector <forward_backward> *forward_all;
  vector <forward_backward> *backward_all;
public:
  void init(hmm *m, vector <sequence *> seq){
    model = m; data = seq;
  }
  void destroy(){
    for(int i = 0; i < data.size(); i++){
      (data.at(i))->destroy();
    }
    model->destroy();
    return;
  }
  void dump(){ /* job classの内容を表示 */
    model->dump();
    for(int i = 0; i < data.size(); i++){
      (data.at(i))->dump();
    }
  }
  void show_seq_head(int i){
    cout << data.at(i) -> header() << endl;
  }
  void set_forward_backward(vector <forward_backward> f,
			    vector <forward_backward> b){
    forward_all = &f;  backward_all = &b;  return;
  }
  hmm * get_model(){ return model; }
  sequence get_data(int i){return *(data.at(i));}
  int num_of_data(){ return data.size(); }
  vector <forward_backward> *get_backward_all(){
    return backward_all; }
  vector <forward_backward> *get_forward_all(){
    return forward_all; }
  forward_backward get_forward(int i){
    return (forward_all -> at(i)); }
  forward_backward get_backward(int i){
    return (backward_all -> at(i)); }
  long double lp_x(int j){ return get_forward(j) . lp_x(); }
  friend int viterbi_body(job &j, int i);
  friend void forward_backward::init(const job j, int i);
};

vector<string> split(const string &str, char delim);
istream &getline_wocomment(char c, istream &is, string &str);

void baum_welch_Mstep(job &myjob, long double **lakl, long double **lekb);
long double calc_lAkl(job &myjob, int k, int l);
long double calc_lEkb(job &myjob, int k, int b);
void baum_welch_Estep(job &myjob, long double **lakl, long double **lekb);
void baum_welch(job &myjob);


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
  int comment_start;
  while((comment_start= str.find(c)) == 0){
    getline(is, str);
  }
  if(comment_start > 0)
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
  boost::to_upper(alphabet);
  
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

sequence *seq_init(string header, int *ary, int length){
  sequence *seq = new sequence;
  seq -> len = length;  seq -> head = header;
  seq -> array = new int [length];
  for(int i = 0; i < length; i++){(seq -> array)[i] = ary[i];}
  return seq;
}

vector <sequence *> hmm::get_data(ifstream &data_fs){
  /* データファイルに含まれる配列をsequence構造体のリストとして返す */
  vector<sequence *> data;
  string head, data_str, buf;

  while(getline(data_fs, buf)){
    if(buf[0] == '>'){
      /* 新しい配列の情報が始まった */
      if(data_str.length() > 0){
	/* 一つ前の配列の情報があるとき，それをvectorに格納する */
	data.push_back(seq_init(head, 
				data_convert(data_str), 
				data_str.length()));
      }
      head.erase();  head = buf; /* 次のヘッダー */
      data_str.erase(); /* 配列情報を破棄 */
    }else{ /* 前の配列の続き */
      data_str += buf;
    }
  }
  if(data_str.length() > 0){
    /* 一つ前の配列の情報があるとき，それをvectorに格納する */
    data.push_back(seq_init(head, data_convert(data_str), 
			    data_str.length()));
  }
  return data;
}

int prepare(job &myjob, char *param_file, char *data_file){
  /* FILE stream を開き，パラメータファイルを読み込む 
   * job構造体に情報を書き込む */

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
  vector <sequence *> seq_data = model -> get_data(data_fs);
  data_fs.close();

  /* job構造体をつくる */
  myjob.init(model, seq_data);

  return EXIT_SUCCESS;
}

int viterbi_body(job &j, int i){
  viterbi vit;
  vit.init((j.data).at(i)->length(), j.model);

  /* viterbi変数(log)の初期化 */
  vit.lv[0] = logl(1);
  for(int s = 1; s < vit.s_size; s++){
    vit.lv[s] = logl(0);
  }

  /* アルゴリズム本体を回す */
  for(int t = 0; t < vit.length; t++){
    vit.repeat(t, ((j.data).at(i)->ary())[t]);
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

void forward_backward::init(job j, int i){
  len = (j.get_data(i)).length();
  s_size = j.model -> s_size();
  tbl = new long double * [len + 1];
  for(int t = 0; t < len + 1; t++){
    tbl[t] = new long double [s_size];
    for(int s = 0; s  < s_size ; s++){ tbl[t][s] = 0; }
  }
  scale = new long double [len + 1];
  for(int t = 0; t < len + 1; t++){ scale[t] = 0; }
  model = j.model;
  return;
}

int forward_backward::forward(const sequence data, int t_len){
  /* forward変数(\forall t, \sum_s tbl[t][s] = 1 と規格化)の初期化 */
  tbl[0][0] = 1; scale[0] = logl(1);
  for(int s = 1; s < s_size; s++){ tbl[0][s] = 0; }
  /* アルゴリズム本体を回す */
  for(int t = 1; t <= t_len; t++){
    int c = data.array[t - 1]; // 1文字目はary[0]に入っている
    for(int l = 1; l < s_size; l++){
      long double sum = 0;
      for(int k = 0; k < s_size; k++){
	sum +=  tbl[t - 1][k] * (model -> get_trans(k, l));
      }
      tbl[t][l] = model -> get_emit(l, c) * sum;
    }
    /* \sum_s tbl[t][s] = 1 と規格化する */
    for(int s = 0; s < s_size; s++){ scale[t] += tbl[t][s]; }
    for(int s = 0; s < s_size; s++){ tbl[t][s] /= scale[t]; }

    /* scaling factor(log)の和を計算しておく
     * scale[t] の中身 は, t = 1..t の各scalig factor の積の対数
     * 前向きアルゴリズムの終了条件から, これはP(x_1..t) の対数に等しい */
    scale[t] = logl(scale[t]);  scale[t] += scale[t - 1]; 
    //printf("P(x_1..%d) = %Lf\n", t, expl(scale[t]));
  }
  return 0;
}

int forward_backward::backward(const sequence data, int t_len){
  /* backward変数(\forall t, \sum_s tbl[t][s] = 1 と規格化)の初期化 */
  tbl[t_len][0] = 0;  scale[t_len] = logl(s_size - 1);
  for(int s = 1; s < s_size; s++){ tbl[t_len][s] = 1.0 / (s_size - 1); }

  /* アルゴリズム本体を回す */
  for(int t = t_len - 1; t >= 0; t--){
    int c = data.array[t]; // 1文字目はary[0]に入っている
    for(int k = 0; k < s_size; k++){
      long double sum = 0;
      for(int l = 0; l < s_size; l++){
	sum += expl( (model -> get_ltrans(k, l)) +
		     (model -> get_lemit(l,c)) +
		     logl(tbl[t + 1][l]) );
      }
      tbl[t][k] = sum;
    }
    /* \sum_s tbl[t][s] = 1 と規格化する */
    for(int s = 0; s < s_size; s++){ scale[t] += tbl[t][s];}
    for(int s = 0; s < s_size; s++){ tbl[t][s] /= scale[t]; }

    /* scaling factor(log)の和を計算しておく
     * scale[t] の中身 は, t = t..T の各scalig factor の積の対数
     * 後ろ向きアルゴリズムの終了条件から, 
     * P(x_1..t) は e^(scale[t]) * b_0(t) となる */
    scale[t] = logl(scale[t]);  scale[t] += scale[t + 1];
    //printf("P(x_%d..%d) = %Lf\n", t, t_len, expl(scale[t]) * tbl[t][0] );
  }
  return 0;
}

long double calc_lAkl(job &myjob, int k, int l){
  /* Baum Welch で log a_kl (transition prob.) を学習データから推定 */
  long double sum_j = 0;
  for(int j = 0; j < myjob.num_of_data(); j++){
    long double sum_t = 0;
    for(int t = 0; t < myjob.get_data(j).length() - 1 ; t++){
      long double sum_log = 0;
      sum_log += myjob.get_forward(j)  . lf(t, k);
      sum_log += myjob.get_model()->get_ltrans( k, l );
      sum_log += myjob.get_model()->get_lemit(l, myjob.get_data(j).array[t]);
      sum_log += myjob.get_backward(j) . lb(t + 1, l);
      sum_t = expl(sum_log);
    }
    sum_t -= myjob.lp_x(j);
    sum_j += sum_t;
  }
  return sum_j; /* log(A_kl) に相当 */
}


long double calc_lEkb(job &myjob, int k, int b){
  /* Baum Welch で log e_k(b) (emittion prob.) を学習データから推定 */
  long double sum_j = 0;
  for(int j = 0; j < myjob.num_of_data(); j++){
    long double sum_t = 0;
    for(int t = 0; t < myjob.get_data(j).length(); t++){
      if(myjob.get_data(j).array[t] == b){
	sum_t += expl( myjob.get_forward(j) . lf(t, k) +
		       myjob.get_backward(j). lb(t, k)  );
      }
    }
    sum_t -= myjob.lp_x(j);
    sum_j += sum_t;
  }
  return sum_j; /* log(E_k(b)) に相当 */
}

void baum_welch_Mstep(job &myjob, long double **lAkl, long double **lEkb){
  /* パラメータの最尤推定 */
  int state_size = myjob.get_model()->s_size();
  int alphabet_size = myjob.get_model()->a_size();
  for(int k = 0; k < state_size; k++){
    for(int l = 0; l < state_size; l++){
      myjob.get_model()->replace_trans(k, l,
				       lAkl[k][l] - lAkl[k][state_size]);
    }
    for(int b = 0; b < alphabet_size; b++){
      myjob.get_model()->replace_emit(k, b,
				      lEkb[k][b] - lEkb[k][alphabet_size]);
    }
  }
  return;
}

void baum_welch_Estep(job &myjob, long double **lAkl, long double **lEkb){
  /* 期待値計算 */
  int state_size = myjob.get_model()->s_size();
  int alphabet_size = myjob.get_model()->a_size();
  for(int k = 0; k < state_size; k++){
    for(int l = 0; l < state_size; l++){
      lAkl[k][l] = calc_lAkl(myjob, k, l);
      lAkl[k][state_size] = expl(lAkl[k][l]);
    }
    lAkl[k][state_size] = logl(lAkl[k][state_size]); /* 行の和 */
    for(int b = 0; b < alphabet_size; b++){
      lEkb[k][b] = calc_lEkb(myjob, k, b);
      lEkb[k][alphabet_size] = expl(lEkb[k][b]);
    }
    lEkb[k][alphabet_size] = logl(lEkb[k][alphabet_size]); /* 行の和 */
  }
  return; /* Mステップで利用するために行の和を計算した */
}

long double calc_lpx_sum(job &myjob){
  long double logsum = 0;
  for(int j = 0; j < myjob.num_of_data(); j++){
    logsum += myjob.lp_x(j);
  }
  return logsum;
}

void baum_welch(job &myjob){
  long double **lAkl = new long double * [myjob.get_model()->s_size()];
  long double **lEkb = new long double * [myjob.get_model()->s_size()];
  for(int k = 0; k < myjob.get_model()->s_size(); k++){
    lAkl[k] = new long double [myjob.get_model()->s_size() + 1];
    lEkb[k] = new long double [myjob.get_model()->a_size() + 1];
    /* (+ 1) してあるのは, 行成分の和を計算して格納するため */
  }

  /* EM algorithm */
  baum_welch_Estep(myjob, lAkl, lEkb);
  baum_welch_Mstep(myjob, lAkl, lEkb);

  for(int k = 0; k < myjob.get_model()->s_size(); k++){
    delete [] lAkl[k]; delete [] lEkb[k];
  }
  delete [] lAkl; delete [] lEkb;

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

    /* 前向き・後ろ向きの実行準備 */
    vector <forward_backward> forward_all;
    vector <forward_backward> backward_all;
    for(int i = 0; i < myjob.num_of_data(); i++){
      forward_backward forward;
      forward.init(myjob, i);
      forward_all.push_back(forward);
      
      forward_backward backward;
      backward.init(myjob, i);
      backward_all.push_back(backward);
    }

    /* jobにforward, backward へのポインタを入れる */
    myjob.set_forward_backward(forward_all, backward_all);

#if 0
    for(int i = 0; i < myjob.num_of_data(); i++){
      /* 前向きアルゴリズムの適用 */
      forward_all.at(i)  . forward(  myjob.get_data(i),
				     myjob.get_data(i).length());
      /* 後ろ向きアルゴリズムの適用 */
      backward_all.at(i) . backward( myjob.get_data(i),
				     myjob.get_data(i).length());
    }
    long double lpxsum = calc_lpx_sum(myjob);
    cout << "init :" << lpxsum << endl;
    long double lpxsum_before ;
    for(int loop = 0; loop < 10; loop++){
      
      baum_welch(myjob);

    for(int i = 0; i < myjob.num_of_data(); i++){
      /* 前向きアルゴリズムの適用 */
      forward_all.at(i)  . forward(  myjob.get_data(i),
				     myjob.get_data(i).length());
      /* 後ろ向きアルゴリズムの適用 */
      backward_all.at(i) . backward( myjob.get_data(i),
				     myjob.get_data(i).length());
    }

      lpxsum_before = lpxsum;
      lpxsum = calc_lpx_sum(myjob);
      cout << "loop (" << loop << ") :" << lpxsum << endl;
    }

#endif



    // 下記にあるのはBaum-Welchを入れる前の状態
#if 1
    for(int i = 0; i < myjob.num_of_data(); i++){
      /* 前向きアルゴリズムの適用 */
      forward_all.at(i)  . forward(  myjob.get_data(i),
				     myjob.get_data(i).length());
      //forward_all.at(i)  . forward_chk();

      /* 後ろ向きアルゴリズムの適用 */
      backward_all.at(i) . backward( myjob.get_data(i),
				     myjob.get_data(i).length());
      //backward_all.at(i) . backward_chk();
    }

    //    baum_welch(myjob);
#endif
    for(int i = 0; i < myjob.num_of_data(); i++){
      myjob.show_seq_head(i);
      /* Viterbi アルゴリズムの適用 */
      viterbi_body(myjob, i);
    }

    myjob.destroy();
    return 0;
  }
}
