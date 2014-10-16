#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cfloat>

using namespace std;

class job;
class hmm;
class sequence;
string str_delete_space(string);
istream &getline_wocomment(char, istream &, string &);
string str_delete_space(string);

class job {
  hmm *model; vector <sequence *> data;
public:
  void init(hmm *m, vector <sequence *> seq){
    model = m; data = seq;
  }
  void dump();
  void viterbi();
  void viterbi_dump();
  void forback_prep();
  void forward();
  void backward();
  void lpx_dump();
  long double lpx_sum();
  int BaumWelch();
  void Estep(long double **, long double **);
  void Mstep(long double **, long double **);
  long double Akl(int, int);
  long double Ekc(int, int);
  friend job *read_from_input_file(char *, char *);
};

class sequence{
  string header;  int length;  int *ary;  int *vit;
  long double **f;  long double **b;
public:
  void dump();
  void viterbi(hmm *);
  void viterbi_dump();
  void forback_prep(hmm *);
  void forward(hmm *);
  void backward(hmm *);
  long double lpx(int);
  long double lpx2(int);
  long double forward_tk(int t, int k){  return f[t][k]; }
  long double backward_tk(int t, int k){ return b[t][k]; }
  long double Akl(hmm *, int, int);
  long double Akl_sub(hmm *, int);
  long double Ekc(hmm *, int, int);
  long double Ekc_sub(hmm *, int);
#if 0
  int len()     { return length; }
  string head() { return header; }
  int x(int i)  { return ary[i]; }
#endif
  friend sequence *seq_init(string header, int *ary, int length);
};

class hmm{
  long double ** trans;   long double ** emit;
  long double ** ltrans;  long double ** lemit;
  int a_size;  int s_size;  char *alph;
  int alph_to_digit(char c){
    /* hmm の中に含まれるアルファベットの文字を受け取って、
     * 何番目の文字であるかindexを返す*/
    for(int i = 0; i < a_size; i++){
      if(alph[i] == toupper(c)){ return i; }
    }
    /* HMMのアルファベットテーブルの中に見つからない */
    cerr << "data file conteins unknown alphabet " << c << endl;
    exit(EXIT_FAILURE);
  }
public:
  void dump();
  void init(int _a_size, int _s_size);
  void init();
  vector <sequence *> get_data(ifstream &data_fs);
  int *data_convert(const string str){ /* str のデータをint *に変換 */
    int *data = new int [str.length()];
    for(int i = 0; i < str.length(); i++){
      data[i] = this -> alph_to_digit(str[i]);
    }
    return data;
  }
  int get_s_size(){ return s_size; }
  int get_a_size(){ return a_size; }
  void set_trans(int, int, long double);
  void set_emit(int, int, long double);
  friend job *read_from_input_file(char *, char *);
  friend void sequence::viterbi(hmm *);
  friend void sequence::forback_prep(hmm *);
  friend void sequence::forward(hmm *);
  friend void sequence::backward(hmm *);
  friend void job::lpx_dump();
  friend long double sequence::Akl(hmm *, int, int);
  friend long double sequence::Ekc(hmm *, int, int);
};

vector<string> split(const string &str, char delim){
  /* strを delim でsplitして vector<string> として返す */
  istringstream iss(str);
  string tmp;
  vector<string> res;
  while(getline(iss, tmp, delim)){
    res.push_back(tmp);
  }
  return res;
}

istream &getline_wocomment(char c, istream &is, string &str){
  /*  read from stream without comments, * 
   *  'c' will be treated as a delimiter */
  getline(is, str);  int comment_start;
  while((comment_start= str.find(c)) == 0){getline(is, str);}
  if(comment_start > 0){ str.erase(comment_start); }
  return is;
}

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
T ** allocat_mat(T ** matrix, int n, int m){
  /* 2D 配列にメモリを割り当てる */
  matrix = new T * [n];
  for(int i = 0; i < n; i++){
    matrix [i] = new T [m];
  }
  return matrix;
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

string str_delete_space(string buf){
  string newstr;
  for(int i = 0; i < buf.length(); i++){
    if(buf[i] != ' '){ newstr += buf[i]; }
  }
  return newstr;
}

sequence *seq_init(string header, int *ary, int length){
  sequence *seq = new sequence;
  seq -> length = length;  seq -> header = header;
  seq -> ary = new int [length];
  for(int i = 0; i < length; i++){(seq -> ary)[i] = ary[i];}
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

job *read_from_input_file(char *param_file, char *data_file){
  /* open file stream */
  ifstream param_fs(param_file);
  if (param_fs.fail()){
    cerr << "cannot open params file" << endl;
    exit(1);
  }
  ifstream data_fs(data_file);
  if (data_fs.fail()){
    cerr << "cannot open data file" << endl;
    exit(1);
  }

  /* create a hmm model from input file*/
  hmm *model = new hmm;
  string buf; char comment_ch = '%';

  getline_wocomment(comment_ch, param_fs, buf);
  sscanf(buf.c_str(), "%d", &(model -> a_size));
  
  model -> alph = new char [model -> a_size];
  getline_wocomment(comment_ch, param_fs, buf);
  for(int i = 0; i < buf.length(); i++){
    model -> alph[i] =  toupper(str_delete_space(buf) [i]);
  }

  getline_wocomment(comment_ch, param_fs, buf);
  sscanf(buf.c_str(), "%d", &(model -> s_size));

  model -> init();

  /* transittion probability */
  for(int i = 0; i < model -> s_size; i++){
    getline_wocomment(comment_ch, param_fs, buf);
    vector<string> input_str = split(buf, ' ');
    long double input;
    for(int j = 0; j < model -> s_size; j++){
      sscanf(input_str[j].c_str(), "%Lf", &input);
      model ->  trans[i][j] = input;
      model -> ltrans[i][j] = logl(input);
    }
  }

  /* emittion probability */
  for(int i = 1; i < model -> s_size; i++){
    getline_wocomment(comment_ch, param_fs, buf);
    vector<string> input_str = split(buf, ' ');
    long double input;
    for(int j = 0; j < model -> a_size; j++){
      sscanf(input_str[j].c_str(), "%Lf", &input);
      model ->  emit[i][j] = input;
      model -> lemit[i][j] = logl(input);
    }
  }
  for(int j = 0; j < model -> a_size; j++){
    model ->  emit[0][j] = 0;
    model -> lemit[0][j] = logl(0);
  }
  job *jb = new job;
  jb -> init(model, model -> get_data(data_fs));
  return jb;
}

void job::dump(){
  model -> dump();
  for(int i = 0; i < data.size(); i++){ (data.at(i))->dump(); }
}

void hmm::dump(){
      /* HMM構造体の内容を表示する */
  cout << "number of alphabet :" << a_size << endl;
  cout << "alphabet are       :" << alph << endl;
  cout << "number of states   :" << s_size << endl;
  cout << "transition probability matrix is as follows: "
       << endl;
  show_matrix(trans,  "%10Lf", s_size, s_size);
  cout << "log transition probability matrix is as follows: "
       << endl;
  show_matrix(ltrans, "%10Lf", s_size, s_size);
  cout << "emittion probabiliry matrix is as follows: " 
       << endl;
  show_matrix(emit,   "%10Lf", s_size, a_size);
  cout << "log emittion probabiliry matrix is as follows: " 
       << endl;
  show_matrix(lemit,  "%10Lf", s_size, a_size);
  return;
}

void sequence::dump(){
  cout << header << endl << "length : " << length << endl;
  for(int i = 0; i < length; i++){ cout << ary[i];}
  cout << endl;
}

void hmm::init(int _a_size, int _s_size){
  a_size = _a_size;  s_size = _s_size;
  return this -> init();
}

void hmm::init(){
  trans  = allocat_mat(trans, s_size, s_size);
  ltrans = allocat_mat(trans, s_size, s_size);
  emit   = allocat_mat(trans, s_size, a_size);
  lemit  = allocat_mat(trans, s_size, a_size);
  for(int j = 0; j < a_size; j++){
    emit[0][j] = 0;  lemit[0][j] = logl(0);
  }
  return;
}

void job::viterbi(){
  for(int i = 0; i < data.size(); i++){ (data.at(i))->viterbi(model); }
}

void job::forback_prep(){
  for(int i = 0; i < data.size(); i++){ (data.at(i))->forback_prep(model); }
}

void job::forward(){
  for(int i = 0; i < data.size(); i++){ (data.at(i))->forward(model); }
}

void job::backward(){
  for(int i = 0; i < data.size(); i++){ (data.at(i))->backward(model); }
}

void job::viterbi_dump(){
  for(int i = 0; i < data.size(); i++){
    (data.at(i))->viterbi_dump(); 
  }
}

void job::lpx_dump(){
  /* 前向き・後ろ向きアルゴリズムの実装が正しいかを確認するため全確率を計算 */
  for(int i = 0; i < data.size(); i++){
    cout << " [forward]  " << "log( P(x^" <<i << ") ) = " 
	 << (data.at(i))->lpx(model -> s_size) << endl;
    cout << " [backward] " << "log( P(x^" <<i << ") ) = " 
	 << (data.at(i))->lpx2(model -> s_size) << endl;
  }
}

long double job::lpx_sum(){
  long double sum = 0;
  for(int j = 0; j < data.size(); j++){
    sum += (data.at(j)) -> lpx(model -> get_s_size()); 
  }
  return sum;
}

void sequence::viterbi_dump(){
  cout << header << endl;
  for(int t = 0; t < length; t++){ cout << vit[t]; }
  cout << endl; return;
}

void sequence::viterbi(hmm *model){
  vit = new int [length];
  long double *lv1, *lv2, **trbk;
  lv1 = new long double [length];
  lv2 = new long double [length];
  trbk = allocat_mat(trbk, length, model -> s_size);
  
  /* Viterbi変数の初期化 */
  lv1[0] = logl(1);
  for(int s = 1; s < model -> s_size; s++){ lv1[s] = logl(0); }

  /* Viterbiアルゴリズムの繰り返しステップ */
  long double *lv_before, *lv; /* viterbi 変数は表2行を交互に使う*/
  for(int t = 0; t < length; t++){
    int c = ary[t];
    if(t % 2 == 0){ lv_before = lv1; lv = lv2; }
    else          { lv_before = lv2; lv = lv1; }
    for(int l = 0; l < model -> s_size; l++){ /* 状態が k => l に遷移した */
      long double max = -1 * INFINITY; int max_index = -1;
      for(int k = 0; k < model -> s_size; k++){
	long double ltrans_kl = model -> ltrans[k][l];
	long double temp = lv_before[k] + ltrans_kl;
	if(temp > max){ max = temp; max_index = k; }
      }
      lv[l] = model->lemit[l][c] + max; /* viterbi変数を計算 */
      trbk[t][l] = max_index; /* trace back ポインタをセット */
    }
  }
  /* Viterbiアルゴリズムのトレースバック */
  vit[length - 1] = find_max_index(lv, model -> s_size);
  for(int t = length - 1; t > 0; t--){
    vit[t - 1] = trbk[t][vit[t]];
  }
  delete [] lv1;  delete [] lv2;
}

void sequence::forback_prep(hmm *model){
  /* 前向き後ろ向きの準備を行う。 メモリを割り当てて，値をセットする。 */

  /* 前向きアルゴリズムの準備 */
  f = allocat_mat(f, length + 1, model -> s_size + 1);
  for(int t = 0; t <= length; t++){ f[t][model -> s_size] = 0; }
  /* 時刻tでのscaling factor は f[t][model -> s_size]に格納する */
  f[0][0] = 1; f[0][model -> s_size] = 1;
  for(int s = 1; s < model -> s_size; s++){ f[0][s] = 0;}
  
  /* 後ろ向きアルゴリズムの準備 */
  b = allocat_mat(f, length + 1, model -> s_size + 1);
  for(int t = 0; t <= length; t++){ b[t][model -> s_size] = 0; }
  /* 時刻tでのscaling factor は f[t][model -> s_size]に格納する */
  b[length][0] = 0; b[length][model -> s_size] = 1;
  for(int s = 1; s < model -> s_size; s++){ 
    b[length][s] = 1.0;
  }
}

void sequence::forward(hmm *model){
  for(int t = 1; t <= length; t++){
    int c = ary[t - 1];
    /* まずは愚直に計算する */
    for(int l = 1; l < model -> get_s_size(); l++){
      long double sum = 0;
      for(int k = 0; k < model -> get_s_size(); k++){
	sum += f[t - 1][k] * model -> trans[k][l];
      }
      f[t][l] = model -> emit[l][c] * sum;
      //cout << "f[" << t << "][" << l << "] = " << f[t][l] << endl;
    }

    /* スケーリングを行う \sum_s f[t][s] = 1 と規格化 */
    for(int s = 0; s < model -> get_s_size(); s++){ 
      f[t][model -> get_s_size()] += f[t][s];
    }
#if 0 /* 表が正しく計算されているか表示して Debug する */
    printf(" t =%3d : ", t);
    for(int s = 0; s < model -> get_s_size() + 1; s++){ 
      printf("  %Lf", f[t][s]);
    }
    printf("\n");
#endif
    for(int s = 0; s < model -> get_s_size(); s++){ 
      f[t][s] /= f[t][model -> get_s_size()];
    }
  }
  return;
}

void sequence::backward(hmm *model){
  for(int t = length - 1; t >= 0; t--){
    int c = ary[t]; /* (t + 1) - 1 */
    for(int k = 0; k < model -> get_s_size(); k++){
      b[t][k] = 0;
      for(int l = 1; l < model -> get_s_size(); l++){
	b[t][k] += model -> trans[k][l] * model -> emit[l][c] * b[t + 1][l];
      }
    }
    /* スケーリングを行う \sum_s f[t][s] = 1 と規格化 */
    for(int s = 0; s < model -> get_s_size(); s++){ 
      b[t][model -> get_s_size()] += b[t][s];
    }
#if 0 /* 表が正しく計算されているか表示して Debug する */
    printf(" t =%3d : ", t);
    for(int s = 0; s < model -> get_s_size() + 1; s++){ 
      printf("  %Lf", b[t][s]);
    }
    printf("\n");
#endif
    for(int s = 0; s < model -> get_s_size(); s++){ 
      b[t][s] /= b[t][model -> get_s_size()];
    }
  }
  return;
}


inline long double sequence::lpx(int s_size){
  /* 前向きアルゴリズムの結果を用いて計算する */
  long double lpx = 0;   /* log ( P(x) ) */
  for(int t = 0; t <= length; t++){ lpx += logl(f[t][s_size]); }
  return lpx;
}

inline long double sequence::lpx2(int s_size){
  /* 後ろ向きアルゴリズムの結果を用いて計算する */
  long double lpx = logl(b[0][0]);
  for(int t = length - 1; t >= 0; t--){
    lpx += logl(b[t][s_size]); 
  }
  return lpx;
}

int job::BaumWelch(){
  long double **A = 
    allocat_mat(A, model -> get_s_size(), model -> get_s_size() + 1);
  long double **E = 
    allocat_mat(E, model -> get_s_size(), model -> get_a_size() + 1);

  long double lpx = -1 * DBL_MAX, lpx_before = 0;

  forward(); backward();
  int t = 0;
  while(1){
    Estep(A, E); Mstep(A, E); t++;
    forward(); backward(); lpx_before = lpx; lpx = lpx_sum();
    //cout << lpx_before << "\t" << lpx - lpx_before << endl;
    if(lpx - lpx_before < 0.5 || t > 40) break;
  }
  return t;
}

inline void hmm::set_trans(int k, int l, long double Akl){
  trans[k][l] = Akl; ltrans[k][l] = logl(Akl); return;
}

inline void hmm::set_emit(int k, int c, long double Ekc){
  emit[k][c] = Ekc;  lemit[k][c] = logl(Ekc);  return;
}


void job::Mstep(long double **A, long double **E){
  for(int k = 0; k < model -> get_s_size(); k++){
    for(int l = 0; l < model -> get_s_size(); l++){
      model -> set_trans(k, l, A[k][l] / A[k][model -> get_s_size()]);
    }
    for(int c = 0; c < model -> get_a_size(); c++){
      model -> set_emit(k, c, E[k][c] / E[k][model -> get_a_size()]);
    }
  }
}

void job::Estep(long double **A, long double **E){
  /* Aklの表とEkcの表を埋める*/
  for(int k = 0; k < model -> get_s_size(); k++){
    A[k][model -> get_s_size()] = 0; /* 行kのlに関する和をここに入れる */
    for(int l = 0; l < model -> get_s_size(); l++){
      A[k][l] = Akl(k, l); A[k][model -> get_s_size()] += A[k][l];
    }
    E[k][model -> get_a_size()] = 0; /* 行kのcに関する和をここに入れる */
    for(int c = 0; c < model -> get_a_size(); c++){
      E[k][c] = Ekc(k, c); E[k][model -> get_a_size()] += E[k][c];
    }
  }
  return;
}

inline long double job::Akl(int k, int l){
  long double sum = 0;
  for(int j = 0; j < data.size(); j++){
    sum += (data.at(j))->Akl(model, k, l);
  }
  return sum;
}

inline long double sequence::Akl_sub(hmm *model, int t){
  long double sum = 0;
  for(int k = 0; k < model -> get_s_size(); k++){
    sum += f[t][k] * b[t][k];
  }
  return sum;
}

inline long double sequence::Akl(hmm *model, int k, int l){
  long double sum = 0;
  for(int t = 0; t < length - 1; t++){
    sum += f[t][k] * b[t + 1][l] 
      / b[t][model -> get_s_size()] 
      * (model -> trans[k][l]) 
      * (model -> emit[l][ary[t - 1]])
      / Akl_sub(model, t);
  }
  return sum;
}

inline long double job::Ekc(int k, int c){
  long double sum = 0;
  for(int j = 0; j < data.size(); j++){
    sum += (data.at(j))->Ekc(model, k, c);
  }
  return sum;
}

inline long double sequence::Ekc_sub(hmm *model, int t){
  return Akl_sub(model, t);
}

inline long double sequence::Ekc(hmm *model, int k, int c){
  long double sum_cond = 0;
  for(int t = 0; t < length - 1; t++){
    long double temp = Ekc_sub(model, t);
    if( ary[t - 1] == c ){
      sum_cond += f[t][k] * b[t][k] / temp;
    }
  }
  return sum_cond;
}

#if 0
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
#endif

