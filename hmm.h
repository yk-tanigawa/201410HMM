#ifndef __HMM__
#define __HMM__

class hmm{
  /* HMM の構造体 */

public:
  long double **trans;  /* 遷移確率 */
  long double **emit;   /* 出力確率 */
  long double **ltrans; /* 遷移確率のlog */
  long double **lemit;  /* 出力確率のlog */

  int state_size;  /* 状態数 */
  int alph_size;   /* アルファベットの数 */
  char *alph;      /* 出力アルファベット記号 */
  hmm(const int a_size, const int s_size);
  hmm(){}
  int alph_to_digit(char c);
  void dump();
  hmm *init(const int a_size, const int s_size);
  void destroy();
  void set_trans(int i, int j, long double ld){
    if(i < 0 || state_size <= i || j < 0 || state_size <= j){
      perror("set_trans"); exit(EXIT_FAILURE);
    }else{
      trans[i][j] = ld;
      ltrans[i][j] = logl(ld);
    }
    return;
  }
  long double get_trans(int i, int j){
    if(i < 0 || state_size <= i || j < 0 || state_size <= j){
      perror("get_trans"); exit(EXIT_FAILURE);
    }else{
      return trans[i][j];
    }
  }
  void set_emit(int i, int j, long double ld){
    if(i < 0 || state_size <= i || j < 0 || alph_size <= j){
      perror("set_trans"); exit(EXIT_FAILURE);
    }else{
      emit[i][j] = ld;
      lemit[i][j] = logl(ld);
    }
    return;
  }
  void friend viterbi_repeat(hmm *hmm, FILE *data, 
			     long double *lv, long double *lv_before);
};

hmm *read_params(FILE *);


#endif
