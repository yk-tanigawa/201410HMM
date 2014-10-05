#ifndef __HMM__
#define __HMM__

class hmm{
public:
  /* HMM の構造体 */
  int state_size;  /* 状態数 */
  int alph_size;   /* アルファベットの数 */
  char *alph;      /* 出力アルファベット記号 */
  long double **trans;  /* 遷移確率 */
  long double **emit;   /* 出力確率 */
  long double **ltrans; /* 遷移確率のlog */
  long double **lemit;  /* 出力確率のlog */
  hmm(const int a_size, const int s_size);
  hmm(){}
  hmm *init(const int a_size, const int s_size);
  int alph_to_digit(char c);
  void dump();
};

void hmm_destroy(hmm *);
hmm *read_params(FILE *);


#endif
