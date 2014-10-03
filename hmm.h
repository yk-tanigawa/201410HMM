#ifndef __HMM__
#define __HMM__

typedef struct hmm_{
  /* HMM の構造体 */
  int state_size;  /* 状態数 */
  int alph_size;   /* アルファベットの数 */
  char *alph;      /* 出力アルファベット記号 */
  long double **trans;  /* 遷移確率 */
  long double **emit;   /* 出力確率 */
  long double **ltrans; /* 遷移確率のlog */
  long double **lemit;  /* 出力確率のlog */
} hmm;


int hmm_alph_to_digit(hmm *, char);
hmm *hmm_init(const int, const int);
void hmm_destroy(hmm *);
void hmm_dump(const hmm *);
hmm *read_params(FILE *);


#endif
