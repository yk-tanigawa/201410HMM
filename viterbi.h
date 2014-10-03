#include "hmm.h"

#ifndef __VITERBI__
#define __VITERBI__

typedef struct tracebk_{
  /* ある時刻t でのtrace back pointerを表す構造体
   * 連結リストとしてtrack backのための情報を保持 */
  int t;   /* time t */
  int *tr; /* trace back pointer 長さ(状態数)の配列*/
  struct tracebk_ *next;   /* linked list */
} tracebk;


void results_show(int *, int);
void double_ary_cpy(long double *, const long double *,const int);
void int_ary_dump(int *, int);
int max_long_double(long double *, int);

tracebk *tracebk_push(tracebk *, int *, int);
void viterbi_tracebk_dump(tracebk *tr, int state_size);
void viterbi_dump(const long double *, const int);
void viterbi_dump_expl(const long double *, const int);
void tracebk_destroy(tracebk *);

int *viterbi_traceback(tracebk *, hmm *, int);
void viterbi_repeat(hmm *, FILE *, long double *, long double *);
int viterbi_repeat_init(hmm *, FILE *);
int viterbi_file_open(const char *, const char*);
int viterbi(FILE *, FILE *);

void viterbi_repeat_debug1(int, int, long double, long double, 
			      long double, long double *);
void viterbi_repeat_debug2(hmm *, int, int, long double, long double *);
void tracebk_debug(void);

#endif

