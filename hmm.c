#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "hmm.h"
#include "viterbi.h"

inline int hmm_alph_to_digit(hmm *hmm, char c){
  /* 
   * hmm の中に含まれるアルファベットの文字を受け取って、
   * 何番目の文字であるかindexを返す
   */
  int i = 0;
  for(i = 0; i < hmm -> alph_size; i++){
    if(((hmm -> alph)[i]) == c){
      return i;
    }
  }
  /* HMMのアルファベットテーブルの中に見つからない */
  fprintf(stderr, "data file conteins unknown alphabet '%c'", c);
  exit(EXIT_FAILURE);
}

hmm *hmm_init(const int alph_size, const int state_size){
  /*
   * 構造体hmmのメモリ領域を確保する。
   * 引数には状態数と、アルファベットの数を与える
   * この関数では、callocを行うだけで、
   * パラメータのセットは別の関数で行う。
   */
  int i;
  hmm *new = calloc(1, sizeof(hmm));
  new -> alph   = calloc(alph_size,  sizeof(char));
  new -> trans  = calloc(state_size, sizeof(long double *));
  new -> ltrans  = calloc(state_size, sizeof(long double *));
  for(i = 0; i < state_size; i++){
    (new -> trans)[i] = calloc(state_size, sizeof(long double));
    (new -> ltrans)[i] = calloc(state_size, sizeof(long double));
  }
  new -> emit   = calloc(state_size,  sizeof(long double *));
  new -> lemit   = calloc(state_size,  sizeof(long double *));
  for(i = 0; i < state_size; i++){
    (new -> emit)[i]  = calloc(alph_size, sizeof(long double));
    (new -> lemit)[i]  = calloc(alph_size, sizeof(long double));
  }
  new -> state_size = state_size;
  new -> alph_size = alph_size;
  return new;
}

void hmm_destroy(hmm *hmm){
  /* HMMの構造体をfreeする */
  int i;
  for(i = 0; i < (hmm-> state_size); i++){
    free((hmm -> trans)[i]);
  }
  for(i = 0; i < (hmm-> state_size); i++){
    free((hmm -> emit)[i]);
  }
  free(hmm -> trans);
  free(hmm -> emit);
  free(hmm -> alph);
  free(hmm);
  return;
}

void hmm_dump(const hmm *data){
  /* HMM構造体の内容を表示する */
  int i, j;
  fprintf(stdout, "number of states   : %d\n", data -> state_size);
  fprintf(stdout, "number of alphabet : %d\n", data -> alph_size);
  fprintf(stdout, "alphabet are       : %s\n", data -> alph);
  fprintf(stdout, "transition probability matrix is as follows:\n");
  for(i = 0; i < data -> state_size; i++){
    for(j = 0; j < data -> state_size; j++){
      fprintf(stdout, "%10Lf", (data -> trans)[i][j]);
    }
    fprintf(stdout, "\n");
  }
  fprintf(stdout, "log transition probability matrix is as follows:\n");
  for(i = 0; i < data -> state_size; i++){
    for(j = 0; j < data -> state_size; j++){
      fprintf(stdout, "%10Lf", (data -> ltrans)[i][j]);
    }
    fprintf(stdout, "\n");
  }
  fprintf(stdout, "emittion probabiliry matrix is as follows:\n");
  for(i = 0; i < data -> state_size; i++){
    for(j = 0; j < data -> alph_size; j++){
      fprintf(stdout, "%10Lf", (data -> emit)[i][j]);
    }
    fprintf(stdout, "\n");
  }
  fprintf(stdout, "log emittion probabiliry matrix is as follows:\n");
  for(i = 0; i < data -> state_size; i++){
    for(j = 0; j < data -> alph_size; j++){
      fprintf(stdout, "%10Lf", (data -> lemit)[i][j]);
    }
    fprintf(stdout, "\n");
  }
  return;
}

hmm *read_params(FILE *fp_params){
  /*
   * FILE *fp_params が指す先のファイルから，
   * HMM のパラメータを読み込み，HMM構造体にセットして，
   * そのポインタを返す。
   */
  int i = 0, j = 0; /* for ループ用の変数 */
  int alph_size, state_size; /* 読み込んだ入力を格納する変数 */

  fscanf(fp_params, "%d", &alph_size);
    
  fflush(stdin);
  char *alph = calloc(alph_size, sizeof(char));
  for(i = 0; i < alph_size; i++){
    fscanf(fp_params, "%1s", &(alph[i]));
    fflush(stdin);
  }
  fscanf(fp_params, "%d", &state_size);
  fflush(stdin);
    
  hmm *model = hmm_init(alph_size, state_size);
  strcpy(model -> alph, alph);
  /* アルファベット集合を構造体にコピーする */
  
  /* 状態遷移確率を読み込む */
  for(i = 0; i < state_size; i++){
    for(j = 0; j < state_size; j++){
      fscanf(fp_params, "%Lf", &((model -> trans)[i][j]));
      (model -> ltrans)[i][j] = log((model -> trans)[i][j]);
      fflush(stdin);
    }
  }
    
  /* 出力確率を読み込む */
  for(j = 0; j < alph_size; j++){
    (model -> lemit)[0][j] = log(0);
  }
  /* s0の出力確率は0であることに注意。s1から読み込む */
  for(i = 1; i < state_size; i++){
    for(j = 0; j < alph_size; j++){
      fscanf(fp_params, "%Lf", &((model -> emit)[i][j]));
      (model -> lemit)[i][j] = log((model -> emit)[i][j]);
      fflush(stdin);
    }
  }
  return model;
}
