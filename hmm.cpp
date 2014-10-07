#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "hmm.h"
#include "viterbi.h"
#include <iostream>
using namespace std;

int hmm::alph_to_digit(char c){
  /* 
   * hmm の中に含まれるアルファベットの文字を受け取って、
   * 何番目の文字であるかindexを返す
   */
  for(int i = 0; i < alph_size; i++){
    if(alph[i] == c){
      return i;
    }
  }
  /* HMMのアルファベットテーブルの中に見つからない */
  fprintf(stderr, "data file conteins unknown alphabet '%c'", c);
  exit(EXIT_FAILURE);
}

hmm::hmm(const int a_size, const int s_size){
  /*
   * 構造体hmmのメモリ領域を確保する。
   * 引数には状態数と、アルファベットの数を与える
   * この関数では、callocを行うだけで、
   * パラメータのセットは別の関数で行う。
   */
#if 1
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
#endif
}

#if 1
hmm *hmm::init(const int a_size, const int s_size){
  /*
   * 本当はコンストラクタをつかって行いたいから将来的に
   * この関数は廃止される。
   */
  /*
   * 構造体hmmのメモリ領域を確保する。
   * 引数には状態数と、アルファベットの数を与える
   * この関数では、callocを行うだけで、
   * パラメータのセットは別の関数で行う。
   */
  hmm *new_hmm = new hmm();
  new_hmm -> alph   = new char [a_size];
  new_hmm -> trans   = new long double * [s_size];
  new_hmm -> ltrans  = new long double * [s_size];
  new_hmm -> emit   = new long double * [s_size];
  new_hmm -> lemit  = new long double * [s_size];
  for(int i = 0; i < s_size; i++){
    (new_hmm -> trans)[i]  = new long double [s_size];
    (new_hmm -> ltrans)[i] = new long double [s_size];
    (new_hmm -> emit)[i]  = new long double [a_size];
    (new_hmm -> lemit)[i] = new long double [a_size];
  }
  new_hmm -> state_size = s_size;
  new_hmm -> alph_size = a_size;
  return new_hmm;
}
#endif

void hmm::destroy(){
  /* HMMの構造体をfreeする */
  for(int i = 0; i < state_size; i++){
    delete [] trans[i];
    delete [] ltrans[i];
    delete [] emit[i];
    delete [] lemit[i];
  }
  delete [] trans;
  delete [] ltrans;
  delete [] emit;
  delete [] lemit;
  delete [] alph;
  delete [] trans;
  delete this;
  return;
}

void hmm::dump(){
  /* HMM構造体の内容を表示する */
  fprintf(stdout, "number of states   : %d\n", state_size);
  fprintf(stdout, "number of alphabet : %d\n", alph_size);
  fprintf(stdout, "alphabet are       : %s\n", alph);
  fprintf(stdout, "transition probability matrix is as follows:\n");
  for(int i = 0; i < state_size; i++){
    for(int j = 0; j < state_size; j++){
      fprintf(stdout, "%10Lf", trans[i][j]);
    }
    fprintf(stdout, "\n");
  }
  fprintf(stdout, "log transition probability matrix is as follows:\n");
  for(int i = 0; i < state_size; i++){
    for(int j = 0; j < state_size; j++){
      fprintf(stdout, "%10Lf", ltrans[i][j]);
    }
    fprintf(stdout, "\n");
  }
  fprintf(stdout, "emittion probabiliry matrix is as follows:\n");
  for(int i = 0; i < state_size; i++){
    for(int j = 0; j < alph_size; j++){
      fprintf(stdout, "%10Lf", emit[i][j]);
    }
    fprintf(stdout, "\n");
  }
  fprintf(stdout, "log emittion probabiliry matrix is as follows:\n");
  for(int i = 0; i < state_size; i++){
    for(int j = 0; j < alph_size; j++){
      fprintf(stdout, "%10Lf", lemit[i][j]);
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
  int alph_size, state_size; /* 読み込んだ入力を格納する変数 */

  fscanf(fp_params, "%d", &alph_size);
  fflush(stdin);

  char *alph = new char [alph_size];
  for(int i = 0; i < alph_size; i++){
    fscanf(fp_params, "%1s", &(alph[i]));
    fflush(stdin);
  }
  fscanf(fp_params, "%d", &state_size);
  fflush(stdin);

  /* アルファベット集合を構造体にコピーする */  
  /* class object をポインタとして宣言したときのconstructerの
   * 使い方がよくわかっていないから下記のようにした */
  hmm temp;
  hmm *model = temp.init(alph_size, state_size);
  for(int i = 0; i < alph_size; i++){
    (model -> alph)[i] = alph[i];
  }
  delete alph;

  /* 状態遷移確率を読み込む */
  for(int i = 0; i < state_size; i++){
    for(int j = 0; j < state_size; j++){
      long double input;
      fscanf(fp_params, "%Lf", &input);
      model -> set_trans(i, j, input);
      fflush(stdin);
    }
  }
  /* 出力確率を読み込む */
  for(int j = 0; j < alph_size; j++){
    model -> set_emit(0, j, 0);
  }
  /* s0の出力確率は0であることに注意。s1から読み込む */
  for(int i = 1; i < state_size; i++){
    for(int j = 0; j < alph_size; j++){
      long double input;
      fscanf(fp_params, "%Lf", &input);
      model -> set_emit(i, j, input);
      fflush(stdin);
    }
  }
  return model;
}
