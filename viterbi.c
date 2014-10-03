#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "hmm.h"
#include "viterbi.h"

#define BUF_SIZE 256

void double_ary_cpy(long double *target, const long double *source,
		    const int length){
  /* 配列のコピーを行う:
   * target に sourceの内容を書き込む */
  int i;
  for(i = 0; i < length; i++){
    target[i] = source[i];
  }
  return;
}

void results_show(int *ary, int length){
  /* viterbiアルゴリズムの結果を表示する。 */
  int i;
  for(i = 0; i < length; i++){
    fprintf(stdout, "%d", ary[i]);
  }
  fprintf(stdout, "\n");
  fflush(stdout);
  return; 
}

void int_ary_dump(int *ary, int length){
  /* int型の長さlengthの配列の内容を表示する。*/
  int i;
  //puts("int_ary_dump");
  for(i = 0; i < length; i++){
    fprintf(stdout, "%3d", ary[i]);
  }
  fprintf(stdout, "\n");
  fflush(stdout);
  return; 
}

int max_index_long_double(long double *ary, int length){
  /* long double型の配列*ary（配列長length）の中で
   * 最大要素のindexを返す*/
  int i, n = 0;
  long double max = ary[0];
  for(i = 1; i < length; i++){
    if(ary[i] > max){ max = ary[i]; n = i; }
  }
  return n;
}


tracebk *tracebk_push(tracebk *list, int *ary, int length){
  /* trace back pointerをスタックに積む 
   * スタックは連結リストとして実装されている */
  int i;
  tracebk *new = calloc(1, sizeof(tracebk));
  new -> tr = calloc(length, sizeof(int));
  for(i = 0; i < length; i++){
    (new -> tr)[i] = ary[i];
  }
  new -> next = list;
  if(list != NULL){
    (new -> t) = ((list -> t) + 1);
  }else{
    (new -> t) = 1;
  }
  //viterbi_tracebk_dump(new, length);
  return new;
}


void viterbi_tracebk_dump(tracebk *tr, int state_size){
  /* tracebk *の連結リストの内容を表示 */
  if(tr == NULL){
    return;
  }else{
    //puts("viterbi_tracebk_dump >");
    printf("%2d: ", tr -> t);
    int_ary_dump(tr -> tr, state_size);
    viterbi_tracebk_dump(tr -> next, state_size);
    return;
  }
}

void viterbi_dump(const long double *ary, const int length){
  int i;
  for(i = 0; i < length; i++){
    fprintf(stdout, "%10Lf", ary[i]);
  }
  fprintf(stdout, "\n");
  fflush(stdout);
  return;
}

void viterbi_dump_expl(const long double *ary, const int length){
  int i;
  for(i = 0; i < length; i++){
    fprintf(stdout, "%10Lf", expl(ary[i]));
  }
  fprintf(stdout, "\n");
  fflush(stdout);
  return;
}

void viterbi_repeat(hmm *hmm, FILE *data, 
		       long double *lv, long double *lv_before){
  /* viterbiの繰り返しアルゴリズムの本体 */

  /* trace back pointerの宣言 */
  tracebk *tr = NULL;
  char c;
  int l = 0, k = 0, c_index;
  while(fscanf( data , "%1s" , &c ) != EOF){
    //puts("start of loop");
    //viterbi_tracebk_dump(tr, hmm ->  state_size);
    if(c != '\n'){
      /* Viterbi変数はひとつ前のものだけわかればよい
       * よって，lv_beforeにlvを上書きする
       * lv_beforeに入っていた内容は忘れる */
      double_ary_cpy(lv_before, lv, hmm -> state_size);
      /* 文字cを数値に変換 */
      c_index = hmm_alph_to_digit(hmm, c);
      /* traceback用の配列を確保*/
      int *tracebk = calloc(hmm -> state_size, sizeof(int));
      //printf("tracebk[] allocated: ");
      //int_ary_dump(tracebk, hmm -> state_size);
      for(l = 0; l < (hmm -> state_size); l++){
	/* l ... 時刻Tでありえる全ての隠れ状態についてのループ */
	long double max = -1 * INFINITY;
	int max_index = -1;
	for(k = 0; k < (hmm -> state_size); k++){
	  /* k ... 時刻 T - 1でありえる全ての隠れ状態について */
	  long double ltrans_kl = (hmm -> ltrans)[k][l];
	  long double temp = lv_before[k] + ltrans_kl;
	  if(temp > max){max = temp;  max_index = k;}
          /* viterbi_repeat_debug1(k, l, temp, max,
	     ltrans_kl, lv_before);*/
	}
	/* Viterbi変数(log) を計算して，トレースバックをセット*/
	lv[l] = (hmm -> lemit)[l][c_index] + max;
	tracebk[l] = max_index;
        /* viterbi_repeat_debug2(hmm, l, c_index, max, lv); */
      }     
#if 0
      printf("tracebk[] computed : ");
      int_ary_dump(tracebk, hmm -> state_size);
#endif
      tr = tracebk_push(tr, tracebk, hmm -> state_size);
      free(tracebk);
#if 0
      fprintf(stdout, "%c", c);      
      viterbi_dump_expl(lv, hmm -> state_size);
      viterbi_tracebk_dump(tr, hmm ->  state_size);
      puts("end of loop");
#endif
    }
  }

  /* tracebackして結果を表示する */
  int *results = 
    viterbi_traceback(tr, hmm, max_index_long_double(lv, hmm -> state_size));
  results_show(results, tr -> t);
  free(results);
  return;
}

int *viterbi_traceback(tracebk *tr, hmm *hmm,
		       int tracebk_start){
  /* viterbiのトレースバックを行い，結果を入れた配列を返す */
  //  printf("length :%d tracebk from %d \n", length, tracebk_start);
  int i, length = tr -> t;
  int *results = calloc(length, sizeof(int));
  for(i = length; i > 0; i--){
    results[i - 1] = tracebk_start;
    //printf("%d: %d\n", i, tracebk_start);
    tracebk_start  = (tr -> tr)[tracebk_start];
    tr = tr -> next;
  }
  return results;
}

int viterbi_repeat_init(hmm *hmm, FILE *data){
  /* Viterbiの繰り返しステップ開始処理を行い，
   * ループを行う関数に投げる */

  /* FASTA のヘッダー行を読み込む */
  char first_line[BUF_SIZE];
  fscanf( data, "%s", first_line);

  /* Viterbi変数を宣言 */
  long double *lv_before = calloc((hmm -> state_size), sizeof(long double));
  /* lv_before: t = t - 1 での Viterbi変数 の log をとったもの*/
  long double *lv        = calloc((hmm -> state_size), sizeof(long double));
  /* lv       : t = t     での Viterbi変数 の log をとったもの*/
  lv[0] = logl(1);
  int i;
  for(i = 1; i < hmm -> state_size; i++){
    lv[i] = logl(0);
  }

  viterbi_repeat(hmm, data, lv, lv_before);

  free(lv_before);
  free(lv);
  return EXIT_SUCCESS;
}

inline int viterbi_file_open(const char *params, const char *data){
  /* HMM のパラメータファイルと，データファイルの名前を
   * 受け取ってviterbiアルゴリズムを適用する。
   * ファイルが存在しない場合は中止 */

  FILE *fp_params, *fp_data;
  if((fp_params = fopen(params, "r")) == NULL){
    perror("cannot open params file");
    exit(EXIT_FAILURE);
  }else{
    /* HMMのパラメータファイルは存在 */
    if((fp_data = fopen(data, "r")) == NULL){
      perror("cannot open data file");
      exit(EXIT_FAILURE);
    }else{
      /* データのファイルも存在 */
      viterbi(fp_params, fp_data);
      fclose(fp_data);
      fclose(fp_params);
      return EXIT_SUCCESS;
    }
  }
}

int viterbi(FILE *params, FILE *data){
  /* 入力を受け取ってviterbiを回す */
  hmm *hmm = read_params(params);
  //hmm_dump(hmm);
  viterbi_repeat_init(hmm, data);
  hmm_destroy(hmm);
  return EXIT_SUCCESS;
}  

/*
 * 以下debug用関数
 */

void viterbi_repeat_debug1(int k, int l,
			     long double temp, long double max,
			     long double ltrans_kl, long double *lv_before){
  fprintf(stdout, "t(%d -> %d) = %Lf ",k, l, expl(ltrans_kl));
  fprintf(stdout, "v_bef[k] = %Lf ",expl(lv_before[k]));
  fprintf(stdout, "temp = %Lf, max = %Lf\n",expl(temp), expl(max));
  return;
}

void viterbi_repeat_debug2(hmm *hmm, int l, int c_index,
			      long double max, long double *lv){
  fprintf(stdout, "l = %d ", l);
  fprintf(stdout, "emittion prob = %Lf ", expl((hmm -> lemit)[l][c_index]));
  fprintf(stdout, "v = %Lf\n",
	  expl((hmm -> lemit)[l][c_index] + max));
  fprintf(stdout, "v[l] = %Lf\n", expl(lv[l]));
}

void tracebk_debug(void){
  int i;
  tracebk *temp = NULL;  
  for(i = 0; i < 4; i++){
    int *ary = calloc(4, sizeof(int));
    ary[0] = i;
    temp = tracebk_push(temp, ary, 4);
    free(ary);
  }
  return;
}


int main(int argc, char *argv[]){
  //tracebk_debug();
  if(argc < 2){
    fprintf(stderr, "usage: $%s <parameter file> <FASTA file>\n", argv[0]);
    return EXIT_FAILURE;
    //    return viterbi_file_open("params.txt", "sample-RNA.fa");
  }else{
    return viterbi_file_open(argv[1], argv[2]);
  }
}
