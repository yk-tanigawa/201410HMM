20141001 HMM の基本アルゴリズムの実装
=========

バグ
 * 前向きアルゴリズムでセグメンテーションフォルトになる

実装が終わったもの
 * Viterbiアルゴリズム

作業途中のもの
 * 前向き アルゴリズム
 * 後ろ向き アルゴリズム
 * Baum-Welch アルゴリズム

今後実装予定のもの
 * MEA アラインメント
 * viterbi, forward, backward, baum-welch, MEA, それぞれについて
   個別の実行ファイルを生成するようにする

その他
 * 過去のバージョンにCで実装したものがある
 * 現在のバージョンはC++で再度書き直したもの
 * 現在は, sample-RNA.fa の配列をViterbiと前向き・後ろ向きのテストデータとして使うことを想定している。
 * FASTAフォーマットをサポート
  - 入力ファイル中に複数の配列があるときに全て読み込む
  - 配列データとヘッダー行の両方を保存している
  - HMM のパラメータファイルは、配列とconsistentなように作製されている必要がある
   + たとえば、配列中のあらわれる文字の数や種類を正しく表記するなど
 * 入力ファイル中の大文字と小文字を区別しない仕様
 * C++ で再度書きなおし。クラス設計の変更
 * 
