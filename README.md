20141001 HMM の基本アルゴリズムの実装
=========

実装が終わったもの
 * Viterbiアルゴリズム
 * 前向き アルゴリズム
 * 後ろ向き アルゴリズム

今後実装予定のもの
 * FASTAフォーマットのサポート
  - 入力ファイル中に複数の配列があるときに正しく全て読み込む
  - 複数の配列を格納するためのデータ構造を追加する予定
 * Baum-Welch アルゴリズム
 * MEA アラインメント
 * viterbi, forward, backward, baum-welch, MEA, それぞれについて
   個別の実行ファイルを生成するようにする

その他
 * 過去のバージョンにCで実装したものがある
 * 現在のバージョンはC++で再度書き直したもの
 * 現在は, sample-RNA.fa の配列をテストデータとして使うことを想定している。

