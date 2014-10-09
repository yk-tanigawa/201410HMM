make clean
echo "+---------- 課題1 Viterbi アルゴリズム -----------+"
echo "+------- 課題2 前向き後ろ向き アルゴリズム -------+"
echo " * buildして実行します。しばらくお待ちください..."
make hmm
./hmm params.txt sample-RNA.fa > "kadai_out.txt"
make clean
echo " * 実行結果は次の通りです。"
cat kadai_out.txt
echo " * 結果は kadai_out.txt にも出力されました。"
