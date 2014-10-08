make clean
echo "+---------- 課題1 Viterbi アルゴリズム ----------+"
echo " * buildして実行します。しばらくお待ちください..."
make hmm
./hmm params.txt sample-RNA.fa > "kadai1_out.txt"
make clean
echo " * 実行結果は次の通りです。"
cat kadai1_out.txt
echo " * 結果は kadai1_out.txt にも出力されました。"
