make clean

echo "+---------- 課題1 Viterbi アルゴリズム -----------+"
g++ kadai1main.cpp
./a.out input/params.txt input/sample-RNA.fa > "kadai1_out.txt"
echo " * 実行結果は次の通りです。"
cat kadai1_out.txt

echo "+------- 課題2 前向き後ろ向き アルゴリズム -------+"
g++ kadai2main.cpp
./a.out input/params.txt input/sample-RNA.fa > "kadai2_out.txt"
echo " * 実行結果は次の通りです。"
cat kadai2_out.txt

echo "+------- 課題3 前向き後ろ向き アルゴリズム -------+"
g++ kadai3main.cpp
./a.out input/params3.txt input/trna-seq.txt > "kadai3_out.txt"
echo " * 実行結果は次の通りです。"
cat kadai3_out.txt

rm a.out
