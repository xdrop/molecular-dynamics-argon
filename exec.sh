g++ -O3 argon.cpp -o argon -std=c++11
cd remout
foldr=$(date +%Y%m%d_%H%M%S)
mkdir $foldr
cd $foldr
ln -s ../../argon
./argon
