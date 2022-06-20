#include<bits/stdc++.h>
using namespace std;
using F = double;
const F inf = 1e18;
pair<F,vector<F>> Simplex_Method(int N,int M,vector<F> C,vector<vector<F>> A,vector<F> B){
  //SUM(j=1,N) X[j]*C[j]を最大化する非負数列Xを構築します
  //各i(i=1,M)において SUM(j=1,N) A[i][j]*X[j] <= B[i] を満たすようにする
  //Fはdoubleだったり分数だったり

  assert(C.size() == N && A.size() == M && B.size() == M);
  for(int i=0;i<M;i++) assert(A[i].size() == N);
  F Ans = 0;//最大値
  vector<F> X;//最大値を取る数列x

  //スラック変数を導入する
  //M個の条件に対し1つずつスラック変数を用意するので、スラック変数をM個追加
  X.resize(N+M);
  C.resize(N+M);
  for(int i=0;i<M;i++) A[i].resize(N+M);

  //スラック変数を初期化 X[N+i]=B[i],A[i][N+i]=1
  for(int i=0;i<M;i++) X[N+i] = B[i];
  for(int i=0;i<M;i++) A[i][N+i] = 1;

  //スラック変数であるならば条件のindex,ないならば-1
  vector<int> Slack(N+M,-1);
  for(int i=0;i<M;i++) Slack[N+i] = i;

  //i番目の条件のスラック変数のindex
  vector<int> A_Slack(M);
  for(int i=0;i<M;i++) A_Slack[i] = N+i;

  //操作ができるまで
  while(true){
    //正の係数をとる非スラック変数を探す(最小添字規則適用)
    int target = -1;
    for(int i=0;i<N+M;i++){
      if(Slack[i] == -1 && C[i] > 0){
        target = i;break;
      }
    }

    //もしすべての係数が正でなかったら、操作を終了
    if(target == -1) break;

    //スラック変数の余裕がある限りX[target]を増やす
    //各条件で増やせる量のminを取れば良い
    //targetと、minを取る条件のスラック変数をswapするので、その条件のindexも記録する(最小添字規則適用)
    //いくらでも増やせる場合、infを返し操作を終了

    F minval = inf;
    int minitr = -1;
    for(int i=0;i<M;i++){
      if(A[i][target] <= 0) continue;
      //2つ目の条件式で最小添字規則適用(小数誤差が難しいところですが...)
      if(B[i]/A[i][target] < minval || (B[i]/A[i][target] == minval && A_Slack[i] < A_Slack[minitr])){
        minval = B[i]/A[i][target];
        minitr = i;
      }
    }
    if(minval == inf){
      Ans = inf;break;
    }

    //A[minitr]以外の条件式に対して変更を行う
    for(int i=0;i<M;i++){
      if(i == minitr || A[i][target] == 0) continue;
      F propo = A[i][target]/A[minitr][target];
      //B[i]の値を変更
      B[i] -= B[minitr]*propo;
      //スラック変数の値を変更
      X[A_Slack[i]] = B[i]/A[i][A_Slack[i]];
      //A[i]の各値を変更
      for(int j=0;j<N+M;j++){
        A[i][j] -= A[minitr][j]*propo;
      }
    }

    //最大化式も変更
    F propo = C[target]/A[minitr][target];
    Ans += B[minitr]*propo;
    for(int j=0;j<N+M;j++){
      C[j] -= A[minitr][j]*propo;
    }

    //最後にスラック変数をSwap
    X[target] = B[minitr]/A[minitr][target];
    X[A_Slack[minitr]] = 0;
    
    Slack[target] = minitr;
    Slack[A_Slack[minitr]] = -1;
    A_Slack[minitr] = target;
  }
  
  //出力
  vector<F> SubX;
  copy(X.begin(),X.begin()+N,back_inserter(SubX));
  return make_pair(Ans,SubX);
}
