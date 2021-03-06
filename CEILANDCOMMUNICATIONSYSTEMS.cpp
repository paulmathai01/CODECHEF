#include<stdio.h>
#include<ctime>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#define REP(i,a,b) for(i=a;i<b;i++)
#define rep(i,n) REP(i,0,n)
 
#define ll long long
 
ll ab(ll x){ if(x < 0) return -x; return x; }
 
int intArrayBinarySearchSmallElement(int d[],int size,int n){
int c = size/2; int i;
if(size==1){if(d[0]<n) return 0; return -1;}
if(size<=0) return -1;
if(d[c] < n) return intArrayBinarySearchSmallElement(d+c,size-c,n) + c;
return intArrayBinarySearchSmallElement(d,c,n);
}
 
int N, Q, M;
int X[100100]; ll sumX[100100];
int U[51][51];
int towerY[100000], towerR[100000], towerC[100000], towerT[100000], towerSize;
 
int treeNode, dataSize;
int treeLeft[500000], treeRight[500000];
char treeFixed[500000];
ll treeSum[500000]; int treeLT[500000], treeRT[500000];
 
void treeRecalc(int node){
int n1 = node * 2 + 1, n2 = node * 2 + 2;
int sz = treeRight[node] - treeLeft[node];
int t;
 
if(sz == 1) return;
if(treeFixed[node] == 0) return;
 
treeFixed[n1] = treeFixed[n2] = 1;
treeFixed[node] = 0;
 
t = treeLT[node];
treeLT[n1] = treeRT[n1] = treeLT[n2] = treeRT[n2] = t;
treeSum[n1] = 2 * towerC[t] * ((sumX[treeRight[n1]] - sumX[treeLeft[n1]]) - (ll)towerY[t] * (sz/2));
if(treeSum[n1] < 0) treeSum[n1] *= -1;
treeSum[n2] = 2 * towerC[t] * ((sumX[treeRight[n2]] - sumX[treeLeft[n2]]) - (ll)towerY[t] * (sz/2));
if(treeSum[n2] < 0) treeSum[n2] *= -1;
}
 
void treeGetSum(int A, int B, ll *sum, int *LT, int *RT, int node){
int n1 = node * 2 + 1, n2 = node * 2 + 2;
int sz = treeRight[node] - treeLeft[node];
ll sum1, sum2; int LT1, LT2, RT1, RT2;
 
if(A < treeLeft[node]) A = treeLeft[node];
if(B > treeRight[node]) B = treeRight[node];
 
if(A == treeLeft[node] && B == treeRight[node]){
*sum = treeSum[node];
*LT = treeLT[node];
*RT = treeRT[node];
return;
}
 
treeRecalc(node);
 
if(B <= treeRight[n1]){
treeGetSum(A, B, sum, LT, RT, n1);
return;
}
 
if(treeLeft[n2] <= A){
treeGetSum(A, B, sum, LT, RT, n2);
return;
}
 
treeGetSum(A, B, &sum1, &LT1, &RT1, n1);
treeGetSum(A, B, &sum2, &LT2, &RT2, n2);
if(sum1==-1 || sum2==-1){ *sum = -1; return; }
*sum = sum1 + sum2;
*LT = LT1;
*RT = RT2;
if(LT2 != RT1) *sum += U[towerT[LT2]][towerT[RT1]];
}
 
void treeSetTower(int A, int B, int t, int node){
int n1 = node * 2 + 1, n2 = node * 2 + 2;
int sz = treeRight[node] - treeLeft[node];
ll sum1, sum2; int LT1, LT2, RT1, RT2;
 
if(A < treeLeft[node]) A = treeLeft[node];
if(B > treeRight[node]) B = treeRight[node];
if(A >= B) return;
 
if(A == treeLeft[node] && B == treeRight[node]){
treeFixed[node] = 1;
treeSum[node] = 2 * towerC[t] * ((sumX[B] - sumX[A]) - (ll)towerY[t] * sz);
if(treeSum[node] < 0) treeSum[node] *= -1;
treeLT[node] = treeRT[node] = t;
return;
}
 
treeRecalc(node);
 
treeSetTower(A, B, t, n1);
treeSetTower(A, B, t, n2);
sum1 = treeSum[n1]; sum2 = treeSum[n2];
LT1 = treeLT[n1]; LT2 = treeLT[n2];
RT1 = treeRT[n1]; RT2 = treeRT[n2];
 
if(sum1==-1 || sum2==-1){ treeSum[node] = -1; return; }
treeSum[node] = sum1 + sum2;
treeLT[node] = LT1;
treeRT[node] = RT2;
if(LT2 != RT1) treeSum[node] += U[towerT[LT2]][towerT[RT1]];
}
 
 
int main(){
double start = clock();
int i,j,k,l;
int A, B, Y, R, C, T, t;
int t1, t2;
int mode;
ll res; int LT, RT;
 
scanf("%d%d%d",&N,&M,&Q);
rep(i,N) scanf("%d",X+i);
rep(i,M) rep(j,M) scanf("%d",U[i]+j);
 
sumX[0] = 0;
rep(i,N) sumX[i+1] = sumX[i] + X[i];
towerSize = 0;
 
dataSize = 1;
while(dataSize < N) dataSize *= 2;
treeNode = dataSize * 2 - 1;
rep(i,treeNode) treeFixed[i] = 0, treeSum[i] = treeLT[i] = treeRT[i] = -1;
treeLeft[0] = 0, treeRight[0] = dataSize;
rep(i,treeNode/2){
int n1 = i * 2 + 1, n2 = i * 2 + 2;
int sz = treeRight[i] - treeLeft[i];
treeLeft[n1] = treeLeft[i];
treeRight[n1] = treeLeft[i] + sz / 2;
treeLeft[n2] = treeLeft[i] + sz / 2;
treeRight[n2] = treeRight[i];
}
 
while(Q--){
scanf("%d",&mode);
if(mode == 1){
scanf("%d%d",&A,&B); A--; B--;
treeGetSum(A, B+1, &res, &LT, &RT, 0);
if(res < 0){
puts("impossible");
} else {
res -= ab(X[A] - towerY[LT]) * towerC[LT];
res -= ab(X[B] - towerY[RT]) * towerC[RT];
printf("%lld\n",res);
}
} else {
scanf("%d%d%d%d",&Y,&R,&C,&T); T--;
t = towerSize++;
towerY[t] = Y; towerR[t] = R; towerC[t] = C; towerT[t] = T;
 
A = intArrayBinarySearchSmallElement(X, N, Y - R) + 1;
B = intArrayBinarySearchSmallElement(X, N, Y + 1);
if(A <= B) treeSetTower(A, B+1, t, 0);
 
A = intArrayBinarySearchSmallElement(X, N, Y + 1) + 1;
B = intArrayBinarySearchSmallElement(X, N, Y + R + 1);
if(A <= B) treeSetTower(A, B+1, t, 0);
}
}
 
if ((double)(clock()-start)/CLOCKS_PER_SEC > 2.0) puts("Time limit exceeded");
 
return 0;
} 


