#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <algorithm>
#define ll long long
using namespace std;
 
const ll Inf = 100000000, Lim = 1000000;
const int MaxN = 210, MaxM = MaxN * MaxN * 6;
int P[MaxN], V[MaxN][MaxN];
int N, M, Rx[MaxM], Ry[MaxM], Rl[MaxM], Rs[MaxM], Rt[MaxM];

class MinCostFlow_Template {
public:
	int En[MaxN], Next[MaxM], Point[MaxM], F[MaxM], Tot, S, T, H[MaxM * 10], Pre[MaxN], Cnt[MaxN];
	ll V[MaxM], Dist[MaxN];
	bool In[MaxN];
	ll Spfa() {
		memset(Dist, 63, sizeof(Dist));
		memset(Cnt, 0, sizeof(Cnt));
		int L = 0, R = 1;
		H[R] = S;
		In[S] = 1;
		Dist[S] = 0;
		while (L++ != R) {
			int U = H[L];
			if (Cnt[U]++ >= N * 2) return -Inf * Inf;
			for (int i = En[U]; i; i = Next[i])
				if (F[i] && (Dist[Point[i]] > Dist[U] + V[i])) {
					Dist[Point[i]] = Dist[U] + V[i];
					Pre[Point[i]] = i;
					if (!In[Point[i]]) {
						In[Point[i]] = 1;
						if ((Dist[Point[i]] > Dist[U]) || !L) H[++R] = Point[i];
						else H[L--] = Point[i];
					}
				}
			In[U] = 0;
		}
		return Dist[T];
	}
public:
	void Init(int _S, int _T) {
		memset(En, 0, sizeof(En));
		memset(In, 0, sizeof(In));
		Tot = 1;
		S = _S;	T = _T;
	}
	void Add(int X, int Y, int Flow, ll Cost) {
		Next[++Tot] = En[X];	En[X] = Tot;	Point[Tot] = Y;	F[Tot] = Flow;	V[Tot] = Cost;
		Next[++Tot] = En[Y];	En[Y] = Tot;	Point[Tot] = X;	F[Tot] = 0;	V[Tot] = -Cost;
	}
	ll MinCostFlow() {
		ll Ret = 0, Tmp;
		while ((Tmp = Spfa()) < Inf) {
			if (Tmp == -Inf * Inf) return -Inf * Inf;
			int Flow = Inf;
			for (int X = T; X != S; X = Point[Pre[X] ^ 1])
				Flow = min(Flow, F[Pre[X]]);
			Ret += Tmp * Flow;
			for (int X = T; X != S; X = Point[Pre[X] ^ 1]) {
				F[Pre[X]] -= Flow;
				F[Pre[X] ^ 1] += Flow;
			}
		}
		return Ret;
	}
}	G;
 
void GetAnswer() {
	for (int i = 1; i <= N * 2; i++) V[i][0] = 0, V[0][i] = Lim;
	for (int i = 1; i <= N * 2; i++)
		for (int j = G.En[i]; j; j = G.Next[j])
			if (!(j & 1) && (G.Point[j] <= N * 2) && (G.F[j ^ 1])) {
				//printf("%d %d %d\n", i, G.Point[j], int(G.V[j]));
				V[G.Point[j]][i] = G.V[j];
				V[i][G.Point[j]] = -G.V[j];
			}
	memset(P, 63, sizeof(P));

	static bool In[MaxN];
	static int Q[MaxM * 10];
	int L = 0, R = 1;
	memset(In, 0, sizeof(In));
	Q[R] = 0;
	In[0] = 1;
	P[0] = 0;
	while (L++ != R) {
		for (int i = 0; i <= N * 2; i++)
			if (P[i] > P[Q[L]] + V[Q[L]][i]) {
				P[i] = P[Q[L]] + V[Q[L]][i];
				if (!In[i]) Q[++R] = i;
			}
		In[Q[L]] = 0;
	}
}
 
void Solve() {
	memset(V, 63, sizeof(V));
	int Deg0[MaxN], Deg1[MaxN];
	ll Sum = 0;
	memset(Deg0, 0, sizeof(Deg0));
	memset(Deg1, 0, sizeof(Deg1));
	memset(P, 0, sizeof(P));
	scanf("%d%d", &N, &M);
	int S = N * 2 + 1, T = S + 1;
	G.Init(S, T);
	for (int i = 1; i <= M; i++) {
		int X, Y, L, Ls, Lt;
		scanf("%d%d%d%d%d", &X, &Y, &L, &Ls, &Lt);
		Rx[i] = X;	Ry[i] = Y;
		Rl[i] = L;	Rs[i] = Ls;	Rt[i] = Lt;

		G.Add(X, Y + N, Inf, Lt - L);
		G.Add(Y + N, X, Inf, L - Ls);
		
		Sum += L;
		Deg0[X]++;
		Deg1[Y]++;
		V[Y + N][X] = Lt - L;
		V[X][Y + N] = L - Ls;
	}
	
	for (int i = 1; i <= N; i++) {
		G.Add(S, i, Deg0[i], -Inf);
		Sum += Inf * Deg0[i];
		G.Add(i + N, T, Deg1[i], -Inf);
		Sum += Inf * Deg1[i];
		
		G.Add(i, T, Inf, 0);
		G.Add(i + N, T, Inf, 0);
	}
	
	Sum += G.MinCostFlow();
	if (Sum < -Inf * Inf / 2) puts("Unlike");
	else {
		GetAnswer();
		printf("%lld\n", Sum);
		for (int i = 1; i < N; i++) printf("%d ", P[i]); printf("%d\n", P[N]);
		for (int i = 1; i < N; i++) printf("%d ", P[i + N]); printf("%d", P[N + N]);
		printf("\n");
	}
}

int main()
{ int T;
  scanf("%d",&T);
  for(inti=1;i<=T;i++) Solve();
  return 0;
}

