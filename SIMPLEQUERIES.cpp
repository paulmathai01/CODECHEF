#include <set>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#define ll long long
#define rint register int
using namespace std;

const int Mod = 1000000007, MaxN = 200010, K = 8;
int cmod(int X) {
	return X >= Mod ? X - Mod : X;
}
struct Data_Struct {
	int C1, C2, C3;
	Data_Struct() {}
	Data_Struct(int D) {
		C1 = D;
		C2 = C3 = 0;
	}
	Data_Struct operator += (const Data_Struct &A) {
		C3 = cmod(C3 + A.C3);
		C2 = cmod(C2 + A.C2);
		C1 = cmod(C1 + A.C1);
		return *this;
	}
	Data_Struct operator - (const Data_Struct &A) const {
		Data_Struct Ret;
		Ret.C1 = cmod(C1 - A.C1 + Mod);
		Ret.C2 = cmod(C2 - A.C2 + Mod);
		Ret.C3 = cmod(C3 - A.C3 + Mod);
		return Ret;
	}
}	null;

int Qry[MaxN][3], Q, N, Val[MaxN], Pre[MaxN];

char c; int fret;
int fastin() {
	for (c = getchar(); (c < '0') || (c > '9'); c = getchar());
	for (fret = c - '0', c = getchar(); (c >= '0') && (c <= '9'); c = getchar()) fret = fret * 10 + c - '0';
	return fret;
}

namespace Initiation {
	int Ord[MaxN], Dfn;
	struct Tree_Node {
		int Son[2], F, Size, Cnt;
	}	T[MaxN];
	int Root, Total, Left;

	inline void Update(int Now) {
		T[Now].Size = T[Now].Cnt;
		if (T[Now].Son[0]) T[Now].Size += T[T[Now].Son[0]].Size;
		if (T[Now].Son[1]) T[Now].Size += T[T[Now].Son[1]].Size;
	}
	void Rotate(int U) {
		rint F = T[U].F, P = (T[F].Son[1] == U), Q = (T[T[F].F].Son[1] == F);
		T[F].Son[P] = T[U].Son[!P];	if (T[F].Son[P]) T[T[F].Son[P]].F = F;
		T[U].F = T[F].F;	T[T[F].F].Son[Q] = U;
		T[U].Son[!P] = F;	T[F].F = U;
		Update(F);	Update(U);
	}
	void Splay(int U, int Root) {
		for (; T[U].F != Root; Rotate(U)) {
			rint F = T[U].F, P = (T[F].Son[1] == U), Q = (T[T[F].F].Son[1] == F);
			if (T[F].F == Root) return Rotate(U);
			Rotate(P == Q ? F : U);
		}
	}

	int Find(int Now, int Ind) {
		if (T[T[Now].Son[0]].Size >= Ind) return Find(T[Now].Son[0], Ind);
		Ind -= T[T[Now].Son[0]].Size;
		if (Ind == T[Now].Cnt) return Now;
		return Find(T[Now].Son[1], Ind - T[Now].Cnt);
	}
	int FindLeft(int Now) {
		if (T[Now].Son[0]) return FindLeft(T[Now].Son[0]);
		else return Now;
	}
	void Dfs(int Now) {
		if (Now == 0) return;
		Dfs(T[Now].Son[0]);
		Ord[Now] = ++Dfn;
		Dfs(T[Now].Son[1]);
	}

	void Insert(int i, int Now) {
		T[Now].Cnt = 1;
		int F = Find(T[Root].Son[0], i);
		Splay(F, Root);
		T[Now].Son[1] = T[F].Son[1];
		T[T[Now].Son[1]].F = Now;
		T[F].Son[1] = Now;
		T[Now].F = F;
		Update(Now);	Update(F);
	}

	pair<int, int> Main() {
		Root = ++Total;
		Left = ++Total;
		T[Left].Size = T[Left].Cnt = 1;
		T[Left].F = Root;
		T[Root].Son[0] = Left;
		Update(Root);

		int N, Q, QTot = 0;
		N = fastin(); Q = fastin();
		for (int i = 1; i <= N; i++) {
			int X = fastin(), Now = ++Total;
			T[Now].Cnt = T[Now].Size = 1;
			Qry[++QTot][1] = Now;
			Qry[QTot][2] = X;
			Insert(i, Now);
		}
		for (int i = 1; i <= Q; i++) {
			++QTot;
			Qry[QTot][0] = fastin();
			if (Qry[QTot][0] == 1) {
				Qry[QTot][1] = fastin(), Qry[QTot][2] = fastin();
				Qry[QTot][1] = Find(T[Root].Son[0], Qry[QTot][1] + 1);
				Splay(Qry[QTot][1], Root);
				Qry[QTot][2] = Find(T[Root].Son[0], Qry[QTot][2] + 1);
				Splay(Qry[QTot][2], Root);
			} else
			if (Qry[QTot][0] == 2) {
				Qry[QTot][0] = 0;
				int X = fastin();
				Qry[QTot][1] = Find(T[Root].Son[0], X + 1);
				Splay(Qry[QTot][1], Root);
				Qry[QTot][2] = fastin();
			}else
			if (Qry[QTot][0] == 3) {
				Qry[QTot][0] = 0;
				int X = fastin();
				Qry[QTot][1] = Find(T[Root].Son[0], X + 1);
				Splay(Qry[QTot][1], Root);
				T[Qry[QTot][1]].Cnt = 0;
				T[Qry[QTot][1]].Size--;
			}else
			if (Qry[QTot][0] == 4) {
				Qry[QTot][0] = 0;
				int key, Now, ind;
				ind = fastin(); key = fastin();
				int F = Find(T[Root].Son[0], ind + 1);
				Splay(F, Root);
				Now = FindLeft(T[F].Son[0]);
				if ((Now == 0) || (T[Now].Cnt == 1)) {
					Now = ++Total;
					Insert(ind + 1, Now);
				} else {
					Splay(Now, Root);
					T[Now].Cnt = 1;
					Update(Now);
				}
				Qry[QTot][1] = Now;
				Qry[QTot][2] = key;
			}else
			if (Qry[QTot][0] == 5) {
				Qry[QTot][1] = fastin(), Qry[QTot][2] = fastin();
				Qry[QTot][1] = Find(T[Root].Son[0], Qry[QTot][1] + 1);
				Splay(Qry[QTot][1], Root);
				Qry[QTot][2] = Find(T[Root].Son[0], Qry[QTot][2] + 1);
				Splay(Qry[QTot][2], Root);
			}
		}

		bool Flag = 0;
		Dfs(T[Root].Son[0]);
		for (int i = 1; i <= QTot; i++) {
			Qry[i][1] = Ord[Qry[i][1]] - 1;
			if ((Qry[i][0] == 1) || (Qry[i][0] == 5)) {
				Qry[i][2] = Ord[Qry[i][2]] - 1;
				Flag = 1;
			}
		}
		if (Flag == 0) exit(0);
		return make_pair(Total - 2, QTot);
	}
}

set <ll> Set;
typedef set <ll> :: iterator LLI;

namespace Tree {
	const int MaxT = (1 << 23) + 23;
	Data_Struct Ans, Data[MaxT];
	int AnsNum;
	int Root[MaxN], Lenth[MaxN], Total, Son[MaxT][K], Size[MaxT];

	int P, D, D2, D3;
	bool Type;
	inline void Excute(Data_Struct &Now) {
		if (Type) {
			Now.C1 = cmod(Now.C1 + D);
			Now.C2 = cmod(Now.C2 + D2);
			Now.C3 = cmod(Now.C3 + D3);
		}
		else {
			Now.C1 = cmod(Now.C1 - D + Mod);
			Now.C2 = cmod(Now.C2 - D2 + Mod);
			Now.C3 = cmod(Now.C3 - D3 + Mod);
		}
	}
	void TModify(int &Now, int L, int R) {
		if (Now == 0) Now = ++Total;
		if (!Type) Size[Now]--;
		else Size[Now]++;
		if (Size[Now] == 0) { Now = 0; return; }
		if (L == R) return Excute(Data[Now]);
		rint Step = Lenth[R - L + 1];
		for (int i = 0, r = L + Step; i < K; i++, r += Step)
			if (P < r) {
				TModify(Son[Now][i], r - Step, r - 1);
				if (Size[Now] == Size[Son[Now][i]]) Data[Now] = Data[Son[Now][i]];
				else Excute(Data[Now]);
				return;
			}
		return;
	}
	void TQuery(int Now, int L, int R) {
		if (!Size[Now] || (L > P)) return;
		rint Step = Lenth[R - L + 1];
		for (int i = 0, r = L + Step; i < K; i++, r += Step)
			if (r - 1 <= P)  {
				if (Size[Son[Now][i]]) Ans += Data[Son[Now][i]];
			} else {
				TQuery(Son[Now][i], r - Step, r - 1);
				return;
			}
	}
	void TQueryNum(int Now, int L, int R) {
		if (!Size[Now] || (L > P)) return;
		rint Step = Lenth[R - L + 1];
		for (int i = 0, r = L + Step; i < K; i++, r += Step)
			if (r - 1 <= P) AnsNum += Size[Son[Now][i]];
			else {
				TQueryNum(Son[Now][i], r - Step, r - 1);
				return;
			}
	}

	inline void Modify(int ind, int _Type) {
		P = Pre[ind], D = Val[ind], D2 = (ll)D * D % Mod, D3 = (ll)D2 * D % Mod;
		Type = (_Type > 0);
		for (int i = ind; i <= N; i += (i & -i))
			TModify(Root[i], 0, N - 1);
	}
	inline Data_Struct Query(int ind) {
		Ans = null;
		for (int i = ind; i; i -= (i & -i))
			TQuery(Root[i], 0, N - 1);
		return Ans;
	}
	inline int QueryNum(int ind) {
		AnsNum = 0;
		for (int i = ind; i; i -= (i & -i))
			TQueryNum(Root[i], 0, N - 1);
		return AnsNum;
	}
	inline Data_Struct Query(int L, int R, int _Lim) {
		P = _Lim;
		Data_Struct Ret = Query(R) - Query(L - 1);
		return Ret;
	}
	inline int QueryNum(int L, int R, int _Lim) {
		P = _Lim;
		return QueryNum(R) - QueryNum(L - 1);
	}
}

int main()
{
	ll inv6 = 166666668;
	pair<int, int> Init = Initiation::Main();
	N = Init.first, Q = Init.second;

	for (int i = 1; i <= N; i++) Tree::Lenth[i] = i / K + 1;

	for (int i = 1; i <= Q; i++) {
		if (Qry[i][0] == 0) {
			int U = Qry[i][1];
			if (Val[U] != 0) {
				ll Now = (ll)Val[U] * Mod + U;
				Set.erase(Now);
				LLI Nxt = Set.upper_bound(Now);
				if ((Nxt != Set.end()) && ((*Nxt) / Mod == Val[U])) {
					int NPos = (*Nxt) % Mod, NPre;
					if (Nxt == Set.begin()) NPre = 0;
					else {
						Nxt--;
						if ((*Nxt) / Mod == Val[U]) NPre = (*Nxt) % Mod;
						else NPre = 0;
					}
					Tree::Modify(NPos, -1);
					Pre[NPos] = NPre;
					Tree::Modify(NPos, 1);
				}
				Tree::Modify(U, -1);
			}

			Val[U] = Qry[i][2];
			if (Val[U] != 0) {
				ll Now = (ll)Val[U] * Mod + U;
				Set.insert(Now);
				LLI Nxt = Set.upper_bound(Now);
				if ((Nxt != Set.end()) && ((*Nxt) / Mod == Val[U])) {
					int NPos = (*Nxt) % Mod;
					Tree::Modify(NPos, -1);
					Pre[NPos] = U;
					Tree::Modify(NPos, 1);
				}
				int NPre = 0;
				Nxt--;
				if (Nxt != Set.begin())  {
					Nxt--;
					if ((*Nxt) / Mod == Val[U]) NPre = (*Nxt) % Mod;
				}
				Pre[U] = NPre;
				Tree::Modify(U, 1);
			}
		} else {
			if (Qry[i][0] == 5) printf("%d\n", Tree::QueryNum(Qry[i][1], Qry[i][2], Qry[i][1] - 1));
			else {
				Data_Struct Ans = Tree::Query(Qry[i][1], Qry[i][2], Qry[i][1] - 1);
				ll O = (ll)Ans.C1 * Ans.C1 % Mod * Ans.C1 % Mod - 3ll * Ans.C1 * Ans.C2 % Mod + 2ll * Ans.C3;
				O = O % Mod * inv6 % Mod;
				if (O < 0) O += Mod;
				printf("%lld\n", O);
			}
		}
	}
	return 0;
}
