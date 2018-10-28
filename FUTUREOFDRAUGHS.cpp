#ifdef _MSC_VER
#define _CRT_SECURE_NO_DEPRECATE
#pragma comment(linker, "/STACK:66777216")
#else
#pragma GCC optimize("O3")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx")
#endif

#include <bits/stdc++.h>

using namespace std;

const int MAX_N = 2000;

int n;

int size[MAX_N];
short g[MAX_N][MAX_N];
int bbase[MAX_N];
int base[MAX_N];
int match[MAX_N];
int tu;
int used[MAX_N];
int tb;
int blossom[MAX_N];

#define gc getchar_unlocked

int read_int() {
char c = gc();
while (!isdigit(c)) c = gc();
int ret = 0;
while (isdigit(c)) {
ret = (ret<<1) + (ret<<3) + (c&15);
c = gc();
}
return ret;
}

int tp;
int lp[MAX_N];
int p[MAX_N];
int tim;
int used_lca[MAX_N];

int lca(int a, int b) {
++tim;
do {
used_lca[a] = tim;
if (match[a] == -1) break;
a = base[p[match[a]]];
}
while (true);
while (used_lca[b] < tim) {
b = base[p[match[b]]];
}
return b;
}

int qEnd;
int q[MAX_N];

void mark_path(int v, int b, int children) {
while (base[v] != b) {
if (used[v] < tu) {
used[v] = tu;
q[qEnd++] = v;
}
if (used[match[v]] < tu) {
used[match[v]] = tu;
q[qEnd++] = match[v];
}
blossom[base[v]] = blossom[base[match[v]]] = tb;
p[v] = children;
children = match[v];
v = p[children];
}
}

bool print;
int vt;
int avt[MAX_N];

int find_path(int root) {
++tp;
	memcpy(base, bbase, n<<2);
used[root] = ++tu;
int qBeg = 0;
qEnd = 1;
q[0] = root;
while (qBeg < qEnd) {
int v = q[qBeg++];	
if (print) {
if (avt[v] < vt) 
avt[v] = vt;
else
continue;
}
for (int i = 0; i < size[v]; ++i) {
int to = g[v][i];
if (base[v] == base[to] || match[v] == to) continue;
if (to == root || match[to] != -1 && lp[match[to]] == tp) {
int cur_base = lca(base[v], base[to]);
++tb;
mark_path(v, cur_base, to);
mark_path(to, cur_base, v);
for (int i = 0; i < n; ++i) {
if (blossom[base[i]] == tb) {
base[i] = cur_base;
}
}
}
else if (lp[to] < tp) {
lp[to] = tp;
p[to] = v;
if (match[to] == -1) return to;
to = match[to];
used[to] = tu;
q[qEnd++] = to;
}
}
}
return -1;
}

int main() {
#ifndef ONLINE_JUDGE
	freopen("sample.in", "r", stdin);
#endif
tim = 0;
	memset(used_lca, 0, MAX_N<<2);
tu = 0;
	memset(used, 0, MAX_N<<2);
tb = 0;
	memset(blossom, 0, MAX_N<<2);
tp = 0;
	memset(lp, 0, MAX_N<<2);
vt = 0;
	memset(avt, 0, MAX_N<<2);
for (int i = 0; i < MAX_N; ++i) {
bbase[i] = i;
}
int nTests = read_int();
while (nTests--) {
n = read_int(); 
int m = read_int();
		memset(size, 0, n<<2);
for (int edge = 0; edge < m; ++edge) {
int x = read_int()-1;
int y = read_int()-1;
g[x][size[x]++] = y;
g[y][size[y]++] = x;
}
		memset(match, -1, n<<2);

print = false;
for (int i = n-1; i >= 0; --i) {
if (match[i] == -1) {
int v = find_path(i);
while (v != -1) {
int pv = p[v];
int ppv = match[pv];
match[v] = pv;
match[pv] = v;
v = ppv;
}
}
}
print = true;
++vt;
for (int i = 0; i < n; ++i) {
if (match[i] == -1) {
find_path(i);
}
}
int answer = 0;
for (int i = 0; i < n; ++i) {
answer += avt[i] == vt;
}
		printf("%d\n", answer);
}
return 0;
} 

#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <set>
#include <map>
#include <queue>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <cstring>
#include <cctype>
#include <cassert>
#include <limits>
#include <functional>
#define rep(i,n) for(int (i)=0;(i)<(int)(n);++(i))
#define rer(i,l,u) for(int (i)=(int)(l);(i)<=(int)(u);++(i))
#define reu(i,l,u) for(int (i)=(int)(l);(i)<(int)(u);++(i))
#if defined(_MSC_VER) || __cplusplus > 199711L
#define aut(r,v) auto r = (v)
#else
#define aut(r,v) __typeof(v) r = (v)
#endif
#define each(it,o) for(aut(it, (o).begin()); it != (o).end(); ++ it)
#define all(o) (o).begin(), (o).end()
#define pb(x) push_back(x)
#define mp(x,y) make_pair((x),(y))
#define mset(m,v) memset(m,v,sizeof(m))
#define INF 0x3f3f3f3f
#define INFL 0x3f3f3f3f3f3f3f3fLL
using namespace std;
typedef vector<int> vi; typedef pair<int,int> pii; typedef vector<pair<int,int> > vpii; typedef long long ll;
template<typename T, typename U> inline void amin(T &x, U y) { if(y < x) x = y; }
template<typename T, typename U> inline void amax(T &x, U y) { if(x < y) x = y; }
 
template<int MOD>
struct ModInt {
	static const int Mod = MOD;
	unsigned x;
	ModInt(): x(0) { }
	ModInt(signed sig) { int sigt = sig % MOD; if(sigt < 0) sigt += MOD; x = sigt; }
	ModInt(signed long long sig) { int sigt = sig % MOD; if(sigt < 0) sigt += MOD; x = sigt; }
	int get() const { return (int)x; }
	
	ModInt &operator+=(ModInt that) { if((x += that.x) >= MOD) x -= MOD; return *this; }
	ModInt &operator-=(ModInt that) { if((x += MOD - that.x) >= MOD) x -= MOD; return *this; }
	ModInt &operator*=(ModInt that) { x = (unsigned long long)x * that.x % MOD; return *this; }
	ModInt &operator/=(ModInt that) { return *this *= that.inverse(); }
	
	ModInt operator+(ModInt that) const { return ModInt(*this) += that; }
	ModInt operator-(ModInt that) const { return ModInt(*this) -= that; }
	ModInt operator*(ModInt that) const { return ModInt(*this) *= that; }
	ModInt operator/(ModInt that) const { return ModInt(*this) /= that; }
	
	ModInt inverse() const {
		signed a = x, b = MOD, u = 1, v = 0;
		while(b) {
			signed t = a / b;
			a -= t * b; std::swap(a, b);
			u -= t * v; std::swap(u, v);
		}
		if(u < 0) u += Mod;
		ModInt res; res.x = (unsigned)u;
		return res;
	}
 
	ModInt operator-() const { ModInt t; t.x = x == 0 ? 0 : Mod - x; return t; }
};
typedef ModInt<1000000007> mint;
 
void truncateLeadingZerosOfPolynomial(const vector<mint> &p, vector<mint> &res) {
	int n = (int)p.size();
	while(n > 0 && p[n-1].x == 0) -- n;
	res.assign(p.begin(), p.begin() + n);
}
 
int extendedEuclideanAlgorithm(const vector<mint> &f, const vector<mint> &g, int k, vector<mint> &out_r, vector<mint> &out_t) {
	if((!f.empty() && f.back().x == 0) || (!g.empty() && g.back().x == 0)) {
		vector<mint> nf, ng;
		truncateLeadingZerosOfPolynomial(f, nf);
		truncateLeadingZerosOfPolynomial(g, ng);
		return extendedEuclideanAlgorithm(nf, ng, k, out_r, out_t);
	}
	if(f.size() < g.size())
		return extendedEuclideanAlgorithm(g, f, k, out_r, out_t);
	assert(k >= 0);
 
	int n[2];
	vector<mint> r[2], t[2];//, s[2];
	int nf = (int)f.size() - 1, ng = (int)g.size() - 1;
	n[0] = nf, n[1] = ng;
	int bufsize = n[0] == -1 ? 1 : n[0] + 1;
	rep(i, 2) {
		r[i].resize(bufsize);
		t[i].resize(bufsize);
	}
	rer(j, 0, n[0]) r[0][j] = f[j];
	rer(j, 0, n[1]) r[1][j] = g[j];
	t[1][0] = mint(1);
 
	if(n[0] < k) {
		out_r.swap(r[0]); out_r.resize(n[0] + 1);
		out_t.swap(t[0]); out_t.resize(1);
		return 0;
	}else if(n[1] < k) {
		out_r.swap(r[1]); out_r.resize(n[0] + 1);
		out_t.swap(t[1]); out_t.resize(1);
		return 1;
	}
 
	for(int i = 1; n[1] >= 0; ++ i) {
		swap(n[0], n[1]);
		r[0].swap(r[1]);
		t[0].swap(t[1]);
		//r[0] = r_i; r[1] = r_{i-1}, r_{i+1}
 
		mint ilc = r[0][n[0]].inverse();
		//q_i     <- r_{i-1} quo r_i
		//r_{i+1} <- r_{i-1} - q_i r_i
		//s_{i+1} <- s_{i-1} - q_i s_i
		//t_{i+1} <- t_{i-1} - q_i t_i
		int nr = n[0], ns = ng - n[1], nt = nf - n[1];
 
		for(int j = n[1] - nr; j >= 0; -- j) {
			mint c = r[1][nr + j] * ilc;
			for(int l = 0; l <= nr; ++ l)
				r[1][j + l] -= c * r[0][l];
			for(int l = 0; l <= nt; ++ l)
				t[1][j + l] -= c * t[0][l];
		}
		while(n[1] >= 0 && r[1][n[1]].x == 0)
			-- n[1];
		mint ip = n[1] == -1 ? mint(1) : r[1][n[1]].inverse();
 
		if(n[1] < k) {
			out_r.swap(r[1]); out_r.resize(n[1] + 1);
			out_t.swap(t[1]); out_t.resize(nf - n[0] + 1);
			return i+1;
		}
	}
	assert(false);
	return -1;
}
 
//�?¤�?½ï¿½�?¤�?¸ï¿½�?£ï¿½�?®�?¦ï¿½�?°�?¥ï¿½ï¿½�?£ï¿½�?®minimum polynomial�?£ï¿½ï¿½�?¨�?¨ï¿½�?§�?®ï¿½�?£ï¿½ï¿½�?£ï¿½ï¿½�?£ï¿½ï¿½
//�?¥ï¿½�?¥�?¥ï¿½ï¿½�?£ï¿½ï¿½linearly recurrent�?£ï¿½�?§�?£ï¿½�?ª�?£ï¿½ï¿½�?¥ �?´�?¥ï¿½ï¿½�?£ï¿½�?¯�?£ï¿½ï¿½res.back() == 0 �?£ï¿½�?�"�?£ï¿½�?ª�?£ï¿½�?£�?£ï¿½�?¦�?£ï¿½ï¿½
//res �?£ï¿½�?®�?¦ï¿½�?�"�?¥�?°�?¾�?£ï¿½�?®0�?£ï¿½�?®�?¦ï¿½�?°�?£ï¿½�?¯�?£ï¿½ï¿½�?£ï¿½ï¿½�?£ï¿½�?®�?¦ï¿½�?°�?£ï¿½ �?£ï¿½ï¿½ a �?£ï¿½�?®�?¦ï¿½�?�"�?¥�?°�?¾�?£ï¿½ï¿½�?¥ï¿½ï¿½�?£ï¿½�?£�?£ï¿½ï¿½�?£ï¿½ï¿½�?£ï¿½�?®�?£ï¿½ï¿½ a �?£ï¿½�?® linearly recurrent �?£ï¿½�?ª�?¦ï¿½ï¿½�?©ï¿½�?·�?£ï¿½�?® prefix �?£ï¿½�?§�?£ï¿½ï¿½�?£ï¿½ï¿½�?£ï¿½ï¿½�?£ï¿½�?¨�?£ï¿½ï¿½�?¨�?¡�?¨�?£ï¿½ï¿½�?£ï¿½ï¿½
void computeMinimumPolynomialForLinearlyRecurrentSequence(const vector<mint> &a, vector<mint> &res) {
	int n2 = (int)a.size(), n = n2 / 2;
	assert(n2 % 2 == 0);
	vector<mint> x(n2 + 1), s, t;
	x[n2] = 1;
	int i = extendedEuclideanAlgorithm(x, a, n, s, t);
	if(s.size() == n2 + 1) {	//a is the zero sequence
		res.assign(1, mint());
		res[0] = mint(1);
		return;
	}
 
	int d = max((int)s.size(), (int)t.size() - 1);
	int e = 0;
	while(t[e].x == 0) ++ e;
	mint c = t[e].inverse();
	res.assign(d + 1, mint());
	reu(i, e, t.size())
		res[d - i] = t[i] * c;
}
 
void checkLinearRecurrence(const vector<mint> &a, int n) {
	if((int)a.size() <= n * 2)
		return;
	vector<mint> v(a.begin(), a.begin() + n * 2), m;
	computeMinimumPolynomialForLinearlyRecurrentSequence(v, m);
	int d = (int)m.size() - 1;
	reu(i, d, a.size()) {
		mint x;
		rer(j, 0, d) x += m[j] * a[i - d + j];
		if(x.x != 0) {
			cerr << "it isn't linearly recurrent of order <= " << n << endl;
			return;
		}
	}
	cerr << "it is linearly recurrent of order " << d << endl;
}
 
vector<mint> fact, factinv;
void nCr_computeFactinv(int N) {
	N = min(N, mint::Mod - 1);
	fact.resize(N+1); factinv.resize(N+1);
	fact[0] = 1;
	rer(i, 1, N) fact[i] = fact[i-1] * i;
	factinv[N] = fact[N].inverse();
	for(int i = N; i >= 1; i --) factinv[i-1] = factinv[i] * i;
}
mint nCr(int n, int r) {
	if(n >= mint::Mod)
		return nCr(n % mint::Mod, r % mint::Mod) * nCr(n / mint::Mod, r / mint::Mod);
	return r > n ? 0 : fact[n] * factinv[n-r] * factinv[r];
}
 
#if defined(WIN32) && !defined(_WINDOWS_)
struct FILETIME {
	unsigned dwLowDateTime, dwHighDateTime;
};
extern "C" void* __stdcall GetCurrentProcess(void);
extern "C" int __stdcall GetProcessTimes(void *hProcess, FILETIME *lpCreationTime, FILETIME *lpExitTime, FILETIME *lpKernelTime, FILETIME *lpUserTime);
#endif
#ifndef WIN32
#include <sys/time.h>
#include <sys/resource.h>
#endif
 
void getCPUTime(double &userTime, double &systemTime) {
#ifdef WIN32
	void *handle = GetCurrentProcess();
	FILETIME dummy1, dummy2, kernel, user;
	GetProcessTimes(handle, &dummy1, &dummy2, &kernel, &user);
	userTime = user.dwHighDateTime * 429.4967296 + user.dwLowDateTime * 1e-7;
	systemTime = kernel.dwHighDateTime * 429.4967296 + kernel.dwLowDateTime * 1e-7;
#else
	struct rusage ru;
	getrusage(RUSAGE_SELF, &ru);
	userTime = ru.ru_utime.tv_sec + ru.ru_utime.tv_usec * 1e-6;
	systemTime = ru.ru_stime.tv_sec + ru.ru_stime.tv_usec * 1e-6;
#endif
}
double getCPUTime() {
	double user, sys;
	getCPUTime(user, sys);
	return user + sys;
}
 
struct CPUTimeIt {
	double user, sys;
	const char *msg;
 
	CPUTimeIt(const char *msg_): msg(msg_) { getCPUTime(user, sys); }
	~CPUTimeIt() {
		double userEnd, sysEnd;
		getCPUTime(userEnd, sysEnd);
		fprintf(stderr, "%s: user %.6fs / sys %.6fs\n", msg, userEnd - user, sysEnd - sys);
	}
 
	operator bool() { return false; }
};
#define CPUTIMEIT(s) if(CPUTimeIt cputimeit_##__LINE__ = s); else
 
 
int main() {
	nCr_computeFactinv(10000);
	const int NN = 50, KK = 10000;
 
	int T;
	while(~scanf("%d", &T)) {
		vector<vector<vi> > graphs(T);
		rep(i, T) {
			int N, M;
			scanf("%d%d", &N, &M);
//			N=NN,M=N*(N-1)/2;
			vector<vi> g(N);
			rep(i, M) {
				int x, y;
				scanf("%d%d", &x, &y), -- x, -- y;
//				x=rand()%N,y=rand()%N;
				g[x].push_back(y);
				g[y].push_back(x);
			}
			graphs[i] = g;
		}
		int Q;
		scanf("%d", &Q);
		vector<vpii> queries((T+1) * (T+1));
		int MaxK = 1;
		rep(i, Q) {
			int L, R, K;
			scanf("%d%d%d", &L, &R, &K), -- L;
//			L=rand()%T,R=rand()%T;if(L>R)swap(L,R);K=rand()%KK+1;
			queries[L * (T+1) + R].push_back(mp(K, i));
			amax(MaxK, K);
		}
		vector<vector<mint> > nums(T, vector<mint>(MaxK+1));
//		CPUTIMEIT("calculate nums")
		rep(i, T) {
			const vector<vi> &g = graphs[i];
			int N = g.size();
			int K = min(MaxK, N * 2 - 1);
			vector<mint> dp, ndp;
			rep(s, N) {
				nums[i][0] += 1;
				ndp.assign(N, mint());
				ndp[s] = 1;
				rep(k, K) {
					ndp.swap(dp);
					ndp.assign(N, mint());
					rep(i, N) {
						mint x = dp[i];
						if(x.x == 0) continue;
						ndp[i] += x;
						each(j, g[i])
							ndp[*j] += x;
					}
					nums[i][k+1] += ndp[s];
				}
			}
			if(K < MaxK) {
				vector<mint> a(nums[i].begin(), nums[i].begin() + (K + 1)), m;
				computeMinimumPolynomialForLinearlyRecurrentSequence(a, m);
				assert(m.back().x == 1);
				int d = (int)m.size() - 1;
				rer(k, K + 1, MaxK) {
					mint x;
					rep(j, d)
						x += m[j] * nums[i][k - d + j];
					nums[i][k] = -x;
				}
			}
		}
 
		const int XX = 8;
 
		vector<vector<mint> > coefs(MaxK+1);
		rer(i, 0, MaxK) {
			coefs[i].resize((i+1) + XX);
			coefs[i][0] = (i + 1) % 2;
			reu(j, 1, i)
				coefs[i][j] = coefs[i-1][j-1] - coefs[i-1][j];
			coefs[i][i] = 1;
		}
		rer(i, 0, MaxK)
			coefs[i][0] -= 1;
 
		vector<mint> ans(Q);
//		CPUTIMEIT("calculate ans")
		rep(L, T) {
			vector<mint> prods((MaxK+1) + XX, mint(1));
			rer(R, L+1, T) {
				rer(k, 0, MaxK)
					prods[k] *= nums[R-1][k];
 
				vpii &v = queries[L * (T+1) + R];
				sort(all(v));
 
				rep(i, v.size()) {
					int K = v[i].first;
					if(i > 0 && v[i-1].first == K) {
						ans[v[i].second] = ans[v[i-1].second];
						continue;
					}
					const unsigned *pcoefs = reinterpret_cast<const unsigned*>(&coefs[K][0]);
					const unsigned *pprods = reinterpret_cast<const unsigned*>(&prods[0]);
//					mint sum;
//					rer(k, 0, K)
//						sum += coefs[K][k] * prods[k];
					unsigned sum = 0;
					int k;
					for(k = 0; k <= K; k += XX) {
						unsigned long long t = 0ULL;
						t += (unsigned long long)pcoefs[k] * pprods[k];
						t += (unsigned long long)pcoefs[k+1] * pprods[k+1];
						t += (unsigned long long)pcoefs[k+2] * pprods[k+2];
						t += (unsigned long long)pcoefs[k+3] * pprods[k+3];
						t += (unsigned long long)pcoefs[k+4] * pprods[k+4];
						t += (unsigned long long)pcoefs[k+5] * pprods[k+5];
						t += (unsigned long long)pcoefs[k+6] * pprods[k+6];
						t += (unsigned long long)pcoefs[k+7] * pprods[k+7];
						unsigned u = t % mint::Mod;
						if((sum += u) >= mint::Mod) sum -= mint::Mod;
					}
					ans[v[i].second].x = sum;
				}
			}
		}
		rep(i, Q)
			printf("%d\n", ans[i].get());
	}
	return 0;
}
 


