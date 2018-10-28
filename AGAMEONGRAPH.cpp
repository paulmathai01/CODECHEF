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
 

