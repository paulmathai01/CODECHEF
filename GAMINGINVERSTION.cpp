#include <iostream>
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <utility>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <cstring>
#include <cassert>
#define rf freopen("in.in", "r", stdin)
#define wf freopen("out.out", "w", stdout)
#define rep(i, s, n) for(int i=int(s); i<=int(n); ++i)
using namespace std;
const int mx = 1e5+10, mod = 1e9+7;
 
int a[mx], fa[mx], mrk[mx], ord[mx], qlen[mx];
int n, q, fidx, oidx=-1;
long long CNT;
 
//------------ Treaps + BIT -------//
void up(int i, int v)
{
	for(; i<mx; i+= i&-i)
		mrk[i] += v;
}
int qu(int i)
{
	int ret = 0;
	for(; i; i-= i&-i)
		ret += mrk[i];
	return ret;
}
int search(int val)
{
	int s = 1, e = mx-1;
	while(s<=e)
	{
		int m = (s+e)>>1;
		if(qu(m)>=val) e=m-1;
		else s=m+1;
	}
	return s;
}
 
struct node
{
	int pr, sz, val, i;
	node* l, *r;
 
	node(node* put, int ins, int priority, int s)
	{
		l=r=put;
		val=ins; sz=s; i=++oidx;
		pr = priority;
	}
 
}; node* dummy, *root;
 
node* rotate_right(node* cur)
{
	node *x = cur->l, *y = x->r;
	x->r=cur, cur->l = y;
	
	cur->sz = cur->l->sz + cur->r->sz + 1;
	x->sz = x->l->sz + x->r->sz + 1;
	
	return x;
}
node* rotate_left(node* cur)
{
	node *x = cur->r, *y = x->l;
	x->l=cur, cur->r = y;
	
	cur->sz = cur->l->sz + cur->r->sz + 1;
	x->sz = x->l->sz + x->r->sz + 1;
	
	return x;
}
 
void insert(node* &cur, int pos, int vins)
{
	if(cur == dummy)
	{
		int priority = rand()%int(1e9) + 1;
		return void(cur = new node(dummy, vins, priority, 1));
	}
 
	if(cur->l->sz+1 >= pos)
	{
		insert(cur->l, pos, vins); cur->sz+=1;
		if(cur->l->pr < cur->pr) cur = rotate_right(cur);
	}
	else
	{
		pos -= cur->l->sz+1, insert(cur->r, pos, vins); cur->sz+=1;
		if(cur->r->pr < cur->pr) cur = rotate_left(cur);
	}
}
 
void inorder(node* &cur)
{
	if(cur == dummy) return;
 
	inorder(cur->l);
	fa[++fidx] = cur->val, ord[cur->i] = fidx;
	inorder(cur->r);
}
 
//---------------------------------//
 
//-----------Segment tree+BIT------//
vector<int> st[mx*4];
int *sb[mx*4];
 
void build(int x, int y, int r)
{
	if(x==y)
	{
		st[r].push_back(a[fa[x]]);
		sb[r] = new int[2];
		return void(sb[r][0]=sb[r][1]=0);
	}
	
	int m=(x+y)>>1;
	build(x, m, 2*r); build(m+1, y, 2*r+1);
 
	merge(st[2*r].begin(), st[2*r].end(), st[2*r+1].begin(), st[2*r+1].end(), back_inserter(st[r]));
	sb[r] = new int[y-x+2];
	rep(i, 0, y-x+1) sb[r][i]=0;
}
 
void up(int i, int* bitv, int sz)
{
	for(; i<sz; i+= i&-i)
		bitv[i] += 1;
}
int qu(int i, int* bitv)
{
	int ret = 0;
	for(; i; i-= i&-i)
		ret += bitv[i];
	return ret;
}
 
void query(int x, int y, int r, int k, int val)
{
	if(x == y)
	{
		up(1, sb[r], y-x+2);
		return;
	}
 
	int m = (x+y)>>1;
	if(k>m)
	{
		int idx = upper_bound(st[2*r].begin(), st[2*r].end(), val)-st[2*r].begin();
		CNT += qu(st[2*r].size(), sb[2*r]) - qu(idx, sb[2*r]);
 
		query(m+1, y, 2*r+1, k, val);
	}
	else
	{
		int idx = lower_bound(st[2*r+1].begin(), st[2*r+1].end(), val)-st[2*r+1].begin();
		CNT += qu(idx, sb[2*r+1]);
 
		query(x, m, 2*r, k, val);
	}
 
	int upidx = lower_bound(st[r].begin(), st[r].end(), val)-st[r].begin();
	up(upidx+1, sb[r], y-x+2);
}
//---------------------------------//
 
int main()
{
	//rf; wf;
 
	dummy = new node(NULL, 0, -mx, 0);
	root = dummy;
 
	scanf("%d %d", &n, &q);
	rep(i, 1, n)
	{
		scanf("%d", &a[i]);
		up(i, 1);
	}
 
	rep(i, 1, q)
	{
		int x, y, k;
		scanf("%d %d %d", &x, &y, &k);
		qlen[i] = y-x+1;
 
		rep(j, x, y)
		{
			int idx = search(x);
			up(idx, -1);
			insert(root, k, idx); ++k;
		}
	}
	inorder(root);
	build(1, n, 1);
 
	rep(i, 1, q)
	{
		rep(j, 1, qlen[i])
		{
			int k = ord[qlen[i-1]+j], val = a[fa[k]];
			query(1, n, 1, k, val);
		}
		
		qlen[i]+=qlen[i-1];
		printf("%lld\n", CNT);
	}
 
	return 0;
} 