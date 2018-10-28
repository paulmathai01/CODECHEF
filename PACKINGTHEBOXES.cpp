#include <cstdio>
#include <algorithm>
#include <vector>
#include <set>
using namespace std;
const long long p = 1000000000;
 
inline long long mul(long long x,long long y)
{
    return (x*y)%p;
}
 
int main()
{
    int k;
    scanf("%d",&k);
    set<pair<long long,int> > s;
    for(int i=0;i<k;i++)
    {
        long long x;
        scanf("%lld",&x);
        s.insert(make_pair(x,i));
    }
    vector<pair<int,int> > e;
    for(int i=k;i<2*k-1;i++)
    {
        pair<long long,int> f1 = *s.begin();
        s.erase(s.begin());
        pair<long long,int> f2 = *s.begin();
        s.erase(s.begin());
        s.insert(make_pair(f1.first+f2.first, i));
        e.push_back(make_pair(i,f1.second));
        e.push_back(make_pair(i,f2.second));
    }
    int root = 2*k-2;
    vector<int> deg(root+1);
    sort(e.begin(),e.end());
    for(int i=e.size()-1;i>=0;i--)
        deg[e[i].second] = deg[e[i].first]+1;
    vector<int> cnt(k);
    for(int i=0;i+1<deg.size();i++) //root-a nie liczymy
        cnt[deg[i]]++;
   vector<long long> tmp(k/2+1);
   tmp[0]=1;
   for(int i=1;i<=k/2;i++)
   {
       tmp[i] = mul(tmp[i-1],(2*i-1));
   }
   long long ans = 1;
   for(int i=0;i<k;i++)
       ans = mul(ans,tmp[cnt[i]/2]);
   printf("%lld %lld\n",s.begin()->first,ans);
}
 