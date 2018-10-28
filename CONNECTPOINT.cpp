#include<algorithm>
#include<cctype>
#include<climits>
#include<cmath>
#include<complex>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<ctime>
#include<deque>
#include<fstream>
#include<functional>
#include<iomanip>
#include<iostream>
#include<list>
#include<map>
#include<numeric>
#include<queue>
#include<set>
#include<sstream>
#include<stack>
#include<string>
#include<utility>
#include<vector>
#define sf scanf
#define si(a)scanf("%d",&a)
#define pf printf
#define pi(a)printf("%d",a)
#define fr freopen
#define ps for(;;)
#define pb push_back
#define mp make_pair
#define lp(i,a,b)for(int i=a;i<=int(b);++i)
#define rp(i,a,b)for(int i=int(a);i>=b;--i)
#define vp(i,a)for(int i=0;i<int(a.size());++i)
#define wp(i,a)for(int i=int(a.size())-1;i>=0;--i)
#define mc int T;scanf("%d",&T);for(int I=1;I<=T;++I)
using namespace std;
typedef long long ll;
typedef long double ld;
typedef vector<int>vi;
typedef pair<int,int>pii;
class MaximalPlanarGraph{
private:
    int n,m;
    std::vector<std::set<int> >to2;
    std::vector<std::vector<int> >to;
    std::vector<int>dec,removed,mark,invc,rt;
    std::vector<std::list<int>::iterator>dpos,pos;
    bool order(const int&v1,const int&v2,const int&vn){
        rt[0]=v1;
        rt[1]=v2;
        rt[n-1]=vn;
        std::fill(invc.begin(),invc.end(),0);
        invc[v1]=1;
        invc[v2]=1;
        invc[vn]=1;
        std::list<int>deg;
        dpos[vn]=deg.insert(deg.begin(),vn);
        std::fill(dec.begin(),dec.end(),0);
        dec[v1]=2;
        dec[v2]=2;
        dec[vn]=2;
        for(int i=n-1;i>=2;--i){
            if(deg.empty())
                return false;
            int v=*deg.begin();
            deg.erase(deg.begin());
            invc[v]=-1;
            rt[i]=v;
            for(int u:to[v]){
                if(invc[u]==1){
                    if(u!=v1&&u!=v2&&dec[u]==2)
                        deg.erase(dpos[u]);
                    --dec[u];
                    if(u!=v1&&u!=v2&&dec[u]==2)
                        dpos[u]=deg.insert(deg.begin(),u);
                }else if(invc[u]==0){
                    invc[u]=2;
                }
            }
            for(int u:to[v]){
                if(invc[u]==2)
                    for(int w:to[u]){
                        if(invc[w]==1){
                            if(w!=v1&&w!=v2&&dec[w]==2)
                                deg.erase(dpos[w]);
                            ++dec[w];
                            if(w!=v1&&w!=v2&&dec[w]==2)
                                 dpos[w]=deg.insert(deg.begin(),w);
                            ++dec[u];
                        }else if(invc[w]==2)
                            ++dec[u];
                    }
            }
            for(int u:to[v]){
                if(invc[u]==2){
                    invc[u]=1;
                    if(dec[u]==2)
                        dpos[u]=deg.insert(deg.begin(),u);
                }
            }
        }
        return true;
    }
    bool embed(){
        std::list<int>ext;
        std::fill(mark.begin(),mark.end(),0);
        int marker=0;
        pos[rt[1]]=ext.insert(ext.begin(),rt[1]);
        pos[rt[2]]=ext.insert(ext.begin(),rt[2]);
        pos[rt[0]]=ext.insert(ext.begin(),rt[0]);
        std::fill(removed.begin(),removed.end(),0);
        removed[rt[1]]=1;
        removed[rt[2]]=1;
        removed[rt[0]]=1;
        for(int i=3;i<n;++i){
            int v=rt[i];
            removed[v]=1;
            std::vector<int>can;
            ++marker;
            for(int u:to[v])
                if(removed[u]){
                    mark[u]=marker;
                    can.push_back(u);
                }
            int start=-1,end=-1;
            for(int u:can){
                std::list<int>::iterator it=pos[u];
                if(it==std::list<int>::iterator())
                    return false;
                if(it==ext.begin()){
                    if(start!=-1)
                        return false;
                    start=u;
                }else{
                    std::list<int>::iterator tmp=it;
                    if(mark[*(--tmp)]!=marker){
                        if(start!=-1)
                            return false;
                        start=u;
                    }
                }
                std::list<int>::iterator tmp=it;
                ++tmp;
                if(tmp==ext.end()){
                    if(end!=-1)
                        return false;
                    end=u;
                }else{
                    if(mark[*tmp]!=marker){
                        if(end!=-1)
                            return false;
                        end=u;
                    }
                }
            }
            if(start==-1||end==-1)
                return false;
            for(int u:can){
                if(u!=start&&u!=end){
                    ext.erase(pos[u]);
                    pos[u]=std::list<int>::iterator();
                }
            }
            pos[v]=ext.insert(pos[end],v);
        }
        return true;
    }
    bool istri(const int&u,const int&v,const int&w){
        return to2[u].count(v)&&to2[v].count(w)&&to2[w].count(u);
    }
 
public:
    MaximalPlanarGraph(const int&_n):
        n(_n),to(n),to2(n),m(0),rt(n),invc(n),dec(n),dpos(n),pos(n),removed(n),mark(n){
    }
    void AddEdge(const int&u,const int&v){
        to[u-1].push_back(v-1);
        to[v-1].push_back(u-1);
        to2[u-1].insert(v-1);
        to2[v-1].insert(u-1);
        ++m;
    }
    bool Recognize(){
        if(n==1&&m==0)
            return true;
        if(n==2&&m==1)
            return true;
        if(n==3&&m==3)
            return true;
        if(n<=3)
            return false;
        if(m!=3*n-6)
            return false;
        int v1;
        for(v1=0;v1<n;++v1)
            if(to[v1].size()<3)
                return false;
        for(v1=0;v1<n;++v1)
            if(to[v1].size()<=5)
                break;
        if(v1>=n)
            return false;
        int v2=to[v1].back();
        for(int i=0;i+1<to[v1].size();++i){
            int vn=to[v1][i];
            if(istri(v1,v2,vn)){
                if(!order(v1,v2,vn))
                    continue;
                if(!embed())
                    continue;
                return true;
            }
        }
        return false;
    }
};
int main(){
    mc{
        int n,m;
        si(n),si(m);
        MaximalPlanarGraph g(n);
        lp(i,1,m){
            int u,v;
            si(u),si(v);
            g.AddEdge(u,v);
        }
        pf("%d\n",int(g.Recognize()));
    }
    return 0;
}


