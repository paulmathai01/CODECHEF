#include <cstdio>
#include <cmath>
#include <cctype>
#include <cassert>
#include <algorithm>
#include <map>
#include <queue>
using namespace std;
 
int s[5], len;
 
long long solve(){
	if(len==1)
		return s[0];
	int B=1;
	for(int i=0; i<len; i++)
		B=max(B, s[i]+1);
	long long ans=0;
 
	for(;; B++){
		long long num=s[0];
		for(int i=1; i<len; i++)
			num=num*B+s[i];
		if(ans && num>=ans)
			break;
		long long check=num;
		for(int p=2; p<B; p++)
			while(!(check%p))
				check/=p;
		if(check!=1)
			continue;
		
		map<long long, long long> m;
		priority_queue<long long> q;
		m[num]=0;
		for(q.push(num); !q.empty(); q.pop()){
			long long n=q.top();
			for(int p=2; p<B; p++){
				if(!(n%p)){
					long long d=n/p;
					long long cost=p+B*m[n];
					if(ans && ans<=cost)
						continue;
					if(!m[d]){
						m[d]=cost;
						q.push(d);
					}else if(m[d]>cost)
						m[d]=cost;
				}
			}
		}
		if(m[1] && (!ans || ans>m[1])){
			ans=m[1];
		}
	}
	return ans;
}
 
int main(){
	int T;
	for(scanf("%d", &T); T--; ){
		char S[6];
		scanf("%s", S);
		for(len=0; S[len]; len++){
			s[len]=isdigit(S[len]) ? S[len]-'0' : S[len]-'A'+10;
		}
		printf("%lld\n", solve());
	}
	return 0;
}