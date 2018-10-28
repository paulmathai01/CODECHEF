#include <iostream>
#include<stdio.h>
#include<math.h>
#define mod 1000000007
using namespace std;
struct node
{
	double val;
	int pos;
};
int main() {
	// your code goes here
 long long int  i,j,ans[100005],a[100005],k,n,m;
 int min[100005];
 struct node pq[100005];
 scanf("%lld%lld",&n,&k);
for(i=1;i<=n;i++)
scanf("%lld",&a[i]);
ans[0]=1;
min[1]=0;
pq[1].val=log(a[1]);
pq[1].pos=1;
int start=1,end=1;
for(i=2;i<=n;i++)
{
	if(i-pq[start].pos>k)
	start++;
	min[i]=pq[start].pos;
	double temp=log(a[i])+pq[start].val;
	while(temp<pq[end].val && end>=start)
	end--;
	end++;
	pq[end].val=temp;
	pq[end].pos=i;
}
i=n;
ans[n]=1;
while(i>0)
{
	ans[n]=(ans[n]*a[i])%mod;
	i=min[i];
}
printf("%lld\n",ans[n]%mod);
	return 0;
}   