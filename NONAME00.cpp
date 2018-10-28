#include <iostream.h>
int main()
{int a[50];
 int n,i,j,k,size,count;
 k=0;
 cin>>n;
 cin>>size;
 while(k<n)
 {for (i=0;i<size;i++)
  {for(j=0;j<size;j++)
   {if (i==j)
    continue;
    else if(a[i]==a[j])
    count++;
   }
  }
 cout<<count;
 }
