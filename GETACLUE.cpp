#include <vector>
#include <list>
#include <cassert>
#include <sstream>
#include <map>
#include <set>
#include <climits>
#include <deque>
#include <fstream>
#include <stack>
#include <bitset>
#include <stack>
#include <queue>
#include <algorithm>
#include <functional>
#include <numeric>
#include <cstring>
#include <utility>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
using namespace std;
 
template<class A, class B> A cvt(B x) {stringstream s;s<<x;A r;s>>r;return r;}
 
#define FOR(i,a,b) for(int i= (int )a ; i < (int )b ; ++i) 
#define REV(i,a,b) for(int i= (int )a ; i >= (int)b ; --i)
#define REP(i,n) FOR(i,0,n)
#define DEP(i,n) REV(i,n,0)
#define PB push_back
#define PP pop()
#define EM empty()
#define INF 1000000000
#define PF push_front
#define ALL(x) x.begin(),x.end()
#define SORT(x) sort(ALL(x))
#define V(x) vector< x >
#define Debug false
#define PRINT(x)        cout << #x << " " << x << endl
#define LET(x,a) 	    __typeof(a) x(a)
#define IFOR(i,a,b) 	for(LET(i,a);i!=(b);++i)
#define EACH(it,v)  	IFOR(it,v.begin(),v.end())
#define PRESENT(c,x) 	((c).find(x) != (c).end())
#define SZ(x) 		x.size()
#define CPRESENT(c,x) 	(find(c.begin(),c.end(),x) != (c).end())
#define D(N) 		int N
#define S(N)		scanf("%d",&N)
#define FASTIO          1
 
typedef pair<int,int>   PI;
typedef pair<int,PI>    TRI;
typedef V( int )        VI;
typedef V( PI  )        VII;
typedef V( string )     VS;
typedef long long       LL;
typedef long double     LD;
 
/* FastIO, generally required these days ;) */
 
#ifndef FASTIO
char *ipos, *opos, InpFile[20000000], OutFile[20000000], DIP[20];
inline int input(int flag=0) {
 
   while(*ipos <= 32) ++ipos;
   if ( flag  ) return (*ipos++ - '0'); /* For getting Boolean Characters */
   int x=0, neg = 0;char c;
   while( true ) {
      c=*ipos++; if(c == '-') neg = 1;
      else {
	 if (c<=32) return neg?-x:x;
	 x=(x<<1)+(x<<3)+c-'0';
      }
   }
}
inline void output(int x,int flag) {
   int y,dig=0;
   while (x||!dig) { y=x/10;DIP[dig++]=x-((y << 3) + (y << 1))+'0';x=y;}
   while (dig--) *opos++=DIP[dig];
   *opos++= flag ? '\n' : ' ';
}
inline void InitFASTIO() {
   ipos = InpFile; opos = OutFile;
   fread_unlocked(InpFile,20000000,1,stdin);
}
inline void FlushFASTIO() {
   fwrite_unlocked(OutFile,opos-OutFile,1,stdout);	
}
#endif
 
/* Main Code Starts from here */
 
#define Max 100
const int N = 20;
 
int sum[2][Max][Max][3], size[2][Max], rank[2][N][N][N];
int Input[10][3], Given[10], Iorders[10];
//int list[10][Max*10], start[10], end[10];
 
int numSolutions[10][Max], totalSolutions[10][Max][Max];
set<int> forwhich[10][Max];
bool isStillValid[10][Max];
int validSolutions[10][Max*10], ssize[10];
 
 
int tmpSolution[3], visited[Max];
 
void init() {
   memset(visited,0,sizeof visited);
}
bool isEmpty(int j) {
   return (ssize[j] == 0);
}
bool allEmpty() {
   REP(j,10) if(!isEmpty(j)) return false;
   return true;
}
void makeEmpty() {
   memset(ssize,0,sizeof ssize);
}
bool isSet(VI& present) {
   REP(i,3) if(visited[present[i]]) return true;
   return false;
}
void setVal(VI& present,int which,int total,int index)  {
   REP(l,3) present[l] = sum[which][total][index][l];
}
void Set(VI& present) {
   REP(i,3) visited[present[i]] = 1;
}
void unSet(VI& present) {
   REP(i,3) visited[present[i]] = 0;
}
 
 
void pre() {
 
   int t = 0, size1, size2;
   FOR(i,1,19) FOR(j,i+1,19) FOR(k,j+1,19) {
 
      t = i + j + k;
      size1 = size[0][t]; size2 = size[1][t];
 
      sum[0][t][size1][0] = i; 
      sum[0][t][size1][1] = j; 
      sum[0][t][size1][2] = k;
      rank[0][i][j][k] = size[0][t]++;
 
      if(i < 7 && j > 6 && j < 13 && k > 12 && k < 19) {
	 sum[1][t][size2][0] = i; 
	 sum[1][t][size2][1] = j; 
	 sum[1][t][size2][2] = k;
	 rank[1][i][j][k] = size[1][t]++;
      }
   }
}
 
 
void undoit(int index,int oldI, int sf, VI& track, int order) {
   if(index > 5 || order < 0) return;
   if(index == 5) {
      REP(i,5) {
	 totalSolutions[i][track[i]][order] -= 1;
	 if(!totalSolutions[i][track[i]][order]) {
	    // This means we have eliminated one of the possible combination for that player and combination.
	    numSolutions[i][track[i]] -= 1;
	    if(numSolutions[i][track[i]] == 1) {
	       //forwhich[i][track[i]].erase(order); // This tell us that pair of [Player,Combination] has a unique solution [5,order].
	       validSolutions[i][ssize[i]++] = track[i]; // Put it into the ith player's valid regime.
	       if(Debug) printf("I inserted this %dth index for player %d\n",track[i],i);
	    }
	    forwhich[i][track[i]].erase(order); // This tell us that pair of [Player,Combination] has a unique solution [5,order].
	 }
      }
   }
   else {
      assert(index <= 4);VI present(5,0);
      int start = index == oldI ? sf : 0, end = index == oldI ? sf+1 : size[0][Given[index]];
      FOR(i,start,end) {
	 track[index] = i; setVal(present,0,Given[index],i); 
	 if(!isSet(present) && isStillValid[index][i]) {
	    Set(present);undoit(index+1,oldI,sf,track,order);
	    unSet(present);
	 }
      }
   }
}
void doit(int index,VI& track, int order) {
   if(index > 5 || order < 0) return;
   if(index == 5) {
      REP(i,5) {
	 if(!totalSolutions[i][track[i]][order]) {
	    numSolutions[i][track[i]] += 1; // Number of Solution for same index, track[index] but different Combination of Given[5].
	    forwhich[i][track[i]].insert(order);
	    // If At the end, this is one, then we can add it to possible solutions.
	    // It means that if the Player has this combination then, there can be only one possible solution. And hence he knows the Solution, and hence wins.
	 }
	 totalSolutions[i][track[i]][order] += 1; // Number of solutions for same combination of Given[5].
      }
   }
   else {
      assert(index <= 4); VI present(5,0);
      REP(i,size[0][Given[index]]) {
	 track[index] = i; setVal(present,0,Given[index],i);
	 //if(Debug) REP(j,3) printf("i == %d, index == %d, %d%c",i,index,present[j],j == 2 ? 10 : 32);
	 if(!isSet(present)) {
	    Set(present);doit(index+1,track,order);
	    unSet(present);
	 }
      }
   }
}
 
void generateInput() {
 
   FILE *f = freopen("in.txt","w",stdout);
 
   srand(time(NULL));
   int kases = 100; cout << kases << endl;
   set<int> t;
   while(kases--) {
 
      t.clear();
      int x = rand()%6+1, y = rand()%6+7, z = rand()%6+13;
      t.insert(x); t.insert(y); t.insert(z);
 
      vector<TRI> tmp;
      REP(i,5) {
	 do { x = rand()%18 + 1; } while(t.count(x)); t.insert(x);
	 do { y = rand()%18 + 1; } while(t.count(y)); t.insert(y);
	 do { z = rand()%18 + 1; } while(t.count(z)); t.insert(z);
	 tmp.PB(TRI(x,PI(y,z)));
      }
      int a[] = {0,1,2,3,4}; random_shuffle(a,a+5);
      REP(i,5) {
	 printf("%d %d %d\n",tmp[a[i]].first,tmp[a[i]].second.first, tmp[a[i]].second.second);
      }
      puts("");
   }
 
   fclose(f);
}
 
 
int main() {
 
 
   //double stime = clock();
 
   //generateInput();
   //freopen("in.txt","r",stdin);
   //freopen("out.txt","w",stdout);
 
 
   pre();
   int kases; scanf("%d",&kases);
   while(kases--) {
 
      int total = 0, winner = -1;
      REP(i,5) {
	 Given[i] = 0;
	 REP(j,3) {
	    scanf("%d",&Input[i][j]); 
	    total += Input[i][j]; Given[i] += Input[i][j];
	    // Given[i] --> total sum for the ith player.
	 }
	 sort(Input[i],Input[i]+3);
	 Iorders[i] = rank[0][Input[i][0]][Input[i][1]][Input[i][2]]; // Represents the rank of the given input combination in the precomputed one.
      }
      Given[5] = 171 - total; // This is the sum of the three cards kept before starting the game.
 
      if(Debug) {
	 REP(i,5) printf("The Ranks -- for %d it is %d\n",i+1,Iorders[i]);
	 REP(i,6) printf("Initial sum for %d is -- %d\n",i+1,Given[i]);
      }
 
      // Get all possible combinations for Given[5] == Required Sum.
      // For every Player, and a corresponsing triplet, get how many solutions are possible
      // If only one when we need to consider that else not.
      // This will form the initial construction part.
 
      memset(visited,0,sizeof visited); // For, browsing through valid possibilites.
      memset(numSolutions,0, sizeof numSolutions); // numSolutions --> Number of Possibilites [of three cards kept] for the pair {Player, Index}.
      memset(totalSolutions,0, sizeof totalSolutions); // total number of possibilities for triplet {Player, Index1, [Three cards's index, Index2]}.
      memset(isStillValid,-1,sizeof isStillValid); // Is [Player, Index] valid ? That is can Player have Index Combination ?
      //memset(forwhich,-1,sizeof forwhich); // For which [Three Card's Index] is numSolutions[player][index] = 1.
      REP(i,10) REP(j,Max) forwhich[i][j].clear();
 
 
      VI track(5,0), present(5,0); // Keeps Track of the index of a particular player while browsing through the solution.
      REP(i,size[1][Given[5]]) {
	 // For all Possible Combinations of Sum Given[5].
	 setVal(present,1,Given[5],i); Set(present);
	 doit(0,track,i); 
	 unSet(present);
      }
 
      // Now, we have what everyone knows at first hand.
      init();
 
      // What are the possible solutions ?
      makeEmpty(); // validSolutions[i] // Contians the valid solutions for each and every person.
 
      REP(i,5) {
	 if(Debug) {
	    printf("The initial valid indexes for player %d are - ", i+1);
	 }
	 REP(j,size[0][Given[i]]) {
	    if(numSolutions[i][j] == 1) {
	       // Means if one of these is the combination, then we have a unique answer.
	       // These are the only ones which are valid.
	       validSolutions[i][ssize[i]++] = j;
	       if(Debug) {
		  printf("%d ",j);
	       }
	    }
	 }
	 if(Debug) puts("");
      }
 
 
      while(!allEmpty() && winner == -1) {
 
	 // We start the game from here.
	 REP(player,5) {
	    // Now, Start the Game.
	    while(!isEmpty(player)) {
	       // We have some to invalidate, or something to verify on.
	       int top = validSolutions[player][--ssize[player]]; 
 
	       // Do we have the combination top ?
	       if(top == Iorders[player]) {
		  // So This player wins.
		  winner = player; makeEmpty();
		  break;
	       }
 
	       // Else, invalidate various other solutions, add if required.
	       // Now Validate, and Unvalidate things.
	       
	       int order = forwhich[player][top].size() ? *(forwhich[player][top].begin()) : -1; // For which combination of the before kept three cards ?
 
	       if(Debug) {
		  printf("Top is %d, and Order is %d\n",top,order);
	       }
	       setVal(present,1,Given[5],order); Set(present);
	       undoit(0,player,top,track,order);
	       unSet(present);
 
	       isStillValid[player][top] = false; // No, player can't have the combination with index top.
	    }
	 }
      }
 
      if(winner != -1) ++winner;
      printf("%d\n",winner);
   }
 
   //double etime = clock();
   //cerr << (etime-stime)/(CLOCKS_PER_SEC) << endl;
   
   return 0;
}