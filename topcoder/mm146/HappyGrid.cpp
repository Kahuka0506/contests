#include <bits/stdc++.h>
#define rep(i,n) for(int i=0; i<n; i++)
#define reps(i,s,n) for(int i=s; i<n; i++)
#define per(i,n) for(int i=n-1; i>=0; i--)
#define pers(i,n,s) for(int i=n-1; i>=s; i--)
#define all(v) v.begin(),v.end()
#define fi first
#define se second
#define pb push_back
#define si(v) int(v.size())
#define lb(v,n) lower_bound(all(v),n)
#define lbi(v,n) int(lower_bound(all(v),n) - v.begin())
#define ub(v,n) upper_bound(all(v),n)
#define ubi(v,n) int(upper_bound(all(v),n) - v.begin())
 
#define mod 1000000007
#define infi 1010000000
#define infl 1100000000000000000
 
#define outve(v) {for(auto i : v) cout << i << " ";cout << endl;}
#define outmat(v) for(auto i : v){for(auto j : i) cout << j << " ";cout << endl;}
#define errve(v) {for(auto i : v) cerr << i << " ";cerr << endl;}
#define errmat(v) for(auto i : v){for(auto j : i) cerr << j << " ";cerr << endl;}
#define in(n,v) for(int i=0; i<(n); i++){cin >> v[i];}
#define IN(n,m,v) rep(i,n) rep(j,m){cin >> v[i][j];}
#define cyes cout << "Yes" << endl
#define cno cout << "No" << endl
#define cYES cout << "YES" << endl
#define cNO cout << "NO" << endl
#define csp << " " <<
#define outset(n) cout << fixed << setprecision(n);
 
using namespace std;
using ll = long long;
using ull = unsigned long long;
using uint = unsigned int;
using ld = long double;
using vi = vector<int>;
using vvi = vector<vector<int>>;
using vd = vector<double>;
using vvd = vector<vector<double>>;
using vl = vector<ll>;
using vvl = vector<vector<ll>>;
using vs = vector<string>;
using pii = pair<int,int>;
using pll = pair<ll,ll>;
template<typename T> using ve = vector<T>;
template<typename T> using vv = vector<vector<T>>;
template<typename T> using pq2 = priority_queue<T>;
template<typename T> using pq1 = priority_queue<T,vector<T>,greater<T>>;
template<typename T> bool chmax(T &a, T b) {if(a < b) {a = b;return 1;}return 0;}
template<typename T> bool chmin(T &a, T b) {if(a > b) {a = b;return 1;}return 0;}

int popcnt(uint n) {return __builtin_popcount(n);}
int popcntl(ull n) {return __builtin_popcountll(n);}
int bsr(uint n) {return 31 - __builtin_clz(n);}
int bsrl(ull n) {return 63 - __builtin_clzll(n);}
int bsf(uint n) {return __builtin_ctz(n);}
int bsfl(ull n) {return __builtin_ctzll(n);}


unsigned int xorshift() {
    static unsigned int tx = 123456789, ty=362436069, tz=521288629, tw=88675123;
    unsigned int tt = (tx^(tx<<11));
    tx = ty; ty = tz; tz = tw;
    return ( tw=(tw^(tw>>19))^(tt^(tt>>8)) );
}

class TimerChrono {
    std::chrono::system_clock::time_point ti;
    double time;
public:
    TimerChrono(){}
    void start(){ti = std::chrono::system_clock::now();}
    inline double get_time(){
        return time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now() - ti).count() / 1000.0);
    }
};

const double TIME_LIMIT = 9000;
int N,C,K,P;
const int N_MAX = 31;
int grid[N_MAX][N_MAX];
const int dh[4] = {0,1,0,-1};
const int dw[4] = {1,0,-1,0};


int score_comp_K = 0;
int score_penalty = 0;

int cal_score1(int (&grid)[N_MAX][N_MAX], int h0, int w0, int h1, int w1, int turn=0){
    int vis[N] = {};
    int target_cell[10];
    int cn_target = 0;
    target_cell[cn_target] = h0*N+w0;
    cn_target++;
    target_cell[cn_target] = h1*N+w1;
    cn_target++;
    rep(k,4){
        int nh = h0 + dh[k], nw = w0 + dw[k];
        if(nh < 0 || nh >= N || nw < 0 || nw >= N || grid[nh][nw] == -1) continue;
        target_cell[cn_target] = nh*N+nw;
        cn_target++;
    }
    rep(k,4){
        int nh = h1 + dh[k], nw = w1 + dw[k];
        if(nh < 0 || nh >= N || nw < 0 || nw >= N || grid[nh][nw] == -1) continue;
        target_cell[cn_target] = nh*N+nw;
        cn_target++;
    }
    
    queue<int> que;
    rep(i,cn_target){
        int h_s = target_cell[i]/N;
        int w_s = target_cell[i]%N;
        if(vis[h_s] & 1 << w_s) continue;
        while(!que.empty()) que.pop();
        que.push(h_s*N+w_s);
        vis[h_s] |= 1<<w_s;
        int nn = 0;
        while (!que.empty()) {
            int hw = que.front();
            que.pop();
            int h = hw/N, w = hw%N;
            nn++;
            rep(k,4){
                int nh = h + dh[k], nw = w + dw[k];
                if(nh < 0 || nh >= N || nw < 0 || nw >= N) continue;
                if((vis[nh] & 1<<nw) || grid[nh][nw] != grid[h_s][w_s]) continue;
                vis[nh] |= 1<<nw;
                que.push(nh*N+nw);
            }
        }
        if(nn == K) score_comp_K-=P;
        else score_penalty -= (nn-K)*(nn-K);
    }
    swap(grid[h0][w0],grid[h1][w1]);
    
    rep(i,N) vis[i] = 0;
    rep(i,cn_target){
        int h_s = target_cell[i]/N;
        int w_s = target_cell[i]%N;
        if(vis[h_s] & 1<<w_s) continue;
        while(!que.empty()) que.pop();
        que.push(h_s*N+w_s);
        vis[h_s] |= 1<<w_s;
        int nn = 0;
        while (!que.empty()) {
            int hw = que.front();
            que.pop();
            int h = hw/N, w = hw%N;
            nn++;
            rep(k,4){
                int nh = h + dh[k], nw = w + dw[k];
                if(nh < 0 || nh >= N || nw < 0 || nw >= N) continue;
                if((vis[nh] & 1<<nw) || grid[nh][nw] != grid[h_s][w_s]) continue;
                vis[nh] |= 1<<nw;
                que.push(nh*N+nw);
            }
        }
        if(nn == K) score_comp_K+=P;
        else score_penalty += (nn-K)*(nn-K);
    }
    //cerr << 1 csp score_penalty csp score_comp_K << endl;
    return (score_comp_K-turn)*100-score_penalty*10;
}

int cal_score(int (&grid)[N_MAX][N_MAX], int turn=0, int real_score=0){
    score_comp_K = score_penalty = 0;
    
    int res = 0;
    int vis[N][N] = {};
    int cn = 1;
    int penalty = 0;
    queue<int> que;
    rep(i,N) rep(j,N){
        if(vis[i][j] > 0 || grid[i][j] == -1) continue;
        vis[i][j] = cn;
        while(!que.empty()) que.pop();
        que.push(i*N+j);
        int color = grid[i][j];
        int num = 0;
        while (!que.empty()) {
            int hw = que.front();
            que.pop();
            int h = hw/N, w = hw%N;
            num++;
            rep(k,4){
                int nh = h + dh[k], nw = w + dw[k];
                if(nh < 0 || nh >= N || nw < 0 || nw >= N) continue;
                if(vis[nh][nw] > 0 || grid[nh][nw] != color) continue;
                vis[nh][nw] = cn;
                que.push(nh*N+nw);
            }
        }
        if(num == K) res+=P, score_comp_K+=P;
        else penalty += (num-K)*(num-K), score_penalty += (num-K)*(num-K);
    }
    //cerr << score_penalty csp score_comp_K << endl;
    if(real_score == 1) return res-turn;
    return (res-turn)*100-penalty*10;
}
int cal_score_real(int (&grid)[N_MAX][N_MAX], int turn=0){
    return cal_score(grid,turn,1);
}



ve<pair<pii,pii>> create_action(int (&grid_init)[N_MAX][N_MAX],int (&grid_target)[N_MAX][N_MAX]){
    
    ve<pair<pii,pii>> res;
    int confirm[N][N] = {};
    
    int i = 0, j = 0;
    int direction = 0;
    int vis[N][N] = {};
    int from[N][N] = {};
    rep(i,N) rep(j,N) from[i][j] = -1, vis[i][j] = -1;
    queue<int> que;
    rep(ii,N*N) {
        if(grid_init[i][j] == grid_target[i][j]) {
            confirm[i][j] = 1;
        }else{
            while (!que.empty()) que.pop();
            que.push(i*N+j);
            from[i][j] = -1;
            vis[i][j] = ii;
            int target_hw = -1;
            while (!que.empty()) {
                int hw = que.front();
                que.pop();
                int h = hw/N, w = hw%N;
                rep(kk,4){
                    int k = (kk+direction+1)%4;
                    int nh = h + dh[k], nw = w + dw[k];
                    if(nh < 0 || nh >= N || nw < 0 || nw >= N) continue;
                    if(vis[nh][nw] == ii || grid_init[nh][nw] == -1) continue;
                    vis[nh][nw] = ii;
                    from[nh][nw] = h*N+w;
                    if(grid_init[nh][nw] == grid_target[i][j] && confirm[nh][nw] == 0){
                        target_hw = nh*N+nw;
                        break;
                    }
                    que.push(nh*N+nw);
                }
                if(target_hw != -1) break;
            }
            int s = target_hw;
            ve<pair<pii,pii>> an;
            int cn_confirm = 0;
            int cc = 0;
            while (from[s/N][s%N] != -1) {
                an.pb(make_pair(make_pair(s/N,s%N),make_pair(from[s/N][s%N]/N,from[s/N][s%N]%N)));
                s = from[s/N][s%N];
            }
            rep(k,si(an)) {
                res.pb(an[k]);
                swap(grid_init[an[k].fi.fi][an[k].fi.se],grid_init[an[k].se.fi][an[k].se.se]);
                if(confirm[an[k].fi.fi][an[k].fi.se] == 1) {
                    if(grid_init[an[k].fi.fi][an[k].fi.se] != grid_target[an[k].fi.fi][an[k].fi.se]) cn_confirm++;
                    cc++;
                }
            }
            if(cn_confirm > 0){
                for (int k = si(an)-2; k >= 0; k--) {
                    res.pb(an[k]);
                    swap(grid_init[an[k].fi.fi][an[k].fi.se],grid_init[an[k].se.fi][an[k].se.se]);
                    if(confirm[an[k].se.fi][an[k].se.se] == 1) {
                        cc--;
                        if(cc == 0) break;
                    }
                }
            }
            confirm[i][j] = 1;
        }
        
        //next i,j
        int ni = i + dh[direction], nj = j + dw[direction];
        if(ni < 0 || ni >= N || nj < 0 || nj >= N || confirm[ni][nj] == 1) {
            direction = (direction+1)%4;
            ni = i + dh[direction], nj = j + dw[direction];
            if(ni < 0 || ni >= N || nj < 0 || nj >= N || confirm[ni][nj] == 1) break;
        }
        i = ni, j = nj;
    }
    return res;
}

int create_action_cnt(int (&grid_init)[N_MAX][N_MAX], int (&grid_target)[N_MAX][N_MAX]){
    int res = 0;
    int confirm[N] = {};
    
    int i = 0, j = 0;
    int direction = 0;
    int vis[N][N] = {};
    int from[N][N] = {};
    queue<int> que;
    int an[N*N] = {};
    rep(ii,N*N) {
        if(grid_init[i][j] == grid_target[i][j]) {
            confirm[i] |= 1<<j;
        }else{
            while (!que.empty()) que.pop();
            que.push(i*N+j);
            from[i][j] = -1;
            vis[i][j] = ii+2;
            int target_hw = -1;
            while (!que.empty()) {
                int hw = que.front();
                que.pop();
                int h = hw/N, w = hw%N;
                rep(kk,4){
                    int k = (kk+direction+1)%4;
                    int nh = h + dh[k], nw = w + dw[k];
                    if(nh < 0 || nh >= N || nw < 0 || nw >= N) continue;
                    if(vis[nh][nw] == ii+2 || grid_init[nh][nw] == -1) continue;
                    vis[nh][nw] = ii+2;
                    from[nh][nw] = h*N+w;
                    if(grid_init[nh][nw] == grid_target[i][j] && !(confirm[nh]&1<<nw)){
                        target_hw = nh*N+nw;
                        break;
                    }
                    que.push(nh*N+nw);
                }
                if(target_hw != -1) break;
            }
            
            int s = target_hw;
            int cn_confirm = 0;
            int d = 0;
            int cc = 0;
            while (from[s/N][s%N] != -1) {
                an[d] = s*10000 + from[s/N][s%N];
                s = from[s/N][s%N];
                d++;
            }
            rep(k,d) {
                res++;
                int h = (an[k]/10000)/N, w = (an[k]/10000)%N;
                swap(grid_init[(an[k]/10000)/N][(an[k]/10000)%N],grid_init[(an[k]%10000)/N][(an[k]%10000)%N]);
                if((confirm[h]&1<<w)) {
                    if(grid_init[h][w] != grid_target[h][w]) cn_confirm++;
                    cc++;
                }
            }
            if(cn_confirm > 0){
                for (int k = d-2; k >= 0; k--) {
                    res++;
                    swap(grid_init[(an[k]/10000)/N][(an[k]/10000)%N],grid_init[(an[k]%10000)/N][(an[k]%10000)%N]);
                    if((confirm[(an[k]%10000)/N]&1<<((an[k]%10000)%N))) {
                        cc--;
                        if(cc == 0) break;
                    }
                }
            }
            confirm[i] |= 1<<j;
        }
        
        //next i,j
        int ni = i + dh[direction], nj = j + dw[direction];
        if(ni < 0 || ni >= N || nj < 0 || nj >= N || (confirm[ni]&1<<nj)) {
            direction = (direction+1)%4;
            ni = i + dh[direction], nj = j + dw[direction];
            if(ni < 0 || ni >= N || nj < 0 || nj >= N || (confirm[ni]&1<<nj)) break;
        }
        i = ni, j = nj;
    }
    return res;
}

ve<pair<pii,pii>> optimize_action(int (grid_init)[N_MAX][N_MAX], int (&grid)[N_MAX][N_MAX], double use_time){
    
    ve<pair<pii,pii>> ans_first;
    int score_swap = 0, score_swap0 = 0;

    TimerChrono timer1;
    timer1.start();
    
    int grid_init_[N_MAX][N_MAX];
    rep(i,N) rep(j,N) grid_init_[i][j] = grid_init[i][j];
    score_swap = create_action_cnt(grid_init_, grid);
    score_swap0 = score_swap;
    
    int loop_ = 0;
    while (timer1.get_time() < TIME_LIMIT*use_time) {
        loop_++;
        
        rep(i,N) rep(j,N) grid_init_[i][j] = grid_init[i][j];
        
        int n = -1;
        if(xorshift()%2==0 && si(ans_first) > 0) n = xorshift()%si(ans_first);
        pair<pii,pii> ans_first_n0;
        
        int h = xorshift()%(N-1), w = xorshift()%(N-1), k = xorshift()%2;
        int nh = h + dh[k], nw = w + dw[k];
        while (grid_init_[h][w] == -1 || grid_init_[h+dh[k]][w+dw[k]] == -1 || grid_init_[h][w]==grid_init_[nh][nw]) {
            h = xorshift()%(N-1), w = xorshift()%(N-1), k = xorshift()%2;
            nh = h + dh[k], nw = w + dw[k];
        }
        if(n == -1) ans_first.pb(make_pair(make_pair(h,w),make_pair(nh,nw)));
        else {
            ans_first_n0 = ans_first[n];
            ans_first[n] = make_pair(make_pair(h,w),make_pair(nh,nw));
        }
        
        rep(i,si(ans_first)) swap(grid_init_[ans_first[i].fi.fi][ans_first[i].fi.se],grid_init_[ans_first[i].se.fi][ans_first[i].se.se]);
        
        
        int score_swap1 = create_action_cnt(grid_init_, grid)+si(ans_first);
        if(score_swap1 < score_swap) {
            score_swap = score_swap1;
        }else{
            if(n == -1) ans_first.pop_back();
            else ans_first[n] = ans_first_n0;
        }
    }
    //cerr << loop_ csp score_swap0 csp score_swap csp si(ans_first) << endl;
    if(score_swap0 <= score_swap) ans_first.clear();
    
    ve<pair<pii,pii>> ans;
    rep(i,N) rep(j,N) grid_init_[i][j] = grid_init[i][j];
    rep(i,si(ans_first)){
        swap(grid_init_[ans_first[i].fi.fi][ans_first[i].fi.se],grid_init_[ans_first[i].se.fi][ans_first[i].se.se]);
        ans.pb(ans_first[i]);
    }
    ve<pair<pii,pii>> ans_second = create_action(grid_init_, grid);
    rep(i,si(ans_second)) ans.pb(ans_second[i]);
    
    return ans;
}

struct Component {
    int color;
    int cells[8];
    Component(int color_, int (&cells_)[8]):color(color_){
        rep(i,K) cells[i] = cells_[i];
    }
};

void swap_component_color(int (grid_init)[N_MAX][N_MAX], int (&grid)[N_MAX][N_MAX], double use_time){
    TimerChrono timer;
    timer.start();
    
    ve<Component> components;
    
    int vis[N_MAX][N_MAX] = {};
    queue<int> que;
    rep(i,N) rep(j,N){
        if(vis[i][j] > 0 || grid[i][j] == -1) continue;
        vis[i][j] = 1;
        while(!que.empty()) que.pop();
        que.push(i*N+j);
        int color = grid[i][j];
        int num = 0;
        vi cells_;
        while (!que.empty()) {
            int hw = que.front();
            que.pop();
            cells_.pb(hw);
            int h = hw/N, w = hw%N;
            num++;
            rep(k,4){
                int nh = h + dh[k], nw = w + dw[k];
                if(nh < 0 || nh >= N || nw < 0 || nw >= N) continue;
                if(vis[nh][nw] > 0 || grid[nh][nw] != color) continue;
                vis[nh][nw] = 1;
                que.push(nh*N+nw);
            }
        }
        if(num == K) {
            int cells[8] = {};
            rep(i,K) cells[i] = cells_[i];
            components.pb(Component(color, cells));
        }
    }
    
    
    int score = cal_score(grid,0);
    int grid_init_[N_MAX][N_MAX];
    rep(i,N) rep(j,N) grid_init_[i][j] = grid_init[i][j];
    score -= create_action_cnt(grid_init_, grid)*100;
    //cerr << score << endl;
    
    int loop = 0;
    while (timer.get_time() < TIME_LIMIT*use_time) {
        loop++;
        
        int n = xorshift()%si(components);
        int m = xorshift()%(si(components)-1);
        if(n == m) m++;
        while (components[n].color == components[m].color) {
            n = xorshift()%si(components);
            m = xorshift()%(si(components)-1);
            if(n == m) m++;
        }
        
        rep(i,K){
            int a = components[n].cells[i], b = components[m].cells[i];
            swap(grid[a/N][a%N],grid[b/N][b%N]);
        }
        
        int score1 = cal_score(grid,0);
        rep(i,N) rep(j,N) grid_init_[i][j] = grid_init[i][j];
        score1 -= create_action_cnt(grid_init_, grid)*100;
        
        if(score < score1){
            score = score1;
        }else{
            rep(i,K){
                int a = components[n].cells[i], b = components[m].cells[i];
                swap(grid[a/N][a%N],grid[b/N][b%N]);
            }
        }
        //if(loop%100 == 0) cerr << loop csp score csp score1 csp timer.get_time() << endl;
    }
    //cerr << score << endl;
    return;
}

void solve(){
    
    TimerChrono timer;
    timer.start();
    
    cin >> N >> C >> K >> P;
    int wall_num = 0;
    rep(i,N) rep(j,N){
        int x;
        cin >> x;
        x--;
        grid[i][j] = x;
        if(x == -1) wall_num++;
    }
    int grid_init[N_MAX][N_MAX];
    int grid_init_[N_MAX][N_MAX];
    rep(i,N) rep(j,N) grid_init[i][j] = grid[i][j];
    double W = (double)wall_num/(double)(N*N);
    
    int best_score = 0;
    int score = 0;
    int grid_best[N_MAX][N_MAX];
    ve<pair<pii,pii>> ans;
    
    //FILE* fp;
    //fp = fopen("data_sa.txt","w");
    
    int loop = 0;
    double T1 = 200, T0 = 0;
    if(N <= 13) T1 = 500;
    else if(N <= 20) T1 = 600;
    else if(N <= 25) T1 = 200;
    double T = T1;
    
    
    double P_logis = 1;
    //if((double)P >= -3.7*(double)N+135.0 || (K==8&&N>=25) || (W<=0.06&&N>=25)) P_logis = 0;
    if((N >= 11&&11<=(P/K)) || (N >= 23&&7<=(P/K)) || (N >= 15&&9<=(P/K)) || (K==8&&N>=25) || (W<=0.06&&N>=25)) P_logis = 0;
    
    score = cal_score(grid,0);
    if(P_logis >= 0.5){
        rep(i,N) rep(j,N) grid_init_[i][j] = grid_init[i][j];
        score -= create_action_cnt(grid_init_, grid)*100;
    }
    
    double time = timer.get_time();
    while (time < TIME_LIMIT*0.85) {
        loop++;
        
        int score_comp_K0 = score_comp_K;
        int score_penalty0 = score_penalty;
        
        int h = xorshift()%(N-1), w = xorshift()%(N-1);
        int d_h = xorshift()%4*(xorshift()%2?-1:1);
        int d_w = xorshift()%4*(xorshift()%2?-1:1);
        int nh = h + d_h, nw = w + d_w;
        chmax(nh,0);chmax(nw,0);chmin(nh,N-1);chmin(nw,N-1);
        while ((h==nh&&w==nw) || grid[h][w]==-1 || grid[nh][nw] == -1 || grid[h][w]==grid[nh][nw]) {
            h = xorshift()%N, w = xorshift()%N;
            d_h = xorshift()%4*(xorshift()%2?-1:1);
            d_w = xorshift()%4*(xorshift()%2?-1:1);
            nh = h + d_h, nw = w + d_w;
            chmax(nh,0);chmax(nw,0);chmin(nh,N-1);chmin(nw,N-1);
        }
        
        int score1 = cal_score1(grid,h,w,nh,nw,0);

        if(P_logis >= 0.5){
            rep(i,N) rep(j,N) grid_init_[i][j] = grid_init[i][j];
            score1 -= create_action_cnt(grid_init_, grid)*100;
        }
        
        if(score < score1){
            score = score1;
        }else{
            //T = 200.0/(1.0 + exp(-(0.40-timer.get_time()/(TIME_LIMIT*0.85))*8.0));
            T = T1 + (T0-T1)*time/(TIME_LIMIT*0.85);
            
            if((double)(xorshift()%infi)/(double)infi < exp((score1-score)/T)){
                score = score1;
            }else{
                score_comp_K = score_comp_K0;
                score_penalty = score_penalty0;
                swap(grid[h][w],grid[nh][nw]);
                
            }
        }
        if(chmax(best_score,score)){
            rep(i,N) rep(j,N) grid_best[i][j] = grid[i][j];
        }
        time = timer.get_time();
        //if(loop%100 == 0) cerr << loop csp score csp score1 csp timer.get_time() csp T << endl;
        //if(loop%100 == 0) fprintf(fp,"%d %.5f %d %d %.5f \n",loop,timer.get_time(),score,score1,T);
    }
    //fclose(fp);
    
    swap_component_color(grid_init,grid_best, 0.1);
    ans = optimize_action(grid_init, grid_best, 0.1);
    
    cerr << N csp C csp K csp P csp wall_num csp loop csp best_score-si(ans) csp cal_score_real(grid_best,si(ans)) csp si(ans) csp timer.get_time() << endl;
    
    cout << si(ans) << endl;
    rep(i,si(ans)) cout << ans[i].fi.fi csp ans[i].fi.se csp ans[i].se.fi csp ans[i].se.se << endl;
    
}



int main(){
    
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    solve();

    return 0;
}

