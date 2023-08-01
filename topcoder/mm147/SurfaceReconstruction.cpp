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


double TIME_LIMIT = 9800.0;
const int LOCAL = 0;
const int TESTER = 0;
const int FILE_OUT = 0;


class SVGWriter {
    string filename;
    FILE *fp;
    int Heigth, Width;
    int cnt = 0;
    int fin = 1;
public:
    SVGWriter(string filename_, int H_=600, int W_=600):filename(filename_),Heigth(H_),Width(W_){
        fp = fopen(filename.c_str(),"w");
    }
    void newSVG(){
        fin = 0;
        fprintf(fp,"<svg height=\"%d\" id=\"f%d\" viewBox=\"-5 -5 %d %d\" width=\"%d\" xmlns=\"http://www.w3.org/2000/svg\">\n",Heigth,cnt,Heigth,Width,Width);
        fprintf(fp,"<rect fill=\"white\" height=\"810\" width=\"810\" x=\"-5\" y=\"-5\"/>\n");
    }
    void clear(){
        fprintf(fp,"</svg>");
        fin = 1;
        cnt++;
    }
    void path(int x0, int y0, int x1, int y1, int width=3, string color="lightgray"){
        fprintf(fp,"<path d=\"M%d,%d L%d,%d\" stroke-width=\"%d\" stroke=\"%s\"/>\n",x0,y0,x1,y1,width,color.c_str());
    }
    void circle(int x, int y, int r, string color="lightgray"){
        fprintf(fp,"<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\"/>\n", x,y,r, color.c_str());
    }
    void rect(int x, int y, int height, int width, string color="lightgray"){
        fprintf(fp,"<rect x=\"%d\" y=\"%d\" height=\"%d\" width=\"%d\" fill=\"%s\" />",x,y,height,width,color.c_str());
    }
    void text(int x, int y, int font_size, string s, string color="black"){
        fprintf(fp,"<text x=\"%d\" y=\"%d\" font-family=\"Verdana\" font-size=\"%d\" fill=\"%s\">", x,y,font_size,color.c_str());
        fprintf(fp, "%s", s.c_str());
        fprintf(fp, "</text>");
    }
    void save(){
        if(fin == 0) fprintf(fp,"</svg>");
        fclose(fp);
    }
};






const int dh[4] = {-1,1,0,0};
const int dw[4] = {0,0,-1,1};
const int dh_9[9] = {-1,-1,-1,0,0,0,1,1,1};
const int dw_9[9] = {-1,0,1,-1,0,1,-1,0,1};
const int N_MAX = 60;
const int NN_MAX = 3600;
int N,E;
double A;

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


double exp_table[730000];
void init_exp_table(){
    rep(i,730000) exp_table[i] = exp(-((double)i)/1000);
}
double sin_table[500000];
const double PI=3.14159265358979323846;
void init_sin_table(){
    rep(i,500000) sin_table[i] = sin((double)(i-250000)/10000.0);
}




struct Param {
    int type;
    int offh, offw;
    double s1, s2;
    double amp;
    Param();
    Param(int type_, int offh_, int offw_, double s1_=0.0, double s2_=0.0, double amp_=0.0):
    type(type_),offh(offh_),offw(offw_),s1(s1_),s2(s2_),amp(amp_){}
};

void apply(double (&grid)[NN_MAX], Param &p){
    rep(h,N) rep(w,N){
        double dz = 0;
        if(p.type == 0) dz = (h-p.offh)*p.s1 - (w-p.offw)*p.s2;
        else if(p.type == 1) dz = p.amp*(sin((w-p.offw)*p.s1) + cos((h-p.offh)*p.s2));
        else if(p.type == 2) dz = (h-p.offh)^(w-p.offw);
        else dz = p.amp*exp(-((h-p.offh)*(h-p.offh)*p.s1 + (w-p.offw)*(w-p.offw)*p.s2));
        grid[h*N+w] += dz;
    }
}
void normalize(double (&grid)[NN_MAX], int (&res)[NN_MAX],double grid_max=-infi, double grid_min=infi){
    if(grid_max < -100000000.0){
        rep(i,N*N) {
            chmin(grid_min,grid[i]);
            chmax(grid_max,grid[i]);
            res[i] = 0;
        }
    }
    grid_max -= grid_min;
    rep(i,N*N) res[i] = int(255*(grid[i]-grid_min)/grid_max);
}
Param rand_param(mt19937 &engine);
void init(mt19937 &engine, int &N, double &A, int &E, int (&grid)[NN_MAX]);
void print_svg_grid(int (&grid)[NN_MAX], int f, int (&checked)[NN_MAX], string name);


int cal_score(int (&pred)[NN_MAX], int (&exact)[NN_MAX], int query){
    double res1 = 0;
    rep(i,N) rep(j,N) res1 += (pred[i*N+j]-exact[i*N+j])*(pred[i*N+j]-exact[i*N+j]);
    res1 /= (double)(N*N);
    res1 /= (64.0*64.0);
    double res2 = (double)query/(N*N/9);
    return int(100000.0*(A*res1 + (1.0-A)*res2));
}





void updata_grid(ve<Param> &params, int (&grid)[NN_MAX]){
    double grid0[NN_MAX] = {};
    rep(i,si(params)) apply(grid0,params[i]);
    normalize(grid0,grid);
}

int cal_mse_score(vi &checked_h, vi &checked_w, vi &checked_grid, int (&grid)[NN_MAX]){
    int res = 0;
    int h,w;
    rep(i,si(checked_h)){
        h = checked_h[i], w = checked_w[i];
        res += (grid[h*N+w]-checked_grid[i])*(grid[h*N+w]-checked_grid[i]);
    }
    return res;
}


void change_param_neigbor(Param &param){
    if(param.type == 0){//slope
        int p = xorshift()%100;
        if(p <= 20){
            param.offh += xorshift()%9-4;
            param.offw += xorshift()%9-4;
            chmax(param.offh,0),chmin(param.offh,N-1);
            chmax(param.offw,0),chmin(param.offw,N-1);
        }else{
            param.s1 += (double)(xorshift())/(double)(UINT_MAX)*2.0-1.0;
            param.s2 += (double)(xorshift())/(double)(UINT_MAX)*2.0-1.0;
            chmax(param.s1,-14.0),chmin(param.s1,14.0);
            chmax(param.s2,-14.0),chmin(param.s2,14.0);
        }
        
    }else if(param.type == 1){//sincos
        int p = xorshift()%100;
        if(p <= 15){
            param.offh += xorshift()%9-4;
            param.offw += xorshift()%9-4;
            chmax(param.offh,0),chmin(param.offh,N-1);
            chmax(param.offw,0),chmin(param.offw,N-1);
        }else if(p <= 60){
            param.s1 += (double)(xorshift())/(double)(UINT_MAX)*0.08-0.04;
            param.s2 += (double)(xorshift())/(double)(UINT_MAX)*0.08-0.04;
            chmax(param.s1,-0.4),chmin(param.s1,0.4);
            chmax(param.s2,-0.4),chmin(param.s2,0.4);
        }else{
            param.amp += (double)(xorshift())/(double)(UINT_MAX)*20.0-10.0;
            chmax(param.amp,0.0),chmin(param.amp,100.0);
        }
        
    }else if(param.type == 2){//xor
        param.offh += xorshift()%9-4;
        param.offw += xorshift()%9-4;
        chmax(param.offh,0),chmin(param.offh,N-1);
        chmax(param.offw,0),chmin(param.offw,N-1);

    }else{//xy^2
        int p = xorshift()%100;
        if(p <= 15){
            param.offh += xorshift()%9-4;
            param.offw += xorshift()%9-4;
            chmax(param.offh,0),chmin(param.offh,N-1);
            chmax(param.offw,0),chmin(param.offw,N-1);
        }else if(p <= 60){
            param.s1 += (double)(xorshift())/(double)(UINT_MAX)*0.02-0.01;
            param.s2 += (double)(xorshift())/(double)(UINT_MAX)*0.02-0.01;
            chmax(param.s1,0.001),chmin(param.s1,0.1);
            chmax(param.s2,0.001),chmin(param.s2,0.1);
        }else{
            param.amp += (double)(xorshift())/(double)(UINT_MAX)*20.0-10.0;
            chmax(param.amp,0.0),chmin(param.amp,100.0);
        }
    }
}

void SA(vi &checked_h, vi &checked_w, vi &checked_grid, int (&surface)[NN_MAX], double time_limit=TIME_LIMIT){
    TimerChrono timer;
    timer.start();
    
    int grid[NN_MAX] = {};
    ve<Param> params;
    rep(j,3){
        rep(i,10) {
            Param p(min(i,3),xorshift()%N,xorshift()%N,0,0,0);
            params.pb(p);
        }
    }
    
    double grid_orig[N*N] = {};
    rep(i,si(params)) {
        Param &p = params[i];
        rep(h,N) rep(w,N){
            double dz = 0;
            if(p.type == 0) dz = (h-p.offh)*p.s1 - (w-p.offw)*p.s2;
            else if(p.type == 1) dz = p.amp*(sin((w-p.offw)*p.s1) + cos((h-p.offh)*p.s2));
            else if(p.type == 2) dz = (h-p.offh)^(w-p.offw);
            else dz = p.amp*exp(-((h-p.offh)*(h-p.offh)*p.s1 + (w-p.offw)*(w-p.offw)*p.s2));
            grid_orig[h*N+w] += dz;
        }
    }
    int score = infi;
    int best_score = infi;
    ve<Param> best_params;
    
    FILE* fp;
    if(LOCAL == 1 && TESTER == 0) fp = fopen("data_sa.txt","w");
    int loop = 0;
    double T1 = 600, T0 = 80;
    double T = T1;
    double time = timer.get_time();
    
    Param param0(0,0,0,0,0);
    double grid_orig0[N*N];
    while(time < time_limit){
        loop++;
        int n = xorshift()%30;
        
        param0 = params[n];
        rep(i,N*N) grid_orig0[i] = grid_orig[i];
        
        change_param_neigbor(params[n]);
        
        double min_grid = infi, max_grid = -infi;
        Param &p_new = params[n];
        Param &p_old = param0;
        {
            max_grid = -infi, min_grid = infi;
            int c = 0;
            rep(h,N) rep(w,N){
                double dz = 0;
                if(p_new.type == 0) {
                    dz = ((h-p_new.offh)*p_new.s1 - (w-p_new.offw)*p_new.s2) - ((h-p_old.offh)*p_old.s1 - (w-p_old.offw)*p_old.s2);
                }else if(p_new.type == 1) {
                    dz = p_new.amp*(sin_table[int(((w-p_new.offw)*p_new.s1)*10000)+250000]+ sin_table[int(((h-p_new.offh)*p_new.s2+PI/2)*10000)+250000])
                    - p_old.amp*(sin_table[int((w-p_old.offw)*p_old.s1*10000)+250000]+ sin_table[int(((h-p_old.offh)*p_old.s2+PI/2)*10000)+250000]);
                    
                }else if(p_new.type == 2) dz = ((h-p_new.offh)^(w-p_new.offw)) - ((h-p_old.offh)^(w-p_old.offw));
                else {
                    dz = p_new.amp*exp_table[int(((h-p_new.offh)*(h-p_new.offh)*p_new.s1 + (w-p_new.offw)*(w-p_new.offw)*p_new.s2)*1000)] - p_old.amp*exp_table[int(((h-p_old.offh)*(h-p_old.offh)*p_old.s1 + (w-p_old.offw)*(w-p_old.offw)*p_old.s2)*1000)];
                }
                grid_orig[c] += dz;
                chmin(min_grid,grid_orig[c]);
                chmax(max_grid,grid_orig[c]);
                c++;
            }
        }
        int score1 = 0;
        {
           int h,w,d;
           rep(i,si(checked_h)){
               h = checked_h[i], w = checked_w[i];
               d = int((grid_orig[h*N+w]-min_grid)/(max_grid-min_grid)*255)-checked_grid[i];
               score1 += d*d;
           }
        }
        if(score1 <= score){
            score = score1;
        }else{
            if((double)(xorshift())/(double)UINT_MAX < exp((score-score1)/T)){
                score = score1;
            }else{
                rep(i,N*N) grid_orig[i] = grid_orig0[i];
                params[n] = param0;
            }
        }
        
        time = timer.get_time();
        //T = T1 + (T0-T1)*time/time_limit;
        //T = pow(T1,(1.0-time/time_limit)) + pow(T0,(time/time_limit));
        T = 700.0/(1.0 + exp(-(0.70-time/time_limit)*31.0));
        
        if(chmin(best_score,score)) best_params = params;
        //if(TESTER == 0 && LOCAL == 1 && loop % 1000 == 0) cerr << loop csp score csp score1 csp time csp T << endl;
        if(TESTER == 0 && LOCAL == 1 && FILE_OUT == 1 && loop % 10 == 0) fprintf(fp,"%d %.6lf %d %d %.5lf\n",loop,time,score,score1,T);
        
    }
    if(LOCAL == 1 && TESTER == 0) cerr << " (" << loop << "|" << score << ") ";
    if(TESTER == 0 && LOCAL == 1) fclose(fp);
    updata_grid(best_params,grid);
    rep(i,N) rep(j,N) surface[i*N+j] = grid[i*N+j];
    
}



int solve(int seed, int N_IN=-1, double A_IN=-1){
    TimerChrono timer;
    timer.start();
    mt19937 engine(seed);
    init_exp_table();
    init_sin_table();
    
    int grid_exact[NN_MAX] = {};
    if(LOCAL == 1) {
        N = N_IN;
        A = A_IN;
        init(engine,N,A,E,grid_exact);
    }else if(LOCAL == 0) {
        cin >> A >> N;
    }
    
    int surface[NN_MAX] = {};
    int checked[NN_MAX] = {};
    memset(surface,-1,sizeof(surface));
    
    int num_d = 9;
    
    if(N <= 25){
        if(A >= 0.7) num_d = 4;
        else if(A >= 0.4) num_d = 2;
        else num_d = 2;
    }else if(N <= 40) {
        if(A >= 0.7) num_d = 6;
        else if(A >= 0.4) num_d = 4;
        else num_d = 2;
    }else if(N <= 50){
        if(A >= 0.7) num_d = 6;
        else if(A >= 0.4) num_d = 6;
        else num_d = 4;
    }else {
        if(A >= 0.7) num_d = 8;
        else if(A >= 0.4) num_d = 6;
        else num_d = 4;
    }
    num_d++;
    
    int dist = N/num_d;
    int h_check = dist/2+(N-dist*num_d)/2, w_check = dist/2+(N-dist*num_d)/2;
    
    int query = 0;
    vi checked_h,checked_w,checked_grid;
    int ccc = 0;
    rep(i,num_d) {
        if(h_check+1 >= N) break;
        rep(j,num_d){
            if(w_check+1 >= N) break;
            ccc++;
            if(ccc % 2 == 1){
                query++;
                if(LOCAL == 0) cout << h_check csp w_check << endl;
                rep(k,9) {
                    int nh = h_check + dh_9[k], nw = w_check + dw_9[k];
                    checked[nh*N+nw]++;
                    int res;
                    if(LOCAL == 0) cin >> res;
                    else if(LOCAL == 1) res = grid_exact[nh*N+nw];
                    surface[nh*N+nw] = res;
                    
                    checked_h.pb(nh);
                    checked_w.pb(nw);
                    checked_grid.pb(surface[nh*N+nw]);
                }
                if(LOCAL == 0){
                    int timevalue;
                    cin >> timevalue;
                }
            }
            w_check += dist;
        }
        h_check += dist;
        w_check = dist/2+(N-dist*num_d)/2;
        
        chmax(w_check,1);
    }
    if(LOCAL == 0) cout << "done" << endl;
    
    
    int sim_num = 4;
    if(N < 40) sim_num = 8;
    int surface_sum[NN_MAX] = {};
    rep(sim,sim_num){
        int surface_i[NN_MAX];
        memcpy(surface_i,surface,sizeof(surface));
        SA(checked_h,checked_w,checked_grid,surface_i,(TIME_LIMIT-timer.get_time())/(sim_num-sim));
        rep(i,N) rep(j,N) surface_sum[i*N+j] += surface_i[i*N+j];
        if(LOCAL == 1 && TESTER == 0) cerr << cal_score(surface_i, grid_exact, query) << " ";
    }
    rep(i,N*N) surface[i] = surface_sum[i]/sim_num;
    
    if(LOCAL == 1){
        if(TESTER == 0) cerr << endl;
        int score = cal_score(surface, grid_exact, query);
        if(TESTER == 0) cerr << score csp score-int(100000.0*(1.0-A)*(double)query/(N*N/9)) csp int(100000.0*(1.0-A)*(double)query/(N*N/9)) csp N csp A csp E csp timer.get_time() csp seed << endl;
        
        if(TESTER == 0){
            int diff_gird[NN_MAX] = {};
            rep(i,N) rep(j,N) diff_gird[i*N+j] = surface[i*N+j]-grid_exact[i*N+j];
            print_svg_grid(surface, 1, checked, "fig_gird_pred.svg");
            print_svg_grid(diff_gird, 10, checked, "fig_gird_diff.svg");
            print_svg_grid(grid_exact, 0, checked, "fig_gird_exact.svg");
        }
        
        return score;
        
    }else if(LOCAL == 0){
        rep(i,N) rep(j,N) cout << surface[i*N+j] << endl;
    }
    
    
    return 1;
}



int main(){
        
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    srand((unsigned)time(NULL));
    
    int seed = rand();
    if(TESTER == 0){
        solve(1237640781);
    }else if(TESTER == 1){
        ll score_sum = 0;
        rep(t,10){
            int N_IN = xorshift()%10+20;
            double A_IN = (double)(xorshift()%infi)/(double)infi*0.1 + 0.2;
            int score = solve(t,N_IN,A_IN);
            score_sum += score;
        }
        cout << score_sum << endl;
    }
    
    return 0;
}








Param rand_param(mt19937 &engine){
    int type = engine()%11;
    chmin(type,3);
    int offh = engine()%N, offw = engine()%N;
    double s1=0,s2=0,amp=0;
    if(type == 0) {
        s1 = (double)(engine()%INT_MAX)/(double)(INT_MAX)*10.0*(engine()%2==0?-1:1);
        s2 = (double)(engine()%INT_MAX)/(double)(INT_MAX)*10.0*(engine()%2==0?-1:1);
    }else if(type == 1){
        s1 = (double)(engine()%INT_MAX)/(double)(INT_MAX)*0.4*(engine()%2==0?-1:1);
        s2 = (double)(engine()%INT_MAX)/(double)(INT_MAX)*0.4*(engine()%2==0?-1:1);
        amp = (double)(engine()%INT_MAX)/(double)(INT_MAX)*90.0 + 10.0;
    }else if(type == 2){
    }else{
        s1 = (double)(engine()%INT_MAX)/(double)(INT_MAX)*0.099 + 0.001;
        s2 = (double)(engine()%INT_MAX)/(double)(INT_MAX)*0.099 + 0.001;
        amp = (double)(engine()%INT_MAX)/(double)(INT_MAX)*90.0 + 10.0;
    }
    Param P(type,offh,offw,s1,s2,amp);
    return P;
}

void init(mt19937 &engine, int &N, double &A, int &E, int (&grid)[NN_MAX]){
    if(N < 0) N = engine()%41+20;
    if(A < 0) A = (double)(engine()%INT_MAX)/(double)(INT_MAX)*0.6 + 0.2;
    E = engine()%29 + 2;
    ve<Param> params_exact;
    rep(i,E){
        Param P = rand_param(engine);
        params_exact.pb(P);
    }
    double grid0[NN_MAX] = {};
    rep(i,E) apply(grid0,params_exact[i]);
    normalize(grid0,grid);
}

void print_svg_grid(int (&grid)[NN_MAX], int f, int (&checked)[NN_MAX], string name="-1"){
    if(name == "-1") name = "fig_grid.svg";
    SVGWriter svg(name);
    svg.newSVG();
    rep(i,N) rep(j,N){
        string c = "rgb(" + to_string(grid[i*N+j])+ "," + to_string(grid[i*N+j]) + "," + to_string(grid[i*N+j]) + ")";
        if(f == 10) {
            if(grid[i*N+j] >= 0) c = "rgb(255," + to_string(255-grid[i*N+j]) + ", " + to_string(255-grid[i*N+j])+ ")";
            else c = "rgb(" + to_string(255+grid[i*N+j]) + ", " + to_string(255+grid[i*N+j]) + ",255)";
        }else{
            if(grid[i*N+j] <= 255/2) c = "rgb(255," + to_string(grid[i*N+j]*2) + ", " + to_string(grid[i*N+j]*2)+ ")";
            else c = "rgb(" + to_string(255-2*(grid[i*N+j]-255/2)) + ", " + to_string(255-2*(grid[i*N+j]-255/2)) + ",255)";
        }
        svg.rect(600.0*(double)i/N,600.0*(double)j/N,600.0/(double)N,600.0/(double)N,c);
        if(checked[i*N+j] != 0) svg.circle(600.0*(double)i/N+600.0/(double)N/2,600.0*(double)j/N+600.0/(double)N/2,600.0/(double)N/7,"red");
    }
    svg.save();
}

