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

const int numberOfCellsMAX = 105;
const int numberOfBasesMAX = 2;

int ti;
double time_,timeRap;
std::chrono::system_clock::time_point start, end_;
std::time_t time_stamp;



int basesIdx[2][numberOfBasesMAX];
int numberOfCells,numberOfBases;
int initTotalNumberOfEggs = 0, initTotalNumberOfCrystals = 0;
int DistanceCell[numberOfCellsMAX][numberOfCellsMAX];

struct Cell {
    int idx, type, resources;
    int ants[2];
    int neighbors[6];

    Cell(){}
    Cell(int idx_, int type_, int resources_, int neighbors_[6], int myAnts_ = 0, int oppAnts_ = 0) {
        idx = idx_;
        type = type_, resources = resources_;
        rep(i,6) neighbors[i] = neighbors_[i];
        ants[0] = myAnts_, ants[1] = oppAnts_;
    }
};

struct State {
    Cell cells[numberOfCellsMAX];
    int beacon[2][numberOfCellsMAX];
    int totalNumberOfEggs = 0, totalNumberOfCrystals = 0;
    int turn = 0;
    int acquiredAnts[2] = {0,0};
    int acquiredCrystals[2] = {0,0};
    
    State(Cell (&cells_)[numberOfCellsMAX]):cells(cells_){
        rep(p,2) rep(i,numberOfCellsMAX) beacon[p][i] = 0;
        rep(i,numberOfCells){
            Cell &c = cells[i];
            if(c.type == 1) totalNumberOfEggs += c.resources;
            else if(c.type == 2) totalNumberOfCrystals += c.resources;
        }
        initTotalNumberOfEggs = totalNumberOfEggs, initTotalNumberOfCrystals = totalNumberOfCrystals;
    }
    void update_state_from_input(int is_input){
        totalNumberOfEggs = totalNumberOfCrystals = 0;
        acquiredAnts[0] = acquiredAnts[1] = 0;
        if(is_input == 1){
            cin >> acquiredCrystals[0] >> acquiredCrystals[1];
            for (int i = 0; i < numberOfCells; i++) {
                int resources, myAnts, oppAnts;
                cin >> resources >> myAnts >> oppAnts;
                cells[i].resources = resources;
                cells[i].ants[0] = myAnts;
                cells[i].ants[1] = oppAnts;
                acquiredAnts[0] += myAnts;
                acquiredAnts[1] += oppAnts;
                if(cells[i].type == 1) totalNumberOfEggs += resources;
                else if(cells[i].type == 2) totalNumberOfCrystals += resources;
            }
        }
        rep(i,numberOfCells) if(cells[i].resources == 0) cells[i].type = 0;
    }
    
    void printInfo(){
        cerr << "myAnts:" << acquiredAnts[0] << " oppAnts:" << acquiredAnts[1] << endl;
        cerr << "myCrystals:" << acquiredCrystals[0] << " oppCrystals:" << acquiredCrystals[1] << endl;
    }
    void printAction(){
        int f = 0;
        rep(i,numberOfCells)if(beacon[0][i] > 0){
            cout << "BEACON" csp i csp beacon[0][i] << ";";
            f++;
        }
        if(f == 0) cout << "WAIT;";
        cout << endl;
    }
    
    
    void advance(vi newBeacon, int player){
        rep(i,numberOfCells) {
            beacon[player][i] = newBeacon[i];
        }
    }
    
    void antAllocationAndMove(int player=0){
        //allocation
        int antSum = 0, beaconSum = 0;
        rep(i,numberOfCells){
            antSum += cells[i].ants[player];
            beaconSum += beacon[player][i];
        }
        
        double scalingFactor = (double)antSum / (double)beaconSum;
        int beaconCell[numberOfCells], wiggleRoom[numberOfCells];
        rep(i,numberOfCells) {
            beaconCell[i] = 0;
            int highBeaconValue = (int)ceil((double)beacon[player][i] * scalingFactor);
            int lowBeaconValue = (int)(beacon[player][i] * scalingFactor);
            beaconCell[i] = max(1, lowBeaconValue);
            wiggleRoom[i] = highBeaconValue - beaconCell[i];
        }
        ve<ll> allPairs;
        rep(i,numberOfCells) if(cells[i].ants[player] > 0){
            rep(j,numberOfCells) if(beacon[player][j] != 0) allPairs.pb(DistanceCell[i][j]*1000000ll+i*1000ll+j);
        }
        sort(all(allPairs));
        
        ve<ll> allocations;
        bool stragglers = false;
        int cn_fin = 0;
        int used[si(allPairs)];
        rep(i,si(allPairs)) used[i] = 0;
        while (cn_fin < si(allPairs)) {
            int ii = -1;
            for (ll dist_ant_beacon : allPairs) {
                ii++;
                int antIdx = (dist_ant_beacon%1000000ll)/1000ll;
                int beaconIdx = dist_ant_beacon%1000ll;
                if(used[ii] != 0) continue;
                
                int maxAlloc = (int)(stragglers ? min(cells[antIdx].ants[player], beaconCell[beaconIdx] + wiggleRoom[beaconIdx]) : min(cells[antIdx].ants[player], beaconCell[beaconIdx]));
                if (maxAlloc > 0) {
                    allocations.pb(antIdx*1000000ll + beaconIdx*1000ll + maxAlloc);
                    cells[antIdx].ants[player] -= maxAlloc;
                    if (!stragglers) {
                        beaconCell[beaconIdx] -= maxAlloc;
                    } else {
                        beaconCell[beaconIdx] -= (maxAlloc - wiggleRoom[beaconIdx]);
                        wiggleRoom[beaconIdx] = 0;
                    }
                }
                
                if(cells[antIdx].ants[player] <= 0) used[ii] = 1, cn_fin++;
            }
            stragglers = true;
        }
        
        //move
        for (ll ant_beacon_maxAlloc : allocations) {
            int antIdx = (ant_beacon_maxAlloc/1000000ll);
            int beaconIdx = (ant_beacon_maxAlloc%1000000ll)/1000ll;
            int maxAlloc = ant_beacon_maxAlloc%1000ll;
            if(antIdx == beaconIdx) {
                cells[antIdx].ants[player] += maxAlloc;
                //cerr << "antmove : " << antIdx csp antIdx csp maxAlloc << endl;
                continue;
            }
            int idx = 0;
            ll dd = -infl;
            int ddd = infi;
            rep(i,6) if(cells[antIdx].neighbors[i] != -1){
                int v = cells[antIdx].neighbors[i];
                chmin(ddd,DistanceCell[beaconIdx][v]);
            }
            rep(i,6) if(cells[antIdx].neighbors[i] != -1){
                int v = cells[antIdx].neighbors[i];
                if(DistanceCell[beaconIdx][v] != ddd) continue;
                if(chmax(dd,(ll)cells[v].ants[player]*10000000+beacon[player][beaconIdx])) idx = v;
            }
            //cells[antIdx].ants[player] -= maxAlloc;
            cells[idx].ants[player] += maxAlloc;
            //cerr << "antmove : " << antIdx csp idx csp maxAlloc << endl;
        }
        
        //cerr << player << "::ant:";rep(i,numberOfCells) {cerr << cells[i].ants[player] << " ";}cerr << endl;
    }
    
    
    pii harvest(int harvestType){
        int chanes[2][numberOfCells];
        int harvestAnt[2] = {0,0}, harvestCrystal[2] = {0,0};
        rep(p,2){
            rep(i,numberOfCells) chanes[p][i] = 0;
            pq2<int> pq;
            rep(i,numberOfBases) {
                chanes[p][basesIdx[p][i]] = cells[basesIdx[p][i]].ants[p];
                pq.push(chanes[p][basesIdx[p][i]]*1000+basesIdx[p][i]);
            }
            while (!pq.empty()) {
                int c_u = pq.top();
                int c = c_u/1000, u = c_u%1000;
                pq.pop();
                if(chanes[p][u] > c) continue;
                rep(neigh,6){
                    int v = cells[u].neighbors[neigh];
                    if(v == -1) continue;
                    if(chmax(chanes[p][v], min(cells[v].ants[p],c))) pq.push(chanes[p][v]*1000+v);
                }
            }
        }
        
        int harvestResource[numberOfCells];
        rep(i,numberOfCells) harvestResource[i] = 0;
        rep(p,2){
            int canMove[numberOfCells];
            rep(i,numberOfCells) {
                if(chanes[p][i] >= chanes[(p+1)%2][i]) canMove[i] = 1;
                else canMove[i] = 0;
            }
            queue<int> que;
            rep(i,numberOfBases) if(chanes[p][basesIdx[p][i]] > 0) {
                canMove[basesIdx[p][i]] = 2;
                que.push(basesIdx[p][i]);
            }
            while (!que.empty()) {
                int u = que.front();
                que.pop();
                rep(neigh,6){
                    int v = cells[u].neighbors[neigh];
                    if(v == -1) continue;
                    if(canMove[v] == 1 && chanes[p][v] > 0){
                        canMove[v] = 2;
                        que.push(v);
                    }
                }
            }
            rep(i,numberOfCells) if(canMove[i] == 2 && cells[i].type == harvestType){
                if(cells[i].type == 1){
                    int r = min(cells[i].resources, chanes[p][i]);
                    harvestAnt[p] += r;
                    harvestResource[i] += r;
                    acquiredAnts[p] += r*numberOfBases;
                    rep(base,numberOfBases) cells[basesIdx[p][base]].ants[p] += r;
                }else if(cells[i].type == 2){
                    int r = min(cells[i].resources, chanes[p][i]);
                    harvestCrystal[p] += r;
                    harvestResource[i] += r;
                    acquiredCrystals[p] += r;
                }
            }
        }
        
        rep(i,numberOfCells) cells[i].resources = max(0,cells[i].resources-harvestResource[i]);
        
        if(harvestType == 1) return make_pair(harvestAnt[0],harvestAnt[1]);
        else return make_pair(harvestCrystal[0],harvestCrystal[1]);
    }
};







//---------------------------------------------------------
//action
vi simpleAction(State &state, int player=0){
    
    vi bestBeacon(numberOfCells,0);
    int bestHarvestScore = -infi;
    
    int policyChangeEggCrystal = 0;
    if(state.totalNumberOfCrystals > initTotalNumberOfCrystals*3/4) policyChangeEggCrystal = 1;
    if(state.acquiredCrystals[0]+state.acquiredCrystals[1] >= int((double)initTotalNumberOfCrystals*0.7)
       || state.acquiredCrystals[1] > int((double)initTotalNumberOfCrystals*0.5*0.5)) policyChangeEggCrystal = 2;
    
    
    int simulationLength = 4;
    if(state.acquiredCrystals[0]+state.acquiredCrystals[1] >= int((double)initTotalNumberOfCrystals*0.8)) simulationLength = 2;
    if(policyChangeEggCrystal == 1) simulationLength = 6;
    
    
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    vi newBeacon(numberOfCells,0);
    int usedNumberOfBeacon = 0;
    vi d(numberOfCells,infi);
    vi from(numberOfCells,-1);
    vi visitOrder;
    pq1<int> pq;
    rep(j,numberOfBases){
        pq.push(0*1000+basesIdx[player][j]);
        d[basesIdx[player][j]] = 0;
    }
    
    rep(i,numberOfCells) if(state.cells[i].type > 0){
        while (!pq.empty()) {
            ll c_u = pq.top();
            int c = c_u/1000, u = c_u%1000;
            pq.pop();
            if(d[u] < c) continue;
            Cell &cNow = state.cells[u];
            rep(neigh, 6) if(cNow.neighbors[neigh] != -1){
                int v = cNow.neighbors[neigh];
                int w = 1;
                if(state.cells[v].type != 0) w = 0;
                if(state.cells[v].ants[player] <= 0) w++;
                
                
                if(chmin(d[v], d[u]+w)) {
                    pq.push(d[v]*1000+v);
                    from[v] = u;
                }
            }
        }
        
        int idx = -1, min_d = infi;
        rep(j,numberOfCells) if(state.cells[j].type > 0){
            if(policyChangeEggCrystal == 1){
                if(state.cells[j].type != 1) continue;
            }else if(policyChangeEggCrystal == 2){
                if(state.cells[j].type != 2) continue;
            }
            if(newBeacon[j] != 0) continue;
            if(chmin(min_d,d[j])) idx = j;
        }
        if(idx == -1) {
            rep(j,numberOfCells) if(state.cells[j].type > 0){
                if(newBeacon[j] != 0) continue;
                if(chmin(min_d,d[j])) idx = j;
            }
        }
        
        int s = idx;
        int cn = 0;
        while (s != -1) {
            s = from[s];
            if(newBeacon[s] != 0) cn += 1;
        }
    
        if(usedNumberOfBeacon+cn > state.acquiredAnts[player]) continue;
        usedNumberOfBeacon += cn;
        s = idx;
        while (s != -1) {
            if(newBeacon[s] == 0) {
                d[s] = 0;
                pq.push(s);
            }
            newBeacon[s] = 1;
            s = from[s];
        }
        
        
        State newState = state;
        int futureHarvest = 0;
        double r = 1.0;
        rep(t,simulationLength){
            newState.advance(newBeacon,0);
            newState.antAllocationAndMove(0);
            pii harvestAnt = newState.harvest(1);
            pii harvestCrystal = newState.harvest(2);
            //futureHarvest += (harvestAnt.fi + harvestCrystal.fi) - (harvestAnt.se + harvestCrystal.se);
            futureHarvest += int(100*(harvestAnt.fi*(policyChangeEggCrystal==1?7:1) + (policyChangeEggCrystal==2?3:1)*harvestCrystal.fi)*r);
            if(harvestCrystal.se*2 >= initTotalNumberOfCrystals) futureHarvest -= 1000;
            //r *= 0.9;
            time_ = static_cast<double>(chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now() - start).count() / 1000.0);
            if(time_ > 90) break;
        }
        //cerr << bestHarvestScore csp futureHarvest << endl;
        if(chmax(bestHarvestScore, futureHarvest)) {
            bestBeacon = newBeacon;
            visitOrder.pb(idx);
        }
        time_ = static_cast<double>(chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now() - start).count() / 1000.0);
        if(time_ > 93) break;
       
    }
    if(si(bestBeacon) <= 0) bestBeacon = newBeacon;
    
    
    {
        newBeacon.assign(numberOfCells,0);
        usedNumberOfBeacon = 0;
        d.assign(numberOfCells,infi);
        from.assign(numberOfCells,-1);
        while (!pq.empty()) pq.pop();
        rep(j,numberOfBases){
            pq.push(0*1000+basesIdx[player][j]);
            d[basesIdx[player][j]] = 0;
        }
        State state1 = state;
        state1.advance(newBeacon,0);
        state1.antAllocationAndMove(0);
        rep(i,numberOfCells) if(state1.cells[i].type > 0){
            while (!pq.empty()) {
                ll c_u = pq.top();
                int c = c_u/1000, u = c_u%1000;
                pq.pop();
                if(d[u] < c) continue;
                Cell &cNow = state1.cells[u];
                rep(neigh, 6) if(cNow.neighbors[neigh] != -1){
                    int v = cNow.neighbors[neigh];
                    int w = 1;
                    if(state1.cells[v].type != 0) w = 0;
                    if(state1.cells[v].ants[player] <= 0) w++;
                    
                    
                    if(chmin(d[v], d[u]+w)) {
                        pq.push(d[v]*1000+v);
                        from[v] = u;
                    }
                }
            }
            
            int idx = -1, min_d = infi;
            rep(j,numberOfCells) if(state1.cells[j].type > 0){
                if(policyChangeEggCrystal == 1){
                    if(state1.cells[j].type != 1) continue;
                }else{
                    if(state1.cells[j].type != 2) continue;
                }
                if(newBeacon[j] != 0) continue;
                if(chmin(min_d,d[j])) idx = j;
            }
            if(idx == -1) {
                rep(j,numberOfCells) if(state1.cells[j].type > 0){
                    if(newBeacon[j] != 0) continue;
                    if(chmin(min_d,d[j])) idx = j;
                }
            }
            
            int s = idx;
            int cn = 0;
            while (s != -1) {
                s = from[s];
                if(newBeacon[s] != 0) cn += 1;
            }
        
            if(usedNumberOfBeacon+cn > state1.acquiredAnts[player]) continue;
            usedNumberOfBeacon += cn;
            s = idx;
            while (s != -1) {
                if(newBeacon[s] == 0) {
                    d[s] = 0;
                    pq.push(s);
                }
                newBeacon[s] = 1;
                s = from[s];
            }
            
            
            State newState = state;
            int futureHarvest = 0;
            double r = 1.0;
            rep(t,simulationLength){
                newState.advance(newBeacon,0);
                newState.antAllocationAndMove(0);
                pii harvestAnt = newState.harvest(1);
                pii harvestCrystal = newState.harvest(2);
                //futureHarvest += (harvestAnt.fi + harvestCrystal.fi) - (harvestAnt.se + harvestCrystal.se);
                futureHarvest += int(100*(harvestAnt.fi*(policyChangeEggCrystal==1?7:1) + (policyChangeEggCrystal==2?3:1)*harvestCrystal.fi)*r);
                if(harvestCrystal.se*2 >= initTotalNumberOfCrystals) futureHarvest -= 1000;
                //r *= 0.9;
                time_ = static_cast<double>(chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now() - start).count() / 1000.0);
                if(time_ > 90) break;
            }
            //cerr << bestHarvestScore csp futureHarvest << endl;
            if(chmax(bestHarvestScore, futureHarvest)) {
                bestBeacon = newBeacon;
                visitOrder.pb(idx);
            }
            time_ = static_cast<double>(chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now() - start).count() / 1000.0);
            if(time_ > 85) break;
        }
    }
    
    time_ = static_cast<double>(chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now() - start).count() / 1000.0);
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    while (time_ < 80) {
        newBeacon.assign(numberOfCells,0);
        usedNumberOfBeacon = 0;
        d.assign(numberOfCells,infi);
        from.assign(numberOfCells,-1);
        while (!pq.empty()) pq.pop();
        rep(j,numberOfBases){
            pq.push(0*1000+basesIdx[player][j]);
            d[basesIdx[player][j]] = 0;
        }

        int n1 = rand()%si(visitOrder);
        int n2 = rand()%si(visitOrder);
        swap(visitOrder[n1],visitOrder[n2]);
    
        time_ = static_cast<double>(chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now() - start).count() / 1000.0);
        
        for(int i : visitOrder){
            while (!pq.empty()) {
                ll c_u = pq.top();
                int c = c_u/1000;
                int u = c_u%1000;
                pq.pop();
                if(d[u] < c) continue;
                Cell &cNow = state.cells[u];
                rep(neigh, 6) if(cNow.neighbors[neigh] != -1){
                    int v = cNow.neighbors[neigh];
                    int w = 1;
                    if(state.cells[v].type != 0) w = 0;
                    if(state.cells[v].ants[player] <= 0) w++;
                    
                    
                    if(chmin(d[v], d[u]+w)) {
                        pq.push(d[v]*1000+v);
                        from[v] = u;
                    }
                }
            }
            
            time_ = static_cast<double>(chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now() - start).count() / 1000.0);
            if(time_ > 85) break;
            
            int idx = i;
            int s = idx;
            int cn = 0;
            while (s != -1) {
                s = from[s];
                if(newBeacon[s] != 0) cn += 1;
            }
        
            usedNumberOfBeacon += cn;
            s = idx;
            while (s != -1) {
                if(newBeacon[s] == 0) {
                    d[s] = 0;
                    pq.push(s);
                }
                newBeacon[s] = 1;
                s = from[s];
            }
            
            State newState = state;
            int futureHarvest = 0;
            double r = 1.0;
            rep(t,simulationLength){
                newState.advance(newBeacon,0);
                newState.antAllocationAndMove(0);
                pii harvestAnt = newState.harvest(1);
                pii harvestCrystal = newState.harvest(2);
                //futureHarvest += (harvestAnt.fi + harvestCrystal.fi) - (harvestAnt.se + harvestCrystal.se);
                futureHarvest += int(100*(harvestAnt.fi*(policyChangeEggCrystal==1?7:1) + (policyChangeEggCrystal==2?3:1)*harvestCrystal.fi)*r);
                if(harvestCrystal.se*2 >= initTotalNumberOfCrystals) futureHarvest -= 1000;
                //r *= 0.9;
                time_ = static_cast<double>(chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now() - start).count() / 1000.0);
                if(time_ > 85) break;
            }
            //cerr << bestHarvestScore csp futureHarvest << endl;
            if(chmax(bestHarvestScore, futureHarvest)) {
                bestBeacon = newBeacon;
            }else{
                swap(visitOrder[n1],visitOrder[n2]);
            }
            time_ = static_cast<double>(chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now() - start).count() / 1000.0);
            if(time_ > 85) break;
        }
    }
    
    
    
    
    
    return bestBeacon;
    
}


vi greedyAction(State &state, int player=0){
    
    rep(i,numberOfCells) if(state.cells[i].ants[(player+1)%2] > 0) state.cells[i].ants[(player+1)%2]++;
    
    vi bestBeacon(numberOfCells,0);
    int bestHarvestScore = -infi;
    
    int policyChangeEggCrystal = 0;
    if(state.totalNumberOfCrystals > initTotalNumberOfCrystals*3/4) policyChangeEggCrystal = 1;
    if(state.acquiredCrystals[0]+state.acquiredCrystals[1] >= int((double)initTotalNumberOfCrystals*0.7)
       || state.acquiredCrystals[1] > int((double)initTotalNumberOfCrystals*0.5*0.5)) policyChangeEggCrystal = 2;
    
    int simulationLength = 4;
    if(state.acquiredCrystals[0]+state.acquiredCrystals[1] >= int((double)initTotalNumberOfCrystals*0.8)) simulationLength = 2;
    if(policyChangeEggCrystal == 1) simulationLength = 6;
    
    
    
    vi newBeacon(numberOfCells,0);
    int usedNumberOfBeacon = 0;
    vi d(numberOfCells,infi);
    vi from(numberOfCells,-1);
    vi visitOrder;
    pq1<int> pq;
    rep(j,numberOfBases){
        pq.push(0*1000+basesIdx[player][j]);
        d[basesIdx[player][j]] = 0;
    }
    rep(i,numberOfCells) if(state.cells[i].type > 0){
        while (!pq.empty()) {
            ll c_u = pq.top();
            int c = c_u/1000, u = c_u%1000;
            pq.pop();
            if(d[u] < c) continue;
            Cell &cNow = state.cells[u];
            rep(neigh, 6) if(cNow.neighbors[neigh] != -1){
                int v = cNow.neighbors[neigh];
                int w = 1;
                if(state.cells[v].type != 0) w = 0;
                if(state.cells[v].ants[player] <= 0) w++;
                
                
                if(chmin(d[v], d[u]+w)) {
                    pq.push(d[v]*1000+v);
                    from[v] = u;
                }
            }
        }
        
        int idx = -1, min_d = infi;
        rep(j,numberOfCells) if(state.cells[j].type > 0){
            if(policyChangeEggCrystal == 1){
                if(state.cells[j].type != 1) continue;
            }else if(policyChangeEggCrystal == 2){
                if(state.cells[j].type != 2) continue;
            }
            if(newBeacon[j] != 0) continue;
            if(chmin(min_d,d[j])) idx = j;
        }
        if(idx == -1) {
            rep(j,numberOfCells) if(state.cells[j].type > 0){
                if(newBeacon[j] != 0) continue;
                if(chmin(min_d,d[j])) idx = j;
            }
        }
        
        int s = idx;
        int cn = 0;
        while (s != -1) {
            s = from[s];
            if(newBeacon[s] != 0) cn += 1;
        }
    
        if(usedNumberOfBeacon+cn > state.acquiredAnts[player]) continue;
        usedNumberOfBeacon += cn;
        s = idx;
        while (s != -1) {
            if(newBeacon[s] == 0) {
                d[s] = 0;
                pq.push(s);
            }
            newBeacon[s] = 1;
            s = from[s];
        }
        
        
        State newState = state;
        int futureHarvest = 0;
        double r = 1.0;
        rep(t,simulationLength){
            newState.advance(newBeacon,0);
            newState.antAllocationAndMove(0);
            pii harvestAnt = newState.harvest(1);
            pii harvestCrystal = newState.harvest(2);
            //futureHarvest += (harvestAnt.fi + harvestCrystal.fi) - (harvestAnt.se + harvestCrystal.se);
            futureHarvest += int(100*(harvestAnt.fi*(policyChangeEggCrystal==1?7:1) + (policyChangeEggCrystal==2?3:1)*harvestCrystal.fi)*r);
            //r *= 0.9;
            time_ = static_cast<double>(chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now() - start).count() / 1000.0);
            if(time_ > 90) break;
        }
        //cerr << bestHarvestScore csp futureHarvest << endl;
        if(chmax(bestHarvestScore, futureHarvest)) {
            bestBeacon = newBeacon;
            visitOrder.pb(idx);
        }
        time_ = static_cast<double>(chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now() - start).count() / 1000.0);
        if(time_ > 93) break;
       
    }
    if(si(bestBeacon) <= 0) bestBeacon = newBeacon;
    newBeacon = bestBeacon;
    
    int bestHarvestScore_back = 0;
    vi bestBeacon_back = bestBeacon;
    {
        State newState = state;
        int futureHarvest = 0;
        rep(t,2){
            newState.advance(newBeacon,0);
            newState.antAllocationAndMove(0);
            pii harvestAnt = newState.harvest(1);
            pii harvestCrystal = newState.harvest(2);
            //futureHarvest += (harvestAnt.fi + harvestCrystal.fi) - (harvestAnt.se + harvestCrystal.se);
            futureHarvest += int(100*(harvestAnt.fi*(policyChangeEggCrystal==1?7:1) + (policyChangeEggCrystal==2?3:1)*harvestCrystal.fi));
        }
        bestHarvestScore_back = futureHarvest;
    }
    
    
    
    bestBeacon.assign(numberOfCells,0);
    vi nowBeacon(numberOfCells,0);
    rep(i,numberOfCells){
        bestBeacon[i] = newBeacon[i], nowBeacon[i] = 0;
        if(newBeacon[i] > 0 && state.beacon[player][i] > 0) nowBeacon[i] = 1;
    }
    bestHarvestScore = 0;
    queue<int> que;
    rep(i,numberOfBases) que.push(basesIdx[player][i]);
    
    int futureHarvest = 0;
    State newState = state;
    rep(tt,2){
        newState.advance(nowBeacon,0);
        newState.antAllocationAndMove(0);
        pii harvestAnt = newState.harvest(1);
        pii harvestCrystal = newState.harvest(2);
        futureHarvest += harvestAnt.fi*7+harvestCrystal.fi;
    }
    if(chmax(bestHarvestScore, futureHarvest)) bestBeacon = nowBeacon;
    
    vi beacon_cells;
    rep(i,numberOfCells) if(nowBeacon[i] > 0) beacon_cells.pb(i);
    while (time_ < 90) {
        
        int max_harvest = -infi;
        int id = 0;
        for (int i : beacon_cells) {
            nowBeacon[i]++;
            
            int futureHarvest = 0;
            State newState = state;
            rep(tt,2){
                newState.advance(nowBeacon,0);
                newState.antAllocationAndMove(0);
                pii harvestAnt = newState.harvest(1);
                pii harvestCrystal = newState.harvest(2);
                futureHarvest += int(100*(harvestAnt.fi*(policyChangeEggCrystal==1?7:1) + (policyChangeEggCrystal==2?3:1)*harvestCrystal.fi));
            }
            if(chmax(max_harvest, futureHarvest)) id = i;
            
            nowBeacon[i]--;
            time_ = static_cast<double>(chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now() - start).count() / 1000.0);
            if(time_ > 90) break;
        }
        
        nowBeacon[id]++;
        int futureHarvest = 0;
        State newState = state;
        rep(tt,2){
            newState.advance(nowBeacon,0);
            newState.antAllocationAndMove(0);
            pii harvestAnt = newState.harvest(1);
            pii harvestCrystal = newState.harvest(2);
            futureHarvest += int(100*(harvestAnt.fi*(policyChangeEggCrystal==1?7:1) + (policyChangeEggCrystal==2?3:1)*harvestCrystal.fi));
        }
        if(chmax(bestHarvestScore, max_harvest)) bestBeacon = nowBeacon;
        
        time_ = static_cast<double>(chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now() - start).count() / 1000.0);
        if(time_ > 90) break;
    }

    if(bestHarvestScore_back > bestHarvestScore) bestBeacon = bestBeacon_back;
    
    
    vi res(numberOfCells,0);
    rep(i,numberOfCells) res[i] = bestBeacon[i];
    
    return res;
}


//-------------------------------------------------------------------------------------
void playGame(State &state){
    int loop = 0;
    while (1) {
        loop++;
        timeRap = 1.0 * (clock() - ti) / CLOCKS_PER_SEC;
        state.update_state_from_input(1);
        
        start = std::chrono::system_clock::now();
        
        int cn = 0;
        rep(i,numberOfCells) if(state.cells[i].type == 2) cn++;
        
        vi action;
        if(loop <= 4 || cn == 1) action = greedyAction(state,0);
        else action = simpleAction(state,0);
        state.advance(action,0);
        
        //vi action_opp = simpleAction(state,1);
        //state.advance(action_opp,1);
        
        state.antAllocationAndMove(0);
        //state.antAllocationAndMove(1);
        
        //pii harvestAnt = state.harvest(1);
        //pii harvestCrystal = state.harvest(2);
        //cerr << harvestAnt.fi csp harvestAnt.se << " | " << harvestCrystal.fi csp harvestCrystal.se << endl;
        
        state.printAction();
        state.printInfo();
        time_ = 1.0 * (clock() - ti) / CLOCKS_PER_SEC;
        cerr << "time :" csp time_-timeRap  << endl;
        
        time_ = static_cast<double>(chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now() - start).count() / 1000.0);
        cerr << "time :" csp time_ << endl;
    }
}



void solve(){
    
    ti = clock();
    
    cin >> numberOfCells;
    Cell cells[numberOfCellsMAX];
    for (int i = 0; i < numberOfCells; i++) {
        int type, initialResources;
        int neigh[6];
        cin >> type >> initialResources;
        rep(j,6){
            int a;
            cin >> a;
            neigh[j] = a;
        }
        Cell cell(i, type, initialResources, neigh, 0, 0);
        cells[i] = cell;
    }
    
    cin >> numberOfBases;
    for (int i = 0; i < numberOfBases; i++) {
        int myBaseIndex;
        cin >> myBaseIndex;
        basesIdx[0][i] = myBaseIndex;
    }
    for (int i = 0; i < numberOfBases; i++) {
        int oppBaseIndex;
        cin >> oppBaseIndex;
        basesIdx[1][i] = oppBaseIndex;
    }
    rep(i,numberOfCells) {
        rep(j,numberOfCells) DistanceCell[i][j] = infi;
        DistanceCell[i][i] = 0;
        pq1<pii> pq;
        pq.push(make_pair(0,i));
        while (!pq.empty()) {
            auto [c,u] = pq.top();
            pq.pop();
            if(DistanceCell[i][u] < c) continue;
            rep(neigh,6){
                int v = cells[u].neighbors[neigh];
                if(v == -1) continue;
                if(chmin(DistanceCell[i][v], c+1)) pq.push(make_pair(c+1,v));
            }
        }
    }

    State state(cells);
    playGame(state);
    
    
}


int main(){
    
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    solve();

    return 0;
}
