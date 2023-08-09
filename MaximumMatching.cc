#include <bits/stdc++.h>

using namespace std;

#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define all(x) begin(x), end(x)
#define sz(x) (int)(x).size()
typedef long long ll;
typedef pair<int, int> pii;
typedef vector<int> vi;

mt19937  rng;
const ll mod = 1000000007; // faster if const


ll modpow(ll b, ll e) {
	ll ans = 1;
	for (; e; b = b * b % mod, e /= 2)
		if (e & 1) ans = ans * b % mod;
	return ans;
}


int matInv(vector<vector<ll>>& A) {
	int n = sz(A); vi col(n);
	vector<vector<ll>> tmp(n, vector<ll>(n));
	rep(i,0,n) tmp[i][i] = 1, col[i] = i;

	rep(i,0,n) {
		int r = i, c = i;
		rep(j,i,n) rep(k,i,n) if (A[j][k]) {
			r = j; c = k; goto found;
		}
		return i;
found:
		A[i].swap(A[r]); tmp[i].swap(tmp[r]);
		rep(j,0,n) swap(A[j][i], A[j][c]), swap(tmp[j][i], tmp[j][c]);
		swap(col[i], col[c]);
		ll v = modpow(A[i][i], mod - 2);
		rep(j,i+1,n) {
			ll f = A[j][i] * v % mod;
			A[j][i] = 0;
			rep(k,i+1,n) A[j][k] = (A[j][k] - f*A[i][k]) % mod;
			rep(k,0,n) tmp[j][k] = (tmp[j][k] - f*tmp[i][k]) % mod;
		}
		rep(j,i+1,n) A[i][j] = A[i][j] * v % mod;
		rep(j,0,n) tmp[i][j] = tmp[i][j] * v % mod;
		A[i][i] = 1;
	}

	for (int i = n-1; i > 0; --i) rep(j,0,i) {
		ll v = A[j][i];
		rep(k,0,n) tmp[j][k] = (tmp[j][k] - v*tmp[i][k]) % mod;
	}

	rep(i,0,n) rep(j,0,n)
		A[col[i]][col[j]] = tmp[i][j] % mod + (tmp[i][j] < 0 ? mod : 0);
	return n;
}

vector<pii> generalMatching(int N, vector<pii>& ed) {
	vector<vector<ll>> mat(N, vector<ll>(N)), A;
	for (pii pa : ed) {
		int a = pa.first, b = pa.second, r = rand() % mod;
		mat[a][b] = r, mat[b][a] = (mod - r) % mod;
	}

	int r = matInv(A = mat), M = 2*N - r, fi, fj;
	assert(r % 2 == 0);

	if (M != N) do {
		mat.resize(M, vector<ll>(M));
		rep(i,0,N) {
			mat[i].resize(M);
			rep(j,N,M) {
				int r = rand() % mod;
				mat[i][j] = r, mat[j][i] = (mod - r) % mod;
			}
		}
	} while (matInv(A = mat) != M);

	vi has(M, 1); vector<pii> ret;
	rep(it,0,M/2) {
		rep(i,0,M) if (has[i])
			rep(j,i+1,M) if (A[i][j] && mat[i][j]) {
				fi = i; fj = j; goto done;
		} assert(0); done:
		if (fj < N) ret.emplace_back(fi, fj);
		has[fi] = has[fj] = 0;
		rep(sw,0,2) {
			ll a = modpow(A[fi][fj], mod-2);
			rep(i,0,M) if (has[i] && A[i][fj]) {
				ll b = A[i][fj] * a % mod;
				rep(j,0,M) A[i][j] = (A[i][j] - A[fi][j] * b) % mod;
			}
			swap(fi,fj);
		}
	}
	return ret;
}

//p_v denotes probability of preserving a vertex, p_e denotes probability of preserving an edge
vector<pii> stochasticMatching(int N, vector<pii> ed, vector<long double> &p_v, vector<vector<long double>> &p_e) {
    vector<int> sampled_vertices(N);
    vector<pii> realization_ed;
    
    rep(it, 0, N)
        if (uniform_real_distribution<long double>(0, 1)(rng) < p_v[it])
            sampled_vertices[it] = 1;
    
    rep(it, 0, (int)ed.size())
        if (sampled_vertices[ed[it].first] && sampled_vertices[ed[it].second] && uniform_real_distribution<long double>(0, 1)(rng) < p_e[ed[it].first][ed[it].second])
            realization_ed.push_back(ed[it]);
    
    return generalMatching(N, realization_ed);
}

long double stochasticMatching_avg(int N, vector<pii> ed, vector<long double> &p_v, vector<vector<long double>> &p_e) {
    long double sum = 0;
    int T = 5;
    
    rep(it, 0, T)
        sum += (long double)stochasticMatching(N, ed, p_v, p_e).size();
        
    return sum / (long double)T;
}

vector<pii> stochastic_inc(int N, vector<pii> ed, vector<long double> &p_v, vector<vector<long double>> &p_e, int R) {
    cout << "Est of expected matching in the solution after consecutive iterations:\n";
    vector<pii> sol;
    vector<long double> percent;
    rep(it, 0, R) {
	    cout << stochasticMatching_avg(N, sol, p_v, p_e) << " ";
	    percent.push_back((long double)sol.size() / (long double)ed.size());
	    vector<pii> inc = stochasticMatching(N, ed, p_v, p_e);
	    
	    rep(it1, 0, (int)inc.size()) {
	        sol.push_back(inc[it1]);
	        sol.push_back({inc[it1].second, inc[it1].first});
	    }
        sort(sol.begin(), sol.end());
        sol.erase(unique(sol.begin(), sol.end()), sol.end());
	}
	cout << "\n";
	cout << "Percent of edges taken to the solution after consecutive iterations:\n";
	rep(it, 0, R)
	    cout << percent[it] << " ";
	cout << "\n";
	return sol;
}


int main() {
    cin.tie(0)->sync_with_stdio(0);
	//cin.exceptions(cin.failbit);
	rng = mt19937(chrono::steady_clock::now().time_since_epoch().count());
	
	//input format:
	//first line: n - number of vertices; m - number of edges
	//next m lines: u v w - an edge from u to v with weight w (always 1)
	//
	//next 2n lines: description of probabilities of realization
	//line 2i: prob. of realizing vertex i
	//line 2i + 1: list of prob. of realizing edges incident to i
	
    int n, m;
    int a, b, weight;
    long double prob;
    
	vector<pii> ed;
	cin >> n >> m;
	vector<vector<int>> nb(n);
	rep(it, 0, m) {
	    cin >> a >> b >> weight;
	    ed.push_back({a, b}); 
	    nb[a].push_back(b);
    }
	vector<long double> p_v;
	vector<vector<long double>> p_e(n);
	rep(it, 0, n) {
	    cin >> prob;
	    p_v.push_back(1 - prob);
	    p_e[it].resize(n);
	    rep(it1, 0, (int)nb[it].size()) {
	        cin >> prob;
	        p_e[it][nb[it][it1]] = 1 - prob;
	    }
	}
	
	int R = 80; //number of iterations
	
	cout << "Matching size of the original graph: " << generalMatching(n, ed).size() << "\n";
	cout << "Est expected matching of the realized graph: " << stochasticMatching_avg(n, ed, p_v, p_e) << "\n";
	
	stochastic_inc(n, ed, p_v, p_e, R);
	return 0;
}
