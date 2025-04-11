### given array a find its prefix sum array by convoluting it with an array of the same size full of ones

### why?

### because each index will correspond to itself and all after it when performing the multiplication

### now say we want to find the kth prefix sum array 

### what do I mean

### For example,

If  a0 = {1, 2, 1, 3}

a1  = {1, 3, 4, 7}
    
a2 = {1, 4, 8, 15}

a3 = {1, 5, 13, 28}  

### all we want to do is raise the array full of ones to the power k and multiply it by a

### implementation:

```
//بسم الله الرحمن الرحيم
#include <bits/stdc++.h>
#define ll int

#define ordered_set tree<ll, null_type,less_equal<ll>, rb_tree_tag,tree_order_statistics_node_update>
using namespace std;
#include <tr2/dynamic_bitset>
using dynamic_bitset = tr2::dynamic_bitset<>;
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
#include <ext/pb_ds/assoc_container.hpp>
using namespace __gnu_pbds;
using namespace std;
#define ld long double
#define fr first;
#define sc second;
#define pii pair<ll, ll>
#define all(x) x.begin(), x.end()
#define FIO ios_base::sync_with_stdio(0), cin.tie(0), cout.tie(0);
#define is(V) (V < M and is_prime[V])
#define LF (x*2+1)
#define RT (x*2+2)
#define md ((lx+rx) >> 1)
#define vi vector<ll>
#define vv vector<vector<ll>>
#define vvii vector<vector<pii>>
#define vb vector<bool>
#define mem(x) memset(x, -1, sizeof x)
#define db dynamic_bitset
#define print(v) for (auto& i : v) cout << i << " ";
const int N = 3e5 + 9, mod = 7340033;

struct base {
    double x, y;
    base() { x = y = 0; }
    base(double x, double y): x(x), y(y) { }
};
inline base operator + (base a, base b) { return base(a.x + b.x, a.y + b.y); }
inline base operator - (base a, base b) { return base(a.x - b.x, a.y - b.y); }
inline base operator * (base a, base b) { return base(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x); }
inline base conj(base a) { return base(a.x, -a.y); }
int lim = 1;
vector<base> roots = {{0, 0}, {1, 0}};
vector<int> rev = {0, 1};
const double PI = acosl(- 1.0);
void ensure_base(int p) {
    if(p <= lim) return;
    rev.resize(1 << p);
    for(int i = 0; i < (1 << p); i++) rev[i] = (rev[i >> 1] >> 1) + ((i & 1)  <<  (p - 1));
    roots.resize(1 << p);
    while(lim < p) {
        double angle = 2 * PI / (1 << (lim + 1));
        for(int i = 1 << (lim - 1); i < (1 << lim); i++) {
            roots[i << 1] = roots[i];
            double angle_i = angle * (2 * i + 1 - (1 << lim));
            roots[(i << 1) + 1] = base(cos(angle_i), sin(angle_i));
        }
        lim++;
    }
}
void fft(vector<base> &a, int n = -1) {
    if(n == -1) n = a.size();
    assert((n & (n - 1)) == 0);
    int zeros = __builtin_ctz(n);
    ensure_base(zeros);
    int shift = lim - zeros;
    for(int i = 0; i < n; i++) if(i < (rev[i] >> shift)) swap(a[i], a[rev[i] >> shift]);
    for(int k = 1; k < n; k <<= 1) {
        for(int i = 0; i < n; i += 2 * k) {
            for(int j = 0; j < k; j++) {
                base z = a[i + j + k] * roots[j + k];
                a[i + j + k] = a[i + j] - z;
                a[i + j] = a[i + j] + z;
            }
        }
    }
}
//eq = 0: 4 FFTs in total
//eq = 1: 3 FFTs in total
vector<int> multiply(vector<int> &a, vector<int> &b, ll lim, int eq = 0) {
    int need = a.size() + b.size() - 1;
    int p = 0;
    while((1 << p) < need) p++;
    ensure_base(p);
    int sz = 1 << p;
    vector<base> A, B;
    if(sz > (int)A.size()) A.resize(sz);
    for(int i = 0; i < (int)a.size(); i++) {
        int x = (a[i] % mod + mod) % mod;
        A[i] = base(x & ((1 << 15) - 1), x >> 15);
    }
    fill(A.begin() + a.size(), A.begin() + sz, base{0, 0});
    fft(A, sz);
    if(sz > (int)B.size()) B.resize(sz);
    if(eq) copy(A.begin(), A.begin() + sz, B.begin());
    else {
        for(int i = 0; i < (int)b.size(); i++) {
            int x = (b[i] % mod + mod) % mod;
            B[i] = base(x & ((1 << 15) - 1), x >> 15);
        }
        fill(B.begin() + b.size(), B.begin() + sz, base{0, 0});
        fft(B, sz);
    }
    double ratio = 0.25 / sz;
    base r2(0,  - 1), r3(ratio, 0), r4(0,  - ratio), r5(0, 1);
    for(int i = 0; i <= (sz >> 1); i++) {
        int j = (sz - i) & (sz - 1);
        base a1 = (A[i] + conj(A[j])), a2 = (A[i] - conj(A[j])) * r2;
        base b1 = (B[i] + conj(B[j])) * r3, b2 = (B[i] - conj(B[j])) * r4;
        if(i != j) {
            base c1 = (A[j] + conj(A[i])), c2 = (A[j] - conj(A[i])) * r2;
            base d1 = (B[j] + conj(B[i])) * r3, d2 = (B[j] - conj(B[i])) * r4;
            A[i] = c1 * d1 + c2 * d2 * r5;
            B[i] = c1 * d2 + c2 * d1;
        }
        A[j] = a1 * b1 + a2 * b2 * r5;
        B[j] = a1 * b2 + a2 * b1;
    }
    fft(A, sz); fft(B, sz);
    vector<int> res(min(need,lim));
    for(int i = 0; i < min(need,lim); i++) {
        long long aa = A[i].x + 0.5;
        long long bb = B[i].x + 0.5;
        long long cc = A[i].y + 0.5;
        res[i] = (aa + ((bb % mod) << 15) + ((cc % mod) << 30))%mod;
    }
    return res;
}

vector<int> pow(vector<int>& a, int p) {
    vector<int> res;
    res.emplace_back(1);
    while(p) {
        if(p & 1) res = multiply(res, a, a.size());
        a = multiply(a, a, a.size());
        p >>= 1;
    }
    return res;
}

void solve(){
    ll n, k;
    cin >> n >> k;
    vi a(n);
    vi b(n,1);
    for (auto& i : a) cin >> i;
    vi p = pow(b,k);
    vi ans = multiply(a,p,a.size());
    for (ll i =0; i < n; i++) cout << ans[i] << " ";
    cout << "\n";
}

int main() {
    ll t = 1;
    //cin >> t;
    while (t--) {
        solve();
    }
}

```



