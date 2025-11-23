#ifndef RMQ_FAST_GUARD
#define RMQ_FAST_GUARD

#include <bits/stdc++.h>

#include "rmq_support.hpp"
#include "int_vector.hpp"
#include "bits.hpp"
#include "util.hpp"

using namespace std;
#define fwd(i, a, n) for (int i = (a); i < (n); i++)
#define rep(i, n) fwd(i, 0, n)
#define all(X) X.begin(), X.end()
#define sz(X) int((X).size())
#define pb push_back
#define eb emplace_back
#define st first
#define nd second
using pii = pair<int, int>; using vi = vector<int>;
using ll = long long; using ld = long double;
#ifdef LOC
auto SS = signal(6, [](int) { *(int *)0 = 0; });
#define DTP(x, y) auto operator << (auto &o, auto a) -> decltype(y, o) { o << "("; x; return o << ")"; }
DTP(o << a.st << ", " << a.nd, a.nd);
DTP(for (auto i : a) o << i << ", ", all(a));
void dump(auto... x) { (( cerr << x << ", " ), ...) << '\n'; }
#define deb(x...) cerr << setw(4) << __LINE__ << ":[" #x "]: ", dump(x)
#else
#define deb(...) 0
#endif


namespace sdsl
{
using sparse_table = rmq_support_sparse_table<true,false>;

struct RMQ_Fast
{
    typedef typename bit_vector::size_type size_type;
	static constexpr int B = 32; // not larger!
	sparse_table s;
	int_vector<> m;
	int_vector<> a, c, c_indexes;
	int_vector<> block_minimums;
	int_vector<> block_minumums_indices;

	RMQ_Fast(int_vector<>* A) : m(sz(*A)), a(*A), c(sz(*A)), c_indexes(sz(*A))
	{
		block_minimums = int_vector<>(sz(a) / B + 1, 0);
		uint32_t mi = 0;
		rep(i, sz(a)) {
			if (i % B || a[i] < block_minimums[i / B])
			{
				block_minimums[i / B] = a[i];
				block_minumums_indices[i / B] = i; // TODO: keep only in-block space, and use util::bit_compress, to save space
			}
			mi <<= 1;
			while (mi && a[i] < a[i - __builtin_ctz(mi)])
			{
				mi ^= (1u << __builtin_ctz(mi));
			}
			m[i] = mi ^= 1;
			c[i] = a[i - __lg((uint32_t)m[i])];
			c_indexes[i] = i - __lg((uint32_t)m[i]);
		}
		s = sparse_table(&block_minimums);
	}
	size_type operator()(size_type l, size_type r) const
	{
		if (r - l + 1 < B)
		{
			return r - __lg(m[r] & ((1u << (r - l + 1)) - 1));
		}
		size_type local_min, local_idx;
		if (c[l + B - 1] <= c[r])
		{
			local_min = c[l + B- 1];
			local_idx = c_indexes[l + B - 1];
		}
		else
		{
			local_min = c[r];
			local_idx = c_indexes[r];
		}
		l = (l + B - 1) / B;
		r = r / B - 1;
		size_type sparse_min_idx = s(l, r);
		size_type sparse_min = block_minimums[sparse_min_idx];
		sparse_min_idx = block_minumums_indices[sparse_min_idx];
		if (sparse_min < local_min or (sparse_min == local_min and sparse_min_idx < local_idx))
		{
			return sparse_min_idx;
		}
		else
		{
			return local_idx;
		}
	}

    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="") const 
    {
        return 1;
    }
};
} // namespace sdsl

#endif // RMQ_FAST_GUARD