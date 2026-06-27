#ifndef RMQ_SDSL_FAST_GUARD
#define RMQ_SDSL_FAST_GUARD

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

template<uint32_t t_sparseTable_block_size, uint32_t t_bitmask_size, uint32_t... t_recurisve_sizes>
class RMQ_SDSL_Fast
{
	using recursive_rmq = RMQ_SDSL_Fast<t_sparseTable_block_size, t_recurisve_sizes..., 0>;
public:
    typedef typename bit_vector::size_type size_type;
private:
	static constexpr int B = 32; // not larger!
	std::unique_ptr<sparse_table> m_sparse_table;
	std::unique_ptr<recursive_rmq> m_recursive_rmq;
	int_vector<> m;
	int_vector<> a, c, c_indexes;
	int_vector<> block_minimums;
	int_vector<> block_minumums_indices;
	int_vector<> m_sample_idx, m_sample_val;

	template<class t_rac>
	void build_sparse_table(const t_rac* v) {
		size_type n = v->size();
		size_type block_size = t_sparseTable_block_size;
		size_type sample_size = n/block_size + ((n % block_size) != 0);
		m_sample_idx = int_vector<>(sample_size);
		m_sample_val = int_vector<>(sample_size);
		for(size_type i = 0; i < sample_size; ++i) {
			size_type min_idx = i * block_size;
			for(size_type j = i*block_size; j < std::min((i+1)*block_size,n); ++j) {
				if((*v)[j] < (*v)[min_idx]) min_idx = j;
			}
			m_sample_val[i] = (*v)[min_idx];
			m_sample_idx[i] = min_idx - i * block_size;
		}
		util::bit_compress(m_sample_idx);
		util::bit_compress(m_sample_val);
		m_sparse_table = std::make_unique<sparse_table>(&m_sample_val);
	}

public:
	template<class t_rac>
	RMQ_SDSL_Fast(const t_rac* A=nullptr) : m_sparse_table(nullptr), m_recursive_rmq(nullptr)
	{
		if (A == nullptr)
		{
			return;
		}
		if (t_bitmask_size == 0)
		{
			block_minimums = int_vector<>(A->size(),0);
			for (size_t i = 0; i < A->size(); ++i) 
				block_minimums[i] = (*A)[i];
			m_sparse_table = std::make_unique<sparse_table>(A);
			return;
		}
		m = int_vector<>(sz(*A), 0);
		a = int_vector<>(*A);
		c = int_vector<>(sz(*A), 0);
		c_indexes = int_vector<>(sz(*A), 0);
		block_minimums = int_vector<>(sz(*A) / B + 1, 0);
		block_minumums_indices = int_vector<>(sz(*A) / B + 1, 0);
		uint32_t mi = 0;
		rep(i, sz(a)) {
			if (!(i % B) or a[i] < block_minimums[i / B])
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
		util::bit_compress(block_minimums);
		util::bit_compress(block_minumums_indices);
		m_recursive_rmq = std::make_unique<recursive_rmq>(&block_minimums);
		if (t_sparseTable_block_size)
		{
			build_sparse_table(A);
		}
	}
	size_type operator()(size_type l, size_type r) const
	{
		if (t_bitmask_size == 0)
		{
			return (*m_sparse_table)(l, r);
		}
		if (t_sparseTable_block_size > 0)
		{
			const size_type block_size = t_sparseTable_block_size;
			size_type i = l / block_size;
			size_type j = r / block_size;
			size_type min_block = (*m_sparse_table)(i,j);
			size_type min_idx = m_sample_idx[min_block] + min_block * block_size;
			if(l <= min_idx and min_idx <= r)
			{
				return min_idx;
			}
		}
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
		size_type recursive_min_idx = (*m_recursive_rmq)(l, r);
		size_type recursive_min = block_minimums[recursive_min_idx];
		recursive_min_idx = block_minumums_indices[recursive_min_idx];
		if (recursive_min < local_min or (recursive_min == local_min and recursive_min_idx < local_idx))
		{
			return recursive_min_idx;
		}
		else
		{
			return local_idx;
		}
	}

	int getFunctionAnswered() {
		return 0;
	}

    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="") const 
    {
        structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
		size_type written_bytes = 0;
		if(t_bitmask_size) {
			written_bytes += m.serialize(out, child, "m");
			written_bytes += c.serialize(out, child, "c");
			written_bytes += c_indexes.serialize(out, child, "c_indexes");
			written_bytes += block_minimums.serialize(out, child, "block_minimums");
			written_bytes += block_minumums_indices.serialize(out, child, "block_minimums_indices");
			if(t_sparseTable_block_size) {
				written_bytes += m_sample_idx.serialize(out, child, "sample_idx");
				written_bytes += m_sample_val.serialize(out, child, "sample_val");
				written_bytes += m_sparse_table->serialize(out, child, "sparse_table");
			}
			written_bytes += m_recursive_rmq->serialize(out, child, "rmq_recursive");
		} else {
			written_bytes += block_minimums.serialize(out, child, "block_minimums");
			written_bytes += m_sparse_table->serialize(out, child, "sparse_rmq");
		}
		structure_tree::add_size(child, written_bytes);
		return written_bytes;
    }
};
} // namespace sdsl

#endif // RMQ_SDSL_FAST_GUARD