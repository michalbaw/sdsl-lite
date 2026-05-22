#ifndef RMQ_SDSL_BITMASKS_COMPRESSED
#define RMQ_SDSL_BITMASKS_COMPRESSED

#include <bits/stdc++.h>
#include <limits>
#include <optional>
#include <utility>
#include <print>

#include "rmq_support.hpp"
#include "rmq_sdsl_fast.hpp"
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

using return_size_type = bit_vector::size_type;

struct rmq_sdsl_bitmasks_return_type
{
	std::optional<return_size_type> left_candidate;
	std::optional<return_size_type> middle_candidate;
	std::optional<return_size_type> right_candidate;
};

// template<bool t_min = true, bool t_use_bitmasks = false, uint32_t t_st_block_size = 0, uint32_t t_super_block_size=1024, uint32_t... t_block_sizes>
// class rmq_succinct_rec_new;

template<uint32_t t_bitmask_size, uint32_t... t_recurisve_sizes>
class RMQ_SDSL_Bitmasks_Compressed
{
	using recursive_rmq = RMQ_SDSL_Fast<0, t_recurisve_sizes...>;
	// using recursive_rmq = rmq_succinct_rec_new<true, true, 0, 1024,128,0>;
public:
    typedef typename bit_vector::size_type size_type;
private:
	static constexpr int B = t_bitmask_size;
	using bitmask_type = uint32_t;
	std::unique_ptr<recursive_rmq> m_recursive_rmq;
	std::unique_ptr<sparse_table> m_sparse_table;
	int_vector<> masks_right;
	int_vector<> masks_left;
	// int_vector<> input_array;
	int_vector<> block_minimums;
	int_vector<> block_minumums_indices;
	int_vector<> right_neighbour_lower_bound;
	int_vector<> left_neighbour_lower_bound;

public:
	template<class t_rac>
	RMQ_SDSL_Bitmasks_Compressed(const t_rac* A=nullptr) : m_recursive_rmq(nullptr)
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
		masks_left = int_vector<>(sz(*A)/B + 1, 0);
		masks_right = int_vector<>(sz(*A)/B + 1, 0);
		right_neighbour_lower_bound = int_vector<>(sz(*A)/B + 1, 0);
		left_neighbour_lower_bound = int_vector<>(sz(*A)/B + 1, 0);
		int_vector<> input_array = int_vector<>(*A);
		block_minimums = int_vector<>(sz(*A) / B + 1, 0);
		block_minumums_indices = int_vector<>(sz(*A) / B + 1, 0);
		bitmask_type mi = 0;
		for (int i = 0; i < sz(input_array); i++)
		{
			if ((i % B) == 0 or i == sz(input_array) - 1 or input_array[i] < block_minimums[i / B])
			{
				block_minimums[i / B] = input_array[i];
				block_minumums_indices[i / B] = i;
			}
			mi <<= 1;
			while (mi && input_array[i] < input_array[i - __builtin_ctzll(mi)])
			{
				mi ^= (1u << __builtin_ctzll(mi));
			}
			mi ^= 1;
			if ((i % B) == B - 1)
			{
				masks_left[i/B] = mi;
			}
		}
		mi = 0;
		for (int i = sz(input_array) - 1; i >= 0; i--)
		{
			mi <<= 1;
			while (mi && input_array[i] <= input_array[i + __builtin_ctzll(mi)])
			{
				mi ^= (1u << __builtin_ctzll(mi));
			}
			mi ^= 1;
			if ((i % B) == 0) {
				masks_right[i/B] = mi;
			}
		}
		size_type neighbour_min = 0;
		for (int i = B; i < sz(input_array); i++)
		{
			if ((i%B) == 0)
			{
				neighbour_min = block_minimums[i/B - 1];
			}
			if (input_array[i] < neighbour_min)
			{
				left_neighbour_lower_bound[i/B] = i;
				i = (i / B) * B + B - 1;
				continue;
			}
			if ((i%B) == B-1)
			{
				left_neighbour_lower_bound[i/B] = i;
			}
		}
		for (int i = (sz(input_array) / B) * B - 1; i >= 0; i--)
		{
			if ((i%B) == B - 1)
			{
				neighbour_min = block_minimums[i/B + 1];
			}
			if (input_array[i] <= neighbour_min)
			{
				right_neighbour_lower_bound[i/B] = i;
				i = (i / B) * B;
				continue;
			}
			if ((i%B) == 0)
			{
				right_neighbour_lower_bound[i/B] = i;
			}
		}
		util::bit_compress(block_minimums);
		util::bit_compress(block_minumums_indices);
		m_recursive_rmq = std::make_unique<recursive_rmq>(&block_minimums);
	}

	rmq_sdsl_bitmasks_return_type operator()(size_type l, size_type r) const
	{
		if (t_bitmask_size == 0)
		{
			return {.middle_candidate = 0};
		}
		if (r - l + 1 < 3 * B)
		{
			return {.middle_candidate = 0};
		}

		size_type left_block = l / B;
		size_type right_block = r / B;

		{
			auto recursive_min_idx = m_recursive_rmq->operator()(left_block, right_block);
			recursive_min_idx = block_minumums_indices[recursive_min_idx];
			if (recursive_min_idx >= l and recursive_min_idx <= r)
			{
				return {.middle_candidate = recursive_min_idx};
			}
		}

		size_type right_edge = (r / B) * B;
		size_type right_idx = right_edge + __lg(static_cast<bitmask_type>(masks_right[right_edge / B]) & ((1lu << (r - right_edge + 1)) - 1));
		// size_type right_val = input_array[right_idx];

		size_type left_edge = ((l + B) / B) * B - 1;
		size_type left_idx = left_edge - __lg(static_cast<bitmask_type>(masks_left[left_edge / B]) & ((1lu << (left_edge - l + 1)) - 1));
		// size_type left_val = input_array[left_idx];

		size_type left_inside_bound = (l + B - 1) / B;
		size_type right_inside_bound = (r + 1) / B - 1;
		rmq_sdsl_bitmasks_return_type returned_value;

		if (left_idx <= right_neighbour_lower_bound[left_block])
		{
			returned_value.left_candidate = left_idx;
		}

		if (right_idx >= left_neighbour_lower_bound[right_block])
		{
			returned_value.right_candidate = right_idx;
		}

		size_type recursive_min_idx = m_recursive_rmq->operator()(left_inside_bound, right_inside_bound);
		// size_type recursive_min = block_minimums[recursive_min_idx];
		recursive_min_idx = block_minumums_indices[recursive_min_idx];
		returned_value.middle_candidate = recursive_min_idx;

		// size_type true_idx, true_val = numeric_limits<unsigned long>::max();
		// for (int i = l; i <= r; i++)
		// {
		// 	if (input_array[i] < true_val)
		// 	{
		// 		true_idx = i;
		// 		true_val = input_array[i];
		// 	}
		// }

		// if (true_idx == left_idx and not returned_value.left_candidate.has_value() or true_idx == right_idx and not returned_value.right_candidate.has_value())
		// {
		// 	std::println(std::cerr, "Oops");
		// 	std::println(std::cerr, "l: {}, r: {}, true: {}->{}, left: {}->{}, mid: {}->{}, right: {}->{}",
		// 		l, r, true_idx, true_val, left_idx, left_val, recursive_min_idx, recursive_min, right_idx, right_val);
		// 	std::println(std::cerr, "Left is present: {}, right: {}", returned_value.left_candidate.has_value(), returned_value.right_candidate.has_value());
		// 	std::println(std::cerr, "Left cutoff idx: {}, right cutoff idx: {}", right_neighbour_lower_bound[left_block], left_neighbour_lower_bound[right_block]);
		// }

		return returned_value;
	}

    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="") const 
    {
        structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
		size_type written_bytes = 0;
		if(t_bitmask_size) {
			written_bytes += masks_left.serialize(out, child, "masks_left");
			written_bytes += masks_right.serialize(out, child, "masks_right");
			// written_bytes += input_array.serialize(out, child, "a");
			written_bytes += block_minimums.serialize(out, child, "block_minimums");
			// written_bytes += block_minumums_indices.serialize(out, child, "block_minimums_indices");
			written_bytes += m_recursive_rmq->serialize(out, child, "rmq_recursive");
		} else {
			written_bytes += block_minimums.serialize(out, child, "block_minimums");
		}
		structure_tree::add_size(child, written_bytes);
		return written_bytes;
    }
};
} // namespace sdsl

#endif // RMQ_SDSL_BITMASKS_COMPRESSED