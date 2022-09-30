import java.util.*;

public class RangeQueries {
	/*
	 * Zero-indexed Fenwick Tree to compute range sums.
	 * You need ceil(log2(N)) bits to store values N values from the range [0, N).
	 * Commonly used bits values:
	 *   - For N = 10^3, need bits = 10.
	 *   - For N = 10^4, need bits = 14.
	 *   - For N = 10^5, need bits = 17.
	 *   - For N = 10^6, need bits = 20.
	 */
	public static class LongFenwickTreeRangeSum {
		private long[] partial;

		public LongFenwickTreeRangeSum(int bits) {
			partial = new long[1 << bits];
		}

		public void insert(int index, long value) {
			long delta = value - get(index, index);
			increment(index, delta);
		}

		public void increment(int index, long value) {
			int curr = index;
			while (curr < partial.length) {
				partial[curr] += value;
				curr += Integer.lowestOneBit(~curr);
			}
		}

		public long get(int index) {
			return get(index, index);
		}

		public long get(int loInclusive, int hiInclusive) {
			long sum = prefixSum(hiInclusive);
			if (loInclusive > 0) {
				sum -= prefixSum(loInclusive - 1);
			}
			return sum;
		}

		private long prefixSum(int hiInclusive) {
			long sum = 0;
			int curr = hiInclusive;
			while (curr >= 0) {
				sum += partial[curr];
				curr -= Integer.lowestOneBit(~curr);
			}
			return sum;
		}

		public static LongFenwickTreeRangeSum newWithSize(int size) {
			return new LongFenwickTreeRangeSum(Integer.SIZE - Integer.numberOfLeadingZeros(size));
		}
	}

	/*
	 * Zero-indexed Fenwick Tree to compute range sums.
	 * You need ceil(log2(N)) bits to store values N values from the range [0, N).
	 * Commonly used bits values:
	 *   - For N = 10^3, need bits = 10.
	 *   - For N = 10^4, need bits = 14.
	 *   - For N = 10^5, need bits = 17.
	 *   - For N = 10^6, need bits = 20.
	 */
	public static class IntFenwickTreeRangeSum {
		private int[] partial;

		public IntFenwickTreeRangeSum(int bits) {
			partial = new int[1 << bits];
		}

		public void insert(int index, int value) {
			int delta = value - get(index, index);
			increment(index, delta);
		}

		public void increment(int index, int value) {
			int curr = index;
			while (curr < partial.length) {
				partial[curr] += value;
				curr += Integer.lowestOneBit(~curr);
			}
		}

		public int get(int index) {
			return get(index, index);
		}

		public int get(int loInclusive, int hiInclusive) {
			int sum = prefixSum(hiInclusive);
			if (loInclusive > 0) {
				sum -= prefixSum(loInclusive - 1);
			}
			return sum;
		}

		private int prefixSum(int hiInclusive) {
			int sum = 0;
			int curr = hiInclusive;
			while (curr >= 0) {
				sum += partial[curr];
				curr -= Integer.lowestOneBit(~curr);
			}
			return sum;
		}

		public static IntFenwickTreeRangeSum newWithSize(int size) {
			return new IntFenwickTreeRangeSum(Integer.SIZE - Integer.numberOfLeadingZeros(size));
		}
	}

	/**
	 * Here is a generic Segment Tree implementation.
	 * It requires a Combiner to define how to merge values.
	 * It takes an optional DefaultProvider to replace null values; otherwise, the Combiner will need to manually handle null values.
	 *
	 * NOTE: Slightly prefer ArraySegmentTree for performance!
	 */
	public static class GenericSegmentTree<T> {
		public ArrayList<SegmentTreeNode> leaves;
		public SegmentTreeNode root;
		public Combiner<T> combiner;
		public DefaultProvider<T> defaultProvider;

		public GenericSegmentTree(int n, Combiner<T> cmb) {
			this(n, cmb, null);
		}

		public GenericSegmentTree(int n, Combiner<T> cmb, DefaultProvider<T> defProv) {
			this.combiner = cmb;
			this.defaultProvider = defProv;
			this.leaves = new ArrayList<SegmentTreeNode>(n);
			for (int i = 0; i < n; ++i) {
				this.leaves.add(null);
			}
			this.root = new SegmentTreeNode(null, 0, n - 1);
		}

		public GenericSegmentTree(T[] vals, Combiner<T> cmb, DefaultProvider<T> defProv) {
			this(vals.length, cmb, defProv);
			for (int i = 0; i < vals.length; i++) {
				this.insert(i, vals[i]);
			}
		}

		public void insert(int idx, T v) {
			this.leaves.get(idx).setAndUpdate(v);
		}

		public T get(int idx) {
			return this.leaves.get(idx).val;
		}

		public T get(int lowerInclusive, int upperInclusive) {
			return this.root.getRange(lowerInclusive, upperInclusive);
		}

		public static interface Combiner<T> {
			public T combine(T lhs, T rhs);
		}

		public static interface DefaultProvider<T> {
			public T getDefault();
		}

		private class SegmentTreeNode {
			public int L;
			public int R;

			public T val;

			public SegmentTreeNode parent;
			public SegmentTreeNode left;
			public SegmentTreeNode rite;

			public SegmentTreeNode(SegmentTreeNode p, int lower, int upper) {
				this.parent = p;
				this.L = lower;
				this.R = upper;

				if (lower == upper) {
					// access outer class GenericSegmentTree
					leaves.set(lower, this);
				} else {
					int mid = (lower + upper) / 2;
					this.left = new SegmentTreeNode(this, lower, mid);
					this.rite = new SegmentTreeNode(this, mid + 1, upper);
				}
			}

			public void setAndUpdate(T v) {
				this.val = v;
				this.update();
			}

			public void update() {
				if (this.left != null && this.rite != null) {
					// access outer class GenericSegmentTree
					this.val = combiner.combine(this.left.getValueOrDefault(), this.rite.getValueOrDefault());
				} else if (this.left != null) {
					this.val = this.left.getValueOrDefault();
				} else if (this.rite != null) {
					this.val = this.rite.getValueOrDefault();
				}

				if (this.parent != null) {
					this.parent.update();
				}
			}

			public T getRange(int lower, int upper) {
				if (this.L >= lower && this.R <= upper) {
					return getValueOrDefault();
				} else if (this.L > upper || this.R < lower) {
					// access outer class GenericSegmentTree
					return defaultProvider.getDefault();
				} else {
					// access outer class GenericSegmentTree
					return combiner.combine(this.left.getRange(lower, upper), this.rite.getRange(lower, upper));
				}
			}
			
			private T getValueOrDefault() {
				if (val != null) {
					return val;
				}
				if (defaultProvider != null) {
					return defaultProvider.getDefault();
				}
				return null;
			}
		}
	}

	public static class ArraySegmentTree<T> {
		private int bits;
		private ArrayList<T> values = new ArrayList<>();
		
		private Merger<T> merger;
		private T defaultValue;
		
		public ArraySegmentTree(int bits, Merger<T> merger, T defaultValue) {
			this.bits = bits;
			this.merger = merger;
			this.defaultValue = defaultValue;
			
			int nodeCount = 1 << (bits + 1);
			while (values.size() < nodeCount) {
				values.add(defaultValue);
			}
		}
		
		public void insert(int index, T data) {
			int curr = (values.size() >> 1) + index;
			values.set(curr, data);
			curr >>= 1;
			while (curr > 0) {
				int left = curr << 1;
				int right = left + 1;
				values.set(curr, merger.merge(values.get(left), values.get(right)));
				curr >>= 1;
			}
		}
		
		public T get(int loInclusive, int hiInclusive) {
			if (loInclusive > hiInclusive) {
				return defaultValue;
			}
			
			int curr = 1;
			for (int d = 0; d < bits ; ++d) {
				int shift = bits - d;
				int LL = (curr << shift) - (values.size() >> 1);
				int LR = LL + (1 << (shift - 1)) - 1;
				int RL = LR + 1;
				int lCurr = curr << 1;
				int rCurr = lCurr + 1;
				if (hiInclusive <= LR) {
					curr = lCurr;
				} else if (loInclusive >= RL) {
					curr = rCurr;
				} else {
					T leftValue = getGTE(loInclusive, lCurr, d + 1);
					T rightValue = getLTE(hiInclusive, rCurr, d + 1);
					return merger.merge(leftValue, rightValue);
				}
			}
			return values.get(curr);
		}
		
		private T getGTE(int loInclusive, int curr, int dStart) {
			T ans = defaultValue;
			for (int d = dStart; d < bits; ++d) {
				int shift = bits - d;
				int LL = (curr << shift) - (values.size() >> 1);
				int LR = LL + (1 << (shift - 1)) - 1;
				int RL = LR + 1;
				int rCurr = (curr << 1) + 1;
				if (loInclusive <= LL) {
					break;
				}
				if (loInclusive >= RL) {
					curr = rCurr;
				} else {
					ans = merger.merge(ans, values.get(rCurr));
					curr <<= 1;
				}
			}
			return merger.merge(ans, values.get(curr));
		}
		
		private T getLTE(int hiInclusive, int curr, int dStart) {
			T ans = defaultValue;
			for (int d = dStart; d < bits; ++d) {
				int shift = bits - d;
				int LL = (curr << shift) - (values.size() >> 1);
				int LR = LL + (1 << (shift - 1)) - 1;
				int RR = LL + (1 << shift) - 1;
				int lCurr = curr << 1;
				if (hiInclusive >= RR) {
					break;
				}
				if (hiInclusive <= LR) {
					curr = lCurr;
				} else {
					ans = merger.merge(ans, values.get(lCurr));
					curr = lCurr + 1;
				}
			}
			return merger.merge(ans, values.get(curr));
		}

		public static <T> ArraySegmentTree<T> newWithSize(int size, Merger<T> merger, T defaultValue) {
			return new ArraySegmentTree<T>(Integer.SIZE - Integer.numberOfLeadingZeros(size), merger, defaultValue);
		}
		
		public static interface Merger<T> {
			public T merge(T a, T b);
		}
	}

	/**
	 * AVLTreeRangeSum performs arbitrary-length range sums in O(log N) time.
	 */
	public static class AVLTreeRangeSum {
		private AVLTreeNode root = new AVLTreeNode(0);

		public void insert(long k, long v) {
			root.insert(k, v);
			root = AVLTreeNode.rebalance(root);
		}

		public void increment(long k, long v) {
			root.increment(k, v);
			root = AVLTreeNode.rebalance(root);
		}

		public long get(long k) {
			return root.get(k);
		}

		public long get(long lo, long hi) {
			return root.get(lo, hi);
		}
		
		public static class AVLTreeNode {
			private long key;
			private long val;
			private long sum;
			private int ht;

			private AVLTreeNode left;
			private AVLTreeNode right;

			public AVLTreeNode(long k) {
				this.key = k;
			}

			public void insert(long k, long v) {
				if (k < key) {
					left = getOrCreate(left, k);
					left.insert(k, v);
				} else if (k > key) {
					right = getOrCreate(right, k);
					right.insert(k, v);
				} else {
					val = v;
				}
				left = rebalance(left);
				right = rebalance(right);
				update();
			}

			public void increment(long k, long v) {
				if (k < key) {
					left = getOrCreate(left, k);
					left.increment(k, v);
				} else if (k > key) {
					right = getOrCreate(right, k);
					right.increment(k, v);
				} else {
					val += v;
				}
				left = rebalance(left);
				right = rebalance(right);
				update();
			}

			public long get(long k) {
				return get(k, k);
			}

			public long get(long lo, long hi) {
				AVLTreeNode curr = this;
				while (curr != null) {
					if (lo <= curr.key && curr.key <= hi) {
						break;
					}
					if (hi < curr.key) {
						curr = curr.left;
					} else if (lo > curr.key) {
						curr = curr.right;
					}
				}
				if (curr == null) {
					return 0;
				}

				long ans = curr.val;
				if (curr.left != null) {
					ans += curr.left.getSumGTE(lo);
				}
				if (curr.right != null) {
					ans += curr.right.getSumLTE(hi);
				}
				return ans;
			}

			private static AVLTreeNode rotateRight(AVLTreeNode root) {
				AVLTreeNode pivot = root.left;
				root.left = pivot.right;
				pivot.right = root;
				root.update();
				pivot.update();
				return pivot;
			}

			private static AVLTreeNode rotateLeft(AVLTreeNode root) {
				AVLTreeNode pivot = root.right;
				root.right = pivot.left;
				pivot.left = root;
				root.update();
				pivot.update();
				return pivot;
			}

			private void update() {
				computeHeight();
				computeSum();
			}

			private int computeHeight() {
				ht = 1 + Math.max(getHeight(left), getHeight(right));
				return ht;
			}

			private long computeSum() {
				sum = val + getSum(left) + getSum(right);
				return sum;
			}

			private long getSumLTE(long k) {
				AVLTreeNode curr = this;
				long sum = 0;
				while (curr != null) {
					if (k < curr.key) {
						curr = curr.left;
					} else {
						sum += curr.val + getSum(curr.left);
						curr = curr.right;
					}
				}
				return sum;
			}

			private long getSumGTE(long k) {
				AVLTreeNode curr = this;
				long sum = 0;
				while (curr != null) {
					if (k > curr.key) {
						curr = curr.right;
					} else {
						sum += curr.val + getSum(curr.right);
						curr = curr.left;
					}
				}
				return sum;
			}

			private static AVLTreeNode rebalance(AVLTreeNode node) {
				if (node == null) {
					return null;
				}
				int bf = balanceFactor(node);
				if (bf > 1) {
					return rotateRight(node);
				} else if (bf < -1) {
					return rotateLeft(node);
				} else {
					return node;
				}
			}

			private static int balanceFactor(AVLTreeNode root) {
				return getHeight(root.left) - getHeight(root.right);
			}

			private static AVLTreeNode getOrCreate(AVLTreeNode node, long k) {
				if (node != null) {
					return node;
				}
				return new AVLTreeNode(k);
			}

			private static int getHeight(AVLTreeNode node) {
				if (node == null) {
					return 0;
				}
				return node.ht;
			}

			private static long getSum(AVLTreeNode node) {
				if (node == null) {
					return 0;
				}
				return node.sum;
			}
		}
	}

	/**
	 * AVLSegmentTree does generic range operations in O(log N) time.
	 * See AVLSegmentTree.LongSum for example usage.
	 * 
	 * NOTE: If possible, use non-generic versions, since they run significantly faster.
	 */
	public static class AVLSegmentTree<KT, VT> {
		private AVLTreeNode root;
		
		private Comparator<KT> comparer;
		private Merger<VT> merger;
		private DefaultProvider<VT> defaultProvider;
		
		public AVLSegmentTree(KT initKey, Comparator<KT> comparer, Merger<VT> merger, DefaultProvider<VT> defaultProvider) {
			this.root = new AVLTreeNode(initKey);
			this.comparer = comparer;
			this.merger = merger;
			this.defaultProvider = defaultProvider;
		}

		public void insert(KT key, VT value) {
			root.insert(key, value);
			root = rebalance(root);
		}

		public VT get(KT k) {
			return root.get(k);
		}

		public VT get(KT lo, KT hi) {
			return root.get(lo, hi);
		}

		private class AVLTreeNode {
			private KT key;
			private VT val;
			private VT sum;  // more of an accumulator than a true "sum"
			private int ht;

			private AVLTreeNode left;
			private AVLTreeNode right;

			public AVLTreeNode(KT key) {
				this.key = key;
			}

			public void insert(KT k, VT v) {
				int c = comparer.compare(k, key);
				if (c < 0) {
					left = getOrCreate(left, k);
					left.insert(k, v);
				} else if (c > 0) {
					right = getOrCreate(right, k);
					right.insert(k, v);
				} else {
					val = v;
				}
				left = rebalance(left);
				right = rebalance(right);
				update();
			}

			public VT get(KT k) {
				return get(k, k);
			}

			public VT get(KT lo, KT hi) {
				AVLTreeNode curr = this;
				while (curr != null) {
					if (comparer.compare(lo, curr.key) <= 0 && comparer.compare(curr.key, hi) <= 0) {
						break;
					}
					if (comparer.compare(hi, curr.key) < 0) {
						curr = curr.left;
					} else if (comparer.compare(lo, curr.key) > 0) {
						curr = curr.right;
					}
				}
				if (curr == null) {
					return defaultProvider.getDefault();
				}
				
				VT ans = curr.getVal();
				if (curr.left != null) {
					ans = merger.merge(ans, curr.left.getSumGTE(lo));
				}
				if (curr.right != null) {
					ans = merger.merge(ans, curr.right.getSumLTE(hi));
				}
				return ans;
			}

			private void update() {
				computeHeight();
				computeSum();
			}

			private int computeHeight() {
				ht = 1 + Math.max(getHeight(left), getHeight(right));
				return ht;
			}

			private VT computeSum() {
				sum = merger.merge(getVal(), merger.merge(getSum(left), getSum(right)));
				return sum;
			}
			
			private VT getVal() {
				if (val == null) {
					return defaultProvider.getDefault();
				}
				return val;
			}

			private VT getSumLTE(KT k) {
				AVLTreeNode curr = this;
				VT sum = defaultProvider.getDefault();
				while (curr != null) {
					if (comparer.compare(k, curr.key) < 0) {
						curr = curr.left;
					} else {
						sum = merger.merge(sum, merger.merge(curr.getVal(), getSum(curr.left)));
						curr = curr.right;
					}
				}
				return sum;
			}
			
			private VT getSumGTE(KT k) {
				AVLTreeNode curr = this;
				VT sum = defaultProvider.getDefault();
				while (curr != null) {
					if (comparer.compare(k, curr.key) > 0) {
						curr = curr.right;
					} else {
						sum = merger.merge(sum, merger.merge(curr.getVal(), getSum(curr.right)));
						curr = curr.left;
					}
				}
				return sum;
			}

			private AVLTreeNode getOrCreate(AVLTreeNode node, KT k) {
				if (node != null) {
					return node;
				}
				return new AVLTreeNode(k);
			}

			private VT getSum(AVLTreeNode node) {
				if (node == null) {
					return defaultProvider.getDefault();
				}
				return node.sum;
			}
		}

		private int getHeight(AVLTreeNode node) {
			if (node == null) {
				return 0;
			}
			return node.ht;
		}

		private int balanceFactor(AVLTreeNode root) {
			return getHeight(root.left) - getHeight(root.right);
		}

		private AVLTreeNode rebalance(AVLTreeNode node) {
			if (node == null) {
				return null;
			}
			int bf = balanceFactor(node);
			if (bf > 1) {
				return rotateRight(node);
			} else if (bf < -1) {
				return rotateLeft(node);
			} else {
				return node;
			}
		}

		private AVLTreeNode rotateRight(AVLTreeNode root) {
			AVLTreeNode pivot = root.left;
			root.left = pivot.right;
			pivot.right = root;
			root.update();
			pivot.update();
			return pivot;
		}

		private AVLTreeNode rotateLeft(AVLTreeNode root) {
			AVLTreeNode pivot = root.right;
			root.right = pivot.left;
			pivot.left = root;
			root.update();
			pivot.update();
			return pivot;
		}
		
		public static interface Merger<T> {
			public T merge(T a, T b);
		}
		
		public static interface DefaultProvider<T> {
			public T getDefault();
		}
		
		/*
		 * Helper class to generate a segment tree that does sum of long values, using long key.
		 */
		public static class LongSum extends AVLSegmentTree<Long, Long> {
			public LongSum() {
				super(
					0L,
					new Comparator<Long>() {
						@Override
						public int compare(Long a, Long b) {
							return Long.compare(a, b);
						}
					},
					new Merger<Long>() {
						@Override
						public Long merge(Long a, Long b) {
							return a + b;
						}
					},
					new DefaultProvider<Long>() {
						@Override
						public Long getDefault() {
							return 0L;
						}
					}
				);
			}
			
			public void increment(long k, long v) {
				insert(k, get(k) + v);
			}
		}
	}

	/**
	 * PrimitiveAVLSegmentTree does generic range operations in O(log N) time.
	 * Use your editor's find-and-replace to rename the types into primitives, since Java doesn't support generics of primitives.
	 */
	public static class PrimitiveAVLSegmentTree {
		// Placeholder types so that the class will compile.
		// Delete these after doing find-and-replace.
		private static final class PrimitiveKeyType {}
		private static final class PrimitiveValueType {}
		
		// Implement the key comparison function, return type should be same as Comparator.compare.
		private static int compareKey(PrimitiveKeyType a, PrimitiveKeyType b) {
			throw new UnsupportedOperationException("Not implemented yet.");
			
			// Example implementation for natural ordering of long keys.
			// return Long.compare(a, b);
		}
		
		// Implement the value merge function.
		private static PrimitiveValueType mergeValue(PrimitiveValueType a, PrimitiveValueType b) {
			throw new UnsupportedOperationException("Not implemented yet.");
			
			// Example implementation for sum.
			// return a + b;
		}
		
		// The default value if nothing in range.
		private static PrimitiveValueType DEFAULT_VALUE = null;
		
		private AVLTreeNode root;

		public void insert(PrimitiveKeyType key, PrimitiveValueType value) {
			if (root == null) {
				root = new AVLTreeNode(key);
			}
			root.insert(key, value);
			root = rebalance(root);
		}

		public PrimitiveValueType get(PrimitiveKeyType k) {
			if (root == null) {
				return DEFAULT_VALUE;
			}
			return root.get(k);
		}

		public PrimitiveValueType get(PrimitiveKeyType lo, PrimitiveKeyType hi) {
			if (root == null) {
				return DEFAULT_VALUE;
			}
			return root.get(lo, hi);
		}

		private static class AVLTreeNode {
			private PrimitiveKeyType key;
			private PrimitiveValueType val;
			private PrimitiveValueType sum;  // more of an accumulator than a true "sum"
			private int ht;

			private AVLTreeNode left;
			private AVLTreeNode right;

			public AVLTreeNode(PrimitiveKeyType key) {
				this.key = key;
				this.val = DEFAULT_VALUE;
			}

			public void insert(PrimitiveKeyType k, PrimitiveValueType v) {
				int c = compareKey(k, key);
				if (c < 0) {
					left = getOrCreate(left, k);
					left.insert(k, v);
				} else if (c > 0) {
					right = getOrCreate(right, k);
					right.insert(k, v);
				} else {
					val = v;
				}
				left = rebalance(left);
				right = rebalance(right);
				update();
			}

			public PrimitiveValueType get(PrimitiveKeyType k) {
				return get(k, k);
			}

			public PrimitiveValueType get(PrimitiveKeyType lo, PrimitiveKeyType hi) {
				AVLTreeNode curr = this;
				while (curr != null) {
					if (compareKey(lo, curr.key) <= 0 && compareKey(curr.key, hi) <= 0) {
						break;
					}
					if (compareKey(hi, curr.key) < 0) {
						curr = curr.left;
					} else if (compareKey(lo, curr.key) > 0) {
						curr = curr.right;
					}
				}
				if (curr == null) {
					return DEFAULT_VALUE;
				}
				
				PrimitiveValueType ans = curr.val;
				if (curr.left != null) {
					ans = mergeValue(ans, curr.left.getSumGTE(lo));
				}
				if (curr.right != null) {
					ans = mergeValue(ans, curr.right.getSumLTE(hi));
				}
				return ans;
			}

			private void update() {
				computeHeight();
				computeSum();
			}

			private int computeHeight() {
				ht = 1 + Math.max(getHeight(left), getHeight(right));
				return ht;
			}

			private PrimitiveValueType computeSum() {
				sum = mergeValue(val, mergeValue(getSum(left), getSum(right)));
				return sum;
			}

			private PrimitiveValueType getSumLTE(PrimitiveKeyType k) {
				AVLTreeNode curr = this;
				PrimitiveValueType sum = DEFAULT_VALUE;
				while (curr != null) {
					if (compareKey(k, curr.key) < 0) {
						curr = curr.left;
					} else {
						sum = mergeValue(sum, mergeValue(curr.val, getSum(curr.left)));
						curr = curr.right;
					}
				}
				return sum;
			}
			
			private PrimitiveValueType getSumGTE(PrimitiveKeyType k) {
				AVLTreeNode curr = this;
				PrimitiveValueType sum = DEFAULT_VALUE;
				while (curr != null) {
					if (compareKey(k, curr.key) > 0) {
						curr = curr.right;
					} else {
						sum = mergeValue(sum, mergeValue(curr.val, getSum(curr.right)));
						curr = curr.left;
					}
				}
				return sum;
			}

			private AVLTreeNode getOrCreate(AVLTreeNode node, PrimitiveKeyType k) {
				if (node != null) {
					return node;
				}
				return new AVLTreeNode(k);
			}

			private PrimitiveValueType getSum(AVLTreeNode node) {
				if (node == null) {
					return DEFAULT_VALUE;
				}
				return node.sum;
			}
		}

		private static int getHeight(AVLTreeNode node) {
			if (node == null) {
				return 0;
			}
			return node.ht;
		}

		private static int balanceFactor(AVLTreeNode root) {
			return getHeight(root.left) - getHeight(root.right);
		}

		private static AVLTreeNode rebalance(AVLTreeNode node) {
			if (node == null) {
				return null;
			}
			int bf = balanceFactor(node);
			if (bf > 1) {
				return rotateRight(node);
			} else if (bf < -1) {
				return rotateLeft(node);
			} else {
				return node;
			}
		}

		private static AVLTreeNode rotateRight(AVLTreeNode root) {
			AVLTreeNode pivot = root.left;
			root.left = pivot.right;
			pivot.right = root;
			root.update();
			pivot.update();
			return pivot;
		}

		private static AVLTreeNode rotateLeft(AVLTreeNode root) {
			AVLTreeNode pivot = root.right;
			root.right = pivot.left;
			pivot.left = root;
			root.update();
			pivot.update();
			return pivot;
		}
	}
}
