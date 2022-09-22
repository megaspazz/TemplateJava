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

	/*
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
				AVLTreeNode[] path = getLeafToRootPath(k);
				path[0].val = v;
				updatePath(path);
			}

			public void increment(long k, long v) {
				AVLTreeNode[] path = getLeafToRootPath(k);
				path[0].val += v;
				updatePath(path);
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

			private AVLTreeNode[] getLeafToRootPath(long k) {
				ArrayList<AVLTreeNode> lst = new ArrayList<>();
				AVLTreeNode curr = this;
				lst.add(curr);
				while (curr.key != k) {
					if (k < curr.key) {
						curr = curr.left = getOrCreate(curr.left, k);
					} else {
						curr = curr.right = getOrCreate(curr.right, k);
					}
					lst.add(curr);
				}
				Collections.reverse(lst);
				return lst.toArray(new AVLTreeNode[0]);
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

			private static void updatePath(AVLTreeNode[] path) {
				for (AVLTreeNode node : path) {
					node.left = rebalance(node.left);
					node.right = rebalance(node.right);
					node.update();
				}
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
}
