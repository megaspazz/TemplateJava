import java.math.*;
import java.util.*;

public class Methods {

	/**
	 * Method for computing combinations.  Returns an array of int[], which contains k elements from { 0, 1, ..., n - 1 }.
	 * Each int[] will be sorted in ascending order, and they array of int[] will be sorted ascending as well.
	 */
	public static int[][] getCombinations(int n, int k) {
		ArrayList<int[]> lst = new ArrayList<int[]>();
		int[] arr = new int[k];
		combinationHelper(n, k, 0, 0, arr, lst);
		return lst.toArray(new int[0][]);
	}

	/**
	 * Helper method for computing combinations.  Used by GetCombinations(...).
	 */
	private static void combinationHelper(int n, int k, int i, int t, int[] curr, ArrayList<int[]> lst) {
		if (t == k) {
			lst.add(curr);
			return;
		}

		curr[t] = i;

		if (k - t > n - i)
			return;

		int[] next = new int[k];
		System.arraycopy(curr, 0, next, 0, k);

		combinationHelper(n, k, i + 1, t + 1, next, lst);
		combinationHelper(n, k, i + 1, t, curr, lst);
	}

	/**
	 * Method for computing permutations.  Returns an array of int[], which contains k elements from { 0, 1, ..., n - 1 }.
	 * Each int[] will be sorted in ascending order, and they array of int[] will be sorted ascending as well.
	 */
	public static int[][] getPermutations(int n, int k) {
		ArrayList<int[]> lst = new ArrayList<int[]>();
		int[] arr = new int[k];
		boolean[] used = new boolean[n];
		permutationHelper(n, k, 0, arr, used, lst);
		return lst.toArray(new int[0][]);
	}

	/**
	 * Helper method for computing permutations.  Used by GetPermutations(...).
	 */
	private static void permutationHelper(int n, int k, int t, int[] curr, boolean[] used, ArrayList<int[]> lst) {
		if (t == k) {
			int[] copy = new int[k];
			System.arraycopy(curr, 0, copy, 0, k);
			lst.add(copy);
			return;
		}

		for (int i = 0; i < n; i++) {
			if (used[i])
				continue;

			curr[t] = i;
			used[i] = true;
			permutationHelper(n, k, t + 1, curr, used, lst);
			used[i] = false;
		}
	}

	/**
	 * Method for computing permutations with repeats.  Returns an array of int[], which contains k elements from { 0, 1, ..., n - 1 }.
	 * Each of the int[] elements will be sorted in ascending order.
	 */
	public static int[][] getAll(int n, int k) {
		ArrayList<int[]> lst = new ArrayList<int[]>();
		int[] arr = new int[k];
		allHelper(n, k, 0, arr, lst);
		return lst.toArray(new int[0][]);
	}

	/**
	 * Helper method for computing permutations with repeats.  Used by GetAll(...).
	 */
	private static void allHelper(int n, int k, int t, int[] curr, ArrayList<int[]> lst) {
		if (t == k) {
			int[] copy = new int[k];
			System.arraycopy(curr, 0, copy, 0, k);
			lst.add(copy);
			return;
		}

		for (int i = 0; i < n; i++) {
			curr[t] = i;
			allHelper(n, k, t + 1, curr, lst);
		}
	}

	/**
	 * Swaps the elements at the two indices in the given array.
	 */
	public static <T> void swap(T[] arr, int i, int j) {
		T tmp = arr[i];
		arr[i] = arr[j];
		arr[j] = tmp;
	}

	/**
	 * Reverses the entire array in-place.
	 */
	public static <T> void reverse(T[] arr) {
		reverse(arr, 0, arr.length - 1);
	}

	/**
	 * Reverses the elements in the range [lo, hi).
	 * Note that the upper index is exclusive.
	 */
	public static <T> void reverse(T[] arr, int lo, int hi) {
		for (int i = lo, j = hi - 1; i < j; i++, j--) {
			swap(arr, i, j);
		}
	}

	/**
	 * Mutates the array parameter to be the next permutation.
	 */
	public static <T extends Comparable<T>> boolean nextPermutation(T[] arr) {
		int p, q;
		for (p = arr.length - 2; p >= 0 && arr[p].compareTo(arr[p + 1]) >= 0; p--);
		if (p < 0) {
			return false;
		}
		for (q = arr.length - 1; arr[q].compareTo(arr[p]) <= 0; q--);
		swap(arr, p, q);
		reverse(arr, p + 1, arr.length);
		return true;
	}

	/**
	 * Computes the GCD (greatest common denominator) between two numbers.
	 */
	public static int gcd(int a, int b) {
		if (a < b)
			return gcd(b, a);

		int r = a % b;
		if (r == 0)
			return b;

		return gcd(b, r);
	}

	/**
	 * Computes the LCM (least common multiple) between all the numbers.
	 */
	public static int lcm(int... arr) {
		int lcm = arr[0];
		for (int i = 1; i < arr.length; i++)
			lcm = lcm * arr[i] / gcd(lcm, arr[i]);

		return lcm;
	}

	/**
	 * Computes the GCD (greatest common denominator) between two numbers.
	 */
	public static long gcd(long a, long b) {
		if (a < b)
			return gcd(b, a);

		long r = a % b;
		if (r == 0)
			return b;

		return gcd(b, r);
	}

	/**
	 * Computes the LCM (least common multiple) between all the numbers.
	 */
	public static long lcm(long... arr) {
		long lcm = arr[0];
		for (int i = 1; i < arr.length; i++)
			lcm = lcm * arr[i] / gcd(lcm, arr[i]);

		return lcm;
	}

	/**
	 * Computes the modular inverse, such that: ak % m = 1, for some k.
	 * See this page for details:  http://rosettacode.org/wiki/Modular_inverse
	 */
	public static long modInverse(long a, long b) {
		return BigInteger.valueOf(a).modInverse(BigInteger.valueOf(b)).longValue();
	}

	/**
	 * Computes the value of (b ^ e) % m.
	 */
	public static long modPow(long b, long e, long m) {
		return BigInteger.valueOf(b).modPow(BigInteger.valueOf(e), BigInteger.valueOf(m)).longValue();
	}

	/**
	 * Gets the specified bit (0 or 1) in the number.
	 */
	public static int getBit(int n, int b) {
		return (1 & (n >> b));
	}

	/**
	 * Sets the specified bit to one in the number.
	 */
	public static int setBit(int n, int b) {
		return ((1 << b) | n);
	}

	/**
	 * Sets the specified bit to zero in the number.
	 */
	public static int unsetBit(int n, int b) {
		return (~(1 << b) & n);
	}

	/**
	 * Sets the specified bit to the given value (0 or 1) in the number.
	 */
	public static int setBitValue(int n, int b, int v) {
		if (v == 0) {
			return unsetBit(n, b);
		} else {
			return setBit(n, b);
		}
	}

	/**
	 * Checks if the number is a power of two.
	 */
	public static boolean powerOfTwo(long n) {
		return ((n & (n - 1)) == 0);
	}

	/**
	 * Computes all the primes up to the specified number.
	 */
	public static Integer[] primesTo(int n) {
		ArrayList<Integer> primes = new ArrayList<Integer>();

		boolean[] prime = new boolean[n + 1];
		for (int i = 0; i < prime.length; i++)
			prime[i] = true;

		prime[0] = false;
		prime[1] = false;
		for (int i = 2; i < prime.length; i++) {
			if (prime[i]) {
				primes.add(i);
				if ((long) i * i <= n) {
					for (int j = i * i; j < prime.length; j += i) {
						prime[j] = false;
					}
				}
			}
		}

		return primes.toArray(new Integer[0]);
	}

	/**
	 * Returns a frequency table representing the prime factorization of the given number.
	 * For example, freq[a] = b means that primes[a] is used b times in the prime factorization.
	 */
	public static int[] primeFactors(int n, int[] primes) {
		int[] freq = new int[primes.length];
		int index = 0;
		while (index < primes.length && primes[index] <= n) {
			int prime = primes[index];
			if (n % prime == 0) {
				freq[index]++;
				n /= prime;
			} else {
				index++;
			}
		}
		return freq;
	}

	/**
	 * Computes all the combinations (i.e. Pascal's triangle) for the given levels.
	 * For example, combos[n][k] = (n choose k).
	 */
	public static long[][] computeCombinations(int n) {
		long[][] combos = new long[n + 1][n + 1];
		combos[0][0] = 1;
		for (int r = 1; r <= n; r++) {
			combos[r][0] = 1;
			for (int c = 1; c <= r; c++) {
				combos[r][c] = combos[r - 1][c - 1] + combos[r - 1][c];
			}
		}
		return combos;
	}

	/**
	 * Computes all the combinations (i.e. Pascal's triangle) for the given levels, modulo some number.
	 * For example, combos[n][k] = (n choose k) % mod.
	 */
	public static long[][] computeCombinations(int n, int mod) {
		long[][] combos = new long[n + 1][n + 1];
		combos[0][0] = 1;
		for (int r = 1; r <= n; r++) {
			combos[r][0] = 1;
			for (int c = 1; c <= r; c++) {
				combos[r][c] = (combos[r - 1][c - 1] + combos[r - 1][c]) % mod;
			}
		}
		return combos;
	}

	/**
	 * Returns true if there is a path from the source 's' to the sink 't' in the residual graph.
	 * Also fills parent[] to store the path.
	 * See here for more info:  http://www.geeksforgeeks.org/ford-fulkerson-algorithm-for-maximum-flow-problem/
	 */
	public static boolean fordFulkersonHelper(int[][] resid, int s, int t, int[] parent) {
		int V = resid.length;
		boolean[] visited = new boolean[V];
		LinkedList<Integer> q = new LinkedList<Integer>();
		q.push(s);
		visited[s] = true;
		parent[s] = -1;

		while (!q.isEmpty()) {
			int u = q.pop();
			for (int v = 0; v < V; v++) {
				if (!visited[v] && resid[u][v] > 0) {
					q.push(v);
					parent[v] = u;
					visited[v] = true;
				}
			}
		}

		return visited[t];
	}

	/**
	 * Returns the maximum flow from 's' to 't' in the given graph.
	 * See here for more info:  http://www.geeksforgeeks.org/ford-fulkerson-algorithm-for-maximum-flow-problem/
	 */
	public static int fordFulkerson(int[][] graph, int s, int t) {
		int V = graph.length;
		int[][] resid = new int[V][V];
		int[] parent = new int[V];
		int maxFlow = 0;

		for (int u = 0; u < V; u++) {
			for (int v = 0; v < V; v++) {
				resid[u][v] = graph[u][v];
			}
		}

		while (fordFulkersonHelper(resid, s, t, parent)) {
			int pathFlow = Integer.MAX_VALUE;
			for (int v = t; v != s; v = parent[v]) {
				int u = parent[v];
				pathFlow = Math.min(pathFlow, resid[u][v]);
			}
			for (int v = t; v != s; v = parent[v]) {
				int u = parent[v];
				resid[u][v] -= pathFlow;
				resid[v][u] += pathFlow;
			}
			maxFlow += pathFlow;
		}

		return maxFlow;
	}

	/**
	 * Returns true if a matching for vertex 'u' is possible.
	 * See here for more info:  http://www.geeksforgeeks.org/maximum-bipartite-matching/
	 */
	public static boolean bipartiteMatchingHelper(boolean[][] bpGraph, int u, boolean[] seen, int[] matchR) {
		int N = bpGraph[0].length;
		for (int v = 0; v < N; v++) {
			if (bpGraph[u][v] && !seen[v]) {
				seen[v] = true;
				if (matchR[v] < 0 || bipartiteMatchingHelper(bpGraph, matchR[v], seen, matchR)) {
					matchR[v] = u;
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * Returns the maximum bipartite matching from an an adjacency matrix.
	 * Note:  bpGraph[i][j] = true if there is an edge from i to j.
	 * Note:  matchIJ (array of length M) is an output variable containing the matchings, such that matchIJ[i] = j means that there is a match from i to j.
	 * Note:  matchJI (array of length N) is an output variable containing the matchings, such that matchJI[j] = i means that there is a match from i to j.
	 * See here for more info:  http://www.geeksforgeeks.org/maximum-bipartite-matching/
	 */
	public static int bipartiteMatching(boolean[][] bpGraph, int[] matchIJ, int[] matchJI) {
		int ans = bipartiteMatching(bpGraph, matchJI);

		for (int i = 0; i < matchJI.length; i++) {
			matchIJ[i] = -1;
		}

		for (int j = 0; j < matchJI.length; j++) {
			int i = matchJI[j];
			if (i >= 0) {
				matchIJ[i] = j;
			}
		}

		return ans;
	}

	/**
	 * Returns the maximum bipartite matching from an an adjacency matrix.
	 * Note:  bpGraph[i][j] = true if there is an edge from i to j.
	 * Note:  matchJI (array of length N) is an output variable containing the matchings, such that matchJI[j] = i means that there is a match from i to j.
	 * See here for more info:  http://www.geeksforgeeks.org/maximum-bipartite-matching/
	 */
	public static int bipartiteMatching(boolean[][] bpGraph, int[] matchJI) {
		int M = bpGraph.length;
		int N = bpGraph[0].length;

		for (int i = 0; i < N; i++) {
			matchJI[i] = -1;
		}

		int ans = 0;
		for (int u = 0; u < M; u++) {
			boolean[] seen = new boolean[N];
			if (bipartiteMatchingHelper(bpGraph, u, seen, matchJI)) {
				ans++;
			}
		}

		return ans;
	}

	/**
	 * Returns the maximum bipartite matching from an an adjacency matrix.
	 * Overload of the bipartiteMatching function without output parameters.
	 * See here for more info:  http://www.geeksforgeeks.org/maximum-bipartite-matching/
	 */
	public static int bipartiteMatching(boolean[][] bpGraph) {
		int N = bpGraph[0].length;
		int[] matchJI = new int[N];
		return bipartiteMatching(bpGraph, matchJI);
	}

	/**
	 * Overload of the bipartiteMatching function taking an adjacency matrix of int[][] instead of boolean[][].
	 */
	public static int bipartiteMatching(int[][] intGraph) {
		boolean[][] bpGraph = intToBooleanAdjMat(intGraph);
		return bipartiteMatching(bpGraph);
	}

	/**
	 * Overload of the bipartiteMatching function taking an adjacency matrix of int[][] instead of boolean[][].
	 * Note:  matchJI (array of length N) is an output variable containing the matchings, such that matchJI[j] = i means that there is a match from i to j.
	 */
	public static int bipartiteMatching(int[][] intGraph, int[] matchJI) {
		boolean[][] bpGraph = intToBooleanAdjMat(intGraph);
		return bipartiteMatching(bpGraph, matchJI);
	}

	/**
	 * Overload of the bipartiteMatching function taking an adjacency matrix of int[][] instead of boolean[][].
	 * Note:  matchIJ (array of length M) is an output variable containing the matchings, such that matchIJ[i] = j means that there is a match from i to j.
	 * Note:  matchJI (array of length N) is an output variable containing the matchings, such that matchJI[j] = i means that there is a match from i to j.
	 */
	public static int bipartiteMatching(int[][] intGraph, int[] matchIJ, int[] matchJI) {
		boolean[][] bpGraph = intToBooleanAdjMat(intGraph);
		return bipartiteMatching(bpGraph, matchIJ, matchJI);
	}

	/**
	 * Implementation of Hopcroft-Karp algorithm for finding maximum bipartite matching in O(sqrt(V) * E) time.
	 * See this page for implementation details:  http://www.sanfoundry.com/java-program-hopcroft-karp-algorithm/
	 */
	public static class HopcroftKarp {
		private static final int NIL = 0;
		private static final int INF = Integer.MAX_VALUE;

		private IntList[] adj;
		private int cx, cy;

		public HopcroftKarp(int xCnt, int yCnt) {
			this.cx = xCnt;
			this.cy = yCnt;
			this.adj = new IntList[cx + cy + 1];
			for (int i = 0; i < adj.length; ++i) {
				this.adj[i] = new IntList();
			}
		}

		public void addEdge(int u, int v) {
			this.adj[u].add(cx + v);
			this.adj[cx + v].add(u);
		}

		public int maxBipartiteMatching() {
			int[] pair = new int[cx + cy + 1];
			int[] dist = new int[cx + cy + 1];
			int matching = 0;
			while (BFS(pair, dist)) {
				for (int v = 1; v <= cx; ++v) {
					if (pair[v] == NIL) {
						if (DFS(v, pair, dist)) {
							matching = matching + 1;
						}
					}
				}
			}
			return matching;
		}

		private boolean BFS(int[] pair, int[] dist) {
			LinkedList<Integer> queue = new LinkedList<Integer>();
			for (int v = 1; v <= cx; ++v) {
				if (pair[v] == NIL) {
					dist[v] = 0;
					queue.add(v);
				} else {
					dist[v] = INF;
				}
			}
			dist[NIL] = INF;
			while (!queue.isEmpty()) {
				int v = queue.poll();
				if (dist[v] < dist[NIL]) {
					for (int u : adj[v]) {
						if (dist[pair[u]] == INF) {
							dist[pair[u]] = dist[v] + 1;
							queue.add(pair[u]);
						}
					}
				}
			}
			return dist[NIL] != INF;
		}

		private boolean DFS(int v, int[] pair, int[] dist) {
			if (v != NIL) {
				for (int u : adj[v]) {
					if (dist[pair[u]] == dist[v] + 1) {
						if (DFS(pair[u], pair, dist)) {
							pair[u] = v;
							pair[v] = u;
							return true;
						}
					}
				}
				dist[v] = INF;
				return false;
			}
			return true;
		}

		public static class IntList extends ArrayList<Integer> {
			private static final long serialVersionUID = 1542278293693820951L;
		}
	}

	/**
	 * Converts an integer adjacency matrix of 1's and 0's to a boolean adjacency matrix.
	 * Useful with bipartiteMatching, which takes adjancency matrix of boolean[][] as input (instead of int[][]).
	 */
	public static boolean[][] intToBooleanAdjMat(int[][] mat) {
		int M = mat.length;
		int N = mat[0].length;
		boolean[][] bMat = new boolean[M][N];
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < N; j++) {
				bMat[i][j] = (mat[i][j] != 0);
			}
		}
		return bMat;
	}

	/**
	 * Here is a Segment Tree implementation.  It currently does sums of longs.
	 * See the comments to find where to change the data-type and function of this Segment Tree.
	 */
	public static class SegmentTree {
		public SegmentTreeNode[] leaves;
		public SegmentTreeNode root;

		public SegmentTree(int n) {
			this.leaves = new SegmentTreeNode[n];
			this.root = new SegmentTreeNode(this, null, 0, n - 1);
		}

		// modify the data-type of this segment tree
		public SegmentTree(long[] vals) {
			this(vals.length);
			for (int i = 0; i < vals.length; i++) {
				this.insert(i, vals[i]);
			}
		}

		// modify the data-type of this segment tree
		public void insert(int idx, long v) {
			this.leaves[idx].setAndUpdate(v);
		}

		// modify the data-type of this segment tree
		public long get(int idx) {
			return this.leaves[idx].val;
		}

		// modify the data-type of this segment tree
		public long get(int lower, int upper) {
			return this.root.getRange(lower, upper);
		}

		private static class SegmentTreeNode {
			public int L;
			public int R;

			// modify the data-type of this segment tree
			public long val;

			public SegmentTree tree;
			public SegmentTreeNode parent;
			public SegmentTreeNode left;
			public SegmentTreeNode rite;

			public SegmentTreeNode(SegmentTree t, SegmentTreeNode p, int lower, int upper) {
				this.tree = t;
				this.parent = p;
				this.L = lower;
				this.R = upper;

				if (lower == upper) {
					this.tree.leaves[lower] = this;
				} else {
					int mid = (lower + upper) / 2;
					this.left = new SegmentTreeNode(tree, this, lower, mid);
					this.rite = new SegmentTreeNode(tree, this, mid + 1, upper);
				}
			}

			// modify the data-type of this segment tree
			public void setAndUpdate(long v) {
				this.val = v;
				this.update();
			}

			// modify the function (i.e. max, min, sum, etc.) in this method
			public void update() {
				if (this.left != null && this.rite != null) {
					this.val = this.left.val + this.rite.val;
				} else if (this.left != null) {
					this.val = this.left.val;
				} else if (this.rite != null) {
					this.val = this.rite.val;
				}

				if (this.parent != null) {
					this.parent.update();
				}
			}

			// modify the function & the default return value in this method
			public long getRange(int lower, int upper) {
				if (this.L >= lower && this.R <= upper) {
					return this.val;
				} else if (this.L > upper || this.R < lower) {
					// modify the default value here (if it's not found)
					return 0;
				} else {
					// modify the function that this Segment Tree uses
					return this.left.getRange(lower, upper) + this.rite.getRange(lower, upper);
				}
			}
		}
	}

	/**
	 * Here is a generic Segment Tree implementation.
	 * It requires a Combiner and DefaultProvider for the custom datatype.
	 */
	public static class GenericSegmentTree<T> {
		public ArrayList<SegmentTreeNode> leaves;
		public SegmentTreeNode root;
		public Combiner<T> combiner;
		public DefaultProvider<T> defaultProvider;

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

		public T get(int lower, int upper) {
			return this.root.getRange(lower, upper);
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
					this.val = combiner.combine(this.left.val, this.rite.val);
				} else if (this.left != null) {
					this.val = this.left.val;
				} else if (this.rite != null) {
					this.val = this.rite.val;
				}

				if (this.parent != null) {
					this.parent.update();
				}
			}

			public T getRange(int lower, int upper) {
				if (this.L >= lower && this.R <= upper) {
					return this.val;
				} else if (this.L > upper || this.R < lower) {
					// access outer class GenericSegmentTree
					return defaultProvider.getDefault();
				} else {
					// access outer class GenericSegmentTree
					return combiner.combine(this.left.getRange(lower, upper), this.rite.getRange(lower, upper));
				}
			}
		}
	}

	/**
	 * Here is a Segment Tree implementation.  It currently does sums of longs.
	 * It contains multiple values in the leaves.
	 */
	public static class LeafArraySegmentTree {
		private long L;
		private long R;

		// modify the data-type of this segment tree
		private long val;
		private long[] leaf;
		private int width;

		private LeafArraySegmentTree parent;
		private LeafArraySegmentTree left;
		private LeafArraySegmentTree rite;

		public LeafArraySegmentTree(long lo, long hi, int w) {
			this(lo, hi, w, null);
		}

		private LeafArraySegmentTree(long lo, long hi, int w, LeafArraySegmentTree p) {
			this.L = lo;
			this.R = hi;
			this.width = w;
			this.parent = p;
			if (hi - lo + 1 <= width) {
				int size = (int) (hi - lo + 1);
				this.leaf = new long[size];
			}
		}

		public LeafArraySegmentTree getLeaf(long k) {
			if (leaf != null) {
				return this;
			}
			long M = (L + R) >> 1;
			if (L <= k && k <= M) {
				if (left == null) {
					left = new LeafArraySegmentTree(L, M, width, this);
				}
				return left.getLeaf(k);
			} else {
				if (rite == null) {
					rite = new LeafArraySegmentTree(M + 1, R, width, this);
				}
				return rite.getLeaf(k);
			}
		}

		public void increment(long k, long v) {
			LeafArraySegmentTree ast = getLeaf(k);
			int offset = (int) (k - ast.L);
			ast.leaf[offset] += v;
			ast.val += v;
			ast = ast.parent;
			while (ast != null) {
				ast.val = valueOf(ast.left) + valueOf(ast.rite);
				ast = ast.parent;
			}
		}

		public long get(long k) {
			return get(k, k);
		}

		public long get(long lo, long hi) {
			if (L > hi || R < lo) {
				return defaultValue();
			}
			if (L >= lo && R <= hi) {
				return val;
			}
			if (leaf != null) {
				long ans = 0;
				for (int i = 0; i < leaf.length; ++i) {
					if (lo <= L + i && L + i <= hi) {
						ans += leaf[i];
					}
				}
				return ans;
			}
			return tryGet(left, lo, hi) + tryGet(rite, lo, hi);
		}

		private static long defaultValue() {
			return 0;
		}

		private static long valueOf(LeafArraySegmentTree ast) {
			if (ast == null) {
				return defaultValue();
			}
			return ast.val;
		}

		private static long tryGet(LeafArraySegmentTree ast, long lo, long hi) {
			if (ast == null) {
				return defaultValue();
			}
			return ast.get(lo, hi);
		}
	}

	/**
	 * Here is a Segment Tree implementation.  It currently does sums of longs.
	 * It lazily initializes the child nodes as they are needed.
	 */
	public static class LazySegmentTree {
		private long L;
		private long R;

		// modify the data-type of this segment tree
		private long val;

		private LazySegmentTree parent;
		private LazySegmentTree left;
		private LazySegmentTree rite;

		public LazySegmentTree(long lo, long hi) {
			this(lo, hi, null);
		}

		private LazySegmentTree(long lo, long hi, LazySegmentTree p) {
			this.L = lo;
			this.R = hi;
			this.parent = p;
		}

		public LazySegmentTree getLeaf(long k) {
			if (L == R) {
				return this;
			}
			long M = (L + R) >> 1;
			if (L <= k && k <= M) {
				if (left == null) {
					left = new LazySegmentTree(L, M, this);
				}
				return left.getLeaf(k);
			} else {
				if (rite == null) {
					rite = new LazySegmentTree(M + 1, R, this);
				}
				return rite.getLeaf(k);
			}
		}

		public void set(long k, long v) {
			LazySegmentTree ast = getLeaf(k);
			ast.val = v;
			ast = ast.parent;
			while (ast != null) {
				ast.val = valueOf(ast.left) + valueOf(ast.rite);
				ast = ast.parent;
			}
		}

		public void increment(long k, long v) {
			LazySegmentTree ast = getLeaf(k);
			ast.val += v;
			ast = ast.parent;
			while (ast != null) {
				ast.val = valueOf(ast.left) + valueOf(ast.rite);
				ast = ast.parent;
			}
		}

		public long get(long k) {
			return get(k, k);
		}

		public long get(long lo, long hi) {
			if (L > hi || R < lo) {
				return defaultValue();
			}
			if (L >= lo && R <= hi) {
				return val;
			}
			return tryGet(left, lo, hi) + tryGet(rite, lo, hi);
		}

		private static long defaultValue() {
			return 0;
		}

		private static long valueOf(LazySegmentTree ast) {
			if (ast == null) {
				return defaultValue();
			}
			return ast.val;
		}

		private static long tryGet(LazySegmentTree ast, long lo, long hi) {
			if (ast == null) {
				return defaultValue();
			}
			return ast.get(lo, hi);
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
				return getSumLTE(hi) - getSumLTE(lo - 1);
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
	 * Here is a Segment Tree implementation for range sums.
	 * It uses lazy propagation for greater efficiency.
	 */
	public static class SegmentTreeRangeSum {
		public SegmentTreeNodeRS[] leaves;
		public SegmentTreeNodeRS root;

		public SegmentTreeRangeSum(int n) {
			this.leaves = new SegmentTreeNodeRS[n];
			this.root = new SegmentTreeNodeRS(this, null, 0, n - 1);
		}

		public void addRange(int lower, int upper, long v) {
			this.root.addRange(lower, upper, v);
		}

		public long getRange(int lower, int upper) {
			return this.root.getRange(lower, upper);
		}

		private static class SegmentTreeNodeRS {
			public int L;
			public int R;

			public long val;

			public boolean update;
			public long change;

			public SegmentTreeRangeSum tree;
			public SegmentTreeNodeRS left;
			public SegmentTreeNodeRS rite;

			public SegmentTreeNodeRS(SegmentTreeRangeSum t, SegmentTreeNodeRS p, int lower, int upper) {
				this.tree = t;
				this.L = lower;
				this.R = upper;

				if (lower == upper) {
					this.tree.leaves[lower] = this;
				} else {
					int mid = (lower + upper) / 2;
					this.left = new SegmentTreeNodeRS(tree, this, lower, mid);
					this.rite = new SegmentTreeNodeRS(tree, this, mid + 1, upper);
				}
			}

			public void addRange(int lower, int upper, long v) {
				if (lower > upper) {
					return;
				}

				this.val += v * (upper - lower + 1);

				if (lower == this.L && upper == this.R) {
					trySetUpdate(this.left, v);
					trySetUpdate(this.rite, v);
					return;
				}

				tryAddRange(this.left, lower, upper, v);
				tryAddRange(this.rite, lower, upper, v);
			}

			public long getRange(int lower, int upper) {
				if (lower > upper) {
					return 0;
				}

				if (this.update) {
					this.val += this.change * (this.R - this.L + 1);
					trySetUpdate(this.left, this.change);
					trySetUpdate(this.rite, this.change);
					this.update = false;
					this.change = 0;
				}

				if (this.L == lower && this.R == upper) {
					return this.val;
				}

				long sum = 0;
				if (this.left != null) {
					int endL = Math.min(this.left.R, upper);
					sum += this.left.getRange(lower, endL);
				}
				if (this.rite != null) {
					int startR = Math.max(this.rite.L, lower);
					sum += this.rite.getRange(startR, upper);
				}
				return sum;
			}

			public static void tryAddRange(SegmentTreeNodeRS node, int lower, int upper, long v) {
				if (node != null) {
					int start = Math.max(node.L, lower);
					int end = Math.min(node.R, upper);
					node.addRange(start, end, v);
				}
			}

			private static void trySetUpdate(SegmentTreeNodeRS node, long v) {
				if (node != null) {
					node.update = true;
					node.change += v;
				}
			}
		}
	}

	/**
	 * Implementation of Segment Tree based on perfect binary tree indexed on an array.
	 * NOTE:  This might not actually work, so use at your own risk!
	 */
	public static class SegmentTreeArray {
		public int[] left;
		public int[] rite;
		public long[] values;

		private final int LAYERS;

		public SegmentTreeArray(int layers) {
			int len = (1 << (layers + 1)) - 1;
			left = new int[len];
			rite = new int[len];
			values = new long[len];
			LAYERS = layers;
			setLayer(0, 0, (1 << layers) - 1);
		}

		private void setLayer(int idx, int lower, int upper) {
			this.left[idx] = lower;
			this.rite[idx] = upper;
			if (lower != upper) {
				int mid = (lower + upper) / 2;
				setLayer(getLeftChild(idx), lower, mid);
				setLayer(getRiteChild(idx), mid + 1, upper);
			}
		}

		public void insert(int idx, long val) {
			int pos = getActualIdx(this.LAYERS, idx);
			this.values[pos] = val;
			if (idx > 0) {
				this.update(getParentIdx(pos));
			}
		}

		private void update(int idx) {
			int idxL = getLeftChild(idx);
			int idxR = getRiteChild(idx);
			this.values[idx] = this.values[idxL] + this.values[idxR];
			if (idx > 0) {
				this.update(getParentIdx(idx));
			}
		}

		public long get(int idx) {
			int pos = getActualIdx(this.LAYERS, idx);
			return this.values[pos];
		}

		public long getRange(int lower, int upper) {
			return getRangeHelper(0, lower, upper);
		}

		private long getRangeHelper(int idx, int lower, int upper) {
			if (this.left[idx] >= lower && this.rite[idx] <= upper) {
				return this.values[idx];
			} else if (this.left[idx] > upper || this.rite[idx] < lower) {
				return 0;
			} else {
				return getRangeHelper(getLeftChild(idx), lower, upper) + getRangeHelper(getRiteChild(idx), lower, upper);
			}
		}

		public static int getActualIdx(int layer, int idx) {
			return (1 << layer) - 1 + idx;
		}

		private static int getParentIdx(int idx) {
			return (idx - 1) / 2;
		}

		private static int getLeftChild(int idx) {
			return 2 * idx + 1;
		}

		private static int getRiteChild(int idx) {
			return 2 * idx + 2;
		}
	}

	/**
	 * Data structure for calculating sub-array sums in O(1) time.
	 * This inserts values in O(blockCount + blockSize) time.
	 * Usually, for an array of size N, set blockCount = blockSize = sqrt(N) for maximum efficiency.
	 */
	public static class BlockArray {
		private final int blockCount;
		private final int blockSize;

		long[] values;
		long[] blockSums;
		long[][] subBlocks;

		public BlockArray(int sz, int cnt) {
			this.values = new long[cnt * sz];
			this.blockSums = new long[cnt];
			this.subBlocks = new long[cnt][sz];

			this.blockSize = sz;
			this.blockCount = cnt;
		}

		public void add(int idx, long val) {
			this.values[idx] += val;
			int blockNum = idx / this.blockSize;
			int blockIdx = idx % this.blockSize;
			for (int i = blockIdx; i < this.blockSize; i++) {
				this.subBlocks[blockNum][i] += val;
			}
			for (int i = blockNum; i < blockCount; i++) {
				this.blockSums[i] += val;
			}
		}

		public long get(int idx) {
			return this.values[idx];
		}

		public long getSum(int idx) {
			if (idx < 0) {
				return 0;
			}
			int blockNum = idx / this.blockSize;
			int blockIdx = idx % this.blockSize;
			long blockPart = (blockNum > 0) ? this.blockSums[blockNum - 1] : 0;
			long subPart = this.subBlocks[blockNum][blockIdx];
			return blockPart + subPart;
		}

		public long getSum(int lower, int upper) {
			return this.getSum(upper) - this.getSum(lower - 1);
		}
	}

	/**
	 * Node used for finding Eulerian paths.
	 */
	public static class EulerNode {
		public String val;
		public LinkedList<EulerNode> nexts;
		public int inDeg;

		public EulerNode(String v) {
			this.val = v;
			this.nexts = new LinkedList<EulerNode>();
			this.inDeg = 0;
		}
	}

	/**
	 * Finds an Eulerian path if one exists.
	 * Otherwise returns null.
	 * This will destructively modify the nodes!
	 */
	public static EulerNode[] findEulerianPath(EulerNode[] nodes) {
		EulerNode start = null;
		EulerNode end = null;
		for (EulerNode n : nodes) {
			int outDeg = n.nexts.size();
			if (n.inDeg + 1 == outDeg) {
				if (start == null) {
					start = n;
				} else {
					return null;
				}
			} else if (n.inDeg - 1 == outDeg) {
				if (end == null) {
					end = n;
				} else {
					return null;
				}
			} else if (Math.abs(n.inDeg - outDeg) > 1) {
				return null;
			}
		}

		if ((start == null) ^ (end == null)) {
			return null;
		}

		if (start != null) {
			end.nexts.addLast(start);
			start.inDeg++;
		} else {
			for (int i = 0; i < nodes.length; i++) {
				if (!nodes[i].nexts.isEmpty() && nodes[i].inDeg > 0) {
					start = nodes[i];
					end = nodes[i];
					break;
				}
			}
		}

		LinkedList<EulerNode> procQ = new LinkedList<EulerNode>();
		ArrayList<EulerNode> done = new ArrayList<EulerNode>();

		procQ.add(start);
		while (!procQ.isEmpty()) {
			EulerNode curr = procQ.removeFirst();
			if (curr.nexts.isEmpty()) {
				done.add(curr);
				continue;
			}
			ArrayList<EulerNode> cycle = getAnyCycleFrom(curr);
			for (int i = cycle.size() - 1; i >= 0; i--) {
				procQ.addFirst(cycle.get(i));
			}
		}
		if (start != end) {
			done.remove(done.size() - 1);
		}

		for (EulerNode n : nodes) {
			if (!n.nexts.isEmpty() || n.inDeg != 0) {
				return null;
			}
		}

		int startIdx = 0;
		int sz = done.size();
		for (int i = 1; i < sz; i++) {
			if (done.get(i) == start && done.get(i - 1) == end) {
				startIdx = i;
				break;
			}
		}

		EulerNode[] arr = new EulerNode[sz];
		for (int i = 0; i < sz; i++) {
			arr[i] = done.get((i + startIdx) % sz);
		}
		return arr;
	}

	/**
	 * Gets any cycle starting with the given node.
	 * This will destructively modify the nodes.
	 */
	public static ArrayList<EulerNode> getAnyCycleFrom(EulerNode start) {
		ArrayList<EulerNode> lst = new ArrayList<EulerNode>();
		lst.add(start);
		EulerNode curr = start;
		do {
			curr = curr.nexts.removeFirst();
			curr.inDeg--;
			lst.add(curr);
		} while (curr != start);
		return lst;
	}

	/**
	 * Creates a sparse table for looking up the minimum between two indices.
	 * Data should not be modified once the table is created.
	 * It takes O(n log(n)) time to create the table.
	 */
	public static class SparseTable {
		private int[][] table;

		public SparseTable(int[] a) {
			int n = a.length;
			int p = Integer.SIZE - Integer.numberOfLeadingZeros(n - 1);

			this.table = new int[n][p];

			for (int i = 0; i < n; i++) {
				this.table[i][0] = a[i];
			}
			for (int j = 1; j < p; j++) {
				for (int i = 0; i < n; i++) {
					int off = (1 << (j - 1));
					int next = (i + off >= n) ? Integer.MAX_VALUE : this.table[i + off][j - 1];
					this.table[i][j] = Math.min(this.table[i][j - 1], next);
				}
			}
		}

		public int getMin(int lower, int upper) {
			int k = (int) (Math.log(upper - lower + 1) / Math.log(2));
			int off = (1 << k);
			return Math.min(this.table[lower][k], this.table[upper - off + 1][k]);
		}
	}

	/**
	 * Gets the component of the graph that the input node is part of.
	 */
	public static HashSet<InOutNode> getComponent(InOutNode seed) {
		HashSet<InOutNode> visited = new HashSet<InOutNode>();
		LinkedList<InOutNode> q = new LinkedList<InOutNode>();
		q.addLast(seed);
		while (!q.isEmpty()) {
			InOutNode curr = q.removeFirst();
			if (visited.contains(curr)) {
				continue;
			}
			visited.add(curr);
			for (InOutNode n : curr.nexts) {
				q.addLast(n);
			}
			for (InOutNode n : curr.prevs) {
				q.addLast(n);
			}
		}
		return visited;
	}

	/**
	 * This method returns a list of the nodes after they have been topologically sorted.
	 * It will return null if it was impossible to do so.
	 * NOTE: This will destructively modify the list!
	 */
	public static ArrayList<InOutNode> toplogicalSort(HashSet<InOutNode> nodes) {
		ArrayList<InOutNode> lst = new ArrayList<InOutNode>();
		HashSet<InOutNode> starts = new HashSet<InOutNode>();
		LinkedList<InOutNode> q = new LinkedList<InOutNode>();
		for (InOutNode n : nodes) {
			if (n.prevs.size() == 0) {
				starts.add(n);
				q.addLast(n);
			}
		}
		while (!q.isEmpty()) {
			InOutNode n = q.removeFirst();
			if (!starts.contains(n)) {
				continue;
			}
			lst.add(n);
			for (InOutNode m : n.nexts) {
				m.prevs.remove(n);
				if (m.prevs.size() == 0) {
					starts.add(m);
					q.addLast(m);
				}
			}
			n.nexts.clear();
		}
		for (InOutNode n : nodes) {
			if ((n.prevs.size() > 0) || (n.nexts.size() > 0)) {
				return null;
			}
		}
		return lst;
	}

	/**
	 * A Node containing information about incoming and outgoing edges.
	 */
	public static class InOutNode {
		public HashSet<InOutNode> nexts;
		public HashSet<InOutNode> prevs;

		public InOutNode() {
			this.nexts = new HashSet<InOutNode>();
			this.prevs = new HashSet<InOutNode>();
		}
	}

	/**
	 * Implementation of Manacher's algorithm to find the longest palindromic substring in O(N) time.
	 * See this page for more details:  http://en.wikipedia.org/wiki/Longest_palindromic_substring
	 */
	public static class ManachersAlgorithm {
		public static String findLongestPalindrome(String s) {
			if (s == null || s.length() == 0)
				return "";

			char[] s2 = addBoundaries(s.toCharArray());
			int[] p = findPalindromeLengths(s2);
			int len = 0;
			int c = 0;
			for (int i = 1; i < s2.length; i++) {
				if (len < p[i]) {
					len = p[i];
					c = i;
				}
			}
			char[] ss = Arrays.copyOfRange(s2, c - len, c + len + 1);
			return String.valueOf(removeBoundaries(ss));
		}

		/*
		 * Returns an array containing the longest palindromic substring length at each index.
		 * Input is a string with the boundaries already added.
		 */
		public static int[] findPalindromeLengths(char[] s2) {
			int[] p = new int[s2.length];
			int c = 0, r = 0; // Here the first element in s2 has been processed.
			int m = 0, n = 0; // The walking indices to compare if two elements are the same
			for (int i = 1; i < s2.length; i++) {
				if (i > r) {
					p[i] = 0;
					m = i - 1;
					n = i + 1;
				} else {
					int i2 = c * 2 - i;
					if (p[i2] < (r - i)) {
						p[i] = p[i2];
						m = -1; // This signals bypassing the while loop below.
					} else {
						p[i] = r - i;
						n = r + 1;
						m = i * 2 - n;
					}
				}
				while (m >= 0 && n < s2.length && s2[m] == s2[n]) {
					p[i]++;
					m--;
					n++;
				}
				if ((i + p[i]) > r) {
					c = i;
					r = i + p[i];
				}
			}
			return p;
		}

		/*
		 * Returns an array containing the longest palindromic substring length at each index.
		 * Input is the original string.
		 */
		public static int[] findPalindromeLengths(String s) {
			if (s == null || s.length() == 0)
				return null;

			char[] s2 = addBoundaries(s.toCharArray());
			int[] p = findPalindromeLengths(s2);
			return p;
		}

		private static char[] addBoundaries(char[] cs) {
			if (cs == null || cs.length == 0)
				return "||".toCharArray();

			char[] cs2 = new char[cs.length * 2 + 1];
			for (int i = 0; i < (cs2.length - 1); i = i + 2) {
				cs2[i] = '|';
				cs2[i + 1] = cs[i / 2];
			}
			cs2[cs2.length - 1] = '|';
			return cs2;
		}

		private static char[] removeBoundaries(char[] cs) {
			if (cs == null || cs.length < 3)
				return "".toCharArray();

			char[] cs2 = new char[(cs.length - 1) / 2];
			for (int i = 0; i < cs2.length; i++) {
				cs2[i] = cs[i * 2 + 1];
			}
			return cs2;
		}
	}

	public static class DeterministicMillerRabin {
		private static final BigInteger[] TWO_POW = new BigInteger[100];
		static {
			TWO_POW[0] = BigInteger.ONE;
			for (int i = 1; i < TWO_POW.length; i++) {
				TWO_POW[i] = TWO_POW[i - 1].add(TWO_POW[i - 1]);
			}
		}

		private static final BigInteger TWO_BI = BigInteger.valueOf(2);

		public static long modPow(long b, BigInteger e_bi, long m) {
			long fac = b;
			long curr = 1;
			while (e_bi.compareTo(BigInteger.ZERO) > 0) {
				if (e_bi.testBit(0)) {
					curr = (curr * fac) % m;
				}
				fac = (fac * fac) % m;
				e_bi = e_bi.divide(TWO_BI);
			}
			return curr;
		}

		public static boolean isPrime(long n) {
			if ((n % 2 == 0 && n != 2) || (n < 2) || (n % 3 == 0 && n != 3)) {
				return false;
			}
			if (n <= 3) {
				return true;
			}

			long d = n / 2;
			long s = 1;
			while (d % 2 == 0) {
				d /= 2;
				s++;
			}

			BigInteger d_bi = BigInteger.valueOf(d);
			BigInteger n_bi = BigInteger.valueOf(n);

			int[] primesToTest = getPrimesToTest(n);
			for (int a : primesToTest) {
				if (try_composite(a, d_bi, n_bi, s)) {
					return false;
				}
			}
			return true;
		}

		private static boolean try_composite(int a, BigInteger d_bi, BigInteger n_bi, long s) {
			BigInteger a_bi = BigInteger.valueOf(a);
			if (a_bi.modPow(d_bi, n_bi).equals(BigInteger.ONE)) {
				return false;
			}
			for (int i = 0; i < s; i++) {
				if (a_bi.modPow(TWO_POW[i].multiply(d_bi), n_bi).equals(n_bi.subtract(BigInteger.ONE))) {
					return false;
				}
			}
			return true;
		}

		private static int[] getPrimesToTest(long n) {
			if (n < 1373653) {
				return new int[] {2, 3};
			} else if (n < 9080191) {
				return new int[] {31, 73};
			} else if (n < 4759123141L) {
				return new int[] {2, 7, 61};
			} else if (n < 1122004669633L) {
				return new int[] {2, 13, 23, 1662803};
			} else if (n < 2152302898747L) {
				return new int[] {2, 3, 5, 7, 11};
			} else if (n < 3474749660383L) {
				return new int[] {2, 3, 5, 7, 11, 13};
			} else {
				return new int[] {2, 3, 5, 7, 11, 13, 17};
			}
		}
	}

	/**
	 * Returns the index that starts the largest lexicographical rotation of the input String.
	 */
	public static int largestLexicographicalRotation(String S) {
		int N = S.length();
		char[] arr = (S + S).toCharArray();
		int[] f = new int[2 * N];
		Arrays.fill(f, -1);
		int k = 0;
		for (int j = 1; j < 2 * N; j++) {
			int i = f[j - k - 1];
			while (i != -1 && arr[j] != arr[k + i + 1]) {
				if (arr[j] > arr[k + i + 1]) {
					k = j - i - 1;
				}
				i = f[i];
			}
			if (i == -1 && arr[j] != arr[k + i + 1]) {
				if (arr[j] > arr[k + i + 1]) {
					k = j;
				}
				f[j - k] = -1;
			} else {
				f[j - k] = i + 1;
			}
		}
		return k;
	}

	/**
	 * Node used in Tarjan's algorithm from finding strongly connected components.
	 * The field "low" is null if not yet visited, and is set to Integer.MAX_VALUE once it is part of a strong connected component.
	 */
	public static class NodeSCC {
		public NodeSetSCC neighbors;
		public Integer low;

		public NodeSCC(int p) {
			this.neighbors = new NodeSetSCC();
		}
	}

	/**
	 * HashSet of NodeSCC used in Tarjan's algorithm of finding strongly connected components.
	 */
	public static class NodeSetSCC extends HashSet<NodeSCC> {
		private static final long serialVersionUID = 5191526706334303089L;
	}

	/**
	 * Implementation of Tarjan's algorithm for finding strongly connected components.
	 * The parameter for the constructor is an array of nodes representing the graph.  These nodes will be modified!
	 * After construction, the strongly connected components will be in the public field called "scc" where it can be accessed.
	 * Run-time: O(|V| + |E|)
	 */
	public static class TarjanSCC {
		private int index;
		private LinkedList<NodeSCC> stack;

		public ArrayList<NodeSetSCC> scc;

		public TarjanSCC(NodeSCC[] nodes) {
			this.index = 0;
			this.stack = new LinkedList<NodeSCC>();
			this.scc = new ArrayList<NodeSetSCC>();

			for (int i = 0; i < nodes.length; i++) {
				if (nodes[i].low == null) {
					strongConnect(nodes[i]);
				}
			}
		}

		private void strongConnect(NodeSCC curr) {
			curr.low = index++;
			stack.push(curr);
			boolean root = true;
			for (NodeSCC next : curr.neighbors) {
				if (next.low == null) {
					strongConnect(next);
				}
				if (curr.low > next.low) {
					curr.low = next.low;
					root = false;
				}
			}
			if (root) {
				NodeSetSCC set = new NodeSetSCC();
				NodeSCC n;
				do {
					n = stack.pop();
					n.low = Integer.MAX_VALUE;
					set.add(n);
				} while (curr != n);
				scc.add(set);
			}
		}
	}

	/**
	 * Computes a power by exponentiation by squaring method.
	 * Consider using Math.pow(...) instead.
	 */
	public static double fastPow(double b, int e) {
		double v = 1;
		double p = b;
		while (e > 0) {
			if (e % 2 == 1) {
				v *= p;
			}
			p *= p;
			e /= 2;
		}
		return v;
	}

	/**
	 * Computes the number of inversions in an array.
	 * See for more details:  http://algs4.cs.princeton.edu/22mergesort/Inversions.java.html
	 */
	public static class Inversion {
		// merge and count
		private static <T extends Comparable<T>> long merge(T[] a, T[] aux, int lo, int mid, int hi) {
			long inversions = 0;

			// copy to aux[]
			for (int k = lo; k <= hi; k++) {
				aux[k] = a[k];
			}

			// merge back to a[]
			int i = lo, j = mid + 1;
			for (int k = lo; k <= hi; k++) {
				if (i > mid)
					a[k] = aux[j++];
				else if (j > hi)
					a[k] = aux[i++];
				else if (less(aux[j], aux[i])) {
					a[k] = aux[j++];
					inversions += (mid - i + 1);
				} else
					a[k] = aux[i++];
			}
			return inversions;
		}

		// return the number of inversions in the subarray b[lo..hi]
		// side effect b[lo..hi] is rearranged in ascending order
		private static <T extends Comparable<T>> long count(T[] a, T[] b, T[] aux, int lo, int hi) {
			long inversions = 0;
			if (hi <= lo)
				return 0;
			int mid = lo + (hi - lo) / 2;
			inversions += count(a, b, aux, lo, mid);
			inversions += count(a, b, aux, mid + 1, hi);
			inversions += merge(b, aux, lo, mid, hi);
			return inversions;
		}


		// count number of inversions in the array a[] - do not overwrite a[]
		public static <T extends Comparable<T>> long count(T[] a) {
			T[] b = Arrays.copyOf(a, a.length);
			T[] aux = Arrays.copyOf(a, a.length);
			for (int i = 0; i < a.length; i++)
				b[i] = a[i];
			long inversions = count(a, b, aux, 0, a.length - 1);
			return inversions;
		}

		// is v < w ?
		private static <T extends Comparable<T>> boolean less(T v, T w) {
			return (v.compareTo(w) < 0);
		}
	}

	/**
	 * Implementation of Suffix Array (which also contains Least Common Prefix array)
	 * See the editorial for implmentation details:  https://www.hackerearth.com/may-hem-15/algorithm/rotations-and-inversions/
	 */
	public static class SuffixArray {
		public int[] index, LCP;

		public SuffixArray(int[] arr) {
			this.index = calcSuffixArray(arr);
			this.LCP = calcLCP(this.index, arr);
		}

		private static int[] calcSuffixArray(final int[] str) {
			int n = str.length;
			Integer[] order = new Integer[n];
			for (int i = 0; i < n; i++)
				order[i] = n - 1 - i;
			Arrays.sort(order, new Comparator<Integer>() {
				public int compare(Integer o1, Integer o2) {
					return str[o1] - str[o2];
				}
			});
			// sa[i] - suffix on i'th position after sorting by first len characters
			// rank[i] - position of the i'th suffix after sorting by first len characters
			int[] sa = new int[n];
			int[] rank = new int[n];
			for (int i = 0; i < n; i++) {
				sa[i] = order[i];
				rank[i] = str[i];
			}
			for (int len = 1; len < n; len *= 2) {
				int[] r = rank.clone();
				for (int i = 0; i < n; i++) {
					// condition s1 + len < n simulates 0-symbol at the end of the string
					// a separate class is created for each suffix followed by 0-symbol
					rank[sa[i]] = i > 0 && r[sa[i - 1]] == r[sa[i]] && sa[i - 1] + len < n && r[sa[i - 1] + len / 2] == r[sa[i] + len / 2] ? rank[sa[i - 1]] : i;
				}
				// Suffixes are already sorted by first len characters
				// Now sort suffixes by first len * 2 characters
				int[] cnt = new int[n];
				for (int i = 0; i < n; i++)
					cnt[i] = i;
				int[] s = sa.clone();
				for (int i = 0; i < n; i++) {
					// s[i] - order of suffixes sorted by first len characters
					// (s[i] - len) - order of suffixes sorted only by second len characters
					int s1 = s[i] - len;
					// sort only suffixes of length > len, others are already sorted
					if (s1 >= 0)
						sa[cnt[rank[s1]]++] = s1;
				}
			}
			return sa;
		}

		// longest common prefixes array in O(n)
		private static int[] calcLCP(int[] sa, int[] s) {
			int n = sa.length;
			int[] rank = new int[n];
			for (int i = 0; i < n; i++)
				rank[sa[i]] = i;
			int[] lcp = new int[n];
			for (int i = 0, h = 0; i < n; i++) {
				if (rank[i] < n - 1) {
					int j = sa[rank[i] + 1];
					while (Math.max(i, j) + h < s.length && s[i + h] == s[j + h]) {
						++h;
					}
					lcp[rank[i]] = h;
					if (h > 0)
						--h;
				}
			}
			return lcp;
		}
	}

	/**
	 * Computes the prefix array for the KMP algorithm.
	 */
	public static int[] getPrefixKMP(String W) {
		int n = W.length();
		int[] T = new int[n];
		T[0] = -1;
		T[1] = 0;
		int pos = 2;
		int cnd = 0;
		while (pos < n) {
			if (W.charAt(pos - 1) == W.charAt(cnd)) {
				cnd++;
				T[pos] = cnd;
				pos++;
			} else if (cnd > 0) {
				cnd = T[cnd];
			} else {
				T[pos] = 0;
				pos++;
			}
		}
		return T;
	}

	/**
	 * Finds the indices of all occurrences of W in S.
	 * This includes all overlapping substrings.
	 */
	public static Integer[] getMatchesKMP(String S, String W) {
		int m = 0;
		int i = 0;
		int[] T = getPrefixKMP(S);
		ArrayList<Integer> lst = new ArrayList<Integer>();
		while (m + i < S.length()) {
			if (W.charAt(i) == S.charAt(m + i)) {
				if (i == W.length() - 1) {
					lst.add(m);
					i = 0;
					m++;
				} else {
					i++;
				}
			} else {
				if (T[i] > -1) {
					m += i - T[i];
					i = T[i];
				} else {
					i = 0;
					m++;
				}
			}
		}
		return lst.toArray(new Integer[0]);
	}

	// Simplex solver for linear programming
	// See this page for details:  http://algs4.cs.princeton.edu/65reductions/Simplex.java.html
	public static class Simplex {
		private static final double EPSILON = 1.0E-10;
		private double[][] a; // tableaux
		private int M; // number of constraints
		private int N; // number of original variables

		private int[] basis; // basis[i] = basic variable corresponding to row i
		// only needed to print out solution, not book

		// sets up the simplex tableaux
		public Simplex(double[][] A, double[] b, double[] c) {
			M = b.length;
			N = c.length;
			a = new double[M + 1][N + M + 1];
			for (int i = 0; i < M; i++)
				for (int j = 0; j < N; j++)
					a[i][j] = A[i][j];
			for (int i = 0; i < M; i++)
				a[i][N + i] = 1.0;
			for (int j = 0; j < N; j++)
				a[M][j] = c[j];
			for (int i = 0; i < M; i++)
				a[i][M + N] = b[i];

			basis = new int[M];
			for (int i = 0; i < M; i++)
				basis[i] = N + i;

			solve();

			// check optimality conditions
			assert check(A, b, c);
		}

		// run simplex algorithm starting from initial BFS
		private void solve() {
			while (true) {

				// find entering column q
				int q = bland();
				if (q == -1)
					break; // optimal

				// find leaving row p
				int p = minRatioRule(q);
				if (p == -1)
					throw new ArithmeticException("Linear program is unbounded");

				// pivot
				pivot(p, q);

				// update basis
				basis[p] = q;
			}
		}

		// lowest index of a non-basic column with a positive cost
		private int bland() {
			for (int j = 0; j < M + N; j++)
				if (a[M][j] > 0)
					return j;
			return -1; // optimal
		}

		// find row p using min ratio rule (-1 if no such row)
		private int minRatioRule(int q) {
			int p = -1;
			for (int i = 0; i < M; i++) {
				if (a[i][q] <= 0)
					continue;
				else if (p == -1)
					p = i;
				else if ((a[i][M + N] / a[i][q]) < (a[p][M + N] / a[p][q]))
					p = i;
			}
			return p;
		}

		// pivot on entry (p, q) using Gauss-Jordan elimination
		private void pivot(int p, int q) {

			// everything but row p and column q
			for (int i = 0; i <= M; i++)
				for (int j = 0; j <= M + N; j++)
					if (i != p && j != q)
						a[i][j] -= a[p][j] * a[i][q] / a[p][q];

			// zero out column q
			for (int i = 0; i <= M; i++)
				if (i != p)
					a[i][q] = 0.0;

			// scale row p
			for (int j = 0; j <= M + N; j++)
				if (j != q)
					a[p][j] /= a[p][q];
			a[p][q] = 1.0;
		}

		// return optimal objective value
		public double value() {
			return -a[M][M + N];
		}

		// return primal solution vector
		public double[] primal() {
			double[] x = new double[N];
			for (int i = 0; i < M; i++)
				if (basis[i] < N)
					x[basis[i]] = a[i][M + N];
			return x;
		}

		// return dual solution vector
		public double[] dual() {
			double[] y = new double[M];
			for (int i = 0; i < M; i++)
				y[i] = -a[M][N + i];
			return y;
		}


		// is the solution primal feasible?
		private boolean isPrimalFeasible(double[][] A, double[] b) {
			double[] x = primal();

			// check that x >= 0
			for (int j = 0; j < x.length; j++) {
				if (x[j] < 0.0) {
					return false;
				}
			}

			// check that Ax <= b
			for (int i = 0; i < M; i++) {
				double sum = 0.0;
				for (int j = 0; j < N; j++) {
					sum += A[i][j] * x[j];
				}
				if (sum > b[i] + EPSILON) {
					return false;
				}
			}
			return true;
		}

		// is the solution dual feasible?
		private boolean isDualFeasible(double[][] A, double[] c) {
			double[] y = dual();

			// check that y >= 0
			for (int i = 0; i < y.length; i++) {
				if (y[i] < 0.0) {
					return false;
				}
			}

			// check that yA >= c
			for (int j = 0; j < N; j++) {
				double sum = 0.0;
				for (int i = 0; i < M; i++) {
					sum += A[i][j] * y[i];
				}
				if (sum < c[j] - EPSILON) {
					return false;
				}
			}
			return true;
		}

		// check that optimal value = cx = yb
		private boolean isOptimal(double[] b, double[] c) {
			double[] x = primal();
			double[] y = dual();
			double value = value();

			// check that value = cx = yb
			double value1 = 0.0;
			for (int j = 0; j < x.length; j++)
				value1 += c[j] * x[j];
			double value2 = 0.0;
			for (int i = 0; i < y.length; i++)
				value2 += y[i] * b[i];
			if (Math.abs(value - value1) > EPSILON || Math.abs(value - value2) > EPSILON) {
				return false;
			}

			return true;
		}

		private boolean check(double[][] A, double[] b, double[] c) {
			return isPrimalFeasible(A, b) && isDualFeasible(A, c) && isOptimal(b, c);
		}
	}

	/**
	 * Implementation of a Priority Queue using a TreeSet, which allows for arbitrary removal in O(log n) time.
	 */
	public static class MyPriorityQueue<T extends Comparable<T>> {
		private static final long MIN_INCR = 0;
		private static final long MAX_INCR = 1000000000000000L;

		private long incr = MIN_INCR;
		private TreeSet<WrapElt<T>> ts;

		public MyPriorityQueue() {
			this.ts = new TreeSet<WrapElt<T>>();
		}

		public boolean isEmpty() {
			return ts.isEmpty();
		}

		public T poll() {
			WrapElt<T> elt = ts.last();
			if (elt == null) {
				return null;
			} else {
				ts.remove(elt);
				return elt.value;
			}
		}

		public T peek() {
			WrapElt<T> elt = ts.last();
			if (elt == null) {
				return null;
			} else {
				return elt.value;
			}
		}

		public void offer(T val) {
			WrapElt<T> elt = new WrapElt<T>(val, incr++);
			ts.add(elt);
		}

		public boolean removeFirst(T val) {
			WrapElt<T> w = new WrapElt<T>(val, MIN_INCR - 1);
			WrapElt<T> ret = ts.higher(w);
			return tryRemoveMatch(val, ret);
		}

		public boolean removeLast(T val) {
			WrapElt<T> w = new WrapElt<T>(val, MAX_INCR);
			WrapElt<T> ret = ts.lower(w);
			return tryRemoveMatch(val, ret);
		}

		private boolean tryRemoveMatch(T val, WrapElt<T> ret) {
			if (ret == null) {
				return false;
			}
			if (ret.value.compareTo(val) == 0) {
				ts.remove(ret);
				return true;
			} else {
				return false;
			}
		}

		public static class WrapElt<T extends Comparable<T>> implements Comparable<WrapElt<T>> {
			public long id;
			public T value;

			public WrapElt(T val, long i) {
				this.value = val;
				this.id = i;
			}

			public int compareTo(WrapElt<T> w) {
				int dv = this.value.compareTo(w.value);
				if (dv != 0) {
					return dv;
				} else {
					return Long.signum(this.id - w.id);
				}
			}
		}
	}

	/**
	 * Compressed binary trie that inserts and finds matches from least significant to most significant bits.
	 */
	public static class PatriciaTrie {
		public static final int BITS = 30;

		public Edge left, rite;

		public void insert(int num) {
			PatriciaTrie curr = this;
			int bits = BITS;
			while (bits > 0) {
				int b = num & 1;
				if (b == 0) {
					if (curr.left == null) {
						PatriciaTrie next = new PatriciaTrie();
						Edge e = new Edge(num, bits, next);
						curr.left = e;
						return;
					} else {
						int p = curr.left.prefixLen(num);
						if (p < curr.left.len) {
							curr.left.splitAt(p);
						}
						bits -= curr.left.len;
						num >>= curr.left.len;
						curr = curr.left.dest;
					}
				} else {
					if (curr.rite == null) {
						PatriciaTrie next = new PatriciaTrie();
						Edge e = new Edge(num, bits, next);
						curr.rite = e;
						return;
					} else {
						int p = curr.rite.prefixLen(num);
						if (p < curr.rite.len) {
							curr.rite.splitAt(p);
						}
						bits -= curr.rite.len;
						num >>= curr.rite.len;
						curr = curr.rite.dest;
					}
				}
			}
		}

		public int getClosest(int num) {
			PatriciaTrie curr = this;
			int bits = 0;
			int sum = 0;
			while (bits < BITS) {
				int b = (num >> bits) & 1;
				if ((b == 0 && curr.left != null) || (b == 1 && curr.rite == null)) {
					sum |= curr.left.val << bits;
					bits += curr.left.len;
					curr = curr.left.dest;
				} else {
					sum |= curr.rite.val << bits;
					bits += curr.rite.len;
					curr = curr.rite.dest;
				}
			}
			return sum;
		}

		public static class Edge {
			public int val;
			public int len;
			public PatriciaTrie dest;

			public Edge(int v, int ln, PatriciaTrie d) {
				this.val = v;
				this.len = ln;
				this.dest = d;
			}

			public int prefixLen(int num) {
				int cnt = 0;
				for (int i = 0; i < len; i++) {
					if (((num >> i) & 1) == ((this.val >> i) & 1)) {
						cnt++;
					} else {
						break;
					}
				}
				return cnt;
			}

			public void splitAt(int idx) {
				PatriciaTrie mid = new PatriciaTrie();
				int frontMask = 0;
				for (int i = 0; i < idx; i++) {
					frontMask |= 1 << i;
				}
				int backMask = 0;
				for (int i = idx; i < this.len; i++) {
					backMask |= 1 << i;
				}
				int front = this.val & frontMask;
				int back = (this.val & backMask) >> idx;
				Edge tmp = new Edge(back, this.len - idx, this.dest);
				int b = back & 1;
				if (b == 0) {
					mid.left = tmp;
				} else {
					mid.rite = tmp;
				}
				this.val = front;
				this.len = idx;
				this.dest = mid;
			}
		}
	}

	/**
	 * Binary trie that inserts and finds matches from most significant to least significant bits.
	 * Consider using the PatriciaTrie instead, along with the reverseBits(...) method.
	 */
	public static class BinaryTrie {
		private static final int BITS = 30;

		public BinaryTrie left, rite;

		public void insert(int n) {
			BinaryTrie curr = this;
			for (int i = BITS - 1; i >= 0; i--) {
				int b = (n >> i) & 1;
				if (b == 0) {
					if (curr.left == null) {
						curr.left = new BinaryTrie();
					}
					curr = curr.left;
				} else {
					if (curr.rite == null) {
						curr.rite = new BinaryTrie();
					}
					curr = curr.rite;
				}
			}
		}

		public int getBestMatch(int n) {
			BinaryTrie curr = this;
			int sum = 0;
			for (int i = BITS - 1; i >= 0; i--) {
				int b = (n >> i) & 1;
				if ((b == 0 && curr.left != null) || (b == 1 && curr.rite == null)) {
					curr = curr.left;
				} else {
					sum |= 1 << i;
					curr = curr.rite;
				}
			}
			return sum;
		}
	}

	/**
	 * Returns a number consisting of the reverse of the specified number of the least significant bits in the given number.
	 */
	public static int reverseBits(int num, int bits) {
		int sum = 0;
		for (int i = 0; i < bits; i++) {
			sum = (sum << 1) | (num & 1);
			num >>= 1;
		}
		return sum;
	}

	/**
	 * Returns a longest increasing subsequence in the given array.
	 */
	public int[] longestIncreasingSubsequence(int[] X) {
		final int N = X.length;
		int[] P = new int[N];
		int[] M = new int[N + 1];
		int L = 0;
		for (int i = 0; i < N; i++) {
			int lo = 1, hi = L;
			while (lo <= hi) {
				int mid = (lo + hi + 1) / 2;
				if (X[M[mid]] < X[i]) {
					lo = mid + 1;
				} else {
					hi = mid - 1;
				}
			}
			int newL = lo;
			P[i] = M[newL - 1];
			M[newL] = i;
			if (newL > L) {
				L = newL;
			}
		}
		int[] S = new int[L];
		int k = M[L];
		for (int i = L - 1; i >= 0; i--) {
			S[i] = X[k];
			k = P[k];
		}
		return S;
	}

	/**
	 * Returns the shortest distance from a point to the plane defined by three points.
	 * If the point that is the shortest distance lies outside the triangle defined by the three points, it returns infinity.
	 */
	public static double planeDist(double[] VD, double[] pa, double[] pb, double[] pc) {
		double[] u = subtract(pb, pa);
		double[] v = subtract(pc, pa);
		double[] n = cross(u, v);
		double[] w = subtract(VD, pa);
		double nn = Math.pow(norm(n), 2);
		double gamma = dot(cross(u, w), n) / nn;
		double beta = dot(cross(w, v), n) / nn;
		double alpha = 1 - gamma - beta;
		if (inRange(alpha, 0, 1) && inRange(beta, 0, 1) && inRange(gamma, 0, 1)) {
			double[] ip = add(add(scalar(pa, alpha), scalar(pb, beta)), scalar(pc, gamma));
			return norm(subtract(VD, ip));
		} else {
			return Double.POSITIVE_INFINITY;
		}
	}

	/**
	 * Returns the shortest distance from a point to a line defined by two points.
	 * If the point that is the shortest distance lies outside the line segment defined by the two points, it returns infinity.
	 */
	public static double lineDist(double[] VD, double[] pa, double[] pb) {
		double[] ab = subtract(pb, pa);
		double[] av = subtract(VD, pa);
		double[] bv = subtract(VD, pb);
		if (leq(dot(av, ab), 0) || geq(dot(bv, ab), 0)) {
			return Double.POSITIVE_INFINITY;
		} else {
			return norm(cross(ab, av)) / norm(ab);
		}
	}

	/**
	 * Returns the distance between two points.
	 */
	public static double pointDist(double[] VD, double[] p) {
		return norm(subtract(p, VD));
	}

	/**
	 * Returns the maximum of all of the input arguments.
	 */
	public static double max(double a, double... arr) {
		for (double v : arr) {
			a = Math.max(a, v);
		}
		return a;
	}

	/**
	 * Returns the minimum of all of the input arguments.
	 */
	public static double min(double a, double... arr) {
		for (double v : arr) {
			a = Math.min(a, v);
		}
		return a;
	}

	/**
	 * Returns whether the two doubles are approximately equal.
	 * You can change the value of epsilon that is the error tolerance.
	 */
	public static boolean eq(double d1, double d2) {
		return Math.abs(d1 - d2) <= 1e-12;
	}

	/**
	 * Returns whether the second double is approximately greater than or equal to the first double.
	 * This depends on the value of epsilon that is the error tolerance.
	 */
	public static boolean geq(double d1, double d2) {
		return (d1 > d2) || eq(d1, d2);
	}

	/**
	 * Returns whether the second double is approximately less than or equal to the first double.
	 * This depends on the value of epsilon that is the error tolerance.
	 */
	public static boolean leq(double d1, double d2) {
		return (d1 < d2) || eq(d1, d2);
	}

	/**
	 * Returns whether the second double is approximately greater than the first double.
	 * This depends on the value of epsilon that is the error tolerance.
	 */
	public static boolean lt(double d1, double d2) {
		return !geq(d1, d2);
	}

	/**
	 * Returns the square of a number.
	 */
	public static double sq(double d) {
		return d * d;
	}

	/**
	 * Returns whether a value lies approximately within a given range.
	 * This depends on the value of epsilon that is the error tolerance.
	 */
	public static boolean inRange(double val, double lo, double hi) {
		return geq(val, lo) && leq(val, hi);
	}

	/**
	 * Returns the sum of two 3D vectors.
	 */
	public static double[] add(double[] u, double[] v) {
		double[] r = new double[3];
		r[0] = u[0] + v[0];
		r[1] = u[1] + v[1];
		r[2] = u[2] + v[2];
		return r;
	}

	/**
	 * Returns the difference between two 3D vectors.
	 */
	public static double[] subtract(double[] u, double[] v) {
		double[] r = new double[3];
		r[0] = u[0] - v[0];
		r[1] = u[1] - v[1];
		r[2] = u[2] - v[2];
		return r;
	}

	/**
	 * Returns the cross product of two 3D vectors.
	 */
	public static double[] cross(double[] u, double[] v) {
		double[] r = new double[3];
		r[0] = (u[1] * v[2]) - (u[2] * v[1]);
		r[1] = (u[2] * v[0]) - (u[0] * v[2]);
		r[2] = (u[0] * v[1]) - (u[1] * v[0]);
		return r;
	}

	/**
	 * Returns the dot product of two 3D vectors.
	 */
	public static double dot(double[] u, double[] v) {
		return (u[0] * v[0]) + (u[1] * v[1]) + (u[2] * v[2]);
	}

	/**
	 * Returns the norm of a 3D vector.
	 */
	public static double norm(double[] v) {
		return Math.sqrt((v[0] * v[0]) + (v[1] * v[1]) + (v[2] * v[2]));
	}

	/**
	 * Returns a 3D vector multiplied by a scalar.
	 */
	public static double[] scalar(double[] v, double s) {
		double[] r = new double[3];
		r[0] = v[0] * s;
		r[1] = v[1] * s;
		r[2] = v[2] * s;
		return r;
	}

	/*
	 * DisjointSet does union-find in approximately linear time.
	 */
	public static class DisjointSet {
		private int[] rank;
		private int[] parent;

		public DisjointSet(int n) {
			this.rank = new int[n];
			this.parent = new int[n];

			for (int i = 0; i < n; ++i) {
				parent[i] = i;
			}
		}

		public int find(int x) {
			if (parent[x] != x) {
				parent[x] = find(parent[x]);
			}
			return parent[x];
		}

		public boolean union(int x, int y) {
			int xr = find(x);
			int yr = find(y);

			if (xr == yr) {
				return false;
			}

			if (rank[xr] < rank[yr]) {
				parent[xr] = yr;
			} else if (rank[xr] > rank[yr]) {
				parent[yr] = xr;
			} else {
				parent[xr] = yr;
				++rank[yr];
			}
			return true;
		}
	}

	/*
	 * QuadTree is basically a 2D segment tree for compute range sums.
	 */
	public static class QuadTree {
		public long xl, xr, yl, yr;
		public long val;

		private QuadTree parent;
		private QuadTree topLeft;
		private QuadTree topRight;
		private QuadTree botLeft;
		private QuadTree botRight;

		public QuadTree(long xlo, long xhi, long ylo, long yhi, QuadTree p) {
			this.xl = xlo;
			this.xr = xhi;
			this.yl = ylo;
			this.yr = yhi;
			this.parent = p;
		}

		public QuadTree getLeaf(long x, long y) {
			if (x == xl && x == xr && y == yl && y == yr) {
				return this;
			}
			long xm = (xl + xr) >> 1;
			long ym = (yl + yr) >> 1;
			if (x <= xm) {
				if (y <= ym) {
					if (botLeft == null) {
						botLeft = new QuadTree(xl, xm, yl, ym, this);
					}
					return botLeft.getLeaf(x, y);
				} else {
					if (topLeft == null) {
						topLeft = new QuadTree(xl, xm, ym + 1, yr, this);
					}
					return topLeft.getLeaf(x, y);
				}
			} else {
				if (y <= ym) {
					if (botRight == null) {
						botRight = new QuadTree(xm + 1, xr, yl, ym, this);
					}
					return botRight.getLeaf(x, y);
				} else {
					if (topRight == null) {
						topRight = new QuadTree(xm + 1, xr, ym + 1, yr, this);
					}
					return topRight.getLeaf(x, y);
				}
			}
		}

		public void set(long x, long y, long v) {
			QuadTree qt = getLeaf(x, y);
			qt.val = v;
			qt = qt.parent;
			while (qt != null) {
				qt.val = valueOf(qt.topLeft) + valueOf(qt.topRight) + valueOf(qt.botLeft) + valueOf(qt.botRight);
				qt = qt.parent;
			}
		}

		public long get(long xlo, long xhi, long ylo, long yhi) {
			if (xlo > xr || xhi < xl || ylo > yr || yhi < yl) {
				return defaultValue();
			}
			if (xl >= xlo && xr <= xhi && yl >= ylo && yr <= yhi) {
				return val;
			}
			return tryGet(topLeft, xlo, xhi, ylo, yhi) + tryGet(topRight, xlo, xhi, ylo, yhi) + tryGet(botLeft, xlo, xhi, ylo, yhi) + tryGet(botRight, xlo, xhi, ylo, yhi);
		}

		private static long defaultValue() {
			return 0;
		}

		private static long valueOf(QuadTree qt) {
			if (qt == null) {
				return defaultValue();
			}
			return qt.val;
		}

		private static long tryGet(QuadTree qt, long xlo, long xhi, long ylo, long yhi) {
			if (qt == null) {
				return defaultValue();
			}
			return qt.get(xlo, xhi, ylo, yhi);
		}
	}

	/**
	 * Generic implementation of the Quickselect algorithm for finding the k-th smallest element (zero-indexed).
	 * Requires a Comparator for the custom datatype.
	 * A side-effect is:
	 *   - all elements to the left of the k-th index will be less than or equal to the k-th element.
	 *   - all elements to the right of the k-th index will be greater than or equal to the k-th element.
	 */
	public static class Quickselect {
		public static <T extends Comparable<T>> T get(T[] A, int k) {
			return get(A, k, new Comparator<T>() {
				@Override
				public int compare(T lhs, T rhs) {
					return lhs.compareTo(rhs);
				}
			});
		}

		public static <T> T get(T[] A, int k, Comparator<T> cmp) {
			final int N = A.length;
			int lowerBound = 0;
			int upperBound = N - 1;
			while (lowerBound < upperBound) {
				int L = lowerBound;
				int R = upperBound;
				int M = lowerBound;
				int w = upperBound - lowerBound + 1;
				int p = lowerBound + ((((677 * lowerBound + 132241 * upperBound) % w) + w) % w);
				swap(A, p, L);
				while (M < R) {
					int compareAgainstPivot = cmp.compare(A[M + 1], A[M]);
					if (compareAgainstPivot < 0) {
						swap(A, M + 1, L);
						++L;
						++M;
					} else if (compareAgainstPivot > 0) {
						swap(A, M + 1, R);
						--R;
					} else {
						++M;
					}
				}
				if (L <= k && k <= M) {
					break;
				}
				if (k < L) {
					upperBound = L - 1;
				} else {
					lowerBound = M + 1;
				}
			}
			return A[k];
		}

		private static <T> void swap(T[] A, int i, int j) {
			T tmp = A[i];
			A[i] = A[j];
			A[j] = tmp;
		}
	}

	/**
	 * Uses rolling hash to quickly determine substring equality.
	 */
	public static class SubHash {
		private static final long P = 2147483647;
		private static final int MAXLEN = 100001;

		private static int UPTO = 1;
		private static long[] POW26 = new long[MAXLEN];
		private static long[] INV26 = new long[MAXLEN];
		static {
			POW26[0] = 1;
			POW26[1] = 26;
			INV26[0] = 1;
			INV26[1] = modInverse(26, P);
		}

		private static void loadPows(int upper) {
			for (int i = UPTO + 1; i <= upper; ++i) {
				POW26[i] = (POW26[i - 1] * POW26[1]) % P;
				INV26[i] = (INV26[i - 1] * INV26[1]) % P;
			}
			UPTO = Math.max(UPTO, upper);
		}

		private int[] S;
		private long[] H;

		public SubHash(String x) {
			loadPows(x.length());
			S = new int[x.length()];
			for (int i = 0; i < x.length(); ++i) {
				S[i] = x.charAt(i) - 'a';
			}
			H = new long[S.length + 1];
			for (int i = 0; i < S.length; ++i) {
				H[i + 1] = (H[i] + (S[i] * POW26[i])) % P;
			}
		}

		public long sub(int i, int j) {
			int len = j - i;
			long diff = (H[j] - H[i] + P) % P;
			long hash = (diff * INV26[i]) % P;
			return hash;
		}

		public int length() {
			return S.length;
		}

		public boolean subEqual(int a, int b, int len) {
			for (int i = 0; i < len; ++i) {
				if (S[a + i] != S[b + i]) {
					return false;
				}
			}
			return true;
		}
	}
}
