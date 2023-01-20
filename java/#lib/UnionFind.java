public class UnionFind {
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

	public static class DisjointSetMerge<T> {
		private int[] rank;
		private int[] parent;
		private ArrayList<T> data;
		private Merger<T> merger;

		public DisjointSetMerge(int n, Merger<T> m) {
			this(n, null, m);
		}

		public DisjointSetMerge(int n, T defaultValue, Merger<T> m) {
			this.rank = new int[n];
			this.parent = new int[n];
			this.data = new ArrayList<>();
			while (this.data.size() < n) {
				this.data.add(defaultValue);
			}
			this.merger = m;

			for (int i = 0; i < n; ++i) {
				parent[i] = i;
			}
		}

		public T getData(int x) {
			return data.get(find(x));
		}

		public void setData(int x, T d) {
			data.set(find(x), d);
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

			data.set(find(xr), merger.merge(data.get(xr), data.get(yr)));

			return true;
		}

		public static interface Merger<T> {
			public T merge(T a, T b);
		}
	}

	/**
	 * Solves offline dynamic connectivity problem in O(K log K), where K is number of operations.
	 * 
	 * Returns an array of query results, where the i-th result is the number of connected components at the i-th count query.
	 * Note that the results array is zero-indexed, and the length of the array is the number of count queries.
	 * 
	 * See for more details:
	 *     https://codeforces.com/gym/100551/problem/A
	 * 
	 * Code is heavily adapted from C++ solution:
	 *     https://codeforces.com/gym/100551/submission/157416340
	 * 
	 * EXAMPLE USAGE:
	 *     DynamicConnectivityOffline.Queries queries = new DynamicConnectivityOffline.Queries();
	 *     queries.appendAddQuery(1, 2);
	 *     queries.appendCountQuery();
	 *     queries.appendRemoveQuery(1, 2);
	 *     ... add more queries ...
	 *     int[] result = DynamicConnectivityOffline.solve(N, queries);
	 *     
	 * NOTE:  Removing an edge that does not exist is NO-OP.
	 */
	public static class DynamicConnectivityOffline {
		public static int[] solve(int n, Queries q) {
			Solver s = new Solver(n, q.sz, q.tp, q.u, q.v);
			s.fullSolve();
			return Arrays.copyOf(s.out, s.oi);
		}

		public static int[] solve(int n, int k, char[] tp, int[] u, int[] v) {
			Solver s = new Solver(n, k, tp, u, v);
			s.fullSolve();
			return Arrays.copyOf(s.out, s.oi);
		}

		public static class Queries {
			private char[] tp;
			private int[] u, v;
			private int sz;

			public Queries() {
				this(2);
			}

			public Queries(int k) {
				this.tp = new char[k];
				this.u = new int[k];
				this.v = new int[k];
			}

			public void appendAddQuery(int a, int b) {
				append(ADD, a, b);
			}

			public void appendRemoveQuery(int a, int b) {
				append(REMOVE, a, b);
			}

			public void appendCountQuery() {
				append(COUNT, -1, -1);
			}

			private void append(char c, int a, int b) {
				if (sz >= tp.length) {
					tp = Arrays.copyOf(tp, tp.length << 1);
					u = Arrays.copyOf(u, u.length << 1);
					v = Arrays.copyOf(v, v.length << 1);
				}
				tp[sz] = c;
				u[sz] = a;
				v[sz] = b;
				++sz;
			}
		}

		private static final char ADD = '+';
		private static final char REMOVE = '-';
		private static final char COUNT = '?';

		private static class Solver {
			private final int n;
			private final int k;

			private final int[] u;
			private final int[] v;
			private final char[] tp;

			private final int[] nxt;
			private final int[] sz;

			private int res = 0;

			private final int[] stk;
			private int si = 0;

			private final int[] inv;

			private final int[] out;
			private int oi = 0;

			private final Map<IntPair, Integer> mp = new HashMap<>();

			public Solver(int n, int k, char[] tp, int[] u, int[] v) {
				this.n = n;
				this.k = k;

				this.nxt = new int[n];
				this.sz = new int[n];

				this.u = u;
				this.v = v;
				this.tp = tp;

				this.inv = new int[k];
				Arrays.fill(inv, -1);

				this.stk = new int[k];
				this.out = new int[k];
			}

			int root(int u) {
				while (nxt[u] != u) {
					u = nxt[u];
				}
				return u;
			}

			void join(int u, int v) {
				u = root(u);
				v = root(v);

				if (u == v) {
					return;
				}

				if (sz[u] < sz[v]) {
					int tmp = u;
					u = v;
					v = tmp;
				}

				stk[si++] = v;
				sz[u] += sz[v];
				--res;
				nxt[v] = u;
			}

			void undo(int newSz) {
				while (si > newSz) {
					int v = stk[--si];
					sz[nxt[v]] -= sz[v];
					nxt[v] = v;
					++res;
				}
			}

			void solve(int l, int r) {
				boolean ok = false;
				for (int i = l; i <= r; ++i) {
					ok |= (tp[i] == COUNT);
				}
				if (!ok) {
					return;
				}

				if (l == r) {
					if (tp[l] == COUNT) {
						out[oi++] = res;
					}
					return;
				}

				int mid = (l + r) >> 1;
				int init = si;

				for (int i = mid + 1; i <= r; ++i) {
					if (tp[i] != COUNT && inv[i] < l && inv[i] >= 0) {
						join(u[i], v[i]);
					}
				}
				solve(l, mid);
				undo(init);

				for (int i = l; i <= mid; ++i) {
					if (tp[i] != COUNT && (inv[i] > r || inv[i] < 0)) {
						join(u[i], v[i]);
					}
				}
				solve(mid + 1, r);
				undo(init);
			}

			void fullSolve() {
				for (int i = 0; i < n; ++i) {
					nxt[i] = i;
					sz[i] = 1;
				}
				for (int i = 0; i < k; ++i) {
					if (tp[i] == COUNT) {
						continue;
					}

					if (u[i] > v[i]) {
						int tmp = u[i];
						u[i] = v[i];
						v[i] = tmp;
					}

					IntPair key = new IntPair(u[i], v[i]);
					if (tp[i] == REMOVE) {
						inv[i] = mp.get(key);
						inv[inv[i]] = i;
						mp.remove(key);
					} else {
						mp.put(key, i);
					}
				}
				res = n;
				solve(0, k - 1);
			}
		}

		private static class IntPair implements Comparable<IntPair> {
			public int first, second;

			public IntPair(int first, int second) {
				this.first = first;
				this.second = second;
			}

			@Override
			public int hashCode() {
				return first * 131071 + second * 31;
			}

			@Override
			public boolean equals(Object obj) {
				IntPair ip = (IntPair) obj;
				return first == ip.first && second == ip.second;
			}

			@Override
			public int compareTo(IntPair ip) {
				int df = Integer.compare(first, ip.first);
				if (df != 0) {
					return df;
				}
				return Integer.compare(second, ip.second);
			}
		}
	}
}
