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
}
