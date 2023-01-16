import java.util.*;

public class Graphs {
	/**
	 * Finds an Eulerian path in a graph, if one exists.
	 */
	// EXAMPLE USAGE:
	/*
		EulerianPath.Graph g = new EulerianPath.Graph(numNodes);
		g.addDirectedEdge(a, b);
		g.addDirectedEdge(u, v);
		g.addDirectedEdge(x, y);
		int[] order = EulerianPath.getNodeOrder(g);
	 */
	public static class EulerianPath {
		public static int[] getNodeOrder(Graph g) {
			ArrayDeque<Node> toProcess = new ArrayDeque<>();
			
			final int N = g.nodes.length;
			int[] indegRem = new int[g.nodes.length];
			int[] selfLoops = new int[N];
			for (Node u : g.nodes) {
				for (Edge e : u.edges) {
					++indegRem[e.dest.id];
					if (e.dest.id == u.id) {
						++selfLoops[u.id];
					}
				}
			}

			Node start = null;
			for (Node u : g.nodes) {
				if (u.edges.size() > indegRem[u.id]) {
					start = u;
					break;
				}
			}
			if (start == null) {
				for (Node u : g.nodes) {
					if (((u.edges.size() - selfLoops[u.id]) & 1) != 0) {
						start = u;
						break;
					}
				}
			}
			if (start == null) {
				start = g.nodes[0];
			}

			int expectedLength = g.numEdges + 1;
			int[] ans = new int[expectedLength];
			int ansIdx = 0;
			boolean[] used = new boolean[g.numEdges];
			toProcess.push(start);
			while (!toProcess.isEmpty()) {
				Node curr = toProcess.peek();

				Edge e = getNextEdge(curr, indegRem, used);
				if (e == null) {
					toProcess.poll();
					ans[ansIdx++] = curr.id;
					continue;
				}

				used[e.id] = true;
				toProcess.push(e.dest);
			}

			if (ansIdx != expectedLength) {
				return null;
			}
			return ans;
		}

		private static Edge getNextEdge(Node u, int[] indegRem, boolean[] used) {
			while (indegRem[u.id] > 0) {
				Edge e = u.edges.get(--indegRem[u.id]);
				if (!used[e.id]) {
					return e;
				}
			}
			return null;
		}

		public static class Node {
			public int id;
			public ArrayList<Edge> edges = new ArrayList<>();

			public Node(int id) {
				this.id = id;
			}
		}

		public static class Edge {
			public int id;
			public Node dest;

			public Edge(int id, Node dest) {
				this.id = id;
				this.dest = dest;
			}
		}

		public static class Graph {
			private Node[] nodes;
			private int numEdges = 0;

			public Graph(int numNodes) {
				nodes = new Node[numNodes];
				for (int i = 0; i < numNodes; ++i) {
					nodes[i] = new Node(i);
				}
			}

			public void addDirectedEdge(int u, int v) {
				nodes[u].edges.add(new Edge(numEdges++, nodes[v]));
			}

			public void addUndirectedEdge(int u, int v) {
				nodes[u].edges.add(new Edge(numEdges, nodes[v]));
				nodes[v].edges.add(new Edge(numEdges, nodes[u]));
				++numEdges;
			}
		}
	}

	/**
	 * Computes shortest path between all pairs of vertices.
	 */
	public static class FloydWarshall {
		/**
		 * NOTE:  Using INF to represent non-existent edges must be in the range [-2^30, 2^30).
		 * NOTE:  The longest shortest-path must be less than INF.
		 */
		public static int[][] ints(int[][] dist) {
			final int N = dist.length;

			int[][] best = new int[N][];
			for (int i = 0; i < N; ++i) {
				best[i] = Arrays.copyOf(dist[i], N);
			}
			for (int k = 0; k < N; ++k) {
				for (int i = 0; i < N; ++i) {
					for (int j = 0; j < N; ++j) {
						int d = best[i][k] + best[k][j];
						if (d < best[i][j]) {
							best[i][j] = d;
						}
					}
				}
			}
			return best;
		}

		/**
		 * NOTE:  Using INF to represent non-existent edges must be in the range [-2^62, 2^62).
		 * NOTE:  The longest shortest-path must be less than INF.
		 */
		public static long[][] longs(long[][] dist) {
			final int N = dist.length;

			long[][] best = new long[N][];
			for (int i = 0; i < N; ++i) {
				best[i] = Arrays.copyOf(dist[i], N);
			}
			for (int k = 0; k < N; ++k) {
				for (int i = 0; i < N; ++i) {
					for (int j = 0; j < N; ++j) {
						long d = best[i][k] + best[k][j];
						if (d < best[i][j]) {
							best[i][j] = d;
						}
					}
				}
			}
			return best;
		}
	}

	/**
	 * Computes lowest common ancestor of nodes in a tree.
	 *
	 * Usage:
	 *   1)  Build bidirectional graph on LCA.Node, populating only "index" and "neighbors" fields.
	 *   2)  With root node, initialize LCA in O(N log N) time:
	 *           LCA lca = new LCA(root);
	 *   3)  To find lowest common ancestor of two nodes in O(1) time:
	 *           LCA.Node ancestor = lca.lowestCommonAncestor(u, v);
	 */
	public static class LCA {
		private ArrayList<Integer> idx = new ArrayList<>();
		private SparseTable tbl;

		public LCA(Node root) {
			root.initRootTree();
			ArrayList<Node> lst = traverse(root);
			tbl = new SparseTable(lst.toArray(new Node[0]));
		}

		public Node lowestCommonAncestor(Node u, Node v) {
			return lowestCommonAncestor(u.index, v.index);
		}

		public Node lowestCommonAncestor(int u, int v) {
			int left = Math.min(idx.get(u), idx.get(v));
			int rite = Math.max(idx.get(u), idx.get(v));
			return tbl.getMin(left, rite);
		}

		private ArrayList<Node> traverse(Node start) {
			ArrayList<Node> lst = new ArrayList<>();
			Stack<Node> nodes = new Stack<>();
			Stack<Integer> indices = new Stack<>();
			nodes.push(start);
			indices.push(0);
			while (!nodes.isEmpty()) {
				Node u = nodes.pop();
				int i = indices.pop();
				while (u.index >= idx.size()) {
					idx.add(null);
				}
				if (idx.get(u.index) == null) {
					idx.set(u.index, lst.size());
				}
				lst.add(u);
				if (i < u.neighbors.size()) {
					nodes.push(u);
					indices.push(i + 1);
					nodes.push(u.neighbors.get(i));
					indices.push(0);
					++i;
				}
			}
			return lst;
		}

		public static class Node {
			public ArrayList<Node> neighbors = new ArrayList<>();
			public int index;
			public int level;

			public Node(int idx) {
				this.index = idx;
			}

			public void initRootTree() {
				Stack<Node> nodes = new Stack<>();
				Stack<Node> parents = new Stack<>();
				Stack<Integer> levels = new Stack<>();
				nodes.push(this);
				parents.push(null);
				levels.push(0);
				while (!nodes.isEmpty()) {
					Node u = nodes.pop();
					Node p = parents.pop();
					int lvl = levels.pop();
					u.neighbors.remove(p);
					u.level = lvl;
					for (Node v : u.neighbors) {
						nodes.push(v);
						parents.push(u);
						levels.push(lvl + 1);
					}
				}
			}
		}

		public static Node nodeMin(Node a, Node b) {
			if (a == null && b == null) {
				return null;
			} else if (a == null) {
				return b;
			} else if (b == null) {
				return a;
			} else {
				return (a.level < b.level) ? a : b;
			}
		}

		public static class SparseTable {
			public Node[] arr;
			public Node[][] table;
			public int num;
			public int pow;

			public SparseTable(Node[] a) {
				int n = a.length;
				int p = Math.max(1, Integer.SIZE - Integer.numberOfLeadingZeros(n));

				this.arr = a;
				this.num = n;
				this.pow = p;
				this.table = new Node[n][p];

				for (int i = 0; i < n; i++) {
					this.table[i][0] = a[i];
				}
				for (int j = 1; j < p; j++) {
					for (int i = 0; i < n; i++) {
						int off = (1 << (j - 1));
						Node next = (i + off >= n) ? null : this.table[i + off][j - 1];
						this.table[i][j] = nodeMin(this.table[i][j - 1], next);
					}
				}
			}

			public Node getMin(int lower, int upper) {
				int k = 31 - Integer.numberOfLeadingZeros(upper - lower + 1);
				int off = (1 << k);
				return nodeMin(this.table[lower][k], this.table[upper - off + 1][k]);
			}
		}
	}

	/**
	 * Performs topological sort based on constructed input graph.
	 * Returns any valid order if possible.
	 * Returns null if impossible, e.g. has a cycle in the graph.
	 *
	 * EXAMPLE USAGE, SEQUENTIAL INT KEY (faster):
	 *   TopologicalSort.IntGraph g = new TopologicalSort.IntGraph(N_MAX);
	 *   g.addEdge(1, 2);
	 *   g.addEdge(1, 3);
	 *   ... add more edges ...
	 *   int[] order = TopologicalSort.sort(g);
	 *   
	 *
	 * EXAMPLE USAGE, GENERIC (slower):
	 *   TopologicalSort.Graph<String, Object> g = new TopologicalSort.Graph();
	 *   g.addEdge("abc", "xyz");
	 *   g.addEdge("abc", "123456");
	 *   ... add more edges ...
	 *   List<TopologicalSort.NodeData<Integer, Object>> order = TopologicalSort.sort(g);
	 */
	public static class TopologicalSort {
		public static <K, V> List<NodeData<K, V>> sort(Graph<K, V> g) {
			HashMap<K, Integer> indegRem = new HashMap<>();
			Queue<Node<K, V>> q = new ArrayDeque<>();
			for (Node<K, V> u : g.nodes.values()) {
				if (u.indegree == 0) {
					q.offer(u);
				} else {
					indegRem.put(u.data.id, u.indegree);
				}
			}
			List<NodeData<K, V>> order = new ArrayList<>();
			while (!q.isEmpty()) {
				Node<K, V> u = q.poll();
				order.add(u.data);
				for (Node<K, V> v : u.next) {
					int vIndegree = indegRem.get(v.data.id);
					if (vIndegree == 1) {
						indegRem.remove(v.data.id);
						q.offer(v);
					} else {
						indegRem.put(v.data.id, vIndegree - 1);
					}
				}
			}
			if (order.size() != g.nodes.size()) {
				return null;
			}
			return order;
		}

		public static int[] sort(IntGraph g) {
			final int N = g.next.length;

			IntDeque q = new IntDeque();
			int[] indegRem = new int[N];
			int wantSize = 0;
			for (int i = 0; i < N; ++i) {
				if (g.next[i] == null) {
					continue;
				}

				++wantSize;
				if (g.indegree[i] == 0) {
					q.offer(i);
				} else {
					indegRem[i] = g.indegree[i];
				}
			}

			IntDeque order = new IntDeque();
			while (!q.isEmpty()) {
				int u = q.poll();
				order.add(u);
				for (int v : g.next[u].toArray()) {
					if (--indegRem[v] == 0) {
						q.offer(v);
					}
				}
			}
			if (order.size() != wantSize) {
				return null;
			}
			return order.toArray();
		}

		public static class NodeData<K, V> {
			public K id;
			public V value;

			public NodeData(K i, V val) {
				id = i;
				value = val;
			}

			public NodeData(K i) {
				id = i;
			}
		}

		public static class Node<K, V> {
			public NodeData<K, V> data;

			private ArrayList<Node<K, V>> next = new ArrayList<>();
			private int indegree = 0;

			public Node(K i, V val) {
				data = new NodeData<>(i, val);
			}

			public Node(K i) {
				data = new NodeData<>(i);
			}
		}

		public static class Graph<K, V> {
			public HashMap<K, Node<K, V>> nodes = new HashMap<>();

			public Node<K, V> getOrCreate(K key) {
				Node<K, V> u = nodes.getOrDefault(key, null);
				if (u == null) {
					u = new Node<>(key);
					nodes.put(key, u);
				}
				return u;
			}

			public Node<K, V> setData(K key, V val) {
				Node<K, V> u = getOrCreate(key);
				u.data.value = val;
				return u;
			}

			public void addDirectedEdge(K kFrom, K kTo) {
				Node<K, V> to = getOrCreate(kTo);
				++to.indegree;
				getOrCreate(kFrom).next.add(to);
			}

			public static <V> Graph<Integer, V> withIntKeys(int loInclusive, int hiExclusive) {
				Graph<Integer, V> g = new Graph<>();
				for (int i = loInclusive; i < hiExclusive; ++i) {
					g.getOrCreate(i);
				}
				return g;
			}
		}

		public static class IntGraph {
			private IntDeque[] next;
			private int[] indegree;

			public IntGraph(int n) {
				this.next = new IntDeque[n];
				this.indegree = new int[n];
			}

			public IntDeque getOrCreate(int id) {
				IntDeque u = next[id];
				if (u == null) {
					u = new IntDeque();
					next[id] = u;
				}
				return u;
			}

			public void addDirectedEdge(int idFrom, int idTo) {
				getOrCreate(idTo);
				++indegree[idTo];
				getOrCreate(idFrom).add(idTo);
			}
		}

		/*
		 * Trimmed down version of IntDeque for convenience.
		 * 
		 * See Lists.java for full implementation:
		 * https://github.com/megaspazz/TemplateJava/blob/master/java/%23lib/Lists.java
		 */
		private static class IntDeque {
			private int[] arr;
			private int off;
			private int len;

			public IntDeque() {
				this(2);
			}

			public IntDeque(int capacity) {
				this.arr = new int[capacity];
			}

			public void addLast(int x) {
				if (len == arr.length) {
					increaseCapacity();
				}
				int idx = index(off + len);
				arr[idx] = x;
				++len;
			}

			public int peekFirst() {
				return arr[off];
			}

			public int removeFirst() {
				int ans = peekFirst();
				off = index(off + 1);
				--len;
				return ans;
			}

			public void add(int x) {
				addLast(x);
			}

			public void offer(int x) {
				addLast(x);
			}

			public int poll() {
				return removeFirst();
			}

			public int size() {
				return len;
			}

			public boolean isEmpty() {
				return size() == 0;
			}

			public int[] toArray() {
				if (len == 0) {
					return new int[0];
				}
				int idx = index(off + len);
				if (idx > off) {
					return Arrays.copyOfRange(arr, off, idx);
				}
				int[] A = new int[len];
				int endLen = arr.length - off;
				int startLen = len - endLen;
				System.arraycopy(arr, off, A, 0, endLen);
				System.arraycopy(arr, 0, A, endLen, startLen);
				return A;
			}

			private void increaseCapacity() {
				int[] next = new int[arr.length << 1];
				int endLen = arr.length - off;
				System.arraycopy(arr, off, next, 0, endLen);
				System.arraycopy(arr, 0, next, endLen, off);
				arr = next;
				off = 0;
			}

			private int index(int i) {
				if (i >= arr.length) {
					i -= arr.length;
				} else if (i < 0) {
					i += arr.length;
				}
				return i;
			}
		}
	}
}
