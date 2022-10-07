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
}
