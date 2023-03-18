import java.util.*;

public class SquareRootDecomposition {
	/**
	 * Generic class to apply Mo's algorithm to answer offline queries.
	 * 
	 * Let:
	 *     N = upper bound of R values, e.g. all queries lie within [0, N].
	 *     K = bucket size, sqrt(N) is a good value to use.
	 *     Q = number of queries.
	 *     T = runtime of `include` or `exclude`.
	 *     R = runtime of `getResult`.
	 * 
	 * The total runtime will be T*(N*K + N*N/K) + R*Q.
	 * For optimal values of K, this is T*(N*sqrt(N)) + R*Q.
	 * 
	 * For typical usage, recommend using `answerQueries` methods, which defaults to using K = sqrt(N).
	 * 
	 * It may be possible to get slightly better performance by fine-tuning K.
	 * - Constants may receive compile-time optimizations on arithmetic division.
	 * - Powers of two may be even faster, if compiler can optimize them to right-shifts.
	 */
	public static class MosAlgorithm {
		private static final int K = 450;

		public static <T> ArrayList<T> answerQueries(Query[] queries, State<T> state) {
			return answerSortedQueries(sortQueries(queries), state);
		}

		public static <T> ArrayList<T> answerQueriesK(Query[] queries, State<T> state) {
			return answerSortedQueries(sortQueriesK(queries), state);
		}

		public static int[] answerQueries(Query[] queries, IntState state) {
			return answerSortedQueries(sortQueries(queries), state);
		}

		public static int[] answerQueriesK(Query[] queries, IntState state) {
			return answerSortedQueries(sortQueriesK(queries), state);
		}

		public static long[] answerQueries(Query[] queries, LongState state) {
			return answerSortedQueries(sortQueries(queries), state);
		}

		public static long[] answerQueriesK(Query[] queries, LongState state) {
			return answerSortedQueries(sortQueriesK(queries), state);
		}

		private static int[] answerSortedQueries(Query[] queries, IntState state) {
			final int Q = queries.length;

			int[] ans = new int[Q];
			int L = 0;
			int R = 0;
			for (Query q : queries) {
				while (R < q.R) {
					state.include(++R);
				}
				while (L > q.L) {
					state.include(--L);
				}
				while (R > q.R) {
					state.exclude(R--);
				}
				while (L < q.L) {
					state.exclude(L++);
				}
				ans[q.id] = state.getResult();
			}
			return ans;
		}

		private static long[] answerSortedQueries(Query[] queries, LongState state) {
			final int Q = queries.length;

			long[] ans = new long[Q];
			int L = 0;
			int R = 0;
			for (Query q : queries) {
				while (R < q.R) {
					state.include(++R);
				}
				while (L > q.L) {
					state.include(--L);
				}
				while (R > q.R) {
					state.exclude(R--);
				}
				while (L < q.L) {
					state.exclude(L++);
				}
				ans[q.id] = state.getResult();
			}
			return ans;
		}

		private static <T> ArrayList<T> answerSortedQueries(Query[] queries, State<T> state) {
			final int Q = queries.length;

			ArrayList<T> ans = new ArrayList<T>(Q);
			for (int i = 0; i < Q; ++i) {
				ans.add(null);
			}

			int L = 0;
			int R = 0;
			for (Query q : queries) {
				while (R < q.R) {
					state.include(++R);
				}
				while (L > q.L) {
					state.include(--L);
				}
				while (R > q.R) {
					state.exclude(R--);
				}
				while (L < q.L) {
					state.exclude(L++);
				}
				ans.set(q.id, state.getResult());
			}
			return ans;
		}

		private static Query[] sortQueries(Query[] queries) {
			int rMax = 0;
			for (Query q : queries) {
				rMax = Math.max(rMax, q.R);
			}

			Arrays.sort(queries, Query.byBucket((int) Math.sqrt(rMax + 1)));
			return queries;
		}

		private static Query[] sortQueriesK(Query[] queries) {
			Arrays.sort(queries, Query.BY_BUCKET_K);
			return queries;
		}

		public static interface Transition {
			public void include(int idx);
			public void exclude(int idx);
		}

		public static interface State<T> extends Transition {
			public T getResult();
		}

		public static interface IntState extends Transition {
			public int getResult();
		}

		public static interface LongState extends Transition {
			public long getResult();
		}

		public static class Query {
			public int id;
			public int L, R;

			public Query(int id, int L, int R) {
				this.id = id;
				this.L = L;
				this.R = R;
			}

			public static final Comparator<Query> byBucket(int bucketSize) {
				return new Comparator<Query>() {
					@Override
					public int compare(Query a, Query b) {
						int db = Integer.compare(a.L / bucketSize, b.L / bucketSize);
						if (db != 0) {
							return db;
						}
						return Integer.compare(a.R, b.R);
					}
				};
			}

			private static final Comparator<Query> BY_BUCKET_K = new Comparator<Query>() {
				@Override
				public int compare(Query a, Query b) {
					int db = Integer.compare(a.L / K, b.L / K);
					if (db != 0) {
						return db;
					}
					return Integer.compare(a.R, b.R);
				}
			};
		}
	}
}
