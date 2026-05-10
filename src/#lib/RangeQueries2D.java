public class RangeQueries2D {
	/**
	 * Performs static 2D range sums.
	 * Note that arguments are inclusive-lower exclusive-upper for function sumRange(int r1, int c1, int r2, int c2).
	 */
	public static class SubmatrixSumInt {
		private FlatMatrixInt mat;
		private int N, M;

		public SubmatrixSumInt(int[][] matrix) {
			N = matrix.length;
			M = matrix[0].length;
			mat = new FlatMatrixInt(N + 1, M + 1);
			for (int i = 0; i < N; ++i) {
				for (int j = 0; j < M; ++j) {
					mat.set(i + 1, j + 1, mat.get(i + 1, j) + mat.get(i, j + 1) - mat.get(i, j) + matrix[i][j]);
				}
			}
		}

		public int sumRange(int r1, int c1, int r2, int c2) {
			return mat.get(r2, c2) - mat.get(r2, c1) - mat.get(r1, c2) + mat.get(r1, c1);
		}

		public int sumRange(Range r) {
			return sumRange(r.r1, r.c1, r.r2, r.c2);
		}

		public static class Range {
			public int r1, c1, r2, c2;

			public Range(int r1, int c1, int r2, int c2) {
				this.r1 = r1;
				this.c1 = c1;
				this.r2 = r2;
				this.c2 = c2;
			}
		}

		private static class FlatMatrixInt {
			private final int M;
			private final int[] A;

			public FlatMatrixInt(int N, int M) {
				this.M = M;
				this.A = new int[N * M];
			}

			private int index(int i, int j) {
				return i * M + j;
			}

			public void set(int i, int j, int x) {
				A[index(i, j)] = x;
			}

			public int get(int i, int j) {
				return A[index(i, j)];
			}
		}
	}

	/**
	 * Performs static 2D range sums.
	 * Note that arguments are inclusive-lower exclusive-upper for function sumRange(int r1, int c1, int r2, int c2).
	 */
	public static class SubmatrixSum {
		private long[][] sum;
		private int N, M;

		public SubmatrixSum(int[][] matrix) {
			this(toLongMatrix(matrix));
		}

		public SubmatrixSum(long[][] matrix) {
			if (matrix.length == 0 || matrix[0].length == 0) {
				return;
			}
			N = matrix.length;
			M = matrix[0].length;
			sum = new long[N + 1][M + 1];
			for (int i = 0; i < N; ++i) {
				for (int j = 0; j < M; ++j) {
					sum[i + 1][j + 1] = sum[i + 1][j] + sum[i][j + 1] - sum[i][j] + matrix[i][j];
				}
			}
		}

		public long sumRange(int r1, int c1, int r2, int c2) {
			if (sum == null) {
				return 0;
			}
			return sum[r2][c2] - sum[r2][c1] - sum[r1][c2] + sum[r1][c1];
		}

		public long sumRange(Range r) {
			return sumRange(r.r1, r.c1, r.r2, r.c2);
		}

		public static class Range {
			public int r1, c1, r2, c2;

			public Range(int r1, int c1, int r2, int c2) {
				this.r1 = r1;
				this.c1 = c1;
				this.r2 = r2;
				this.c2 = c2;
			}
		}

		private static long[][] toLongMatrix(int[][] mat) {
			long[][] ans = new long[mat.length][];
			for (int i = 0; i < mat.length; ++i) {
				ans[i] = new long[mat[i].length];
				for (int j = 0; j < mat[i].length; ++j) {
					ans[i][j] = mat[i][j];
				}
			}
			return ans;
		}
	}
}
