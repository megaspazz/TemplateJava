public class RangeQueries2D {
	/**
	 * Performs static 2D range sums.
	 * Note that arguments are inclusive-lower exclusive-upper for function sumRange(int r1, int c1, int r2, int c2).
	 */
	public static class SubmatrixSumInt {
		private int[][] sum;
		private int N, M;

		public SubmatrixSumInt(int[][] matrix) {
			if (matrix.length == 0 || matrix[0].length == 0) {
				return;
			}
			N = matrix.length;
			M = matrix[0].length;
			sum = new int[N + 1][M + 1];
			for (int i = 0; i < N; ++i) {
				for (int j = 0; j < M; ++j) {
					sum[i + 1][j + 1] = sum[i + 1][j] + sum[i][j + 1] - sum[i][j] + matrix[i][j];
				}
			}
		}

		public int sumRange(int r1, int c1, int r2, int c2) {
			if (sum == null) {
				return 0;
			}
			return sum[r2][c2] - sum[r2][c1] - sum[r1][c2] + sum[r1][c1];
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
