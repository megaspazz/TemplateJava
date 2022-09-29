import java.util.*;

public class ArraysAndStrings {
	public static class Quickselect {
		private static final Random RNG = new Random();

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
				int p = lowerBound + RNG.nextInt(upperBound - lowerBound + 1);
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

		public static int get(int[] A, int k) {
			final int N = A.length;
			int lowerBound = 0;
			int upperBound = N - 1;
			while (lowerBound < upperBound) {
				int L = lowerBound;
				int R = upperBound;
				int M = lowerBound;
				int p = lowerBound + RNG.nextInt(upperBound - lowerBound + 1);
				
				swap(A, p, L);
				while (M < R) {
					if (A[M + 1] < A[M]) {
						swap(A, M + 1, L);
						++L;
						++M;
					} else if (A[M + 1] > A[M]) {
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

		private static void swap(int[] A, int i, int j) {
			int tmp = A[i];
			A[i] = A[j];
			A[j] = tmp;
		}
	}

	public static class Sort {
		public static int[] ints(int[] A) {
			return mergesortInts(A);
		}

		public static long[] longs(long[] A) {
			final int N = A.length;

			Long[] arr = new Long[N];
			for (int i = 0; i < N; ++i) {
				arr[i] = A[i];
			}

			Arrays.sort(arr);

			for (int i = 0; i < N; ++i) {
				A[i] = arr[i];
			}
			return A;
		}

		public static double[] doubles(double[] A) {
			final int N = A.length;

			Double[] arr = new Double[N];
			for (int i = 0; i < N; ++i) {
				arr[i] = A[i];
			}

			Arrays.sort(arr);

			for (int i = 0; i < N; ++i) {
				A[i] = arr[i];
			}
			return A;
		}
		
		private static int[] mergesortInts(int[] A) {
			final int N = A.length;
			
			int[] L = new int[N >> 1];
			int[] R = new int[N >> 1];
			for (int w = 2; w < N; w <<= 1) {
				int hw = w << 1;
				for (int i = 0; i + hw < N; i += w) {
					int j = i + hw;
					int k = Math.min(i + w, N);
					
					System.arraycopy(A, i, L, 0, hw);
					System.arraycopy(A, j, R, 0, k - j);
					
					int p = i;
					int a = 0;
					int b = 0;
					while (a < L.length && b < R.length) {
						if (L[a] < R[b]) {
							A[p++] = L[a++];
						} else {
							A[p++] = R[b++];
						}
					}
					if (a < L.length) {
						System.arraycopy(L, a, A, p, L.length - a);
					} else {
						System.arraycopy(R, b, A, p, R.length - b);
					}
				}
			}
			return A;
		}
	}

	/**
	 * Generic binary search to find the first or last value resulting in a matching condition.
	 */
	// EXAMPLE USAGE (find insertion index in sorted array `A`):
	/*
		int insertionIndex = BinarySearch.firstThat(0, A.length, new BinarySearch.IntCheck() {
			@Override
			public boolean valid(int index) {
				return A[index] > mid;
			}
		});
	 */
	public static class BinarySearch {
		// Finds the left-most value that satisfies the IntCheck in the range [L, R).
		// It will return R if the nothing in the range satisfies the check.
		public static int firstThat(int L, int R, IntCheck check) {
			while (L < R) {
				int M = (L >> 1) + (R >> 1) + (L & R & 1);
				if (check.valid(M)) {
					R = M;
				} else {
					L = M + 1;
				}
			}
			return L;
		}

		// Finds the right-most value that satisfies the IntCheck in the range [L, R).
		// It will return L - 1 if nothing in the range satisfies the check.
		public static int lastThat(int L, int R, IntCheck check) {
			int firstValue = firstThat(L, R, new IntCheck() {
				@Override
				public boolean valid(int value) {
					return !check.valid(value);
				}
			});
			return firstValue - 1;
		}

		// Finds the left-most value that satisfies the LongCheck in the range [L, R).
		public static long firstThat(long L, long R, LongCheck check) {
			while (L < R) {
				long M = (L >> 1) + (R >> 1) + (L & R & 1);
				if (check.valid(M)) {
					R = M;
				} else {
					L = M + 1;
				}
			}
			return L;
		}

		// Finds the right-most value that satisfies the IntCheck in the range [L, R).
		// It will return L - 1 if nothing in the range satisfies the check.
		public static long lastThat(long L, long R, LongCheck check) {
			long firstValue = firstThat(L, R, new LongCheck() {
				@Override
				public boolean valid(long value) {
					return !check.valid(value);
				}
			});
			return firstValue - 1;
		}

		public static interface LongCheck {
			public boolean valid(long value);
		}

		public static interface IntCheck {
			public boolean valid(int value);
		}
	}

	private static class Shuffle {
		private static final ThreadLocalRandom RNG = ThreadLocalRandom.current();
		
		public static int[] ints(int[] A) {
			for (int i = A.length - 1; i > 0; --i) {
				int j = RNG.nextInt(0, i + 1);
				swapInts(A, i, j);
			}
			return A;
		}
		
		private static void swapInts(int[] A, int i, int j) {
			int tmp = A[i];
			A[i] = A[j];
			A[j] = tmp;
		}
	}
}
