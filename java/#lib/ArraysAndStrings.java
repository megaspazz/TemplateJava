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
			Integer[] arr = Arrays.stream(A).mapToObj(Integer::valueOf).toArray(Integer[]::new);
			Arrays.sort(arr);
			for (int i = 0; i < A.length; ++i) {
				A[i] = arr[i];
			}
			return A;
		}

		public static long[] longs(long[] A) {
			Long[] arr = Arrays.stream(A).mapToObj(Long::valueOf).toArray(Long[]::new);
			Arrays.sort(arr);
			for (int i = 0; i < A.length; ++i) {
				A[i] = arr[i];
			}
			return A;
		}

		public static double[] doubles(double[] A) {
			Double[] arr = Arrays.stream(A).mapToObj(Double::valueOf).toArray(Double[]::new);
			Arrays.sort(arr);
			for (int i = 0; i < A.length; ++i) {
				A[i] = arr[i];
			}
			return A;
		}
	}
}
