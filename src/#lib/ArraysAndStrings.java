import java.util.*;
import java.util.concurrent.*;

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

		public static long get(long[] A, int k) {
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

		public static long getReversed(long[] A, int k) {
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
					if (A[M + 1] > A[M]) {
						swap(A, M + 1, L);
						++L;
						++M;
					} else if (A[M + 1] < A[M]) {
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

		private static void swap(long[] A, int i, int j) {
			long tmp = A[i];
			A[i] = A[j];
			A[j] = tmp;
		}
	}

	public static class Sort {
		public static int[] ints(int[] A) {
			final int N = A.length;

			int[] buf = new int[N];
			for (int lw = 1; lw < N; lw <<= 1) {
				int w = lw << 1;
				for (int i = 0; i + lw <= N; i += w) {
					int k = Math.min(i + w, N);

					int segLen = k - i;
					System.arraycopy(A, i, buf, 0, segLen);

					int p = i;
					int a = 0;
					int b = lw;
					while (a < lw && b < segLen) {
						if (buf[a] < buf[b]) {
							A[p++] = buf[a++];
						} else {
							A[p++] = buf[b++];
						}
					}
					if (a < lw) {
						System.arraycopy(buf, a, A, p, lw - a);
					} else {
						System.arraycopy(buf, b, A, p, segLen - b);
					}
				}
			}
			return A;
		}

		public static long[] longs(long[] A) {
			final int N = A.length;

			long[] buf = new long[N];
			for (int lw = 1; lw < N; lw <<= 1) {
				int w = lw << 1;
				for (int i = 0; i + lw <= N; i += w) {
					int k = Math.min(i + w, N);

					int segLen = k - i;
					System.arraycopy(A, i, buf, 0, segLen);

					int p = i;
					int a = 0;
					int b = lw;
					while (a < lw && b < segLen) {
						if (buf[a] < buf[b]) {
							A[p++] = buf[a++];
						} else {
							A[p++] = buf[b++];
						}
					}
					if (a < lw) {
						System.arraycopy(buf, a, A, p, lw - a);
					} else {
						System.arraycopy(buf, b, A, p, segLen - b);
					}
				}
			}
			return A;
		}

		public static double[] doubles(double[] A) {
			final int N = A.length;

			double[] buf = new double[N];
			for (int lw = 1; lw < N; lw <<= 1) {
				int w = lw << 1;
				for (int i = 0; i + lw <= N; i += w) {
					int k = Math.min(i + w, N);

					int segLen = k - i;
					System.arraycopy(A, i, buf, 0, segLen);

					int p = i;
					int a = 0;
					int b = lw;
					while (a < lw && b < segLen) {
						if (buf[a] < buf[b]) {
							A[p++] = buf[a++];
						} else {
							A[p++] = buf[b++];
						}
					}
					if (a < lw) {
						System.arraycopy(buf, a, A, p, lw - a);
					} else {
						System.arraycopy(buf, b, A, p, segLen - b);
					}
				}
			}
			return A;
		}

		public static class Reversed {
			public static int[] ints(int[] A) {
				Sort.ints(A);
				for (int i = 0, j = A.length - 1; i < j; ++i, --j) {
					int tmp = A[i];
					A[i] = A[j];
					A[j] = tmp;
				}
				return A;
			}

			public static long[] longs(long[] A) {
				Sort.longs(A);
				for (int i = 0, j = A.length - 1; i < j; ++i, --j) {
					long tmp = A[i];
					A[i] = A[j];
					A[j] = tmp;
				}
				return A;
			}

			public static double[] doubles(double[] A) {
				Sort.doubles(A);
				for (int i = 0, j = A.length - 1; i < j; ++i, --j) {
					double tmp = A[i];
					A[i] = A[j];
					A[j] = tmp;
				}
				return A;
			}
		}
	}

	public static class BucketSort {
		public static <T> T[] sortReversed(T[] items, IntConverter<T> converter) {
			sort(items, converter);
			for (int i = 0, j = items.length - 1; i < j; ++i, --j) {
				T tmp = items[i];
				items[i] = items[j];
				items[j] = tmp;
			}
			return items;
		}

		public static <T> T[] sort(T[] items, IntConverter<T> converter) {
			final int N = items.length;

			int[] values = new int[N];
			int minVal = Integer.MAX_VALUE;
			int maxVal = Integer.MIN_VALUE;
			for (int i = 0; i < items.length; ++i) {
				values[i] = converter.toInt(items[i]);
				minVal = Math.min(minVal, values[i]);
				maxVal = Math.max(maxVal, values[i]);
			}

			int capacity = maxVal - minVal + 1;
			ArrayList<ArrayList<T>> buckets = new ArrayList<>(capacity);
			for (int i = 0; i < capacity; ++i) {
				buckets.add(null);
			}

			for (int i = 0; i < items.length; ++i) {
				int bucketIndex = values[i] - minVal;
				ArrayList<T> lst = buckets.get(bucketIndex);
				if (lst == null) {
					lst = new ArrayList<>();
					buckets.set(bucketIndex, lst);
				}
				lst.add(items[i]);
			}

			int p = 0;
			for (ArrayList<T> lst : buckets) {
				if (lst == null) {
					continue;
				}
				for (T item : lst) {
					items[p++] = item;
				}
			}
			return items;
		}

		public static int[] sortIntsReversed(int[] items) {
			sortInts(items);
			for (int i = 0, j = items.length - 1; i < j; ++i, --j) {
				int tmp = items[i];
				items[i] = items[j];
				items[j] = tmp;
			}
			return items;
		}

		public static int[] sortInts(int[] items) {
			final int N = items.length;

			int minVal = Integer.MAX_VALUE;
			int maxVal = Integer.MIN_VALUE;
			for (int i = 0; i < items.length; ++i) {
				minVal = Math.min(minVal, items[i]);
				maxVal = Math.max(maxVal, items[i]);
			}

			int capacity = maxVal - minVal + 1;
			int[] counts = new int[capacity];

			for (int i = 0; i < N; ++i) {
				int bucketIndex = items[i] - minVal;
				++counts[bucketIndex];
			}

			int p = 0;
			for (int i = 0; i < capacity; ++i) {
				int origVal = i + minVal;
				for (int j = 0; j < counts[i]; ++j) {
					items[p++] = origVal;
				}
			}
			return items;
		}

		public static char[] sortCharsReversed(char[] items) {
			sortChars(items);
			for (int i = 0, j = items.length - 1; i < j; ++i, --j) {
				char tmp = items[i];
				items[i] = items[j];
				items[j] = tmp;
			}
			return items;
		}

		public static char[] sortChars(char[] items) {
			final int N = items.length;

			int minVal = Integer.MAX_VALUE;
			int maxVal = Integer.MIN_VALUE;
			for (int i = 0; i < items.length; ++i) {
				minVal = Math.min(minVal, items[i]);
				maxVal = Math.max(maxVal, items[i]);
			}

			int capacity = maxVal - minVal + 1;
			int[] counts = new int[capacity];

			for (int i = 0; i < N; ++i) {
				int bucketIndex = items[i] - minVal;
				++counts[bucketIndex];
			}

			int p = 0;
			for (int i = 0; i < capacity; ++i) {
				char origVal = (char) (i + minVal);
				for (int j = 0; j < counts[i]; ++j) {
					items[p++] = origVal;

				}
			}
			return items;
		}

		public static interface IntConverter<T> {
			public int toInt(T item);
		}
	}

	public static class RadixSort {
		/**
		 * Sorts the array "in-place", such that the returned array is the input array with elements in sorted order.
		 *
		 * Note that the following must hold true:
		 *     max(arr) - min(arr) <= Integer.MAX_VALUE
		 */
		public static int[] ints(int[] arr) {
			final int[] origArr = arr;
			final int N = arr.length;
			final int radix = Math.max(4, N);

			int minValue = Integer.MAX_VALUE;
			for (int x : arr) {
				minValue = Math.min(minValue, x);
			}

			int maxValue = Integer.MIN_VALUE;
			for (int i = 0; i < N; ++i) {
				arr[i] -= minValue;
				maxValue = Math.max(maxValue, arr[i]);
			}

			int exp = 1;
			int[] aux = new int[N];
			int[] idx = new int[N];
			int[] count = new int[radix];
			while (true) {
				Arrays.fill(count, 0);

				for (int i = 0; i < N; ++i) {
					idx[i] = arr[i] / exp % radix;
					++count[idx[i]];
				}

				for (int i = 1; i < radix; ++i) {
					count[i] += count[i - 1];
				}

				for (int i = N - 1; i >= 0; --i) {
					aux[--count[idx[i]]] = arr[i];
				}

				int[] tmp = aux;
				aux = arr;
				arr = tmp;

				if (exp > maxValue / radix) {
					break;
				}

				exp *= radix;
			}

			for (int i = 0; i < N; ++i) {
				arr[i] += minValue;
			}

			if (arr != origArr) {
				System.arraycopy(arr, 0, origArr, 0, N);
			}

			return origArr;
		}
	}

	/**
	 * Performs RadixSort with a fixed radix that is a power of two to take advantage of bitwise operators for efficiency.
	 * If all calls are on large arrays, e.g. at least 50 elements or so, prefer this over normal RadixSort.
	 */
	public static class FixedRadixSort {
		private static final int RADIX_BITS_FOR_INT = 11;
		private static final int MASK_FOR_INT = (1 << RADIX_BITS_FOR_INT) - 1;

		/**
		 * Sorts the array "in-place", such that the returned array is the input array with elements in sorted order.
		 *
		 * Note that the following must hold true:
		 *     max(arr) - min(arr) <= Integer.MAX_VALUE
		 */
		public static int[] ints(int[] arr) {
			final int[] origArr = arr;
			final int N = arr.length;

			int minValue = Integer.MAX_VALUE;
			for (int x : arr) {
				minValue = Math.min(minValue, x);
			}

			int maxValue = Integer.MIN_VALUE;
			for (int i = 0; i < N; ++i) {
				arr[i] -= minValue;
				maxValue = Math.max(maxValue, arr[i]);
			}

			int[] aux = new int[N];
			int[] count = new int[1 << RADIX_BITS_FOR_INT];
			int shift = 0;
			for (shift = 0; shift < Integer.SIZE && (maxValue >> shift) > 0; shift += RADIX_BITS_FOR_INT) {
				Arrays.fill(count, 0);

				for (int i = 0; i < N; ++i) {
					final int idx = (arr[i] >> shift) & MASK_FOR_INT;
					++count[idx];
				}

				for (int i = 1; i < count.length; ++i) {
					count[i] += count[i - 1];
				}

				for (int i = N - 1; i >= 0; --i) {
					final int idx = (arr[i] >> shift) & MASK_FOR_INT;
					aux[--count[idx]] = arr[i];
				}

				int[] tmp = aux;
				aux = arr;
				arr = tmp;
			}

			for (int i = 0; i < N; ++i) {
				arr[i] += minValue;
			}

			if (arr != origArr) {
				System.arraycopy(arr, 0, origArr, 0, N);
			}

			return origArr;
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

	/**
	 * Wrapper class for Knuth-Morris-Pratt (KMP) string search algorithm.
	 *   - Use `KMP.buildTable(...)` if only the failure table containing prefix matches is needed.
	 *   - Use `KMP.search(...)` to find all matching indices.
	 */
	public static class KMP {
		public int[] W;
		public int[] T;

		public KMP(String needle) {
			this(stringToIntArray(needle));
		}

		public KMP(char[] needle) {
			this(charArrayToIntArray(needle));
		}

		public KMP(int[] needle) {
			this.W = needle;
			this.T = buildTable(this.W);
		}

		public int[] search(String S) {
			return search(stringToIntArray(S));
		}

		public int[] search(char[] S) {
			return search(charArrayToIntArray(S));
		}

		public int[] search(int[] S) {
			final int N = S.length;
			final int M = W.length;

			if (N == 0) {
				if (M == 0) {
					return new int[1];
				} else {
					return new int[0];
				}
			}
			if (M == 0) {
				int[] arr = new int[N];
				for (int i = 0; i < N; ++i) {
					arr[i] = i;
				}
				return arr;
			}

			int j = 0;
			int k = 0;
			int[] found = new int[N];
			int p = 0;
			while (j < S.length) {
				if (W[k] == S[j]) {
					++j;
					++k;
					if (k == M) {
						found[p++] = j - k;
						k = T[k];
					}
				} else {
					k = T[k];
					if (k < 0) {
						++j;
						++k;
					}
				}
			}
			return Arrays.copyOf(found, p);
		}

		public static int[] buildTable(String W) {
			return buildTable(stringToIntArray(W));
		}

		public static int[] buildTable(char[] W) {
			return buildTable(charArrayToIntArray(W));
		}

		public static int[] buildTable(int[] W) {
			final int N = W.length;
			int[] T = new int[N + 1];
			T[0] = -1;
			int pos = 1;
			int cnd = 0;
			while (pos < N) {
				if (W[pos] == W[cnd]) {
					T[pos] = cnd;
				} else {
					T[pos] = cnd;
					while (cnd >= 0 && W[pos] != W[cnd]) {
						cnd = T[cnd];
					}
				}
				++pos;
				++cnd;
			}
			T[pos] = cnd;
			return T;
		}

		public static int[] search(String S, String W) {
			return new KMP(W).search(S);
		}

		public static int[] search(char[] S, char[] W) {
			return new KMP(W).search(S);
		}

		public static int[] search(int[] S, int[] W) {
			return new KMP(W).search(S);
		}

		private static int[] stringToIntArray(String str) {
			return charArrayToIntArray(str.toCharArray());
		}

		private static int[] charArrayToIntArray(char[] S) {
			final int N = S.length;

			int[] A = new int[N];
			for (int i = 0; i < N; ++i) {
				A[i] = S[i];
			}
			return A;
		}
	}

	/**
	 * Computes a polynomial hash to do fast equality checks.
	 * Excluding collisions, two subarrays are considered equal if their hashes by the `sub` function are equal.
	 * 
	 * NOTE:  If collisions are a concern, use SubHashMulti instead.
	 * NOTE:  It is NOT hack-resistant!
	 */
	public static class SubHash {
		private static final int P = 2147483647;
		private static final int K = 104723;
		private static final int OFFSET = 7;

		private static final int MAX_LEN = 2_000_002;

		private static int UPTO = 1;
		private static long[] POW26 = new long[MAX_LEN];
		private static long[] INV26 = new long[MAX_LEN];
		static {
			POW26[0] = 1;
			POW26[1] = K;
			INV26[0] = 1;
			INV26[1] = modInverse(K, P);
		}

		private static void loadPows(int upper) {
			for (int i = UPTO + 1; i <= upper; ++i) {
				POW26[i] = (POW26[i - 1] * POW26[1]) % P;
				INV26[i] = (INV26[i - 1] * INV26[1]) % P;
			}
			UPTO = Math.max(UPTO, upper);
		}

		private final int[] S;
		private final long[] H;

		public SubHash(String x) {
			loadPows(x.length());
			S = new int[x.length()];
			for (int i = 0; i < x.length(); ++i) {
				S[i] = x.charAt(i) - 'a' + OFFSET;
			}
			H = new long[S.length + 1];
			for (int i = 0; i < S.length; ++i) {
				H[i + 1] = (H[i] + (S[i] * POW26[i])) % P;
			}
		}

		public long sub(int i, int j) {
			long diff = (H[j] - H[i] + P) % P;
			long hash = (diff * INV26[i]) % P;
			return hash;
		}

		public int length() {
			return S.length;
		}

		/**
		 * Computes the value of (b ^ e) % mod.
		 */
		private static long modPow(long b, long e, long mod) {
			long p = b;
			long ans = 1;
			while (e > 0) {
				if ((e & 1) == 1) {
					ans = ans * p % mod;
				}
				p = p * p % mod;
				e >>= 1;
			}
			return ans;
		}

		/**
		 * Computes the modular inverse, such that: ak % MOD = 1, for some k.
		 * See this page for details:  http://rosettacode.org/wiki/Modular_inverse
		 */
		public static long modInverse(long a, long mod) {
			return modPow(a, mod - 2, mod);
		}
	}

	/**
	 * Computes multiple polynomial hashes to do fast equality checks.
	 * Excluding collisions, two subarrays are considered equal if their hashes by the `sub` function are equal.
	 * 
	 * To reduce collisions, add new distinct primes to `P` (modulo) and `K` (base), but it will cause a performance penalty.
	 * Note that the choices of P and K must not exceed 2147483647, e.g. it must fit within a signed 32-bit integer.
	 * 
	 * Some additional pairs to consider:
	 *   - P = 2122331213, K = 104717
	 *   - P = 2124749677, K = 104711
	 * 
	 * NOTE:  If only a single value in `P` is needed, consider using SubHash instead.
	 * NOTE:  It is NOT hack-resistant!
	 */
	public static class SubHashMulti {
		private static final int[] P = {2131131137, 2147483647};
		private static final int[] K = {104723, 104729};

		private static final int HASHES = P.length;

		private static final int MAX_LEN = 2_000_002;

		private static int UPTO = 1;
		private static long[][] POW = new long[HASHES][MAX_LEN];
		private static long[][] INV = new long[HASHES][MAX_LEN];
		static {
			for (int j = 0; j < HASHES; ++j) {
				POW[j][0] = 1;
				POW[j][1] = K[j];
				INV[j][0] = 1;
				INV[j][1] = modInverse(K[j], P[j]);
			}
		}

		private static void loadPows(int upper) {
			for (int j = 0; j < HASHES; ++j) {
				for (int i = UPTO + 1; i <= upper; ++i) {
					POW[j][i] = POW[j][i - 1] * POW[j][1] % P[j];
					INV[j][i] = INV[j][i - 1] * INV[j][1] % P[j];
				}
			}
			UPTO = Math.max(UPTO, upper);
		}

		private final long[] S;
		private final long[][] H;

		public SubHashMulti(long[] x) {
			loadPows(x.length);

			S = x;
			H = new long[HASHES][S.length + 1];
			for (int j = 0; j < HASHES; ++j) {
				for (int i = 0; i < S.length; ++i) {
					H[j][i + 1] = (H[j][i] + S[i] * POW[j][i]) % P[j];
				}
			}
		}

		public SubHashMulti(int[] x) {
			this(toLongArray(x));
		}

		public SubHashMulti(char[] x) {
			this(toLongArray(x));
		}

		public SubHashMulti(String x) {
			this(x.toCharArray());
		}

		public long[] sub(int loInclusive, int hiExclusive) {
			long[] hash = new long[HASHES];
			for (int j = 0; j < HASHES; ++j) {
				hash[j] = (H[j][hiExclusive] + P[j] - H[j][loInclusive]) * INV[j][loInclusive] % P[j];
			}
			return hash;
		}

		public int length() {
			return S.length;
		}

		private static long[] toLongArray(char[] x) {
			long[] arr = new long[x.length];
			for (int i = 0; i < x.length; ++i) {
				arr[i] = x[i];
			}
			return arr;
		}

		private static long[] toLongArray(int[] x) {
			long[] arr = new long[x.length];
			for (int i = 0; i < x.length; ++i) {
				arr[i] = x[i];
			}
			return arr;
		}

		/**
		 * Computes the value of (b ^ e) % mod.
		 */
		private static long modPow(long b, long e, long mod) {
			long p = b;
			long ans = 1;
			while (e > 0) {
				if ((e & 1) == 1) {
					ans = ans * p % mod;
				}
				p = p * p % mod;
				e >>= 1;
			}
			return ans;
		}

		/**
		 * Computes the modular inverse, such that: ak % MOD = 1, for some k.
		 * See this page for details:  http://rosettacode.org/wiki/Modular_inverse
		 */
		public static long modInverse(long a, long mod) {
			return modPow(a, mod - 2, mod);
		}
	}
}
