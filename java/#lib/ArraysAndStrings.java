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
	 * To reduce collisions, increase `P`, but it comes with a performance penalty.
	 * Some reasonable primes for P:
	 *   - 2147483647
	 *   - 888888877777777
	 *   - 1111111111111111111
	 */
	public static class SubHash {
		private static final long P = 2147483647;
		private static final long K = 104729;

		private static final int MAXLEN = 2_000_002;

		private static final LongModMath LMM = new LongModMath(P);

		private static int UPTO = 1;
		private static long[] POW = new long[MAXLEN];
		private static long[] INV = new long[MAXLEN];
		static {
			POW[0] = 1;
			POW[1] = K;
			INV[0] = 1;
			INV[1] = LMM.modInverse(K);
		}

		private static void loadPows(int upper) {
			for (int i = UPTO + 1; i <= upper; ++i) {
				POW[i] = LMM.multiply(POW[i - 1], POW[1]);
				INV[i] = LMM.multiply(INV[i - 1], INV[1]);
			}
			UPTO = Math.max(UPTO, upper);
		}

		private final long[] S;
		private final long[] H;

		public SubHash(long[] x) {
			loadPows(x.length);

			S = x;
			H = new long[S.length + 1];
			for (int i = 0; i < S.length; ++i) {
				H[i + 1] = LMM.add(H[i], LMM.multiply(S[i], POW[i]));
			}
		}

		public SubHash(int[] x) {
			this(toLongArray(x));
		}

		public SubHash(char[] x) {
			this(toLongArray(x));
		}

		public SubHash(String x) {
			this(x.toCharArray());
		}

		public long sub(int loInclusive, int j) {
			long diff = LMM.subtract(H[j], H[loInclusive]);
			long hash = LMM.multiply(diff, INV[loInclusive]);
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

		// Included as nested class to make copy-pasting easier.
		// Refer to MathUtils.LongModMath for the canonical version.
		public static class LongModMath {
			private static final long RAW_MULTIPLY_MAX = 3037000499L;

			private final long mod;

			private final int chunkSize;
			private final long chunkMask;

			public LongModMath(long mod) {
				this.mod = mod;
				this.chunkSize = Long.SIZE - Long.numberOfLeadingZeros(Long.MAX_VALUE / (mod + 1) + 1) - 1;
				this.chunkMask = (1L << chunkSize) - 1;
			}

			public long multiply(long a, long b) {
				if (mod <= RAW_MULTIPLY_MAX) {
					return a * b % mod;
				}
				return multiplyInternal(a, b);
			}

			public long multiply(long... arr) {
				long ans = 1;
				for (long x : arr) {
					ans = multiply(ans, x);
				}
				return ans;
			}

			public long add(long a, long b) {
				long ans = a + b;
				if (ans >= mod) {
					ans -= mod;
				}
				return ans;
			}

			public long add(long... arr) {
				long ans = 0;
				for (long x : arr) {
					ans = add(ans, x);
				}
				return ans;
			}

			public long subtract(long a, long b) {
				return add(a, mod - b);
			}

			/**
			 * Computes the value of (b ^ e) % MOD.
			 */
			public long modPow(long b, long e) {
				long p = b;
				long ans = 1;
				while (e > 0) {
					if ((e & 1) == 1) {
						ans = multiply(ans, p);
					}
					p = multiply(p, p);
					e >>= 1;
				}
				return ans;
			}

			/**
			 * Computes the modular inverse, such that: ak % MOD = 1, for some k.
			 * See this page for details:  http://rosettacode.org/wiki/Modular_inverse
			 */
			public long modInverse(long a) {
				return modPow(a, mod - 2);
			}

			private long multiplyInternal(long a, long b) {
				if (a > b) {
					return multiplyInternal(b, a);
				}
				if (a == 0) {
					return 0;
				}

				long ans = 0;
				while (a > 0) {
					long mask = a & chunkMask;
					if (mask > 0) {
						ans = add(ans, (mask * b) % mod);
					}
					b = (b << chunkSize) % mod;
					a >>= chunkSize;
				}
				return ans;
			}
		}
	}
}
