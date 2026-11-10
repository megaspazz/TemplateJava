import java.util.concurrent.*;

public class SubstringHashing {
	/**
	 * Computes a polynomial hash to do fast equality checks.
	 * Excluding collisions, two subarrays are considered equal if their hashes by the `get` function are equal.
	 *
	 * NOTE:  If collisions are a concern, use `SubHashTwo` or `SubHashMulti` instead.
	 * NOTE:  If collisions are not a concern, consider using `SubHashMersenne`, which is slightly faster but has higher collision probability.
	 */
	public static class SubHash {
		private static final long P = 3_013_622_539L;
		private static final long K = ThreadLocalRandom.current().nextLong(2, P - 1);

		private static final int MAX_LEN = 2_000_002;

		private static int UPTO = 0;
		private static long[] POW = new long[MAX_LEN];
		static {
			POW[0] = 1;
		}

		private static void loadPows(int upper) {
			for (int i = UPTO + 1; i <= upper; ++i) {
				POW[i] = (POW[i - 1] * K) % P;
			}
			UPTO = Math.max(UPTO, upper);
		}

		private final long[] S;
		private final long[] H;

		public SubHash(long[] x, boolean copyInput) {
			loadPows(x.length);

			this.S = x.clone();
			this.H = new long[S.length + 1];
			for (int i = 0; i < S.length; ++i) {
				H[i + 1] = add(H[i] * K % P, S[i]);
			}
		}

		public SubHash(long[] x) {
			this(x, true);
		}

		public SubHash(int[] x) {
			this(toLongArray(x), false);
		}

		public SubHash(char[] x) {
			this(toLongArray(x), false);
		}

		public SubHash(String x) {
			this(x.toCharArray());
		}

		public long get(int loInclusive, int hiExclusive) {
			return subtract(H[hiExclusive], H[loInclusive] * POW[hiExclusive - loInclusive] % P);
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

		private static long add(long a, long b) {
			long res = a + b;
			return res >= P ? res - P : res;
		}

		private static long subtract(long a, long b) {
			long res = a - b;
			return res < 0 ? res + P : res;
		}
	}

	/**
	 * Computes a polynomial hash to do fast equality checks.
	 * Excluding collisions, two subarrays are considered equal if their hashes by the `get` function are equal.
	 * 
	 * NOTE:  If collisions are a concern, use `SubHashMulti` instead.
	 */
	public static class SubHashTwo {
		private static final long P = 3_037_000_493L;
		private static final long KA = ThreadLocalRandom.current().nextLong(2, P - 1);
		private static final long KB = ThreadLocalRandom.current().nextLong(2, P - 1);

		private static final int MAX_LEN = 2_000_002;

		private static int UPTO = 0;
		private static long[] POW_A = new long[MAX_LEN];
		private static long[] POW_B = new long[MAX_LEN];
		static {
			POW_A[0] = POW_B[0] = 1;
		}

		private static void loadPows(int upper) {
			for (int i = UPTO + 1; i <= upper; ++i) {
				POW_A[i] = POW_A[i - 1] * KA % P;
				POW_B[i] = POW_B[i - 1] * KB % P;
			}
			UPTO = Math.max(UPTO, upper);
		}

		private final long[] S;
		private final long[] HA;
		private final long[] HB;

		public SubHashTwo(long[] x, boolean copyInput) {
			final int N = x.length;

			loadPows(N);

			this.S = copyInput ? x.clone() : x;
			this.HA = new long[N + 1];
			this.HB = new long[N + 1];
			for (int i = 0; i < N; ++i) {
				HA[i + 1] = add(HA[i] * KA % P, S[i]);
				HB[i + 1] = add(HB[i] * KB % P, S[i]);
			}
		}

		public SubHashTwo(long[] x) {
			this(x, true);
		}

		public SubHashTwo(int[] x) {
			this(toLongArray(x), false);
		}

		public SubHashTwo(char[] x) {
			this(toLongArray(x), false);
		}

		public SubHashTwo(String x) {
			this(toLongArray(x), false);
		}

		public long get(int loInclusive, int hiExclusive) {
			int shiftPow = hiExclusive - loInclusive;
			long hashA = subtract(HA[hiExclusive], HA[loInclusive] * POW_A[shiftPow] % P);
			long hashB = subtract(HB[hiExclusive], HB[loInclusive] * POW_B[shiftPow] % P);
			return hashA << 32 | hashB;
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

		private static long[] toLongArray(String x) {
			final int N = x.length();

			long[] arr = new long[N];
			for (int i = 0; i < N; ++i) {
				arr[i] = x.charAt(i);
			}
			return arr;
		}

		private static long add(long a, long b) {
			long res = a + b;
			return res >= P ? res - P : res;
		}

		private static long subtract(long a, long b) {
			long res = a - b;
			return res < 0 ? res + P : res;
		}
	}

	/**
	 * Computes multiple polynomial hashes to do fast equality checks.
	 * Excluding collisions, two subarrays are considered equal if their hashes by the `get` function are equal.
	 * 
	 * To reduce collisions, increase `HASHES`.
	 * For performance tuning, adjust `MAX_LEN` to be the maximum length of a subarray for the constraints of the problem.
	 * 
	 * NOTE:  If only a single value in `P` is needed, consider using `SubHash` instead.
	 * NOTE:  If only two values in `P` are needed, consider using `SubHashTwo` instead.
	 * NOTE:  Consider using `SubHashMersenneMulti`, which is faster, but has a higher chance of collision for the same number of hashes.
	 */
	public static class SubHashMulti {
		private static final int MAX_LEN = 2_000_002;

		private static final int HASHES = 4;

		private static final long P = 3_037_000_493L;

		private static final long[] K = new long[HASHES];
		static {
			for (int i = 0; i < HASHES; ++i) {
				K[i] = ThreadLocalRandom.current().nextLong(2, P - 1);
			}
		}

		private static int UPTO = 0;
		private static long[][] POW = new long[HASHES][MAX_LEN + 1];
		static {
			for (int j = 0; j < HASHES; ++j) {
				POW[j][0] = 1;
			}
		}

		private static void loadPows(int upper) {
			for (int j = 0; j < HASHES; ++j) {
				for (int i = UPTO + 1; i <= upper; ++i) {
					POW[j][i] = POW[j][i - 1] * K[j] % P;
				}
			}
			UPTO = Math.max(UPTO, upper);
		}

		private final long[] S;
		private final long[][] H;

		public SubHashMulti(long[] x, boolean copyInput) {
			final int N = x.length;

			loadPows(N);

			this.S = copyInput ? x.clone() : x;
			this.H = new long[HASHES][N + 1];
			for (int j = 0; j < HASHES; ++j) {
				for (int i = 0; i < N; ++i) {
					H[j][i + 1] = add(H[j][i] * K[j] % P, S[i]);
				}
			}
		}

		public SubHashMulti(long[] x) {
			this(x, true);
		}

		public SubHashMulti(int[] x) {
			this(toLongArray(x), false);
		}

		public SubHashMulti(char[] x) {
			this(toLongArray(x), false);
		}

		public SubHashMulti(String x) {
			this(toLongArray(x), false);
		}

		public long[] get(int loInclusive, int hiExclusive) {
			int shiftPow = hiExclusive - loInclusive;
			long[] hash = new long[HASHES];
			for (int j = 0; j < HASHES; ++j) {
				hash[j] = subtract(H[j][hiExclusive], H[j][loInclusive] * POW[j][shiftPow] % P);
			}
			return hash;
		}

		public int length() {
			return S.length;
		}

		private static long[] toLongArray(char[] x) {
			final int N = x.length;

			long[] arr = new long[N];
			for (int i = 0; i < x.length; ++i) {
				arr[i] = x[i];
			}
			return arr;
		}

		private static long[] toLongArray(int[] x) {
			final int N = x.length;

			long[] arr = new long[N];
			for (int i = 0; i < N; ++i) {
				arr[i] = x[i];
			}
			return arr;
		}

		private static long[] toLongArray(String x) {
			final int N = x.length();

			long[] arr = new long[N];
			for (int i = 0; i < N; ++i) {
				arr[i] = x.charAt(i);
			}
			return arr;
		}

		private static long add(long a, long b) {
			long res = a + b;
			return res >= P ? res - P : res;
		}

		private static long subtract(long a, long b) {
			long res = a - b;
			return res < 0 ? res + P : res;
		}
	}

	/**
	 * Computes a polynomial hash to do fast equality checks.
	 * Excluding collisions, two subarrays are considered equal if their hashes by the `get` function are equal.
	 * It is slightly faster than regular `SubHash`, since it uses a Mersenne prime for `P`, so we can implement modulo using bit operations.
	 */
	public static class SubHashMersenne {
		private static final int E = 31;
		private static final long P = (1L << E) - 1;
		private static final long K = ThreadLocalRandom.current().nextLong(2, P - 1);

		private static final int MAX_LEN = 2_000_002;

		private static int UPTO = 0;
		private static long[] POW = new long[MAX_LEN];
		static {
			POW[0] = 1;
		}

		private static void loadPows(int upper) {
			for (int i = UPTO + 1; i <= upper; ++i) {
				POW[i] = mersenneMod(POW[i - 1] * K);
			}
			UPTO = Math.max(UPTO, upper);
		}

		private final long[] S;
		private final long[] H;

		public SubHashMersenne(long[] x) {
			loadPows(x.length);

			this.S = x.clone();
			this.H = new long[S.length + 1];
			for (int i = 0; i < S.length; ++i) {
				H[i + 1] = mersenneMod(mersenneMod(H[i] * K) + S[i]);
			}
		}

		public SubHashMersenne(int[] x) {
			this(toLongArray(x));
		}

		public SubHashMersenne(char[] x) {
			this(toLongArray(x));
		}

		public SubHashMersenne(String x) {
			this(x.toCharArray());
		}

		public long get(int loInclusive, int hiExclusive) {
			return mersenneMod(H[hiExclusive] - mersenneMod(H[loInclusive] * POW[hiExclusive - loInclusive]) + P);
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

		private static long mersenneMod(long x) {
			long r = (x & P) + (x >> E);
			return r >= P ? r - P : r;
		}
	}

	/**
	 * Computes a polynomial hash to do fast equality checks.
	 * Excluding collisions, two subarrays are considered equal if their hashes by the `get` function are equal.
	 * It is slightly faster than regular `SubHash`, since it uses a Mersenne prime for `P`, so we can implement modulo using bit operations.
	 */
	public static class SubHashMersenneTwo {
		private static final int E = 31;
		private static final long P = (1L << E) - 1;
		private static final long K = ThreadLocalRandom.current().nextLong(2, P - 1);

		private static final int MAX_LEN = 2_000_002;

		private static int UPTO = 0;
		private static long[] POW = new long[MAX_LEN];
		static {
			POW[0] = 1;
		}

		private static void loadPows(int upper) {
			for (int i = UPTO + 1; i <= upper; ++i) {
				POW[i] = mersenneMod(POW[i - 1] * K);
			}
			UPTO = Math.max(UPTO, upper);
		}

		private final long[] S;
		private final long[] HA;
		private final long[] HB;

		public SubHashMersenneTwo(long[] x) {
			loadPows(x.length);

			this.S = x.clone();
			this.HA = new long[S.length + 1];
			this.HB = new long[S.length + 1];
			for (int i = 0; i < S.length; ++i) {
				HA[i + 1] = mersenneMod(mersenneMod(HA[i] * K) + S[i]);
				HB[i + 1] = mersenneMod(mersenneMod(HB[i] * K) + S[i]);
			}
		}

		public SubHashMersenneTwo(int[] x) {
			this(toLongArray(x));
		}

		public SubHashMersenneTwo(char[] x) {
			this(toLongArray(x));
		}

		public SubHashMersenneTwo(String x) {
			this(x.toCharArray());
		}

		public long get(int loInclusive, int hiExclusive) {
			int shiftPow = hiExclusive - loInclusive;
			long hashA = mersenneMod(HA[hiExclusive] - mersenneMod(HA[loInclusive] * POW[shiftPow]) + P);
			long hashB = mersenneMod(HB[hiExclusive] - mersenneMod(HB[loInclusive] * POW[shiftPow]) + P);
			return hashA << 32 | hashB;
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

		private static long mersenneMod(long x) {
			long r = (x & P) + (x >> E);
			return r >= P ? r - P : r;
		}
	}

	/**
	 * Computes a polynomial hash to do fast equality checks.
	 * Excluding collisions, two subarrays are considered equal if their hashes by the `get` function are equal.
	 * 
	 * It is slightly faster than `SubHashTwo`, since it uses a Mersenne prime for `P`, so we can implement modulo using bit operations.
	 * However, it does have a higher chance of collisions.
	 */
	public static class SubHashMersenneLong {
		private static final int E = 61;
		private static final long F = Long.SIZE - E;
		private static final long P = (1L << E) - 1;
		private static final long K = ThreadLocalRandom.current().nextLong(2, P - 1);

		private static final int MAX_LEN = 2_000_002;

		private static int UPTO = 0;
		private static long[] POW = new long[MAX_LEN];
		static {
			POW[0] = 1;
		}

		private static void loadPows(int upper) {
			for (int i = UPTO + 1; i <= upper; ++i) {
				POW[i] = mersenneMultiplyMod(POW[i - 1], K);
			}
			UPTO = Math.max(UPTO, upper);
		}

		private final long[] S;
		private final long[] H;

		public SubHashMersenneLong(long[] x) {
			loadPows(x.length);

			this.S = x.clone();
			this.H = new long[S.length + 1];
			for (int i = 0; i < S.length; ++i) {
				H[i + 1] = mersenneMod(mersenneMultiplyMod(H[i], K) + S[i]);
			}
		}

		public SubHashMersenneLong(int[] x) {
			this(toLongArray(x));
		}

		public SubHashMersenneLong(char[] x) {
			this(toLongArray(x));
		}

		public SubHashMersenneLong(String x) {
			this(x.toCharArray());
		}

		public long get(int loInclusive, int hiExclusive) {
			return mersenneMod(H[hiExclusive] - mersenneMultiplyMod(H[loInclusive], POW[hiExclusive - loInclusive]) + P);
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

		private static long mersenneMod(long x) {
			long r = (x & P) + (x >> E);
			return r >= P ? r - P : r;
		}

		private static long mersenneMultiplyMod(long a, long b) {
			long hi = Math.multiplyHigh(a, b);
			long lo = a * b;
			long r = (lo & P) + (lo >>> E | hi << F);
			return r >= P ? r - P : r;
		}
	}

	/**
	 * Computes multiple polynomial hashes to do fast equality checks.
	 * Excluding collisions, two subarrays are considered equal if their hashes by the `get` function are equal.
	 * 
	 * To reduce collisions, increase `HASHES`.
	 * For performance tuning, adjust `MAX_LEN` to be the maximum length of a subarray for the constraints of the problem.
	 */
	public static class SubHashMersenneMulti {
		private static final int MAX_LEN = 2_000_002;

		private static final int HASHES = 4;

		private static final int E = 31;
		private static final long P = (1L << E) - 1;

		private static final long[] K = new long[HASHES];
		static {
			for (int i = 0; i < HASHES; ++i) {
				K[i] = ThreadLocalRandom.current().nextLong(2, P - 1);
			}
		}

		private static int UPTO = 0;
		private static long[][] POW = new long[HASHES][MAX_LEN + 1];
		static {
			for (int j = 0; j < HASHES; ++j) {
				POW[j][0] = 1;
			}
		}

		private static void loadPows(int upper) {
			for (int j = 0; j < HASHES; ++j) {
				for (int i = UPTO + 1; i <= upper; ++i) {
					POW[j][i] = mersenneMod(POW[j][i - 1] * K[j]);
				}
			}
			UPTO = Math.max(UPTO, upper);
		}

		private final long[] S;
		private final long[][] H;

		public SubHashMersenneMulti(long[] x, boolean copyInput) {
			final int N = x.length;

			loadPows(N);

			this.S = copyInput ? x.clone() : x;
			this.H = new long[HASHES][N + 1];
			for (int j = 0; j < HASHES; ++j) {
				for (int i = 0; i < N; ++i) {
					H[j][i + 1] = add(mersenneMod(H[j][i] * K[j]), S[i]);
				}
			}
		}

		public SubHashMersenneMulti(long[] x) {
			this(x, true);
		}

		public SubHashMersenneMulti(int[] x) {
			this(toLongArray(x), false);
		}

		public SubHashMersenneMulti(char[] x) {
			this(toLongArray(x), false);
		}

		public SubHashMersenneMulti(String x) {
			this(toLongArray(x), false);
		}

		public long[] get(int loInclusive, int hiExclusive) {
			int shiftPow = hiExclusive - loInclusive;
			long[] hash = new long[HASHES];
			for (int j = 0; j < HASHES; ++j) {
				hash[j] = subtract(H[j][hiExclusive], mersenneMod(H[j][loInclusive] * POW[j][shiftPow]));
			}
			return hash;
		}

		public int length() {
			return S.length;
		}

		private static long[] toLongArray(char[] x) {
			final int N = x.length;

			long[] arr = new long[N];
			for (int i = 0; i < x.length; ++i) {
				arr[i] = x[i];
			}
			return arr;
		}

		private static long[] toLongArray(int[] x) {
			final int N = x.length;

			long[] arr = new long[N];
			for (int i = 0; i < N; ++i) {
				arr[i] = x[i];
			}
			return arr;
		}

		private static long[] toLongArray(String x) {
			final int N = x.length();

			long[] arr = new long[N];
			for (int i = 0; i < N; ++i) {
				arr[i] = x.charAt(i);
			}
			return arr;
		}

		private static long add(long a, long b) {
			long res = a + b;
			return res >= P ? res - P : res;
		}

		private static long subtract(long a, long b) {
			long res = a - b;
			return res < 0 ? res + P : res;
		}

		private static long mersenneMod(long x) {
			long r = (x & P) + (x >> E);
			return r >= P ? r - P : r;
		}
	}

	/**
	 * Computes multiple polynomial hashes to do fast equality checks.
	 * Excluding collisions, two subarrays are considered equal if their hashes by the `get` function are equal.
	 * 
	 * To reduce collisions, increase `HASHES`.
	 * For performance tuning, adjust `MAX_LEN` to be the maximum length of a subarray for the constraints of the problem.
	 */
	public static class SubHashMersenneLongMulti {
		private static final int MAX_LEN = 2_000_002;

		private static final int HASHES = 2;

		private static final int E = 61;
		private static final long F = Long.SIZE - E;
		private static final long P = (1L << E) - 1;
		private static final long[] K = new long[HASHES];

		static {
			for (int i = 0; i < HASHES; ++i) {
				K[i] = ThreadLocalRandom.current().nextLong(2, P - 1);
			}
		}

		private static int UPTO = 0;
		private static long[][] POW = new long[HASHES][MAX_LEN + 1];
		static {
			for (int j = 0; j < HASHES; ++j) {
				POW[j][0] = 1;
			}
		}

		private static void loadPows(int upper) {
			for (int j = 0; j < HASHES; ++j) {
				for (int i = UPTO + 1; i <= upper; ++i) {
					POW[j][i] = mersenneMultiplyMod(POW[j][i - 1], K[j]);
				}
			}
			UPTO = Math.max(UPTO, upper);
		}

		private final long[] S;
		private final long[][] H;

		public SubHashMersenneLongMulti(long[] x, boolean copyInput) {
			final int N = x.length;

			loadPows(N);

			this.S = copyInput ? x.clone() : x;
			this.H = new long[HASHES][N + 1];
			for (int j = 0; j < HASHES; ++j) {
				for (int i = 0; i < N; ++i) {
					H[j][i + 1] = add(mersenneMultiplyMod(H[j][i], K[j]), S[i]);
				}
			}
		}

		public SubHashMersenneLongMulti(long[] x) {
			this(x, true);
		}

		public SubHashMersenneLongMulti(int[] x) {
			this(toLongArray(x), false);
		}

		public SubHashMersenneLongMulti(char[] x) {
			this(toLongArray(x), false);
		}

		public SubHashMersenneLongMulti(String x) {
			this(toLongArray(x), false);
		}

		public long[] get(int loInclusive, int hiExclusive) {
			int shiftPow = hiExclusive - loInclusive;
			long[] hash = new long[HASHES];
			for (int j = 0; j < HASHES; ++j) {
				hash[j] = subtract(H[j][hiExclusive], mersenneMultiplyMod(H[j][loInclusive], POW[j][shiftPow]));
			}
			return hash;
		}

		public int length() {
			return S.length;
		}

		private static long[] toLongArray(char[] x) {
			final int N = x.length;

			long[] arr = new long[N];
			for (int i = 0; i < x.length; ++i) {
				arr[i] = x[i];
			}
			return arr;
		}

		private static long[] toLongArray(int[] x) {
			final int N = x.length;

			long[] arr = new long[N];
			for (int i = 0; i < N; ++i) {
				arr[i] = x[i];
			}
			return arr;
		}

		private static long[] toLongArray(String x) {
			final int N = x.length();

			long[] arr = new long[N];
			for (int i = 0; i < N; ++i) {
				arr[i] = x.charAt(i);
			}
			return arr;
		}

		private static long add(long a, long b) {
			long res = a + b;
			return res >= P ? res - P : res;
		}

		private static long subtract(long a, long b) {
			long res = a - b;
			return res < 0 ? res + P : res;
		}

		private static long mersenneMultiplyMod(long a, long b) {
			long hi = Math.multiplyHigh(a, b);
			long lo = a * b;
			long r = (lo & P) + (lo >>> E | hi << F);
			return r >= P ? r - P : r;
		}
	}
}
