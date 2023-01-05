public class MathUtils {
	private static final int MOD = 1_000_000_007;

	/**
	 * Pre-compute factorial function and mod-inverse of factorial function in linear time.
	 */
	private static final long[] F = new long[400001];
	private static final long[] INV = new long[F.length];
	private static final long[] FI = new long[F.length];
	static {
		INV[1] = 1;
		for (int i = 2; i < INV.length; ++i) {
			INV[i] = MOD - (MOD / i) * INV[MOD % i] % MOD;
		}

		F[0] = FI[0] = 1;
		for (int i = 1; i < F.length; ++i) {
			F[i] = (i * F[i - 1]) % MOD;
			FI[i] = (INV[i] * FI[i - 1]) % MOD;
		}
	}

	/**
	 * Computes the modulo result of (n choose k).
	 * Requires pre-computing the factorial function and mod-inverse of factorial function.
	 */
	public static long C(int n, int k) {
		return F[n] * FI[k] % MOD * FI[n - k] % MOD;
	}

	/**
	 * Computes the GCD (greatest common denominator) between all the numbers.
	 * NOTE: Only works with non-negative numbers!
	 */
	public static int gcd(int... arr) {
		int g = 0;
		for (int x : arr) {
			g = gcd(g, x);
		}
		return g;
	}

	/**
	 * Computes the GCD (greatest common denominator) between two numbers.
	 * NOTE: Only works with non-negative numbers!
	 */
	public static int gcd(int a, int b) {
		if (a < b) {
			return gcd(b, a);
		}

		if (b == 0) {
			return a;
		}

		int r = a % b;
		if (r == 0) {
			return b;
		}

		return gcd(b, r);
	}

	/**
	 * Computes the LCM (least common multiple) between all the numbers.
	 * NOTE:  Only works with positive numbers!
	 */
	public static int lcm(int... arr) {
		int ans = 1;
		for (int x : arr) {
			ans *= x / gcd(ans, x);
		}
		return ans;
	}

	/**
	 * Computes the GCD (greatest common denominator) between all the numbers.
	 * NOTE: Only works with non-negative numbers!
	 */
	public static long gcd(long... arr) {
		long g = 0;
		for (long x : arr) {
			g = gcd(g, x);
		}
		return g;
	}

	/**
	 * Computes the GCD (greatest common denominator) between two numbers.
	 * NOTE: Only works with non-negative numbers!
	 */
	public static long gcd(long a, long b) {
		if (a < b) {
			return gcd(b, a);
		}

		if (b == 0) {
			return a;
		}

		long r = a % b;
		if (r == 0) {
			return b;
		}

		return gcd(b, r);
	}

	/**
	 * Computes the LCM (least common multiple) between all the numbers.
	 * NOTE:  Only works with positive numbers!
	 */
	public static long lcm(long... arr) {
		long ans = arr[0];
		for (long x : arr) {
			ans *= x / gcd(ans, x);
		}
		return ans;
	}

	/**
	 * Computes arithmetic operations modulo some constant MOD, should be a prime.
	 * It is guaranteed to not overflow, as long as all operands are within the range [0, MOD).
	 * 
	 * NOTE: It's recommended to copy the code in this class and inline it within your other functions,
	 *       since it's cumbersome to write `LongModMath.multiply(a, LongModMath.add(b, c)`,
	 *       compared to simply writing `multiply(a, add(b, c))`.
	 * 
	 * NOTE:  If you need multiple instances of this class, e.g. one modulo for computing an answer,
	 *        and another modulo for polynomial string hashing, copy and rename a copy of the static class.
	 *        In each copy, create a class-local static variable named `MOD` for the different modulos.
	 *        It's recommended to rename the class something shorter to make the calling code less cumbersome.
	 *        If these classes are not static, there will be performance issues of unknown causes.
	 */
	public static class LongModMath {
		private static final long RAW_MULTIPLY_MAX = 3037000499L;

		private static final int CHUNK_SIZE = Long.SIZE - Long.numberOfLeadingZeros(Long.MAX_VALUE / (MOD + 1) + 1) - 1;
		private static final long CHUNK_MASK = (1L << CHUNK_SIZE) - 1;

		@SuppressWarnings("unused")
		public static long multiply(long a, long b) {
			if (MOD <= RAW_MULTIPLY_MAX) {
				return a * b % MOD;
			}
			return multiplyInternal(a, b);
		}

		public static long multiply(long... arr) {
			long ans = 1;
			for (long x : arr) {
				ans = multiply(ans, x);
			}
			return ans;
		}

		public static long add(long a, long b) {
			long ans = a + b;
			if (ans >= MOD) {
				ans -= MOD;
			}
			return ans;
		}

		public static long add(long... arr) {
			long ans = 0;
			for (long x : arr) {
				ans = add(ans, x);
			}
			return ans;
		}

		public static long subtract(long a, long b) {
			return add(a, MOD - b);
		}

		/**
		 * Computes the value of (b ^ e) % MOD.
		 */
		public static long modPow(long b, long e) {
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
		public static long modInverse(long a) {
			return modPow(a, MOD - 2);
		}

		private static long multiplyInternal(long a, long b) {
			if (a > b) {
				return multiplyInternal(b, a);
			}
			if (a == 0) {
				return 0;
			}

			long ans = 0;
			while (a > 0) {
				long mask = a & CHUNK_MASK;
				if (mask > 0) {
					ans = add(ans, (mask * b) % MOD);
				}
				b = (b << CHUNK_SIZE) % MOD;
				a >>= CHUNK_SIZE;
			}
			return ans;
		}
	}

	/**
	 * Computes all the primes up to the specified number.
	 */
	public static int[] primesTo(int n) {
		ArrayList<Integer> primes = new ArrayList<Integer>();

		boolean[] prime = new boolean[n + 1];
		Arrays.fill(prime, 2, n + 1, true);
		for (int i = 2; i < prime.length; i++) {
			if (prime[i]) {
				primes.add(i);
				if ((long) i * i <= n) {
					for (int j = i * i; j < prime.length; j += i) {
						prime[j] = false;
					}
				}
			}
		}

		int[] ans = new int[primes.size()];
		for (int i = 0; i < ans.length; ++i) {
			ans[i] = primes.get(i);
		}
		return ans;
	}
	
	/**
	 * Represent fractions in the form: n / d.
	 * 
	 * NOTE:  The value of n * d must not exceed the limits of what a 64-bit integer can hold.
	 * NOTE:  Do not use the private constructor; instead use Fraction.canonical(...).  Otherwise, comparisons won't work! 
	 */
	private static class Fraction implements Comparable<Fraction> {
		public long n, d;

		private Fraction(long n, long d) {
			this.n = n;
			this.d = d;
		}

		@Override
		public int compareTo(Fraction obj) {
			Fraction f = (Fraction) obj;
			return Long.compare(n * f.d, f.n * d);
		}

		public static Fraction min(Fraction a, Fraction b) {
			if (a.compareTo(b) < 0) {
				return a;
			} else {
				return b;
			}
		}

		public static Fraction max(Fraction a, Fraction b) {
			if (a.compareTo(b) > 0) {
				return a;
			} else {
				return b;
			}
		}

		public static Fraction canonical(int n, int d) {
			if (n == 0) {
				return Fraction.ZERO;
			}

			if (d < 0) {
				n = -n;
				d = -d;
			}
			return new Fraction(n, d);
		}

		public static final Fraction ZERO = new Fraction(0, 1);
	}

	private static int[] primeDivisorsTo(int n) {
		int[] div = new int[n + 1];
		div[1] = 1;
		for (int i = 2; i < div.length; i++) {
			if (div[i] == 0) {
				div[i] = i;
				if (i < div.length / i) {
					for (int j = i * i; j < div.length; j += i) {
						div[j] = i;
					}
				}
			}
		}
		return div;
	}
}
