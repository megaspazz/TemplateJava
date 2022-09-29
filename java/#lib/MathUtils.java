public class Combinatorics {
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
	 * NOTE:  Only works with non-negative numbers!
	 */
	public static int lcm(int... arr) {
		int ans = 1;
		for (int x : arr) {
			ans = ans * x / gcd(ans, x);
		}
		return ans;
	}

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
	 * NOTE:  Only works with non-negative numbers!
	 */
	public static long lcm(long... arr) {
		long ans = arr[0];
		for (long x : arr) {
			ans = ans * x / gcd(ans, x);
		}
		return ans;
	}
	
	/**
	 * Computes arithmetic operations modulo some constant MOD, should be a prime.
	 * It is guaranteed to not overflow, as long as the following requirements are satisfied.
	 *   - MOD must be less than or equal to Long.MAX_VALUE / 2^CHUNK_SIZE.
	 *   - All operands must be within the range [0, MOD).
	 * 
	 * NOTE: It's recommended to copy the code in this class and inline it within your other functions,
	 *       since it's cumbersome to write `LongModMath.multiply(a, LongModMath.add(b, c)`,
	 *       compared to simply writing `multiply(a, add(b, c))`.
	 */
	public static class LongModMath {
		private static final long MOD = 1111111111111111111L;

		private static final int CHUNK_SIZE = 3;
		private static final int CHUNK_MASK = (1 << CHUNK_SIZE) - 1;

		private static long multiply(long a, long b) {
			if (a > b) {
				return multiply(b, a);
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
		
		private static long add(long a, long b) {
			long ans = a + b;
			if (ans >= MOD) {
				ans -= MOD;
			}
			return ans;
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
}
