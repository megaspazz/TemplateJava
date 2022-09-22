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
	private static long C(int n, int k) {
		return F[n] * FI[k] % MOD * FI[n - k] % MOD;
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

	/**
	 * Computes the GCD (greatest common denominator) between two numbers.
	 * NOTE:  Only works with non-negative numbers!
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
}
