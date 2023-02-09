public class Hashing {
	public static class IntKey {
		public int value;

		private IntKey(int value) {
			this.value = value;
		}

		@Override
		public int hashCode() {
			return mulberry32(value + ADD_MIX);
		}

		@Override
		public boolean equals(Object obj) {
			IntKey other = (IntKey) obj;
			return value == other.value;
		}

		@Override
		public String toString() {
			return Integer.toString(value);
		}

		private static final int mulberry32(int x) {
			int z = (x + 0x6D2B79F5);
			z = (z ^ (z >>> 15)) * (z | 1);
			z ^= z + (z ^ (z >>> 7)) * (z | 61);
			return z ^ (z >>> 14);
		}

		public static IntKey of(int value) {
			if (CACHE_MIN <= value && value < CACHE_MAX) {
				return CACHE[value - CACHE_MIN];
			}
			return new IntKey(value);
		}

		private static final int ADD_MIX = mulberry32((int) System.nanoTime());

		private static final int CACHE_MIN = -256;
		private static final int CACHE_MAX = 256;
		private static final IntKey[] CACHE = new IntKey[CACHE_MAX - CACHE_MIN];
		static {
			for (int i = CACHE_MIN; i < CACHE_MAX; ++i) {
				CACHE[i - CACHE_MIN] = new IntKey(i);
			}
		}
	}

	public static class LongKey {
		public long value;

		private LongKey(long value) {
			this.value = value;
		}

		@Override
		public int hashCode() {
			return Long.hashCode(splitmix64(value + ADD_MIX));
		}

		@Override
		public boolean equals(Object obj) {
			LongKey other = (LongKey) obj;
			return value == other.value;
		}

		@Override
		public String toString() {
			return Long.toString(value);
		}

		private static final long splitmix64(long x) {
			long z = x + 0x9E3779B97F4A7C15L;
			z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
			z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
			return z ^ (z >>> 31);
		}

		public static LongKey of(long value) {
			if (CACHE_MIN <= value && value < CACHE_MAX) {
				return CACHE[(int) (value - CACHE_MIN)];
			}
			return new LongKey(value);
		}

		private static final long ADD_MIX = splitmix64(System.nanoTime());

		private static final int CACHE_MIN = -256;
		private static final int CACHE_MAX = 256;
		private static final LongKey[] CACHE = new LongKey[CACHE_MAX - CACHE_MIN];
		static {
			for (int i = CACHE_MIN; i < CACHE_MAX; ++i) {
				CACHE[i - CACHE_MIN] = new LongKey(i);
			}
		}
	}
}
