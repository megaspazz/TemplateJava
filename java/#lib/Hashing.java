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
			IntKey other = (IntKey)obj;
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

		private static final int ADD_MIX = mulberry32((int)System.nanoTime());

		private static final int CACHE_MIN = -256;
		private static final int CACHE_MAX = 256;
		private static final IntKey[] CACHE = new IntKey[CACHE_MAX - CACHE_MIN];
		static {
			for (int i = CACHE_MIN; i < CACHE_MAX; ++i) {
				CACHE[i - CACHE_MIN] = new IntKey(i);
			}
		}
	}
}
