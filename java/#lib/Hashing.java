public class Hashing {
	public static class IntKey {
		private static final int ADD_MIX = mulberry32((int)System.nanoTime());
		
		public int value;

		public IntKey(int value) {
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
	}
}
