import java.util.*;

public class CommonUtils {
	/**
	 * Counts the frequency of objects.
	 * Change to extend TreeMap instead, if ordering of objects is required.
	 * 
	 * NOTE:  If `total` is needed, only use `increment(...)` for updates.
	 */
	public static class CountMapInt<T> extends HashMap<T, Integer> {
		private static final long serialVersionUID = -1501598139835601959L;

		public int total;

		public int getCount(T k) {
			return getOrDefault(k, 0);
		}

		public void increment(T k, int v) {
			total += v;
			int next = getCount(k) + v;
			if (next == 0) {
				remove(k);
			} else {
				put(k, next);
			}
		}

		public static <T> CountMapInt<T> fromArray(T[] A) {
			CountMapInt<T> cm = new CountMapInt<>();
			for (T x : A) {
				cm.increment(x, 1);
			}
			return cm;
		}
	}

	/**
	 * Counts the frequency of objects.
	 * Change to extend TreeMap instead, if ordering of objects is required.
	 * 
	 * NOTE:  If `total` is needed, only use `increment(...)` for updates.
	 */
	public static class CountMapLong<T> extends HashMap<T, Long> {
		private static final long serialVersionUID = -9079906779955923767L;

		public long total;

		public long getCount(T k) {
			return getOrDefault(k, 0L);
		}

		public void increment(T k, long v) {
			total += v;
			long next = getCount(k) + v;
			if (next == 0) {
				remove(k);
			} else {
				put(k, next);
			}
		}

		public static <T> CountMapLong<T> fromArray(T[] A) {
			CountMapLong<T> cm = new CountMapLong<>();
			for (T x : A) {
				cm.increment(x, 1);
			}
			return cm;
		}
	}

	/**
	 * Keeps track of the top two elements inserted.
	 * If `first` and `second` are initialized, `count` will only include elements that were strictly greater than the initial `second` value.
	 */
	public static class TopTwoInt {
		public int count;
		public int first;
		public int second;

		public TopTwoInt() {
			this(Integer.MIN_VALUE);
		}

		public TopTwoInt(int init) {
			this(init, init);
		}

		public TopTwoInt(int first, int second) {
			this.first = Math.max(first, second);
			this.second = Math.min(first, second);
		}

		public void add(int x) {
			if (x < second) {
				return;
			}

			if (x > first) {
				second = first;
				first = x;
			} else {
				second = x;
			}
			count = Math.min(2, count + 1);
		}
		
		@Override
		public String toString() {
			return "[" + first + ", " + second + "]";
		}
	}

	/**
	 * Keeps track of the top two elements inserted.
	 * If `first` and `second` are initialized, `count` will only include elements that were strictly greater than the initial `second` value.
	 */
	public static class TopTwoLong {
		public int count;
		public long first;
		public long second;

		public TopTwoLong() {
			this(Long.MIN_VALUE);
		}

		public TopTwoLong(long init) {
			this(init, init);
		}

		public TopTwoLong(long first, long second) {
			this.first = Math.max(first, second);
			this.second = Math.min(first, second);
		}

		public void add(long x) {
			if (x < second) {
				return;
			}

			if (x > first) {
				second = first;
				first = x;
			} else {
				second = x;
			}
			count = Math.min(2, count + 1);
		}
	}

	/**
	 * Keeps track of the bottom two elements inserted.
	 * If `first` and `second` are initialized, `count` will only include elements that were strictly less than the initial `second` value.
	 */
	public static class BottomTwoLong {
		public int count;
		public long first;
		public long second;

		public BottomTwoLong() {
			this(Long.MAX_VALUE);
		}

		public BottomTwoLong(long init) {
			this(init, init);
		}

		public BottomTwoLong(long first, long second) {
			this.first = Math.max(first, second);
			this.second = Math.min(first, second);
		}
		
		public BottomTwoLong(long[] arr) {
			this();
			for (long x : arr) {
				add(x);
			}
		}

		public void add(long x) {
			if (x > second) {
				return;
			}

			if (x < first) {
				second = first;
				first = x;
			} else {
				second = x;
			}
			count = Math.min(2, count + 1);
		}
	}
	
	public static class IntMultiSet {
		public final int offset;
		public final long[] count;
		
		private int uniq;
		private long total;
		
		public IntMultiSet(int hiExclusive) {
			this(0, hiExclusive);
		}
		
		public IntMultiSet(int loInclusive, int hiExclusive) {
			this.offset = loInclusive;
			this.count = new long[hiExclusive - loInclusive];
		}
		
		public int uniqueCount() {
			return uniq;
		}
		
		public long size() {
			return total;
		}
		
		public void increment(int k) {
			increment(k, 1);
		}
		
		public void decrement(int k) {
			increment(k, -1);
		}
		
		public void increment(int k, long v) {
			final int key = k - offset;
			
			if (count[key] == 0) {
				++uniq;
			}
			
			count[key] += v;
			total += v;
			
			if (count[key] == 0) {
				--uniq;
			}
		}
	}
}
