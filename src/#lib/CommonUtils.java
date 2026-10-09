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

		public int decrement(T k) {
			return decrement(k, 1);
		}

		public int decrement(T k, int v) {
			return increment(k, -v);
		}

		public int increment(T k) {
			return increment(k, 1);
		}

		public int increment(T k, int v) {
			total += v;
			int next = getCount(k) + v;
			if (next == 0) {
				remove(k);
			} else {
				put(k, next);
			}
			return next;
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

		public long decrement(T k) {
			return decrement(k, 1);
		}

		public long decrement(T k, long v) {
			return increment(k, -v);
		}

		public long increment(T k) {
			return increment(k, 1);
		}

		public long increment(T k, long v) {
			total += v;
			long next = getCount(k) + v;
			if (next == 0) {
				remove(k);
			} else {
				put(k, next);
			}
			return next;
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

	/**
	 * Flat int-to-long multi-set designed to be re-used across multiple test cases.
	 * For example, initialize it as a static variable, or pass it around everywhere it's needed.
	 * Instead of instantiating a new one, to re-use an existing one, just call the `clear()` method.
	 */
	public static class ReusableIntMultiSet {
		private final int N;

		public final int offset;
		public final long[] count;

		private int uniq;
		private long total;

		private final int[] updates;
		private int updatesIdx;

		public ReusableIntMultiSet(int hiExclusive) {
			this(0, hiExclusive);
		}

		public ReusableIntMultiSet(int loInclusive, int hiExclusive) {
			this.N = hiExclusive - loInclusive;
			this.offset = loInclusive;
			this.count = new long[N];
			this.updates = new int[N];
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

		public long getCount(int k) {
			return count[getIndex(k)];
		}

		public void increment(int k, long v) {
			final int key = getIndex(k);

			if (count[key] == 0) {
				++uniq;
			}

			count[key] += v;
			total += v;

			if (count[key] == 0) {
				--uniq;
			} else if (updatesIdx < N) {
				updates[updatesIdx++] = key;
			}
		}

		public void clear() {
			if (updatesIdx >= N) {
				Arrays.fill(count, 0);
			} else {
				for (int i = 0; i < updatesIdx; ++i) {
					count[updates[i]] = 0;
				}
			}
			updatesIdx = 0;
		}

		private int getIndex(int k) {
			return k - offset;
		}
	}

	/**
	 * Flat int-to-int map designed to be re-used across multiple test cases.
	 * For example, initialize it as a static variable, or pass it around everywhere it's needed.
	 * Instead of instantiating a new one, to re-use an existing one, just call the `clear()` method.
	 */
	public static class ResuableFixedIntToIntMap {
		private final int[] A;
		private final boolean[] seen;
		private final int[] insertions;
		private int insertionsIdx = 0;

		public ResuableFixedIntToIntMap(int N) {
			this.A = new int[N];
			this.seen = new boolean[N];
			this.insertions = new int[N];
		}

		public void put(int k, int v) {
			A[k] = v;
			seen[k] = true;

			if (insertionsIdx < A.length) {
				insertions[insertionsIdx++] = k;
			}
		}

		public boolean containsKey(int k) {
			return seen[k];
		}

		public int get(int k) {
			return A[k];
		}

		public int getOrDefault(int k, int defaultValue) {
			return containsKey(k) ? get(k) : defaultValue;
		}

		public void clear() {
			if (insertionsIdx >= A.length) {
				Arrays.fill(seen, false);
			} else {
				for (int i = 0; i < insertionsIdx; ++i) {
					seen[insertions[i]] = false;
				}
			}
			insertionsIdx = 0;
		}
	}
}
