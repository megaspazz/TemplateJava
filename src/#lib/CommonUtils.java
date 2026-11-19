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

		public boolean isEmpty() {
			return size() == 0;
		}

		public int uniqueCount() {
			return uniq;
		}

		public int size() {
			return uniq;
		}

		public long totalCount() {
			return total;
		}

		public void increment(int k) {
			increment(k, 1);
		}

		public void decrement(int k) {
			increment(k, -1);
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

		public long getCount(int k) {
			return count[getIndex(k)];
		}

		public void clear() {
			uniq = 0;
			total = 0;

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
	public static class ReusableFixedIntToIntMap {
		private final int offset;
		private final int N;
		private final int[] A;
		private final boolean[] seen;
		private final int[] insertions;
		private int insertionsIdx = 0;

		public ReusableFixedIntToIntMap(int loInclusive, int hiExclusive) {
			this.offset = loInclusive;
			this.N = hiExclusive - loInclusive;
			this.A = new int[N];
			this.seen = new boolean[N];
			this.insertions = new int[N];
		}

		public ReusableFixedIntToIntMap(int N) {
			this(0, N);
		}

		public void put(int k, int v) {
			int i = getIndex(k);
			A[i] = v;
			seen[i] = true;
			if (insertionsIdx < A.length) {
				insertions[insertionsIdx++] = i;
			}
		}

		public boolean containsKey(int k) {
			int i = getIndex(k);
			return seen[i];
		}

		public int get(int k) {
			int i = getIndex(k);
			return A[i];
		}

		public int getOrDefault(int k, int defaultValue) {
			return containsKey(k) ? get(k) : defaultValue;
		}

		public void clear() {
			if (insertionsIdx >= N) {
				Arrays.fill(seen, false);
			} else {
				for (int i = 0; i < insertionsIdx; ++i) {
					seen[insertions[i]] = false;
				}
			}
			insertionsIdx = 0;
		}

		private int getIndex(int k) {
			return k - offset;
		}
	}

	/**
	 * Flat generic int-to-any map designed to be re-used across multiple test cases.
	 * For example, initialize it as a static variable, or pass it around everywhere it's needed.
	 * Instead of instantiating a new one, to re-use an existing one, just call the `clear()` method.
	 * 
	 * TODO:  Implement offset from ReusableFixedIntToIntMap.
	 */
	public static class ReusableFixedIntMap<T> {
		private final T[] A;
		private final boolean[] seen;
		private final int[] insertions;
		private int insertionsIdx = 0;

		public ReusableFixedIntMap(int N, T[] template) {
			this.A = Arrays.copyOf(template, N);
			Arrays.fill(A, null);

			this.seen = new boolean[N];
			this.insertions = new int[N];
		}

		public void put(int k, T v) {
			A[k] = v;
			seen[k] = true;

			if (insertionsIdx < A.length) {
				insertions[insertionsIdx++] = k;
			}
		}

		public boolean containsKey(int k) {
			return seen[k];
		}

		public T get(int k) {
			return getOrDefault(k, null);
		}

		public T getOrDefault(int k, T defaultValue) {
			return containsKey(k) ? A[k] : defaultValue;
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

	/**
	 * Helper methods for packing two int values into a single long value.
	 * Useful for efficiently representing a class that would have two int values, particularly in sorting.
	 * Copy the functions without the class definition for less verbose code.
	 */
	public static class IntPair {
		public static long make(int first, int second) {
			return ((long) first << 32) | second;
		}

		public static int first(long pair) {
			return (int) (pair >>> 32);
		}

		public static int second(long pair) {
			return (int) (pair & ~0);
		}
	}
}