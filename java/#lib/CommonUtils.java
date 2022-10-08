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
}
