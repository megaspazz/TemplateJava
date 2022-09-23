import java.util.*;

public class CommonUtils {
	/**
	 * Counts the frequency of objects.
	 * Change to extend TreeMap instead, if ordering of objects is required.
	 */
	public static class CountMap<T> extends HashMap<T, Long> {
		private static final long serialVersionUID = -9079906779955923767L;

		public long getCount(T k) {
			return getOrDefault(k, 0L);
		}

		public void increment(T k, long v) {
			long next = getCount(k) + v;
			if (next == 0) {
				remove(k);
			} else {
				put(k, next);
			}
		}

		public static <T> CountMap<T> fromArray(T[] A) {
			CountMap<T> cm = new CountMap<>();
			for (T x : A) {
				cm.increment(x, 1);
			}
			return cm;
		}
	}
}
