import java.util.*;

public class TopAndBottom {
	/**
	 * Keeps track of the top two elements inserted.
	 * If `first` and `second` are initialized, `count` will only include elements that were strictly greater than the initial `second` value.
	 * 
	 * NOTE:  It doesn't work with Integer.MIN_VALUE.
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
	 * 
	 * NOTE:  It doesn't work with Long.MIN_VALUE.
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

		@Override
		public String toString() {
			return "[" + first + ", " + second + "]";
		}
	}

	/**
	 * Keeps track of the bottom two elements inserted.
	 * If `first` and `second` are initialized, `count` will only include elements that were strictly less than the initial `second` value.
	 * 
	 * NOTE:  It doesn't work with Integer.MAX_VALUE.
	 */
	public static class BottomTwoInt {
		public int count;
		public int first;
		public int second;

		public BottomTwoInt() {
			this(Integer.MAX_VALUE);
		}

		public BottomTwoInt(int init) {
			this(init, init);
		}

		public BottomTwoInt(int first, int second) {
			this.first = Math.max(first, second);
			this.second = Math.min(first, second);
		}

		public BottomTwoInt(int[] arr) {
			this();
			for (int x : arr) {
				add(x);
			}
		}

		public void add(int x) {
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

		@Override
		public String toString() {
			return "[" + first + ", " + second + "]";
		}
	}

	/**
	 * Keeps track of the bottom two elements inserted.
	 * If `first` and `second` are initialized, `count` will only include elements that were strictly less than the initial `second` value.
	 * 
	 * NOTE:  It doesn't work with Long.MAX_VALUE.
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

		@Override
		public String toString() {
			return "[" + first + ", " + second + "]";
		}
	}

	/**
	 * Keeps track of the top three elements inserted.
	 * 
	 * NOTE:  It doesn't work with Integer.MIN_VALUE.
	 */
	public static class TopThreeInt {
		public int count;
		public int first;
		public int second;
		public int third;

		public TopThreeInt() {
			this.first = this.second = this.third = Integer.MIN_VALUE;
		}

		public void add(int x) {
			if (x < third) {
				return;
			}

			if (x > first) {
				third = second;
				second = first;
				first = x;
			} else if (x > second) {
				third = second;
				second = x;
			} else {
				third = x;
			}
			count = Math.min(2, count + 1);
		}

		@Override
		public String toString() {
			return "[" + first + ", " + second + ", " + third + "]";
		}
	}

	/**
	 * Keeps track of the bottom three elements inserted.
	 * 
	 * NOTE:  It doesn't work with Integer.MAX_VALUE.
	 */
	public static class BottomThreeInt {
		public int count;
		public int first;
		public int second;
		public int third;

		public BottomThreeInt() {
			this.first = this.second = this.third = Integer.MAX_VALUE;
		}

		public void add(int x) {
			if (x > third) {
				return;
			}

			if (x < first) {
				third = second;
				second = first;
				first = x;
			} else if (x < second) {
				third = second;
				second = x;
			} else {
				third = x;
			}
			count = Math.min(2, count + 1);
		}

		@Override
		public String toString() {
			return "[" + first + ", " + second + ", " + third + "]";
		}
	}

	/**
	 * For a multi-set of (N) long values, it keeps (K) largest values on top and remaining (N - K) smallest values on bottom, with ability to balance on (K).
	 * It keeps the total and sum of values in both the top and bottom.
	 */
	public static class SimpleHourglass {
		public final CountMapInt top = new CountMapInt();
		public final CountMapInt bot = new CountMapInt();

		public int totalCount() {
			return top.total + bot.total;
		}

		public void increment(long key) {
			incrementPreferringTop(key);
		}

		public void incrementPreferringTop(long key) {
			if (bot.isEmpty() || key > bot.lastKey()) {
				top.increment(key);
			} else {
				bot.increment(key);
			}
		}

		public void incrementPreferringBot(long key) {
			if (top.isEmpty() || key < top.firstKey()) {
				bot.increment(key);
			} else {
				top.increment(key);
			}
		}

		public void decrement(long key) {
			decrementTopFirst(key);
		}

		public void decrementTopFirst(long key) {
			if (top.containsKey(key)) {
				top.decrement(key);
			} else {
				bot.decrement(key);
			}
		}

		public void decrementBotFirst(long key) {
			if (bot.containsKey(key)) {
				bot.decrement(key);
			} else {
				top.decrement(key);
			}
		}

		public void balanceTop(int targetSize) {
			balanceBot(Math.max(0, top.total + bot.total - targetSize));
		}

		public void balanceBot(int targetSize) {
			while (bot.total < targetSize && !top.isEmpty()) {
				Map.Entry<Long, Integer> lowestKvp = top.firstEntry();
				long key = lowestKvp.getKey();
				int val = lowestKvp.getValue();
				int rem = Math.min(val, targetSize - bot.total);
				top.decrement(key, rem);
				bot.increment(key, rem);
			}
			while (bot.total > targetSize) {
				Map.Entry<Long, Integer> highestKvp = bot.lastEntry();
				long key = highestKvp.getKey();
				int val = highestKvp.getValue();
				int rem = Math.min(val, bot.total - targetSize);
				bot.decrement(key, rem);
				top.increment(key, rem);
			}
		}

		public static class CountMapInt extends TreeMap<Long, Integer> {
			private static final long serialVersionUID = -9164076262452209132L;

			public int total;
			public long sum;

			public int getCount(Long k) {
				return getOrDefault(k, 0);
			}

			public int decrement(Long k) {
				return decrement(k, 1);
			}

			public int decrement(Long k, int v) {
				return increment(k, -v);
			}

			public int increment(Long k) {
				return increment(k, 1);
			}

			public int increment(Long k, int v) {
				total += v;
				sum += k * v;
				int next = getCount(k) + v;
				if (next == 0) {
					remove(k);
				} else {
					put(k, next);
				}
				return next;
			}
		}
	}
}