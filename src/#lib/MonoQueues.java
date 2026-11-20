public class MonoQueues {
	/**
	 * Queue that keeps track of the minimum value within it.
	 * Useful for finding the min in a sliding window.
	 */
	public static class MinQueueInt {
		private final IntDeque q = new IntDeque();
		private final IntDeque mins = new IntDeque();

		public void offer(int x) {
			while (!mins.isEmpty() && x < mins.peekLast()) {
				mins.removeLast();
			}
			mins.addLast(x);
			q.offer(x);
		}

		public int poll() {
			int x = q.poll();
			if (mins.peekFirst() == x) {
				mins.removeFirst();
			}
			return x;
		}

		public int min() {
			return mins.peekFirst();
		}

		public int size() {
			return q.size();
		}

		private static class IntDeque {
			private int[] arr;
			private int off;
			private int len;

			public IntDeque() {
				this(2);
			}

			public IntDeque(int capacity) {
				this.arr = new int[capacity];
			}

			public void addLast(int x) {
				if (len == arr.length) {
					increaseCapacity();
				}
				int idx = index(off + len);
				arr[idx] = x;
				++len;
			}

			public int peekFirst() {
				return arr[off];
			}

			public int peekLast() {
				int idx = index(off + len - 1);
				return arr[idx];
			}

			public int removeFirst() {
				int ans = peekFirst();
				off = index(off + 1);
				--len;
				return ans;
			}

			public int removeLast() {
				int ans = peekLast();
				--len;
				return ans;
			}

			public void offer(int x) {
				addLast(x);
			}

			public int poll() {
				return removeFirst();
			}

			public int size() {
				return len;
			}

			public boolean isEmpty() {
				return size() == 0;
			}

			private void increaseCapacity() {
				int[] next = new int[arr.length << 1];
				int endLen = arr.length - off;
				System.arraycopy(arr, off, next, 0, endLen);
				System.arraycopy(arr, 0, next, endLen, off);
				arr = next;
				off = 0;
			}

			private int index(int i) {
				if (i >= arr.length) {
					i -= arr.length;
				} else if (i < 0) {
					i += arr.length;
				}
				return i;
			}
		}
	}

	/**
	 * Queue that keeps track of the maximum value within it.
	 * Useful for finding the max in a sliding window.
	 */
	public static class MinQueueLong {
		private final LongDeque q = new LongDeque();
		private final LongDeque mins = new LongDeque();

		public void offer(long x) {
			while (!mins.isEmpty() && x < mins.peekLast()) {
				mins.removeLast();
			}
			mins.addLast(x);
			q.offer(x);
		}

		public long poll() {
			long x = q.poll();
			if (mins.peekFirst() == x) {
				mins.removeFirst();
			}
			return x;
		}

		public long min() {
			return mins.peekFirst();
		}

		public int size() {
			return q.size();
		}

		private static class LongDeque {
			private long[] arr;
			private int off;
			private int len;

			public LongDeque() {
				this(2);
			}

			public LongDeque(int capacity) {
				this.arr = new long[capacity];
			}

			public void addLast(long x) {
				if (len == arr.length) {
					increaseCapacity();
				}
				int idx = index(off + len);
				arr[idx] = x;
				++len;
			}

			public long peekFirst() {
				return arr[off];
			}

			public long peekLast() {
				int idx = index(off + len - 1);
				return arr[idx];
			}

			public long removeFirst() {
				long ans = peekFirst();
				off = index(off + 1);
				--len;
				return ans;
			}

			public long removeLast() {
				long ans = peekLast();
				--len;
				return ans;
			}

			public void offer(long x) {
				addLast(x);
			}

			public long poll() {
				return removeFirst();
			}

			public int size() {
				return len;
			}

			public boolean isEmpty() {
				return size() == 0;
			}

			private void increaseCapacity() {
				long[] next = new long[arr.length << 1];
				int endLen = arr.length - off;
				System.arraycopy(arr, off, next, 0, endLen);
				System.arraycopy(arr, 0, next, endLen, off);
				arr = next;
				off = 0;
			}

			private int index(int i) {
				if (i >= arr.length) {
					i -= arr.length;
				} else if (i < 0) {
					i += arr.length;
				}
				return i;
			}
		}
	}

	/**
	 * Queue that keeps track of the maximum value within it.
	 * Useful for finding the max in a sliding window.
	 */
	public static class MaxQueueInt {
		private final IntDeque q = new IntDeque();
		private final IntDeque maxes = new IntDeque();

		public void offer(int x) {
			while (!maxes.isEmpty() && x > maxes.peekLast()) {
				maxes.removeLast();
			}
			maxes.addLast(x);
			q.offer(x);
		}

		public int poll() {
			int x = q.poll();
			if (maxes.peekFirst() == x) {
				maxes.removeFirst();
			}
			return x;
		}

		public int max() {
			return maxes.peekFirst();
		}

		public int size() {
			return q.size();
		}

		private static class IntDeque {
			private int[] arr;
			private int off;
			private int len;

			public IntDeque() {
				this(2);
			}

			public IntDeque(int capacity) {
				this.arr = new int[capacity];
			}

			public void addLast(int x) {
				if (len == arr.length) {
					increaseCapacity();
				}
				int idx = index(off + len);
				arr[idx] = x;
				++len;
			}

			public int peekFirst() {
				return arr[off];
			}

			public int peekLast() {
				int idx = index(off + len - 1);
				return arr[idx];
			}

			public int removeFirst() {
				int ans = peekFirst();
				off = index(off + 1);
				--len;
				return ans;
			}

			public int removeLast() {
				int ans = peekLast();
				--len;
				return ans;
			}

			public void offer(int x) {
				addLast(x);
			}

			public int poll() {
				return removeFirst();
			}

			public int size() {
				return len;
			}

			public boolean isEmpty() {
				return size() == 0;
			}

			private void increaseCapacity() {
				int[] next = new int[arr.length << 1];
				int endLen = arr.length - off;
				System.arraycopy(arr, off, next, 0, endLen);
				System.arraycopy(arr, 0, next, endLen, off);
				arr = next;
				off = 0;
			}

			private int index(int i) {
				if (i >= arr.length) {
					i -= arr.length;
				} else if (i < 0) {
					i += arr.length;
				}
				return i;
			}
		}
	}

	/**
	 * Queue that keeps track of the maximum value within it.
	 * Useful for finding the max in a sliding window.
	 */
	public static class MaxQueueLong {
		private final LongDeque q = new LongDeque();
		private final LongDeque maxes = new LongDeque();

		public void offer(long x) {
			while (!maxes.isEmpty() && x > maxes.peekLast()) {
				maxes.removeLast();
			}
			maxes.add(x);
			q.offer(x);
		}

		public long poll() {
			long x = q.poll();
			if (maxes.peekFirst() == x) {
				maxes.removeFirst();
			}
			return x;
		}

		public int size() {
			return q.size();
		}

		private static class LongDeque {
			private long[] arr;
			private int off;
			private int len;

			public LongDeque() {
				this(2);
			}

			public LongDeque(int capacity) {
				this.arr = new long[capacity];
			}

			public void addLast(long x) {
				if (len == arr.length) {
					increaseCapacity();
				}
				int idx = index(off + len);
				arr[idx] = x;
				++len;
			}

			public long peekFirst() {
				return arr[off];
			}

			public long peekLast() {
				int idx = index(off + len - 1);
				return arr[idx];
			}

			public long removeFirst() {
				long ans = peekFirst();
				off = index(off + 1);
				--len;
				return ans;
			}

			public long removeLast() {
				long ans = peekLast();
				--len;
				return ans;
			}

			public void add(long x) {
				addLast(x);
			}

			public void offer(long x) {
				addLast(x);
			}

			public long poll() {
				return removeFirst();
			}

			public int size() {
				return len;
			}

			public boolean isEmpty() {
				return size() == 0;
			}

			private void increaseCapacity() {
				long[] next = new long[arr.length << 1];
				int endLen = arr.length - off;
				System.arraycopy(arr, off, next, 0, endLen);
				System.arraycopy(arr, 0, next, endLen, off);
				arr = next;
				off = 0;
			}

			private int index(int i) {
				if (i >= arr.length) {
					i -= arr.length;
				} else if (i < 0) {
					i += arr.length;
				}
				return i;
			}
		}
	}

	/**
	 * Queue that keeps track of the minimum and maximum values within it.
	 * Useful for finding the min and max in a sliding window.
	 */
	public static class MonoQueueInt {
		private final IntDeque q = new IntDeque();
		private final IntDeque mins = new IntDeque();
		private final IntDeque maxes = new IntDeque();

		public void offer(int x) {
			while (!mins.isEmpty() && x < mins.peekLast()) {
				mins.removeLast();
			}
			mins.addLast(x);

			while (!maxes.isEmpty() && x > maxes.peekLast()) {
				maxes.removeLast();
			}
			maxes.addLast(x);

			q.offer(x);
		}

		public int poll() {
			int x = q.poll();
			if (mins.peekFirst() == x) {
				mins.removeFirst();
			}
			if (maxes.peekFirst() == x) {
				maxes.removeFirst();
			}
			return x;
		}

		public int min() {
			return mins.peekFirst();
		}

		public int max() {
			return maxes.peekFirst();
		}

		public int size() {
			return q.size();
		}

		private static class IntDeque {
			private int[] arr;
			private int off;
			private int len;

			public IntDeque() {
				this(2);
			}

			public IntDeque(int capacity) {
				this.arr = new int[capacity];
			}

			public void addLast(int x) {
				if (len == arr.length) {
					increaseCapacity();
				}
				int idx = index(off + len);
				arr[idx] = x;
				++len;
			}

			public int peekFirst() {
				return arr[off];
			}

			public int peekLast() {
				int idx = index(off + len - 1);
				return arr[idx];
			}

			public int removeFirst() {
				int ans = peekFirst();
				off = index(off + 1);
				--len;
				return ans;
			}

			public int removeLast() {
				int ans = peekLast();
				--len;
				return ans;
			}

			public void offer(int x) {
				addLast(x);
			}

			public int poll() {
				return removeFirst();
			}

			public int size() {
				return len;
			}

			public boolean isEmpty() {
				return size() == 0;
			}

			private void increaseCapacity() {
				int[] next = new int[arr.length << 1];
				int endLen = arr.length - off;
				System.arraycopy(arr, off, next, 0, endLen);
				System.arraycopy(arr, 0, next, endLen, off);
				arr = next;
				off = 0;
			}

			private int index(int i) {
				if (i >= arr.length) {
					i -= arr.length;
				} else if (i < 0) {
					i += arr.length;
				}
				return i;
			}
		}
	}

	/**
	 * Queue that keeps track of the minimum and maximum values within it.
	 * Useful for finding the min and max in a sliding window.
	 */
	public static class MonoQueueLong {
		private final LongDeque q = new LongDeque();
		private final LongDeque mins = new LongDeque();
		private final LongDeque maxes = new LongDeque();

		public void offer(int x) {
			while (!mins.isEmpty() && x < mins.peekLast()) {
				mins.removeLast();
			}
			mins.addLast(x);

			while (!maxes.isEmpty() && x > maxes.peekLast()) {
				maxes.removeLast();
			}
			maxes.addLast(x);

			q.offer(x);
		}

		public long poll() {
			long x = q.poll();
			if (mins.peekFirst() == x) {
				mins.removeFirst();
			}
			if (maxes.peekFirst() == x) {
				maxes.removeFirst();
			}
			return x;
		}

		public long min() {
			return mins.peekFirst();
		}

		public long max() {
			return maxes.peekFirst();
		}

		public int size() {
			return q.size();
		}

		private static class LongDeque {
			private long[] arr;
			private int off;
			private int len;

			public LongDeque() {
				this(2);
			}

			public LongDeque(int capacity) {
				this.arr = new long[capacity];
			}

			public void addLast(long x) {
				if (len == arr.length) {
					increaseCapacity();
				}
				int idx = index(off + len);
				arr[idx] = x;
				++len;
			}

			public long peekFirst() {
				return arr[off];
			}

			public long peekLast() {
				int idx = index(off + len - 1);
				return arr[idx];
			}

			public long removeFirst() {
				long ans = peekFirst();
				off = index(off + 1);
				--len;
				return ans;
			}

			public long removeLast() {
				long ans = peekLast();
				--len;
				return ans;
			}

			public void offer(long x) {
				addLast(x);
			}

			public long poll() {
				return removeFirst();
			}

			public int size() {
				return len;
			}

			public boolean isEmpty() {
				return size() == 0;
			}

			private void increaseCapacity() {
				long[] next = new long[arr.length << 1];
				int endLen = arr.length - off;
				System.arraycopy(arr, off, next, 0, endLen);
				System.arraycopy(arr, 0, next, endLen, off);
				arr = next;
				off = 0;
			}

			private int index(int i) {
				if (i >= arr.length) {
					i -= arr.length;
				} else if (i < 0) {
					i += arr.length;
				}
				return i;
			}
		}
	}
}