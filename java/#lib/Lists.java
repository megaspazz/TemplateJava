import java.util.function.*;

public class Lists {
	public static class IntList {
		public int[] arr;
		public int pos;

		public IntList(int capacity) {
			this.arr = new int[capacity];
		}

		public IntList() {
			this(10);
		}

		public IntList(IntList other) {
			this.arr = Arrays.copyOfRange(other.arr, 0, other.pos);
			this.pos = other.pos;
		}

		public void add(int x) {
			if (pos >= arr.length) {
				resize(arr.length << 1);
			}
			arr[pos++] = x;
		}

		public void insert(int i, int x) {
			if (pos >= arr.length) {
				resize(arr.length << 1);
			}
			System.arraycopy(arr, i, arr, i + 1, pos++ - i);
			arr[i] = x;
		}

		public int get(int i) {
			return arr[i];
		}

		public void set(int i, int x) {
			arr[i] = x;
		}

		public void clear() {
			pos = 0;
		}

		public int size() {
			return pos;
		}

		public boolean isEmpty() {
			return pos == 0;
		}

		public int last() {
			return arr[pos - 1];
		}

		public int removeLast() {
			return arr[--pos];
		}

		public void push(int x) {
			add(x);
		}

		public int pop() {
			return removeLast();
		}

		public boolean remove(int x) {
			for (int i = 0; i < pos; ++i) {
				if (arr[i] == x) {
					removeAt(i);
					return true;
				}
			}
			return false;
		}

		public void removeAt(int i) {
			System.arraycopy(arr, i + 1, arr, i, --pos - i);
		}

		public IntList subList(int fromIndex, int toIndexExclusive) {
			IntList lst = new IntList();
			lst.arr = Arrays.copyOfRange(arr, fromIndex, toIndexExclusive);
			lst.pos = toIndexExclusive - fromIndex;
			return lst;
		}

		public void forEach(IntConsumer consumer) {
			for (int i = 0; i < pos; ++i) {
				consumer.accept(arr[i]);
			}
		}

		public int[] toArray() {
			return Arrays.copyOf(arr, pos);
		}

		public static IntList of(int... items) {
			IntList lst = new IntList(items.length);
			System.arraycopy(items, 0, lst.arr, 0, items.length);
			lst.pos = items.length;
			return lst;
		}

		public String join(CharSequence sep) {
			StringBuilder sb = new StringBuilder();
			joinToBuffer(sb, sep);
			return sb.toString();
		}

		@Override
		public String toString() {
			StringBuilder sb = new StringBuilder();
			sb.append('[');
			joinToBuffer(sb, ", ");
			sb.append(']');
			return sb.toString();
		}

		private void resize(int newCapacity) {
			arr = Arrays.copyOf(arr, newCapacity);
		}

		private void joinToBuffer(StringBuilder sb, CharSequence sep) {
			for (int i = 0; i < pos; ++i) {
				if (i > 0) {
					sb.append(sep);
				}
				sb.append(arr[i]);
			}
		}
	}

	public static class LongList {
		public long[] arr;
		public int pos;

		public LongList(int capacity) {
			this.arr = new long[capacity];
		}

		public LongList() {
			this(10);
		}

		public LongList(LongList other) {
			this.arr = Arrays.copyOfRange(other.arr, 0, other.pos);
			this.pos = other.pos;
		}

		public void add(long x) {
			if (pos >= arr.length) {
				resize(arr.length << 1);
			}
			arr[pos++] = x;
		}

		public void insert(int i, long x) {
			if (pos >= arr.length) {
				resize(arr.length << 1);
			}
			System.arraycopy(arr, i, arr, i + 1, pos++ - i);
			arr[i] = x;
		}

		public long get(int i) {
			return arr[i];
		}

		public void set(int i, long x) {
			arr[i] = x;
		}

		public void clear() {
			pos = 0;
		}

		public int size() {
			return pos;
		}

		public boolean isEmpty() {
			return pos == 0;
		}

		public long last() {
			return arr[pos - 1];
		}

		public long removeLast() {
			return arr[--pos];
		}

		public void push(int x) {
			add(x);
		}

		public long pop() {
			return removeLast();
		}

		public boolean remove(long x) {
			for (int i = 0; i < pos; ++i) {
				if (arr[i] == x) {
					removeAt(i);
					return true;
				}
			}
			return false;
		}

		public void removeAt(int i) {
			System.arraycopy(arr, i + 1, arr, i, --pos - i);
		}

		public LongList subList(int fromIndex, int toIndexExclusive) {
			LongList lst = new LongList();
			lst.arr = Arrays.copyOfRange(arr, fromIndex, toIndexExclusive);
			lst.pos = toIndexExclusive - fromIndex;
			return lst;
		}

		public void forEach(LongConsumer consumer) {
			for (int i = 0; i < pos; ++i) {
				consumer.accept(arr[i]);
			}
		}

		public long[] toArray() {
			return Arrays.copyOf(arr, pos);
		}

		public static LongList of(int... items) {
			LongList lst = new LongList(items.length);
			System.arraycopy(items, 0, lst.arr, 0, items.length);
			lst.pos = items.length;
			return lst;
		}

		public String join(CharSequence sep) {
			StringBuilder sb = new StringBuilder();
			joinToBuffer(sb, sep);
			return sb.toString();
		}

		@Override
		public String toString() {
			StringBuilder sb = new StringBuilder();
			sb.append('[');
			joinToBuffer(sb, ", ");
			sb.append(']');
			return sb.toString();
		}

		private void resize(int newCapacity) {
			arr = Arrays.copyOf(arr, newCapacity);
		}

		private void joinToBuffer(StringBuilder sb, CharSequence sep) {
			for (int i = 0; i < pos; ++i) {
				if (i > 0) {
					sb.append(sep);
				}
				sb.append(arr[i]);
			}
		}
	}

	public static class IntDeque {
		private int[] arr;
		private int off;
		private int len;
		
		public IntDeque() {
			this(12);
		}
		
		public IntDeque(int capacity) {
			this.arr = new int[capacity];
		}
		
		public void addFirst(int x) {
			if (len == arr.length) {
				increaseCapacity();
			}
			if (off == 0) {
				off = arr.length;
			}
			arr[--off] = x;
			++len;
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
			int idx = index(off + len);
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
		
		public void add(int x) {
			addLast(x);
		}
		
		public void offer(int x) {
			addLast(x);
		}
		
		public int poll() {
			return removeFirst();
		}
		
		public int peekQueue() {
			return peekFirst();
		}
		
		public void push(int x) {
			addLast(x);
		}
		
		public int pop() {
			return removeLast();
		}
		
		public int peekStack() {
			return peekLast();
		}
		
		public int peek() {
			return peekFirst();
		}
		
		public int get(int i) {
			if (i >= len) {
				throw new ArrayIndexOutOfBoundsException(String.format("index %d out of range [0, %d)", i, len));
			}
			int idx = index(i + off);
			return arr[idx];
		}
		
		public void set(int i, int x) {
			if (i >= len) {
				throw new ArrayIndexOutOfBoundsException(String.format("index %d out of range [0, %d)", i, len));
			}
			int idx = index(i + off);
			arr[idx] = x;
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
			}
			return i;
		}
	}

	/**
	 * This is an priority queue containing int values in descending order.
	 *   - All operations run in amortized O(1) time.
	 *   - However, it can only offer values that are greater-than-or-equal-to the last value of poll().
	 *   - Before the first poll() call, any values can be inserted.
	 */
	public static class AscendingDiscretePriorityQueue<T> {
		private int idx;
		private ArrayList<ArrayDeque<T>> queues;
		private IntConverter<T> converter;

		public AscendingDiscretePriorityQueue(IntConverter<T> converter) {
			this.queues = new ArrayList<>();
			this.converter = converter;
		}

		public boolean offer(T item) {
			int x = converter.toInt(item);
			if (x < idx) {
				return false;
			}
			getOrCreate(x).offer(item);
			return true;
		}

		public T peek() {
			return nextQueue().peek();
		}

		public T poll() {
			return nextQueue().poll();
		}

		public boolean isEmpty() {
			return nextQueue() == null;
		}

		private ArrayDeque<T> nextQueue() {
			while (idx < queues.size()) {
				ArrayDeque<T> q = queues.get(idx);
				if (q != null && !q.isEmpty()) {
					return q;
				}
				++idx;
			}
			return null;
		}

		private ArrayDeque<T> getOrCreate(int i) {
			while (queues.size() <= i) {
				queues.add(null);
			}
			ArrayDeque<T> q = queues.get(i);
			if (q == null) {
				q = new ArrayDeque<>();
				queues.set(i, q);
			}
			return q;
		}

		public static interface IntConverter<T> {
			public int toInt(T item);
		}
	}
}
