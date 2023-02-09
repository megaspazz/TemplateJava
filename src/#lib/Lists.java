import java.util.*;
import java.util.function.*;

public class Lists {
	/**
	 * Circular buffer of int values, can be used as:
	 *   - ArrayList: values are added to end.
	 *   - Queue: values are added to end and removed from front.
	 *   - Stack: values are added to and removed from front.
	 */
	public static class IntDeque {
		private int[] arr;
		private int off;
		private int len;

		public IntDeque() {
			this(2);
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

		public void add(int x) {
			addLast(x);
		}

		public void offer(int x) {
			addLast(x);
		}

		public int poll() {
			return removeFirst();
		}

		public void push(int x) {
			addFirst(x);
		}

		public int pop() {
			return removeFirst();
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

		public int[] toArray() {
			if (len == 0) {
				return new int[0];
			}
			int idx = index(off + len);
			if (idx > off) {
				return Arrays.copyOfRange(arr, off, idx);
			}
			int[] A = new int[len];
			int endLen = arr.length - off;
			int startLen = len - endLen;
			System.arraycopy(arr, off, A, 0, endLen);
			System.arraycopy(arr, 0, A, endLen, startLen);
			return A;
		}

		@Override
		public String toString() {
			StringBuilder sb = new StringBuilder();
			sb.append('[');
			printToBuffer(sb, ", ");
			sb.append(']');
			return sb.toString();
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

		private void printToBuffer(StringBuilder sb, CharSequence sep) {
			for (int i = 0; i < len; ++i) {
				if (i > 0) {
					sb.append(sep);
				}
				sb.append(get(i));
			}
		}

		public static IntDeque of(int... arr) {
			IntDeque deq = new IntDeque();
			for (int x : arr) {
				deq.add(x);
			}
			return deq;
		}
	}

	/**
	 * Circular buffer of long values, can be used as:
	 *   - ArrayList: values are added to end.
	 *   - Queue: values are added to end and removed from front.
	 *   - Stack: values are added to and removed from front.
	 */
	public static class LongDeque {
		private long[] arr;
		private int off;
		private int len;

		public LongDeque() {
			this(2);
		}

		public LongDeque(int capacity) {
			this.arr = new long[capacity];
		}

		public void addFirst(long x) {
			if (len == arr.length) {
				increaseCapacity();
			}
			if (off == 0) {
				off = arr.length;
			}
			arr[--off] = x;
			++len;
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

		public void push(long x) {
			addFirst(x);
		}

		public long pop() {
			return removeFirst();
		}

		public long peek() {
			return peekFirst();
		}

		public long get(int i) {
			if (i >= len) {
				throw new ArrayIndexOutOfBoundsException(String.format("index %d out of range [0, %d)", i, len));
			}
			int idx = index(i + off);
			return arr[idx];
		}

		public void set(int i, long x) {
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

		public long[] toArray() {
			if (len == 0) {
				return new long[0];
			}
			int idx = index(off + len);
			if (idx > off) {
				return Arrays.copyOfRange(arr, off, idx);
			}
			long[] A = new long[len];
			int endLen = arr.length - off;
			int startLen = len - endLen;
			System.arraycopy(arr, off, A, 0, endLen);
			System.arraycopy(arr, 0, A, endLen, startLen);
			return A;
		}

		@Override
		public String toString() {
			StringBuilder sb = new StringBuilder();
			sb.append('[');
			printToBuffer(sb, ", ");
			sb.append(']');
			return sb.toString();
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

		private void printToBuffer(StringBuilder sb, CharSequence sep) {
			for (int i = 0; i < len; ++i) {
				if (i > 0) {
					sb.append(sep);
				}
				sb.append(get(i));
			}
		}

		public static LongDeque of(long... arr) {
			LongDeque deq = new LongDeque();
			for (long x : arr) {
				deq.add(x);
			}
			return deq;
		}
	}

	public static class IntList {
		public int[] arr;
		public int pos;

		public IntList(int capacity) {
			this.arr = new int[capacity];
		}

		public IntList() {
			this(2);
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
			this(2);
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

		public static LongList of(long... items) {
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
}
