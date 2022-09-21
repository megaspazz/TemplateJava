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
		
		public int last() {
			return arr[pos - 1];
		}
		
		public int removeLast() {
			return arr[--pos];
		}

		public boolean remove(int x) {
			for (int i = 0; i < pos; ++i) {
				if (arr[i] == x) {
					System.arraycopy(arr, i + 1, arr, i, --pos - i);
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

		private void resize(int newCapacity) {
			arr = Arrays.copyOf(arr, newCapacity);
		}
		
		public static IntList of(int... items) {
			IntList lst = new IntList(items.length);
			System.arraycopy(items, 0, lst.arr, 0, items.length);
			lst.pos = items.length;
			return lst;
		}
	}
}
