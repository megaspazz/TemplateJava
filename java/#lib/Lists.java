import java.util.functional.*;

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

		public void add(int x) {
			if (pos >= arr.length) {
				resize(arr.length << 1);
			}
			arr[pos++] = x;
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

		public void forEach(IntConsumer consumer) {
			for (int i = 0; i < pos; ++i) {
				consumer.accept(arr[i]);
			}
		}

		private void resize(int newCapacity) {
			arr = Arrays.copyOf(arr, newCapacity);
		}
	}
}
