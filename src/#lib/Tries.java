public class Tries {
	public static class BinaryTrie {
		private Node root = new Node();
		private int bits;

		public BinaryTrie(int bits) {
			this.bits = bits;
		}

		public void insert(int x) {
			Node curr = root;
			for (int b = bits - 1; b >= 0; --b) {
				++curr.countBelow;
				int idx = (x >> b) & 1;
				curr = curr.getOrCreate(idx);
			}
			++curr.countBelow;
		}

		public void remove(int x) {
			Node curr = root;
			for (int b = bits - 1; b >= 0; --b) {
				--curr.countBelow;
				int idx = (x >> b) & 1;
				curr = curr.get(idx);
			}
			--curr.countBelow;
		}

		public int getMaxMatch(int x) {
			Node curr = root;
			int ans = 0;
			for (int b = bits - 1; b >= 0; --b) {
				int idx = (x >> b) & 1;
				if (!curr.has(idx)) {
					idx ^= 1;
				}
				curr = curr.get(idx);
				ans |= idx << b;
			}
			return ans;
		}

		private class Node {
			public Node next0, next1;
			public int countBelow;

			public boolean has(int b) {
				Node u = get(b);
				return u != null && u.countBelow > 0;
			}

			public Node get(int b) {
				if (b == 0) {
					return next0;
				} else {
					return next1;
				}
			}

			public Node getOrCreate(int b) {
				if (b == 0) {
					if (next0 == null) {
						next0 = new Node();
					}
					return next0;
				} else {
					if (next1 == null) {
						next1 = new Node();
					}
					return next1;
				}
			}
		}
	}
}
