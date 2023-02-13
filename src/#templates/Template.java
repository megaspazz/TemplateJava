import java.io.*;
import java.util.*;

public class Template {
	public static void solveCase(FastIO io, int testCase) {

	}

	// +---------------+ //
	// | CONFIGURATION | //
	// +---------------+ //

	private static final TestType TEST_TYPE = TestType.SINGLE;
	private static final int START_TEST_CASE = 1;

	private static final boolean USE_THREAD = false;
	private static final long THREAD_STACK_SIZE = 1L << 28;

	private static final String INPUT_FILE = null;
	private static final String OUTPUT_FILE = null;

	// +---------------+ //
	// | TEMPLATE CODE | //
	// +---------------+ //

	private static void solve(FastIO io) {
		switch (TEST_TYPE) {
			case SINGLE: {
				solveCase(io, START_TEST_CASE);
				break;
			}
			case MULTIPLE: {
				final int T = io.nextInt();
				for (int t = 0; t < T; ++t) {
					solveCase(io, START_TEST_CASE + t);
				}
				break;
			}
		}
	}

	private static enum TestType {
		SINGLE,
		MULTIPLE,
	}

	public static class FastIO {
		private InputStream reader;
		private PrintWriter writer;

		private byte[] buf = new byte[1024];
		private int curChar;
		private int numChars;

		public FastIO(InputStream r, OutputStream w) {
			reader = r;
			writer = new PrintWriter(new BufferedWriter(new OutputStreamWriter(w)));
		}

		public int read() {
			if (numChars == -1)
				throw new InputMismatchException();
			if (curChar >= numChars) {
				curChar = 0;
				try {
					numChars = reader.read(buf);
				} catch (IOException e) {
					throw new InputMismatchException();
				}
				if (numChars <= 0)
					return -1;
			}
			return buf[curChar++];
		}

		public String readToEnd() {
			StringBuilder sb = new StringBuilder();
			int c = read();
			while (c >= 0) {
				sb.append((char) c);
				c = read();
			}
			return sb.toString();
		}

		public String nextLine() {
			int c = read();
			while (isSpaceChar(c))
				c = read();
			StringBuilder res = new StringBuilder();
			do {
				res.appendCodePoint(c);
				c = read();
			} while (!isEndOfLine(c));
			return res.toString();
		}

		public String nextString() {
			int c = read();
			while (isSpaceChar(c))
				c = read();
			StringBuilder res = new StringBuilder();
			do {
				res.appendCodePoint(c);
				c = read();
			} while (!isSpaceChar(c));
			return res.toString();
		}

		public long nextLong() {
			int c = read();
			while (isSpaceChar(c))
				c = read();
			int sgn = 1;
			if (c == '-') {
				sgn = -1;
				c = read();
			}
			long res = 0;
			do {
				if (c < '0' || c > '9')
					throw new InputMismatchException();
				res *= 10;
				res += c - '0';
				c = read();
			} while (!isSpaceChar(c));
			return res * sgn;
		}

		public int nextInt() {
			int c = read();
			while (isSpaceChar(c))
				c = read();
			int sgn = 1;
			if (c == '-') {
				sgn = -1;
				c = read();
			}
			int res = 0;
			do {
				if (c < '0' || c > '9')
					throw new InputMismatchException();
				res *= 10;
				res += c - '0';
				c = read();
			} while (!isSpaceChar(c));
			return res * sgn;
		}

		// TODO: read this byte-by-byte like the other read functions.
		public double nextDouble() {
			return Double.parseDouble(nextString());
		}

		public int[] nextIntArray(int n) {
			return nextIntArray(n, 0);
		}

		public int[] nextIntArray(int n, int off) {
			int[] arr = new int[n + off];
			for (int i = 0; i < n; i++) {
				arr[i + off] = nextInt();
			}
			return arr;
		}

		public long[] nextLongArray(int n) {
			return nextLongArray(n, 0);
		}

		public long[] nextLongArray(int n, int off) {
			long[] arr = new long[n + off];
			for (int i = 0; i < n; i++) {
				arr[i + off] = nextLong();
			}
			return arr;
		}

		private boolean isSpaceChar(int c) {
			return c == ' ' || c == '\n' || c == '\r' || c == '\t' || c == -1;
		}

		private boolean isEndOfLine(int c) {
			return c == '\n' || c == '\r' || c == -1;
		}

		public void print(Object... arr) {
			for (int i = 0; i < arr.length; i++) {
				if (i != 0) {
					writer.print(' ');
				}
				writer.print(arr[i]);
			}
		}

        public void print(int... arr) {
			for (int i = 0; i < arr.length; i++) {
				if (i != 0) {
					writer.print(' ');
				}
				writer.print(arr[i]);
			}
        }

        public void print(long... arr) {
			for (int i = 0; i < arr.length; i++) {
				if (i != 0) {
					writer.print(' ');
				}
				writer.print(arr[i]);
			}
        }

		public void println(Object... arr) {
			print(arr);
			writer.println();
		}

        public void println(int... arr) {
            print(arr);
            writer.println();
        }

        public void println(long... arr) {
            print(arr);
            writer.println();
        }

		public void printf(String format, Object... args) {
			print(String.format(format, args));
		}

		public void flush() {
			writer.flush();
		}
	}

	private static void runSolution(FastIO io) {
		solve(io);
		io.flush();
	}

	private static class ThreadedSolution implements Runnable {
		private FastIO io;

		public ThreadedSolution(FastIO io) {
			this.io = io;
		}

		@Override
		public void run() {
			runSolution(io);
		}
	}

	public static void main(String[] args) throws FileNotFoundException, InterruptedException {
		InputStream inStream;
		if (INPUT_FILE == null) {
			inStream = System.in;
		} else {
			inStream = new FileInputStream(INPUT_FILE);
		}

		OutputStream outStream;
		if (OUTPUT_FILE == null) {
			outStream = System.out;
		} else {
			outStream = new FileOutputStream(OUTPUT_FILE);
		}

		FastIO io = new FastIO(inStream, outStream);

		if (USE_THREAD) {
			Thread t = new Thread(null, new ThreadedSolution(io), "ThreadedSolution", THREAD_STACK_SIZE);
			t.start();
			t.join();
		} else {
			runSolution(io);
		}
	}
}