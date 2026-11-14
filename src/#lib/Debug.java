import java.util.*;
import java.util.function.*;

/**
 * Helpers for printing debug statements that can be easily toggled on and off, and will automatically format arrays.
 * To remove them, an easy way is to delete the debug definitions, and then delete all the lines that have new compilation errors.
 * 
 * In general order from slowest to fastest for `int x` type and `int[] y` types:
 *     1.  debugf("x = %d, expensiveFn(y) = %s%n", () -> new Object[] {x, expensiveFn(y)})
 *     2.  debugf("x = %d, expensiveFn(y) = %s%n", () -> x, () -> expensiveFunction(y))
 *     3.  debugf("x = %d, expensiveFn(y) = %s%n, x, lazy(() -> expensiveFunction(y))
 *     4.  debugf("x = %d, y = %s%n, x, y)
 *     5.  if (DEBUG) debugf(...)
 */
public class Debug {
	public static final boolean DEBUG = true;

	public static void debugf(String format, Object... args) {
		if (DEBUG) {
			Object[] convertedArgs = new Object[args.length];
			for (int i = 0; i < args.length; ++i) {
				convertedArgs[i] = convert(args[i]);
			}
			System.out.format(format, convertedArgs);
		}
	}

	public static void debugf(String format, Supplier<?>... argSuppliers) {
		debugf(format, (Object[]) argSuppliers);
	}

	public static void debugf(String format, Supplier<Object[]> argsSupplier) {
		debugf(format, argsSupplier.get());
	}

	public static Supplier<?> lazy(Supplier<?> s) {
		return s;
	}

	private static Object convert(Object obj) {
		Object val = (obj instanceof Supplier<?> s) ? s.get() : obj;
		return switch (val) {
			case null -> null;
			case int[] arr -> Arrays.toString(arr);
			case long[] arr -> Arrays.toString(arr);
			case char[] arr -> Arrays.toString(arr);
			case double[] arr -> Arrays.toString(arr);
			case boolean[] arr -> Arrays.toString(arr);
			case byte[] arr -> Arrays.toString(arr);
			case short[] arr -> Arrays.toString(arr);
			case float[] arr -> Arrays.toString(arr);
			case Object[] arr -> Arrays.deepToString(arr);
			default -> val;
		};
	}
}
