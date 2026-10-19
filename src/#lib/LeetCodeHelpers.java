import java.util.*;

public class LeetCodeHelpers {
	public static int[] toIntArray(List<Integer> lst) {
		final int N = lst.size();

		int[] A = new int[N];
		int idx = 0;
		for (int x : lst) {
			A[idx++] = x;
		}
		return A;
	}

	public static int[][] toIntMatrix(List<List<Integer>> grid) {
		final int N = grid.size();

		int[][] A = new int[N][];
		int idx = 0;
		for (List<Integer> row : grid) {
			A[idx++] = toIntArray(row);
		}
		return A;
	}
}
