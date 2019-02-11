
import java.util.Random;

import Plotter.Plotter;

public class Sorting {

	final static int SELECTION_VS_QUICK_LENGTH = 12;
	final static int MERGE_VS_QUICK_LENGTH = 15;
	final static int MERGE_VS_QUICK_SORTED_LENGTH = 12;
	final static int SELECT_VS_MERGE_LENGTH = 16;
	final static double T = 600.0;
	
	private static void print(double[] arr){
		for(int i = 0; i < arr.length ; i++){
			System.out.print(arr[i] + "  ");
		}
		System.out.println("");
	}
	
	
	/**
	 * Sorts a given array using the quick sort algorithm.
	 * At each stage the pivot is chosen to be the rightmost element of the subarray.
	 * 
	 * Should run in average complexity of O(nlog(n)), and worst case complexity of O(n^2)
	 * 
	 * @param arr - the array to be sorted
	 */
	public static void quickSort(double[] arr){
		quickSort(arr, 0, arr.length - 1);
	}
	
	private static void quickSort(double[] arr, int first, int last){
		if(first < last){
			int q = partition(arr, first, last);
			quickSort(arr, first, q - 1);
			quickSort(arr, q + 1 , last);
		}
	}
	
	private static int partition(double[] arr, int first, int last){
		double pivot = arr[last];
		int wall = first;
		for(int j = first ; j < last ; j++){
			if(arr[j] <= pivot){
				wall = wall + 1;
				double temp = arr[j];
				arr[j] = arr[wall - 1];
				arr[wall - 1] = temp;
			}
		}
		arr[last] = arr[wall];
		arr[wall] = pivot;
		return wall;
	}
	
	/**
	 * Sorts a given array using the merge sort algorithm.
	 * 
	 * Should run in complexity O(nlog(n)) in the worst case.
	 * 
	 * @param arr - the array to be sorted
	 */
	public static void mergeSort(double[] arr){
		mergeSort(arr, 0, arr.length - 1);
	}
	
	private static void mergeSort(double[] arr, int first, int last){
		if(first < last){
			int split = (first + last)/2;
			mergeSort(arr, first, split);
			mergeSort(arr, split + 1, last);
			merge(arr, first, split, last);
			
		}
	}
	
	private static void merge(double[] arr, int first, int split, int last){
		//print(arr);
		int sizeFirst = split - first + 1;
		int sizeSec = last - split;
		double[] left = new double[sizeFirst + 1];
		double[] right = new double[sizeSec + 1];
		for(int i = 0 ; i < left.length - 1 ; i++){
			left[i] = arr[first + i];
		}
		for(int j = 0 ; j < right.length - 1 ; j++){
			right[j] = arr[split + j + 1];
		}
		left[left.length - 1] = Double.POSITIVE_INFINITY;
		right[right.length - 1] = Double.POSITIVE_INFINITY;
		int i = 0;
		int j = 0;
		for(int k = first ; k < last + 1 ; k++){
			if(left[i] <= right[j]){
				arr[k] = left[i];
				i++;
			}else{
				arr[k] = right[j];
				j++;
			}
		}
	}
	

	/**
	 * finds the i'th order statistic of a given array.
	 * 
	 * Should run in complexity O(n) in the worst case.
	 * 
	 * @param arr - the array.
	 * @param i - a number between 0 and arr.length - 1.
	 * @return the number which would be at index i, if the array was to be sorted
	 */
	public static double select (double[] arr, int i){
		return select(arr, 0, arr.length - 1, i);
	}
	
	private static double select(double[] arr, int start, int end, int i){
		if(start == end) return arr[start];
		double[] medianArr = new double[(end -start) / 5 + 1];
		int res = (end - start + 1) % 5;
		int index;
		if(res == 0){
			for(int j = 0; j < medianArr.length; j++){
			index = j*5 + start;
			medianArr[j] = median(arr[index], arr[index + 1], arr[index + 2], arr[index + 3], arr[index + 4]);
			}
		}else{
			double[] lastGroup = new double[res];
			for(int k = 0 ; k < res ; k++){
				lastGroup[k] = arr[end - k];
			}
			quickSort(lastGroup);
			medianArr[medianArr.length - 1] = lastGroup[res / 2];
			for(int s = 0; s < medianArr.length - 1; s++){
				index = s*5 + start;
				medianArr[s] = median(arr[index], arr[index + 1], arr[index + 2], arr[index + 3], arr[index + 4]);
				}
		}
		
		double pivot = select(medianArr, 0, medianArr.length - 1, medianArr.length/2);
		int q = partition(arr, start, end, pivot);
		if(q == i) return arr[q];
		else if(q < i) return select(arr, q+1, end, i);
		else return select(arr, start, q-1, i);
		
	}
	
	private static double median(double a, double b, double c, double d, double e){
		double[] temp = {a, b, c, d, e};
		quickSort(temp);
		return temp[2];
		}
	
	private static int partition(double[] arr, int first, int last, double pivot){
		int index;
		for(index = 0; index < arr.length ; index++){
			if(arr[index] == pivot){
				break;
			}
		}
		double temp = arr[index];
		arr[index] = arr[last];
		arr[last] = temp;
		
		return partition(arr, first, last);
	}
		
	
	/**
	 * Sorts a given array using the selection sort algorithm.
	 * 
	 * Should run in complexity O(n^2) in the worst case.
	 * 
	 * @param arr - the array to be sorted
	 */
	public static void selectionSort(double[] arr){
		for(int i = 0 ; i < arr.length - 1; i++){
			int minIndex = i;
			for(int j = i + 1 ; j < arr.length ; j++){
				if(arr[j] < arr[minIndex]){
					minIndex = j;
				}
			}
		if(minIndex != i){
			double temp = arr[i];
			arr[i] = arr[minIndex];
			arr[minIndex] = temp;	
		}
		}
	}
	
	public static void main(String[] args) {
		//selectionVsQuick();
		//mergeVsQuick();
		//mergeVsQuickOnSortedArray();
		//selectVsMerge();
		double[] arr = {3, 1, 6, 5, 10 ,7 ,2};
		//quickSort(arr);
		//mergeSort(arr);
		//selectionSort(arr);
		System.out.println(select(arr, 4));
		//int q = partition(arr, 0, 7);
		//System.out.println(q);
		//print(arr);
		
	}
	
	/**
	 * Compares the selection sort algorithm against quick sort on random arrays
	 */
	public static void selectionVsQuick(){
		double[] quickTimes = new double[SELECTION_VS_QUICK_LENGTH];
		double[] selectionTimes = new double[SELECTION_VS_QUICK_LENGTH];
		long startTime, endTime;
		Random r = new Random();
		for (int i = 0; i < SELECTION_VS_QUICK_LENGTH; i++) {
			long sumQuick = 0;
			long sumSelection = 0;
			for(int k = 0; k < T; k++){
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = r.nextGaussian() * 5000;
					b[j] = a[j];
				}
				startTime = System.currentTimeMillis();
				quickSort(a);
				endTime = System.currentTimeMillis();
				sumQuick += endTime - startTime;
				startTime = System.currentTimeMillis();
				selectionSort(b);
				endTime = System.currentTimeMillis();
				sumSelection += endTime - startTime;
			}
			quickTimes[i] = sumQuick/T;
			selectionTimes[i] = sumSelection/T;
		}
		Plotter.plot("quick sort", quickTimes, "selection sort", selectionTimes);
	}
	
	/**
	 * Compares the merge sort algorithm against quick sort on random arrays
	 */
	public static void mergeVsQuick(){
		double[] quickTimes = new double[MERGE_VS_QUICK_LENGTH];
		double[] mergeTimes = new double[MERGE_VS_QUICK_LENGTH];
		long startTime, endTime;
		Random r = new Random();
		for (int i = 0; i < MERGE_VS_QUICK_LENGTH; i++) {
			long sumQuick = 0;
			long sumMerge = 0;
			for (int k = 0; k < T; k++) {
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = r.nextGaussian() * 5000;
					b[j] = a[j];
				}
				startTime = System.currentTimeMillis();
				quickSort(a);
				endTime = System.currentTimeMillis();
				sumQuick += endTime - startTime;
				startTime = System.currentTimeMillis();
				mergeSort(b);
				endTime = System.currentTimeMillis();
				sumMerge += endTime - startTime;
			}
			quickTimes[i] = sumQuick/T;
			mergeTimes[i] = sumMerge/T;
		}
		Plotter.plot("quick sort", quickTimes, "merge sort", mergeTimes);
	}

	/**
	 * Compares the merge sort algorithm against quick sort on pre-sorted arrays
	 */
	public static void mergeVsQuickOnSortedArray(){
		double[] quickTimes = new double[MERGE_VS_QUICK_SORTED_LENGTH];
		double[] mergeTimes = new double[MERGE_VS_QUICK_SORTED_LENGTH];
		long startTime, endTime;
		Random r = new Random();
		for (int i = 0; i < MERGE_VS_QUICK_SORTED_LENGTH; i++) {
			long sumQuick = 0;
			long sumMerge = 0;
			for (int k = 0; k < T; k++) {
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = j;
					b[j] = j;
				}
				startTime = System.currentTimeMillis();
				quickSort(a);
				endTime = System.currentTimeMillis();
				sumQuick += endTime - startTime;
				startTime = System.currentTimeMillis();
				mergeSort(b);
				endTime = System.currentTimeMillis();
				sumMerge += endTime - startTime;
			}
			quickTimes[i] = sumQuick/T;
			mergeTimes[i] = sumMerge/T;
		}
		Plotter.plot("quick sort on sorted array", quickTimes, "merge sort on sorted array", mergeTimes);
	}

	/**
	 * Compares the select algorithm against sorting an array.
	 */
	public static void selectVsMerge(){
		double[] mergeTimes = new double[MERGE_VS_QUICK_LENGTH];
		double[] selectTimes = new double[MERGE_VS_QUICK_LENGTH];
		long startTime, endTime;
		double x;
		Random r = new Random();
		for (int i = 0; i < MERGE_VS_QUICK_LENGTH; i++) {
			long sumMerge = 0;
			long sumSelect = 0;
			for (int k = 0; k < T; k++) {
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = r.nextGaussian() * 5000;
					b[j] = a[j];
				}
				int index = (int)(Math.random() * size);
				startTime = System.currentTimeMillis();
				mergeSort(a);
				x = a[index];
				endTime = System.currentTimeMillis();
				sumMerge += endTime - startTime;
				startTime = System.currentTimeMillis();
				x = select(b, index);
				endTime = System.currentTimeMillis();
				sumSelect += endTime - startTime;
			}
			mergeTimes[i] = sumMerge/T;
			selectTimes[i] = sumSelect/T;
		}
		Plotter.plot("merge sort and select", mergeTimes, "select", selectTimes);
	}
}
