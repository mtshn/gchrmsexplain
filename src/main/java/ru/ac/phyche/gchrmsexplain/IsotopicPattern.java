package ru.ac.phyche.gchrmsexplain;

import java.util.HashMap;

// https://github.com/emptyport/isotope-abundances/blob/master/ISOTOPES.json
/**
 * A method for generating isotopic distributions. An atom-by-atom sequential
 * merging algorithm is used. If there are isotopic patterns (sets of
 * mz-intensity pairs) for two fragments A and B, then it is easy to calculate
 * the pattern for the combined fragment AB. In this way, element by element,
 * atom by atom can be added. This method is acceptable for gas chromatography.
 * Information on natural isotope distributions for elements is taken from here:
 * https://github.com/emptyport/isotope-abundances/blob/master/ISOTOPES.json
 * https://github.com/emptyport/isotope-abundances It's MIT licensed project,
 * original data from NIST:
 * https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-column-descriptions
 */
public class IsotopicPattern {

	private HashMap<Float, Float> fractionByMZ = new HashMap<Float, Float>();

	/**
	 * all mz in right order from little to large
	 * 
	 * @return
	 */
	public float[] sortedMZ() {
		float[] result = new float[fractionByMZ.size()];
		int j = 0;
		for (float mz : fractionByMZ.keySet()) {
			result[j] = mz;
			j++;
		}
		return result;
	}

	/**
	 * Add new mz-intensity pair
	 * 
	 * @param mz     mz
	 * @param intens intensity
	 */
	public void addIntens(float mz, float intens) {
		if (fractionByMZ.get(mz) == null) {
			fractionByMZ.put(mz, intens);
		} else {
			fractionByMZ.put(mz, fractionByMZ.get(mz) + intens);
		}
	}

	/**
	 * Should be 1, however can differ due to pruning and uncertains
	 * 
	 * @return sum of all probabilities of different isotopomers
	 */
	public float sumProbabilities() {
		float result = 0;
		for (float mz1 : this.fractionByMZ.keySet()) {
			result = result + this.fractionByMZ.get(mz1);
		}
		return result;
	}

	/**
	 * 
	 * @return the isotopic pattern as map. Key - mz (EXACT float values! No
	 *         arithmetics is allowed!), value - intensity. Sum of intensities - 1
	 */
	public HashMap<Float, Float> asMap() {
		return fractionByMZ;
	}

	public IsotopicPattern mergeWith(IsotopicPattern a, float threshold) {
		IsotopicPattern result = new IsotopicPattern();
		for (float mz1 : this.fractionByMZ.keySet()) {
			for (float mz2 : a.fractionByMZ.keySet()) {
				float fraction1 = this.fractionByMZ.get(mz1);
				float fraction2 = a.fractionByMZ.get(mz2);
				if (fraction1 * fraction2 > threshold) {
					result.addIntens(mz1 + mz2, fraction1 * fraction2);
				}
			}
		}
		return result;
	}

	/**
	 * The periodic table
	 * 
	 * @param n element number
	 * @return symbol (6 -> C, 7 -> N etc.)
	 */
	public static String symbolByAtomicNumber(int n) {
		if (n == 1) {
			return ("H");
		}
		if (n == 6) {
			return ("C");
		}
		if (n == 7) {
			return ("N");
		}
		if (n == 8) {
			return ("O");
		}
		if (n == 9) {
			return ("F");
		}
		if (n == 14) {
			return ("Si");
		}
		if (n == 15) {
			return ("P");
		}
		if (n == 16) {
			return ("S");
		}
		if (n == 17) {
			return ("Cl");
		}
		if (n == 33) {
			return ("As");
		}
		if (n == 34) {
			return ("Se");
		}
		if (n == 35) {
			return ("Br");
		}
		if (n == 52) {
			return ("Te");
		}
		if (n == 53) {
			return ("I");
		}
		throw new RuntimeException("Unsupported element! " + n);
	}

	/**
	 * The periodic table
	 * 
	 * @param s element symbol (C - carbon, Br - bromine etc)
	 * @return
	 */
	public static int atomicNumberBySymbol(String s) {
		if (s.equals("H")) {
			return 1;
		}
		if (s.equals("C")) {
			return 6;
		}
		if (s.equals("N")) {
			return 7;
		}
		if (s.equals("O")) {
			return 8;
		}
		if (s.equals("F")) {
			return 9;
		}
		if (s.equals("Si")) {
			return 14;
		}
		if (s.equals("P")) {
			return 15;
		}
		if (s.equals("S")) {
			return 16;
		}
		if (s.equals("Cl")) {
			return 17;
		}
		if (s.equals("As")) {
			return 33;
		}
		if (s.equals("Se")) {
			return 34;
		}
		if (s.equals("Br")) {
			return 35;
		}
		if (s.equals("Te")) {
			return 52;
		}
		if (s.equals("I")) {
			return 53;
		}
		throw new RuntimeException("Unsupported element! " + s);
	}

	/**
	 * Isotopic pattern for single atom of given element. Data are given from this
	 * source: https://github.com/emptyport/isotope-abundances
	 * 
	 * @param el - atomic number of element
	 * @return isotopic pattern for single neutral atom
	 */
	public static IsotopicPattern singleAtomsOfAnElement(int el) {
		IsotopicPattern result = new IsotopicPattern();
		if (el == 1) {// H
			result.fractionByMZ.put(1.00782503223f, 0.999885f);
			result.fractionByMZ.put(2.01410177812f, 0.000115f);
		}
		if (el == 5) {// B
			result.fractionByMZ.put(10.01293695f, 0.199f);
			result.fractionByMZ.put(11.00930536f, 0.801f);
		}

		if (el == 6) {// C
			result.fractionByMZ.put(12f, 0.9893f);
			result.fractionByMZ.put(13.00335483507f, 0.0107f);
		}

		if (el == 7) {// N
			result.fractionByMZ.put(14.00307400443f, 0.99636f);
			result.fractionByMZ.put(15.00010889888f, 0.00364f);
		}

		if (el == 8) {// O
			result.fractionByMZ.put(15.99491461957f, 0.99757f);
			result.fractionByMZ.put(16.9991317565f, 0.00038f);
			result.fractionByMZ.put(17.99915961286f, 0.00205f);
		}

		if (el == 9) {// F
			result.fractionByMZ.put(18.99840316273f, 1f);
		}

		if (el == 14) {// Si
			result.fractionByMZ.put(27.97692653465f, 0.92223f);
			result.fractionByMZ.put(28.9764946649f, 0.04685f);
			result.fractionByMZ.put(29.973770136f, 0.03092f);
		}
		if (el == 15) {// P
			result.fractionByMZ.put(30.97376199842f, 1f);
		}
		if (el == 16) {// S
			result.fractionByMZ.put(31.9720711744f, 0.9499f);
			result.fractionByMZ.put(32.9714589098f, 0.0075f);
			result.fractionByMZ.put(33.967867004f, 0.0425f);
			result.fractionByMZ.put(35.96708071f, 0.0001f);
		}

		if (el == 17) {// Cl
			result.fractionByMZ.put(34.968852682f, 0.7576f);
			result.fractionByMZ.put(36.965902602f, 0.2424f);
		}

		if (el == 33) {// As
			result.fractionByMZ.put(74.92159457f, 1f);
		}

		if (el == 34) {// Se
			result.fractionByMZ.put(73.922475934f, 0.0089f);
			result.fractionByMZ.put(75.919213704f, 0.0937f);
			result.fractionByMZ.put(76.919914154f, 0.0763f);
			result.fractionByMZ.put(77.91730928f, 0.2377f);
			result.fractionByMZ.put(79.9165218f, 0.4961f);
			result.fractionByMZ.put(81.9166995f, 0.0873f);
		}

		if (el == 35) {// Br
			result.fractionByMZ.put(78.9183376f, 0.5069f);
			result.fractionByMZ.put(80.9162897f, 0.4931f);
		}
		if (el == 52) {// Te
			result.fractionByMZ.put(119.9040593f, 0.0009f);
			result.fractionByMZ.put(121.9030435f, 0.0255f);
			result.fractionByMZ.put(122.9042698f, 0.0089f);
			result.fractionByMZ.put(123.9028171f, 0.0474f);
			result.fractionByMZ.put(124.9044299f, 0.0707f);
			result.fractionByMZ.put(125.9033109f, 0.1884f);
			result.fractionByMZ.put(127.90446128f, 0.3174f);
			result.fractionByMZ.put(129.906222748f, 0.3408f);
		}

		if (el == 53) {// I
			result.fractionByMZ.put(126.9044719f, 1f);
		}
		return result;
	}

	/**
	 * Isotopic pattern for molecule X_n, where X - is an element For neutral
	 * molecule
	 * 
	 * @param el        element (atomic number)
	 * @param n         number of atoms
	 * @param threshold Ignore isotopomers whose fraction in the total sum of
	 *                  isotopomers is lower than threshold
	 * @return
	 */
	public static IsotopicPattern multipleAtomsOfAnElement(int el, int n, float threshold) {
		IsotopicPattern result = singleAtomsOfAnElement(el);
		for (int i = 0; i < n - 1; i++) {
			result = result.mergeWith(singleAtomsOfAnElement(el), threshold);
		}
		return result;
	}

	/**
	 * Isotopic pattern for neutral molecule X_n
	 * 
	 * @param elements  atomic numbers of elements
	 * @param numbers   numbers of atoms of each element
	 * @param threshold Ignore isotopomers whose fraction in the total sum of
	 *                  isotopomers is lower than threshold
	 * @return pattern The isotopic distribution is constructed using step-by-step
	 *         addition. In this case, the threshold is applied at each step. That
	 *         is, first we take, for example, the element C, it has two isotopes.
	 *         We add another C to it, we get C2, all isotopomers whose share turned
	 *         out to be below the threshold, we discard, add the next atom. And so
	 *         on for each element: C4, H3, O3, for example. Then we combine the
	 *         element with the element.
	 */
	public static IsotopicPattern isotopicPattern(int[] elements, int[] numbers, float threshold) {
		IsotopicPattern result = multipleAtomsOfAnElement(elements[0], numbers[0], threshold);
		for (int i = 1; i < elements.length; i++) {
			result = result.mergeWith(multipleAtomsOfAnElement(elements[i], numbers[i], threshold), threshold);
		}
		return result;
	}
	/**
	 * 
	 * @return the m/z with the most abundant intensity
	 */
	public float mainMZ() {
		float maxIntens = 0;
		float maxMZ = -1;
		for (float mz : this.fractionByMZ.keySet()) {
			if (fractionByMZ.get(mz) > maxIntens) {
				maxIntens = fractionByMZ.get(mz);
				maxMZ = mz;
			}
		}
		return maxMZ;
	}

	/** 
	 * Absolutely the same as isotopicPattern, but for cation. The difference is mass of electron: 0.0005 Da.
	 * @param elements atomic numbers
	 * @param numbers numbers of atoms
	 * @param threshold (fragment ions with fraction of total current below threshold are excluded)
	 * @return pattern
	 */
	public static IsotopicPattern isotopicPatternSingleCation(int[] elements, int[] numbers, float threshold) {
		IsotopicPattern result = isotopicPattern(elements, numbers, threshold);
		HashMap<Float, Float> r1 = result.fractionByMZ;
		HashMap<Float, Float> r2 = new HashMap<Float, Float>();
		for (Float mz : r1.keySet()) {
			r2.put(mz - 0.0005f, r1.get(mz));
		}
		result.fractionByMZ = r2;
		return result;
	}

	/**
	 * Formula like C4H10
	 * @param elements atomic numbers
	 * @param numbers numbers of atoms
	 * @return formula
	 */
	public static String formulaString(int[] elements, int[] numbers) {
		String result = "";
		for (int j = 0; j < 100; j++) {
			for (int i = 0; i < elements.length; i++) {
				if (elements[i] == j) {
					result += symbolByAtomicNumber(elements[i]) + numbers[i];
				}
			}
		}
		return result.trim();
	}

}
