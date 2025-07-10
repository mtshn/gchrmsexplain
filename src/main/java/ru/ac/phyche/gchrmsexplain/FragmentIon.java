package ru.ac.phyche.gchrmsexplain;

import java.util.HashMap;

/**
 * Fragment ion. A wrapper class around the molecular formula of an ion.
 */
public class FragmentIon {
	public final String formula;
	public final int[] elements;
	public final int[] elementNumbers;

	/**
	 * c
	 * 
	 * @param data Hash table with molecular formula. Key is atomic number of
	 *             element (carbon 6, nitrogen 7, hydrogen 1, etc.), value is number
	 *             of atoms. Hydrogen is also included
	 */
	FragmentIon(HashMap<Integer, Integer> data) {
		HashMap<Integer, Integer> data1 = new HashMap<Integer, Integer>();
		for (Integer i : data.keySet()) {
			if (data.get(i) > 0) {
				data1.put(i, data.get(i));
			}
		}
		int[] elements_ = new int[data1.keySet().size()];
		int[] elementNumbers_ = new int[data1.keySet().size()];
		int j = 0;
		for (Integer i : data1.keySet()) {
			elements_[j] = i;
			elementNumbers_[j] = data1.get(i);
			j++;
		}
		elements = elements_;
		elementNumbers = elementNumbers_;
		formula = IsotopicPattern.formulaString(elements, elementNumbers);
	}

	/**
	 * 
	 * @param elements       atomic numbers of elements
	 * @param elementNumbers numbers of corresponding atoms
	 */
	FragmentIon(int[] elements, int[] elementNumbers) {
		formula = IsotopicPattern.formulaString(elements, elementNumbers);
		this.elements = elements;
		this.elementNumbers = elementNumbers;
	}

	/**
	 * 
	 * @return Hash table with molecular formula. Key is atomic number of element
	 *         (carbon 6, nitrogen 7, hydrogen 1, etc.), value is number of atoms.
	 *         Hydrogen is also included
	 */
	public HashMap<Integer, Integer> asMap() {
		if (elements.length == 0) {
			throw new RuntimeException("Molecular formula with 0 elements");
		}
		HashMap<Integer, Integer> result = new HashMap<Integer, Integer>();
		if (elements.length != elementNumbers.length) {
			throw new RuntimeException("Invalid ion");
		}
		for (int i = 0; i < elements.length; i++) {
			result.put(elements[i], elementNumbers[i]);
		}
		return result;
	}

	/**
	 * For neutral fragment...
	 * @param threshold see IsotopicPattern class
	 * @return isotopic patten
	 */
	public IsotopicPattern isotopicPattern(float threshold) {
		if (elements.length == 0) {
			throw new RuntimeException("Molecular formula with 0 elements");
		}
		return IsotopicPattern.isotopicPattern(elements, elementNumbers, threshold);
	}

	/**
	 * The result of this method differs from method isotopicPattern by the mass of the electron.
	 * @param threshold see IsotopicPattern class
	 * @return isotopic patten
	 */
	public IsotopicPattern isotopicPatternSingleCation(float threshold) {
		if (elements.length == 0) {
			throw new RuntimeException("Molecular formula with 0 elements");
		}
		return IsotopicPattern.isotopicPatternSingleCation(elements, elementNumbers, threshold);
	}
}
