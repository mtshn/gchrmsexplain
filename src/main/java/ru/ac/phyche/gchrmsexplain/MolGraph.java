package ru.ac.phyche.gchrmsexplain;

import java.util.ArrayList;
import java.util.HashMap;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

/**
 * This class is a molecular graph. Atoms are heavy only (all except hydrogen).
 * Only common non-metals. The number of hydrogens is given for each atom, but
 * they are not part of the graph.
 */
public class MolGraph {
	private int[][] edges;
	private int[] elementAtomicNumber;
	private int[] nHydrogens;

	private int[] possibleHChangeFragment;
	private int maxFluorineMigration = 2;

	/**
	 * Molecular formula.
	 * 
	 * @return Hash table with molecular formula. Key is atomic number of element
	 *         (carbon 6, nitrogen 7, hydrogen 1, etc.), value is number of atoms.
	 *         Hydrogen is also included
	 */
	public HashMap<Integer, Integer> molecularFormula() {
		HashMap<Integer, Integer> molecularFormula = new HashMap<Integer, Integer>();
		int nHydrogens = 0;
		for (int i = 0; i < this.elementAtomicNumber.length; i++) {
			nHydrogens = nHydrogens + this.nHydrogens[i];
			int atomicNumber = this.elementAtomicNumber[i];
			if (molecularFormula.get(atomicNumber) == null) {
				molecularFormula.put(atomicNumber, 1);
			} else {
				molecularFormula.put(atomicNumber, 1 + molecularFormula.get(atomicNumber));
			}
		}
		if (nHydrogens != 0) {
			molecularFormula.put(1, nHydrogens);
		}
		return molecularFormula;
	}

	/**
	 * Molecular formula of fragment
	 * 
	 * @param atomsIncluded array have length equal to number of atoms. True if an
	 *                      atom is in fragment. Empty fragment - all false, full
	 *                      mol - all trues. Atom are included (or excluded) with
	 *                      all corresponding implicit hydrogens.
	 * @return Hash table with molecular formula of fragment. Key is atomic number
	 *         of element (carbon 6, nitrogen 7, hydrogen 1, etc.), value is number
	 *         of atoms. Hydrogen is also included
	 */
	public HashMap<Integer, Integer> molecularFormulaFragment(boolean atomsIncluded[]) {
		HashMap<Integer, Integer> molecularFormula = new HashMap<Integer, Integer>();
		int nHydrogens = 0;
		if (atomsIncluded.length != this.nHydrogens.length) {
			throw (new RuntimeException("Unexpected error while fragmenting graph..."));
		}
		for (int i = 0; i < this.elementAtomicNumber.length; i++) {
			if (atomsIncluded[i]) {
				nHydrogens = nHydrogens + this.nHydrogens[i];
				int atomicNumber = this.elementAtomicNumber[i];
				if (molecularFormula.get(atomicNumber) == null) {
					molecularFormula.put(atomicNumber, 1);
				} else {
					molecularFormula.put(atomicNumber, 1 + molecularFormula.get(atomicNumber));
				}
			}
		}
		if (nHydrogens != 0) {
			molecularFormula.put(1, nHydrogens);
		}
		return molecularFormula;
	}

	/**
	 * Hydrogens
	 * 
	 * @return sum of numbers of hydrogens over all atoms. Total number of hydrogens
	 */
	public int nHydrogensAll() {
		int nH = 0;
		for (int i = 0; i < nHydrogens.length; i++) {
			nH += nHydrogens[i];
		}
		if (nHydrogens.length != nAtoms()) {
			throw new RuntimeException("Invalid molecular graph");
		}
		return nH;
	}

	/**
	 * n CH3
	 * 
	 * @param atomsIncluded array have length equal to number of atoms. True if an
	 *                      atom is in fragment. Empty fragment - all false, full
	 *                      mol - all trues. Atom are included (or excluded) with
	 *                      all corresponding implicit hydrogens.
	 * @return number of methyl groups (C with 3 hydrogens) in the fragment
	 */
	public int nMethyls(boolean atomsIncluded[]) {
		int nCH3 = 0;
		for (int i = 0; i < elementAtomicNumber.length; i++) {
			if (nHydrogens[i] == 3) {
				if (atomsIncluded[i]) {
					nCH3++;
				}
			}
		}
		return nCH3;
	}

	/**
	 * n NO2
	 * 
	 * @param atomsIncluded array have length equal to number of atoms. True if an
	 *                      atom is in fragment. Empty fragment - all false, full
	 *                      mol - all trues. Atom are included (or excluded) with
	 *                      all corresponding implicit hydrogens.
	 * @return number of nitro groups in the fragment
	 */
	public int nNitros(boolean atomsIncluded[]) {
		int nNO2 = 0;
		for (int i = 0; i < elementAtomicNumber.length; i++) {
			if (elementAtomicNumber[i] == 7) {
				if (atomsIncluded[i]) {
					int nONeigh = 0;
					for (int j = 0; j < edges.length; j++) {
						if (edges[j][0] == i) {
							if (elementAtomicNumber[edges[j][1]] == 8) {
								if (atomsIncluded[edges[j][1]]) {
									nONeigh++;
								}
							}
						}
						if (edges[j][1] == i) {
							if (elementAtomicNumber[edges[j][0]] == 8) {
								if (atomsIncluded[edges[j][0]]) {
									nONeigh++;
								}
							}
						}
					}
					if (nONeigh == 2) {
						nNO2++;
					}
				}
			}
		}
		return nNO2;
	}

	/**
	 * n CH3
	 * 
	 * @param n             - atomic number of element
	 * @param atomsIncluded array have length equal to number of atoms. True if an
	 *                      atom is in fragment. Empty fragment - all false, full
	 *                      mol - all trues. Atom are included (or excluded) with
	 *                      all corresponding implicit hydrogens.
	 * @return number of occurences of the given element in the fragment
	 */
	public int nElement(int n, boolean atomsIncluded[]) {
		int nE = 0;
		for (int i = 0; i < elementAtomicNumber.length; i++) {
			if (nHydrogens[i] == n) {
				if (atomsIncluded[i]) {
					nE++;
				}
			}
		}
		return nE;
	}

	/**
	 * nodes
	 * 
	 * @return number of non-hyderogen atoms (of nodes, hydrogens are not nodes)
	 */
	public int nAtoms() {
		return elementAtomicNumber.length;
	}

	/**
	 * edges
	 * 
	 * @return number of bonds between non-hyderogen atoms (of edges; hydrogens are
	 *         not nodes)
	 */
	public int nBonds() {
		return edges.length;
	}

	/**
	 * Break the bond.
	 * 
	 * @param edgeNum           number of edges (bonds) to destruct
	 * @param invertedDirection if true node2 will carry "charge", i.e. node2 and
	 *                          corresponding component of connectivity with be
	 *                          marked with true in output array
	 * @return returns array with the same length with number of nodes. True if this
	 *         atom is included in fragment ion. All elements are false for cleavage
	 *         of single edge in ring
	 */
	public boolean[] edgeCleavage(int edgeNum, boolean invertedDirection) {
		int nNodes = this.nAtoms();
		int nEdges = this.nBonds();
		boolean[] result = new boolean[nNodes];
		if (!invertedDirection) {
			result[edges[edgeNum][0]] = true;
		} else {
			result[edges[edgeNum][1]] = true;
		}
		for (int c = 0; c < nEdges; c++) {
			for (int i = 0; i < nEdges; i++) {
				if (i != edgeNum) {
					if (result[this.edges[i][0]]) {
						result[this.edges[i][1]] = true;
					}
					if (result[this.edges[i][1]]) {
						result[this.edges[i][0]] = true;
					}
				}
			}
		}

		int markedNodes = 0;
		for (int i = 0; i < nNodes; i++) {
			if (result[i]) {
				markedNodes++;
			}
		}
		if (markedNodes < nNodes) {
			return result;
		} else {
			return new boolean[nNodes];
		}
	}

	/**
	 * Break two bonds
	 * 
	 * @param edgeNum1           number of edge (bond) to destruct
	 * @param edgeNum2           number of edge (bond) to destruct
	 * @param invertedDirection1 if true node2 of edgeNum1 will carry "charge", i.e.
	 *                           node2 and corresponding component of connectivity
	 *                           with be marked with true in output array
	 * @param invertedDirection2 the same for the second edge
	 * @return returns array with the same length with number of nodes. True if this
	 *         atom is included in fragment ion. All elements are false for cleavage
	 *         of 2 edges that do not cause real fragmentation (e.g. two cleavages
	 *         in two different rings) or the second fragmentation fragments
	 *         "neutral loss" (excluded after the first fragmentation part)
	 */
	public boolean[] twoEdgeCleavages(int edgeNum1, int edgeNum2, boolean invertedDirection1,
			boolean invertedDirection2) {
		int nNodes = this.nAtoms();
		int nEdges = this.nBonds();
		boolean[] marked1 = new boolean[nNodes];
		boolean[] marked2 = new boolean[nNodes];

		if (!invertedDirection1) {
			marked1[edges[edgeNum1][0]] = true;
		} else {
			marked1[edges[edgeNum1][1]] = true;
		}

		if (!invertedDirection2) {
			marked2[edges[edgeNum2][0]] = true;
		} else {
			marked2[edges[edgeNum2][1]] = true;
		}

		for (int c = 0; c < nEdges; c++) {
			for (int i = 0; i < nEdges; i++) {
				if (i != edgeNum1) {
					if (i != edgeNum2) {
						if (marked1[this.edges[i][0]]) {
							marked1[this.edges[i][1]] = true;
						}
						if (marked1[this.edges[i][1]]) {
							marked1[this.edges[i][0]] = true;
						}
						if (marked2[this.edges[i][0]]) {
							marked2[this.edges[i][1]] = true;
						}
						if (marked2[this.edges[i][1]]) {
							marked2[this.edges[i][0]] = true;
						}
					}
				}
			}
		}

		int markedNodesBoth = 0;
		for (int i = 0; i < nNodes; i++) {
			if (marked1[i]) {
				if (marked2[i]) {
					markedNodesBoth++;
				}
			}
		}
		int markedNodes1 = 0;
		for (int i = 0; i < nNodes; i++) {
			if (marked1[i]) {
				markedNodes1++;
			}
		}
		if (markedNodesBoth < nNodes) {
			if (markedNodes1 == markedNodesBoth) { // One ion, 1 or 2 neutral losses
				return marked1;
			}
		}
		return new boolean[nNodes];
	}

	/**
	 * Break theww bonds
	 * 
	 * @param edgeNum1           number of edge (bond) to destruct
	 * @param edgeNum2           number of edge (bond) to destruct
	 * @param edgeNum3           number of edge (bond) to destruct
	 * 
	 * @param invertedDirection1 if true node2 of edgeNum1 will carry "charge", i.e.
	 *                           node2 and corresponding component of connectivity
	 *                           with be marked with true in output array
	 * @param invertedDirection2 the same for the second edge
	 * @param invertedDirection2 the same for the third edge
	 * @return returns array with the same length with number of nodes. True if this
	 *         atom is included in fragment ion. All elements are false for cleavage
	 *         of 2 edges that do not cause real fragmentation (e.g. two cleavages
	 *         in two different rings) or the second fragmentation fragments
	 *         "neutral loss" (excluded after the first fragmentation part)
	 */
	public boolean[] threeEdgeCleavages(int edgeNum1, int edgeNum2, int edgeNum3, boolean invertedDirection1,
			boolean invertedDirection2, boolean invertedDirection3) {
		int nNodes = this.nAtoms();
		int nEdges = this.nBonds();
		boolean[] marked1 = new boolean[nNodes];
		boolean[] marked2 = new boolean[nNodes];
		boolean[] marked3 = new boolean[nNodes];

		if (!invertedDirection1) {
			marked1[edges[edgeNum1][0]] = true;
		} else {
			marked1[edges[edgeNum1][1]] = true;
		}

		if (!invertedDirection2) {
			marked2[edges[edgeNum2][0]] = true;
		} else {
			marked2[edges[edgeNum2][1]] = true;
		}

		if (!invertedDirection3) {
			marked3[edges[edgeNum3][0]] = true;
		} else {
			marked3[edges[edgeNum3][1]] = true;
		}

		for (int c = 0; c < nEdges; c++) {
			for (int i = 0; i < nEdges; i++) {
				if (i != edgeNum1) {
					if (i != edgeNum2) {
						if (i != edgeNum3) {
							if (marked1[this.edges[i][0]]) {
								marked1[this.edges[i][1]] = true;
							}
							if (marked1[this.edges[i][1]]) {
								marked1[this.edges[i][0]] = true;
							}
							if (marked2[this.edges[i][0]]) {
								marked2[this.edges[i][1]] = true;
							}
							if (marked2[this.edges[i][1]]) {
								marked2[this.edges[i][0]] = true;
							}
							if (marked3[this.edges[i][0]]) {
								marked3[this.edges[i][1]] = true;
							}
							if (marked3[this.edges[i][1]]) {
								marked3[this.edges[i][0]] = true;
							}
						}
					}
				}
			}
		}

		int markedNodesThree = 0;
		for (int i = 0; i < nNodes; i++) {
			if (marked1[i]) {
				if (marked2[i]) {
					if (marked3[i]) {
						markedNodesThree++;
					}
				}
			}
		}
		int markedNodes1 = 0;
		for (int i = 0; i < nNodes; i++) {
			if (marked1[i]) {
				markedNodes1++;
			}
		}
		if (markedNodesThree < nNodes) {
			if (markedNodes1 == markedNodesThree) { // One ion, 1 or 2 neutral losses
				if (markedNodes1 != markedNodesThree) {
					throw new RuntimeException("Three node cleavage");
				}
				return marked1;
			}
		}
		return new boolean[nNodes];
	}

	private static HashMap<Integer, Integer> clone(HashMap<Integer, Integer> a) {
		HashMap<Integer, Integer> result = new HashMap<Integer, Integer>();
		for (int i : a.keySet()) {
			result.put(i, (a.get(i) + 0));
		}
		return result;
	}

	// array containing possible variants - how many fluorines could migrate into a
	// given fragment during its formation. It contains both positive and negative
	// numbers (if fluorines migrated FROM the fragment). The algorithm takes into
	// account that only fluorines attached to atoms that included the bond that is
	// being broken migrate over this bond, do "long jumps".
	// all numbers in the array are (abs values) less or equal fluorineMaxMigration
	private int[] fluorinePossibleMigration(boolean[] atomsIncluded, int fluorineMaxMigration) {
		int fluorine = 0;
		for (int i = 0; i < nBonds(); i++) {
			int at = -1;
			if (atomsIncluded[edges[i][1]] && (!atomsIncluded[edges[i][0]])) {
				at = edges[i][0];
			}
			if (atomsIncluded[edges[i][0]] && (!atomsIncluded[edges[i][1]])) {
				at = edges[i][1];
			}
			if (at >= 0) { // atoms that have bonds that are breaking
				for (int j = 0; j < nBonds(); j++) {
					if (j != i) {
						if (edges[j][0] == at) {
							if (elementAtomicNumber[edges[j][1]] == 9) {
								fluorine++;
							}
						}
						if (edges[j][1] == at) {
							if (elementAtomicNumber[edges[j][0]] == 9) {
								fluorine++;
							}
						}
					}
				}
				if (elementAtomicNumber[at] == 9) {
					fluorine++;
				}
			}
		}
		int maxFluorine = fluorine;
		if (maxFluorine > fluorineMaxMigration) {
			maxFluorine = fluorineMaxMigration;
		}
		int[] result = new int[maxFluorine];
		for (int i = 0; i < result.length; i++) {
			result[i] = i + 1;
		}
		return result;
	}

	private int nHAll() { // same as nHydrogensAll
		HashMap<Integer, Integer> f = this.molecularFormulaFragment(this.markAllAtoms());
		int nHAll = f.get(1) == null ? 0 : f.get(1);
		return nHAll;
	}

	private static int get(HashMap<Integer, Integer> f, int i) {
		if (f.get(i) == null) {
			return 0;
		}
		return f.get(i);
	}

	/**
	 * The method checks whether a molecular formula is possible in principle. Each
	 * carbon adds no more than 2 hydrogens (CH2 fragment) when added to an existing
	 * one. The same is true for silicon. Nitrogen no more than 1 (when added to an
	 * existing molecule), the same is true for arsenic, phosphorus. If a molecule
	 * contains both fluorine and one of the following elements: sulfur, phosphorus,
	 * arsenic, tellurium, selenium, then True is always returned, since such atoms
	 * can have an unexpectedly large number of bonds in organofluorines
	 * 
	 * @param f Hash table with molecular formula. Key is atomic number of element
	 *          (carbon 6, nitrogen 7, hydrogen 1, etc.), value is number of atoms.
	 *          Hydrogen is also included
	 * @return bool
	 */
	public boolean checkThatNHydrogensAndHalogensIsPossible(HashMap<Integer, Integer> f) {
		if (((get(f, 16) > 0) && (get(f, 9) > 0)) || ((get(f, 15) > 0) && (get(f, 9) > 0))
				|| ((get(f, 33) > 0) && (get(f, 9) > 0)) || ((get(f, 34) > 0) && (get(f, 9) > 0))
				|| ((get(f, 52) > 0) && (get(f, 9) > 0))) { // S, P, As, Te, Se can have many X-F bonds. No strict
															// limits are posed in such cases.
			return true;
		}

		int n = 2;
		n += 2 * get(f, 6);
		n += get(f, 5) + get(f, 7);
		n += 2 * get(f, 14) + get(f, 15) + get(f, 33);
		n -= get(f, 1);
		n -= get(f, 9);
		n -= get(f, 17);
		n -= get(f, 35);
		n -= get(f, 53);
		return n >= 0;
	}

	/**
	 * The method allows to obtain possible fragment ions for a given fragment
	 * taking into account additional migrations. Without additional rearrangements
	 * - one ion corresponds to a single fragment, see the molecularFormulaFragment
	 * method. However, additional migrations are considered here: loss or
	 * acquisition (migration both INTO the fragment and FROM it) of each element of
	 * possibleHChangeFragment hydrogens, migration of no more than
	 * maxFluorineMigration fluorines through the bond being broken. In this case,
	 * fluorines migrate "locally" (through the bond being broken), and hydrogens
	 * can migrate "globally" so that their total amount does not exceed that of the
	 * original molecule. Additionally, the loss of no more than one halogen
	 * (chlorine bromine iodine) and no more than one methyl group is considered.
	 * Two or more different halogens or halogens and methyl can be lost
	 * simultaneously. For example, methyl, chlorine iodine can "fall off"
	 * simultaneously. But two methyls or two chlorines cannot. Another
	 * rearrangement is the transformation of the nitro group into oxygen. Moreover,
	 * if there are 2 or more nitro groups, each can turn into oxygen. As a result,
	 * the method returns ALL possible options that can be obtained from a given
	 * fragment using the above-described rearrangements.
	 * 
	 * @param atomsIncluded array with the same length with number of nodes. True if
	 *                      this atom is included in fragment ion.
	 * @return list of possible fragments
	 */
	public ArrayList<FragmentIon> formulasForFragment(boolean atomsIncluded[]) {
		boolean empty = true;
		for (int i = 0; i < atomsIncluded.length; i++) {
			if (atomsIncluded[i]) {
				empty = false;
			}
		}
		if (empty) {
			return new ArrayList<FragmentIon>();
		}
		ArrayList<FragmentIon> result = new ArrayList<FragmentIon>();
		HashMap<Integer, Integer> f = this.molecularFormulaFragment(atomsIncluded);
		int nHAll = nHAll();
		result.add(new FragmentIon(f));
		for (int element : new int[] { 17, 35, 53 }) {
			if (f.get(element) != null) {
				if (f.get(element) > 0) {
					HashMap<Integer, Integer> f1 = clone(f);
					f1.put(element, f.get(element) - 1);
					FragmentIon fi = new FragmentIon(f1);
					if (fi.elementNumbers.length != 0) {
						result.add(fi);
					}
				}
			}
		}

		if (this.nMethyls(atomsIncluded) > 0) {
			HashMap<Integer, Integer> f1 = clone(f);
			if ((f1.get(6) < 1) || (f1.get(1) < 3)) {
				throw new RuntimeException("Unknown error. Should be methyl group here!");
			}
			if (!f1.get(6).equals(1)) {
				f1.put(6, f.get(6) - 1);
				f1.put(1, f.get(1) - 3);
			}
			result.add(new FragmentIon(f1));
		}
		int nNitro = this.nNitros(atomsIncluded);
		if (nNitro > 0) {
			HashMap<Integer, Integer> f1 = clone(f);
			if ((f1.get(7) < nNitro) || (f1.get(8) < 2 * nNitro)) {
				throw new RuntimeException("Unknown error. Should be nitro groups here!");
			}
			if (f1.get(6) != null) {
				if (!f1.get(6).equals(0)) {
					for (int j = 1; j < nNitro + 1; j++) {
						f1.put(7, f.get(7) - j);
						f1.put(8, f.get(8) - j);
					}
				}
			}
			result.add(new FragmentIon(f1));
		}
		@SuppressWarnings("unchecked")
		ArrayList<FragmentIon> tmp = (ArrayList<FragmentIon>) result.clone();
		for (FragmentIon a : tmp) {
			HashMap<Integer, Integer> f2 = a.asMap();
			int nH = 0;
			if (f2.get(1) != null) {
				nH = f2.get(1);
			}
			for (int dH : this.possibleHChangeFragment) {
				if ((nH + dH >= 0) && (nH + dH <= nHAll)) {
					HashMap<Integer, Integer> f1 = clone(f2);
					f1.put(1, nH + dH);
					if (checkThatNHydrogensAndHalogensIsPossible(f1)) {
						result.add(new FragmentIon(f1));
					}
				}
			}
		}
		tmp = null;
		int[] possibleMigrationF = this.fluorinePossibleMigration(atomsIncluded, this.maxFluorineMigration);
		@SuppressWarnings("unchecked")
		ArrayList<FragmentIon> tmp1 = (ArrayList<FragmentIon>) result.clone();
		for (FragmentIon a : tmp1) {
			HashMap<Integer, Integer> f2 = a.asMap();
			int nF = 0;
			if (f2.get(9) != null) {
				nF = f2.get(9);
			}
			for (int dF : possibleMigrationF) {
				if (nF + dF >= 0) {
					HashMap<Integer, Integer> f1 = clone(f2);
					f1.put(9, nF + dF);
					if (checkThatNHydrogensAndHalogensIsPossible(f1)) {
						result.add(new FragmentIon(f1));
					}
				}
			}
			for (int dF = -this.maxFluorineMigration; dF < 0; dF++) {
				if (nF + dF >= 0) {
					HashMap<Integer, Integer> f1 = clone(f2);
					f1.put(9, nF + dF);
					FragmentIon fi = new FragmentIon(f1);
					if (fi.elementNumbers.length != 0) {
						result.add(fi);
					}
				}
			}
		}
		return result;
	}

	/**
	 * Full fragment
	 * 
	 * @return array of Trues with same length as number of atoms
	 */
	public boolean[] markAllAtoms() {
		int nNodes = this.nAtoms();
		boolean[] marked = new boolean[nNodes];
		for (int i = 0; i < marked.length; i++) {
			marked[i] = true;
		}
		return marked;
	}

	/**
	 * This method calculates all possible molecular formulas that can theoretically
	 * be formed from the molecular formula of a molecule. It does not take into
	 * account the structure and bond breaks, as well as any laws on the topic of
	 * parity. This is the maximum, most general list.
	 * 
	 * @return fragment ions
	 */
	public ArrayList<FragmentIon> allPossibleFragmentsFormulas() {
		ArrayList<HashMap<Integer, Integer>> result = new ArrayList<HashMap<Integer, Integer>>();
		HashMap<Integer, Integer> f = this.molecularFormulaFragment(this.markAllAtoms());
		for (int i : f.keySet()) {
			if (result.size() == 0) {
				for (int j = 0; j < f.get(i) + 1; j++) {
					HashMap<Integer, Integer> f1 = new HashMap<Integer, Integer>();
					f1.put(i, j);
					result.add(f1);
				}
			} else {
				ArrayList<HashMap<Integer, Integer>> newList = new ArrayList<HashMap<Integer, Integer>>();
				for (int j = 0; j < f.get(i) + 1; j++) {
					for (HashMap<Integer, Integer> x : result) {
						HashMap<Integer, Integer> f1 = clone(x);
						f1.put(i, j);
						newList.add(f1);
					}
				}
				result.addAll(newList);
			}
		}
		ArrayList<FragmentIon> resultClean = new ArrayList<FragmentIon>();
		for (HashMap<Integer, Integer> f1 : result) {
			if (checkThatNHydrogensAndHalogensIsPossible(f1)) {
				FragmentIon fi = new FragmentIon(f1);
				if (fi.elements.length != 0) {
					resultClean.add(fi);
				}
			}
		}
		return resultClean;
	}

	/**
	 * This method generates all possible fragment ions taking into account the
	 * number of possible hydrogen and fluorine migrations and the number of bond
	 * breaks (cleavages). In this case, all possible bond breaks are performed, and
	 * for each resulting fragment, the method formulasForFragment is applied.
	 * Hydrogen migrations are non-local. The total number of hydrogens that can
	 * enter or leave a fragment is specified. Fluorine migrations are local - only
	 * over 1 bond;
	 * 
	 * @param nCleavages           number of bond breaks. 1 2 3 allowed! Other
	 *                             numbers are not allowed
	 * @param maxHDrift            maximum number of hydrogens that can migrate to
	 *                             fragemnt
	 * @param maxHLoss             maximum number of hydrogens that can migrate from
	 *                             fragemnt
	 * @param maxFluorineMigration fluorine migration
	 * @return list of fragments.
	 */
	public ArrayList<FragmentIon> fragmentsFormulas(int nCleavages, int maxHDrift, int maxHLoss,
			int maxFluorineMigration) {
		this.maxFluorineMigration = maxFluorineMigration;
		this.possibleHChangeFragment = new int[maxHDrift + maxHLoss];
		int j1 = 0;
		for (int i = 0 - maxHLoss; i < maxHDrift + 1; i++) {
			if (i != 0) {
				this.possibleHChangeFragment[j1] = i;
				j1++;
			}
		}
		if ((nCleavages < 1) || (nCleavages > 3)) {
			throw new RuntimeException("1-3 bond cleavages are expected");
		}
		ArrayList<FragmentIon> result = new ArrayList<FragmentIon>();
		result.addAll(formulasForFragment(markAllAtoms()));
		for (int i = 0; i < this.edges.length; i++) {
			result.addAll(formulasForFragment(this.edgeCleavage(i, false)));
			result.addAll(formulasForFragment(this.edgeCleavage(i, true)));
		}
		if (nCleavages >= 2) {
			for (int i = 0; i < this.edges.length; i++) {
				for (int j = 0; j < this.edges.length; j++) {
					result.addAll(formulasForFragment(this.twoEdgeCleavages(i, j, false, false)));
					result.addAll(formulasForFragment(this.twoEdgeCleavages(i, j, false, true)));
					result.addAll(formulasForFragment(this.twoEdgeCleavages(i, j, true, false)));
					result.addAll(formulasForFragment(this.twoEdgeCleavages(i, j, true, true)));
				}
			}
		}
		if (nCleavages >= 3) {
			for (int i = 0; i < this.edges.length; i++) {
				for (int j = 0; j < this.edges.length; j++) {
					for (int k = 0; k < this.edges.length; k++) {
						result.addAll(formulasForFragment(this.threeEdgeCleavages(i, j, k, false, false, false)));
						result.addAll(formulasForFragment(this.threeEdgeCleavages(i, j, k, false, false, true)));
						result.addAll(formulasForFragment(this.threeEdgeCleavages(i, j, k, false, true, false)));
						result.addAll(formulasForFragment(this.threeEdgeCleavages(i, j, k, false, true, true)));
						result.addAll(formulasForFragment(this.threeEdgeCleavages(i, j, k, true, false, false)));
						result.addAll(formulasForFragment(this.threeEdgeCleavages(i, j, k, true, false, true)));
						result.addAll(formulasForFragment(this.threeEdgeCleavages(i, j, k, true, true, false)));
						result.addAll(formulasForFragment(this.threeEdgeCleavages(i, j, k, true, true, true)));
					}
				}
			}
		}
		return result;
	}

	private static int nHydrogensNeigh(IAtom a) {
		int count = 0;
		for (IBond b : a.bonds()) {
			if (b.getOther(a).getAtomicNumber().equals(1)) {
				count = count + 1;
			}
		}
		return count + a.getImplicitHydrogenCount();
	}

	/**
	 * From smiles
	 * @param smiles SMILES string, structure
	 */
	public MolGraph(String smiles) {
		IAtomContainer mol = null;
		try {
			SmilesParser parser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
			mol = parser.parseSmiles(smiles.trim());
			Aromaticity arom = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
			arom.apply(mol);
			CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(mol.getBuilder());
			adder.addImplicitHydrogens(mol);
		} catch (CDKException e) {
			e.printStackTrace();
			throw new RuntimeException(e.getMessage());
		}

		ArrayList<Integer> atomsElementAtomicNumber = new ArrayList<Integer>();
		ArrayList<Integer> atomsNHydrogensNeigh = new ArrayList<Integer>();
		HashMap<Integer, Integer> graphNumberByCDKNumber = new HashMap<Integer, Integer>();
		ArrayList<Integer> firstNodesOfBonds = new ArrayList<Integer>();
		ArrayList<Integer> secondNodesOfBonds = new ArrayList<Integer>();

		int heavyAtoms = 0;
		for (IAtom a : mol.atoms()) {
			if (!a.getAtomicNumber().equals(1)) {
				heavyAtoms++;
				int n = a.getIndex();
				graphNumberByCDKNumber.put(n, atomsElementAtomicNumber.size());
				atomsElementAtomicNumber.add(a.getAtomicNumber());
				atomsNHydrogensNeigh.add(nHydrogensNeigh(a));
			}
		}
		if (heavyAtoms != atomsElementAtomicNumber.size()) {
			throw new RuntimeException("Unknown error");
		}
		this.elementAtomicNumber = new int[heavyAtoms];
		this.nHydrogens = new int[heavyAtoms];
		for (int i = 0; i < heavyAtoms; i++) {
			elementAtomicNumber[i] = atomsElementAtomicNumber.get(i);
			nHydrogens[i] = atomsNHydrogensNeigh.get(i);
		}
		int nonHbonds = 0;
		for (IBond b : mol.bonds()) {
			boolean nonH = true;
			for (IAtom a : b.atoms()) {
				if (a.getAtomicNumber().equals(1)) {
					nonH = false;
				}
			}
			nonHbonds = nonH ? nonHbonds + 1 : nonHbonds;
			if (nonH) {
				int a1 = graphNumberByCDKNumber.get(b.getBegin().getIndex());
				int a2 = graphNumberByCDKNumber.get(b.getEnd().getIndex());
				int a11 = a1 < a2 ? a1 : a2;
				int a21 = a1 < a2 ? a2 : a1;
				if (a1 == a2) {
					throw new RuntimeException("Unknown error");
				}
				firstNodesOfBonds.add(a11);
				secondNodesOfBonds.add(a21);
			}
		}
		if (nonHbonds != firstNodesOfBonds.size()) {
			throw new RuntimeException("Unknown error");
		}
		this.edges = new int[nonHbonds][];
		for (int i = 0; i < this.edges.length; i++) {
			this.edges[i] = new int[2];
			this.edges[i][0] = firstNodesOfBonds.get(i);
			this.edges[i][1] = secondNodesOfBonds.get(i);
		}
	}

	/**
	 * 
	 * @return Possible variations in the number of hydrogens during the formation of a fragment ion
	 */
	public int[] getPossibleHChangeFragment() {
		return possibleHChangeFragment;
	}

	/**
	 * 
	 * @param possibleHChangeFragment Possible variations in the number of hydrogens during the formation of a fragment ion
	 */
	public void setPossibleHChangeFragment(int[] possibleHChangeFragment) {
		this.possibleHChangeFragment = possibleHChangeFragment;
	}

	/**
	 * 
	 * @return max number fluorine migration
	 */
	public int getMaxFluorineMigration() {
		return maxFluorineMigration;
	}

	/**
	 * 
	 * @param maxFluorineMigration max number fluorine migration
	 */
	public void setMaxFluorineMigration(int maxFluorineMigration) {
		this.maxFluorineMigration = maxFluorineMigration;
	}

}
