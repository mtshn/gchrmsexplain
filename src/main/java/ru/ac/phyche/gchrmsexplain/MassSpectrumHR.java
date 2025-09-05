package ru.ac.phyche.gchrmsexplain;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import org.apache.commons.io.input.BOMInputStream;
import org.apache.commons.lang3.tuple.Pair;
import org.openscience.cdk.exception.CDKException;

/**
 * High resolution mass spectrum and peak explanation
 */
public class MassSpectrumHR {
	private float[] mzs = null;
	private float[] intensities = null;

	public float thresholdGenerateIsotopic = 1E-3f;

	private static class ComparableF implements Comparable<ComparableF> {
		public float[] x = null;

		@Override
		public int compareTo(ComparableF o) {
			return ((Float) x[0]).compareTo(o.x[0]);
		}

		public ComparableF(float[] x) {
			this.x = x;
		}
	}

	/**
	 * The method loads a spectrum from a CSV (comma separated) file with header
	 * Scan Number,m/z,Intensity,Relative,Segment Number,Flags. Everything above the
	 * header is ignored. m/z and relative intensities are in the corresponding
	 * columns (1 and 3, starting with 0). If the maximum intensity in the original
	 * data is not 999, then the peaks are rescaled so that the maximum peak is 999
	 * 
	 * @param filename                  file name
	 * @param thresholdBy999            Threshold, peaks with intensity below which
	 *                                  are ignored, we always consider that the
	 *                                  maximum peak has intensity 999.
	 * @param thresholdGenerateIsotopic The intensity threshold used when generating
	 *                                  the isotope distribution. Peaks with an
	 *                                  intensity below the threshold are not
	 *                                  generated. If you do not plan to use
	 *                                  non-static methods of this class, which
	 *                                  include consideration/generation of the
	 *                                  isotope distribution, then you can specify
	 *                                  any number. Important! The threshold is
	 *                                  calculated based on the total intensity of
	 *                                  ALL isotope peaks, which is taken to be
	 *                                  equal to 1 (not the base peak, but the
	 *                                  sum!!!)
	 * @return Mass spectrum
	 * @throws IOException io
	 */
	public static MassSpectrumHR fromThermoCSVFile(String filename, float thresholdBy999,
			float thresholdGenerateIsotopic) throws IOException {
		try {
			return fromAnyCSVFile(filename, thresholdBy999, thresholdGenerateIsotopic, 1, 3,
					"Scan Number,m/z,Intensity,Relative,Segment Number,Flags", ',');
		} catch (Throwable e) {
			return fromAnyCSVFile(filename, thresholdBy999, thresholdGenerateIsotopic, 1, 2,
					"Scan Number,m/z,Relative,Segment Number,Flags", ',');
		}
	}

	/**
	 * 
	 * The method loads a spectrum from a CSV (comma separated) file with header
	 * m/z,Intensity. Everything above the header is ignored. m/z and relative
	 * intensities are in the corresponding columns (0 and 1)
	 * 
	 * @param filename                  file name
	 * @param thresholdBy999            Threshold, peaks with intensity below which
	 *                                  are ignored, we always consider that the
	 *                                  maximum peak has intensity 999.
	 * @param thresholdGenerateIsotopic The intensity threshold used when generating
	 *                                  the isotope distribution. Peaks with an
	 *                                  intensity below the threshold are not
	 *                                  generated. If you do not plan to use
	 *                                  non-static methods of this class, which
	 *                                  include consideration/generation of the
	 *                                  isotope distribution, then you can specify
	 *                                  any number. Important! The threshold is
	 *                                  calculated based on the total intensity of
	 *                                  ALL isotope peaks, which is taken to be
	 *                                  equal to 1 (not the base peak, but the
	 *                                  sum!!!)
	 * @return Mass spectrum
	 * @throws IOException io
	 */
	public static MassSpectrumHR fromSimpleCSVFile(String filename, float thresholdBy999,
			float thresholdGenerateIsotopic) throws IOException {
		return fromAnyCSVFile(filename, thresholdBy999, thresholdGenerateIsotopic, 0, 1, "m/z,Intensity", ',');
	}

	private interface MZtester {
		boolean test(float mz);
	}

	/**
	 * This method removes from THIS spectrum all peaks larger than the molecular
	 * ion, except for specially allowed ones. Isotopic peaks of the molecular ion
	 * are also not removed.
	 * 
	 * @param smiles           proposed structure of the molecule (SMILES)
	 * @param permittedAdducts list of m/z's that are allowed and not removed, even
	 *                         if they are greater than the molecular m/z.
	 * @param permittedMZ      List of exact masses that can be added to a molecular
	 *                         ion (and are not removed by this function) For
	 *                         example, if the m/z of a molecular ion is M and OH is
	 *                         added to it (with mass M_OH), then the ion M+M_OH is
	 *                         not removed. M_OH must be in this array! The isotopic
	 *                         distribution of the molecular ion is taken into
	 *                         account
	 * @param mzThreshold      mass determination accuracy in high resolution mass
	 *                         spectrometry
	 */
	public void cutMZAboveMolecularIon(String smiles, float[] permittedAdducts, float[] permittedMZ,
			float mzThreshold) {
		MolGraph m = new MolGraph(smiles);
		HashMap<Integer, Integer> mf = m.molecularFormulaFragment(m.markAllAtoms());
		IsotopicPattern ip = (new FragmentIon(mf)).isotopicPatternSingleCation(thresholdGenerateIsotopic);
		float maxmz = 0;
		for (float mz : ip.asMap().keySet()) {
			if (maxmz < mz) {
				maxmz = mz;
			}
		}
		ArrayList<Float> permittedMZFullList = new ArrayList<Float>();
		for (float mz : permittedMZ) {
			permittedMZFullList.add(mz);
		}
		for (float mz : ip.asMap().keySet()) {
			for (float mz1 : permittedAdducts) {
				permittedMZFullList.add(mz + mz1);
			}
		}
		final float maxmz2 = maxmz;
		MZtester mzTester = new MZtester() {
			public boolean test(float mz) {
				if (mz < (maxmz2 + 1.5f)) {
					return true;
				}
				for (float mz1 : permittedMZFullList) {
					if (Math.abs(mz - mz1) < mzThreshold) {
						return true;
					}
				}
				return false;
			}
		};
		ArrayList<Float> acceptedMZs = new ArrayList<Float>();
		ArrayList<Float> acceptedIntensities = new ArrayList<Float>();
		for (int i = 0; i < this.mzs.length; i++) {
			if (mzTester.test(mzs[i])) {
				acceptedMZs.add(mzs[i]);
				acceptedIntensities.add(intensities[i]);
			}
		}
		this.mzs = new float[acceptedMZs.size()];
		this.intensities = new float[acceptedMZs.size()];
		for (int i = 0; i < this.mzs.length; i++) {
			mzs[i] = acceptedMZs.get(i);
			intensities[i] = acceptedIntensities.get(i);
		}
	}

	/**
	 * 
	 * The method loads a spectrum from a CSV (comma separated) file with header.
	 * Everything above the header is ignored. m/z and relative intensities are in
	 * the corresponding columns (1 and 3, starting with 0). If the maximum
	 * intensity in the original data is not 999, then the peaks are rescaled so
	 * that the maximum peak is 999
	 * 
	 * @param filename                  file name
	 * @param header                    Header string, everything above is ignored.
	 *                                  For example "m/z,Intensity". Header string
	 *                                  is obligatory.
	 * @param thresholdBy999            Threshold, peaks with intensity below which
	 *                                  are ignored, we always consider that the
	 *                                  maximum peak has intensity 999.
	 * @param thresholdGenerateIsotopic See other similar methods. If you do not
	 *                                  plan to use non-static methods of this
	 *                                  class, which include
	 *                                  consideration/generation of the isotope
	 *                                  distribution, then you can specify any
	 *                                  number.
	 * @return Mass spectrum
	 * @throws IOException io
	 */
	public static MassSpectrumHR fromSimpleCSVFile(String filename, String header, float thresholdBy999,
			float thresholdGenerateIsotopic) throws IOException {
		return fromAnyCSVFile(filename, thresholdBy999, thresholdGenerateIsotopic, 0, 1, header, ',');
	}

	/**
	 * @return Mass spectrum as space-separated string
	 */
	public String toString() {
		String result = "";
		for (int i = 0; i < this.mzs.length; i++) {
			result = result + mzs[i] + " " + intensities[i] + " ";
		}
		return result.trim();
	}

	/**
	 * A string containing m/z values, intensities, and no other words or numbers.
	 * The string must contain N m/z-intensity pairs separated by any number of
	 * semicolons, colons, commas, parentheses (i.e., ")", "("), tabs, newlines, and
	 * spaces. For example, "80.1 999 81.0 48", where 999 and 48 are intensities. If
	 * the maximum intensity in the original data is not 999, then the peaks are
	 * rescaled so that the maximum peak is 999
	 * 
	 * @param line                      A string with m/z and intensity pairs
	 * @param thresholdBy999            Threshold, peaks with intensity below which
	 *                                  are ignored, we always consider that the
	 *                                  maximum peak has intensity 999.
	 * @param thresholdGenerateIsotopic See method fromSimpleCSVFile and other
	 * @return mass spectrum
	 */
	public static MassSpectrumHR fromSimpleMZTable(String line, float thresholdBy999, float thresholdGenerateIsotopic) {
		line = line.trim().replace("(", " ").replace(")", " ").replace(";", " ").replace(",", " ").replace(":", " ")
				.trim();
		String[] numbers = line.split("\\s+");
		ArrayList<ComparableF> r1 = new ArrayList<ComparableF>();
		float maxI = 0;
		for (int i = 0; i < numbers.length; i += 2) {
			float mz = Float.parseFloat(numbers[i]);
			float intens = Float.parseFloat(numbers[i + 1]);
			r1.add(new ComparableF(new float[] { mz, intens }));
			if (intens > maxI) {
				maxI = intens;
			}
		}

		ArrayList<ComparableF> r = new ArrayList<ComparableF>();
		for (ComparableF a : r1) {
			if (a.x[1] * 999 / maxI >= thresholdBy999) {
				r.add(a);
			}
		}
		ComparableF[] r2 = r.toArray(new ComparableF[r.size()]);
		Arrays.sort(r2);
		MassSpectrumHR result = new MassSpectrumHR();
		result.mzs = new float[r.size()];
		result.intensities = new float[r.size()];
		int j = 0;
		for (ComparableF a : r2) {
			result.mzs[j] = a.x[0];
			result.intensities[j] = a.x[1] * 999 / maxI;
			j++;
		}
		return result;
	}

	/**
	 * The same as fromSimpleMZTable but the line is loaded from file filename (all
	 * file - is the line, one spectrum in the file)
	 * 
	 * @param filename                  file name
	 * @param thresholdBy999            see description of the method
	 *                                  fromSimpleCSVFile
	 * @param thresholdGenerateIsotopic see description of the method
	 *                                  fromSimpleCSVFile
	 * @return mass spectrum
	 * @throws IOException io
	 */
	public static MassSpectrumHR fromSimpleMZTableFile(String filename, float thresholdBy999,
			float thresholdGenerateIsotopic) throws IOException {
		BufferedReader br = new BufferedReader(
				new InputStreamReader(new BOMInputStream(new FileInputStream(filename))));
		String s = br.readLine();
		String merged = "";
		while (s != null) {
			merged += " " + s;
			s = br.readLine();
		}
		br.close();
		return fromSimpleMZTable(merged, thresholdBy999, thresholdGenerateIsotopic);
	}

	private static String trim(String s) {
		if (s != null) {
			return s.trim();
		} else {
			return null;
		}
	}

	private static String removeKey(String s, char separator) {
		String result = "";
		boolean go = false;
		for (char c : s.toCharArray()) {
			if (go) {
				result += c;
			}
			if (c == separator) {
				go = true;
			}
		}
		return result;
	}

	/**
	 * Loads mass spectrum from MSP file. It is assumed that the file contains a
	 * single spectrum, otherwise the first spectrum is loaded. More detailed
	 * explanation in the fromMSPBlock method.
	 * 
	 * @param filename                  file name
	 * @param thresholdBy999            see description of the method
	 *                                  fromSimpleCSVFile
	 * @param thresholdGenerateIsotopic see description of the method
	 *                                  fromSimpleCSVFile
	 * @return The first member of the returned tuple is the spectrum. The second
	 *         array of strings of length 2 contains the molecule name and the
	 *         SMILES string (structure), in that order.
	 * @throws IOException io
	 */
	public static Pair<MassSpectrumHR, String[]> fromMSPFile(String filename, float thresholdBy999,
			float thresholdGenerateIsotopic) throws IOException {
		BufferedReader br = new BufferedReader(
				new InputStreamReader(new BOMInputStream(new FileInputStream(filename))));
		String s = br.readLine();
		String msp = "";
		boolean go = true;
		while (go) {
			msp += s + "\n";
			s = br.readLine();
			if ((s == null) || (s.toUpperCase().split(":")[0].equals("NAME"))) {
				go = false;
			}
		}
		Pair<MassSpectrumHR, String[]> sp = null;
		sp = MassSpectrumHR.fromMSPBlock(msp, thresholdBy999, thresholdGenerateIsotopic);
		br.close();
		return sp;
	}

	/**
	 * 
	 * Parsing loads mass spectrum from MSP block. Several options are supported for
	 * how the structure can be specified. A separate field containing the word
	 * "SMILES" in the name, for example just "SMILES" or "SMILES_can", the InChI
	 * field, or the Comments field in which (mass bank style) comments in quotes
	 * separated by spaces, for example in the comments field such fragments
	 * "SMILES=CC=C", "SMILES_can=CC=C", "InChI=InChI=1S/C3H6/c1-3-2/h3H,1H2,2H3".
	 * The block must begin with the Name field, the last field before the peak list
	 * must be Num peaks. See also the documentation for NIST mssearch. If the
	 * maximum intensity in the original data is not 999, then the peaks are
	 * rescaled so that the maximum peak is 999
	 * 
	 * @param msp                       msp block (contain newlines)
	 * @param thresholdBy999            Threshold, peaks with intensity below which
	 *                                  are ignored, we always consider that the
	 *                                  maximum peak has intensity 999.
	 * @param thresholdGenerateIsotopic See other similar methods. If you do not
	 *                                  plan to use non-static methods of this
	 *                                  class, which include
	 *                                  consideration/generation of the isotope
	 *                                  distribution, then you can specify any
	 *                                  number.
	 * 
	 * @return The first member of the returned tuple is the spectrum. The second
	 *         array of strings of length 2 contains the molecule name and the
	 *         SMILES string (structure), in that order.
	 * @throws IOException
	 */
	public static Pair<MassSpectrumHR, String[]> fromMSPBlock(String msp, float thresholdBy999,
			float thresholdGenerateIsotopic) throws IOException {

		String[] splt = msp.trim().split("(\\r\\n|\\n|;)");
		HashMap<String, String> fields = new HashMap<String, String>();
		if (!splt[0].toUpperCase().split(":")[0].toLowerCase().equals("name")) {
			throw new RuntimeException(
					"Invalid MSP. First line should be \"Name: ...\" " + splt[0] + " " + splt.length);
		}
		int j = 0;
		while (!splt[j].toLowerCase().split(":")[0].equals("num peaks")) {
			String[] splt2 = splt[j].trim().split(":");
			fields.put(splt2[0].toLowerCase().trim(), removeKey(splt[j], ':'));
			j++;
			if (j == splt.length) {
				throw new RuntimeException("Invalid MSP. No \"Num peaks: ...\" field " + splt[0] + " " + splt.length);
			}
		}
		String[] result2 = new String[2];
		result2[0] = trim(fields.get("name"));
		for (String s2 : fields.keySet()) {
			if (s2.contains("smiles")) {
				result2[1] = fields.get(s2.trim());
			}
		}
		if (result2[1] == null) {// no SMILES, check INCHI
			for (String s2 : fields.keySet()) {
				if (s2.equals("inchi")) {
					try {
						result2[1] = fields.get(s2.trim());
						result2[1] = ParsingNIST23.inchiToSmiles(result2[1], false);
					} catch (CDKException e) {
					}
				}
			}
		}
		if (result2[1] == null) {// no SMILES, check comment
			for (String s2 : fields.keySet()) {
				if (s2.equals("comments")) {
					String comments = fields.get(s2.trim()).trim();

					// Mass bank-style comment in MSP file.
					if ((comments.charAt(0) == '"') && (comments.charAt(comments.length() - 1) == '"')) {
						String[] split3 = ("\" " + comments + " \"").split("\"\\s+\"");
						for (String x : split3) {
							if (x.split("=")[0].toLowerCase().contains("smiles")) {
								result2[1] = removeKey(x, '=');
							}
						}
						if (result2[1] == null) {// no SMILES, search INCHI
							for (String x : split3) {
								if (x.split("=")[0].toLowerCase().equals("inchi")) {
									try {
										result2[1] = removeKey(x, '=');
										result2[1] = ParsingNIST23.inchiToSmiles(result2[1], false);
									} catch (CDKException e) {
									}
								}
							}
						}
					}
				}
			}
		}
		int numpeaks = Integer.parseInt(splt[j].toUpperCase().split(":")[1].trim());
		j++;
		String spectrum = "";
		while (j < splt.length) {
			spectrum = spectrum += " " + splt[j];
			j++;
		}
		spectrum = spectrum.replace("(", "").replace(")", "").replace(",", "").replace(";", "").trim();
		if (spectrum.split("\\s+").length != 2 * numpeaks) {
			throw new RuntimeException("Wrong number of peaks." + msp);
		}
		MassSpectrumHR result1 = fromSimpleMZTable(spectrum, thresholdBy999, thresholdGenerateIsotopic);
		return Pair.of(result1, result2);
	}

	/**
	 * 
	 * The method loads a spectrum from a CSV (comma separated) file with given
	 * header and given numbers of columns in the table with m/z and intensity.
	 * Everything above the header is ignored. m/z and relative intensities are in
	 * the corresponding columns (column numbers start from zero).
	 * 
	 * @param filename                  file name
	 * @param thresholdBy999            Threshold, peaks with intensity below which
	 *                                  are ignored, we always consider that the
	 *                                  maximum peak has intensity 999.
	 * @param thresholdGenerateIsotopic The intensity threshold used when generating
	 *                                  the isotope distribution. Peaks with an
	 *                                  intensity below the threshold are not
	 *                                  generated. If you do not plan to use
	 *                                  non-static methods of this class, which
	 *                                  include consideration/generation of the
	 *                                  isotope distribution, then you can specify
	 *                                  any number. Important! The threshold is
	 *                                  calculated based on the total intensity of
	 *                                  ALL isotope peaks, which is taken to be
	 *                                  equal to 1 (not the base peak, but the
	 *                                  sum!!!)
	 * @param csvColumnMZ               number of column in table with MZs
	 * @param csvColumnIntens           number of column in table with Intensities
	 * @param header                    header
	 * @param separator                 separator (always "," is used)
	 * @return Mass spectrum
	 * @throws IOException io
	 */
	public static MassSpectrumHR fromAnyCSVFile(String filename, float thresholdBy999, float thresholdGenerateIsotopic,
			int csvColumnMZ, int csvColumnIntens, String header, char separator) throws IOException {
		BufferedReader br = new BufferedReader(
				new InputStreamReader(new BOMInputStream(new FileInputStream(filename))));
		String s = br.readLine();
		while (!s.trim().trim().equals(header)) {
			s = br.readLine();
			if (s == null) {
				br.close();
				throw new RuntimeException("Incorrect CSV spectrum. Header not found!");
			}
		}
		s = br.readLine();
		ArrayList<ComparableF> r1 = new ArrayList<ComparableF>();
		float maxI = 0;
		while (s != null) {
			if (!s.trim().equals("")) {
				String[] splt = App.splitCSV(s.trim(), separator);
				float mz = Float.parseFloat(splt[csvColumnMZ]);
				float intens = Float.parseFloat(splt[csvColumnIntens]);
				if (intens > maxI) {
					maxI = intens;
				}
				r1.add(new ComparableF(new float[] { mz, intens }));

			}
			s = br.readLine();
		}
		br.close();
		ArrayList<ComparableF> r = new ArrayList<ComparableF>();
		for (ComparableF a : r1) {
			if (a.x[1] * 999 / maxI >= thresholdBy999) {
				r.add(a);
			}
		}
		ComparableF[] r2 = r.toArray(new ComparableF[r.size()]);
		Arrays.sort(r2);
		MassSpectrumHR result = new MassSpectrumHR();
		result.mzs = new float[r.size()];
		result.intensities = new float[r.size()];
		int j = 0;
		for (ComparableF a : r2) {
			result.mzs[j] = a.x[0];
			result.intensities[j] = a.x[1] * 999 / maxI;
			j++;
		}
		return result;
	}

	/**
	 * Rounding to low resolution. Rounding borders [n-0.38 ... n+ 0.62] ->n See
	 * Khrisanfov M, Samokhin A. A general procedure for rounding m/z values in
	 * low‐resolution mass spectra. Rapid Communications in Mass Spectrometry. 2022
	 * Jun 15;36(11):e9294. 10.1002/rcm.9294
	 * 
	 * @param threshold
	 * @return low res mass spectrum
	 */
	public MassSpectrumLR toLowResolution(float threshold) {
		float[] mzInt = new float[10000];
		for (int i = 0; i < this.mzs.length; i++) {
			for (int j = 0; j < mzInt.length; j++) {
				if (mzs[i] > j - 0.38) {
					if (mzs[i] <= j + 0.62) {
						mzInt[j] += intensities[i];
					}
				}
			}
		}
		float max = -10000;
		for (int i = 0; i < mzInt.length; i++) {
			if (mzInt[i] > max) {
				max = mzInt[i];
			}
		}
		String sp = "";
		for (int i = 0; i < mzInt.length; i++) {
			if (mzInt[i] * 999 / max > threshold) {
				sp += i + " " + mzInt[i] * 999 / max + " ";
			}
		}
		sp = sp.trim();
		return new MassSpectrumLR(sp, false, true);
	}

	/**
	 * Result of mass spectrum peak explanation (taking into account isotopic
	 * peaks). Explanation levels: 0 - exact m/z contradicts any possible fragment;
	 * 1 - exact m/z corresponds to a possible fragment; 2 - isotopic distribution
	 * does not contradict the observed data, but some isotopic peaks merged with
	 * other peaks; 3 - all isotopic peaks are observed perfectly. See also
	 * explainPeaks method.
	 */
	public static class PeakExplained {
		public int level = 0; // explanation level
		public float experimentalMZ = 0;
		public float foundMZ = -1; // exact theoritical mz of proposed formula
		public float isotopicIntensityExplained = 0; // fraction of intensity of isotopic peaks that can be explained
		public String formula = ""; // molecular formula (explanation)
		public IsotopicPattern ip = null; // Isotopic pattern object
		// Explanation levels
		public static final int NOT_EXPLAINED = 0;
		public static final int ONLYMZ = 1;
		public static final int MZ_ISOTOPIC_NOT_PERFECT = 2;
		public static final int COMPLETE = 3;

		public String toString() {
			String result = "";
			result += level;
			if (foundMZ > 0) {
				result += "," + formula + "," + foundMZ + "," + Math.abs(experimentalMZ - foundMZ) + "#";
				// # is require in order to be replaced for additional info
				result += isotopicIntensityExplained * 100 + ",";
				if (ip != null) {
					for (float mz : ip.sortedMZ()) {
						result += "(," + mz + "," + ip.asMap().get(mz) + ",) ";

					}
				}
			}
			return result.trim();
		}

	}

	/**
	 * Different similarity measures between high-resolution mass spectra when using
	 * compare - THIS - is 1
	 */
	public static class SimilarityHR {
		public float recall;
		public float precision;
		public float f1;

		public float weightedRecall;
		public float weightedPrecision;
		public float weightedF1;

		public float dotProduct;
		public float scaledDotProduct;

		public float lowResolutionCompositeSimilarity;

		public int nPeaks1 = 0; // when using compare - THIS - is 1
		public int nPeaks2 = 0;
		public int nPeaksOnly1 = 0;
		public int nPeaksOnly2 = 0;
		public int nPeaksBoth = 0;

		public float fractionIntensityIn1AndIn2 = 0;
		public float fractionIntensityIn2AndIn1 = 0;

		public String toString() {
			String result = "recall," + recall + ",";
			result += "precision," + precision + ",";
			result += "f1," + f1 + ",";
			result += "weightedRecall," + weightedRecall + ",";
			result += "weightedPrecision," + weightedPrecision + ",";
			result += "weightedF1," + weightedF1 + ",";
			result += "dotProduct," + dotProduct + ",";
			result += "scaledDotProduct," + scaledDotProduct + ",";
			result += "lowResolutionCompositeSimilarity," + lowResolutionCompositeSimilarity + ",";
			result += "nPeaks1," + nPeaks1 + ",";
			result += "nPeaks2," + nPeaks2 + ",";
			result += "nPeaksOnly1," + nPeaksOnly1 + ",";
			result += "nPeaksOnly2," + nPeaksOnly2 + ",";
			result += "nPeaksBoth," + nPeaksBoth;
			return result.trim();
		}
	}

	/**
	 * print mass spectrum
	 */
	public void print() {
		if (mzs.length != intensities.length) {
			throw new RuntimeException();
		}
		for (int i = 0; i < mzs.length; i++) {
			System.out.println(mzs[i] + " " + intensities[i]);
		}
	}

	/**
	 * THIS - is 1 in result
	 * 
	 * @param other       second (2) mass spectrum
	 * @param mzThreshold threshold to consider peaks as the same
	 * @return similarity measures
	 */
	public SimilarityHR compare(MassSpectrumHR other, float mzThreshold) {

		int inBoth = 0;
		int inThisAll = 0;
		int inOtherAll = 0;
		int onlyInThis = 0;
		int onlyInOther = 0;

		float intensityInThisInBoth = 0;
		float intensityInOtherInBoth = 0;
		float intensityInThisAll = 0;
		float intensityInOtherAll = 0;
		float intensityOnlyInThis = 0;
		float intensityOnlyInOther = 0;

		float xx = 0;
		float yy = 0;
		float xy = 0;

		float xxW = 0;
		float yyW = 0;
		float xyW = 0;

		for (int i = 0; i < this.mzs.length; i++) {
			inThisAll++;
			intensityInThisAll += this.intensities[i];
			xx += this.intensities[i] * this.intensities[i];
			xxW += this.intensities[i] * this.mzs[i];

			boolean foundInOther = false;
			for (int j = 0; j < other.mzs.length; j++) {
				if (Math.abs(other.mzs[j] - this.mzs[i]) < mzThreshold) {
					foundInOther = true;
					xy += this.intensities[i] * other.intensities[j];
					xyW += Math.sqrt(this.intensities[i] * other.intensities[j] * this.mzs[i] * other.mzs[j]);
				}
			}
			if (foundInOther) {
				inBoth++;
				intensityInThisInBoth += this.intensities[i];
			} else {
				onlyInThis++;
				intensityOnlyInThis += this.intensities[i];
			}
		}

		for (int i = 0; i < other.mzs.length; i++) {
			inOtherAll++;
			intensityInOtherAll += other.intensities[i];
			yy += other.intensities[i] * other.intensities[i];
			yyW += other.intensities[i] * other.mzs[i];
			boolean foundInThis = false;
			for (int j = 0; j < this.mzs.length; j++) {
				if (Math.abs(this.mzs[j] - other.mzs[i]) < mzThreshold) {
					foundInThis = true;
				}
			}
			if (foundInThis) {
				intensityInOtherInBoth += other.intensities[i];
			} else {
				onlyInOther++;
				intensityOnlyInOther += other.intensities[i];
			}
		}

		SimilarityHR result = new SimilarityHR();
		result.precision = inBoth * 1.0f / inThisAll;
		result.recall = inBoth * 1.0f / inOtherAll;
		result.weightedPrecision = intensityInThisInBoth / intensityInThisAll;
		result.weightedRecall = intensityInOtherInBoth / intensityInOtherAll;
		result.f1 = 2 * (result.precision * result.recall) / (result.precision + result.recall);
		result.weightedF1 = 2 * (result.weightedPrecision * result.weightedRecall)
				/ (result.weightedPrecision + result.weightedRecall);

		result.nPeaks1 = inThisAll;
		result.nPeaks2 = inOtherAll;
		result.nPeaksOnly1 = onlyInThis;
		result.nPeaksOnly2 = onlyInOther;
		result.nPeaksBoth = inBoth;

		result.fractionIntensityIn1AndIn2 = intensityInThisInBoth;
		result.fractionIntensityIn2AndIn1 = intensityInOtherInBoth;

		if (Math.abs(intensityInOtherAll - intensityOnlyInOther - intensityInOtherInBoth) > 1e-2) {
			System.out.println(intensityInOtherAll + " " + intensityOnlyInOther + " " + intensityInOtherInBoth);
			throw new RuntimeException("unknown error");
		}
		if (Math.abs(intensityInThisAll - intensityOnlyInThis - intensityInThisInBoth) > 1e-2) {
			throw new RuntimeException("unknown error");
		}

		result.dotProduct = (xy * xy) / (xx * yy);
		result.scaledDotProduct = (xyW * xyW) / (xxW * yyW);

		result.lowResolutionCompositeSimilarity = (this.toLowResolution(1e-7f)).identity(other.toLowResolution(1e-7f));

		return result;
	}

	/**
	 * This method calls the explainPeaks(FragmentIon[] possibleIons,
	 * HashMap<String, String> properties) method on the given molecule. A more
	 * detailed explanation of the algorithm is given there.
	 * 
	 * 
	 * @param smiles     a SMILES string for the molecule
	 * @param properties Properties as a hash table. For example, mzThreshold can be
	 *                   retrieved as
	 *                   Float.parseFloat(properties.get("mzThreshold"))
	 * @param nCleavages The number of bond breaks that is considered when
	 *                   constructing a list of possible fragments. Can take values
	 *                   ​​1 2 3 4. In this case, a ring break requires at least two
	 *                   bond breaks. Migration of hydrogen atoms, migration of
	 *                   fluorine, and additional loss of methyls and halogens are
	 *                   not included in the number of bond breaks. The value 4
	 *                   (!!!!!) does not mean the actual 4 bond breaks, but a list
	 *                   of all possible combinations of atoms that do not
	 *                   contradict the original molecular formula.
	 * @return
	 */
	public PeakExplained[] explainPeaks(String smiles, HashMap<String, String> properties, int nCleavages) {
		MolGraph m = new MolGraph(smiles);
		int maxHDrift = Integer.parseInt(properties.get("maxHDrift"));
		int maxHLoss = Integer.parseInt(properties.get("maxHLoss"));
		int maxFMigration = Integer.parseInt(properties.get("maxFMigration"));
		ArrayList<FragmentIon> f = null;
		if (nCleavages > 0) {
			f = m.fragmentsFormulas(nCleavages, maxHDrift, maxHLoss, maxFMigration);
		} else {
			f = m.allPossibleFragmentsFormulas();
		}
		HashSet<String> molecularFormulas = new HashSet<String>();
		ArrayList<FragmentIon> f1 = new ArrayList<FragmentIon>();
		for (FragmentIon fi : f) {
			if (!molecularFormulas.contains(fi.formula)) {
				molecularFormulas.add(fi.formula);
				f1.add(fi);
			}
		}
		FragmentIon[] x = f1.toArray(new FragmentIon[f1.size()]);
		PeakExplained[] result = null;
		result = explainPeaks(x, properties);
		return result;
	}

	/**
	 * 
	 * The method tries to explain each peak of the mass spectrum according to the
	 * explanation levels (0-3). If there is a possible molecular formula for a peak
	 * (with an accuracy of mzThreshold) among the fragments listed in possibleIons,
	 * then this peak has an explanation level of 1 or higher. If such a formula is
	 * not found, then it gets a level of 0. Sometimes this can happen when several
	 * peaks overlap, and the centroid is significantly shifted from the real m/z of
	 * each of the peaks. After finding the molecular formula to which the peak
	 * corresponds, an isotopic distribution is considered for this peak, and for
	 * each peak of the isotopic distribution its theoretical intensity is
	 * calculated. If it is less than intensityThreshold, then we do not look for
	 * this isotopic peak. Otherwise, we look for an isotopic peak. For each
	 * isotopic peak with a theoretical intensity higher than intensityThreshold, we
	 * perform the following operations. If such a peak is in the spectrum (with an
	 * accuracy of mzThreshold), then we look at its intensity. If it differs from
	 * the theoretical one by no more than absoluteDifferenceForIsotopic *OR* by no
	 * more than percentDifferenceForIsotopic percent (the theoretical intensity is
	 * in the denominator when calculating the percentage error!) then we mark this
	 * peak as 3 (perfect), otherwise as 2 (found but not perfect). If the isotopic
	 * peak peak_iso is not found (with an accuracy of mzThreshold), then a search
	 * is performed for a more intense (intensity higher than the sought theoretical
	 * intensity of the isotopic peak) peak found_peak in the spectrum, such that
	 * Math.abs(mz_found_peak - mz_peak_iso) / mz_peak_iso < (1 / resolution), i.e.
	 * peak with which the original isotopic peak overlapped, then the isotopic peak
	 * is labeled as 2. As a result, the original peak whose interpretation was
	 * performed is labeled with explanation level 3 if ALL isotopic peaks (with
	 * theoretical intensity above intensityThreshold) are labeled as 3. The peak
	 * whose interpretation was performed is labeled with explanation level 2 if ALL
	 * isotopic peaks (with theoretical intensity above intensityThreshold) are
	 * labeled as 2 or 3. Otherwise, it is labeled as 1.
	 * 
	 * @param possibleIons Possible fragment ions constructed for this molecule
	 * @param properties   Properties as a hash table. For example, mzThreshold can
	 *                     be retrieved as
	 *                     Float.parseFloat(properties.get("mzThreshold"))
	 * @return Explanation for each of the peaks. The returned array has a length
	 *         equal to this.mzs.length
	 */

	public PeakExplained[] explainPeaks(FragmentIon[] possibleIons, HashMap<String, String> properties) {
		int nPeaks = this.mzs.length;
		IsotopicPattern[] patterns = new IsotopicPattern[possibleIons.length];
		for (int i = 0; i < patterns.length; i++) {
			patterns[i] = possibleIons[i].isotopicPatternSingleCation(thresholdGenerateIsotopic);
		}
		PeakExplained[] result = new PeakExplained[this.mzs.length];
		float mzThreshold = Float.parseFloat(properties.get("mzThreshold"));
		float resolution = Float.parseFloat(properties.get("resolution"));
		float percentDifferenceForIsotopic = Float.parseFloat(properties.get("percentDifferenceForIsotopic"));
		float absoluteDifferenceForIsotopic = Float.parseFloat(properties.get("absoluteDifferenceForIsotopic"));
		float intensityThreshold = Float.parseFloat(properties.get("intensityThreshold"));
		for (int i = 0; i < nPeaks; i++) {
			PeakExplained e = new PeakExplained();
			e.level = PeakExplained.NOT_EXPLAINED;
			for (int j = 0; j < patterns.length; j++) {
				for (float mz : patterns[j].asMap().keySet()) {
					if (Math.abs(mz - this.mzs[i]) < mzThreshold) {// Possible match
						int explanation = PeakExplained.ONLYMZ;
						boolean peakInIsotopicNotFound = false;
						boolean peakInIsotopicWrong = false;
						float isotopicIntensityExplained = patterns[j].asMap().get(mz);
						for (float mz2 : patterns[j].asMap().keySet()) {
							float intensityToFind = intensities[i] * patterns[j].asMap().get(mz2)
									/ patterns[j].asMap().get(mz);
							if (intensityToFind > intensityThreshold) {
								boolean found = false;
								boolean foundExactly = false;
								for (int k = 0; k < nPeaks; k++) {
									if (Math.abs(mz2 - this.mzs[k]) < mzThreshold) {
										found = true;
										if ((Math.abs(intensityToFind - this.intensities[k]) / intensityToFind)
												* 100 < percentDifferenceForIsotopic) {
											foundExactly = true;
										}
										if (Math.abs(intensityToFind
												- this.intensities[k]) < absoluteDifferenceForIsotopic) {
											foundExactly = true;
										}
									}
								}
								if (!foundExactly) {
									for (int k = 0; k < nPeaks; k++) {
										if (Math.abs(mz2 - this.mzs[k]) / mz2 < (1 / resolution)) {
											if (this.intensities[k] > intensityToFind) {
												found = true;
											}
										}
									}
								}
								if (!found) {
									peakInIsotopicNotFound = true;
								} else {
									if (!foundExactly) {
										peakInIsotopicWrong = true;
									} else {
										if (mz != mz2) {
											isotopicIntensityExplained += patterns[j].asMap().get(mz2);
										}
									}
								}
							}
						}
						if (!peakInIsotopicNotFound) {// all isotopic peaks present
							if (!peakInIsotopicWrong) {
								explanation = PeakExplained.COMPLETE;
							} else {
								explanation = PeakExplained.MZ_ISOTOPIC_NOT_PERFECT;
							}
						}
						if (explanation > e.level) {
							e.level = explanation;
							e.isotopicIntensityExplained = isotopicIntensityExplained;
							e.foundMZ = mz;
							e.ip = patterns[j];
							e.formula = possibleIons[j].formula;
							e.experimentalMZ = mzs[i];
						}
					}
				}
			}
			result[i] = e;
		}
		return result;
	}

	/**
	 * This method returns 0 if the mass spectrum does not contain peaks in the
	 * range from M-3.5 to M+2.5 and 1 otherwise, where M is the expected molecular
	 * weight. This range takes into account possible hydrogen losses. If it returns
	 * 0, it means that nothing resembling a molecular ion is observed at all.
	 * 
	 * @param smiles SMILES string
	 * @param fw     filewriter (obligatory open and non-zero!) where log will be
	 *               written
	 * @return 0 or 1
	 */
	public int findMolecularIonIntegerMass(String smiles, FileWriter fw) {
		try {
			MolGraph m = new MolGraph(smiles);
			HashMap<Integer, Integer> mf = m.molecularFormulaFragment(m.markAllAtoms());
			float maxMZTarget = (new FragmentIon(mf)).isotopicPatternSingleCation(thresholdGenerateIsotopic).mainMZ();
			boolean found = false;
			for (int i = 0; i < this.mzs.length; i++) {
				if (mzs[i] > maxMZTarget - 3.5) {
					if (mzs[i] < maxMZTarget + 2.5) {
						found = true;
						fw.write("*?Molecular?," + mzs[i] + "," + intensities[i] + "\n");
					}
				}
			}
			if (!found) {
				fw.write("**No molecular ion, even at low resolution!**\n");
			}
			return found ? 1 : 0;
		} catch (IOException e) {
			e.printStackTrace();
			throw new RuntimeException(e.getMessage());
		}
	}

	/**
	 * Fraction of total ion current above expected molecular ion. This value should
	 * be zero, otherwise it is either identification error or noise/impurity.
	 * 
	 * @param smiles SMILES
	 * @return fraction of ion current above molecular ion
	 */
	public float fractionIonsAboveMolecularIon(String smiles) {
		MolGraph m = new MolGraph(smiles);
		HashMap<Integer, Integer> mf = m.molecularFormulaFragment(m.markAllAtoms());
		IsotopicPattern ip = (new FragmentIon(mf)).isotopicPatternSingleCation(thresholdGenerateIsotopic);
		float maxmz = 0;
		for (float mz : ip.asMap().keySet()) {
			if (maxmz < mz) {
				maxmz = mz;
			}
		}
		float full = 0;
		float above = 0;
		for (int i = 0; i < this.mzs.length; i++) {
			full += intensities[i];
			if (mzs[i] > maxmz + 1.5) {
				above += intensities[i];
			}
		}
		return above / full;
	}

	/**
	 * Algorithm for searching for a molecular ion. The peaks M+, [M-H+], [M-2H+]
	 * .... [2-xH+] where x is maxHLostMI, are searched for. For each, if the EXACT
	 * mass is found (with an accuracy of mzThreshold), then the corresponding peak
	 * gets level 1. Unlike the explainPeaks method, where all the peaks of the
	 * spectrum (both the main isotopic and minor, for example 13C) were considered
	 * one by one, here at the first stage only the main, most intense peaks of the
	 * isotopic distribution are taken. Then all the isotopic peaks are searched for
	 * this peak. For each isotopic peak. If such a peak is in the spectrum (with an
	 * accuracy of mzThreshold), then we look at its intensity. If it differs from
	 * the theoretical one by no more than absoluteDifferenceForIsotopic *OR* by no
	 * more than percentDifferenceForIsotopic percent (the theoretical intensity is
	 * in the denominator when calculating the percentage error!) then we mark this
	 * peak as found3 (perfect), otherwise as found2 (found but not perfect). If the
	 * isotopic peak peak_iso is not found (with an accuracy of mzThreshold), then a
	 * search is performed for a peak found_peak in the spectrum, such that
	 * Math.abs(mz_found_peak - mz_peak_iso) / mz_peak_iso < (1 / resolution), i.e.
	 * peak with which the original isotopic peak overlapped, then the isotopic peak
	 * is labeled as found2. Next, for the entire ion, we calculate the fraction of
	 * the ion current of the isotopic distribution labeled as found2 and as found3
	 * (including the intensity of the main peak from which we started, it is always
	 * labeled found3). If the fraction of the intensity labeled as found3 exceeds
	 * fractionIsotopicThreshold, then the corresponding peak is labeled as
	 * explained at level 3. If this is not true, then if the sum of the fractions
	 * of the intensities labeled as found2 and found3 exceeds
	 * fractionIsotopicThreshold, then the corresponding peak is labeled as
	 * explained at level 2.
	 * 
	 * @param smiles     SMILES string
	 * @param properties Properties as a hash table. For example, mzThreshold can be
	 *                   retrieved as
	 *                   Float.parseFloat(properties.get("mzThreshold"))
	 * @param fw         filewriter (obligatory open and non-zero!) where log will
	 *                   be written
	 * @return the method returns an array of length maxHLostMI with explanation
	 *         levels (0, 1, 2, 3) for M, M-H, etc.
	 */
	public int[] findMolecularIon(String smiles, HashMap<String, String> properties, FileWriter fw) {
		int maxHLostMI = Integer.parseInt(properties.get("maxHLostMI"));
		float fractionIsotopicThreshold = Float.parseFloat(properties.get("fractionIsotopicThreshold"));
		float mzThreshold = Float.parseFloat(properties.get("mzThreshold"));
		float resolution = Float.parseFloat(properties.get("resolution"));
		float percentDifferenceForIsotopic = Float.parseFloat(properties.get("percentDifferenceForIsotopic"));
		float absoluteDifferenceForIsotopic = Float.parseFloat(properties.get("absoluteDifferenceForIsotopic"));
		MolGraph m = new MolGraph(smiles);
		HashMap<Integer, Integer> mf = m.molecularFormulaFragment(m.markAllAtoms());
		ArrayList<IsotopicPattern> mi = new ArrayList<IsotopicPattern>();
		if (mf.get(1) == null) {
			mf.put(1, 0);
		}
		for (int i = 0; i < maxHLostMI + 1; i++) {
			if ((mf.get(1) != null) && (mf.get(1) - i >= 0)) {
				@SuppressWarnings("unchecked")
				HashMap<Integer, Integer> mf1 = (HashMap<Integer, Integer>) mf.clone();
				mf1.put(1, mf.get(1) - i);
				mi.add((new FragmentIon(mf1)).isotopicPatternSingleCation(thresholdGenerateIsotopic));
			} else {
				mi.add(null);
			}
		}
		int count = 0;
		int[] result = new int[mi.size()];
		for (IsotopicPattern ion : mi) {
			ArrayList<String> outWrite = new ArrayList<String>();
			if (ion != null) {
				int found = 0;
				float maxMZTarget = ion.mainMZ();
				float maxIntens = -1;
				float maxMZFound = -1;

				for (int i = 0; i < this.mzs.length; i++) {
					float mz2 = this.mzs[i];
					if (Math.abs(mz2 - maxMZTarget) < mzThreshold) {
						found = 1;
						maxIntens = this.intensities[i];
						maxMZFound = this.mzs[i];
					}
				}
				if (maxMZFound > 0) {
					outWrite.add("*Molecular - " + count + " H; m/z theor.," + maxMZTarget + ",isotopic fraction,"
							+ ion.asMap().get(maxMZTarget) + ",m/z found," + maxMZFound + ",intensity," + maxIntens);
				} else {
					outWrite.add("*Molecular - " + count + " H; m/z theor.," + maxMZTarget + ",isotopic fraction,"
							+ ion.asMap().get(maxMZTarget) + ",NOT FOUND");
				}
				if (found > 0) {
					float intensityFoundAtLevel2 = ion.asMap().get(maxMZTarget);
					float intensityFoundAtLevel3 = ion.asMap().get(maxMZTarget);
					for (float mz : ion.asMap().keySet()) {
						boolean found3 = false;
						float foundMZ = -1;
						float foundIntens = -1;
						if (mz != maxMZTarget) {
							for (int i = 0; i < this.mzs.length; i++) {
								float mz2 = this.mzs[i];
								if (Math.abs(mz2 - mz) < mzThreshold) {
									float intensity = ion.asMap().get(mz) * maxIntens / ion.asMap().get(maxMZTarget);
									float intensity2 = this.intensities[i];
									if (((100 * Math.abs(intensity - intensity2)
											/ intensity) < percentDifferenceForIsotopic)
											|| (Math.abs(intensity - intensity2) < absoluteDifferenceForIsotopic)) {
										intensityFoundAtLevel3 += ion.asMap().get(mz);
										intensityFoundAtLevel2 += ion.asMap().get(mz);
										found3 = true;
										foundMZ = mz2;
										foundIntens = intensity2;
									}
								}
							}
							if (!found3) {
								for (int i = 0; i < this.mzs.length; i++) {
									float mz2 = this.mzs[i];
									if (Math.abs(mz - this.mzs[i]) / mz < (1 / resolution)) {
										intensityFoundAtLevel2 += ion.asMap().get(mz);
										foundMZ = mz2;
										foundIntens = this.intensities[i];
									}
								}
							}

							if (foundMZ > 0) {
								outWrite.add("*Molecular - " + count + " H; m/z theor.," + mz + ",isotopic fraction,"
										+ ion.asMap().get(mz) + ",m/z found," + foundMZ + ",intensity," + foundIntens
										+ ",level " + (found3 ? "3" : "2"));
							} else {
								outWrite.add("*Molecular - " + count + " H; m/z theor.," + mz + ",isotopic fraction,"
										+ ion.asMap().get(mz) + ",NOT FOUND");
							}

						}
					}
					if (intensityFoundAtLevel2 > fractionIsotopicThreshold) {
						found = 2;
					}
					if (intensityFoundAtLevel3 > fractionIsotopicThreshold) {
						found = 3;
					}
					outWrite.set(0, outWrite.get(0) + "," + "main isotopic peak," + "level " + found
							+ ",isotopic intensity perfectly found," + intensityFoundAtLevel3);
				}
				result[count] = found;
			} else {
				result[count] = 0;
			}
			count++;
			for (String x : outWrite) {
				try {
					if (fw != null) {
						fw.write(x + "\n");
					}
				} catch (IOException e) {
					e.printStackTrace();
					throw new RuntimeException(e.getMessage());
				}
			}
		}
		return result;
	}

	/**
	 * 
	 * @return array of mz values
	 */
	public float[] getMzs() {
		return mzs;
	}

	/**
	 * 
	 * @param mzs array of mz values
	 */
	public void setMzs(float[] mzs) {
		this.mzs = mzs;
	}

	/**
	 * 
	 * @return intensity values (the same length as mzs)
	 */
	public float[] getIntensities() {
		return intensities;
	}

	/**
	 * 
	 * @param intensities intensity values (the same length as mzs)
	 */
	public void setIntensities(float[] intensities) {
		this.intensities = intensities;
	}
}
