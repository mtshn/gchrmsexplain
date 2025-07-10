package ru.ac.phyche.gchrmsexplain;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import org.apache.commons.lang3.tuple.Pair;

import ru.ac.phyche.gchrmsexplain.MassSpectrumHR.PeakExplained;
import ru.ac.phyche.gchrmsexplain.NISTDataBase.SearchResult;

/**
 * The class interprets the mass spectrum using the explainPeaks method
 * belonging to class MassSpectrumHR.
 */
public class ExplainMassSpectrumHR {

	/**
	 * The method interprets the mass spectrum using the explainPeaks method
	 * belonging to class MassSpectrumHR.
	 * 
	 * @param spectrumNameSmiles The first member of the tuple is the spectrum. The
	 *                           second array of strings of length 2 contains the
	 *                           molecule name and the SMILES string (structure), in
	 *                           that order.
	 * @param name2              Additional name of comound (usually name of file
	 *                           with spectrum)
	 * @param outputFile         filewriter (obligatory open and non-zero!) where
	 *                           log will be written
	 * @param properties         Properties as a hash table. For example,
	 *                           mzThreshold can be retrieved as
	 *                           Float.parseFloat(properties.get("mzThreshold"))
	 * @return The method returns an array of length 9, containing the results of
	 *         interpretation for 2, 3 bond breaks and for all possible molecular
	 *         formulas. For each case, the values ​​of the fraction of the total
	 *         intrinsic current interpreted at level 3, 2, 1 are given (the figure
	 *         for level 1 includes the fraction of peaks interpreted at level 1 and
	 *         ABOVE, also for level 2). The order is: 2 bonds level 3, 2 bonds
	 *         level 2, 2 bonds 2 level 1, 3 bonds ..3 2 1... , all possible
	 *         fragments (3 2 1 level).
	 * @throws IOException io
	 */
	public static float[] explainSpectrum(Pair<MassSpectrumHR, String[]> spectrumNameSmiles, String name2,
			FileWriter outputFile, HashMap<String, String> properties) throws IOException {
		String smiles = spectrumNameSmiles.getRight()[1];
		String name = spectrumNameSmiles.getRight()[0];
		MassSpectrumHR x = spectrumNameSmiles.getLeft();
		// four levels of fragmentation: from single cleavage to all fragments
		PeakExplained[] r3Cleavage = x.explainPeaks(smiles, properties, 3);
		PeakExplained[] r2Cleavage = x.explainPeaks(smiles, properties, 2);
		PeakExplained[] r1Cleavage = x.explainPeaks(smiles, properties, 1);
		PeakExplained[] rAllFragments = x.explainPeaks(smiles, properties, -1);

		float allIntensity = 0;
		float[] explainedIntensity1 = new float[4]; // four levels of fragmentation: from single cleavage to all
													// fragments
		float[] explainedIntensity2 = new float[4];
		float[] explainedIntensity3 = new float[4];
		int[] fragmentingLevel = new int[x.getMzs().length];
		int j = 0;
		boolean[] explained = new boolean[x.getMzs().length];
		String[] explanations = new String[x.getMzs().length];
		for (int i = 0; i < x.getMzs().length; i++) {
			allIntensity += x.getIntensities()[i];
		}
		for (PeakExplained[] r : new PeakExplained[][] { r1Cleavage, r2Cleavage, r3Cleavage, rAllFragments }) {
			j++;
			if (r.length != x.getMzs().length) {
				throw new RuntimeException("Error array length");
			}
			for (int i = 0; i < r.length; i++) {
				if (!explained[i]) {
					explanations[i] = r[i].toString();
				}
				if (r[i].level >= 1) {
					explainedIntensity1[j - 1] += x.getIntensities()[i];
					if (!explained[i]) {
						explained[i] = true;
						fragmentingLevel[i] = j;
					}
				}
				if (r[i].level >= 2) {
					explainedIntensity2[j - 1] += x.getIntensities()[i];
				}
				if (r[i].level >= 3) {
					explainedIntensity3[j - 1] += x.getIntensities()[i];
				}
			}

		}
		outputFile.write(smiles + "\n\"" + name + "\"\n\"" + name2 + "\"\n" + "mz_experimen"
				+ "t,intensity,explanation_level,formula,mz_theory,delta_mz,"
				+ "fragmentation_level,fraction_intensity,fraction_isotopic_explained,isotopic_distribution_theory\n");
		for (int i = 0; i < x.getMzs().length; i++) {
			String fragmentationLevel = fragmentingLevel[i] <= 3 ? (fragmentingLevel[i] + "") : "all_fragments";
			String fractionIntensity = (100 * x.getIntensities()[i] / allIntensity) + "";
			String explanation = explanations[i].replace("#", "," + fragmentationLevel + "," + fractionIntensity + ",");
			outputFile.write(x.getMzs()[i] + "," + x.getIntensities()[i] + "," + explanation + "\n");
		}

		float[] result = new float[9];
		for (int i = 0; i < 4; i++) {
			float e3 = 100 * explainedIntensity3[i] / allIntensity;
			float e2 = 100 * explainedIntensity2[i] / allIntensity;
			float e1 = 100 * explainedIntensity1[i] / allIntensity;
			String header = "Fraction of explained intensity when up to " + (i + 1) + " bonds are broken";
			if (i == 3) {
				header = "Fraction of explained intensity with all possible fragments";
			}
			outputFile.write(header + "\n");
			outputFile.write("Explained at level 3:," + e3);
			outputFile.write(",Explained at level 2:," + e2);
			outputFile.write(",Explained at level 1:," + e1 + "\n");
			if (i > 0) {
				result[(i - 1) * 3] = e3;
				result[(i - 1) * 3 + 1] = e2;
				result[(i - 1) * 3 + 2] = e1;
			}
		}
		return result;
	}

	private static float[] explainCSVFileSpectrum(String filename, String smiles, FileWriter outputFile,
			HashMap<String, String> properties) throws IOException {
		MassSpectrumHR x = MassSpectrumHR.fromThermoCSVFile(filename,
				Float.parseFloat(properties.get("csvLoadIntensityThreshold")),
				Float.parseFloat(properties.get("thresholdGenerateIsotopic")));
		String name = "";
		Pair<MassSpectrumHR, String[]> spectrumNameSmiles = Pair.of(x, new String[] { smiles, name });
		return explainSpectrum(spectrumNameSmiles, filename, outputFile, properties);
	}

	/**
	 * Explain mass spectrum and search in NIST database
	 * 
	 * @param csvFilePath     mass spectrum in "thermo" format
	 * @param compoundName    name
	 * @param smiles          structure
	 * @param properties      properties hash table (see other similar methods)
	 * @param shortTable      filewriter (obligatory open and non-zero!) where log
	 *                        will be written
	 * @param fullExplanation filewriter (obligatory open and non-zero!) where log
	 *                        will be written
	 * @param fwSearchResult  filewriter (obligatory open and non-zero!) where log
	 *                        will be written (NIST search)
	 * @param noNistSearch    disable nist search
	 * @param nist            NIST MS database in own format, see NISTDataBase class
	 */
	private static void explainAndSearch(String csvFilePath, String compoundName, String smiles,
			HashMap<String, String> properties, FileWriter shortTable, FileWriter fullExplanation,
			FileWriter fwSearchResult, boolean noNistSearch, NISTDataBase nist) {
		String name = compoundName;
		String shortFileName = (new File(csvFilePath)).getName();
		System.out.println(shortFileName);
		try {
			float csvLoadTreshold = Float.parseFloat(properties.get("csvLoadIntensityThreshold"));
			boolean nistSearch = (!noNistSearch) && (fwSearchResult != null) && (nist != null);
			float foundMF2 = 0;
			int foundRank = -1;
			if (nistSearch) {
				MassSpectrumLR x = MassSpectrumHR
						.fromThermoCSVFile(csvFilePath, 1E-6f,
								Float.parseFloat(properties.get("thresholdGenerateIsotopic")))
						.toLowResolution(csvLoadTreshold);
				SearchResult[] sr = nist.search(x, Integer.parseInt(properties.get("maxSearchRank")));
				String inchikey = ParsingNIST23.smilesToInchiKey(smiles);
				fwSearchResult.write("\"" + name + "\"," + smiles + ",\"" + shortFileName + "\"\n");
				for (int k = 0; k < sr.length; k++) {
					if (inchikey.equals(sr[k].result.getIds().inchiKeyNist)) {
						foundRank = k;
					}
					fwSearchResult.write(k + "," + sr[k].similarity + "," + sr[k].result.toStringFull() + "\n");
				}
				fwSearchResult.write("\n\n");
				foundMF2 = nist.compoundInDatabaseMatchFactor(smiles, x);
			}
			MassSpectrumHR y = MassSpectrumHR.fromThermoCSVFile(csvFilePath, csvLoadTreshold,
					Float.parseFloat(properties.get("thresholdGenerateIsotopic")));
			float[] result = explainCSVFileSpectrum(csvFilePath, smiles, fullExplanation, properties);
			int[] molecularIon = y.findMolecularIon(smiles, properties, fullExplanation);
			int molecularIonLR = y.findMolecularIonIntegerMass(smiles, fullExplanation);
			float fractionAboveMolecular = y.fractionIonsAboveMolecularIon(smiles);
			int molecularIonBestLevel = Math.max(molecularIon[0], Math.max(molecularIon[1], molecularIon[2]));
			String explanationRates = "";
			for (int k = 0; k < result.length; k++) {
				explanationRates += result[k] + ",";
			}
			String molecularIonLevels = "";
			for (int k = 0; k < molecularIon.length; k++) {
				molecularIonLevels += molecularIon[k] + ",";
			}
			String nistStringResult = nistSearch ? (foundRank + 1) + "," + foundMF2 : "";
			shortTable.write("\"" + name + "\",\"" + shortFileName + "\"," + smiles + "," + nistStringResult + ","
					+ molecularIonBestLevel + "," + molecularIonLR + "," + fractionAboveMolecular + ","
					+ explanationRates + "," + molecularIonLevels + "\n");
		} catch (Throwable e) {
			try {
				shortTable.write("\"" + name + "\",\"" + shortFileName + "\"," + smiles
						+ ",Cannot process spectrum! Exception\n");
			} catch (Exception e1) {
				throw new RuntimeException(e1.getMessage());
			}
			e.printStackTrace();
		}
	}

	/**
	 * For all practical purposes, it is recommended to use the App class.
	 * 
	 * @param args args
	 * @throws Exception e
	 */
	public static void main(String[] args) throws Exception {
		BufferedReader br = new BufferedReader(new FileReader("banochki_highres.txt"));
		FileWriter fw = new FileWriter("banochki_highres_explained.txt");
		FileWriter fw1 = new FileWriter("banochki_highres_explaination_level.txt");
		FileWriter fwSearchResult = new FileWriter("banochki_highres_nist_search.txt");
		FileWriter fwComparison = new FileWriter("banochki_highres_comparison.txt");

		String s = br.readLine();
		NISTDataBase nist = new NISTDataBase("/home/xxx/NIST23/mainlib_gcms_nist23.ms");
		HashMap<String, String> properties = App.loadProperties(args);
		while (s != null) {
			if (!s.trim().equals("")) {
				String[] splt = s.trim().split("\\s+");
				String name = splt[0];
				String smiles = splt[1];
				File[] xx = new File("./Data_new/xxx").listFiles();
				boolean found = false;
				ArrayList<String> pathsSpectra = new ArrayList<String>();
				for (int i = 0; i < xx.length; i++) {
					String[] f = xx[i].getName().split("\\_");
					if (f[0].equals(name) || f[1].equals(name) || f[f.length - 1].equals(name)) {
						explainAndSearch(xx[i].getAbsolutePath(), name, smiles, properties, fw1, fw, fwSearchResult,
								false, nist);
						found = true;
						pathsSpectra.add(xx[i].getAbsolutePath());
					}
				}
				float csvLoadThreshold = Float.parseFloat(properties.get("csvLoadIntensityThreshold"));
				float thresholdIsotopic = Float.parseFloat(properties.get("thresholdGenerateIsotopic"));
				float mzThreshold = Float.parseFloat(properties.get("mzThreshold"));
				fwComparison.write("\"" + name + "\"\n");
				for (String path1 : pathsSpectra) {
					for (String path2 : pathsSpectra) {
						try {
							MassSpectrumHR ms1 = MassSpectrumHR.fromThermoCSVFile(path1, csvLoadThreshold,
									thresholdIsotopic);
							MassSpectrumHR ms2 = MassSpectrumHR.fromThermoCSVFile(path2, csvLoadThreshold,
									thresholdIsotopic);
							String name1 = (new File(path1)).getName();
							String name2 = (new File(path1)).getName();
							fwComparison.write("\"" + name1 + "\",\"" + name2 + "\","
									+ ms1.compare(ms2, mzThreshold).toString() + "\n");
						} catch (Throwable e) {
							e.printStackTrace();
						}
					}
				}
				fwComparison.write("\n\n");
				if (!found) {
					fw1.write("\"" + name + "\",NOT FOUND\n");
				}
			} else {
				fw1.write("\n");
			}
			fw.flush();
			fw1.flush();
			fwSearchResult.flush();
			s = br.readLine();
		}
		br.close();
		fw.close();
		fw1.close();
		fwComparison.close();
		fwSearchResult.close();
	}

}
