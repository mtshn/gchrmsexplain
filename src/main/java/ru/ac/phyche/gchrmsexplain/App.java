package ru.ac.phyche.gchrmsexplain;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.File;
import java.io.IOException;
import java.io.ByteArrayInputStream;
import java.io.InputStream;
import java.net.URLDecoder;
import java.util.ArrayList;
import java.util.HashMap;

import org.apache.commons.lang3.tuple.Pair;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

/**
 * Сommand line interface program for interpretation of high resolution gc-ms mass spectra.
 * Dmitriy D. Matyushin, Anastasia Yu. Sholokhova, 2025
 * See for details comments in MassSpectrumHR and ExplainMassSpectrumHR classes.
 * Run with "--help" key in order to see options.
 */
public class App {

	private static String removeSharpComment(String s) {
		String result = "";
		boolean quotes = false;
		boolean go = true;
		for (char c : s.toCharArray()) {
			if ((!quotes) && (c == '#')) {
				go = false;
			}
			if (c == '"') {
				quotes = !quotes;
			}
			if (go) {
				result += c;
			}
		}
		return result;
	}

	public static HashMap<String, String> loadProperties(HashMap<String, String> defaultProp, String filename)
			throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(filename));
		String s = br.readLine();
		HashMap<String, String> result = new HashMap<String, String>();
		if (defaultProp != null) {
			result = defaultProp;
		}
		while (s != null) {
			s = removeSharpComment(s.trim()).trim();
			if (!s.equals("")) {
				String key = "";
				String property = "";
				boolean goKey = true;
				for (char c : s.toCharArray()) {
					if (!goKey) {
						property = property + c;
					}
					if (c == '=') {
						goKey = false;
					}
					if (goKey) {
						key = key + c;
					}
				}
				key = key.trim();
				property = property.trim();
				if (!key.equals("")) {
					result.put(key, property);
				}
			}
			s = br.readLine();
		}
		br.close();
		return result;
	}

	public static String explainOneSpectrum(HashMap<String, String> properties, String fileNameForLog,
			FileWriter fileWriterResult, Pair<MassSpectrumHR, String[]> spectrumNameSmiles) throws IOException {
		String smiles = spectrumNameSmiles.getRight()[1];
		String name = spectrumNameSmiles.getRight()[0];
		MassSpectrumHR y = spectrumNameSmiles.getLeft();
		float[] result0 = ExplainMassSpectrumHR.explainSpectrum(spectrumNameSmiles, fileNameForLog, fileWriterResult,
				properties);
		int[] molecularIon = y.findMolecularIon(smiles, properties, fileWriterResult);
		int molecularIonLR = y.findMolecularIonIntegerMass(smiles, fileWriterResult);
		float fractionAboveMolecular = y.fractionIonsAboveMolecularIon(smiles);
		int molecularIonBestLevel = Math.max(molecularIon[0], Math.max(molecularIon[1], molecularIon[2]));
		String explanationRates = "";
		for (int k = 0; k < result0.length; k++) {
			explanationRates += result0[k] + ",";
		}
		String molecularIonLevels = "";
		for (int k = 0; k < molecularIon.length; k++) {
			molecularIonLevels += molecularIon[k] + ",";
		}
		String result = ("\"" + name + "\",\"" + fileNameForLog + "\"," + smiles + "," + molecularIonBestLevel + ","
				+ molecularIonLR + "," + fractionAboveMolecular + "," + explanationRates + molecularIonLevels + "\n");
		return result;
	}

	private static IAtomContainer smilesToAtomContainer(String s) {
		try {
			SmilesParser parser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
			IAtomContainer mol = parser.parseSmiles(s.trim());
			Aromaticity arom = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
			arom.apply(mol);
			CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(mol.getBuilder());
			adder.addImplicitHydrogens(mol);
			return mol;
		} catch (Exception e) {
			throw new RuntimeException(e.getMessage());
		}
	}

	private static String atomContainerlToSmiles(IAtomContainer mol, boolean stereochemistry) throws CDKException {
		Aromaticity arom = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
		AtomContainerManipulator.suppressHydrogens(mol);
		arom.apply(mol);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(mol.getBuilder());
		adder.addImplicitHydrogens(mol);
		SmilesGenerator sg = null;
		if (!stereochemistry) {
			sg = new SmilesGenerator((SmiFlavor.Unique | SmiFlavor.UseAromaticSymbols) ^ SmiFlavor.AtomicMass);
		} else {
			sg = new SmilesGenerator((SmiFlavor.Absolute | SmiFlavor.UseAromaticSymbols) ^ SmiFlavor.AtomicMass);
		}
		String smiles = sg.create(mol);
		return smiles;
	}

	public static String[] splitCSV(String s, char separator) {
		ArrayList<String> r = new ArrayList<String>();
		int i = 0;
		boolean inQuotes = false;
		String x = "";
		while (i < s.length()) {
			if (!inQuotes) {
				if (s.charAt(i) == '"') {
					inQuotes = true;
				} else {
					if (s.charAt(i) == separator) {
						r.add(x);
						x = "";
					} else {
						x += s.charAt(i);
					}
				}
			} else {
				if (s.charAt(i) == '"') {
					if ((i < s.length() - 1) && (s.charAt(i + 1) == '"')) {
						x += '"';
						i++;
					} else {
						inQuotes = false;
					}
				} else {
					x += s.charAt(i);
				}
			}
			i++;
		}
		r.add(x);
		return r.toArray(new String[r.size()]);
	}

	public static Pair<MassSpectrumHR, String[]>[] loadSpectraFromManyFiles(String csvFileWithFilenames, String prefix,
			String fileFormat, String csvSpectrumHeader, int mzColumn, int intensColumn, float thresholdBy999,
			float thresholdGenerateIsotopic) throws IOException, CDKException {
		ArrayList<String> al = new ArrayList<String>();
		BufferedReader br = new BufferedReader(new FileReader(csvFileWithFilenames));
		String s = br.readLine();
		while (s != null) {
			if (!s.trim().equals("")) {
				if (splitCSV(s, ',').length >= 2) {
					al.add(s.trim());
				}
			}
			s = br.readLine();
		}
		br.close();
		@SuppressWarnings("unchecked")
		Pair<MassSpectrumHR, String[]>[] result = new Pair[al.size()];
		int i = 0;
		for (String s1 : al) {
			String filename = "";
			String name = "";
			String smiles = "";
			String[] splt = splitCSV(s1, ',');
			filename = new File(prefix, splt[0].trim()).getAbsolutePath();
			smiles = splt[1].trim();
			if (splt.length >= 3) {
				name = splt[2].trim();
			}
			MassSpectrumHR sp = loadSingleSpectrumFromFile(filename, fileFormat, csvSpectrumHeader, mzColumn,
					intensColumn, thresholdBy999, thresholdGenerateIsotopic);
			result[i] = Pair.of(sp, new String[] { name, smiles, filename });
			i++;
		}
		return result;
	}

	@SuppressWarnings("unchecked")
	public static Pair<MassSpectrumHR, String[]>[] loadSpectrumFromSingleFile(String filename, String name,
			String smiles, String fileFormat, String csvSpectrumHeader, int mzColumn, int intensColumn,
			float thresholdBy999, float thresholdGenerateIsotopic) throws IOException, CDKException {
		MassSpectrumHR x = loadSingleSpectrumFromFile(filename, fileFormat, csvSpectrumHeader, mzColumn, intensColumn,
				thresholdBy999, thresholdGenerateIsotopic);
		Pair<MassSpectrumHR, String[]>[] result = null;
		result = new Pair[] { Pair.of(x, new String[] { name, smiles }) };
		return result;
	}

	public static MassSpectrumHR loadSingleSpectrumFromFile(String filename, String fileFormat,
			String csvSpectrumHeader, int mzColumn, int intensColumn, float thresholdBy999,
			float thresholdGenerateIsotopic) throws IOException, CDKException {

		if (fileFormat.equals("MSP")) {
			return MassSpectrumHR.fromMSPFile(filename, thresholdBy999, thresholdGenerateIsotopic).getLeft();
		}

		if (fileFormat.equals("CSV")) {
			return MassSpectrumHR.fromAnyCSVFile(filename, thresholdBy999, thresholdGenerateIsotopic, mzColumn,
					intensColumn, csvSpectrumHeader, ',');
		}

		if (fileFormat.equals("TXT")) {
			return MassSpectrumHR.fromSimpleMZTableFile(filename, thresholdBy999, thresholdGenerateIsotopic);
		}
		throw new RuntimeException("Unknown File Type " + fileFormat);
	}

	public static Pair<MassSpectrumHR, String[]>[] loadMSPWithSMILESorINCHI(String filename, float thresholdBy999,
			float thresholdGenerateIsotopic) throws IOException, CDKException {
		BufferedReader br = new BufferedReader(new FileReader(filename));
		String s = br.readLine();
		ArrayList<Pair<MassSpectrumHR, String[]>> rslt = new ArrayList<Pair<MassSpectrumHR, String[]>>();
		while (s != null) {
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
			rslt.add(sp);
		}
		br.close();
		@SuppressWarnings("unchecked")
		Pair<MassSpectrumHR, String[]>[] result = new Pair[rslt.size()];
		result = rslt.toArray(result.clone());
		return result;
	}

	private static String c1(String s) {
		return s.toLowerCase().replaceAll("\\s+", " ").trim();
	}

	public static Pair<MassSpectrumHR, String[]>[] loadSDFNIST(String filename, float thresholdBy999,
			float thresholdGenerateIsotopic) throws IOException, CDKException {
		BufferedReader br = new BufferedReader(new FileReader(filename));
		String s = br.readLine();
		ArrayList<Pair<MassSpectrumHR, String[]>> rslt = new ArrayList<Pair<MassSpectrumHR, String[]>>();
		while (s != null) {
			// reading sdf/mol2 part
			String sdf = s + "\n";
			s = br.readLine();
			sdf += s + "\n";
			s = br.readLine();
			sdf += s + "\n";
			s = br.readLine();
			sdf += s + "\n";
			s = br.readLine();

			while ((s != null) && (!s.matches(".*>\\s*<.*"))) {
				sdf += s + "\n";
				s = br.readLine();
			}
			if (s != null) {
				InputStream sdfStream = new ByteArrayInputStream(sdf.getBytes());
				MDLV2000Reader reader = new MDLV2000Reader(sdfStream);
				IAtomContainer molecule = (IAtomContainer) reader.read(smilesToAtomContainer(""));
				reader.close();
				sdfStream.close();
				String smiles = "Structure_processing_error";
				try {
					smiles = atomContainerlToSmiles(molecule, true);
				} catch (Throwable e) {

				}
				MassSpectrumHR spectrum = null;
				String name = "";
				String mz = "";
				while ((!c1(s).contains("mass spectral peaks")) && (!c1(s).contains("mass spectrum"))) {
					if (s.toLowerCase().contains("<name>")) {
						s = br.readLine();
						name = s.trim();
					}
					s = br.readLine();
				}
				String s1 = c1(s);
				s = br.readLine();
				while ((!s.contains("$$$$")) && (!s.matches(".*>\\\\s*<.*"))) {
					mz += s + "\n";
					s = br.readLine();
				}
				while ((!s.contains("$$$$"))) {
					s = br.readLine();
				}
				s = br.readLine();
				if (s1.contains("mass spectrum")) {// sdf with msp block
					spectrum = MassSpectrumHR.fromMSPBlock(mz, thresholdBy999, thresholdGenerateIsotopic).getLeft();
				} else {
					spectrum = MassSpectrumHR.fromSimpleMZTable(mz, thresholdBy999, thresholdGenerateIsotopic);
				}
				rslt.add(Pair.of(spectrum, new String[] { name, smiles }));
			}
		}
		br.close();
		@SuppressWarnings("unchecked")
		Pair<MassSpectrumHR, String[]>[] result = new Pair[rslt.size()];
		result = rslt.toArray(result.clone());
		return result;
	}

	public static HashMap<String, String> loadProperties(String args[]) throws IOException {
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("--help") || args[i].equals("-?")
					|| args[i].equals("/?")) {
				printHelp();
				System.exit(0);
			}
		}
		if (args.length % 2 != 0) {
			System.out.println("Incorrect input");
			printHelp();
			System.exit(1);
		} else {
			for (int i = 0; i < args.length; i++) {
				if (((args[i].charAt(0) == '-') && (args[i].charAt(1) == '-')) || (args[i].equals("-O"))
						|| (args[i].equals("-O1"))) {
					if (i % 2 != 0) {
						System.out.println("Incorrect input");
						printHelp();
						System.exit(1);
					}
				} else {
					if (i % 2 == 0) {
						System.out.println("Incorrect input");
						printHelp();
						System.exit(1);
					}
				}
			}
		}
		String defaultPropFile = "./properties.txt";
		File jarPath = (new File(URLDecoder.decode(App.class.getProtectionDomain().getCodeSource().getLocation().getPath(),"UTF-8"))).getParentFile();
		defaultPropFile = (new File(jarPath,defaultPropFile)).getAbsolutePath();
		HashMap<String, String> properties = loadProperties(null, defaultPropFile);
		String propFile = null;
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("--propFile")) {
				if (i != args.length - 1) {
					propFile = args[i + 1];
				}
			}
		}

		if (propFile != null) {
			properties = loadProperties(properties, propFile);
		}

		for (int i = 0; i < args.length; i++) {
			if (!args[i].equals("--propFile")) {
				if (((args[i].charAt(0) == '-') && (args[i].charAt(1) == '-')) || (args[i].equals("-O"))
						|| (args[i].equals("-O1"))) {
					if (i != args.length - 1) {
						String key = args[i].replace("-", " ").trim();
						String property = args[i + 1].trim();
						if (!key.equals("")) {
							properties.put(key, property);
						}
					}
				}
			}
		}
		return properties;
	}

	public static void printHelpShort() {
		System.out.println("****************************************************");
		System.out.println("use with option --help in order to view using options");
		System.out.println("usage: --option1 <argument1> --option1 <argument1> ...");
		System.out.println("****************************************************");
	}

	public static void printHelp() {
		System.out.println("usage: --option1 <argument1> --option1 <argument1> ...\n"+"Options:\n"
				+ "--propFile   file with properties. All other options can be also set there.\n"
				+ "             by default properties.txt file is used.\n"
				+ "\n"
				+ "--SMILES     SMILES string (structures can also be provided in input file)\n"
				+ "--name       compound name (structures can also be provided in input file)\n"
				+ "\n"
				+ "-O           output file name, where full interpretation of mass spectra\n"
				+ "             will be written\n"
				+ "-O1          shortened output file: one molecule - one line.\n"
				+ "--fileFormat format of input spectra\n"
				+ "--inputFile  input file name\n"
				+ "--prefix     prefix (path) where the spectrum files listed in the input file are\n"
				+ "             located. Only for those file formats where one spectrum is one file.\n"
				+ "\n"
				+ "# Explanation of mass spectral peaks (high resolution)\n"
				+ "--mzThreshold                      accuracy of mass determination, Da\n"
				+ "--resolution                       HRMS resolution\n"
				+ "--percentDifferenceForIsotopic     the relative (percent) error in the intensity\n"
				+ "                                   of an isotopic peak at which it is considered\n"
				+ "                                   \"perfect\"\n"
				+ "--absoluteDifferenceForIsotopic    the absolute (base peak = 999) error in the\n"
				+ "                                   intensity of an isotopic peak at which it is\n"
				+ "                                   considered \"perfect\"\n"
				+ "If observed intensity = X and theoretical = Y, the isotopic peak is \"perfect\"\n"
				+ "if  100*(X-Y)/Y  <  percentDifferenceForIsotopic *OR* if\n"
				+ "(X-Y) <  absoluteDifferenceForIsotopic\n"
				+ "--intensityThreshold = 5            the intensity of an isotopic peak below\n"
				+ "                                    which we do not attempt to search for an\n"
				+ "                                    isotopic peak (base peak = 999)\n"
				+ "\n"
				+ "# Atom migrations\n"
				+ "--maxHDrift                         maximum number of H-atoms that can migrate\n"
				+ "                                    TO ion\n"
				+ "--maxHLoss                          maximum additional loss of H-atoms\n"
				+ "--maxFMigration                     max number of fluorine atoms that can migrate\n"
				+ "\n"
				+ "# Explaining molecular ion peak\n"
				+ "--maxHLostMI                        consider M, M-H ... M-xH when considering\n"
				+ "                                    molecular ion (usually 2)\n"
				+ "--fractionIsotopicThreshold         The fraction of isotopic peak intensities\n"
				+ "                                    that must be well explained when explaining\n"
				+ "                                    the molecular ion peak\n"
				+ "\n"
				+ "# General mass spectrometry settings\n"
				+ "--csvLoadIntensityThreshold         threshold, peaks with intensity below which\n"
				+ "                                    are discarded during the initial loading of\n"
				+ "                                    the spectrum (base peak (peak with maximum\n"
				+ "                                    intensity) = 999)\n"
				+ "--thresholdGenerateIsotopic         threshold, peaks with intensity below which\n"
				+ "                                    are discarded when calculating the isotopic\n"
				+ "                                    distribution. The intensity is calculated\n"
				+ "                                    from the ENTIRE intensity of the isotopic\n"
				+ "                                    distribution (not from the base peak).\n"
				+ "\n"
				+ "# CSV file load options\n"
				+ "--csvSpectrumHeader                 CSV table header, below which are m/z and\n"
				+ "                                    intensity values. All lines above the header\n"
				+ "                                    are ignored.\n"
				+ "--mzColumn                          The column number in the CSV table that\n"
				+ "                                    contains the m/z values. The FIRST column is\n"
				+ "                                    number 0!\n"
				+ "--intensColumn                      The column number in the CSV table that\n"
				+ "                                    contains the intensity values. The FIRST\n"
				+ "                                    column is number 0!");
	}

	public static void run(HashMap<String, String> properties) throws IOException, CDKException {
		String fileFormat = properties.get("fileFormat");
		String outputFile = properties.get("O");
		String outputFile1 = properties.get("O1");
		String smiles = properties.get("SMILES");
		String name = properties.get("name");

		Pair<MassSpectrumHR, String[]>[] data = null;
		float csvLoadThreshold = Float.parseFloat(properties.get("csvLoadIntensityThreshold"));
		float thresholdIsotopic = Float.parseFloat(properties.get("thresholdGenerateIsotopic"));
		for (String key : properties.keySet()) {
			System.out.println(key + " = " + properties.get(key));
		}
		if ((smiles != null) && (!smiles.trim().equals(""))) {
			if (fileFormat.equals("SDF") || fileFormat.equals("MSP1")) {
				System.out.println(
						"File formats SDF and MSP1 cannot be used with non-empty SMILES string set via properties file or command line");
				System.exit(1);
			}
			int mzColumn = Integer.parseInt(properties.get("mzColumn"));
			int intensColumn = Integer.parseInt(properties.get("intensColumn"));
			data = loadSpectrumFromSingleFile(properties.get("inputFile"), name, smiles, fileFormat,
					properties.get("csvSpectrumHeader"), mzColumn, intensColumn, csvLoadThreshold, thresholdIsotopic);
		} else {
			if (fileFormat.equals("SDF")) {
				data = loadSDFNIST(properties.get("inputFile"), csvLoadThreshold, thresholdIsotopic);
			} else {
				if (fileFormat.equals("MSP1")) {
					data = loadMSPWithSMILESorINCHI(properties.get("inputFile"), csvLoadThreshold, thresholdIsotopic);
				} else {
					int mzColumn = Integer.parseInt(properties.get("mzColumn"));
					int intensColumn = Integer.parseInt(properties.get("intensColumn"));
					data = loadSpectraFromManyFiles(properties.get("inputFile"), properties.get("prefix"), fileFormat,
							properties.get("csvSpectrumHeader"), mzColumn, intensColumn, csvLoadThreshold,
							thresholdIsotopic);
				}
			}
		}
		FileWriter fw = new FileWriter(new File(outputFile));
		FileWriter fw1 = null;
		if (outputFile1 != null) {
			fw1 = new FileWriter(new File(outputFile1));
			String line = "Name,File name,SMILES,Molecular ion level,Molecular ion (or M-xH) presents?"
					+ ",Fraction above molecular,% of total ion current explained,,,,,,,,,Molecular ion explanation level";
			String line2 = ",,,,,,2 bonds broken,,,3 bonds broken,,,All possible fragments,,,,,";
			String line3 = ",,,,,,Level3,Level2,Level1,Level3,Level2,Level1,Level3,Level2,Level1,[M+],[M-H+],[M-2H+]";
			fw1.write(line + "\n" + line2 + "\n" + line3 + "\n");
		}
		for (int i = 0; i < data.length; i++) {
			String fileNameForLog = data[i].getRight().length >= 3 ? data[i].getRight()[2] : "";
			String explanation = (explainOneSpectrum(properties, fileNameForLog, fw, data[i]));
			if (fw1 != null) {
				fw1.write(explanation);
			}
			System.out.println(explanation);
		}
		fw.close();
		if (fw1 != null) {
			fw1.close();
		}
	}

	public static void main(String args[]) throws Exception {
		try {
			HashMap<String, String> properties = loadProperties(args);
			run(properties);
		} catch (Throwable e) {
			System.out.println("Exception " + e.getMessage());
			e.printStackTrace();
			printHelpShort();
		}
	}
}