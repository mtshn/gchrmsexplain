package ru.ac.phyche.gchrmsexplain;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.io.ByteArrayInputStream;
import java.io.InputStream;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.inchi.InChIToStructure;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import io.github.dan2097.jnainchi.InchiStatus;

/**
 * Parse manually extracted MSP and SDF files with mass spectra
 */
public class ParsingNIST23 {

	private static String replaceAlphaBeta(String s) {
		String r = s.replace(".beta.", "a").replace(".alpha.", "a");
		return r;
	}

	private static String convertNameMSPToNameSDF(String s) {
		String r = s.replace(".beta.", "a").replace(".eta.", "u").replace(".alpha.", "a").replace(".pi.", "a")
				.replace(".sigma.", "a").replace(".gamma.", "c").replace(".epsilon.", "i").replace(".delta.", "e")
				.replace(".mu.", "?").replace(".omega.", "e").replace(".+/-.", "n").replace("  #", " #");
		r = r.substring(0, Math.min(r.length(), 80)).trim();
		if (s.equals("(10S*,11R*)-(.+-.)-1.alpha.,4.alpha.,4a.alpha.,4b.beta..,5.alpha.,8.alpha."
				+ ",8b.beta.,9a.alpha.-octahydro-10,11-diphenoxy-1,4:5,8-dimethano-9H-fluoren-9-one")) {
			return ("(10S*,11R*)-(.+-.)-1a,4a,4aa,4ba.,5a,8a,8ba,9aa-octahydro-10,11-diphenoxy-");
		}
		return r;
	}

	public static IAtomContainer smilesToAtomContainer(String s) throws CDKException {
		SmilesParser parser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
		IAtomContainer mol = parser.parseSmiles(s.trim());
		Aromaticity arom = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
		arom.apply(mol);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(mol.getBuilder());
		adder.addImplicitHydrogens(mol);
		return mol;
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

	public static String smilesToInchiKey(String smiles) throws CDKException {
		IAtomContainer mol = smilesToAtomContainer(smiles);
		InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
		InChIGenerator gen = factory.getInChIGenerator(mol);
		String inchi = gen.getInchiKey();
		return inchi;
	}

	public static String smilesToInchi(String smiles) throws CDKException {
		IAtomContainer mol = smilesToAtomContainer(smiles);
		InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
		InChIGenerator gen = factory.getInChIGenerator(mol);
		String inchi = gen.getInchi();
		return inchi;
	}

	public static String inchiToSmiles(String inchi, boolean stereochemistry) throws CDKException {
		InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
		InChIToStructure intostruct = factory.getInChIToStructure(inchi, DefaultChemObjectBuilder.getInstance());
		InchiStatus ret = intostruct.getStatus();
		if (ret.equals(InchiStatus.ERROR)) {
			throw (new CDKException("Inchi status failed!"));
		}
		IAtomContainer mol = intostruct.getAtomContainer();
		return atomContainerlToSmiles(mol, stereochemistry);
	}

	public static String canonical(String smiles, boolean stereochemistry) throws CDKException {
		if ((smiles != null) && (!smiles.trim().equals(""))) {
			return inchiToSmiles(smilesToInchi(smiles), stereochemistry);
		} else {
			return null;
		}
	}

	public static void convertMSPToLines(String mspFile, String outfile) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(mspFile));
		FileWriter fw = new FileWriter(outfile);
		String s1 = br.readLine();
		fw.write("Name|SMILES|Num peaks|Spectrum|NIST ID|DB ID|InChIKey NIST|InChIKey CDK|||||||||Comment" + "\n");
		HashSet<String> compoundsWithZeroIntensity = new HashSet<String>();
		while (s1 != null) {
			boolean invalidFile = false;
			if (s1.indexOf("Name:") != 0) {
				invalidFile = true;
			}
			String name1 = s1.split("Name:")[1].trim();

			s1 = br.readLine();
			if (s1.indexOf("CDK_SMILES:") != 0) {
				invalidFile = true;
			}
			String smiles = s1.split("CDK_SMILES:")[1].trim();

			s1 = br.readLine();
			if (s1.indexOf("CDK_InChIKey:") != 0) {
				invalidFile = true;
			}
			String inchiKeyCDK = s1.split("CDK_InChIKey:")[1].trim();

			s1 = br.readLine();
			if (s1.indexOf("InChIKey:") != 0) {
				invalidFile = true;
			}
			String inchiKeyNIST = s1.split("InChIKey:")[1].trim();
			String nistID = "ERROR!!";
			String dbID = "ERROR!!";
			String comment = "No comments";

			while (s1.indexOf("Num Peaks:") != 0) {
				s1 = br.readLine();
				if ((s1.trim().equals("")) || (s1.indexOf("Name:") == 0)) {
					invalidFile = true;
				}
				if (s1.contains("NIST#:")) {
					nistID = s1.split("NIST#:")[1].trim().split("\\s+")[0];
				}
				if (s1.contains("DB#:")) {
					dbID = s1.split("DB#:")[1].trim().split("\\s+")[0];
				}
				if (s1.indexOf("Comments:") == 0) {
					comment = s1.split("Comments:")[1].trim();
				}
			}
			int numpeaks = Integer.parseInt(s1.split("Num Peaks:")[1].trim());
			String spectrum = "";
			while (!(s1.trim().equals(""))) {
				s1 = br.readLine();
				spectrum += " " + s1;
			}
			String[] spectrumIntsLines = spectrum.replace(';', ' ').trim().split("\\s+");
			int[] spectrumInts = new int[spectrumIntsLines.length];
			for (int i = 0; i < spectrumInts.length; i++) {
				spectrumInts[i] = Integer.parseInt(spectrumIntsLines[i]);
			}
			if ((spectrumInts.length != numpeaks * 2) || (!spectrum.contains("999"))) {
				br.close();
				fw.close();
				throw new RuntimeException(spectrum + name1);
			}

			s1 = br.readLine();
			if (invalidFile) {
				br.close();
				fw.close();
				throw new RuntimeException("Invalid file format " + s1 + " " + name1);
			}
			String spectrumOut = "";
			for (int i = 0; i < numpeaks * 2; i = i + 2) {
				spectrumOut += spectrumInts[i] + " " + 1.0f * spectrumInts[i + 1] + " ";
				if (spectrumInts[i + 1] == 0) {
					compoundsWithZeroIntensity.add(name1 + " " + nistID + " " + inchiKeyNIST);
				}
			}
			spectrumOut = spectrumOut.trim();
			String result = name1.trim();
			result += "|" + smiles;
			result += "|" + numpeaks;
			result += "|" + spectrumOut;
			result += "|" + nistID;
			result += "|" + dbID;
			result += "|" + inchiKeyNIST;
			result += "|" + inchiKeyCDK;
			result += "|||||||||" + comment;
			fw.write(result + "\n");
		}
		System.out.println("Compounds with zero peaks " + compoundsWithZeroIntensity.size());
		br.close();
		fw.close();
	}

	public static void checkFilesAreNotShifted(String mspFile, String sdfFile, String outfile) throws IOException {
		BufferedReader br1 = new BufferedReader(new FileReader(mspFile));
		BufferedReader br2 = new BufferedReader(new FileReader(sdfFile));
		FileWriter fw = new FileWriter(outfile);
		String s1 = br1.readLine();
		String s2 = br2.readLine();
		HashSet<String> compounds = new HashSet<String>();
		int i = 0;
		while ((s1 != null) && (s2 != null)) {
			if (s1.indexOf("Name:") == 0) {
				String name1 = s1.split("Name:")[1].trim();
				s1 = br1.readLine();
				String inchikey1 = s1.split("InChIKey:")[1].trim();
				String nameKey = name1 + " " + inchikey1;
				fw.write(name1);
				fw.write("||");
				if ((!name1.equals("Xenon")) && (!name1.equals("Argon")) && (!name1.equals("Krypton"))) {
					i++;
					compounds.add(nameKey);
					fw.write(s2.trim());
					String name2 = s2.trim();
					while (!s2.trim().equals("$$$$")) {
						s2 = br2.readLine();
					}
					s2 = br2.readLine();
					if (!convertNameMSPToNameSDF(name1).equals(replaceAlphaBeta(name2))) {
						System.out.println(i);
						System.out.println(name1);
						System.out.println(convertNameMSPToNameSDF(name1));
						System.out.println(name2);
						System.out.println(replaceAlphaBeta(name2));
						br1.close();
						br2.close();
						fw.close();
						System.exit(0);
					}
				}
				fw.write("\n");
			}
			s1 = br1.readLine();
		}
		br1.close();
		br2.close();
		fw.close();
	}

	public static void addSmilesToMSP(String mspFileName, String sdfFileName, String outFileName,
			String badMolsFileName, boolean rejectForWhichInChIIsCorrectUpToStereo, boolean noCanonicalization)
			throws IOException {
		BufferedReader sdfFile = new BufferedReader(new FileReader(sdfFileName));
		BufferedReader brNist = new BufferedReader(new FileReader(mspFileName));
		FileWriter fw1 = new FileWriter(outFileName);
		FileWriter fw2 = new FileWriter(badMolsFileName);

		String s1 = brNist.readLine();
		int same = 0;
		int different = 0;
		int stereoSame = 0;
		int wrongMW = 0;
		int i1 = 0;

		while (s1 != null) {
			String result = "";
			boolean invalidFile = false;
			if (s1.indexOf("Name:") != 0) {
				invalidFile = true;
			}
			String name1 = s1.split("Name:")[1].trim();
			result += s1 + "\n" + "XXXXXXXXXXXXXXXX\nXXXYYYYYYYXXXXXX\n";
			s1 = brNist.readLine();
			if (s1.indexOf("InChIKey:") != 0) {
				invalidFile = true;
			}
			String inchiKey2 = s1.split("InChIKey:")[1].trim();
			result += s1 + "\n";
			float mw = 0;
			float mw1 = 0;
			while (!s1.contains("NIST#:")) {
				s1 = brNist.readLine();
				result += s1 + "\n";
				if (s1.indexOf("ExactMass:") == 0) {
					mw = Float.parseFloat(s1.split("ExactMass:")[1].trim());
				}
				if (s1.indexOf("MW:") == 0) {
					mw1 = Float.parseFloat(s1.split("MW:")[1].trim());
				}
				if ((s1.trim().equals("")) || (s1.indexOf("Name:") == 0)) {
					invalidFile = true;
				}
			}
			String nistBase = s1.split("NIST#:")[1].trim().split("\\s+")[0];
			while (!(s1.trim().equals(""))) {
				s1 = brNist.readLine();
				result += s1 + "\n";
			}
			s1 = brNist.readLine();
			if (invalidFile) {
				sdfFile.close();
				brNist.close();
				fw1.close();
				fw2.close();
				throw new RuntimeException("Invalid file format " + s1);
			}

			if ((!name1.equals("Xenon")) && (!name1.equals("Argon")) && (!name1.equals("Krypton"))) {
				String molString = "";
				String s = sdfFile.readLine();
				while (!s.equals("$$$$")) {
					molString += s + "\n";
					s = sdfFile.readLine();
				}
				molString += s;
				// System.out.println(molString);
				InputStream is = new ByteArrayInputStream(molString.getBytes());
				// System.out.println(molecule.getAtomCount());
				String inchiKey1 = "ERROR_PARSING_SDF";
				String smiles = "ERROR";
				float mw2 = -1000;
				boolean badRecord = false;
				try {
					if (noCanonicalization) {
						MDLV2000Reader reader = new MDLV2000Reader(is);
						IAtomContainer molecule = (IAtomContainer) reader.read(smilesToAtomContainer(""));
						smiles = atomContainerlToSmiles(molecule, true);
						inchiKey1 = smilesToInchiKey(smiles);
						reader.close();
						mw2 = (float) MolecularFormulaManipulator.getMass(
								MolecularFormulaManipulator.getMolecularFormula(smilesToAtomContainer(smiles)),
								MolecularFormulaManipulator.MonoIsotopic);	
					} else {
						MDLV2000Reader reader = new MDLV2000Reader(is);
						IAtomContainer molecule = (IAtomContainer) reader.read(smilesToAtomContainer(""));
						smiles = canonical(atomContainerlToSmiles(molecule, true), true);
						inchiKey1 = smilesToInchiKey(smiles);
						reader.close();
						mw2 = (float) MolecularFormulaManipulator.getMass(
								MolecularFormulaManipulator.getMolecularFormula(smilesToAtomContainer(smiles)),
								MolecularFormulaManipulator.MonoIsotopic);
					}
				} catch (Throwable e) {
					badRecord = true;
					fw2.write("CDK_exception_during_conversion | " + name1 + " | NIST# " + nistBase + " | InchiKey: "
							+ inchiKey2 + "\n");
				}

				if (inchiKey1.equals(inchiKey2)) {
					same++;
				} else {
					if (inchiKey1.split("\\-")[0].equals(inchiKey2.split("\\-")[0])) {
						stereoSame++;
						if (rejectForWhichInChIIsCorrectUpToStereo) {
							if (!badRecord) {
								fw2.write("InChIKeys_(from_CDK_and_from_NIST)_are_alsmost_the"
										+ "_same,_but_differ_in_stereoshemistry_or_isotopes | " + name1 + " | NIST# "
										+ nistBase + " | InchiKey: " + inchiKey2 + " | InchiKey_CDK: " + inchiKey1
										+ " | SMILES: " + smiles + "\n");
							}
							badRecord = true;
						}
					} else {
						different++;
						if (!badRecord) {
							fw2.write("InChIKeys_from_CDK_and_from_NIST_are_DIFFERENT" + " | " + name1 + " | NIST# "
									+ nistBase + " | InchiKey: " + inchiKey2 + " | InchiKey_CDK: " + inchiKey1
									+ " | SMILES: " + smiles + "\n");
						}
						badRecord = true;
					}
				}

				float massThreshold = 0.4f;

				if (mw < 0.1) { // no correct exact mass in NIST
					mw = mw1; // Integer mass from NIST
					massThreshold = 1.1f;
				}

				if (Math.abs(mw - mw2) > massThreshold) {
					if (!badRecord) {
						fw2.write("Differetnt_molecular weight_from_NIST_and_CDK" + " | " + name1 + " | NIST# "
								+ nistBase + " | InchiKey: " + inchiKey2 + " | InchiKey_CDK: " + inchiKey1
								+ " | SMILES: " + smiles + " | MW_NIST " + mw + " | MW_CDK " + mw2 + " | MW_THREHOLD "
								+ massThreshold + "\n");
						wrongMW++;
					}
					badRecord = true;
				}

				if (i1 % 100 == 0) {
					System.out.print("same " + same / (0.01f * i1) + " ");
					System.out.print("different " + different / (0.01f * i1) + " ");
					System.out.print("stereoSame " + stereoSame / (0.01f * i1) + " ");
					System.out.println(i1);
				}
				if (!badRecord) {
					fw1.write(result.replace("XXXXXXXXXXXXXXXX", "CDK_SMILES: " + smiles).replace("XXXYYYYYYYXXXXXX",
							"CDK_InChIKey: " + inchiKey1));
				}
			} else {
				fw2.write("ARGON_OR_KYPTON_OR_XENON | " + name1 + " | NIST# " + nistBase + " | InchiKey: " + inchiKey2
						+ "\n");
			}
			i1++;
		}
		sdfFile.close();
		brNist.close();
		fw1.close();
		fw2.close();
		System.out.print("same " + same / (0.01f * i1) + " ");
		System.out.print("different " + different / (0.01f * i1) + " ");
		System.out.print("stereoSame " + stereoSame / (0.01f * i1) + " ");
		System.out.println(i1);
		System.out.print("same " + same + " ");
		System.out.print("different " + different + " ");
		System.out.print("stereoSame " + stereoSame + " ");
		if (!rejectForWhichInChIIsCorrectUpToStereo) {
			System.out.print("wrongMW " + wrongMW + " ");
		}
		System.out.println(i1);
	}

	public static void main(String[] args) throws Exception {
		String prefix = "/home/xxx/NIST23/mainlib_gcms_nist23";
		String mspFileName = prefix + ".msp";
		String sdfFileName = prefix + ".sdf";
		String namesFileName = prefix + "_names.txt";
		String mspWithSMILESFileName = prefix + "_smiles.msp";
		String badRecordsFileName = prefix + "_bad_records.txt";
		String parsedLinesFileName = prefix + ".ms";

		checkFilesAreNotShifted(mspFileName, sdfFileName, namesFileName);
		addSmilesToMSP(mspFileName, sdfFileName, mspWithSMILESFileName, badRecordsFileName, true,false);
		convertMSPToLines(mspWithSMILESFileName, parsedLinesFileName);
	}

}