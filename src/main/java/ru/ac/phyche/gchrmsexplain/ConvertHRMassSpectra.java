package ru.ac.phyche.gchrmsexplain;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringWriter;

import org.apache.commons.lang3.tuple.Pair;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.layout.StructureDiagramGenerator;

public class ConvertHRMassSpectra {

	public static void convertThermoCSVToSimpleCSV(String inputFile, String outputFile, float thresholdBy999)
			throws IOException {
		MassSpectrumHR x = MassSpectrumHR.fromThermoCSVFile(inputFile, thresholdBy999, 0);
		FileWriter fw = new FileWriter("outputFile");
		fw.write("m/z,Intensity");
		for (int i = 0; i < x.getMzs().length; i++) {
			fw.write(x.getMzs()[i] + "," + x.getIntensities()[i] + "\n");
		}
		fw.close();
	}

	public static void convertMultipleCSVToMultipleMSP(String inputCSVWithFilenames, String prefix, String outputPrefix,
			String csvSpectrumHeader, int mzColumn, int intensColumn, float threshold)
			throws IOException, CDKException {
		Pair<MassSpectrumHR, String[]>[] data = App.loadSpectraFromManyFiles(inputCSVWithFilenames, prefix, "CSV",
				csvSpectrumHeader, mzColumn, intensColumn, threshold, 0);
		for (Pair<MassSpectrumHR, String[]> x : data) {
			MassSpectrumHR s = x.getLeft();
			String name = x.getRight()[0];
			String sourceFilename = x.getRight()[2];
			String[] sp = sourceFilename.split("\\.");
			String outfileName = sp[sp.length - 1].toLowerCase().equals("csv")
					? sourceFilename.replace("." + sp[sp.length - 1], ".msp")
					: sourceFilename + ".msp";
			outfileName = new File(outputPrefix, outfileName).getAbsolutePath();
			FileWriter fw = new FileWriter(outfileName);
			fw.write("Name: " + name + "\n");
			fw.write("Num peaks: " + s.getMzs().length + "\n");
			for (int i = 0; i < s.getMzs().length; i++) {
				fw.write(s.getMzs()[i] + " " + s.getIntensities()[i] + "\n");
			}
			fw.close();
		}

	}

	public static void convertMultipleCSVToSingleMSP(String inputCSVWithFilenames, String prefix, String outputFile,
			String csvSpectrumHeader, int mzColumn, int intensColumn, float threshold)
			throws IOException, CDKException {
		Pair<MassSpectrumHR, String[]>[] data = App.loadSpectraFromManyFiles(inputCSVWithFilenames, prefix, "CSV",
				csvSpectrumHeader, mzColumn, intensColumn, threshold, 0);
		FileWriter fw = new FileWriter(outputFile);
		for (Pair<MassSpectrumHR, String[]> x : data) {
			MassSpectrumHR s = x.getLeft();
			String name = x.getRight()[0];
			String smiles = x.getRight()[1];
			fw.write("Name: " + name + "\n");
			fw.write("Comments: \"SMILES=" + smiles + "\"\n");
			fw.write("Num peaks: " + s.getMzs().length + "\n");
			for (int i = 0; i < s.getMzs().length; i++) {
				fw.write(s.getMzs()[i] + " " + s.getIntensities()[i] + "\n");
			}
			fw.write("\n");
		}
		fw.close();
	}

	public static void convertMultipleCSVToSingleSDF(String inputCSVWithFilenames, String prefix, String outputFile,
			String csvSpectrumHeader, int mzColumn, int intensColumn, float threshold, boolean secondFormat)
			throws IOException, CDKException {
		Pair<MassSpectrumHR, String[]>[] data = App.loadSpectraFromManyFiles(inputCSVWithFilenames, prefix, "CSV",
				csvSpectrumHeader, mzColumn, intensColumn, threshold, 0);
		FileWriter fw = new FileWriter(outputFile);
		int n = 0;
		for (Pair<MassSpectrumHR, String[]> x : data) {
			StringWriter sw = new StringWriter();
			MDLV2000Writer writer = new MDLV2000Writer(sw);
			MassSpectrumHR s = x.getLeft();
			String name = x.getRight()[0];
			String smiles = x.getRight()[1];
			IAtomContainer ac = ParsingNIST23.smilesToAtomContainer(smiles);
			StructureDiagramGenerator sdg = new StructureDiagramGenerator();
			sdg.setMolecule(ac);
			sdg.generateCoordinates();
			ac = sdg.getMolecule();
			writer.write(ac);
			writer.close();
			String mol2 = sw.toString();
			sw.close();
			mol2 = name.toUpperCase().substring(0, Math.min(70, name.length())) + mol2;
			mol2 = mol2.replace("M  END\n", "").replace("999 V2000", "");
			fw.write(mol2);
			n = n + 1;
			if (secondFormat) {
				fw.write(">  <NAME>\n" + name + "\n\n");
				fw.write(">  <ID>\n" + n + "\n\n");
				fw.write(">  <NUM PEAKS>\n" + s.getMzs().length + "\n\n");
				fw.write(">  <MASS SPECTRAL PEAKS>\n");
			} else {
				fw.write(">  <MASS SPECTRUM>\n");
				fw.write("Name: " + name + "\n");
				fw.write("DB#: " + n + "\n");
				fw.write("Comments: \"SMILES=" + smiles + "\"\n");
				fw.write("Num peaks: " + s.getMzs().length + "\n");
			}
			for (int i = 0; i < s.getMzs().length; i++) {
				fw.write(s.getMzs()[i] + " " + s.getIntensities()[i] + "\n");
			}
			fw.write("\n$$$$\n");
		}
		fw.close();

	}

	public static void main(String[] args) throws Exception {
		convertMultipleCSVToMultipleMSP("./sample_inputs/CSV/input1.csv", "./sample_inputs/data", "",
				"Scan Number,m/z,Intensity,Relative,Segment Number,Flags", 1, 3, 1);
		convertMultipleCSVToSingleMSP("./sample_inputs/CSV/input1.csv", "./sample_inputs/data",
				"./sample_inputs/data/spectra.msp", "Scan Number,m/z,Intensity,Relative,Segment Number,Flags", 1, 3, 1);
		convertMultipleCSVToSingleSDF("./sample_inputs/CSV/input1.csv", "./sample_inputs/data",
				"./sample_inputs/data/spectra.sdf", "Scan Number,m/z,Intensity,Relative,Segment Number,Flags", 1, 3, 1,
				true);
		convertMultipleCSVToSingleSDF("./sample_inputs/CSV/input1.csv", "./sample_inputs/data",
				"./sample_inputs/data/spectra1.sdf", "Scan Number,m/z,Intensity,Relative,Segment Number,Flags", 1, 3, 1,
				false);

	}

}
