package ru.ac.phyche.gchrmsexplain;

import java.util.HashMap;

public class MassSpectrumLR {//Low resolution

	public static class SpectrumIDs {
		public String name = null;
		public String smiles = null;
		public String idNIST = null;
		public String idDB = null;
		public String inchiKeyNist = null;
	}

	private int[] mzs = null;
	private float[] intensities = null;
	private float[] scaledIntensities = null;
	private float[] sqrtIntensities = null;
	private SpectrumIDs ids;

	public MassSpectrumLR(int[] mzs, float[] intensities) {
		if (mzs.length != intensities.length) {
			throw new RuntimeException("MZs and intensities arrays must have identical length");
		}
		this.mzs = mzs;
		this.intensities = intensities;
	}

	public MassSpectrumLR(String spectrum, boolean fullString, boolean performScale) {
		String spectrumLine = null;
		String[] splt0 = null;
		int numpeaksInFile = -1;
		if (fullString) {
			splt0 = (spectrum.trim()+" ").split("\\|");
			this.ids=new SpectrumIDs();
			ids.name = splt0[0].trim();
			ids.smiles = splt0[1].trim();
			numpeaksInFile = Integer.parseInt(splt0[2].trim());
			ids.idNIST = splt0[4].trim();
			ids.idDB = splt0[5].trim();
			ids.inchiKeyNist = splt0[6].trim();
		}

		if (!fullString) {
			spectrumLine = spectrum;
		} else {
			spectrumLine = splt0[3];
		}
		String[] splt1 = spectrumLine.trim().split("\\s+");
		if (splt1.length % 2 != 0) {
			throw new RuntimeException("Incorrect spectrum line. Odd number of values.");
		}
		int numpeaks = splt1.length / 2;
		if (fullString) {
			if (numpeaks != numpeaksInFile) {
				throw new RuntimeException("Wrong number of peaks. "+numpeaks+" "+numpeaksInFile);
			}
		}
		this.mzs = new int[numpeaks];
		this.intensities = new float[numpeaks];
		for (int i = 0; i < numpeaks; i++) {
			mzs[i] = Integer.parseInt(splt1[2 * i]);
			intensities[i] = Float.parseFloat(splt1[2 * i + 1]);
		}
		if (performScale) {
			this.scaling();
		}
	}

	private static String trim(String s) {
		if (s != null) {
			return s.trim();
		} else {
			return null;
		}
	}

	private static String nullToEmpty(String s) {
		if (s != null) {
			return s;
		} else {
			return "";
		}
	}

	public static MassSpectrumLR fromMSP(String s) {
		String[] splt = s.trim().split("[\\n;]");
		HashMap<String, String> fields = new HashMap<String, String>();
		SpectrumIDs ids = new SpectrumIDs();
		if (!splt[0].toUpperCase().split(":")[0].equals("NAME")) {
			throw new RuntimeException("Invalid MSP. First line should be \"Name: ...\" " + splt[0]);
		}
		int j = 0;
		while (!splt[j].toUpperCase().split(":")[0].equals("NUM PEAKS")) {
			String[] splt2 = splt[j].toUpperCase().trim().split(":");
			fields.put(splt2[0].trim(), splt2[1].trim());
			j++;
			if (j == splt.length) {
				throw new RuntimeException("Invalid MSP. No \"Num peaks: ...\" field" + splt[0]);
			}
		}
		ids.name = trim(fields.get("NAME"));
		ids.inchiKeyNist = trim(fields.get("INCHIKEY"));
		ids.idDB = trim(fields.get("DB#"));
		ids.idNIST = trim(fields.get("NIST#"));
		for (String s2 : fields.keySet()) {
			if (s2.contains("SMILES")) {
				ids.smiles = fields.get(s2.trim());
			}
		}
		int numpeaks = Integer.parseInt(splt[j].toUpperCase().split(":")[1].trim());
		j++;
		String spectrum = "";
		while (j < splt.length) {
			spectrum = spectrum += " " + splt[j];
			j++;
		}
		spectrum = spectrum.trim();
		if (spectrum.split("\\s+").length != 2 * numpeaks) {
			throw new RuntimeException("Wrong number of peaks.");
		}
		MassSpectrumLR result = new MassSpectrumLR(spectrum, false, false);
		result.ids = ids;
		return result;
	}

	public String toString() {
		String result = "";
		for (int i = 0; i < this.numPeaks(); i++) {
			result = result + mzs[i] + " " + intensities[i] + " ";
		}
		return result.trim();
	}

	public String toStringFull() {
		String result = nullToEmpty(this.getName()) + "|" + nullToEmpty(this.getSmiles()) + "|";
		result += this.numPeaks() + "|" + this.toString() + "|";
		result += nullToEmpty(this.ids.idNIST) + "|" + nullToEmpty(this.ids.idDB) + "|";
		result += nullToEmpty(this.ids.inchiKeyNist) + "||||||||||";
		return result.trim();
	}

	public void scaling() {
		sqrtIntensities = new float[this.numPeaks()];
		scaledIntensities = new float[this.numPeaks()];
		for (int i = 0; i < this.numPeaks(); i++) {
			sqrtIntensities[i] = (float) Math.sqrt(intensities[i]);
			scaledIntensities[i] = (float) (sqrtIntensities[i] * Math.sqrt(mzs[i]));
		}
	}

	public int getMZ(int i) {
		return mzs[i];
	}

	public float getIntensity(int i) {
		return intensities[i];
	}

	public float getSqrtIntensity(int i) {
		return sqrtIntensities[i];
	}

	public float getScaledIntensity(int i) {
		return scaledIntensities[i];
	}

	public int[] getMzs() {
		return mzs;
	}

	public void setMzs(int[] mzs) {
		this.mzs = mzs;
	}

	public float[] getIntensities() {
		return intensities;
	}

	public void setIntensities(float[] intensities) {
		this.intensities = intensities;
	}

	public float[] getScaledIntensities() {
		return scaledIntensities;
	}

	public void setScaledIntensities(float[] scaledIntensities) {
		this.scaledIntensities = scaledIntensities;
	}

	public float[] getSqrtIntensities() {
		return sqrtIntensities;
	}

	public void setSqrtIntensities(float[] sqrtIntensities) {
		this.sqrtIntensities = sqrtIntensities;
	}

	public void setIds(SpectrumIDs ids) {
		this.ids = ids;
	}

	// THIS - unknown. Second - library
	public float identity(MassSpectrumLR librarySpectrum) {
		float mf = IdentitySearch.calcMatchFactor(this, librarySpectrum, true, -1, -1, false);
		return mf;
	}

	public SpectrumIDs getIds() {
		return this.ids;
	}

	public String getName() {
		return this.ids.name;
	}

	public String getSmiles() {
		return this.ids.smiles;
	}

	public int numPeaks() {
		if (mzs.length != intensities.length) {
			throw new RuntimeException("MZs and intensities arrays must have identical length");
		}
		return mzs.length;
	}

}
