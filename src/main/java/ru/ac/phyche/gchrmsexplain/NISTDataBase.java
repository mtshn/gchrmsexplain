package ru.ac.phyche.gchrmsexplain;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.io.FileReader;


public class NISTDataBase {
	private MassSpectrumLR[] data = null;

	public static class SearchResult implements Comparable<SearchResult> {
		public float similarity;
		public MassSpectrumLR result;

		@Override
		public int compareTo(SearchResult o) {
			return 0 - ((Float) similarity).compareTo(o.similarity);
		}
	}

	private static int[] ints(int n) {
		int[] a = new int[n];
		for (int i = 0; i < n; i++) {
			a[i] = i;
		}
		return a;
	}

	private static int[] intsrnd(int n) {
		int[] a = ints(n);
		Random rnd = new Random();
		for (int i = 0; i < n; i++) {
			int x = a[i];
			int j = rnd.nextInt(n);
			int y = a[j];
			a[i] = y;
			a[j] = x;
		}
		return a;
	}

	public NISTDataBase(String nistFileName) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(nistFileName));
			String s = br.readLine();
			s = br.readLine();
			ArrayList<String> linesList = new ArrayList<String>();
			while (s != null) {
				linesList.add(s.trim());
				s = br.readLine();
			}
			br.close();
			String[] lines = linesList.toArray(new String[linesList.size()]);
			MassSpectrumLR[] result = new MassSpectrumLR[lines.length];
			Arrays.stream(intsrnd(lines.length)).parallel().forEach(i -> {
				result[i] = new MassSpectrumLR(lines[i], true, true);
			});
			this.data = result;
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public float compoundInDatabaseMatchFactor(String smiles, MassSpectrumLR querySpectrum) {
		try {
			String inchikey = ParsingNIST23.smilesToInchiKey(smiles);
			float result = -1;
			for (int i = 0; i < data.length; i++) {
				if (data[i].getIds().inchiKeyNist.equals(inchikey)) {
					return querySpectrum.identity(data[i]);
				}
			}
			return result;
		} catch (Exception e) {
			e.printStackTrace();
			return -100;
		}
	}

	public SearchResult[] search(MassSpectrumLR query, int n) {
		SearchResult[] x = new SearchResult[data.length];
		Arrays.stream(intsrnd(x.length)).parallel().forEach(i -> {
			float mf = query.identity(data[i]);
			SearchResult r = new SearchResult();
			r.result = data[i];
			r.similarity = mf;
			x[i] = r;
		});
		Arrays.sort(x);
		SearchResult[] result = new SearchResult[n];
		for (int i = 0; i < n; i++) {
			result[i] = x[i];
		}
		return result;
	}

}
