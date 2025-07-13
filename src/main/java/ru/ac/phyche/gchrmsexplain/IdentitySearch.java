package ru.ac.phyche.gchrmsexplain;

/**
 * This class is an implementation of the Identity algorithm for comparing
 * low-resolution mass spectra, which is used in the NIST mssearch software.
 * This class is an implementation of the Identity algorithm for comparing
 * low-resolution mass spectra, which is used in the NIST mssearch software. The
 * implementation is based on a publicly available implementation
 * (https://github.com/AndreySamokhin/mssearchr/**, MIT license), see also the
 * article Samokhin A. The Identity Algorithm: How the Most Popular Electron
 * Ionization Mass Spectral Library Search Engine Actually Works. Journal of the
 * American Society for Mass Spectrometry. 2024 Nov 18;35(12):3178-83.
 * 10.1021/jasms.4c00350
 * 
 */
public class IdentitySearch {

	/**
	 * This method implements algorithm from paper Samokhin A. The Identity
	 * Algorithm: How the Most Popular Electron Ionization Mass Spectral Library
	 * Search Engine Actually Works. Journal of the American Society for Mass
	 * Spectrometry. 2024 Nov 18;35(12):3178-83. 10.1021/jasms.4c00350 The code is
	 * java edition of the following code:
	 * https://github.com/AndreySamokhin/mssearchr/ (MIT licensed)
	 * 
	 * @param msUnknown       mass spectrum on unknown compounds
	 * @param msLibrary       mass spectrum of library compound
	 * @param isIdentity      if True Identity (composite) similarity algorithm is
	 *                        used, else - dot product similarity
	 * @param minMZ           minimum m/z considered (including)
	 * @param maxMZ           minimum m/z considered (excluding)
	 * @param isReverseSearch if reverse searc
	 * @return similarity score
	 */
	public static float calcMatchFactor(MassSpectrumLR msUnknown, MassSpectrumLR msLibrary, boolean isIdentity,
			int minMZ, int maxMZ, boolean isReverseSearch) {

// Rprintf("tol = %e\n", tol);

// counters
		int iU = 0;
		int iL = 0;
		final int nPeaksU = msUnknown.numPeaks();
		final int nPeaksL = msLibrary.numPeaks();

// finding starting m/z value
		if (minMZ > 0) {
			while ((iU < nPeaksU) && (msUnknown.getMZ(iU) < minMZ)) {
				iU++;
			}
			while ((iL < nPeaksL) && (msLibrary.getMZ(iL) < minMZ)) {
				iL++;
			}
		} else {
			if (msUnknown.getMZ(iU) > msLibrary.getMZ(iL)) {
				while ((iL < nPeaksL) && (msUnknown.getMZ(iU) > msLibrary.getMZ(iL))) {
					iL++;
				}
			} else {
				while ((iU < nPeaksU) && (msUnknown.getMZ(iU) < msLibrary.getMZ(iL))) {
					iU++;
				}
			}
		}
		// at least one mass spectrum does not contain peaks in the set m/z range
		if ((iU >= nPeaksU) || (iL >= nPeaksL)) {
			return 0;
		}

// Identity algorithm:
//   mf = (term1 * n1 + term2 * n2) / (n1 + n2)
//   term1 = cos_a^2 = sum_ul^2 / (sum_uu * sum_ll)
//   term2 = r = sum_rm / sum_m

// Similarity algorithm:
//   mf = cos_a^2 = sum_ul^2 / (sum_uu * sum_ll)

		float sum_ul = 0;
		float sum_uu = 0;
		float sum_ll = 0;
		float sum_rm = 0;
		int sum_m = 0;

		int n1 = 0; // the number of peaks in both library and unknown spectrum
		int n2 = 0; // the number of ratios of intensities

		float prev_w2_u = 0;
		float prev_w2_l = 0;

		boolean should_calc_term2 = false;
		boolean should_continue = true;

		final int kBoth = 0; // Calculation type
		final int kUnknown = 1;
		final int kLibrary = 2;

		while (should_continue) {

			int calc_type;
			if (iU >= nPeaksU) {
				calc_type = kLibrary;
			} else {
				if (iL >= nPeaksL) {
					calc_type = kUnknown;
				} else {
					if (msUnknown.getMZ(iU) > msLibrary.getMZ(iL)) {
						calc_type = kLibrary;
					} else {
						if (msUnknown.getMZ(iU) < msLibrary.getMZ(iL)) {
							calc_type = kUnknown;
						} else {
							if (isReverseSearch && (msLibrary.getIntensity(iL) == 0)) {
								calc_type = kUnknown;
							} else {
								calc_type = kBoth;
							}
						}
					}
				}
			}

			if (calc_type == kBoth) {
				if ((msUnknown.getIntensity(iU) > 1) || (msLibrary.getIntensity(iL) > 1)) {
					sum_ul += msUnknown.getScaledIntensity(iU) * msLibrary.getScaledIntensity(iL);
					sum_uu += msUnknown.getScaledIntensity(iU) * msUnknown.getScaledIntensity(iU);
					sum_ll += msLibrary.getScaledIntensity(iL) * msLibrary.getScaledIntensity(iL);
					if (msUnknown.getScaledIntensity(iU) > 0 && msLibrary.getScaledIntensity(iL) > 0) {
						n1++;
					}

					if (isIdentity) {
						if (should_calc_term2) {
							float r_num = msUnknown.getSqrtIntensity(iU) * prev_w2_l;
							float r_denom = prev_w2_u * msLibrary.getSqrtIntensity(iL);
							if ((r_num * r_denom) > 0) {
								if (r_denom > r_num)
									sum_rm += msUnknown.getMZ(iU) * r_num / r_denom;
								else
									sum_rm += msUnknown.getMZ(iU) * r_denom / r_num;
								sum_m += msUnknown.getMZ(iU);
								n2++;
							}
						}
						should_calc_term2 = true;
						prev_w2_u = msUnknown.getSqrtIntensity(iU);
						prev_w2_l = msLibrary.getSqrtIntensity(iL);
					}
				}
				iU++;
				iL++;

			} else if (calc_type == kUnknown) {
				if (msUnknown.getIntensity(iU) > 1) {
					should_calc_term2 = false;
					if (!isReverseSearch) {
						sum_uu += msUnknown.getScaledIntensity(iU) * msUnknown.getScaledIntensity(iU);
					}
				}
				iU++;

			} else { // i.e., 'calc_type == kLibrary'
				if (msLibrary.getIntensity(iL) > 1) {
					should_calc_term2 = false;
					sum_ll += msLibrary.getScaledIntensity(iL) * msLibrary.getScaledIntensity(iL);
				}
				iL++;
			}

			if (iU >= nPeaksU && iL >= nPeaksL) {
				should_continue = false;
			}

			if (maxMZ > 0) {
				if ((iU >= nPeaksU || (msUnknown.getMZ(iU) > maxMZ))
						&& (iL >= nPeaksL || (msLibrary.getMZ(iL) > maxMZ))) {
					should_continue = false;
				}
			}
		}

		if (n1 == 0) {
			return 0;
		}

		float term1 = sum_ul * sum_ul / (sum_uu * sum_ll); // i.e., cos_a^2;
		if (sum_m > 0) { // if 'sum_m = 0', then 'n2 = 0' and 'term2 = 0'
			return 1000.0f * (term1 * n1 + (sum_rm / sum_m) * n2) / (n1 + n2) - 0.5f;
		} else {
			return 1000.0f * term1 - 0.5f;
		}
	}
}
