package ru.ac.phyche.gchrmsexplain;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Image;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashMap;

import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JEditorPane;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.ScrollPaneLayout;

import org.apache.commons.lang3.exception.ExceptionUtils;

/**
 * SWING-based GUI for batch processing of spectra
 */
public class AppBatchGUI {
	public static abstract class ExtractProperties {
		public abstract HashMap<String, String> extract(HashMap<String, String> properties);
	}

	public static String iconPNGPath() {
		return App.fileNameToPathSameDirAsJARLocated("icon.png");
	}

	public static void showEditorWithHTMLContentFromFile(String header, String fileName) {
		try {
			JFrame frame = new JFrame(header);
			Image image = new ImageIcon(iconPNGPath()).getImage();
			frame.setIconImage(image);
			frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			frame.setSize(600, 400);
			JEditorPane p = null;

			p = new JEditorPane();
			p.setContentType("text/html");
			p.setPage(new File(fileName).toURI().toURL());
			p.setEditable(false);

			JScrollPane scrollPane = new JScrollPane(p);
			scrollPane.setLayout(new ScrollPaneLayout());
			scrollPane.setPreferredSize(new Dimension(600, 400));
			scrollPane.add(p);
			frame.add(scrollPane);
			scrollPane.setViewportView(p);
			frame.setVisible(true);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void showEditorPaneWindow(String header, String content, boolean isHTML) {
		JFrame frame = new JFrame(header);
		Image image = new ImageIcon(iconPNGPath()).getImage();
		frame.setIconImage(image);
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		frame.setSize(600, 400);
		Component p = null;
		if (isHTML) {
			p = new JEditorPane();
			((JEditorPane) p).setContentType("text/html");
			((JEditorPane) p).setText(content);
			((JEditorPane) p).setEditable(false);
		} else {
			p = new JTextArea();
			((JTextArea) p).setText(content);
			((JTextArea) p).setEditable(false);
		}
		JScrollPane scrollPane = new JScrollPane(p);
		scrollPane.setLayout(new ScrollPaneLayout());
		scrollPane.setPreferredSize(new Dimension(600, 400));
		scrollPane.add(p);
		frame.add(scrollPane);
		scrollPane.setViewportView(p);
		frame.setVisible(true);
	}

	public static void showHelp() {
		showEditorWithHTMLContentFromFile("Information about software",
				App.fileNameToPathSameDirAsJARLocated("help.html"));
	}

	public static void showHelpFileFormats() {
		showEditorWithHTMLContentFromFile("Information about software",
				App.fileNameToPathSameDirAsJARLocated("fileformats.html"));

	}

	public static void showError(String error) {
		showEditorPaneWindow("Exception", error, false);
	}

	public static ExtractProperties addMSproperties(JPanel panel, int x, int yStart) {
		JLabel labelX1 = new JLabel("Parameters of explanation of mass spectral peaks");
		labelX1.setBounds(x, yStart, 490, 18);
		panel.add(labelX1);

		JLabel label6 = new JLabel("Accuracy of mass determination, Da");
		JTextField tf6 = new JTextField("0.0006");
		label6.setBounds(x, yStart + 30, 260, 18);
		tf6.setBounds(x + 270, yStart + 30, 150, 23);
		panel.add(label6);
		panel.add(tf6);

		JLabel label7 = new JLabel("Resolution");
		JTextField tf7 = new JTextField("30000");
		label7.setBounds(x, yStart + 60, 260, 18);
		tf7.setBounds(x + 270, yStart + 60, 150, 23);
		panel.add(label7);
		panel.add(tf7);

		JLabel label8 = new JLabel("percentDifferenceForIsotopic");
		JTextField tf8 = new JTextField("10");
		tf8.setToolTipText(
				"The relative (percent!) error in the intensity of an isotopic peak at which it is considered \"perfect\"");
		label8.setBounds(x, yStart + 90, 260, 18);
		tf8.setBounds(x + 270, yStart + 90, 150, 23);
		panel.add(label8);
		panel.add(tf8);

		JLabel label9 = new JLabel("absoluteDifferenceForIsotopic");
		JTextField tf9 = new JTextField("15");
		tf9.setToolTipText(
				"The absolute (base peak = 999) error in the intensity of an isotopic peak at which it is considered \"perfect\"");
		label9.setBounds(x, yStart + 120, 260, 18);
		tf9.setBounds(x + 270, yStart + 120, 150, 23);
		panel.add(label9);
		panel.add(tf9);

		JLabel label10 = new JLabel("intensityThreshold");
		JTextField tf10 = new JTextField("5");
		tf10.setToolTipText(
				"The intensity of an isotopic peak below which we do not attempt to search for an isotopic peak (base peak = 999)");
		label10.setBounds(x, yStart + 150, 260, 18);
		tf10.setBounds(x + 270, yStart + 150, 150, 23);
		panel.add(label10);
		panel.add(tf10);

		JLabel labelX2 = new JLabel("Atom migrations");
		labelX2.setBounds(x, yStart + 180, 490, 18);
		panel.add(labelX2);

		JLabel label11 = new JLabel("maxHDrift");
		JTextField tf11 = new JTextField("2");
		tf11.setToolTipText("Maximum number of H-atoms that can migrate TO ion");
		label11.setBounds(x, yStart + 210, 260, 18);
		tf11.setBounds(x + 270, yStart + 210, 150, 23);
		panel.add(label11);
		panel.add(tf11);

		JLabel label12 = new JLabel("maxHLoss");
		JTextField tf12 = new JTextField("3");
		tf12.setToolTipText("Maximum additional loss of H-atoms");
		label12.setBounds(x, yStart + 240, 260, 18);
		tf12.setBounds(x + 270, yStart + 240, 150, 23);
		panel.add(label12);
		panel.add(tf12);

		JLabel label13 = new JLabel("maxFMigration");
		JTextField tf13 = new JTextField("2");
		tf13.setToolTipText("Maximum number of fluorine atoms that can migrate");
		label13.setBounds(x, yStart + 270, 260, 18);
		tf13.setBounds(x + 270, yStart + 270, 150, 23);
		panel.add(label13);
		panel.add(tf13);

		JLabel labelX3 = new JLabel("Explaining molecular ion peak");
		labelX3.setBounds(x, yStart + 300, 535, 18);
		panel.add(labelX3);

		JLabel label14 = new JLabel("maxHLostMI");
		JTextField tf14 = new JTextField("2");
		tf14.setToolTipText("Consider M, M-H, M-2H when considering molecular ion...");
		label14.setBounds(x, yStart + 330, 260, 18);
		tf14.setBounds(x + 270, yStart + 330, 150, 23);
		panel.add(label14);
		panel.add(tf14);

		JLabel label15 = new JLabel("fractionIsotopicThreshold");
		JTextField tf15 = new JTextField("0.96");
		tf15.setToolTipText(
				"The fraction of isotopic peak intensities that must be well explained when explaining the molecular ion peak");
		label15.setBounds(x, yStart + 360, 260, 18);
		tf15.setBounds(x + 270, yStart + 360, 150, 23);
		panel.add(label15);
		panel.add(tf15);

		JLabel labelX4 = new JLabel("General mass spectrometry settings");
		labelX4.setBounds(x, yStart + 390, 535, 18);
		panel.add(labelX4);

		JLabel label16 = new JLabel("csvLoadIntensityThreshold");
		JTextField tf16 = new JTextField("1");
		tf16.setToolTipText(
				"Threshold, peaks with intensity below which are discarded during the initial loading of the spectrum (base peak (peak with maximum intensity) = 999)");
		label16.setBounds(x, yStart + 420, 260, 18);
		tf16.setBounds(x + 270, yStart + 420, 150, 23);
		panel.add(label16);
		panel.add(tf16);

		JLabel label17 = new JLabel("thresholdGenerateIsotopic");
		JTextField tf17 = new JTextField("0.0001");
		tf17.setToolTipText(
				"Threshold, peaks with intensity below which are discarded when calculating the isotopic distribution. The intensity is calculated from the ENTIRE intensity of the isotopic distribution (not from the base peak).");
		label17.setBounds(x, yStart + 450, 260, 18);
		tf17.setBounds(x + 270, yStart + 450, 150, 23);
		panel.add(label17);
		panel.add(tf17);

		ExtractProperties propertiesExtract = new ExtractProperties() {
			@SuppressWarnings("unchecked")
			public HashMap<String, String> extract(HashMap<String, String> properties) {
				HashMap<String, String> properties1;
				properties1 = (HashMap<String, String>) properties.clone();
				properties1.put("mzThreshold", tf6.getText());
				properties1.put("resolution", tf7.getText());
				properties1.put("percentDifferenceForIsotopic", tf8.getText());
				properties1.put("absoluteDifferenceForIsotopic", tf9.getText());
				properties1.put("intensityThreshold", tf10.getText());
				properties1.put("maxHDrift", tf11.getText());
				properties1.put("maxHLoss", tf12.getText());
				properties1.put("maxFMigration", tf13.getText());
				properties1.put("maxHLostMI", tf14.getText());
				properties1.put("fractionIsotopicThreshold", tf15.getText());
				properties1.put("csvLoadIntensityThreshold", tf16.getText());
				properties1.put("thresholdGenerateIsotopic", tf17.getText());
				return properties1;
			}
		};
		return propertiesExtract;
	}

	public static void main(String args[]) throws Exception {
		HashMap<String, String> properties = App.loadProperties(args);

		JFrame frame = new JFrame("Batch GC-HRMS mass spectra interpreter");
		Image image = new ImageIcon(iconPNGPath()).getImage();
		frame.setIconImage(image);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setSize(600, 400);
		JPanel panel = new JPanel();
		panel.setBounds(0, 0, 600, 1000);
		panel.setPreferredSize(new Dimension(550, 1000));
		panel.setLayout(null);
		JScrollPane scrollPane = new JScrollPane(panel);
		scrollPane.setLayout(new ScrollPaneLayout());
		scrollPane.setPreferredSize(new Dimension(500, 400));

		// O
		JLabel label1 = new JLabel("Output file");
		JTextField tf1 = new JTextField("output.csv");
		tf1.setToolTipText("Output file name, where full interpretation of mass spectra will be written");
		JButton btn1 = new JButton("Select");
		label1.setBounds(10, 55, 170, 18);
		tf1.setBounds(280, 55, 150, 23);
		btn1.setBounds(440, 55, 80, 23);
		panel.add(label1);
		panel.add(tf1);
		panel.add(btn1);
		btn1.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				JFileChooser fc = new JFileChooser();
				if ((fc).showSaveDialog(frame) == JFileChooser.APPROVE_OPTION) {
					tf1.setText(fc.getSelectedFile().getAbsolutePath());
				}
			}
		});

		// O1
		JLabel label2 = new JLabel("Output file (short table)");
		JTextField tf2 = new JTextField("out1.csv");
		tf2.setToolTipText("Shortened output file: one molecule - one line");
		JButton btn2 = new JButton("Select");
		label2.setBounds(10, 85, 170, 18);
		tf2.setBounds(280, 85, 150, 23);
		btn2.setBounds(440, 85, 80, 23);
		panel.add(label2);
		panel.add(tf2);
		panel.add(btn2);
		btn2.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				JFileChooser fc = new JFileChooser();
				if ((fc).showSaveDialog(frame) == JFileChooser.APPROVE_OPTION) {
					tf2.setText(fc.getSelectedFile().getAbsolutePath());
				}
			}
		});

		// fileFormat
		JLabel label3 = new JLabel("Mass spectrum format");
		JComboBox<String> cb3 = new JComboBox<String>(new String[] { "CSV", "TXT", "MSP", "MSP1", "SDF" });
		cb3.setToolTipText("Format of input spectra");
		JButton btn3 = new JButton("Help");
		label3.setBounds(10, 115, 170, 18);
		cb3.setBounds(280, 115, 150, 23);
		btn3.setBounds(440, 115, 69, 23);
		panel.add(label3);
		panel.add(cb3);
		panel.add(btn3);
		btn3.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				showHelpFileFormats();
			}
		});
		// inputFile
		JLabel label4 = new JLabel("Input file");
		JTextField tf4 = new JTextField("input.csv");
		tf4.setToolTipText("Input file name");
		JButton btn4 = new JButton("Select");
		label4.setBounds(10, 145, 170, 18);
		tf4.setBounds(280, 145, 150, 23);
		btn4.setBounds(440, 145, 80, 23);
		panel.add(label4);
		panel.add(tf4);
		panel.add(btn4);
		btn4.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				JFileChooser fc = new JFileChooser();
				if ((fc).showOpenDialog(frame) == JFileChooser.APPROVE_OPTION) {
					tf4.setText(fc.getSelectedFile().getAbsolutePath());
				}
			}
		});

		// prefix
		JLabel label5 = new JLabel("Mass spectra folder (many files)");
		JTextField tf5 = new JTextField("");
		tf5.setToolTipText(
				"Prefix (path) where the spectrum files listed in the input file are located. Only for those file formats where one spectrum is one file.");
		JButton btn5 = new JButton("Select");
		label5.setBounds(10, 175, 260, 18);
		tf5.setBounds(280, 175, 150, 23);
		btn5.setBounds(440, 175, 80, 23);
		panel.add(label5);
		panel.add(tf5);
		panel.add(btn5);
		btn5.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				JFileChooser fc = new JFileChooser();
				fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				if ((fc).showDialog(frame, "Select folder") == JFileChooser.APPROVE_OPTION) {
					tf5.setText(fc.getSelectedFile().getAbsolutePath());
				}
			}
		});

		JLabel labelX5 = new JLabel("CSV file loading options");
		labelX5.setBounds(10, 685, 535, 18);
		panel.add(labelX5);

		JLabel label18 = new JLabel("csvSpectrumHeader");
		JTextField tf18 = new JTextField("m/z,Intensity");
		tf18.setToolTipText(
				"CSV table header, below which are m/z and intensity values. All lines above the header are ignored.");
		label18.setBounds(10, 715, 260, 18);
		tf18.setBounds(280, 715, 250, 23);
		panel.add(label18);
		panel.add(tf18);

		JLabel label19 = new JLabel("mzColumn");
		JTextField tf19 = new JTextField("1");
		tf19.setToolTipText(
				"The column number in the CSV table that contains the m/z values. The first column is number 1 (not zero!)!");
		label19.setBounds(10, 745, 260, 18);
		tf19.setBounds(280, 745, 150, 23);
		panel.add(label19);
		panel.add(tf19);

		JLabel label20 = new JLabel("intensColumn");
		JTextField tf20 = new JTextField("2");
		tf20.setToolTipText(
				"The column number in the CSV table that contains the intensity values. The first column is number 1 (not zero!)!  ");
		label20.setBounds(10, 775, 260, 18);
		tf20.setBounds(280, 775, 150, 23);
		panel.add(label20);
		panel.add(tf20);

		ExtractProperties propertiesExtractMS = addMSproperties(panel, 10, 205);

		ExtractProperties propertiesExtractOther = new ExtractProperties() {
			@SuppressWarnings("unchecked")
			public HashMap<String, String> extract(HashMap<String, String> properties) {
				HashMap<String, String> properties1;
				properties1 = (HashMap<String, String>) properties.clone();
				properties1.put("O", tf1.getText());
				properties1.put("O1", tf2.getText());
				properties1.put("fileFormat", (String) cb3.getSelectedItem());
				properties1.put("inputFile", tf4.getText());
				properties1.put("prefix", tf5.getText());
				properties1.put("csvSpectrumHeader", tf18.getText());
				properties1.put("mzColumn", "" + (Integer.parseInt(tf19.getText()) + 1));
				properties1.put("intensColumn", "" + (Integer.parseInt(tf20.getText()) + 1));
				return properties1;
			}
		};

		JButton btnRun = new JButton("<html><font size=+2 color=red>Run</font></html>");
		btnRun.setBackground(Color.CYAN);
		btnRun.setBounds(10, 15, 90, 30);
		panel.add(btnRun);
		btnRun.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				HashMap<String, String> p = propertiesExtractMS.extract(properties);
				p = propertiesExtractOther.extract(p);
				try {
					App.run(p);
					String results = "Brief explanation results \n";
					BufferedReader br = new BufferedReader(new FileReader(p.get("O1")));
					String s1 = br.readLine();
					while (s1 != null) {
						results += s1 + "\n";
						s1 = br.readLine();
					}
					br.close();
					results += "\n\nFull explanation results \n\n";
					BufferedReader br1 = new BufferedReader(new FileReader(p.get("O")));
					String s2 = br1.readLine();
					while (s2 != null) {
						results += s2 + "\n";
						s2 = br1.readLine();
					}
					br1.close();
					showEditorPaneWindow("Explanation results", results, false);
				} catch (Throwable e1) {
					showError(ExceptionUtils.getStackTrace(e1));
					e1.printStackTrace();
				}
			}
		});

		JButton btnHelp = new JButton("<html><font size=+1 color=blue>Help</font></html>");
		btnHelp.setBackground(new Color(220, 255, 255));
		btnHelp.setBounds(120, 15, 90, 30);
		panel.add(btnHelp);
		btnHelp.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				showHelp();
			}
		});

		scrollPane.add(panel);
		frame.add(scrollPane);
		scrollPane.setViewportView(panel);
		frame.setVisible(true);

	}
}
