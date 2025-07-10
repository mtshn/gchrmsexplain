package ru.ac.phyche.gchrmsexplain;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashMap;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.ScrollPaneLayout;

import org.apache.commons.lang3.exception.ExceptionUtils;

import ru.ac.phyche.gchrmsexplain.AppBatchGUI.ExtractProperties;

public class AppSingleSpectrumGUI {

	public static void showHelp() {
		AppBatchGUI.showEditorPaneWindow("Information about software", "", true);
	}

	public static void main(String args[]) throws Exception {
		HashMap<String, String> properties = App.loadProperties(args);

		JFrame frame = new JFrame("GC-HRMS mass spectrum interpreter");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setSize(600, 400);
		JPanel panel = new JPanel();
		panel.setBounds(0, 0, 600, 1000);
		panel.setPreferredSize(new Dimension(550, 1000));
		panel.setLayout(null);
		JScrollPane scrollPane = new JScrollPane(panel);
		scrollPane.setLayout(new ScrollPaneLayout());
		scrollPane.setPreferredSize(new Dimension(500, 400));

		JLabel label1 = new JLabel("SMILES or InChI string (structure of molecule):");
		JTextField tf1 = new JTextField("");
		label1.setBounds(10, 55, 500, 18);
		tf1.setBounds(10, 75, 500, 23);
		panel.add(label1);
		panel.add(tf1);

		JLabel label2 = new JLabel("Copy your mass spectrum (ASCII m/z table) here:");
		label2.setBounds(10, 105, 500, 18);
		panel.add(label2);

		JTextArea input = new JTextArea();
		input.setBounds(10, 125, 550, 200);
		panel.add(input);

		ExtractProperties propertiesExtractMS = AppBatchGUI.addMSproperties(panel, 10, 335);

		JButton btnRun = new JButton("<html><font size=+2 color=red>Run</font></html>");
		btnRun.setBackground(Color.CYAN);
		btnRun.setBounds(10, 15, 90, 30);
		panel.add(btnRun);
		btnRun.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				HashMap<String, String> p = propertiesExtractMS.extract(properties);
				p.put("O", "ms_explain11a4bk6soinc.tmp");
				p.put("O1", "ms_explain14a4bk6s8fjh.tmp");
				p.put("inputFile", "ms_explain14a4bk6input.tmp");
				p.put("fileFormat", "TXT");
				try {
					FileWriter fw = new FileWriter(p.get("inputFile"));
					fw.write(input.getText());
					fw.close();
					p.put("SMILES", tf1.getText());
					String x = p.get("SMILES");
					try {
						ParsingNIST23.smilesToAtomContainer(x);
					} catch (Throwable e1) {
						try {
							p.put("SMILES", ParsingNIST23.inchiToSmiles(x, true));
						} catch (Throwable e2) {
							e1.printStackTrace();
							e2.printStackTrace();
							throw (new RuntimeException("Invalid SMILES or InChI " + x));
						}
					}
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
					AppBatchGUI.showEditorPaneWindow("Explanation results", results, false);
				} catch (Throwable e1) {
					AppBatchGUI.showError(ExceptionUtils.getStackTrace(e1));
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
