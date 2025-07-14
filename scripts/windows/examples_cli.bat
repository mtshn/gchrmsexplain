mkdir out
explain_cli.bat -O1 ./out/out1a.csv  -O ./out/out1b.csv --fileFormat MSP --prefix ../gchrmsexplain/example_spectra --inputFile ../gchrmsexplain/example_input/input.csv
explain_cli.bat -O1 ./out/out2a.csv  -O ./out/out2b.csv --fileFormat CSV --csvSpectrumHeader "Scan Number,m/z,Intensity,Relative,Segment Number,Flags" --mzColumn 1 --intensColumn 3  --prefix ../gchrmsexplain/example_spectra --inputFile ../gchrmsexplain/example_input/input2.csv
explain_cli.bat -O1 ./out/out3a.csv  -O ./out/out3b.csv --fileFormat MSP1 --inputFile ../gchrmsexplain/example_spectra/spectra.msp
explain_cli.bat -O1 ./out/out4a.csv  -O ./out/out4b.csv --fileFormat SDF --inputFile ../gchrmsexplain/example_spectra/spectra.sdf
explain_cli.bat -O1 ./out/out5a.csv  -O ./out/out5b.csv --fileFormat SDF --inputFile ../gchrmsexplain/example_spectra/spectra1.sdf
pause

