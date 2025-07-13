Spectrum examples are provided only to illustrate file formats.

## TXT 
A file containing m/z values, intensities, and no other words or numbers. The string must contain N m/z-intensity pairs separated by any number of semicolons, colons, commas, parentheses (i.e., ")", "("), tabs, newlines, and spaces. For example, "80.1 999 81.0 48", where 999 and 48 are intensities. No headers, additional words etc allowed. Example 
```
79.0291 500.0 
80.0261 24.0
80.0324 2.0
80.0369 999.0
81.0403 48.0
81.034 4.0 
```

## MSP
MSP file (NIST file format). One spectrum per file, structure information is ignored.  The file must begin with the Name field, the last field before the peak list must be Num peaks. See also the documentation for NIST mssearch. The header part of the file consists of tags (separated by semicolons or newlines), each of which consists of a tag name and a value separated by a colon. The Name and Num peaks tags are mandatory, the others are not needed. After the Num peaks tag there are pairs of m/z - intensity usually separated by new lines or semicolons. Tag names are not case sensitive. Example
```
Name: pyrazine
Other tag: is not required
Num peaks: 6
79.0291 500.0 ; 80.0261 24.0
80.0324 2.0
80.0369 999.0
81.0403 48.0
81.034 4.0 
```

## CSV
comma separated table with m/z and intensities. It obligatory contains header line. Everything above the header is ignored. m/z and relative intensities are in the corresponding columns (column numbers start from zero). Example (header = m/z,intensity mzColumn=0 intensColumn =1):
```
this line is above header
m/z,intensity
80.1, 999
81.0, 48
```
Another example (header = something,m/z,something else,intensity mzColumn=1 intensColumn =3)
```
this line is above header
something,m/z,something else,intensity
a,80.1,0,999
b,81.0,1,48
```

## MSP1
MSP file (NIST file format). Multiple spectra per file, structure information is ignored.  The file must begin with the Name field, the last field before the peak list must be "Num peaks". See also the documentation for NIST mssearch. Several options are supported for  how the structure can be specified. A separate field containing the word "SMILES" in the name, for example just "SMILES" or "SMILES\_can", the InChI  field, or the Comments field in which (mass bank style) comments in quotes separated by spaces, for example in the comments field such fragments "SMILES=CC=C", "SMILES\_can=CC=C", "InChI=InChI=1S/C3H6/c1-3-2/h3H,1H2,2H3". Blocks corresponding to one spectrum do not have empty lines. Blocks are separated by empty lines (after the end of the peaks of one spectrum, before the Name: field of the next one). Tag names are not case sensitive. Example

```
Name: 4-Amino-4H-1,2,4-triazole
Comments: "SMILES=Nn1cnnc1"
Num peaks: 14
40.01812 5.197597
41.01337 2.4138587
42.0212 10.6804085
43.029026 12.774507
53.01335 3.232251
54.00864 1.1544144
56.024227 1.5152494
57.032043 1.165888
69.03203 12.647626
84.037895 1.7564867
84.04297 999.0
85.04001 15.896615
85.04631 23.503286
85.051 2.4852006

nAme: 4-Amino-4H-1,2,4-triazole
cOmments: "inchi=InChI=1S/C2H4N4/c3-6-1-4-5-2-6/h1-2H,3H2"
num peaks: 14
40.018127 6.9538255
41.01337 2.3435268
42.021202 12.53235
43.029026 12.601754
53.013363 4.170679
54.00864 1.5058227
56.024227 1.4460394
57.03205 1.1557341
69.032074 19.19098
84.03789 1.857533
84.04297 999.0
85.040016 15.837611
85.04632 23.39
85.050995 2.4600854

Name: 4-Amino-4H-1,2,4-triazole;iNchi: InChI=1S/C2H4N4/c3-6-1-4-5-2-6/h1-2H,3H2;num peaks: 14;40.01812 5.197597;41.01337 2.4138587;42.0212 10.6804085;43.029026 12.774507;53.01335 3.232251;54.00864 1.1544144;56.024227 1.5152494;57.032043 1.165888;69.03203 12.647626;84.037895 1.7564867;84.04297 999.0;85.04001 15.896615;85.04631 23.503286;85.051 2.4852006

name:  4-Amino-4H-1,2,4-triazole
SMILES: Nn1cnnc1
Num peaks: 14
40.018127 6.9538255
41.01337 2.3435268
42.021202 12.53235
43.029026 12.601754
53.013363 4.170679
54.00864 1.5058227
56.024227 1.4460394
57.03205 1.1557341
69.032074 19.19098
84.03789 1.857533
84.04297 999.0
85.040016 15.837611
85.04632 23.39
85.050995 2.4600854

name:  4-Amino-4H-1,2,4-triazole
SMILES_can: Nn1cnnc1
Num peaks: 14
40.018127 6.9538255
41.01337 2.3435268
42.021202 12.53235
43.029026 12.601754
53.013363 4.170679
54.00864 1.5058227
56.024227 1.4460394
57.03205 1.1557341
69.032074 19.19098
84.03789 1.857533
84.04297 999.0
85.040016 15.837611
85.04632 23.39
85.050995 2.4600854
```

## SDF
SDF file format. Spectrum and structure in the same file, many compounds and spectra in the same file. File format equal to that used in NIST software: LIB2NIST, NIST MS Interpreter. Further information is available here:\
https://chemdata.nist.gov/mass-spc/ms-search/Library_conversion_tool.html   \
The sdf section should not contain "extra" tags (like V2000), two-dimensional coordinates are specified. The "oldest" format. At the end of each SDF block there should be a line
```
> <MASS SPECTRUM>
```
and after it (without empty lines) a regular MSP block (see above). It ends with an empty line and the line "$$$$".
Another option contains the lines
```
> <NUM PEAKS>
**

> <MASS SPECTRAL PEAKS>
```
Where instead of ** there is a real number of peaks. Then there is a mass spectrum in TXT format, an empty line and the line "$$$$". There should be **no empty line** between $$$$ and first line of the next SDF block.
Example
```
4-AMINO-4H-1,2,4-TRIAZOLE
  CDK     07062513022D

  6  6  0  0  0  0  0  0  0  0
   -2.4472   -0.1316    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9555   -0.2891    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2055   -1.5881    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2542   -1.2765    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.4200    0.2164    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0442    0.8193    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  2  6  1  0  0  0  0
>  <NAME>
4-Amino-4H-1,2,4-triazole

>  <ID>
7

>  <NUM PEAKS>
14

>  <MASS SPECTRAL PEAKS>
40.018127 6.9538255
41.01337 2.3435268
42.021202 12.53235
43.029026 12.601754
53.013363 4.170679
54.00864 1.5058227
56.024227 1.4460394
57.03205 1.1557341
69.032074 19.19098
84.03789 1.857533
84.04297 999.0
85.040016 15.837611
85.04632 23.39
85.050995 2.4600854

$$$$ 
```
Another example
```
4-AMINO-4H-1,2,4-TRIAZOLE
  CDK     07062513022D

  6  6  0  0  0  0  0  0  0  0
   -2.4472   -0.1316    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9555   -0.2891    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2055   -1.5881    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2542   -1.2765    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.4200    0.2164    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0442    0.8193    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  2  6  1  0  0  0  0
>  <MASS SPECTRUM>
Name: 4-Amino-4H-1,2,4-triazole
DB#: 7
Num peaks: 14
40.018127 6.9538255
41.01337 2.3435268
42.021202 12.53235
43.029026 12.601754
53.013363 4.170679
54.00864 1.5058227
56.024227 1.4460394
57.03205 1.1557341
69.032074 19.19098
84.03789 1.857533
84.04297 999.0
85.040016 15.837611
85.04632 23.39
85.050995 2.4600854

$$$$  
```
