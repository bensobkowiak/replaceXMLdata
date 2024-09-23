# replaceXMLdata
Replicate XML files from BEAUti with new sequence data and dates while retaining parameters to run with [BEAST2](https://www.beast2.org).

First create an XML file using BEAUti with models and prior parameters to be replicated. This will be the base XML file. Run this script to make a new XML file with new data (sequences and dates) with the same. This reduces the need to create new XML files in BEAUti for each dataset, which can be tedious!

The script can be run with Rscript using the following command:

```R
Rscript replaceData_XML.R -x base.xml -f newSequence.fasta -d Dates.txt -p old_prefix -n new_prefix -o new_data.xml -i TRUE
```

| Option  | Description |
| ------------- | ------------- |
| -x, --xml  | XML file built using BEAUti that acts as a base XML for the new file |
| -f, --fasta  | New sequence data in the FASTA format, this must be the same type (e.g., all nucleotide, all amino acid) as the base XML  |
| -d, --dates  | A two-column text file, with header, with the first column as the sequence name and second as the date. Dates must be in the same format as the base XML  |
| -p, --old_prefix | The prefix for output files that can be replaced in the base XML, this is usually the name of the alignment file used to build the base XML |
| -n, --new_prefix  | A new prefix to name all the output files |
| -o, --output  | Name of the new XML file |
| -i, --invariant  | Logical, if TRUE, this will add a line under the sequence data in the XML file to estimate the proportion of invariant sites if you are using variant only sequence data, as detailed [here](https://groups.google.com/g/beast-users/c/DuhdMp9JNcA). This script is configured for the Mtb H37Rv reference strain, but it can be adapted for any taxa. The original XML supplied withe -x option must not have this line present. |


_Note: this script is created for BEAST2 and may not be compatible with BEAST v.1. This has also only been tested with one data partition._
