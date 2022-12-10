# HeatMetap
## A quick tool for metagenomic geographical visualization 


HeatMetap is a python-based tool which takes a list of reads as input tagged with location data, and a reference genome, and outputs a heatmap showing the geographical dsitribution of the organism of interest.


### Inputs
- A CSV filed named "input_data.csv"
 The CSV file must have 4 columns:
  - A 'SRR' column whose name matches the name of the input read file*
  - A 'latitude' column and a 'longitude' column
  - An optional 'basepairs_count'column which contains the number of basepairs in the input read file.
 - A FASTA file which acts as the reference genome for the organism of interest
 - Corresponding FASTA files which match the 'SRR' column

*The program automatically detects paired-end reads, so suffixes like '_1' are not required in the column. The file names must also not included file extension. The script automatically parses the file, and can parse fasta files or fasta files compresse with gzip.

#### Usage
```sh
./HeatMetap.py <Input Folder> <Output Folder> <Reference Genome> <Optional Parameters>
```
Optional parameters:
--cpu: Number of CPU cores to use (default: 12)
--coverage: Minimum coverage  in percentage of the reference genome for a sample to be included (default 80)
-- identity: Minimum percent identity for a mapped read to be included (default 50)
--radius: Output heatmap radius of datapoint (default 60)
--blur: Output heatmap blur value per datapoint (default 40)

#### Example Use
```sh
./HeatMetap.py Input Output Input/Reference/1283077.fasta --cpu 12 --identity 60
```

##### Dependencies:
- pandas
- Folium
- numpy
- matplotlib

