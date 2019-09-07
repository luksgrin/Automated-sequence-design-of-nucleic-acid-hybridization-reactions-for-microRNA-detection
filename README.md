# Automated-sequence-design-of-nucleic-acid-hybridization-reactions-for-microRNA-detection
MicroRNA can be found in a variety of biological samples and then they represent important molecular markers for early diagnostic strategies. This work consists in an automated sequence design algorithm of nucleic acid hybridization reactions for microRNA detection. See https://riunet.upv.es/handle/10251/125058 for full work.

To employ this script, it is crucial to have a Python version 3.4 or higher, and the following modules and packages:

- Plotly for Python (https://plot.ly/python/)
- Nupack suite (http://www.nupack.org/). To use this program in the code, the wrapper NuPACK.py (present in this repository) is essential.
- ViennaRNA Package 2 (https://www.tbi.univie.ac.at/RNA/index.html)

The code can function as an importable module or as a standalone script.

If used as a standalone script, it has the following syntax:

python codeluks [INPUT]

Where INPUT can be a DNA sequence corresponding to a miRNA or a .fasta file containing various miRNA sequences.

The output is one .txt file with the initial sequence designed and the final sequences after evolution and a .html file representing the Nupack simulations of the equilibirums.

If used as an importable module, to perform the sequence generation and evolution, the function codeluks.main(NAME, MIRNA) is sufficient.

NAME will correspond to the circuit name and MIRNA to the input sequence.

The function returns a Circuit object that contains all the components' sequences, the shadow circuit's sequences, the score value and the standarized score value.
