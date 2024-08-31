# Yeastriction webservice

Yeastriction is a tool to find CRISPR/Cas9 target sites in various yeast genomes. This is a Python/Dash implementation of the original NodeJS version. If you're looking for the NodeJS version, please visit the [original Yeastriction-NodeJS repository](https://github.com/hillstub/Yeastriction-NodeJS).

Yeastriction first extracts all possible Cas9 target sequences (20 basepairs followed by NGG) from a specified ORF and from its complementary strand. Subsequently, sequences containing 6 or more Ts are discarded as this can terminate transcription (Braglia *et al.*, 2005; Wang *et al.*, 2008). Target sequences are then tested for off-targets (an off-target is defined as a sequence with either the NGG or NAG PAM sequence and 17 or more nucleotides identical to the original 20 bp target sequence (Hsu *et al.*, 2013)) by matching the sequences against the reference genome using Bowtie (version 1) (Langmead *et al.*, 2009). If any off-target is found the original target sequence is discarded. In a next step, the AT content is calculated for the target sequence. Using the RNAfold library (essentially with the parameters `--MEA --noLP –temp=30.`) (Lorenz *et al.*, 2011) the maximum expected accuracy structure of each RNA molecule is calculated. The target sequence is also searched for the presence of restriction sites based on a default list or a user-defined list. The targets can be ranked based on presence of restriction sites (1 for containing and 0 for lacking a restriction site), AT content (1 having the highest AT-content and 0 for the lowest AT-content) and secondary structure (1 having the lowest amount of pairing nucleotides and 0 for the highest number of nucleotides involved in secondary structures (indicated by brackets)). The range for every parameter is determined per locus and used to normalize the values. Subsequently, the target sequences are ranked by summation of the score for each parameter. These ranking scores should only be used to order the targets from a single locus and not to compare targets for different loci.

## Installation

### Using Poetry

1. Clone the repository:
   ```
   git clone https://github.com/hillstub/Yeastriction.git
   cd Yeastriction
   ```

2. Install Poetry if you haven't already:
   ```
   curl -sSL https://install.python-poetry.org | python3 -
   ```

3. Install the project dependencies:
   ```
   poetry install
   ```

4. Activate the virtual environment:
   ```
   poetry shell
   ```

5. Run the application:
   ```
   python app.py
   ```

The application should now be running on `http://localhost:8050`.

### Using Docker

If you prefer to use Docker, we provide a Dockerfile and docker-compose.yaml for easy setup:

1. Make sure you have Docker and Docker Compose installed on your system.

2. Clone the repository:
   ```
   git clone https://github.com/your-username/yeastriction.git
   cd yeastriction
   ```

3. Build and run the Docker container:
   ```
   docker-compose up --build
   ```

This will build the Docker image and start the container. The application will be accessible at `http://localhost:8050`.

To stop the container, use:
```
docker-compose down
```

## Usage

1. Navigate to the main page of the application.
2. Select a strain from the dropdown menu.
3. Choose a CRISPR system.
4. Select an locus from the list.
5. The application will display potential CRISPR/Cas9 target sites for the selected ORF.

## Importing New Genomes

To import new genomes, you need to run the application with the `ALLOW_IMPORT` environment variable set to `true`:

```
ALLOW_IMPORT=true python app.py
```

Then, navigate to the import page and upload your genome files:

1. A tab-separated file (.tab) containing the ORF information.
2. A FASTA file (.fasta or .fa) containing the genome sequence.

Both files should have the same name (except for the extension) corresponding to the strain name.

## Customization

While Yeastriction was built with _Saccharomyces cerevisiae_ in mind, it can be adapted for use with other organisms. Modifications to the ORF symbol matching and name distinction can be made in the appropriate files.

## Citation

If you use Yeastriction in your research, please cite our paper:

Robert Mans, Harmen M. van Rossum, Melanie Wijsman, Antoon Backx, Niels G.A. Kuijpers, Marcel van den Broek, Pascale Daran-Lapujade, Jack T. Pronk, Antonius J.A. van Maris, Jean-Marc G. Daran (2015) CRISPR/Cas9: a molecular Swiss army knife for simultaneous introduction of multiple genetic modifications in *Saccharomyces cerevisiae*. *FEMS Yeast Research* **15**. [[PubMed](http://www.ncbi.nlm.nih.gov/pubmed/25743786)] [[FEMS Yeast Research](http://femsyr.oxfordjournals.org/content/15/2/fov004)]

## References
  * Braglia P, Percudani R & Dieci G (2005) Sequence context effects on oligo(dT) termination signal recognition by *Saccharomyces cerevisiae* RNA polymerase III. *J Biol Chem* 280: 19551–19562.
  * Hsu PD, Scott D A, Weinstein J A, et al. (2013) DNA targeting specificity of RNA-guided Cas9 nucleases. *Nat Biotechnol* 31: 827–832.
  * Langmead B, Trapnell C, Pop M & Salzberg SL (2009) Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. *Genome Biol* 10: R25.
  * Lorenz R, Bernhart SH, Höner Zu Siederdissen C, Tafer H, Flamm C, Stadler PF & Hofacker IL (2011) ViennaRNA Package 2.0. *Algorithms Mol Biol* 6: 26.
  * Wang Q & Wang L (2008) New methods enabling efficient incorporation of unnatural amino acids in yeast. *J Am Chem Soc* 130: 6066–6067.


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

