# ASM Merger

ASM is an anchor alignment based de novo assembler. ASM Merger is an implementation of ASM designed for merging multiple assemblies or collapsing erroneous duplications caused by divergent haplotypes. It takes advantage of medium to long-range read information (a few hundreds to hundreds of kb) and tolerates high levels of heterozygosity. The method is based on local inexact alignments and can handle some amount of errors smoothly, but it is not aimed at high error rates, at least at a high coverage. Similar to an overlap-layout-consensus strategy, our approach i) computes alignments between all reads, ii) constructs and refines a graph representing them, and finally, iii) outputs a consensus sequence. For details on the assembly pipeline, see the "ASM Methods" document.

## Usage

```
usage: anchor_assembly_config [-h] --json-file <filename> [--debug-mode]
                              [--input-scaffolds <contigs_filename> [<scaffolds_agp> ...]]
                              [--paired-reads <libSize> <anchor_policy: 'front'|'ends'|'all'> <file0> <file1>]
                              [--mate-paired-reads <libSize> <anchor_policy: 'front'|'ends'|'all'> <file0> <file1>]
                              [--single-end-reads <anchor_policy: 'front'|'ends'|'all'> <filename>]
                              --asm-directory <directory> --anchor <integer>
                              [--anchor-spacing <integer>]
                              [--min-anchor <integer>] [--memory <integer>]
                              [--use-small-reads] [--coverage <integer>]
                              [--anchorsxchunk <integer>]
                              [--max-size-of-index-chunk <integer>]
                              --divergence <float_over_1>
                              [--min-mapsxanchor <integer>]
                              [--max-mapsxanchor <integer>]
                              [--extra-mapping-parameters <string>]
                              [--min-index-chunks <integer>]
                              [--max-distance-diff <float_over_1>]
                              [--trim <integer>] [--hub-cardinality <integer>]
                              [--repeat-rounds <integer>]
                              [--repeat-resolution-depth <integer>]
                              [--path-expansion-depth <integer>]
                              [--path-selection-ratio <float_over_1>]
                              [--raw-graph-output] [--output-chunks <integer>]
                              [--min-out-scaffold <integer>]
                              [--consensus-type [majority | longest]]

Create a configuration json file for the assembly pipeline

optional arguments:
  -h, --help            show this help message and exit

configuration file:
  --json-file <filename>
                        Configuration JSON to be generated. Required
  --debug-mode          If present, does not remove intermediate files.
                        Intended for advance and debug usage.

input files:
  --input-scaffolds <contigs_filename> [<scaffolds_agp> ...]
                        Input assembled scaffolds to be merged description:
                        input contigs in fasta or fastq format and optinally,
                        the input agp file must be given
  --paired-reads <libSize> <anchor_policy: 'front'|'ends'|'all'> <file0> <file1>
                        Input pair end library description: the expected
                        insert size, anchor spacing policy and input
                        fasta/fastq files must be given
  --mate-paired-reads <libSize> <anchor_policy: 'front'|'ends'|'all'> <file0> <file1>
                        Input mate pair library description: the expected
                        insert size, anchor spacing policy and input
                        fasta/fastq files must be given
  --single-end-reads <anchor_policy: 'front'|'ends'|'all'> <filename>
                        Input single end library description: anchor spacing
                        policy and input fasta/fastq file must be given
  --asm-directory <directory>
                        Output directory. It should not exist. Required

main parameters:
  --anchor <integer>    Main anchor size. Required
  --anchor-spacing <integer>
                        Spacing between anchors. Default: anchor size
  --min-anchor <integer>
                        Min anchor size. Default: anchor size
  --memory <integer>    Max memory available for use. Default: 48000000000
  --use-small-reads     Consider smaller reads than the base anchor size.
                        Default: not activated
  --coverage <integer>  Minimum coverage to not remove a tip or a singleton.
                        Default: 1

mapping parameters:
  --anchorsxchunk <integer>
                        Number of anchors (i.e reads) per alignment job
  --max-size-of-index-chunk <integer>
                        Maximum size of indexing chunk. Default: 10000000000
  --divergence <float_over_1>
                        Maximum divergence to collapse. Required
  --min-mapsxanchor <integer>
                        Minimum number of allowed mappings per anchor.
                        Default: 0
  --max-mapsxanchor <integer>
                        Maximum number of allowed mappings per anchor.
                        Default: 1000
  --extra-mapping-parameters <string>
                        Extra GEM mapping parameters, unparsed. Default:
                        (None)
  --min-index-chunks <integer>
                        Chunk indexes to do the read selection, this should be
                        the minimal possible. Default: 1

graph parameters:
  --max-distance-diff <float_over_1>
                        Maximum allowed distance difference in percentage to
                        allow them to be flagged as equivalent. Default: 0.05
  --trim <integer>      Maximum trim without considering coverage. Default:
                        minimum anchor size
  --hub-cardinality <integer>
                        Node cardinality to be considered a hub to do
                        aggressive filtering. Value 0 deactivates the
                        filtering. Default: 16
  --repeat-rounds <integer>
                        Number of repeat rounds. Default: 1
  --repeat-resolution-depth <integer>
                        Search depth for repeat resolution. Value 0
                        deactivates repeat resolution. Default: 4
  --path-expansion-depth <integer>
                        Path expansion depth (scaffold). Value 0 deactivates
                        path expansion. Default: 10
  --path-selection-ratio <float_over_1>
                        Minimum coverage ratio to select a path, with respect
                        to another. If it is 0, this simplification is not
                        run. Default: 0

output parameters:
  --raw-graph-output    Outputs the graph without applyting the non-
                        overlapping transformation. Default: store_false
  --output-chunks <integer>
                        Number of output chunks to compute the multiple
                        alignment. Default: 1
  --min-out-scaffold <integer>
                        Minimum output size for scaffolds. Default: 200
  --consensus-type [majority | longest]
                        Consensus type for output: either majority consensus
                        or the longest contig is chosen. Default: majority

```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.
