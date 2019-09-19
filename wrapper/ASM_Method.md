# Algorithm and implementation details

The general pipeline is shown in
Figure [\[fig:Pipeline\]](#fig:Pipeline). We use asymptotic notation to
describe the complexity of each step. If not otherwise stated, the cost
for the worst case is presented. We use the following input parameters:

  - \(r\): number of reads

  - \(R\): total number of characters in the reads

  - \(a\): number of anchors

  - \(A\): total number of alignments

  - \(k\): anchor length

A description of each module follows.

## Setup

The basic internal data structures are created given user-defined
parameters. Mainly:

  - `reads` file: internal binary representation of the fasta/fastq
    file(s).

  - `alignments.orig` file: internal binary representation of the anchor
    alignments.

The reads in the original fasta/fastq file are renamed using ascending
numerical identifiers. Anchor sequences are extracted from the reads
following a user-defined sliding window. They are identified by the
(numerical) read id, plus the offset. Both the renamed set of reads and
the anchors are stored in chunks in fasta/fastq format (the same as the
input).

Each of the read subsets is indexed separately using the `gem indexer`.
Then, each anchor subset \(j\) is mapped to each read subset \(i\) using
the `gem mapper` with the user-defined parameters. The indexing and
mapping jobs are independent between them and can be run distributedly.

The resulting alignments are postprocessed to eliminate likely mapping
artefacts, as deletions at the ends of the anchors. Also, in case of
flagging the expected ploidy level \(l\), alignments are clustered
together according to their CIGAR string and only the first stratum-wise
\(l\) groups are kept. Finally, the alignments are transformed into the
internal binary format, which only holds the read id, offset and biggest
indel free block.

The total cost of this step is dominated by mapping and filtering.
Creating the binary files is done online using \(O(1)\) memory.

The current implementation assumes that the `reads` file can be loaded
as a whole into memory. The `alignments` file is always loaded in
chunks.

## Base graph construction.

A direct graph is built based on the sorted list of anchors projected in
one read.

We perform the following sequence of steps to build the base graph:

#### Pileup.

All the alignments are read and processed online. The information of the
alignments is used to place the reads relative to each other building up
a set of pileups. A given read is assigned only to one pileup and
position. Initially, there are as many pileups as reads. If two reads
are aligned to each other and belong to different pileups, their pileups
are joined. Otherwise, nothing is done. Thus, the alignments of a given
anchor might be spread in several positions inside a pileup. These
positions are output.

This step requires \(O(r)\) memory and \(O(A)\) time.

#### Anchor sorting.

The possible locations in which an anchor might be placed are used to
sort the (sub)anchors and then the alignments. All the anchors are
loaded into memory whilst the alignments are loaded into \(c\) chunks of
size \(O(A/c)\).

The actual sorting needs \(O(a)\) memory and \(O(a\text{log} a)\) time.
The whole alignments file needs to be read for each chunk to place the
alignments accordingly, so \(O(Ac)\) time and \(O(A/c)\) memory is
required.

#### Anchor graph. <span id="anchorgraph" label="anchorgraph">\[anchorgraph\]</span>

The sorted alignments file is processed online to create a graph of
consecutive anchors in reads. Given the processed anchors at a point,
every read stores the last anchor in which is contained. For each anchor
(node in the graph), the previous anchors according to the reads are
collected, and an edge or *link* is output between them.

This step requires \(O(r)\) memory and \(O(A)\) time.

#### Merge links. <span id="mergelinks" label="mergelinks">\[mergelinks\]</span>

An equivalent link between two anchors might have been created in
different positions of the pileup. Here, we eliminate redundant links,
by storing in memory the unique set. The processing of the links is done
by chunks. The distance comparison can be strict or accept some
divergence.

Let \(l\) be the total number of input links, and \(c\) the number of
chunks, this step requires \(O(l c)\) time and O\((l/c)\) space.

#### Block graph. <span id="blockgraph" label="blockgraph">\[blockgraph\]</span>

Unambiguous sequential anchor paths are detected and stored as a pileup.
They are called *blocks* and made up the nodes of the block graph. The
*links* show the same relationship as in the anchor graph, but
implicitly, i.e. two blocks are related if the anchors at their ends are
related. All this information is split into 4 binary files: the list of
anchors in a block (`table`), the list of positions of each anchor in
the block (`pileup`), the length and coverage information of each block
(`blocks-info`) and the links between them (`links`). Besides, the
correspondence from an anchor to a block is also stored
(`anchor2block`).

First, the anchor links are processed online, accumulating for each
anchor end how many other ends are directly related, and storing one if
found. Second, for each anchor we try to build an unambiguous path from
it and flag the blocks according they belong to a path or not. Once all
the anchors have been considered, we build up and output the `pileup`,
`table` and `blocks-info` file. Finally, the anchor links are processed
online again, and for each ambiguous anchor end, we output the
corresponding block links. Also, the correspondence of an anchor to a
block is output.

Let \(l\) be the total number of input links, this step requires
\(O(a)\) space and \(O(l)\) time.

## Tip refinement

A base graph can be iteratively refined by considering smaller anchors
at the tips and mapping against other tips.

We perform the following steps:

#### Remove pure tips.

Terminal short blocks for which exists at least another non-tip
alternative are removed. A terminal block is one with neither input or
output links. The graph is loaded in memory and processed in block order
in two rounds.

Let \(l\) be the total number of input links in the graph, this step
requires \(O(l)\) time and \(O(l)\) space.

#### Block graph.

The pruned graph is compacted using the *Block Graph* module in
Section [\[blockgraph\]](#blockgraph). Note that the resulting `pileup`
and `table` files refer to blocks instead of anchors.

#### Blocks2anchor.<span id="block2anchor" label="block2anchor">\[block2anchor\]</span>

`pileup` and `table` files of a graph \(G'\) described in terms of
blocks of a graph \(G\) are translated into `pileup` and `table` files
made up by anchors using the graph \(G\) composition. First, the anchor
content of \(G'\) is collected processing online \(G\) blocks. Then, the
contents of \(G'\) blocks are output in order.

This transformation requires \(O(a \log a)\) time and \(O(a)\) space.

#### Terminal reads and anchors. <span id="terminal" label="terminal">\[terminal\]</span>

For every terminal block end, the fragments of reads contained up to a
given user defined size are output. Besides, the sequence of new
terminal smaller anchors are output. This is done for each read fragment
rather than the consensus.

First, all the sequence is loaded into memory. The `alignments` file
cannot be loaded into memory as a whole in general. To overcome this
limitation, the algorithm proceeds iteratively only loading the
alignments relevant to a chunk of terminal blocks. First, it collects
the anchor ids in as many blocks as the user-defined memory parameter
can handle. Second, it fetches the anchor contents in file order.
Finally, for each block in the chunk, it builds up the read pileup in
the corresponding user defined region and outputs the reads fragments
and the sequence of the new anchors.

Let \(c\) the number of chunks, \(A_t\) the total number of alignments
in terminal blocks, \(a_t\) the number of anchors in terminal blocks,
this step takes \(O(R + A_t/c)\) space and \(O(c A + a_t \log a_t)\)
time.

#### Updating alignments.

The new anchors are mapped to the reads in the tips.

Then, the alignments are postprocessed analogously as in the Setup (see
Section  [1.1](#setup)) and the corresponding `alignments` file is
created.

Finally, the contents of the latter `alignments` file are merged with
the global `alignments` file. The list of read fragments in terminals is
loaded into memory. Then, the original `alignments` file is read,
suppressing the bigger anchors containing terminal read fragments.
Finally, the new smaller anchors are output.

Let \(r_t\) be the number of reads in terminal reads, merging the files
requires \(O(r_t)\) space and \(O(A + a \log r_t )\) time.

## Graph filtering

A filtered set of anchors is created removing artefacts based on the
shape of the base graph. First overlaps larger than the anchor length
(i.e., anchors inside other anchors) are removed because these are
artefacts of how the graph is built. Secondly, a set of candidate blocks
to be filtered is built such that they do not expand significantly the
anchor length and their degree is high. Then, some are selected trying
to not disconnect the graph.

Let \(b\) the number of blocks of the input graph, \(f\) the size of the
candidate set, this step takes \(O(b)\) memory and \(O(a + f k)\) time.

From the new set of anchors, an anchor graph is built and then,
compacted into a block graph analogously as in
Section [\[anchorgraph\]](#anchorgraph).

## Graph simplification

The filtered graph is simplified in several ways, removing spurious
blocks and links. Some of the transformations can be applied several
times.

We describe them in the following.

#### Clean graph.

Given an input graph, bubbles, tips and spurious direct paths are
removed if matching the user specified criteria, in particular coverage.

First, a filtering based on coverage is done, both on blocks and links.
These are removed if there are other alternatives.

Then, spurious direct links and bubbles are removed iteratively,
electing those with more coverage. A link between two blocks is
considered spurious, if there is at least a path connecting them going
through at least one block and with similar distance. To remove a bubble
or a direct link the source and destination must be overlapping or
connected by library links.

Finally, tips are removed based on coverage and length.

Let \(l\) be the total number of input links in the graph, let \(b\) the
total number of blocks, let \(d\) the biggest node cardinality, let
\(i\) the total number of iterations, let \(\delta\) the search depth
for equivalent paths, this step requires \(O(l)\) space and
\(O(i b d^\delta)\) time. Note that \(d^\delta\) is an upper bound. In
general, most of nodes have a small cardinality and the search finishes
much earlier than \(\delta\) steps for most branches.

#### Repeat resolution and scaffolding.

In case pairing information is provided (either by a paired library or
because the input to the assembly are scaffolds), the corresponding
links between blocks are deduced and the average library length
estimated. This is achieved without remapping, simply the alignments
file is read, and for each read the first and last block in which is
contained is kept. Besides, the correspondence from an anchor to a block
is used. Then, paired links are processed and in case the source and
destination correspond to different blocks, the link is output.

Let \(l\) be the total number of input links, this takes \(O(A + l)\)
time and \(O(r + a)\) space.

Then, redundant pair links between blocks are removed by using the
*Merge Links* module in Section [\[mergelinks\]](#mergelinks) and graph
\(G\) is obtained. Afterwards, the standard deviation is computed by
loading \(G\) into memory and reading the unsimplified links.

Let \(l\) be the total number of input links, \(l_G\) be the total
number of links in G, this step takes \(O(l_G + l)\) time and \(O(l_G)\)
space.

Finally, local repeat resolution and path expansion is applied on \(G\).
The summarized pairing information, together with the alignments of the
block ends are loaded into memory. From these alignments, we can deduce
the reads which expand over several blocks. Besides, \(G\) death ends
are enriched with pair links.

A repeat like block is one with more that one input and output. Besides
we impose relative length constraints on the direct neighbouring blocks
to minimize bubble expansion. The algorithm considers each repeat like
block and tries to find an unambiguous crossing path from one direct
neighbour and with enough support. The path is built by read threading
or by finding a path to a pair related block. When searching a path, at
most \(\delta\) blocks further are considered to avoid exponential
explosion. This is repeated until the graph cannot be refined further or
\(\iota\) iterations have been performed, where \(\iota\) is a user
defined parameter.

Then, the graph is simplified further by building paths from seed nodes.
Seed nodes are those with biggest block size and small cardinality.
Specifically, the algorithm tries to build paths between seed nodes
using read threading and pair information. Similarly as for the repeat
resolution, the path searching is limited to avoid exponential
explosion. Note that the non seed nodes may belong to many paths
connecting seed nodes. If no path is found from a seed node, all the non
seed nodes reachable from it are kept.

Let \(B\) be the total number of alignments at block ends, let \(l_G\)
be the total number of links in G, let \(d\) the biggest node
cardinality, this step require \((B + l)\) space and
\(O(\iota B d^\delta)\) time. Note as with Clean Graph module (see
Section [\[clean\]](#clean)), the latter is an upper bound which in
practice is not achieved.

#### Graph pruning.

An input graph is simplified using a greedy strategy that resembles
sspace scaffolder .

The graph is loaded into memory, together with basic block information.
Then, blocks are considered in descending length order. If the block has
not been visited, a path is built from each end as further as possible.
A path can be expanded while a visited node is not found, and the best
supported alternative has more coverage than the second supported
alternative by a user-defined given ratio. Once all blocks have been
considered, if a given block end could not be expanded, it is linked
back again with those alternatives that could not be expanded either.
Finally, the simplified graph is output.

Let \(l\) be the total number of input links, this step takes \(O(l)\)
time and \(O(l)\) space.

## Non-overlapping graph

Direct neighbours in the block graph can be overlapping up to the anchor
length. If the consensus from overlapping blocks were to be output,
redundant sequences will be output at the ends. If this is to be
avoided, the graph can be transformed first into a non overlapping graph
as follows.

#### Overlapping class graph.

We process the `links` file online in several rounds. First, we
determine for each block end its maximum overlap and allocate as many
classes as overlapping positions. Note that the maximum overlap is
\(k\). Besides, each non overlapping fragment is a class by itself. In
the second round, we make equivalent all the overlapping positions
(classes). Then, we output the unique classes by calculating the
transitivity relationships. Also, we output the links between contiguous
classes and correspondences between positions, classes and original
blocks. Finally, in the third round, we output the links between non
overlapping blocks in terms of their corresponding classes.

Let \(l\) be the total number of input links, \(b\) the number of blocks
of the input graph, this step takes \(O(k b)\) space and \(O(k b + l)\)
time.

Then, the output class graph is simplified, by merging links and
compacting the graph as in Section [\[blockgraph\]](#blockgraph).

Note that some small spurious blocks can be created as an artefact of
the method because equivalence relationships are build at a finer
granularity (position level) than the anchors (alignment level, where
indels are encapsulated). These blocks should be filtered out when
outputting the fasta.

#### Classes to blocks.

A class graph is translated into a *block of blocks graph*. First, the
correspondence from classes to non overlapping blocks are loaded into
memory. Then, each block \(\beta\) of the overlapping graph \(G\) is
processed online. Each of the end overlapping positions of \(\beta\) are
translated into classes, and then, to non overlapping blocks, in which
\(\beta\) is recorded. Finally, the contents of the non overlapping
blocks in terms of \(G\) are output.

Let \(c\) be the total number of classes, let \(b_o\) be the total
number of overlapping blocks in the non overlapping graph, this step
takes \(O(b_o \log b_o)\) time and \(O(c + b_o)\) space.

Finally the block of blocks graph is translated into a block (of
anchors) graph as in Section [\[block2anchor\]](#block2anchor).

## Fasta output

From each block in the graph, a *consensus* sequence is built. This is
done by deducing the pileup of read pieces corresponding to the anchors
that belong to a block.

The `alignments` file is loaded by chunks as in *Terminal reads and
anchors* in Section [\[terminal\]](#terminal). For each chunk, the
algorithm computes the consensus of each of the blocks whose data has
been fetch.

If an anchor is marked as repeated is only used for the consensus if
there is no read crossing the repeat. By default, the consensus is
obtained by concatenating the results of pieces delimited by indels. A
majority consensus is performed in the free indel regions, and a
multiple alignment using `mafft` is performed for solving the indels.
Alternatively, ambiguities can be solved by considering only the longest
reads. This is intended to be used when using as input artificial (i.e.
resulting from an assembly) in order to minimize haplotype changes.

Additionally, the pileup of reads for each block could be provided for
snp analyses, but this is not currently implemented.

Let \(c\) the number of chunks, \(A_c\) the total number of alignments
in a chunk, this step takes \(O(A_c + R)\) space, and
\(O(c A + a \log a + R)\) time.
