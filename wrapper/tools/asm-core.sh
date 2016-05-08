#! /bin/bash

#### 0. Manipulate parameters
# Set path if defined
if [[ -n "${bindir}" ]]; then
 	export PATH=${bindir}:$PATH
fi

# # Set anchorchunk parameter
anchorchunk=${memory}/64
debug=0 # Deactivate debug, it could be added as an external parameter

function check_and_abort {
  val=$?
  if [ $val -ne 0 ];  then 
  	exit $val
  fi
}

function enrich_anchor {
  check_and_abort; mv ${alignments_orig} ${alignments_orig}.old
  check_and_abort; gem-aaa_add-anchors-for-uncovered ${workdir}/reads ${alignments_orig}.old ${alignments_orig}
}

function debug_pileup {
	check_and_abort; gem-aaa_id-table-coherence-test ${workdir}/order.pileup >> ${logname} 
	check_and_abort; gem-aaa_id-table-coherence-test ${workdir}/order.pileup.inverse >> ${logname} 
	check_and_abort; gem-aaa_components-coherence-test ${workdir}/components >> ${logname} 
	check_and_abort; gem-aaa_align-coherence-test ${reads} ${workdir}/alignments >> ${logname} 
	check_and_abort; gem-aaa_anchor-coherence-test ${alignments_orig} ${workdir}/anchors >> ${logname} 
}

function pileup { 
	check_and_abort; gem-aaa_pileup ${reads} ${alignments_orig}  ${workdir}/order.pileup ${workdir}/components ${workdir}/read-positions ${workdir}/anchors ${workdir}/alignments >> ${logname} 
	if [ $debug -eq 1 ] ; then debug_pileup ; fi
}

function debug_dictator {
	check_and_abort; gem-aaa_id-table-coherence-test ${workdir}/order.anchors	 >> ${logname} 
	check_and_abort; gem-aaa_anchor-coherence-test ${alignments_orig} ${workdir}/distances.anchors	 >> ${logname} 
	check_and_abort; gem-aaa_offsets-coherence-test ${workdir}/components.anchors >> ${logname} 
}

function predemux {
  check_and_abort; gem-aaa_dictator ${workdir}/anchors ${workdir}/alignments  ${workdir}/order.pileup.inverse ${workdir}/components ${workdir}/order.anchors ${workdir}/distances.anchors ${workdir}/components.anchors >> ${logname} 
  check_and_abort; gem-aaa_fascist -i ${workdir} -m ${memory} >> ${logname}
  # Check correction of operation before deleting file
  if [ $? -eq 0 ];  then 
  	rm ${workdir}/alignments
  	rm ${workdir}/order.pileup.*
  	rm ${workdir}/order.anchors
  fi
}

function debug_links {
	check_and_abort; gem-aaa_links-coherence-test $1 >> ${logname} 
}

function debug_block_graph {
	check_and_abort; gem-aaa_id-table-coherence-test ${workdir}/$1-table >> ${logname} 
	check_and_abort; gem-aaa_pileup-coherence-test ${reads} ${workdir}/$1-table ${workdir}/$1-pileup >> ${logname} 
	debug_links ${workdir}/$1-block-links >> ${logname} 
}

function compact_graph {
	PREFIX=$1
	MASK=$2
	if [ $debug -eq 1 ]; then debug_links ${workdir}/$PREFIX-anchor-links; fi
	check_and_abort; gem-aaa_merge-links ${workdir}/$PREFIX-anchor-links ${workdir}/$PREFIX-merged-anchor-links  ${anchorchunk} >> ${logname}
	if [ $debug -eq 1 ]; then debug_links ${workdir}/$PREFIX-merged-anchor-links; fi
	check_and_abort; gem-aaa_block1-graph-1-1 ${alignments_orig} ${workdir}/$PREFIX-merged-anchor-links ${workdir}/$PREFIX-table ${workdir}/$PREFIX-pileup ${workdir}/$PREFIX-blocks-info ${workdir}/$PREFIX-block-links ${workdir}/$PREFIX-anchor2block $MASK 0  >> ${logname}
	if [ $debug -eq 1 ]; then debug_block_graph $PREFIX; fi
	check_and_abort; gem-aaa_merge-links ${workdir}/$PREFIX-block-links ${workdir}/$PREFIX-links  ${anchorchunk} >> ${logname}
	last_op=$?
	if [ $debug -eq 1 ]; then debug_links ${workdir}/$PREFIX-links; fi
	# Check correction of operation before deleting file
  	if [ $last_op -eq 0 ];  then 
  		rm ${workdir}/$PREFIX-merged-anchor-links
  		rm ${workdir}/$PREFIX-block-links
  	fi
}

function remove_pure_tips_all {
    PREFIX=$1
    check_and_abort; gem-aaa_remove-tips ${workdir}/$PREFIX-blocks-info ${workdir}/$PREFIX-links ${workdir}/$PREFIX-notips-expanded-links ${workdir}/$PREFIX-blocks2keep 1 ${coverage} ${trim}  >> ${logname}
    check_and_abort; gem-aaa_merge-links ${workdir}/$PREFIX-notips-expanded-links ${workdir}/$PREFIX-notips-expanded-merged-links ${anchorchunk} >> ${logname}
    check_and_abort; gem-aaa_block1-graph-1-1 ${workdir}/$PREFIX-blocks-info ${workdir}/$PREFIX-notips-expanded-merged-links ${workdir}/$PREFIX-notips-block-table ${workdir}/$PREFIX-notips-block-pileup ${workdir}/$PREFIX-notips-blocks-info ${workdir}/$PREFIX-notips-block-links ${workdir}/$PREFIX-notips-block2block ${workdir}/$PREFIX-blocks2keep 1 >> ${logname}
    check_and_abort; gem-aaa_merge-links ${workdir}/$PREFIX-notips-block-links ${workdir}/$PREFIX-notips-links  ${anchorchunk} >> ${logname}
    if [ $2 = 1 ]; then
		check_and_abort; gem-aaa_block2anchor-pileup ${alignments_orig} ${workdir}/$PREFIX-blocks-info ${workdir}/$PREFIX-table ${workdir}/$PREFIX-pileup ${workdir}/$PREFIX-anchor2block ${workdir}/$PREFIX-notips-blocks-info ${workdir}/$PREFIX-notips-block-table ${workdir}/$PREFIX-notips-block-pileup ${workdir}/$PREFIX-notips-block2block ${workdir}/$PREFIX-notips-table ${workdir}/$PREFIX-notips-pileup ${workdir}/$PREFIX-notips-anchor2block >> ${logname}
 		last_op=$?
    else
		check_and_abort; gem-aaa_block2anchor-pileup ${alignments_orig}${workdir}/$PREFIX-blocks-info ${workdir}/$PREFIX-table ${workdir}/$PREFIX-pileup ${workdir}/$PREFIX-notips-blocks-info ${workdir}/$PREFIX-notips-block-table ${workdir}/$PREFIX-notips-block-pileup  ${workdir}/$PREFIX-notips-table ${workdir}/$PREFIX-notips-pileup >> ${logname}
    	last_op=$?
    fi
   	# Check correction of operation before deleting file
   	if [ $last_op -eq 0 ];  then 
   		rm ${workdir}/$PREFIX-notips-expanded-links 
   		rm ${workdir}/$PREFIX-notips-expanded-merged-links
   		rm ${workdir}/$PREFIX-notips-block-table
   		rm ${workdir}/$PREFIX-notips-block-pileup
   		rm ${workdir}/$PREFIX-notips-block-links
   		rm ${workdir}/$PREFIX-notips-block2block
   		rm ${workdir}/$PREFIX-blocks2keep
   	fi
}

function remove_pure_tips {
    remove_pure_tips_all $1 1
}

function stats {
	PREFIX=$1
	check_and_abort; gem-aaa_graph1-stats ${workdir}/${PREFIX}-table ${workdir}/${PREFIX}-blocks-info ${workdir}/${PREFIX}-links ${workdir}/${PREFIX}-in-out-stats ${workdir}/${PREFIX}-degree-stats ${workdir}/${PREFIX}-length-stats ${workdir}/${PREFIX}-degree-length-stats	 ${workdir}/${PREFIX}-coverage-stats ${workdir}/${PREFIX}-all-lengths >> ${logname}
}

function clean_graph_all {
	PREFIX=$1
	remove_pure_tips_all $1 $2
	check_and_abort
	PREFIX=$1"-notips"
	PREFIXOUT=$1"-nobubble"
	LIB=$3
	check_and_abort; gem-aaa_clean-graph ${workdir}/$PREFIX-blocks-info ${workdir}/$PREFIX-links ${workdir}/${PREFIXOUT}0-expanded-links ${workdir}/$PREFIX-blocks2keep 1 $3 1 1 ${max_distance_diff_perc} ${anchor} 0 >> ${logname}
	
	let ARTIFACTDIS=(${anchor}+99)/100
	check_and_abort; gem-aaa_join-equivalent-blocks ${workdir}/$PREFIX-blocks-info ${workdir}/${PREFIXOUT}0-expanded-links ${workdir}/${PREFIXOUT}0-noartifact-links ${workdir}/$PREFIX-noartifact-blocks2keep 1 $ARTIFACTDIS  >> ${logname}
	check_and_abort; gem-aaa_merge-links ${workdir}/${PREFIXOUT}0-noartifact-links ${workdir}/${PREFIXOUT}0-noartifact-merged-links ${anchorchunk} >> ${logname}
	check_and_abort; gem-aaa_merge-mask-table ${workdir}/$PREFIX-blocks2keep  ${workdir}/$PREFIX-noartifact-blocks2keep ${workdir}/$PREFIX-noartifact-all-blocks2keep
#	check_and_abort; cp ${PREFIXOUT}0-expanded-links ${workdir}/${PREFIXOUT}0-noartifact-merged-links
#	check_and_abort; cp ${workdir}/$PREFIX-blocks2keep ${workdir}/$PREFIX-noartifact-all-blocks2keep
	check_and_abort; gem-aaa_block1-graph-1-1 ${workdir}/$PREFIX-blocks-info ${workdir}/${PREFIXOUT}0-noartifact-merged-links ${workdir}/${PREFIXOUT}0-block-table ${workdir}/${PREFIXOUT}0-block-pileup ${workdir}/${PREFIXOUT}0-blocks-info ${workdir}/${PREFIXOUT}0-block-links ${workdir}/${PREFIXOUT}0-block2block ${workdir}/$PREFIX-noartifact-all-blocks2keep 1 >> ${logname}	
	check_and_abort; gem-aaa_merge-links ${workdir}/${PREFIXOUT}0-block-links ${workdir}/${PREFIXOUT}0-links  ${anchorchunk} >> ${logname}
	if [ $2 = 1 ]; then
            check_and_abort; gem-aaa_block2anchor-pileup ${alignments_orig} ${workdir}/$PREFIX-blocks-info ${workdir}/$PREFIX-table ${workdir}/$PREFIX-pileup ${workdir}/$PREFIX-anchor2block ${workdir}/${PREFIXOUT}0-blocks-info ${workdir}/${PREFIXOUT}0-block-table ${workdir}/${PREFIXOUT}0-block-pileup ${workdir}/${PREFIXOUT}0-block2block ${workdir}/${PREFIXOUT}0-table ${workdir}/${PREFIXOUT}0-pileup ${workdir}/${PREFIXOUT}0-anchor2block >> ${logname}
	else
            check_and_abort; gem-aaa_block2anchor-pileup ${alignments_orig} ${workdir}/$PREFIX-blocks-info ${workdir}/$PREFIX-table ${workdir}/$PREFIX-pileup ${workdir}/${PREFIXOUT}0-blocks-info ${workdir}/${PREFIXOUT}0-block-table ${workdir}/${PREFIXOUT}0-block-pileup ${workdir}/${PREFIXOUT}0-table ${workdir}/${PREFIXOUT}0-pileup  >> ${logname}
	fi
	check_and_abort; gem-aaa_clean-graph ${workdir}/${PREFIXOUT}0-blocks-info ${workdir}/${PREFIXOUT}0-links ${workdir}/${PREFIXOUT}-expanded-links ${workdir}/${PREFIXOUT}0-blocks2keep 1 $3  ${coverage} ${scaffolds_as_input} ${max_distance_diff_perc} ${anchor} ${trim} >> ${logname}
	check_and_abort; gem-aaa_block1-graph-1-1 ${workdir}/${PREFIXOUT}0-blocks-info ${workdir}/${PREFIXOUT}-expanded-links ${workdir}/${PREFIXOUT}-block-table ${workdir}/${PREFIXOUT}-block-pileup ${workdir}/${PREFIXOUT}-blocks-info ${workdir}/${PREFIXOUT}-block-links ${workdir}/${PREFIXOUT}-block2block ${workdir}/${PREFIXOUT}0-blocks2keep 1 >> ${logname}	
	check_and_abort; gem-aaa_merge-links ${workdir}/${PREFIXOUT}-block-links ${workdir}/${PREFIXOUT}-links  ${anchorchunk} >> ${logname}
	if [ $2 = 1 ]; then
	    check_and_abort; gem-aaa_block2anchor-pileup ${alignments_orig} ${workdir}/${PREFIXOUT}0-blocks-info ${workdir}/${PREFIXOUT}0-table ${workdir}/${PREFIXOUT}0-pileup ${workdir}/${PREFIXOUT}0-anchor2block ${workdir}/${PREFIXOUT}-blocks-info ${workdir}/${PREFIXOUT}-block-table ${workdir}/${PREFIXOUT}-block-pileup ${workdir}/${PREFIXOUT}-block2block ${workdir}/${PREFIXOUT}-table ${workdir}/${PREFIXOUT}-pileup ${workdir}/${PREFIXOUT}-anchor2block >> ${logname}
		last_op=$?
	else
	    check_and_abort; gem-aaa_block2anchor-pileup ${alignments_orig} ${workdir}/${PREFIXOUT}0-blocks-info ${workdir}/${PREFIXOUT}0-table ${workdir}/${PREFIXOUT}0-pileup  ${workdir}/${PREFIXOUT}-blocks-info ${workdir}/${PREFIXOUT}-block-table ${workdir}/${PREFIXOUT}-block-pileup ${workdir}/${PREFIXOUT}-table ${workdir}/${PREFIXOUT}-pileup  >> ${logname}
		last_op=$?
	fi
	if [ $last_op -eq 0 ];  then 
		rm ${workdir}/${PREFIXOUT}0-*
		rm ${workdir}/${PREFIX}-blocks2keep 
		rm ${workdir}/${PREFIX}-noartifact-*
		rm ${workdir}/${PREFIXOUT}0-noartifact-*
		rm ${workdir}/${PREFIXOUT}0-expanded-links
		rm ${workdir}/${PREFIXOUT}-block-*
		rm ${workdir}/${PREFIXOUT}-block2block 
	fi
}

function clean_graph {
    clean_graph_all $1 "1" $2
}

function basic_demux {
	check_and_abort; gem-aaa_demux ${alignments_orig} ${workdir}/distances.anchors ${workdir}/components.anchors ${workdir}/alignments.sorted ${workdir}/components ${workdir}/base-anchor-links >> ${logname}
    check_and_abort; compact_graph "base" ""
    #clean_graph "base"
}

function filtered_demux {	
	let ANCHORPLUSINDEL=${anchor}+${max_distance_diff_perc}
	BASE=base
	check_and_abort; gem-aaa_create-filtered-anchors  ${reads} ${alignments_orig}  ${workdir}/$BASE-anchor2block ${workdir}/$BASE-blocks-info ${workdir}/$BASE-links ${workdir}/anchors-mask ${max_distance_diff_perc} ${hub_cardinality} ${coverage} 
	check_and_abort; gem-aaa_filtered-demux ${workdir}/anchors-mask ${workdir}/distances.anchors ${workdir}/components.anchors ${workdir}/alignments.sorted ${workdir}/components ${workdir}/filtered-anchor-links >> ${logname}
	if [ $? -eq 0 ];  then 
		rm ${workdir}/distances.anchors 
		rm ${workdir}/components.anchors
		rm ${workdir}/alignments.sorted 
		compact_graph "filtered" ${workdir}/anchors-mask
		check_and_abort; clean_graph "filtered" ""
	fi

}


function demux {
 basic_demux
 check_and_abort; filtered_demux
#  clean_graph "filtered"
}

function non_filtered_demux {
 basic_demux
 maskfile=${workdir}/anchors-mask 
 BASE=base
 check_and_abort; gem-aaa_create-true-mask ${workdir}/$BASE-anchor2block ${workdir}/$BASE-blocks-info ${workdir}/anchors-mask 
 check_and_abort; clean_graph "base" ""
}


function rough_scaffold {
    PREFIX=$1
	check_and_abort; gem-aaa_merge-links ${workdir}/$PREFIX-links ${workdir}/$PREFIX-pair-links ${workdir}/$PREFIX-expanded-with-pairs-links ${workdir}/$PREFIX-libraries-deviation ${coverage} ${anchor} ${anchorchunk}
	check_and_abort; gem-aaa_block1-graph-1-1 ${workdir}/$PREFIX-blocks-info ${workdir}/$PREFIX-expanded-with-pairs-links ${workdir}/$PREFIX-with-pairs-block-table ${workdir}/$PREFIX-with-pairs-block-pileup ${workdir}/$PREFIX-with-pairs-blocks-info ${workdir}/$PREFIX-with-pairs-block-links ${workdir}/$PREFIX-with-pairs-block2block $2 1 >> ${logname}	
	check_and_abort; gem-aaa_merge-links ${workdir}/$PREFIX-with-pairs-block-links ${workdir}/$PREFIX-with-pairs-links ${anchorchunk} >> ${logname}
	check_and_abort; gem-aaa_block2anchor-pileup ${alignments_orig} ${workdir}/$PREFIX-blocks-info ${workdir}/$PREFIX-table ${workdir}/$PREFIX-pileup ${workdir}/$PREFIX-anchor2block ${workdir}/$PREFIX-with-pairs-blocks-info ${workdir}/$PREFIX-with-pairs-block-table ${workdir}/$PREFIX-with-pairs-block-pileup ${workdir}/$PREFIX-with-pairs-block2block ${workdir}/$PREFIX-with-pairs-table ${workdir}/$PREFIX-with-pairs-pileup ${workdir}/$PREFIX-with-pairs-anchor2block >> ${logname}
	#clean bubbles and block
	clean_graph $PREFIX-with-pairs ""
	check_and_abort; cp ${workdir}/*-repeats-anchors ${workdir}/$PREFIX-nobubble-repeats-anchors
	if [ $? -eq 0 ];  then 
		rm ${workdir}/$PREFIX-expanded-with-pairs-links
		rm ${workdir}/$PREFIX-with-pairs-block-*
		rm ${workdir}/$PREFIX-with-pairs-block2block
	fi
}

function pair_links_from_scaffolds {
    PREFIX=$1
    check_and_abort; gem-aaa_agp2pairs ${workdir}/scaffolds.agp ${workdir}/pair-read-links ${min_anchor}
    maskfile=${workdir}/anchors-mask 
    check_and_abort; gem-aaa_add-pair-links ${reads} ${alignments_orig} $maskfile ${workdir}/pair-read-links ${workdir}/$PREFIX-blocks-info   ${workdir}/$PREFIX-links ${workdir}/$PREFIX-anchor2block ${workdir}/$PREFIX-pair-translated-links ${workdir}/$PREFIX-read-pair-translated-links ${workdir}/$PREFIX-libraries-counts  ${anchor} 1 1
    check_and_abort; gem-aaa_merge-links ${workdir}/$PREFIX-pair-translated-links ${workdir}/$PREFIX-pair-links ${anchorchunk}
    if [ $? -eq 0 ];  then
    	rm ${workdir}/pair-read-links 
    fi
}

function add_pair_links_from_scaffolds {
    PREFIX=$1
    #add pairs                                                                                      
    pair_links_from_scaffolds $PREFIX 
	rough_scaffold $PREFIX ""
}

function pair_links_from_libraries {
    PREFIX=$1
    maskfile=${workdir}/anchors-mask 
    check_and_abort; gem-aaa_add-pair-links ${reads}  ${alignments_orig} $maskfile ${workdir}/libraries ${workdir}/$PREFIX-blocks-info ${workdir}/$PREFIX-links ${workdir}/$PREFIX-anchor2block ${workdir}/$PREFIX-pair-translated-links ${workdir}/$PREFIX-read-pair-translated-links ${workdir}/$PREFIX-libraries-counts  ${anchor} 1 0
    check_and_abort; gem-aaa_merge-links ${workdir}/$PREFIX-pair-translated-links ${workdir}/$PREFIX-pair-links ${anchorchunk}
}
  

function add_pair_links_from_libraries {
    PREFIX=$1
    pair_links_from_libraries $PREFIX
    rough_scaffold $PREFIX ""
}


function repeat_resolution {
  PREFIX=$1
  if [ $2 -eq 0 ]; then
    pair_links_from_libraries $PREFIX
  else
    pair_links_from_scaffolds $PREFIX
  fi
  do_path_expansion=1
  if [ ${path_expansion_depth} -eq 0 ]; then
  	do_path_expansion=0
  elif  [ ${repeat_resolution_depth} -eq 0 ]; then
  	do_path_expansion=2
  fi
  check_and_abort; gem-aaa_links-standard-deviation ${workdir}/$PREFIX-pair-translated-links ${workdir}/$PREFIX-pair-links ${workdir}/$PREFIX-libraries-deviation
  check_and_abort; gem-aaa_read-repeat-resolution ${reads} ${workdir}/$PREFIX-table ${workdir}/$PREFIX-pileup ${workdir}/$PREFIX-links ${workdir}/$PREFIX-pair-links ${workdir}/$PREFIX-read-pair-translated-links ${workdir}/$PREFIX-blocks-info ${alignments_orig}  ${workdir}/$PREFIX-unfolded-block-table ${workdir}/$PREFIX-unfolded-block-pileup ${workdir}/$PREFIX-unfolded-links ${workdir}/$PREFIX-unfolded-blocks-info ${workdir}/$PREFIX-unfolded-block2block  ${workdir}/$PREFIX-libraries-deviation ${max_distance_diff_perc} ${coverage} ${coverage} ${min_out_scaffold} ${repeat_resolution_depth} ${path_expansion_depth} ${do_path_expansion} >> ${logname}
  check_and_abort; gem-aaa_block2anchor-pileup ${alignments_orig} ${workdir}/$PREFIX-blocks-info ${workdir}/$PREFIX-table ${workdir}/$PREFIX-pileup ${workdir}/$PREFIX-anchor2block ${workdir}/$PREFIX-unfolded-blocks-info ${workdir}/$PREFIX-unfolded-block-table ${workdir}/$PREFIX-unfolded-block-pileup ${workdir}/$PREFIX-unfolded-block2block ${workdir}/$PREFIX-unfolded-table ${workdir}/$PREFIX-unfolded-pileup ${workdir}/$PREFIX-unfolded-anchor2block ${workdir}/$PREFIX-norepeats-repeats-anchors  >> ${logname}
  check_and_abort; gem-aaa_merge-links ${workdir}/$PREFIX-unfolded-links ${workdir}/$PREFIX-unfolded-merged-links ${anchorchunk} >> ${logname}
  check_and_abort; gem-aaa_block1-graph-1-1 ${workdir}/$PREFIX-unfolded-blocks-info ${workdir}/$PREFIX-unfolded-merged-links ${workdir}/$PREFIX-norepeats-block-table ${workdir}/$PREFIX-norepeats-block-pileup ${workdir}/$PREFIX-norepeats-blocks-info ${workdir}/$PREFIX-norepeats-block-links ${workdir}/$PREFIX-norepeats-block2block 1 >> ${logname}	
  check_and_abort; gem-aaa_merge-links ${workdir}/$PREFIX-norepeats-block-links ${workdir}/$PREFIX-norepeats-lib-links ${anchorchunk} >> ${logname}
  check_and_abort; gem-aaa_block2anchor-pileup ${alignments_orig} ${workdir}/$PREFIX-unfolded-blocks-info ${workdir}/$PREFIX-unfolded-table ${workdir}/$PREFIX-unfolded-pileup ${workdir}/$PREFIX-unfolded-anchor2block ${workdir}/$PREFIX-norepeats-blocks-info ${workdir}/$PREFIX-norepeats-block-table ${workdir}/$PREFIX-norepeats-block-pileup ${workdir}/$PREFIX-norepeats-block2block ${workdir}/$PREFIX-norepeats-table ${workdir}/$PREFIX-norepeats-pileup ${workdir}/$PREFIX-norepeats-anchor2block >> ${logname}
  check_and_abort; cp ${workdir}/$PREFIX-libraries-deviation ${workdir}/$PREFIX-norepeats-nobubble-libraries-deviation
  check_and_abort; cp ${workdir}/$PREFIX-norepeats-repeats-anchors ${workdir}/$PREFIX-norepeats-nobubble-repeats-anchors
  check_and_abort; gem-aaa_compact-library-links ${workdir}/$PREFIX-norepeats-nobubble-libraries-deviation ${workdir}/$PREFIX-norepeats-lib-links   ${workdir}/$PREFIX-norepeats-links ${coverage}
   if [ $? -eq 0 ];  then
  	rm ${workdir}/$PREFIX-unfolded-*
  	rm ${workdir}/$PREFIX-norepeats-block-*
  	rm ${workdir}/$PREFIX-norepeats-block2block
  	rm ${workdir}/$PREFIX-norepeats-lib-links 
  fi
  clean_graph $PREFIX-norepeats ${workdir}/$PREFIX-libraries-deviation
}

function remove_library_links {
  PREFIX=$1
  check_and_abort; cp $PREFIX-links $PREFIX-all-links
  check_and_abort; gem-aaa_clean-library-links $PREFIX-all-links $PREFIX-links
  clean_graph $PREFIX ""
}

function solve {
  PREFIX=$1

  if [ $2 -eq 0 ]; then
    pair_links_from_libraries $PREFIX
  else
    pair_links_from_scaffolds $PREFIX
  fi

  check_and_abort; gem-aaa_select-branch-greedy ${workdir}/$PREFIX-blocks-info ${workdir}/$PREFIX-links ${workdir}/$PREFIX-libraries-counts  ${workdir}/$PREFIX-chosen-links ${path_selection_ratio} >> ${logname}
  check_and_abort; gem-aaa_block1-graph-1-1 ${workdir}/$PREFIX-blocks-info ${workdir}/$PREFIX-chosen-links ${workdir}/$PREFIX-solved-block-table ${workdir}/$PREFIX-solved-block-pileup ${workdir}/$PREFIX-solved-blocks-info ${workdir}/$PREFIX-solved-block-links ${workdir}/$PREFIX-solved-block2block 1 >> ${logname}	
  check_and_abort; gem-aaa_merge-links ${workdir}/$PREFIX-solved-block-links ${workdir}/$PREFIX-solved-links ${anchorchunk} >> ${logname}
  check_and_abort; gem-aaa_block2anchor-pileup ${alignments_orig} ${workdir}/$PREFIX-blocks-info ${workdir}/$PREFIX-table ${workdir}/$PREFIX-pileup ${workdir}/$PREFIX-anchor2block ${workdir}/$PREFIX-solved-blocks-info ${workdir}/$PREFIX-solved-block-table ${workdir}/$PREFIX-solved-block-pileup ${workdir}/$PREFIX-solved-block2block ${workdir}/$PREFIX-solved-table ${workdir}/$PREFIX-solved-pileup ${workdir}/$PREFIX-solved-anchor2block >> ${logname}
  check_and_abort; cp ${workdir}/$PREFIX-repeats-anchors ${workdir}/$PREFIX-solved-repeats-anchors
 # clean_graph $PREFIX-solved ${workdir}/$PREFIX-libraries-deviation
  if [ $? -eq 0 ];  then
  	rm ${workdir}/$PREFIX-chosen-links
  	rm ${workdir}/$PREFIX-solved-block-*
  	rm ${workdir}/$PREFIX-solved-block2block
  fi
}

function cutblocks {
    let MINTERMINAL=${min_out_scaffold} 
	PREFIX=$1
	REPEATSFILE=$2   
	check_and_abort; gem-aaa_blocks-by-position ${workdir}/$PREFIX-blocks-info ${workdir}/$PREFIX-links ${workdir}/$PREFIX-overlapping-ends-offsets ${workdir}/$PREFIX-non-overlapping-lengths  ${workdir}/$PREFIX-classes-links ${workdir}/$PREFIX-positions2classes >> ${logname}
	check_and_abort; gem-aaa_merge-links ${workdir}/$PREFIX-classes-links ${workdir}/$PREFIX-expanded-merged-classes-links ${anchorchunk}	>> ${logname}
#       check_and_abort; gem-aaa_basic-cut ${workdir}/$PREFIX-expanded-merged-classes-links  ${workdir}/$PREFIX-simplified-merged-classes-links
	check_and_abort; mv ${workdir}/$PREFIX-expanded-merged-classes-links ${workdir}/$PREFIX-simplified-merged-classes-links
	check_and_abort; gem-aaa_block1-graph-1-1 ${workdir}/$PREFIX-non-overlapping-lengths  ${workdir}/$PREFIX-simplified-merged-classes-links  ${workdir}/$PREFIX-classes-table ${workdir}/$PREFIX-classes-pileup ${workdir}/$PREFIX-non-overlapping-tmp-blocks-info ${workdir}/$PREFIX-non-overlapping-block-links ${workdir}/$PREFIX-classes2blocks 2  >> ${logname}
	check_and_abort; gem-aaa_merge-links ${workdir}/$PREFIX-non-overlapping-block-links ${workdir}/$PREFIX-non-overlapping-links ${anchorchunk}	>> ${logname}
	check_and_abort; gem-aaa_non-overlapping-blocks ${workdir}/$PREFIX-blocks-info ${workdir}/$PREFIX-overlapping-ends-offsets ${workdir}/$PREFIX-non-overlapping-lengths ${workdir}/$PREFIX-positions2classes ${workdir}/$PREFIX-non-overlapping-links ${workdir}/$PREFIX-classes2blocks  ${workdir}/$PREFIX-non-overlapping-block-table  ${workdir}/$PREFIX-non-overlapping-block-pileup ${workdir}/$PREFIX-non-overlapping-tmp-blocks-info  ${workdir}/$PREFIX-non-overlapping-blocks-info  >> ${logname}
#	check_and_abort; gem-aaa_block2anchor-pileup ${alignments_orig} ${workdir}/$PREFIX-blocks-info ${workdir}/$PREFIX-table ${workdir}/$PREFIX-pileup ${workdir}/$PREFIX-non-overlapping-blocks-info ${workdir}/$PREFIX-non-overlapping-block-table ${workdir}/$PREFIX-non-overlapping-block-pileup ${workdir}/$PREFIX-non-overlapping-table ${workdir}/$PREFIX-non-overlapping-pileup >> ${logname}
	check_and_abort; gem-aaa_block2anchor-pileup ${alignments_orig} ${workdir}/$PREFIX-blocks-info ${workdir}/$PREFIX-table ${workdir}/$PREFIX-pileup ${workdir}/$PREFIX-anchor2block ${workdir}/$PREFIX-non-overlapping-blocks-info ${workdir}/$PREFIX-non-overlapping-block-table ${workdir}/$PREFIX-non-overlapping-block-pileup  ${workdir}/$PREFIX-classes2blocks ${workdir}/$PREFIX-non-overlapping-table ${workdir}/$PREFIX-non-overlapping-pileup ${workdir}/$PREFIX-non-overlapping-anchor2block >> ${logname}
	#check_and_abort; gem-aaa_blocks2fasta ${reads} ${alignments_orig} ${workdir}/$PREFIX-non-overlapping-table ${workdir}/$PREFIX-non-overlapping-pileup ${workdir}/$PREFIX-non-overlapping-blocks-info ${workdir}/$PREFIX-sequence-by-block.fasta ${workdir}/$PREFIX-coverage-by-block.fasta ${workdir}/$PREFIX-reads-by-block ${min_out_scaffold} >> ${logname}
#	clean_graph_all $PREFIX"-non-overlapping" 0
#	PREFIX=$PREFIX"-non-overlapping-nobubble" 
	if [ $? -eq 0 ];  then
  		rm ${workdir}/$PREFIX-overlapping-*
  		rm ${workdir}/$PREFIX-expanded-*
  		rm ${workdir}/$PREFIX-classes*
  		rm ${workdir}/$PREFIX-simplified-*
  		rm ${workdir}/$PREFIX-positions2classes
  		rm ${workdir}/$PREFIX-non-overlapping-block-*
  	fi
	PREFIX=$PREFIX"-non-overlapping"
	CONTEXT=$((${anchor}/100 + 5))
	check_and_abort; gem-aaa_blocks2fasta ${reads} ${alignments_orig} ${workdir}/$PREFIX-table ${workdir}/$PREFIX-pileup ${workdir}/$PREFIX-blocks-info ${output_chunks} ${workdir}/$PREFIX-sequence-by-block.fasta ${workdir}/$PREFIX-coverage-by-block.fasta ${workdir}/$PREFIX-reads-by-block $REPEATSFILE ${min_out_scaffold} $CONTEXT $MINTERMINAL ${consensus_type}  >> ${logname}
	 if [ $? -eq 0 ];  then
	 	rm ${workdir}/$PREFIX-lengths
  	fi
}

function blocks2fasta {
    PREFIX=$1
    REPEATSFILE=$2
    MINTERMINAL=${min_out_scaffold}
    CONTEXT=$((${anchor}/100 + 5))
    check_and_abort; gem-aaa_blocks2fasta ${reads} ${alignments_orig} ${workdir}/$PREFIX-table ${workdir}/$PREFIX-pileup ${workdir}/$PREFIX-blocks-info ${output_chunks} ${workdir}/$PREFIX-sequence-by-block.fasta ${workdir}/$PREFIX-coverage-by-block.fasta ${workdir}/$PREFIX-reads-by-block $REPEATSFILE ${min_out_scaffold} $CONTEXT $MINTERMINAL ${consensus_type} >> ${logname}
}

function terminals {
  PREFIX=$1
  let region=${anchor} - 1;
  let newanchor=${anchor}/2;
  check_and_abort; gem-aaa_terminal-reads2fasta ${reads} ${alignments_orig} ${workdir}/$PREFIX-table ${workdir}/$PREFIX-pileup ${workdir}/$PREFIX-blocks-info ${workdir}/terminal-anchors.fasta ${workdir}/terminal-reads $region $newanchor >> ${logname}	
  if [ $? -eq 0 ];  then
  	rm alignments.sorted
  fi
}

function explorer {
  PREFIX=$1
  check_and_abort; gem-aaa_explorer1 ${reads} ${alignments_orig} ${workdir}/$PREFIX-table ${workdir}/$PREFIX-pileup ${workdir}/$PREFIX-blocks-info ${workdir}/$PREFIX-links
}

# ${prefix} can be base or filtered
function core {
# enrich_anchor
 pileup 
 predemux
 demux
# cutblocks ${prefix}"-nobubble"
 repeat_resolution ${prefix}"-nobubble" ${scaffolds_as_input}
 cutblocks ${prefix}"-nobubble-norepeats-nobubble" ${prefix}"-nobubble-norepeats-nobubble-repeats-anchors"
}

function scaffold {
   if  [ $2 -eq 0 ]; then
    add_pair_links_from_libraries  $1
  else
    add_pair_links_from_scaffolds $1
  fi
}


function core_tips {
 basic_demux
 remove_pure_tips "base"
 if [ $last_op -eq 0 ];  then 
 	rm ${workdir}/alignments.sorted
 	rm ${workdir}/distances.anchors 
	rm ${workdir}/components.anchors
 fi
}




#num_params=$#              # Number of command-line parameters.

#Apply results
rm -rf ${logname}

if [[ -z ${module} ]]; then 
		core
else
  if [ ${module} = "pileup" ]; then
	pileup
	predemux
  elif  [ ${module} = "simplified-graph" ]; then
	demux 
  elif  [ ${module} = "basic-graph" ]; then
	non_filtered_demux 
  elif  [ ${module} = "cutblocks" ]; then
	cutblocks ${prefix} ""
  elif  [ ${module} = "cutblocks-repeats" ]; then
#      blocks2fasta ${prefix} ${prefix}-repeats-anchors
       cutblocks ${prefix} ${prefix}-repeats-anchors
  elif  [ ${module} = "blocks2fasta" ]; then
      blocks2fasta ${prefix} ""
  elif  [ ${module} = "blocks2fasta-repeats" ]; then
  	  echo ${prefix}
      blocks2fasta ${prefix} ${prefix}-repeats-anchors  
  elif  [ ${module} = "explorer" ]; then
	explorer  ${prefix}
  elif  [ ${module} = "terminals" ]; then
	terminals  ${prefix}
  elif  [ ${module} = "repeats" ]; then
	repeat_resolution ${prefix} ${scaffolds_as_input}
  elif [ ${module} = "notips-graph" ]; then
        core_tips  
  elif [ ${module} = "core-pipeline" ]; then
        core 
  elif [ ${module} = "remove-tips" ]; then
        remove_pure_tips ${prefix}
  elif [ ${module} = "scaffold" ]; then
      scaffold ${prefix} ${scaffolds_as_input}
  elif [ ${module} = "clean" ]; then
      clean_graph ${prefix} "" 
  elif [ ${module} = "clean-repeats" ]; then
      clean_graph ${prefix}-norepeats ${prefix}-libraries-deviation
  elif [ ${module} = "solve" ]; then
      solve ${prefix}
  elif [ ${module} = "remove-paired-links" ]; then
      remove_library_links ${prefix}
  elif [ ${module} = "stats" ]; then
 	  stats ${prefix}
  else
  	echo "invalid option"
    exit 1
  fi
fi
