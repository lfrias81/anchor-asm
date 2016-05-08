{
    "env": {
        "PATH": "/project/devel/lfrias/gem:/project/devel/lfrias/asm/src:${__file__|parent}:${__file__|parent}/scripts:${__file__|parent}/tools:${PATH}"
    },
    "jip_path" : "${__file__|parent}:${__file__|parent}/tools",
    "extra": ["-C","AD"],
    "account": "AD",
    "queue": "development",
    "jobs": {
        "index_to_be_selected": {
            "threads": 1, 
            "queue": "limited", 
            "time": "3h" 
        }, 
        "selected_index": {
            "threads": 1, 
            "queue": "limited", 
            "time": "3h" 
        }, 
        "selected_biggest_*": {
            "threads": 1, 
            "queue": "limited", 
             "time": "3h" 
        }, 
        "to_index_*": {
            "threads": 1, 
            "queue": "limited", 
            "time": "3h" 
        }, 
        "make_index_*": {
            "threads": 8, 
            "queue": "limited", 
            "time": "3h" 
        }, 
        "to_map_*": {
            "threads": 1, 
            "queue": "limited", 
            "time": "3h" 
        }, 
        "map_and_select_*": {
            "threads": 8, 
            "queue": "limited", 
            "time": "3h"
        }, 
        "map_anchor_chunk_*": {
            "threads": 8, 
            "queue": "limited", 
            "time": "3h"
        }, 
        "collect_*anchor*": {
            "threads": 8
        }, 
        "collect_general_*": {
            "threads": 1 
        }, 
        "fake_parameters_*": {
            "threads": 1, 
            "queue": "limited", 
            "time": "1h"
        }, 
        "filter_anchors_*": {
            "threads": 8, 
            "queue": "limited" 
        }, 
        "alignments_refinement_*": {
            "threads": 1
        }, 
        "assembly_pileup_*": {
            "threads": 8 
        }, 
        "assembly_tips_*": {
            "threads": 1 
        }, 
        "terminal_anchors_*": {
            "threads": 1 
        }, 
        "assembly_graph_*": {
            "threads": 1 
        }, 
        "assembly_solve_*": {
            "threads": 1 
        }, 
        "assembly_repeats_*": {
            "threads": 8 
        }, 
        "assembly_final_*": {
            "threads": 1 
        }, 
        "multiple_align_*": {
            "threads": 1, 
            "queue": "limited"
        },
        "cat_*": {
            "threads": 1, 
            "queue": "limited", 
            "time": "1h" 
        }, 
        "clean_*": {
            "threads": 1, 
            "queue": "limited", 
            "time": "1h" 
        },
        "cleanup": {
            "threads": 1, 
            "queue": "limited",
            "time": "1h" 
        }
    }
}
