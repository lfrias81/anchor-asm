#!/usr/bin/env python

import os
import json
import argparse 
import sys
import math
import subprocess

# TODO: This should be created outside via a template, and cutted with the program

class CreateConfigurationFile(object):
    def __init__(self):
        #To be filled after running the setup
        self.indexChunks = 0
        self.mappingChunks = 0
        self.anchor_list =[]
        self.scaffolds_as_input= False

        #Dictionaries fpr outputting the json
        self.allParameters = {}
        
        self.configParameters = {}
        self.inputParameters = {}     
        self.mainParameters = {}
        self.mappingParameters = {}        
        self.graphParameters = {}
        self.outputParameters = {}
        
        self.envParameters = {}
    
    # Register parameters
    # General function calling all the rest    
    def register_parameters(self, parser):
        self.register_config(parser)
        self.register_input(parser)       
        self.register_main_parameters(parser)
        self.register_mapping_parameters(parser)      
        self.register_graph_parameters(parser)
        self.register_output_parameters(parser)
        
        
    def register_config(self, parser):
        config_group = parser.add_argument_group('configuration file')
        
        config_group.add_argument('--json-file', dest="jsonFile", metavar="<filename>", type=argparse.FileType('w'), required=True, help='Configuration JSON to be generated. Required')
        config_group.add_argument('--debug-mode', dest="debugMode", action='store_true', help="If present, does not remove intermediate files. Intended for advance and debug usage.")

    def required_length(self):
        class RequiredLength(argparse.Action):
            def __call__(self, parser, args, values, option_string=None):
                if not 1<=len(values)<=2:
                    raise argparse.ArgumentTypeError('input-scaffolds parameter requires the sequence file and optionally, the agp file')
                setattr(args, self.dest, values)
        return RequiredLength

    def register_input(self, parser):
        #TODO: make this flexible to be several files
        input_group = parser.add_argument_group('input files') 

        input_group.add_argument('--input-scaffolds', dest="inputScaffolds", nargs='+', action=self.required_length(), metavar=("<contigs_filename>","<scaffolds_agp>"),  help ="Input assembled scaffolds to be merged description: input contigs in fasta or fastq format and optinally, the input agp file must be given")
        input_group.add_argument('--paired-reads', dest="pairedReads", nargs=4, action="append", metavar=("<libSize>", "<anchor_policy: 'front'|'ends'|'all'>", "<file0>", "<file1>"), help ="Input pair end library description: the expected insert size, anchor spacing policy and input fasta/fastq files must be given")
        input_group.add_argument('--mate-paired-reads', dest="matePairedReads", action="append", nargs=4,  metavar=("<libSize>", "<anchor_policy: 'front'|'ends'|'all'>", "<file0>", "<file1>"), help ="Input mate pair library description: the expected insert size, anchor spacing policy and input fasta/fastq files must be given")
        input_group.add_argument('--single-end-reads', dest="singleEndReads", action="append", nargs=2, metavar=("<anchor_policy: 'front'|'ends'|'all'>","<filename>"), help ="Input single end library description:  anchor spacing policy and input fasta/fastq file must be given")
 
        input_group.add_argument('--asm-directory', dest="asmDirectory", metavar="<directory>", required=True, help="Output directory. It should not exist. Required")
         
    def register_main_parameters(self, parser):
        main_group = parser.add_argument_group('main parameters')
        
        main_group.add_argument('--anchor', dest="anchor", metavar="<integer>", type=int, required=True, help="Main anchor size. Required" )
        
        main_group.add_argument('--anchor-spacing', dest="anchorSpacing", metavar="<integer>", type=int, help="Spacing between anchors. Default: anchor size" )
        
        main_group.add_argument('--min-anchor', dest="minAnchor", metavar="<integer>", type=int, help="Min anchor size. Default: anchor size")
        
        main_group.add_argument('--memory', dest="memory", metavar="<integer>", type=int, default=48000000000, help= "Max memory available for use.  Default: %(default)s")
        
        main_group.add_argument('--use-small-reads', dest="useSmallReads", action='store_true', help="Consider smaller reads than the base anchor size. Default: not activated")
       
        main_group.add_argument('--coverage', dest="coverage", metavar="<integer>", type=int, default=1, help="Minimum coverage to not remove a tip or a singleton. Default: %(default)s")
        
    def register_mapping_parameters(self, parser):
        mapping_group = parser.add_argument_group('mapping parameters')
        
        mapping_group.add_argument('--anchorsxchunk', dest="anchorsXChunk", metavar="<integer>", type=int, help="Number of anchors (i.e reads) per alignment job")
        mapping_group.add_argument('--max-size-of-index-chunk', dest="maxSizeOfIndexChunk", metavar="<integer>", type=int, default=10000000000, help="Maximum size of indexing chunk. Default: %(default)s")
        #TODO add some control that divergence is a number between 0 and 1
        mapping_group.add_argument('--divergence', dest="divergence", metavar="<float_over_1>", type=float, required=True, help="Maximum divergence to collapse. Required" )
        mapping_group.add_argument('--min-mapsxanchor', dest="minMapsXAnchor", metavar="<integer>", type=int, default=0, help="Minimum number of allowed mappings per anchor. Default: %(default)s")
        mapping_group.add_argument('--max-mapsxanchor', dest="maxMapsXAnchor", metavar="<integer>", type=int, default=1000, help="Maximum number of allowed mappings per anchor. Default: %(default)s")
        mapping_group.add_argument('--extra-mapping-parameters', dest="extraMappingParameters", metavar="<string>", default="", help="Extra GEM mapping parameters, unparsed. Default: (None)")
        mapping_group.add_argument('--min-index-chunks', dest="minIndexChunks", metavar="<integer>", default=1, help="Chunk indexes to do the read selection, this should be the minimal possible. Default: %(default)s")
 
    def register_output_parameters(self, parser):
        output_group = parser.add_argument_group('output parameters')
        output_group.add_argument('--raw-graph-output', dest="rawGraphOutput", action='store_true', default='store_false', help="Outputs the graph without applyting the non-overlapping transformation. Default:  %(default)s" )
        output_group.add_argument('--output-chunks', dest="outputChunks", metavar="<integer>", default=1, type=int, help= "Number of output chunks to compute the multiple alignment. Default: %(default)s")

        output_group.add_argument('--min-out-scaffold', dest="minOutScaffold", metavar="<integer>", default=200, type=int, help= "Minimum output size for scaffolds. Default: %(default)s")
        # -s parameter   
        output_group.add_argument('--consensus-type', dest="consensusType", metavar="[majority | longest]", default ="majority", choices=['majority', 'longest'], help="Consensus type for output: either majority consensus or the longest contig is chosen. Default: %(default)s")
 
    def register_graph_parameters(self, parser):
        graph_group = parser.add_argument_group('graph parameters')

        # -g parameter
        graph_group.add_argument('--max-distance-diff', dest="maxDistanceDiff", metavar="<float_over_1>", default=0.05, type=float, help= "Maximum allowed distance difference in percentage to allow them to be flagged as equivalent. Default: %(default)s")
        # -r parameter
        graph_group.add_argument('--trim', dest="trim", metavar="<integer>", type=int, help= "Maximum trim without considering coverage. Default: minimum anchor size")
        graph_group.add_argument('--hub-cardinality', dest="hubCardinality", metavar="<integer>", type=int, default=16, help= " Node cardinality to be considered a hub to do aggressive filtering. Value 0 deactivates the filtering. Default: %(default)s")
  
        graph_group.add_argument('--repeat-rounds', dest="repeatRounds", metavar="<integer>", default=1, type=int, help="Number of repeat rounds. Default: %(default)s")       
   
        # -k -l parameter
        graph_group.add_argument('--repeat-resolution-depth', type=int,  metavar="<integer>", dest="repeatResolutionDepth",  default=4, help = "Search depth for repeat resolution. Value 0 deactivates repeat resolution. Default: %(default)s")         
        # -y parameter
        graph_group.add_argument('--path-expansion-depth', type=int, metavar="<integer>",  dest="pathExpansionDepth", default=10, help = "Path expansion depth (scaffold). Value 0 deactivates path expansion. Default: %(default)s")
        # -j parameter
        graph_group.add_argument('--path-selection-ratio', dest="pathSelectionRatio", metavar="<float_over_1>", default=0, type=float, help = "Minimum coverage ratio to select a path, with respect to another. If it is 0, this simplification is not run. Default: %(default)s")
      ######################################################################
    
    def raise_exception_if_not_exists(self, path):
        if not os.path.exists(path):
            raise IOError("Input " + path + " file does not exist")  
    
    def raise_exception_if_not_valid_anchor_policy(self, policy):
        if not (policy =='front' or policy =='ends' or policy == 'all'):
            raise ValueError("Anchor policy " + policy + " is not valid. It must be either 'front', 'ends' or 'all")     
    # Check global conditions on the parameters
    def check_and_set_default_parameters(self,args):
        ## Set defaults
        args.minAnchor = args.minAnchor or args.anchor 
        args.anchorSpacing = args.anchorSpacing or args.anchor 

        self.anchor_list.append(args.anchor)
        next_anchor = args.anchor/2
        while next_anchor >= args.minAnchor:
            self.anchor_list.append(next_anchor)
            next_anchor = next_anchor/2
                  
        #TODO think of some formula
        if not args.anchorsXChunk:
            heuristic = (args.divergence+0.01)*args.anchor
            args.anchorsXChunk = int(100000000/(heuristic*pow(2, math.sqrt(heuristic))))      
        
        
        if not args.trim:
            args.trim = args.minAnchor
        
        ## Check files
        readsActive = args.pairedReads or args.matePairedReads or args.singleEndReads
        if args.inputScaffolds and readsActive:
            raise SyntaxError("Currently, either scaffolds or reads are supported, but not both")        
        if not args.inputScaffolds and not readsActive:
            raise SyntaxError("Either scaffolds or reads have to be specified as input")
        if args.inputScaffolds:
            self.raise_exception_if_not_exists(args.inputScaffolds[0])
            if len(args.inputScaffolds)==2:
                self.raise_exception_if_not_exists(args.inputScaffolds[1])
            self.scaffolds_as_input = True
        if args.singleEndReads:
            for single in args.singleEndReads:
                self.raise_exception_if_not_valid_anchor_policy(single[0])
                self.raise_exception_if_not_exists(single[1])
        if args.matePairedReads:
            for paired in args.matePairedReads:
                self.raise_exception_if_not_valid_anchor_policy(paired[1])
                self.raise_exception_if_not_exists(paired[2])
                self.raise_exception_if_not_exists(paired[3])
        if args.pairedReads:
            for paired in args.pairedReads:                
                self.raise_exception_if_not_valid_anchor_policy(paired[1])
                self.raise_exception_if_not_exists(paired[2])
                self.raise_exception_if_not_exists(paired[3])
        #TODO check the file exists
    
    def get_complete_policy(self, args, policy):
        if policy == 'all':
            return "long-framed(" + str(args.anchor) + "," + str(args.anchorSpacing) + ")"
        elif policy == "front":
            return "short(" + str(args.anchor) + ")"
        elif policy == "ends":
            return "short-framed(" + str(args.anchor) +  ")"
            
    
    def run_raw_setup(self,args):
        # TODO check 
        command = [ "gem-aaa_setup", #os.path.dirname(os.path.realpath(__file__)) + "/bin/gem-aaa_setup",
                   "-p", "bin",
                   "-x", "(" + str(args.maxSizeOfIndexChunk) +",8,,(0,3,0,0))",
                   "-s", str(args.anchorsXChunk),
                   "-m", "(" + str(args.anchorsXChunk) + ",8,,(3,0,0,0))",
                   "-f", "(" + str(args.minMapsXAnchor) + "," + str(args.maxMapsXAnchor) +")",
                   "-o", args.asmDirectory ]
        if args.inputScaffolds:          
            command.append("-i")
            command.append("assembly:assembled(" + args.inputScaffolds[0] + ",ignore,"+ self.get_complete_policy(args, 'all')+ ")")
        if args.singleEndReads:
            for single in args.singleEndReads:
                command.append("-i")
                command.append("illumina:single-ends(" + single[1] + ",ignore,"+ self.get_complete_policy(args, single[0])+ ")")

        if args.matePairedReads:
            for paired in args.matePairedReads:          
                command.append("-i")
                command.append("illumina:mate-pairs(" + paired[0] +
                           ",(" + paired[2] + ",ignore,"+ self.get_complete_policy(args, paired[1])+ ")" + 
                           ",(" + paired[3] + ",ignore,"+ self.get_complete_policy(args, paired[1])+ ")" + 
                           ")")
        
        if args.pairedReads:
            for paired in args.pairedReads:
                command.append("-i")
                command.append("illumina:paired-ends(" + paired[0] +
                             ",(" + paired[2] + ",ignore,"+ self.get_complete_policy(args, paired[1])+ ")" + 
                           ",(" + paired[3] + ",ignore,"+ self.get_complete_policy(args, paired[1])+ ")" + 
                           ")")
        #print command
        if subprocess.call(command):
            raise RuntimeError("Error in the setup.")  
        #Add scaffolds file
        if args.inputScaffolds:  
            if len(args.inputScaffolds)==2:
                command = ["cp", args.inputScaffolds[1], 
                      args.asmDirectory + "/scaffolds.agp"]
            else: #fake file
                command = ["touch", args.asmDirectory + "/scaffolds.agp"]
            if subprocess.call(command):
                raise RuntimeError("Invalid input file.")  
    
            
        with open(args.asmDirectory + '/parameters', 'r') as f:
            data = f.readlines()

            for line in data:
                words = line.split()
                if words[0] == "indexing-chunks":
                    self.indexChunks = words[1]
                    args.minIndexChunks = args.minIndexChunks or  self.indexChunks
                if words[0] == "mapping-chunks":
                    self.mappingChunks = words[1]  
        
        #create subdirectories
        command = ["jip", os.path.dirname(os.path.realpath(__file__)) + "/tools/make_directories.jip",
                   "--asm-dir", args.asmDirectory, 
                   "-a"]
        for a in self.anchor_list:
            command.append(str(a))
                   
        if subprocess.call(command):
            raise RuntimeError("Error in the setup.")  
    
    ####################################################################
    # Store parameters in json 
    # General function calling all the rest 
    def store_parameters(self,arg):
        self.store_config_parameters(arg)
        self.store_input_parameters(arg)
        self.store_main_parameters(arg)    
        self.store_mapping_parameters(arg)   
        self.store_graph_parameters(arg)   
        self.store_output_parameters(arg)  
         
    def store_config_parameters(self,args):
        self.configParameters["jsonFile"] = args.jsonFile.name # We need a string type for the json
        self.configParameters["debugMode"] = args.debugMode 
        self.allParameters["json"] = self.configParameters
        
    def store_input_parameters(self,args):
        self.inputParameters["inputScaffolds"] = args.inputScaffolds
        self.inputParameters["asmDirectory"] = args.asmDirectory
        self.inputParameters["pairedReads"] = args.pairedReads
        self.inputParameters["matePairedReads"] = args.matePairedReads
        self.inputParameters["singleEndReads"] = args.singleEndReads
         
        self.allParameters["input"] = self.inputParameters

    def store_main_parameters(self,args):
        # Build anchor list
        self.mainParameters["scaffoldsAsInput"] = self.scaffolds_as_input
        self.mainParameters["anchors"] = self.anchor_list
        self.mainParameters["anchorSpacing"] = args.anchorSpacing
        self.mainParameters["memory"] = args.memory
        self.mainParameters["coverage"] = args.coverage
        self.mainParameters["useSmallReads"] = args.useSmallReads
        self.allParameters["mainParameters"] = self.mainParameters
                
    def store_mapping_parameters(self,args):
        self.mappingParameters["indexChunks"] = self.indexChunks
        self.mappingParameters["mappingChunks"] = self.mappingChunks
        self.mappingParameters["anchorsXChunk"] = args.anchorsXChunk
        self.mappingParameters["maxSizeOfIndexChunk"] = args.maxSizeOfIndexChunk
        self.mappingParameters["divergence"] = args.divergence
        self.mappingParameters["minMapsXAnchor"] = args.minMapsXAnchor
        self.mappingParameters["maxMapsXAnchor"] = args.maxMapsXAnchor
        self.mappingParameters["extraMappingParameters"] = args.extraMappingParameters       
        self.mappingParameters["minIndexChunks"] = args.minIndexChunks
        self.allParameters["mappingParameters"] = self.mappingParameters
        
    def store_graph_parameters(self,args):
        self.graphParameters["maxDistanceDiff"] = args.maxDistanceDiff
        self.graphParameters["trim"] = args.trim
        self.graphParameters["hubCardinality"] = args.hubCardinality
        self.graphParameters["pathExpansionDepth"] = args.pathExpansionDepth
        self.graphParameters["pathSelectionRatio"] = args.pathSelectionRatio
        self.graphParameters["repeatRounds"] = args.repeatRounds
        self.graphParameters["repeatResolutionDepth"] = args.repeatResolutionDepth
        self.allParameters["graphParameters"] = self.graphParameters
        
    def store_output_parameters(self,args):
        self.outputParameters["rawGraphOutput"] = args.rawGraphOutput
        self.outputParameters["outputChunks"] = args.outputChunks
        self.outputParameters["minOutScaffold"] = args.minOutScaffold
        self.outputParameters["consensusType"] = args.consensusType
        self.allParameters["outputParameters"] = self.outputParameters
        
        
 
##############################################################################################################3
        
# 1. Create object class Configuration File
configManager = CreateConfigurationFile()

# 2. Parsing
# 2.1. Create parser
parser = argparse.ArgumentParser(prog="anchor_assembly_config", description="Create a configuration json file for the assembly pipeline")

# 2.2 Register parameters
configManager.register_parameters(parser)


# 2.3. Do actual parsing, checking for integrity by individual parameters
args = None

try:
    args = parser.parse_args()
except IOError, msg:
    parser.error(str(msg))
except argparse.ArgumentTypeError, msg:
     parser.error(str(msg))

# 3. Check global conditions on the arguments and set defaults
try:
    configManager.check_and_set_default_parameters(args)    
except SyntaxError, msg:
    parser.error(str(msg))
except SyntaxError, msg:
    parser.error(str(msg))
except IOError, msg:
    parser.error(str(msg))
# 4. Run setup, and create files

try:
    configManager.run_raw_setup(args)
except RuntimeError, msg:
    parser.error(str(msg))

# 5. Generate the actual parameters, once created the files Save arguments to super map structure
configManager.store_parameters(args)

# 6. Store JSON file
json.dump(configManager.allParameters, args.jsonFile, indent=2)
    


    
