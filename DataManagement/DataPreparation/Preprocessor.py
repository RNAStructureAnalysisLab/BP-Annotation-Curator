# AUTHOR: Kristopher Church

import json
import csv
import os
import shutil
import re

# TODO add multithreading and a time to finish estimation
# TODO rejecting residue pairs leads to innacuracies in other pipeline steps
# since Table_Extender might treat these as 'nbp' which is 'Not a Base Pair'.
# However, it is not gauranteed that they actually are not base pair
# interactions. (perhaps store them as '?' to represent ambiguity?)
# TODO re-evaluate the contact type dictionaries to see which should be made
# blank or how they can be converted to leontis-westhof style
# Might want to save them as 'REJECT' still but put them into the accepted
# dictionary, and the Table_Extender pipeline step will then remove rows with 
# REJECT from the final table

class Preprocessor:
    INPUT_DIRECTORY = os.path.join(
        'Data', 'Raw', 'AnnotationTools', 'Annotations'
    )
    OUTPUT_DIRECTORY = os.path.join('Data', 'Preprocessed', 'JSON_Annotations')
    initialize_variables = True # False after _on_start() has been invoked
    
    # Following two will be initialized programmatically by _on_start()
    RENAMING_CONVENTION = {}
    DESCRIPTIONS_TO_REJECT = {}
    
    # ACTION: iterates through every CSV file in INPUT_DIRECTORY, and
    # passes them on individually to the helper function, to_json, for conversion
    # into JSON
    @staticmethod
    def convert_all():
        # Programmatically initialize certain class variables
        if Preprocessor.initialize_variables:
            Preprocessor._on_start()
        
        # Ensure OUTPUT_DIRECTORY exists, make empty if so
        if os.path.exists(Preprocessor.OUTPUT_DIRECTORY):
            shutil.rmtree(Preprocessor.OUTPUT_DIRECTORY)
        os.makedirs(Preprocessor.OUTPUT_DIRECTORY, exist_ok=True)
        
        # Begin converting each CSV into a JSON
        for csv_file_name in os.listdir(Preprocessor.INPUT_DIRECTORY):
            if csv_file_name.endswith('.csv'):
                Preprocessor._to_json(csv_file_name)
           
    # INPUT: the name of a CSV file representating a tool annotation
    # ACTION: converts the CSV file into a reformatted JSON file and stores in
    # a new directory called OUTPUT_DIRECTORY
    @staticmethod
    def _to_json(csv_file_name):
        accepted_output_file_name = os.path.splitext(csv_file_name)[0]
        rejected_output_file_name = (
            os.path.splitext(csv_file_name)[0] + '_rejected'
        )
        
        accepted_annotations, rejected_annotations = (
            Preprocessor._load_annotations(csv_file_name)
        )
                
        accepted_annotations, rejected_annotations = Preprocessor._sort_annotations(
            accepted_annotations, rejected_annotations
        )
            
        Preprocessor._export_json(accepted_output_file_name, accepted_annotations)
        Preprocessor._export_json(rejected_output_file_name, rejected_annotations)
                
    # INPUT:s a string representing the name of a CSV annotation file
    # OUTPUT: returns a dictionary representing the base pairing information in
    # the CSV file. Each key is a residue pair tuple, and each value is a set 
    # consisting of the possible base pairing interactions for a given pair.
    # Uses _ensure_format to make the residue pairs be sorted such that the one
    # with the smallest index is on the left, and such that 'description'
    # is strictly base pairing interactions and not also other residue 
    # interactions
    @staticmethod
    def _load_annotations(csv_file_name):
        input_file_path = os.path.join(
            Preprocessor.INPUT_DIRECTORY, csv_file_name
        )
        is_clarna = 'CL' in csv_file_name.split('_')[1]
        accepted_annotations = {}
        rejected_annotations = {}
        with open(input_file_path, 'r') as annotation_csv_file:
            csv_reader = csv.reader(annotation_csv_file)
            next(csv_reader) # skip the first row with only the labels
            for csv_row in csv_reader:
                if is_clarna:
                    residue1, residue2, nucleotides, weight, description, _ = csv_row
                else:
                    residue1, residue2, nucleotides, description, _ = csv_row
                residue_pair = (
                    f'{residue1}{nucleotides[0]}', f'{residue2}{nucleotides[1]}'
                )
                
                new_residue_pair, new_description = Preprocessor._ensure_format(
                    residue_pair, description
                )
                
                # Add to rejected annotations
                if new_description == 'REJECT':
                    if residue_pair in rejected_annotations:
                        rejected_annotations[residue_pair].add(description)
                    else:
                        rejected_annotations[residue_pair] = {description}
                
                # Add to accepted_annotations otherwise
                elif is_clarna: # Adding procedure for ClaRNA
                    if new_description == '':
                        weighted_description = ''
                    else:
                        weighted_description = f'{new_description} {weight}'
                    
                    if new_residue_pair in accepted_annotations and new_description != '':
                        cl_descriptions = list(
                            accepted_annotations[new_residue_pair]
                        )
                        
                        # Retain only instance with the highest weight
                        # Loop assumes there ar at most 2 elements where one is
                        # ''
                        for cl_description in cl_descriptions:
                            if cl_description != '':
                                stored_weight = cl_description.split(' ')[1]
                                if float(weight) > float(stored_weight):
                                    accepted_annotations[new_residue_pair].remove(cl_description)
                                    accepted_annotations[new_residue_pair].add(weighted_description)
                    else:
                        accepted_annotations[new_residue_pair] = {
                            weighted_description
                        }
                else: # Adding procedure for any other tool
                    if new_residue_pair in accepted_annotations:
                        accepted_annotations[new_residue_pair].add(
                            new_description
                        )
                    else:
                        accepted_annotations[new_residue_pair] = {
                            new_description
                        }   
           
        return accepted_annotations, rejected_annotations
    
    @staticmethod
    def _export_json(output_file_name, annotations):
        try:
            pdb_id, tool = output_file_name.split('_')
        except Exception: # output_file_name has '_rejected' in the name
            pdb_id, tool, rejected_tag = output_file_name.split('_')
            tool = f'{tool}_{rejected_tag}'
        output_file_path = os.path.join(
            Preprocessor.OUTPUT_DIRECTORY, tool + '.json'
        )
        
        # load data currently in the JSON file if any
        try:
            with open(output_file_path, 'r') as json_file:
                data = json.load(json_file)
        except:
            data = {}
        for base_pair, descriptions in annotations.items():
            # Change the data type so it is JSON compatible
            base_pair_1 = str(base_pair[0])
            base_pair_2 = str(base_pair[1])
            descriptions = list(descriptions)
            
            if pdb_id not in data:
                data[pdb_id] = {}
            if base_pair_1 not in data[pdb_id]:
                data[pdb_id][base_pair_1] = {}
            data[pdb_id][base_pair_1][base_pair_2] = descriptions
  
        with open(output_file_path, 'w') as json_file:
            json.dump(data, json_file, indent=4)
        
    # INPUT: a tuple representing a residue pair, and a set representing the
    # possible interactions for that residue
    # ACTION: Ensures that the left residue in 'interacting_residues' has the
    # lowest index. Renames the contents of 'descriptions' to only include
    # base pairing contact types in the style of leontis-westhof nomenclature
    # which is the one used by RNA 3D Motif Atlas
    @staticmethod
    def _ensure_format(residue_pair, description):
        was_reversed = False
        converted_residue_pair = ( # Reformat for quantitative comparison
            Preprocessor._convert_residue_pair(residue_pair)
        )
        
        # Rearrange residues to be properly ordered
        if converted_residue_pair[0] > converted_residue_pair[1]:
            was_reversed = True
            first = residue_pair[1]
            second = residue_pair[0]
            residue_pair = (first, second)
            
        # Change description to a proper notation
        description = Preprocessor._rename(description, was_reversed)
        
        return residue_pair, description
        
    # ACTION: Initializes the class variables RENAMING_CONVENTION and
    # DESCRIPTIONS_TO_REJECT. Sets on_start to False so that this is performed
    # only once
    @staticmethod
    def _on_start():
        # Initialize RENAMING_CONVENTION and DESCRIPTIONS_TO_REMOVE
        edge_to_edge_interactions = ['WW', 'WH', 'WS', 'HW', 'HH', 'HS', 'SW', 'SH', 'SS']
        conformations = ['cis', 'tran']
        
        def replace(edge):
            return edge.replace('H', 'M').replace('S', 'm')
        
        for e in edge_to_edge_interactions:
            for c in conformations:
                r3dma_style = f'{c[0]}{e}'
                
                # ClaRNA style contact types to r3dma style ones
                Preprocessor.RENAMING_CONVENTION[f'{e}_{c}'] = r3dma_style
                Preprocessor.RENAMING_CONVENTION[f'?{e}_{c}'] = r3dma_style
                # Other style contact types to r3dma style ones
                Preprocessor.RENAMING_CONVENTION[f'{c[0]}{replace(e[0])}-{replace(e[1])}'] = r3dma_style
                Preprocessor.RENAMING_CONVENTION[f'{c[0]}{replace(e[0])}+{replace(e[1])}'] = r3dma_style
        
        # Create a set of notations that should be removed since they don't cleanly correspond to leontis-westhof nomenclature
        Preprocessor.DESCRIPTIONS_TO_REJECT = {
            '?H_0BPh', '?W_6BPh', '?SW_2BR', '?W_345BPh', '?H_0BR', 'UNK_SHORT_DESC', 'W_6BR', 
            'tW+.', '?diagonal-nc-ww', 'diagonal-nc-ww', '?W_345BR', 'c.-W', 'W_345BPh', 't.-M', 'cW-.', '?H_789BR', 'diagonal-c', 'c.-M', 
            'H_789BPh', 't.-W', 'cW+.', 'H_789BR', 'SW_2BR', '?SW_2BPh', '?W_6BR', '?diagonal-c', 't.+W', 'tW-.', '?S_1BPh', '?S_1BR', 
            '?H_789BPh', 'SW_2BPh', 'tm+.', 'H_0BR', 'W_345BR', 'W_6BPh', 'c.+M', 't.-m', 'tm-.', 't.+m', 'tM-.', 'cM+.', 't.+M', 'c.+W'
        }
        # Types that will be converted to empty strings (stacking, base
        # phosphate, and base ribose interactions):
        # {
        # 'base-ribose-stacking',
        #'?<>', '', '>>', '<<', '<>', '?><', '><', '--', '?<<', '?>>',
        # }
        
        Preprocessor.initialize_variables = False
        
    # INPUT: a string representing a description, and a boolean on whether it
    # should be reversed
    # OUTPUT: a string representing a contact type in the leontis-westhof
    # nomenclature style. If the description was not a base pairing
    # interaction, then it returns an empty string. If it seems like it might
    # be a contact type, but one that can't be converted to leontist-westhof
    # style, then return 'REJECT'
    @staticmethod
    def _rename(description, was_reversed):
        description = description.strip()
        
        if description in Preprocessor.DESCRIPTIONS_TO_REJECT:
            return 'REJECT'
        elif description not in Preprocessor.RENAMING_CONVENTION:
            return ''
        elif was_reversed:
            # TODO Assuming no DSSR styled description
            # Swaps the placement of two edges, IE: ?WH_tran --> ?HW_tran
            delimiter_index = description.find('_')
            description = (
                description[:delimiter_index - 2] + 
                description[delimiter_index - 1] + 
                description[delimiter_index - 2] + 
                description[delimiter_index:]
            )
            
        return Preprocessor.RENAMING_CONVENTION[description]
       
    #INPUT: two dictionaries with a tuple primary key, and a set secondary key
    #OUTPUT: the sorted contents of the original two dictionaries in increasing
    # order, but not strictly lexicographically. IE: D9U is less than D29U
    # because 9 is less than 29. But D9U is less than E5C because D is before E
    @staticmethod
    def _sort_annotations(accepted_annotations, rejected_annotations):
        sorted_accepted_annotations = {
            k: accepted_annotations[k] for k in sorted(
                accepted_annotations.keys(), 
                key=Preprocessor._convert_residue_pair
            )
        }
        sorted_rejected_annotations = {
            k: rejected_annotations[k] for k in sorted(
                rejected_annotations.keys(), 
                key=Preprocessor._convert_residue_pair
            )
        }
        
        return sorted_accepted_annotations, sorted_rejected_annotations
    
    # INPUT: a tuple of two strings representing a residue pair, IE: (D9U, D19U)
    # OUTPUT: a tuple where the strings have been converted into a meaningful
    # format for quantitative comparison. IE: 'D9U < D19U' would be 'True'
    @staticmethod
    def _convert_residue_pair(residue_pair):
        pattern = r'([a-zA-Z]+|\d[a-zA-Z]+|\d)(-?\d+)([^\d]+)'
        residue1, residue2 = residue_pair
        
        residue1_groups = re.match(pattern, residue1)
        residue2_groups = re.match(pattern, residue2)
        
        return (
            residue1_groups.group(1), int(residue1_groups.group(2)),
            residue1_groups.group(3)
        ), (
            residue2_groups.group(1), int(residue2_groups.group(2)),
            residue2_groups.group(3)
        )