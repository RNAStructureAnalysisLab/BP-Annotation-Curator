# AUTHOR: Kristopher Church

import json
import csv
import os
import shutil
import re
import pandas as pd

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
    
    json_cache = {} # key is the output file name, value is the dictionary
    rejects_df = pd.DataFrame(columns=['Tool', 'PDB', 'Residue IDs', 'Original Description'])
    
    # ACTION: iterates through every PDB CSV file in INPUT_DIRECTORY, and
    # passes them on individually to the helper function, to_json, for conversion
    # into JSON compatible dictionary
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
            tool = os.path.splitext(csv_file_name)[0].split('_')[1]
            if csv_file_name.endswith('.csv') and 'MO' != tool:
                Preprocessor._to_json(csv_file_name)
        
        # Save annoation as a JSON file in OUTPUT_DIRECTORY
        for output_file_name, annotation in Preprocessor.json_cache.items():
            Preprocessor._export_json(output_file_name, annotation)
            
        Preprocessor.rejects_df.to_csv(os.path.join('lazy_dir', 'preprocessor_rejects.csv'))
           
    # INPUT: the name of a CSV file representating a tool annotation
    # ACTION: Reformats the data found in the CSV file and stores it in a 
    # JSON compatible dictionary
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
            
        Preprocessor._to_cache(accepted_output_file_name, accepted_annotations)
        Preprocessor._to_cache(rejected_output_file_name, rejected_annotations)
                
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
        pdb, tool = csv_file_name.split('_')
        tool = tool.split('.')[0]
        is_clarna = 'CL' in csv_file_name.split('_')[1]
        is_dssr = 'DSSR' in csv_file_name.split('_')[1]
        accepted_annotations = {}
        rejected_annotations = {}
        with open(input_file_path, 'r') as annotation_csv_file:
            csv_reader = csv.reader(annotation_csv_file)
            next(csv_reader) # skip the first row with only the labels
            for csv_row in csv_reader:
                if is_clarna:
                    residue1, residue2, nucleotides, weight, description, _ = csv_row
                elif is_dssr:
                    residue1, residue2, nucleotides, description = csv_row
                else:
                    residue1, residue2, nucleotides, description, _ = csv_row
                residue_pair = (
                    f'{residue1}{nucleotides[0]}', f'{residue2}{nucleotides[1]}'
                )
                new_residue_pair, new_description = Preprocessor._ensure_format(
                    residue_pair, description
                )
                
                if new_residue_pair not in accepted_annotations.keys():
                    accepted_annotations[new_residue_pair] = set()
                if new_description == 'nbp':
                    continue # create the key, if the set is empty by the end, then fill with "nbp"
                if is_clarna:
                    # we only add to the set if the weight is greater than or
                    # equal to the weights of the descriptions already in the set
                    # if grater, reset the set
                    clarna_descriptions = list(accepted_annotations[new_residue_pair])
                    if clarna_descriptions: # is not empty
                        stored_weight = clarna_descriptions[0].split(' ')[1]
                        if weight < stored_weight:
                            continue
                        if weight > stored_weight:
                            accepted_annotations[new_residue_pair] = set()
                    new_description = f"{new_description} {weight}"
                        
                accepted_annotations[new_residue_pair].add(new_description)
                if new_description == 'REJECT':
                    new_row = {'Tool': tool, 'PDB': pdb, 'Residue IDs': residue_pair, 'Original Description': description}
                    Preprocessor.rejects_df.loc[len(Preprocessor.rejects_df)] = new_row
        #fill empty sets in the dictionary to contain "nbp"
        for key, description_set in accepted_annotations.items():
            if not description_set: # is in fact empty
                accepted_annotations[key].add("nbp")
        return accepted_annotations, rejected_annotations
    
    @staticmethod
    def _export_json(output_file_name, annotations):
        output_file_path = os.path.join(
            Preprocessor.OUTPUT_DIRECTORY, output_file_name + '.json'
        )
  
        with open(output_file_path, 'w') as json_file:
            json.dump(annotations, json_file, indent=4)
        
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
    # TODO: can make this implementation better/more efficient
    @staticmethod
    def _on_start():
        # Initialize RENAMING_CONVENTION and DESCRIPTIONS_TO_REMOVE
        edge_to_edge_interactions = ['WW', 'WH', 'WS', 'HW', 'HH', 'HS', 'SW', 'SH', 'SS']
        conformations = ['cis', 'tran']
        
        def replace(edge):
            return edge.replace('H', 'M').replace('S', 'm')
        
        for e in edge_to_edge_interactions:
            for c in conformations:
                lw_style = f'{c[0]}{e}'
                
                # ClaRNA style contact types to lw style ones
                Preprocessor.RENAMING_CONVENTION[f'{e}_{c}'] = lw_style
                Preprocessor.RENAMING_CONVENTION[f'?{e}_{c}'] = lw_style
                # DSSR style contact types to lw style ones
                Preprocessor.RENAMING_CONVENTION[f'{c[0]}{replace(e[0])}-{replace(e[1])}'] = lw_style
                Preprocessor.RENAMING_CONVENTION[f'{c[0]}{replace(e[0])}+{replace(e[1])}'] = lw_style
                # R3DMA style contact types to lw style ones
                Preprocessor.RENAMING_CONVENTION[lw_style] = lw_style
                Preprocessor.RENAMING_CONVENTION[f'n{lw_style}'] = lw_style
                Preprocessor.RENAMING_CONVENTION[f'{lw_style}a'] = lw_style
                Preprocessor.RENAMING_CONVENTION[f'n{lw_style}a'] = lw_style
                Preprocessor.RENAMING_CONVENTION[f'{c[0]}{e[0].lower()}{e[1]}'] = lw_style
                Preprocessor.RENAMING_CONVENTION[f'n{c[0]}{e[0].lower()}{e[1]}'] = lw_style
                Preprocessor.RENAMING_CONVENTION[f'{c[0]}{e[0].lower()}{e[1]}a'] = lw_style
                Preprocessor.RENAMING_CONVENTION[f'n{c[0]}{e[0].lower()}{e[1]}a'] = lw_style
                Preprocessor.RENAMING_CONVENTION[f'{c[0]}{e[0]}{e[1].lower()}'] = lw_style
                Preprocessor.RENAMING_CONVENTION[f'n{c[0]}{e[0]}{e[1].lower()}'] = lw_style
                Preprocessor.RENAMING_CONVENTION[f'{c[0]}{e[0]}{e[1].lower()}a'] = lw_style
                Preprocessor.RENAMING_CONVENTION[f'n{c[0]}{e[0]}{e[1].lower()}a'] = lw_style
                
                # unknown edge mappings
                Preprocessor.RENAMING_CONVENTION[f'{c[0]}{replace(e[0])}-.'] = f'{c[0]}{e[0]}.'
                Preprocessor.RENAMING_CONVENTION[f'{c[0]}{replace(e[0])}+.'] = f'{c[0]}{e[0]}.'
                Preprocessor.RENAMING_CONVENTION[f'{c[0]}.-{replace(e[0])}'] = f'{c[0]}.{e[0]}'
                Preprocessor.RENAMING_CONVENTION[f'{c[0]}.+{replace(e[0])}'] = f'{c[0]}.{e[0]}'
        
        # Create a set of notations that should be made blank
        Preprocessor.DESCRIPTIONS_TO_REJECT = {
            'UNK_SHORT_DESC', 
            'cBW', 'tBW', 'ncBW', 'ntBW', 'cWB', 'tWB', 'ncWB', 'ntWB'
        }
        '''
        'tW+.',
        't.-M', 'cW-.', 'c.-M', 't.-W', 'cW+.', 
        't.+W', 'tW-.', 'tm+.', 'c.+M', 't.-m', 'tm-.', 
        't.+m', 'tM-.', 'cM+.', 't.+M', 'c.+W',
        '''
        
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
            return 'nbp'
        elif was_reversed:
            if '_' in description:
                # Swaps the placement of two edges, IE: ?WH_tran --> ?HW_tran
                delimiter_index = description.find('_')
                description = (
                    description[:delimiter_index - 2] + 
                    description[delimiter_index - 1] + 
                    description[delimiter_index - 2] + 
                    description[delimiter_index:]
                )
            elif description != '?diagonal-nc-ww' and description != 'diagonal-nc-ww' and description != '?diagonal-c' and description != 'diagonal-c': # Swap placement for a DSSR styled contact type
                description = (
                    description[0] + description[-1] + description[-2] + 
                    description[-3]
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
            
    # INPUT: a string representing a file name, IE: 7A0S_CL, and its annotation
    # ACTION: Further processes output_file_name to remove the part with the 
    # PDB ID. Then stores this as the key in json_cache with its annotation
    @staticmethod
    def _to_cache(output_file_name, annotations):
        try:
            pdb_id, tool = output_file_name.split('_')
            output_file_name = tool
        except ValueError: # output_file_name has '_rejected' in the name
            pdb_id, tool, rejected_tag = output_file_name.split('_')
            output_file_name = f'{tool}_{rejected_tag}'
        
        if output_file_name not in Preprocessor.json_cache:
            Preprocessor.json_cache[output_file_name] = {}
        data = Preprocessor.json_cache[output_file_name]
        
        if pdb_id not in data:
            data[pdb_id] = {}
        for base_pair, descriptions in annotations.items():
            base_1 = str(base_pair[0])
            base_2 = str(base_pair[1])
            descriptions = list(descriptions)

            if base_1 not in data[pdb_id]:
                data[pdb_id][base_1] = {}
            data[pdb_id][base_1][base_2] = descriptions