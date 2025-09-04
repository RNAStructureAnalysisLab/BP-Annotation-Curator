import os
import pandas as pd
from collections import Counter

class Catalogue:
    extended_tables_directory = os.path.join(
        "Data", "Preprocessed", "Extended_Tables"
    )
    extended_table_names = os.listdir(extended_tables_directory)

    @staticmethod
    def run():
        for table_name in Catalogue.extended_table_names:
            cluster_csv = pd.read_csv(os.path.join(
                Catalogue.extended_tables_directory, table_name
            ))
            formalized_cluster = Catalogue._formalize(cluster_csv)
            formalized_cluster['PDB'] = cluster_csv['PDB']
            print(formalized_cluster)
            ranked_cluster = Catalogue._convert_to_rankings(formalized_cluster)
            print(ranked_cluster)
            input("wait")
            
            
    @staticmethod 
    def _formalize(cluster_csv):
        formalized_cluster = pd.DataFrame()
        for column in cluster_csv.columns[::-1]:
            if '-' not in column:
                break
            column_data = cluster_csv[column].str.split(',')
            column_data = Catalogue._standardize(column_data)
            rf, cf = Catalogue._get_frequencies(column_data)
            column_data = Catalogue._condense(column_data)
            column_data = Catalogue._annotate_contact_type(column_data, rf, cf)
            formalized_cluster[column] = column_data
        return formalized_cluster  
    
    @staticmethod
    def _standardize(column_data: pd.Series) -> pd.Series:
        def transform_segment(segment: str) -> str:
            if segment.startswith("nc"):
                return "nc" + segment[2:].upper()
            elif segment.startswith("c"):
                return "c" + segment[1:].upper()
            return segment

        return column_data.apply(lambda lst: [transform_segment(s) for s in lst])
    
    @staticmethod
    def _get_frequencies(column_data: pd.Series):
        rf, cf = {}, {}

        for lst in column_data:
            if not lst:  # skip empty lists
                continue

            # --- Raw frequency count ---
            for elem in lst:
                rf[elem] = rf.get(elem, 0) + 1

            # --- Mode frequency count ---
            # Count frequencies inside the current list
            local_counts = {}
            for elem in lst:
                local_counts[elem] = local_counts.get(elem, 0) + 1

            # Find mode(s) in this list
            max_count = max(local_counts.values())
            modes = [k for k, v in local_counts.items() if v == max_count]

            # Increment cf for each mode
            for mode in modes:
                cf[mode] = cf.get(mode, 0) + 1

        return rf, cf
    
    @staticmethod
    def _condense(column_data: pd.Series) -> pd.Series:
        def condense_list(lst):
            if not lst:
                return []
            counts = Counter(lst)
            # Preserve the order of first occurrence
            seen = []
            condensed = []
            for item in lst:
                if item not in seen:
                    condensed.append(f"{counts[item]}{item}")
                    seen.append(item)
            return condensed

        return column_data.apply(condense_list)
    
    @staticmethod
    def _annotate_contact_type(column_data, rf, cf):
        def annotate_list(lst):
            new_lst = []
            for item in lst:
                # Split the leading number from the contact type
                count_str = ''
                i = 0
                while i < len(item) and item[i].isdigit():
                    count_str += item[i]
                    i += 1
                contact_type = item[i:]  # the string part, e.g., "cWW"
                
                # Look up rf and cf for this contact type (use 0 if missing)
                r_count = rf.get(contact_type, 0)
                c_count = cf.get(contact_type, 0)
                
                # Construct new string
                new_lst.append(f"{count_str}r{r_count}c{c_count}")
            return new_lst

        return column_data.apply(annotate_list)

    @staticmethod
    def _convert_to_rankings(formalized_cluster: pd.DataFrame) -> pd.DataFrame:
        ranked_cluster = pd.DataFrame()
    
        for column in formalized_cluster.columns:
            if column == "PDB":
                ranked_cluster[column] = formalized_cluster[column]
                continue
    
            def rank_list(lst):
                if not lst:
                    return []
    
                # Parse each entry into (count, r_value, c_value, original_item)
                parsed = []
                for item in lst:
                    # Example: "4r68c12"
                    count_part, rest = '', item
                    i = 0
                    while i < len(item) and item[i].isdigit():
                        count_part += item[i]
                        i += 1
                    count = int(count_part)
                    r_index = rest.find('r')
                    c_index = rest.find('c')
                    r_value = int(rest[r_index+1:c_index])
                    c_value = int(rest[c_index+1:])
                    parsed.append((count, r_value, c_value, item))
    
                # --- Dense ranking for r ---
                r_values = sorted({x[1] for x in parsed}, reverse=True)
                r_rank_map = {val: rank+1 for rank, val in enumerate(r_values)}
    
                # --- Dense ranking for c ---
                c_values = sorted({x[2] for x in parsed}, reverse=True)
                c_rank_map = {val: rank+1 for rank, val in enumerate(c_values)}
    
                # Construct new ranked strings
                ranked_list = [f"{count}r{r_rank_map[r]}c{c_rank_map[c]}" for count, r, c, _ in parsed]
                return ranked_list
    
            ranked_cluster[column] = formalized_cluster[column].apply(rank_list)
    
        return ranked_cluster
