import os
import json
import pandas as pd
from collections import Counter
import re
from collections import defaultdict

class Catalogue:
    extended_tables_directory = os.path.join(
        'Data', 'Preprocessed', 'Extended_Tables'    
    )
    
    @staticmethod
    def run():
        tier1_ties, tier2_ties, tier3_ties = {}, {}, {}
        for table_name in os.listdir(Catalogue.extended_tables_directory):
            cluster_df = pd.read_csv(os.path.join(
                Catalogue.extended_tables_directory, table_name    
            ))
            Catalogue._find_ties(cluster_df, tier1_ties, tier2_ties, tier3_ties, table_name)
        print(tier1_ties)
            
    @staticmethod
    def _find_ties(cluster_df, tier1_ties, tier2_ties, tier3_ties, table_name):
        pdb_ids = cluster_df['PDB']
        for column in cluster_df.columns[::-1]:
            if '-' not in column:
                break
            column_data = cluster_df[column].str.split(',')
            column_data = Catalogue._standardize(column_data)
            column_data = Catalogue._condense(column_data)
            Catalogue._add_tier1_ties(column_data)
            #Catalogue._add_tier2_ties(column_data, tier1_ties, column, table_name, pdb_ids)
            #Catalogue._add_tier3_ties(column_data, tier1_ties, column, table_name, pdb_ids)
            
            
    @staticmethod
    def _standardize(column_data: pd.Series) -> pd.Series:
        def transform_segment(segment: str) -> str:
            if segment.startswith("nc"):
                return "c" + segment[2:].upper()
            elif segment.startswith("c"):
                return "c" + segment[1:].upper()
            elif segment.startswith("nt"):
                return "t" + segment[2:].upper()
            elif segment.startswith("t"):
                return "t" + segment[1:].upper()
            return segment
    
        return column_data.apply(lambda lst: [transform_segment(s) for s in lst])

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
    def _add_tier1_ties(column_data):
        for entry in column_data:
            max_count = 0
            most_frequent_contact_types = set()
            for contact_type in entry:
                count, contact_type = Catalogue._split_contact_type(contact_type)
                if most_frequent_contact_types.empty():
                    most_frequent_contact_types.add(contact_type)
                    max_count = count
                elif count == max_count:
                    most_frequent_contact_types.add(contact_type)
                elif count > max_count:
                    max_count = count
                    most_frequent_contact_types = set()
                    most_frequent_contact_types.add(contact_type)
            consensus = Catalogue._resolve_tie(most_frequent_contact_types)
                    
    @staticmethod
    def _split_contact_type(string):
        number, contact_type = "", ""
        for character in string:
            if character.isdigit():
                number += character
            else:
                contact_type += character
        return (int(number), contact_type)
    
    @staticmethod
    def _resolve_tie(most_frequent_contact_types):
        return None
    '''
    @staticmethod
    def _add_tier1_ties(column_data, tier1_ties, column, table_name, pdb_ids):
        """
        Detects ties in contact types and inserts them into tier1_ties.
        
        Args:
            column_data (pd.Series): each row is a list of contact types (like ['5cWW', '1nbp'])
            tier1_ties (dict): dictionary to update in place
            column (str): the column name being processed
            table_name (str): the table name (cluster_table_name)
            pdb_ids (pd.Series): corresponding pdb ids for each row
        """
        
        for i, entry in enumerate(column_data):
            if not entry:
                continue
            
            # Parse (coeff, contact_type)
            parsed = []
            for ct in entry:
                j = 0
                while j < len(ct) and ct[j].isdigit():
                    j += 1
                coeff = int(ct[:j]) if j > 0 else 0
                parsed.append((coeff, ct))
            
            max_coeff = max(coeff for coeff, _ in parsed)
            winners = [ct for coeff, ct in parsed if coeff == max_coeff]
            
            if len(winners) > 1:  # Tie found
                # Ensure nested dict exists
                if table_name not in tier1_ties:
                    tier1_ties[table_name] = {}
                if column not in tier1_ties[table_name]:
                    tier1_ties[table_name][column] = []
                
                # Append the pdb id associated with this tie
                tier1_ties[table_name][column].append(pdb_ids.iloc[i])
    '''
                
    @staticmethod
    def _add_tier2_ties(column_data, tier1_ties, column, table_name, pdb_ids):
        pass
    
    @staticmethod
    def _add_tier3_ties(column_data, tier3_ties, column, table_name, pdb_ids):
        def parse_contact(contact):
            """Split contact into (coefficient, suffix, caps)."""
            match = re.match(r"(\d+)([A-Za-z]+)", contact)
            if not match:
                return 0, contact, set()
            coeff = int(match.group(1))
            suffix = match.group(2)
            caps = {c for c in suffix if c.isupper()}
            return coeff, suffix, caps

        for entry, pdb in zip(column_data, pdb_ids):
            if not entry:
                continue

            # --- Group contacts by overlapping capitalized letters ---
            groups = []
            for contact in entry:
                coeff, suffix, caps = parse_contact(contact)
                placed = False
                for group in groups:
                    if not group["caps"].isdisjoint(caps):  # overlap → same group
                        group["coeff"] += coeff
                        group["suffixes"].append(suffix)
                        group["caps"].update(caps)
                        placed = True
                        break
                if not placed:
                    groups.append({"coeff": coeff, "suffixes": [suffix], "caps": caps})

            # --- Build combined contact types ---
            combined_contacts = [
                (g["coeff"], "".join(g["suffixes"])) for g in groups
            ]

            # --- Find max coefficient ---
            max_coeff = max(c[0] for c in combined_contacts)
            top_contacts = [c for c in combined_contacts if c[0] == max_coeff]

            # --- If tie, store result ---
            if len(top_contacts) > 1:
                tier3_ties.setdefault(table_name, {}).setdefault(column, []).append(pdb)

'''
class Catalogue:
    extended_tables_directory = os.path.join(
        "Data", "Preprocessed", "Extended_Tables"
    )
    extended_table_names = os.listdir(extended_tables_directory)

    @staticmethod
    def run():
        global_lookup = {}
        for table_name in Catalogue.extended_table_names:
            cluster_csv = pd.read_csv(os.path.join(
                Catalogue.extended_tables_directory, table_name
            ))
            formalized_cluster = Catalogue._formalize(cluster_csv)
            formalized_cluster['PDB'] = cluster_csv['PDB']
            ranked_cluster = Catalogue._convert_to_rankings(formalized_cluster)
            Catalogue._build_dictionary(ranked_cluster, table_name, global_lookup)
        ties = Catalogue._filter_ties(global_lookup)
        #print(ties)
        #input("wait")
        
        #let's look up instances of that one example I thought of:
        #print(ties.keys())
        #input("okay")
        
        print(f"""
              {ties[('2r1c1', '2r2c1', '2r2c1')]}
              {ties[('2r1c1', '2r1c1', '2r1c1')]}
              """)
        input("alright")
            
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

    @staticmethod
    def _build_dictionary(ranked_cluster: pd.DataFrame, table_name: str, lookup: dict):
        for col in ranked_cluster.columns:
            if col == "PDB":
                continue
    
            for idx, row in ranked_cluster.iterrows():
                key = tuple(sorted(row[col]))  # tuple of ranked list
                pdb_id = row["PDB"]
    
                value = [pdb_id, col, table_name]
    
                if key not in lookup:
                    lookup[key] = []
                lookup[key].append(value)
    
        return lookup
    
    @staticmethod
    def _filter_ties(global_lookup: dict) -> dict:
        tie_dict = {}
    
        for key, value in global_lookup.items():
            # Extract the leading coefficient from each item in the tuple
            counts = []
            for item in key:
                num_str = ""
                i = 0
                while i < len(item) and item[i].isdigit():
                    num_str += item[i]
                    i += 1
                counts.append(int(num_str))
    
            # Find the maximum coefficient
            max_count = max(counts)
    
            # Count how many items have that maximum
            if counts.count(max_count) >= 2:
                tie_dict[key] = value  # keep only ties
    
        return tie_dict
'''