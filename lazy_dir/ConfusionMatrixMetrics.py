import os
import pandas as pd
from collections import Counter

class ConfusionMatrixMetrics():
    CONFUSION_MATRIX_DIRECTORY = os.path.join(
        '.', 'lazy_data'
    )
    BENCHMARK_PATH = os.path.join(
        '..', 'Data', 'ResultsAnalysis', 'BasePairingAgreementsSummary', 'base_pairing_agreements.csv'
    )
    OUTPUT_DIRECTORY = os.path.join(
        '.', 'metrics'
    )
    
    @staticmethod
    def run() -> None:
        metrics = {'CL': None, 'FR': None, 'MC': None, 'RV': None, 'DSSR': None}
        performance = metrics.copy()
        for confusion_matrix in os.listdir(ConfusionMatrixMetrics.CONFUSION_MATRIX_DIRECTORY):
            tool = confusion_matrix.split('_')[0]
            path = os.path.join(ConfusionMatrixMetrics.CONFUSION_MATRIX_DIRECTORY, confusion_matrix)
            confusion_matrix_df = (
                pd.read_csv(path, index_col=0).apply(pd.to_numeric, errors='raise')
            )
            confusion_matrix_df = confusion_matrix_df[confusion_matrix_df.index != "?"].drop(columns=["?"])
            metrics_df = ConfusionMatrixMetrics._create_metrics_table()
            ConfusionMatrixMetrics._calculate_metrics(metrics_df, confusion_matrix_df)
            performance_df = ConfusionMatrixMetrics._calculate_performance(metrics_df)
            metrics[tool] = metrics_df
            performance[tool] = performance_df
            
        benchmark_df = pd.read_csv(ConfusionMatrixMetrics.BENCHMARK_PATH)
        ambiguity_metrics_df = ConfusionMatrixMetrics._calculate_ambiguity_metrics(benchmark_df)
        
        ConfusionMatrixMetrics._write_tables(metrics, performance, ambiguity_metrics_df)
    
    @staticmethod
    def _create_metrics_table() -> pd.DataFrame:
        labels = [
            "cWW", "cWH", "cWS", "cHW", "cHH", "cHS", "cSW", "cSH", "cSS",
            "tWW", "tWH", "tWS", "tHW", "tHH", "tHS", "tSW", "tSH", "tSS", "nbp"
        ]
        metrics = ["TP", "FP", "FN", "TN"]
        
        return pd.DataFrame(0.0, index=labels, columns=metrics)
    
    @staticmethod
    def _calculate_metrics(metrics_df: pd.DataFrame, confusion_matrix_df: pd.DataFrame) -> None:
        labels = confusion_matrix_df.index
        TP = confusion_matrix_df.values.diagonal()
        FP = confusion_matrix_df.sum(axis=0).values - TP
        FN = confusion_matrix_df.sum(axis=1).values - TP
        TN = confusion_matrix_df.values.sum() - (FP + FN + TP)
        results = pd.DataFrame({
            'TP': TP, 'FP': FP, 'FN': FN, 'TN': TN
        }, index=labels)
        
        for row_index in metrics_df.index:
            #total = results.loc[row_index].sum()
            metrics_df.loc[row_index, 'TP'] = results.loc[row_index, 'TP'] #/total
            metrics_df.loc[row_index, 'FP'] = results.loc[row_index, 'FP'] #/total
            metrics_df.loc[row_index, 'FN'] = results.loc[row_index, 'FN']#/total
            metrics_df.loc[row_index, 'TN'] = results.loc[row_index, 'TN']#/total
            
    @staticmethod
    def _calculate_ambiguity_metrics(benchmark_df: pd.DataFrame) -> pd.DataFrame:
        df_rows = {}
        five_way_tie_count = 0
        two_way_tie_count = 0
        for all_annotations, consensus in zip(benchmark_df['All Annotations'], benchmark_df['Expected Contact Type']):
            if consensus == "unresolved tie":
                continue
            all_annotations = all_annotations.split(',')
            counts = Counter(all_annotations) 
            mode, max_frequency = counts.most_common(1)[0]
            top_contact_types = [contact_type for contact_type, frequency in counts.items() if frequency == max_frequency]
            
            if len(top_contact_types) == 2:
                ambiguity_type = ConfusionMatrixMetrics._get_type(counts)
                current_counts = df_rows.get(ambiguity_type, {'Frequency': 0, 'CL': [0,0], 'FR': [0,0], 'MC': [0,0], 'RV': [0,0], 'DSSR': [0,0]})
                ConfusionMatrixMetrics._update_counts(current_counts, all_annotations, consensus)
                df_rows[ambiguity_type] = current_counts
                two_way_tie_count += 1
            elif len(top_contact_types) > 2:
                #TODO input(f"TODO handle this later {top_contact_types} {all_annotations}")
                five_way_tie_count += 1
                continue
        print(five_way_tie_count)
        print(two_way_tie_count)
        df = pd.DataFrame.from_dict(df_rows, orient='index')
        tool_columns = ["CL", "FR", "MC", "RV", "DSSR"]
        df[tool_columns] = df[tool_columns].apply(lambda col: col.apply(lambda x: ((x[0] / x[1]) * 100) if x[1] != 0 else 0))
        return df
            
    @staticmethod 
    def _get_type(counts: Counter) -> str:
        type_notation = ""
        if len(counts) != 3:
            raise ValueError(f"Expected exactly 3 distinct contact types, got {len(counts)}: {counts}")
        tie_group_1, tie_group_2, outlier = [item for item, _ in counts.most_common(3)] 
        
        type_notation += ConfusionMatrixMetrics._compare(tie_group_1, tie_group_2) + ":"
        first_term = ConfusionMatrixMetrics._compare(tie_group_2, outlier)
        second_term = ConfusionMatrixMetrics._compare(tie_group_1, outlier)
        
        if ("_" == first_term and "_" == second_term) or "_" == second_term:
            type_notation += first_term + "," + second_term
        elif "_" == first_term:
            type_notation += second_term + "," + first_term
        else:
            first_term, second_term = sorted([first_term, second_term], reverse=True)
            type_notation += first_term + "," + second_term
        return type_notation
            
    @staticmethod
    def _compare(group_1: str, group_2: str) -> str:
        if "nbp" in [group_1, group_2]:
            return "_"
        shared_edges = ConfusionMatrixMetrics._get_shared_edge_count(group_1, group_2)
        orientation = "T" if group_1[0] == group_2[0] else "F"
        if shared_edges == 0:
            position = "-"
        elif group_1[1] == group_2[1] or group_1[2] == group_2[2]:
            position = "T"
        else:
            position = "F"
        return f"{shared_edges}{orientation}{position}"
    
    @staticmethod
    def _get_shared_edge_count(group1: str, group2: str) -> int:
        edges1 = [group1[1], group1[2]]
        edges2 = [group2[1], group2[2]]
        return sum(min(edges1.count(edge), edges2.count(edge)) for edge in "WHS")
    
    @staticmethod
    def _update_counts(current_counts: dict[str, list[int]], all_annotations: list[str], consensus: str) -> None:
        index_to_tool= {0: "CL", 1: "FR", 2: "MC", 3: "RV", 4: "DSSR"}
        current_counts['Frequency'] += 1
        for i in range(len(all_annotations)):
            contact_type = all_annotations[i]
            tool = index_to_tool[i]
            if consensus == contact_type:
                current_counts[tool][0] += 1
                current_counts[tool][1] += 1
            else:
                current_counts[tool][1] += 1
    
    @staticmethod
    def _write_tables(metrics: dict[str, pd.DataFrame], performance: dict[str, pd.DataFrame], ambiguity_metrics_df: pd.DataFrame) -> None:
        os.makedirs(ConfusionMatrixMetrics.OUTPUT_DIRECTORY, exist_ok=True)
        for tool, metrics_df in metrics.items():
            path = os.path.join(ConfusionMatrixMetrics.OUTPUT_DIRECTORY, f"{tool}_confusion_matrix_metrics.csv")
            metrics_df.to_csv(path)
        for tool, performance_df in performance.items():
            path = os.path.join(ConfusionMatrixMetrics.OUTPUT_DIRECTORY, f"{tool}_performance.csv")
            performance_df.to_csv(path)
        path = os.path.join(ConfusionMatrixMetrics.OUTPUT_DIRECTORY, "ambiguity_type_metrics.csv")
        ambiguity_metrics_df.to_csv(path)
        
    @staticmethod
    def _calculate_performance(metrics_df: pd.DataFrame) -> pd.DataFrame:
        tp = metrics_df['TP'].astype(float)
        fp = metrics_df['FP'].astype(float)
        fn = metrics_df['FN'].astype(float)
        precision = tp / (tp + fp)
        recall    = tp / (tp + fn)
    
        f1 = 2 * precision * recall / (precision + recall)
    
        performance_df = pd.DataFrame(
            {'F1': f1, 'P-R': precision - recall},
            index=metrics_df.index
        ).fillna(0.0)
    
        return performance_df
    
ConfusionMatrixMetrics.run()