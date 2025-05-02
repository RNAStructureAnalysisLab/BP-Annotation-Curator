# AUTHOR: Kristopher Church

import os
import shutil
import statistics
import pandas as pd

# TODO organize the code later when have time
class Agreement_Analyzer:
    RESULTS_DIRECTORY = os.path.join('Data', 'ResultsAnalysis')
    OUTPUT_DIRECTORY = os.path.join(
        'Data', 'ResultsAnalysis', 'BasePairingAgreementsSummary'
    )
    AGREEMENT_DF = pd.read_csv(
        os.path.join(
            RESULTS_DIRECTORY, 'base_pairing_agreements.csv'
        )
    )
    TOOLS = ['R3DMA', 'CL', 'FR', 'MC', 'RV']
    
    bp_interaction_counts = [[0,0], [0,0], [0,0], [0,0], [0,0]]
    nbp_interaction_counts = [[0,0], [0,0], [0,0], [0,0], [0,0]]
    total_bp_interactions = 0
    total_nbp_interactions = 0
    
    @staticmethod
    def run():
        # Reset the contents of OUTPUT_DIRECTORY folder to be blank
        if os.path.exists(Agreement_Analyzer.OUTPUT_DIRECTORY):
            shutil.rmtree(Agreement_Analyzer.OUTPUT_DIRECTORY)
        os.makedirs(Agreement_Analyzer.OUTPUT_DIRECTORY)
        
        summary = []
        
        for _, row in Agreement_Analyzer.AGREEMENT_DF.iterrows():
            if 'nbp' == row['Consensus Contact Type']:
                for i in range(len(Agreement_Analyzer.TOOLS)):
                    Agreement_Analyzer.nbp_interaction_counts[i][0] += (
                        int(row.iloc[6 + i])
                    )
                    Agreement_Analyzer.nbp_interaction_counts[i][1] += (
                        int(row.iloc[12 + i])
                    )
                Agreement_Analyzer.total_nbp_interactions += 1
            else:
                for i in range(len(Agreement_Analyzer.TOOLS)):
                    Agreement_Analyzer.bp_interaction_counts[i][0] += (
                        int(row.iloc[6 + i])    
                    )
                    Agreement_Analyzer.bp_interaction_counts[i][1] += (
                        int(row.iloc[12 + i])    
                    )
                Agreement_Analyzer.total_bp_interactions += 1
        
        Agreement_Analyzer._normalize()
        summary.append({
            'Interaction Type': 'Base Pairing',
            'Agreement Interaction (R3DMA)': Agreement_Analyzer.bp_interaction_counts[0][0],
            'Agreement Match (R3DMA)': Agreement_Analyzer.bp_interaction_counts[0][1],
            'Agreement Interaction (CL)': Agreement_Analyzer.bp_interaction_counts[1][0],
            'Agreement Match (CL)': Agreement_Analyzer.bp_interaction_counts[1][1],
            'Agreement Interaction (FR)': Agreement_Analyzer.bp_interaction_counts[2][0],
            'Agreement Match (FR)': Agreement_Analyzer.bp_interaction_counts[2][1],
            'Agreement Interaction (MC)': Agreement_Analyzer.bp_interaction_counts[3][0],
            'Agreement Match (MC)': Agreement_Analyzer.bp_interaction_counts[3][1],
            'Agreement Interaction (RV)': Agreement_Analyzer.bp_interaction_counts[4][0],
            'Agreement Matche (RV)': Agreement_Analyzer.bp_interaction_counts[4][1]
        })
        summary.append({
            'Interaction Type': 'Not Base Pairing',
            'Agreement Interaction (R3DMA)': Agreement_Analyzer.nbp_interaction_counts[0][0],
            'Agreement Match (R3DMA)': Agreement_Analyzer.nbp_interaction_counts[0][1],
            'Agreement Interaction (CL)': Agreement_Analyzer.nbp_interaction_counts[1][0],
            'Agreement Match (CL)': Agreement_Analyzer.nbp_interaction_counts[1][1],
            'Agreement Interaction (FR)': Agreement_Analyzer.nbp_interaction_counts[2][0],
            'Agreement Match (FR)': Agreement_Analyzer.nbp_interaction_counts[2][1],
            'Agreement Interaction (MC)': Agreement_Analyzer.nbp_interaction_counts[3][0],
            'Agreement Match (MC)': Agreement_Analyzer.nbp_interaction_counts[3][1],
            'Agreement Interaction (RV)': Agreement_Analyzer.nbp_interaction_counts[4][0],
            'Agreement Matche (RV)': Agreement_Analyzer.nbp_interaction_counts[4][1]
        })
        df = pd.DataFrame(summary)
        output_path = os.path.join(Agreement_Analyzer.OUTPUT_DIRECTORY, 'tool_interaction_agreement_with_concensus.csv')
        df.to_csv(output_path, index=False)
        
        # Get deviation across rows
        row_deviation = Agreement_Analyzer._get_row_differences(
            Agreement_Analyzer.bp_interaction_counts
        )
        print(row_deviation)
        row_deviation = Agreement_Analyzer._get_deviations(row_deviation)
        summary = []
        summary.append({
            'Deviation (R3DMA)': row_deviation[0],
            'Deviation (CL)': row_deviation[1],
            'Deviation (FR)': row_deviation[2],
            'Deviation (MC)': row_deviation[3],
            'Deviation (RV)': row_deviation[4]
        })
        df = pd.DataFrame(summary)
        output_path = os.path.join(Agreement_Analyzer.OUTPUT_DIRECTORY, 'tool_contact_type_disagreement_with_concensus.csv')
        df.to_csv(output_path, index=False)
        
        bp_agreement = Agreement_Analyzer._get_agreement()
        bp_agreement_deviation = Agreement_Analyzer._get_deviations(bp_agreement)
        summary = []
        summary.append({
            'Deviation (R3DMA)': bp_agreement_deviation[0],
            'Deviation (CL)': bp_agreement_deviation[1],
            'Deviation (FR)': bp_agreement_deviation[2],
            'Deviation (MC)': bp_agreement_deviation[3],
            'Deviation (RV)': bp_agreement_deviation[4]   
        })
        df = pd.DataFrame(summary)
        output_path = os.path.join(Agreement_Analyzer.OUTPUT_DIRECTORY, 'tool_base_pairing_reporting.csv')
        df.to_csv(output_path, index=False)
        
        '''
        # Now get differences across columns (bp - nbp)
        column_deviation = Agreement_Analyzer._get_column_differences()
        column_deviation = Agreement_Analyzer._get_deviations(column_deviation)
        summary = []
        summary.append({
            'Deviation (R3DMA)': column_deviation[0],
            'Deviation (CL)': column_deviation[1],
            'Deviation (FR)': column_deviation[2],
            'Deviation (MC)': column_deviation[3],
            'Deviation (RV)': column_deviation[4]
        })
        df = pd.DataFrame(summary)
        output_path = os.path.join(Agreement_Analyzer.OUTPUT_DIRECTORY, 'tool_')
        '''
    @staticmethod
    def _normalize():
        for i in range(len(Agreement_Analyzer.TOOLS)):
            Agreement_Analyzer.bp_interaction_counts[i][0] /= Agreement_Analyzer.total_bp_interactions
            Agreement_Analyzer.bp_interaction_counts[i][1] /= Agreement_Analyzer.total_bp_interactions
            Agreement_Analyzer.nbp_interaction_counts[i][0] /= Agreement_Analyzer.total_nbp_interactions
            Agreement_Analyzer.nbp_interaction_counts[i][1] /= Agreement_Analyzer.total_nbp_interactions
         
    # ACTION: Gets the difference within each pair of the list. The difference 
    # shows the percentage of rows where the tool agreed with the consensus on 
    # whether there is a bp interaction or not, but still disagreed on what 
    # that interaction was.
    # OUTPUT: A list of floats
    @staticmethod
    def _get_row_differences(counts_list):
        differences = []
        for pair in counts_list:
            differences.append(pair[0] - pair[1])
        return differences
    
    # ACTION: Normalizes the differences with the average to where we can see
    # whether it is disagreeing more often or less often than the other tools
    @staticmethod
    def _get_deviations(bp_differences):
        mean = statistics.mean(bp_differences)
        return [x - mean for x in bp_differences]
    
    # Trying to figure out whether tool has a tendency to over report or under report base pairing interactions
    @staticmethod
    def _get_column_differences():
        differences = []
        for i in range(len(Agreement_Analyzer.TOOLS)):
            bp_agreement = Agreement_Analyzer.bp_interaction_counts[i][0]
            nbp_agreement = Agreement_Analyzer.nbp_interaction_counts[i][0]
            differences.append(bp_agreement - nbp_agreement)
        return differences
    
    @staticmethod
    def _get_agreement():
        bp_agreement = []
        for pair in Agreement_Analyzer.bp_interaction_counts:
            bp_agreement.append(pair[0])
        return bp_agreement