#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 18:11:17 2024

@author: azlannisar
"""

def create_boxplots_with_data(dataframes_list, x_values, contig_positions_dict, enzyme_count_df, log_results_df):
    unique_motifs = dataframes_list[0]['Motif'].unique()

    for motif in unique_motifs:
        if '"A"' not in motif:
            continue

        scores = []
        counts = []
        for df in dataframes_list:
            df_motif = df[df['Motif'] == motif]
            score_values = df_motif['Score1'].values
            scores.append(score_values)
            counts.append(len(score_values))

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.boxplot(scores, labels=x_values, showfliers=False)
        ax.set_xlabel('Isolates')
        ax.set_ylabel('Score')
        ax.set_ylim(0, 110)
        ax.set_title(f'6mA Boxplot for {motif}')
        plt.xticks(rotation=45, ha='right')
        ax.axhline(y=25, color='magenta', linestyle='--', linewidth=1)

        for i, count in enumerate(counts):
            ax.text(i + 1, 10, f'n={count}', ha='center', va='top', fontsize=10, color='black', fontweight='bold')

        collected_data = []
        for isolate in x_values:
            df_isolate = dataframes_list[x_values.index(isolate)]
            motif_data = df_isolate[df_isolate['Motif'] == motif]

            for mod_pos in motif_data['Mod_Pos']:
                if isolate in contig_positions_dict and mod_pos in contig_positions_dict[isolate]:
                    entries = contig_positions_dict[isolate][mod_pos]
                    for entry in entries:
                        if entry['MethylationType'] == '6mA':
                            entry_with_isolate = entry.copy()
                            entry_with_isolate['Isolate'] = isolate
                            collected_data.append(entry_with_isolate)

        if collected_data:
            unique_data = pd.DataFrame(collected_data).drop_duplicates().reset_index(drop=True)

            position_column_index = unique_data.columns.get_loc('Position')
            columns_to_process = unique_data.columns[position_column_index + 1:]

            def update_collected_row(row):
                enzyme = row['Enzyme']
                for identified_column in columns_to_process:
                    if pd.notna(row[identified_column]):
                        filtered_enzyme_count = enzyme_count_df[enzyme_count_df['Isolate'] == identified_column]
                        filtered_row = filtered_enzyme_count[
                            (filtered_enzyme_count['Enzyme'] == enzyme) &
                            (filtered_enzyme_count['Motif'] == motif)
                        ]

                        if len(filtered_row) > 0:
                            count_value = filtered_row.iloc[0]['Count']
                            
                            # Lookup beta coefficient from Log_Results_df
                            beta_row = log_results_df[
                                (log_results_df['Isolate'] == identified_column) &
                                (log_results_df['Enzyme'] == enzyme)
                            ]
                            
                            if len(beta_row) > 0:
                                beta_value = beta_row.iloc[0]['beta_coefficient']
                                row[identified_column] = f"{row[identified_column]:.2e} ({count_value})\nß={beta_value:.2f}"
                            else:
                                row[identified_column] = f"{row[identified_column]:.2e} ({count_value})"

                return row

            unique_data = unique_data.apply(update_collected_row, axis=1)
            unique_data = unique_data.loc[:, ~unique_data.columns.isin(['Isolate'])]
            unique_data = unique_data.drop_duplicates().reset_index(drop=True)

            table_data = unique_data.values
            col_labels = unique_data.columns
            cell_text = []

            for row in table_data:
                formatted_row = ['{:.2e}'.format(value) if isinstance(value, (float, np.float64)) else str(value) for value in row]
                cell_text.append(formatted_row)

            table = plt.table(cellText=cell_text, colLabels=col_labels, cellLoc='center', loc='bottom', bbox=[-0.1, -1.7, 1.1, 1.1])
            table.auto_set_font_size(False)
            table.set_fontsize(8)
            table.scale(1, 1.5)
            plt.subplots_adjust(left=0.2, bottom=0.6)
        else:
            plt.subplots_adjust(left=0.2, bottom=0.2)
        plt.show()