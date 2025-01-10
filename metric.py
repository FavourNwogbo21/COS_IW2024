import pandas as pd
from scipy.stats import spearmanr

def calculate_metrics(benchling_csv, guidescan_csv):
    benchling_df = pd.read_csv(benchling_csv)
    guidescan_df = pd.read_csv(guidescan_csv)
    benchling_df['Specificity Score'] = pd.to_numeric(benchling_df['Specificity Score'], errors='coerce')
    benchling_df['Efficiency Score'] = pd.to_numeric(benchling_df['Efficiency Score'], errors='coerce')
    guidescan_df['Specificity'] = pd.to_numeric(guidescan_df['Specificity'], errors='coerce')
    guidescan_df['Cutting Efficiency'] = pd.to_numeric(guidescan_df['Cutting Efficiency'], errors='coerce')
    benchling_df = benchling_df.dropna(subset=['Specificity Score', 'Efficiency Score'])
    guidescan_df = guidescan_df.dropna(subset=['Specificity', 'Cutting Efficiency'])
    guidescan_df['Specificity'] *= 100
    guidescan_df['Cutting Efficiency'] *= 100

    # Align dataframes by matching guides based on 'Sequence' (Benchling) and 'gRNA-Seq' (GuideScan)
    common_guides = set(benchling_df['Sequence']).intersection(set(guidescan_df['gRNA-Seq']))

    # filtering both tools' dataframes to include only the guides in common
    guidescan_df = guidescan_df[guidescan_df['gRNA-Seq'].isin(common_guides)].reset_index(drop=True)
    benchling_df = benchling_df[benchling_df['Sequence'].isin(common_guides)].reset_index(drop=True)

    # metric formulas
    def compute_metric(S, E, weights):
        return weights[0] * S + weights[1] * E

    metrics = {
        'M1': (0.5, 0.5),
        'M2': (0.75, 0.25),
        'M3': (0.25, 0.75)
    }

    # metrics for GuideScan
    for metric, weights in metrics.items():
        guidescan_df[metric] = compute_metric(guidescan_df['Specificity'], guidescan_df['Cutting Efficiency'], weights)

    # metrics for Benchling
    for metric, weights in metrics.items():
        benchling_df[metric] = compute_metric(benchling_df['Specificity Score'], benchling_df['Efficiency Score'], weights)

    # ranking guides within each metric
    for metric in metrics.keys():
        guidescan_df[f'{metric}_Rank'] = guidescan_df[metric].rank(ascending=False)
        benchling_df[f'{metric}_Rank'] = benchling_df[metric].rank(ascending=False)

    # selecting the same 5 sequences for display
    display_sequences = guidescan_df['gRNA-Seq'].head(5).tolist()
    guidescan_display = guidescan_df[guidescan_df['gRNA-Seq'].isin(display_sequences)]
    benchling_display = benchling_df[benchling_df['Sequence'].isin(display_sequences)]

    # calculating Spearman correlations
    correlations = {}
    p_values = {}
    for metric in metrics.keys():
        x_corr, x_pval = spearmanr(guidescan_df[f'{metric}_Rank'], guidescan_df.index)
        y_corr, y_pval = spearmanr(guidescan_df[f'{metric}_Rank'], benchling_df[f'{metric}_Rank'])
        correlations[metric] = {
            'x': x_corr,
            'y': y_corr,
            'z': (x_corr + y_corr) / 2
        }
        p_values[metric] = {
            'x_pval': x_pval,
            'y_pval': y_pval
        }

    # finding the best performing metric
    best_metric = max(correlations.items(), key=lambda item: item[1]['z'])[0]

    # ranking guides by the best performing metric
    guidescan_df['Best_Metric_Rank'] = guidescan_df[best_metric].rank(ascending=False)
    benchling_df['Best_Metric_Rank'] = benchling_df[best_metric].rank(ascending=False)

    # calculating statistics for the best metric
    guidescan_stats = {
        'mean': guidescan_df[best_metric].mean(),
        'median': guidescan_df[best_metric].median(),
        'std': guidescan_df[best_metric].std(),
        'min': guidescan_df[best_metric].min(),
        'max': guidescan_df[best_metric].max()
    }

    benchling_stats = {
        'mean': benchling_df[best_metric].mean(),
        'median': benchling_df[best_metric].median(),
        'std': benchling_df[best_metric].std(),
        'min': benchling_df[best_metric].min(),
        'max': benchling_df[best_metric].max()
    }

    # determining the best performing tool
    if guidescan_stats['mean'] > benchling_stats['mean']:
        best_tool = 'GuideScan'
    else:
        best_tool = 'Benchling'

    # results
    print("GuideScan Summary:")
    print(guidescan_display[['gRNA-Seq', 'Specificity', 'Cutting Efficiency', 'M1_Rank', 'M2_Rank', 'M3_Rank']])
    print("\nBenchling Summary:")
    print(benchling_display[['Sequence', 'Specificity Score', 'Efficiency Score', 'M1_Rank', 'M2_Rank', 'M3_Rank']])
    print("\nMetric Correlations:")
    for metric, values in correlations.items():
        print(f"{metric}: x = {values['x']:.2f}, y = {values['y']:.2f}, z = {values['z']:.2f}")
        print(f"{metric} p-values: x_pval = {p_values[metric]['x_pval']:.2e}, y_pval = {p_values[metric]['y_pval']:.2e}")
    print(f"\nBest Metric: {best_metric}")

    print("\nGuideScan Best Metric Statistics:")
    print(guidescan_stats)
    print("\nBenchling Best Metric Statistics:")
    print(benchling_stats)
    print(f"\nBest Performing Tool: {best_tool}")

    return guidescan_display, benchling_display, correlations, p_values, best_metric, guidescan_stats, benchling_stats, best_tool

print("TP53 Gene")
benchling_csv = 'Benchling 7_1 TP53.csv'
guidescan_csv = 'GuideScan 7_1 TP53.csv'
results = calculate_metrics(benchling_csv, guidescan_csv)


print("------------------------------------------")
print("PAX1 Gene")
benchling_csv2 = 'BenchlingPAX1.csv'
guidescan_csv2 = 'GScanPAX1.csv'
results = calculate_metrics(benchling_csv2, guidescan_csv2)
