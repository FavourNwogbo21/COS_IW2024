import pandas as pd

def analyze_score_ranges(benchling_csv, guidescan_csv):
    # loading CSV files into dataframes
    benchling_df = pd.read_csv(benchling_csv)
    guidescan_df = pd.read_csv(guidescan_csv)
    
    # now convert the columns to numeric to access the integers (they are currently strings)
    benchling_df['Specificity Score'] = pd.to_numeric(benchling_df['Specificity Score'], errors='coerce')
    benchling_df['Efficiency Score'] = pd.to_numeric(benchling_df['Efficiency Score'], errors='coerce')
    guidescan_df['Specificity'] = pd.to_numeric(guidescan_df['Specificity'], errors='coerce')
    guidescan_df['Cutting Efficiency'] = pd.to_numeric(guidescan_df['Cutting Efficiency'], errors='coerce')
    
    benchling_df = benchling_df.dropna(subset=['Specificity Score', 'Efficiency Score'])
    guidescan_df = guidescan_df.dropna(subset=['Specificity', 'Cutting Efficiency'])
    
    # normalize GuideScan scores to a scale of 100
    guidescan_df['Specificity'] *= 100
    guidescan_df['Cutting Efficiency'] *= 100
    bins = [0, 20, 40, 60, 80, 100]
    labels = ['0-20', '20-40', '40-60', '60-80', '80-100']
    
    # putting the scores into ranges
    benchling_df['Specificity Range'] = pd.cut(benchling_df['Specificity Score'], bins=bins, labels=labels, include_lowest=True)
    benchling_df['Efficiency Range'] = pd.cut(benchling_df['Efficiency Score'], bins=bins, labels=labels, include_lowest=True)
    guidescan_df['Specificity Range'] = pd.cut(guidescan_df['Specificity'], bins=bins, labels=labels, include_lowest=True)
    guidescan_df['Efficiency Range'] = pd.cut(guidescan_df['Cutting Efficiency'], bins=bins, labels=labels, include_lowest=True)
    
    benchling_specificity_counts = benchling_df['Specificity Range'].value_counts().sort_index()
    benchling_efficiency_counts = benchling_df['Efficiency Range'].value_counts().sort_index()
    guidescan_specificity_counts = guidescan_df['Specificity Range'].value_counts().sort_index()
    guidescan_efficiency_counts = guidescan_df['Efficiency Range'].value_counts().sort_index()

    benchling_specificity_percentages = (benchling_specificity_counts / benchling_specificity_counts.sum()) * 100
    benchling_efficiency_percentages = (benchling_efficiency_counts / benchling_efficiency_counts.sum()) * 100
    guidescan_specificity_percentages = (guidescan_specificity_counts / guidescan_specificity_counts.sum()) * 100
    guidescan_efficiency_percentages = (guidescan_efficiency_counts / guidescan_efficiency_counts.sum()) * 100
    
    # Print results
    print("Benchling Specificity Score Ranges:")
    print(benchling_specificity_counts)
    print("Percentages:")
    print(benchling_specificity_percentages)
    print("\nBenchling Efficiency Score Ranges:")
    print(benchling_efficiency_counts)
    print("Percentages:")
    print(benchling_efficiency_percentages)
    print("\nGuideScan Specificity Score Ranges:")
    print(guidescan_specificity_counts)
    print("Percentages:")
    print(guidescan_specificity_percentages)
    print("\nGuideScan Efficiency Score Ranges:")
    print(guidescan_efficiency_counts)
    print("Percentages:")
    print(guidescan_efficiency_percentages)

    return {
        'benchling_specificity_ranges': benchling_specificity_counts,
        'benchling_specificity_percentages': benchling_specificity_percentages,
        'benchling_efficiency_ranges': benchling_efficiency_counts,
        'benchling_efficiency_percentages': benchling_efficiency_percentages,
        'guidescan_specificity_ranges': guidescan_specificity_counts,
        'guidescan_specificity_percentages': guidescan_specificity_percentages,
        'guidescan_efficiency_ranges': guidescan_efficiency_counts,
        'guidescan_efficiency_percentages': guidescan_efficiency_percentages,
    }


print("TP53 Gene")
benchling_csv = 'Benchling 7_1 TP53.csv'
guidescan_csv = 'GuideScan 7_1 TP53.csv'
results = analyze_score_ranges(benchling_csv, guidescan_csv)

print("-----------------------------------")
print("PAX1 Gene")
benchling_csv2 = 'BenchlingPAX1.csv'
guidescan_csv2 = 'GScanPAX1.csv'
results = analyze_score_ranges(benchling_csv2, guidescan_csv2)

