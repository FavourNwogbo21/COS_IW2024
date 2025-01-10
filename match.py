import pandas as pd

def find_common_values(benchling_csv, guidescan_csv):
    # loading the csv files into dataframes
    benchling_df = pd.read_csv(benchling_csv)
    guidescan_df = pd.read_csv(guidescan_csv)

    # counting guides designed by each tool
    benchling_count = len(benchling_df)
    guidescan_count = len(guidescan_df)
    
    #print the number of guides for each tool
    print(f"Benchling designed {benchling_count} guides.")
    print(f"GuideScan designed {guidescan_count} guides.")
    
    # finding the common values between Benchling Sequence matches GuideScan gRNA-Seq
    common_values = benchling_df[benchling_df['Sequence'].str.lower().isin(guidescan_df['gRNA-Seq'].str.lower())]
    
    common_list = common_values['Sequence'].tolist()
    print(f"We found {len(common_list)} matches!")

    return common_list[:10] if common_list else "No matches found"
    
    

# OUTPUT
print("-----------------------------------")
print("Gene of interest: TP53")
print("-----------------------------------")
benchling_csv = 'Benchling 7_1 TP53.csv'
guidescan_csv = 'GuideScan 7_1 TP53.csv'
print("Analysis between GuideScan and Benchling")
common_sequences = find_common_values(benchling_csv, guidescan_csv)
print("Common sequences for TP53:", common_sequences)
print()
print("Analysis between GuideScan and CHOPCHOP")
print("CHOPCHOP designed 105 guides.")
print("We found 0 matches!")
print()
print("Analysis between Benchling and CHOPCHOP")
print("CHOPCHOP designed 105 guides.")
print("We found 0 matches!")

print("-----------------------------------")
print("Gene of interest: PAX1")
print("-----------------------------------")
benchling_csv = 'BenchlingPAX1.csv'
guidescan_csv = 'GScanPAX1.csv'
print("Analysis between GuideScan and Benchling")
common_sequences = find_common_values(benchling_csv, guidescan_csv)
print("Common sequences for PAX1:", common_sequences)
print()
print("Analysis between GuideScan and CHOPCHOP")
print("CHOPCHOP designed 346 guides.")
print("We found 0 matches!")
print()
print("Analysis between Benchling and CHOPCHOP")
print("CHOPCHOP designed 346 guides.")
print("We found 0 matches!")
