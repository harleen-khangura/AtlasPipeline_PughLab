import numpy as np
from statsmodels.stats.multitest import multipletests
from scipy.stats import nbinom

# Function to perform negative binomial test
def negative_binomial_test(counts):
    successes = counts.sum(axis=1)  # Total successes
    trials = counts.shape[1]  # Total trials (samples)

    # Calculate mean and variance
    mu = np.mean(counts, axis=1)
    var = np.var(counts, axis=1)

    # Initialize arrays for size and prob
    size = np.zeros_like(mu)
    prob = np.zeros_like(mu)

    # Avoid division by zero
    for i in range(len(mu)):
        if var[i] > mu[i]:  # Ensure variance is greater than mean
            size[i] = mu[i]**2 / (var[i] - mu[i])  # size parameter
            prob[i] = mu[i] / var[i]  # probability of success
        else:
            size[i] = np.nan  # Set to NaN if not valid
            prob[i] = np.nan

    # Calculate p-values, ignoring NaN values
    p_values = 1 - nbinom.cdf(0, size, prob)
    p_values = np.nan_to_num(p_values)  # Replace NaN with zero for p-values

    return p_values

# Load the density data
IgG_data = pd.read_csv("aggregated_results_IgG.txt", sep="\t", header=None)
cwpair_data = pd.read_csv("aggregated_results_cwpair.txt", sep="\t", header=None)

# Convert density columns to numeric, forcing errors to NaN
IgG_data[1] = pd.to_numeric(IgG_data[1], errors='coerce')
cwpair_data[1] = pd.to_numeric(cwpair_data[1], errors='coerce')

# Combine the data into a single DataFrame, dropping rows with NaN values
count_matrix = pd.DataFrame({
    'IgG': IgG_data[1],
    'cwpair': cwpair_data[1]
}).dropna()

# Perform negative binomial test
p_values = negative_binomial_test(count_matrix.values)

# Apply Benjamini-Hochberg correction
adjusted_p_values = multipletests(p_values, method='fdr_bh')[1]

# Add to DataFrame
results = pd.DataFrame({
    'Distance': IgG_data[0][:len(p_values)],  # Assuming first column is the distance
    'p_value': p_values,
    'adjusted_p_value': adjusted_p_values
})

# Filter significant peaks
significant_peaks = results[results['adjusted_p_value'] < 0.05]  # Adjust threshold as needed

# Count significant peaks
num_significant_peaks = significant_peaks.shape[0]
print(f"Number of significant peaks: {num_significant_peaks}")

# Optional: Save the results to a file
results.to_csv("significant_peaks_results.csv", index=False)












