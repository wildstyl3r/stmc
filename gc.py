import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import sys
import csv

def townsend_model(ep, A, B):
    return A * np.exp(-B / ep)

def find_column_by_prefix(df, prefix):
    # Find all columns that start with the given prefix
    matches = [col for col in df.columns if str(col).strip().startswith(prefix)]
    
    if len(matches) == 0:
        print(f"Error: No column headers start with the prefix '{prefix}'.")
        sys.exit(1)
    elif len(matches) > 1:
        print(f"Error: Multiple columns match the prefix '{prefix}': {matches}. Please be more specific.")
        sys.exit(1)
        
    return matches[0]

def calculate_townsend_nlls(csv_path, ep_prefix, alpha_prefix):
    try:
        # Load CSV file
        df = pd.read_csv(csv_path)
        
        # Match columns based on prefixes
        ep_col = find_column_by_prefix(df, ep_prefix)
        alpha_col = find_column_by_prefix(df, alpha_prefix)
        
        # print(f"Matched '{ep_prefix}' -> Column: '{ep_col}'")
        # print(f"Matched '{alpha_prefix}' -> Column: '{alpha_col}'\n")
        
        # Extract data arrays
        ep = df[ep_col].values
        alpha_p = df[alpha_col].values
        
        # Force numeric types and drop NaNs/infinite values
        ep = pd.to_numeric(ep, errors='coerce')
        alpha_p = pd.to_numeric(alpha_p, errors='coerce')
        
        # Filter out zero, negative, or invalid NaN entries
        valid_idx = (ep > 0) & (alpha_p > 0) & (~np.isnan(ep)) & (~np.isnan(alpha_p))
        ep = ep[valid_idx]
        alpha_p = alpha_p[valid_idx]
        
        if len(ep) < 2:
            print("Error: Not enough valid positive numeric data points to perform fitting.")
            sys.exit(1)

        # Estimate initial guesses (p0) using a quick log-linear fallback
        log_y = np.log(alpha_p)
        inv_x = 1.0 / ep
        slope, intercept = np.polyfit(inv_x, log_y, 1)
        initial_A = np.exp(intercept)
        initial_B = -slope
        
        # Perform Non-Linear Least Squares (NLLS) optimization
        popt, pcov = curve_fit(
            townsend_model, 
            ep, 
            alpha_p, 
            p0=[initial_A, initial_B],
            bounds=(0, np.inf)
        )
        
        A_fit, B_fit = popt
        perr = np.sqrt(np.diag(pcov))
        
        # Calculate R²
        residuals = alpha_p - townsend_model(ep, A_fit, B_fit)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((alpha_p - np.mean(alpha_p))**2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0
        
        # Print optimized parameters
        # print("=== Townsend NLLS Fitting Results ===")
        print(f"Data points used: {len(ep)}")
        print(f"Constant A: {A_fit:.2f} ± {perr[0]:.2f} cm^-1 * Torr^-1")
        print(f"Constant B: {B_fit:.2f} ± {perr[1]:.2f} V * cm^-1 * Torr^-1")
        print(f"R² (Non-linear): {r_squared:.6f}")
        return popt
        
    except FileNotFoundError:
        print(f"Error: File not found at '{csv_path}'")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred during execution: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python townsend_nlls_flex.py <path_to_csv> <E/p_prefix> <alpha/p_prefix>")
        sys.exit(1)
    else:
        csv_file = sys.argv[1]
        ep_pref = sys.argv[2]
        alpha_pref = sys.argv[3]

        data = {
            "diss": [],
            "A": [],
            "B": [],
        }
        for r in ["D_0_0", "D_0_1", "D_0_2", "D_0_3", "D_0_4", "D_0_5", "D_0_6", "D_0_75", "D_0_8", "D_0_9", "D_1_0"]:
            a, b = calculate_townsend_nlls(csv_file.replace("D_0_0", r), ep_pref, alpha_pref)
            data["diss"].append(r)
            data["A"].append(a)
            data["B"].append(b)
        df = pd.DataFrame(data)
        df.to_csv('gc.csv', index=False)
