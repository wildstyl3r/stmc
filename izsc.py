import csv
import sys


def process_csv(input_file, output_file, sub_value, scale_value):
    try:
        with open(input_file, mode="r", newline="") as infile, open(
            output_file, mode="w", newline=""
        ) as outfile:

            reader = csv.reader(infile)
            writer = csv.writer(outfile)

            for row_idx, row in enumerate(reader, start=1):
                # Skip empty rows
                if not row:
                    continue

                # Ensure the row has exactly 2 columns
                if len(row) != 2:
                    print(
                        f"Warning: Skipping row {row_idx} because it does not have 2 columns: {row}"
                    )
                    continue

                try:
                    # Convert strings to floats, apply math, and update
                    col1 = float(row[0]) - sub_value
                    col2 = float(row[1]) * scale_value

                    writer.writerow([col1, col2])

                except ValueError:
                    # Handle headers or non-numeric text gracefully
                    if row_idx == 1:
                        print(
                            f"Skipping header row or non-numeric data at line 1: {row}"
                        )
                        # Optional: writer.writerow(row) # Uncomment to keep headers
                    else:
                        print(
                            f"Error: Non-numeric data found at row {row_idx}: {row}"
                        )

        print(f"Success! Processed data saved to '{output_file}'")

    except FileNotFoundError:
        print(f"Error: The file '{input_file}' was not found.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


# --- Example Usage ---
if __name__ == "__main__":
    # Define your files and modification values here
    INPUT_PATH = sys.argv[1]
    OUTPUT_PATH = sys.argv[2]
    EXCITATION_ENERGY = float(sys.argv[3])
    IONIZATION_ENERGY = float(sys.argv[4])
    SCALE_BY = (IONIZATION_ENERGY/(IONIZATION_ENERGY-EXCITATION_ENERGY))**2

    process_csv(INPUT_PATH, OUTPUT_PATH, EXCITATION_ENERGY, SCALE_BY)
