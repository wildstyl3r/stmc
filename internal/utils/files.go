package utils

import (
	"bufio"
	"fmt"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

func ReadFloatPairs(filename string) ([][]float64, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, fmt.Errorf("error opening file: %w", err)
	}
	defer file.Close()

	var result [][]float64

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		parts := strings.Fields(line)

		// Skip empty lines
		if len(parts) == 0 {
			continue
		}

		// Validate number of columns
		if len(parts) != 2 {
			return nil, fmt.Errorf("invalid format in line: %q - expected 2 numbers, got %d", line, len(parts))
		}

		// Convert to float64
		x, err := strconv.ParseFloat(parts[0], 64)
		if err != nil {
			return nil, fmt.Errorf("error parsing float in line %q: %w", line, err)
		}

		y, err := strconv.ParseFloat(parts[1], 64)
		if err != nil {
			return nil, fmt.Errorf("error parsing float in line %q: %w", line, err)
		}

		result = append(result, []float64{x, y})
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading file: %w", err)
	}

	return result, nil
}

func GetFilename(filePath string) string {
	// Get the base name (removes directory components)
	base := filepath.Base(filePath)

	// Remove the extension (everything after last dot)
	ext := filepath.Ext(base)

	// Trim the extension from base name
	nameWithoutExt := strings.TrimSuffix(base, ext)

	return nameWithoutExt
}

func OpenFile(makeDir bool, outputPath string, fileSuffix, modelName string) (*os.File, error) {
	if makeDir && fileSuffix != "" && fileSuffix != "." {
		os.MkdirAll(outputPath+fileSuffix, 0750)
		return os.Create(outputPath + fileSuffix + "/" + modelName + ".txt")
	} else {
		return os.Create(outputPath + modelName + "_" + fileSuffix + ".txt")
	}
}
