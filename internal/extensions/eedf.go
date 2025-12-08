package extensions

import "github.com/wildstyl3r/stmc/internal/model"

const ElectronEnergyDistributionFunctionKey = "ElectronEnergyDistributionFunction"

func ElectronEnergyDistributionFunction(m *model.Model) ([]string, []any, error) {
	return []string{ElectronEnergyDistributionFunctionKey}, []any{}, nil
}
