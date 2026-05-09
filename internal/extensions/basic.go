package extensions

import (
	"fmt"
	"strings"

	"github.com/wildstyl3r/stmc/internal/config"
	"github.com/wildstyl3r/stmc/internal/messages"
	"github.com/wildstyl3r/stmc/internal/model"
	"github.com/wildstyl3r/stmc/internal/utils"
)

func BasicCalculation(parameters *config.ModelParameters, outputDir, modelName string, logger messages.Logger) (dataRow, _ *utils.SheathResult, finalModel, _ *model.Model) {
	m := model.NewModel(parameters)

	LoadExtensions(m.DataHub)
	// additionStep := 0
	m.Run(func(m *model.Model) int {
		if m.TotalElectronsEmittedOnCathode == 0 {
			return m.Parameters.NElectrons
		} else {
			collisions := m.GetMetrics(SingleElectronDetailedCollisionRateKey).(map[string]utils.GriddedInterval)
			for substring := range m.Parameters.RequireCollisionRelativeMargin {
				for processName, counters := range collisions {
					if strings.Contains(processName, substring) {
						excessLeap := utils.RelativeExcessLeap(counters.Values)
						if excessLeap > 8.5 {
							fmt.Printf("\nrelative excess leap for %s is %f\n", processName, excessLeap)
							fmt.Printf("+ %d electrons\n", m.Parameters.AddByNElectrons)
							return m.Parameters.AddByNElectrons
						}
					}
				}
			}
			return 0
		}
	}, logger)
	cr := m.SheathResult()
	return &cr, nil, m, nil
}
