package extensions

import (
	"github.com/wildstyl3r/stmc/internal/config"
	"github.com/wildstyl3r/stmc/internal/messages"
	"github.com/wildstyl3r/stmc/internal/model"
	"github.com/wildstyl3r/stmc/internal/utils"
)

type IonMobilityResult struct {
	utils.Result
	IonMobility float64 `csv:"IonMobility" units:"Length:2,Voltage:-1,Time:-1"`
}

func IonMobilityCalculation(parameters *config.ModelParameters, outputDir, modelName string, logger messages.Logger) (dataRow, _ utils.ResultInterface, finalModel, _ *model.Model) {
	m := model.NewModel(parameters)

	LoadExtensions(m.DataHub)
	// additionStep := 0
	m.Run(func(m *model.Model) int {
		if m.TotalParticlesEmitted == 0 {
			return m.Parameters.NParticles
		} else {
			return 0
		}
	}, false, logger)
	imr := IonMobilityResult{
		Result:      m.CoreResult(),
		IonMobility: 0,
	}
	return &imr, nil, m, nil
}
