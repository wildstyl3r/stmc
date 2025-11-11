package extensions

import (
	"github.com/wildstyl3r/stmc/internal/model"
	"github.com/wildstyl3r/stmc/internal/utils"
)

const (
	NormalizedWallLossRateKey       = "NormalizedWallLossRate"
	NormalizedWallLossRateMarginKey = "NormalizedWallLossRateMargin"
)

func NormalizedWallLossRate(model *model.Model) ([]string, []any, error) {
	wallLosses := make([]float64, model.NumCells)
	wallLossesMargin := make([]float64, model.NumCells)
	for xIndex := range model.NumCells {
		xWallLossesMean, xWallLossesVariance := model.WallLossAtCell[xIndex].MeanAndVariance()
		wallLosses[xIndex] = xWallLossesMean / (model.XStep * model.Parameters.Pressure)
		wallLossesMargin[xIndex] = utils.NormalMargin(0.95, xWallLossesVariance, model.TotalElectronsEmittedOnCathode) / (model.XStep * model.Parameters.Pressure)
	}
	return []string{NormalizedWallLossRateKey, NormalizedWallLossRateMarginKey}, []any{wallLosses, wallLossesMargin}, nil
}
