package extensions

import (
	"github.com/wildstyl3r/lxgata"
	"github.com/wildstyl3r/stmc/internal/model"
	"github.com/wildstyl3r/stmc/internal/utils"
)

const (
	SingleElectronCollisionRateKey       = "NormalizedCollisionRate"
	SingleElectronCollisionRateMarginKey = "NormalizedCollisionRateMargin"
	TotalCollisionsKey                   = "TotalCollisions"
)

func NormalizedCollisionRate(model *model.Model) ([]string, []any, error) {
	collisions := map[lxgata.CollisionType]utils.GriddedInterval{}
	collisionsMargin := map[lxgata.CollisionType]utils.GriddedInterval{}
	collCounters := map[lxgata.CollisionType]float64{}
	for key, val := range model.CollisionAtCell {
		if model.CollisionAtCell[key] != nil {
			currentCollisions := utils.GriddedInterval{Values: make([]float64, model.NumCells), Step: model.XStep}
			currentMargins := utils.GriddedInterval{Values: make([]float64, model.NumCells), Step: model.XStep}
			for xIndex := range model.NumCells {
				xCollisionsMean, xCollisionsVariance := utils.MeanAndVariance(val[xIndex], true)
				currentCollisions.Values[xIndex] = xCollisionsMean / model.XStep
				currentMargins.Values[xIndex] = utils.StudentedMargin(0.95, xCollisionsVariance, model.TotalElectronsPassed) / model.XStep
				// 1.96 is double-sided quantile for 95% confidence
				collCounters[key] += float64(utils.SumIntSlice(model.CollisionAtCell[key][xIndex]))
			}
			currentCollisions.Values = append(currentCollisions.Values, currentCollisions.Values[len(currentCollisions.Values)-1])
			currentMargins.Values = append(currentMargins.Values, currentMargins.Values[len(currentMargins.Values)-1])
			collisions[key] = currentCollisions
			collisionsMargin[key] = currentMargins
			collCounters[key] /= float64(model.TotalElectronsPassed)
		}
	}
	return []string{SingleElectronCollisionRateKey, SingleElectronCollisionRateMarginKey, TotalCollisionsKey}, []any{collisions, collisionsMargin, collCounters}, nil
}
