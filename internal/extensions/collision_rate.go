package extensions

import (
	"github.com/wildstyl3r/lxgata"
	"github.com/wildstyl3r/stmc/internal/constants"
	"github.com/wildstyl3r/stmc/internal/model"
	"github.com/wildstyl3r/stmc/internal/utils"
)

const (
	NormalizedCollisionRateKey       = "NormalizedCollisionRate"
	NormalizedCollisionRateMarginKey = "NormalizedCollisionRateMargin"
	TotalCollisionsKey               = "TotalCollisions"
)

func NormalizedCollisionRate(model *model.Model) ([]string, []any, error) {
	collisions := map[lxgata.CollisionType][]float64{}
	collisionsMargin := map[lxgata.CollisionType][]float64{}
	collCounters := map[lxgata.CollisionType]float64{}
	for key, val := range model.CollisionAtCell {
		if model.CollisionAtCell[key] != nil {
			collisions[key] = make([]float64, model.NumCells)
			collisionsMargin[key] = make([]float64, model.NumCells)
			for xIndex := range model.NumCells {
				xCollisionsMean, xCollisionsVariance := utils.MeanAndVariance(val[xIndex], true)
				collisions[key][xIndex] = xCollisionsMean / (model.XStep * model.Parameters.Pressure)
				collisionsMargin[key][xIndex] = utils.NormalMargin(constants.Quantile95, xCollisionsVariance, model.Parameters.NElectrons) / (model.XStep * model.Parameters.Pressure)
				// 1.96 is double-sided quantile for 95% confidence
				collCounters[key] += float64(utils.SumIntSlice(model.CollisionAtCell[key][xIndex]))
			}
			collCounters[key] /= float64(model.Parameters.NElectrons)
		}
	}
	return []string{NormalizedCollisionRateKey, NormalizedCollisionRateMarginKey, TotalCollisionsKey}, []any{collisions, collisionsMargin, collCounters}, nil
}
