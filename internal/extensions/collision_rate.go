package extensions

import (
	"math"

	"github.com/wildstyl3r/lxgata"
	"github.com/wildstyl3r/stmc/internal/constants"
	"github.com/wildstyl3r/stmc/internal/model"
	"github.com/wildstyl3r/stmc/internal/utils"
)

const (
	NormalizedCollisionRateKey   = "NormalizedCollisionRate"
	NormalizedCollisionRateCIKey = "NormalizedCollisionRateCI"
	TotalCollisionsKey           = "TotalCollisions"
)

func NormalizedCollisionRate(model *model.Model) ([]string, []any, error) {
	collisions := map[lxgata.CollisionType][]float64{}
	collisionsCI := map[lxgata.CollisionType][]float64{}
	collCounters := map[lxgata.CollisionType]float64{}
	for key, val := range model.CollisionAtCell {
		if model.CollisionAtCell[key] != nil {
			collisions[key] = make([]float64, model.NumCells)
			collisionsCI[key] = make([]float64, model.NumCells)
			for xIndex := range model.NumCells {
				collisions[key][xIndex] = float64(utils.SumIntSlice(val[xIndex])) / (float64(model.Parameters.NElectrons) * model.XStep * model.Parameters.Pressure)
				collisionsCI[key][xIndex] = constants.Quantile95 * math.Sqrt(float64(utils.Variance(val[xIndex], true))) /
					(math.Sqrt(float64(model.Parameters.NElectrons)) * model.XStep * model.Parameters.Pressure)
				// 1.96 is double-sided quantile for 95% confidence
				collCounters[key] += float64(utils.SumIntSlice(model.CollisionAtCell[key][xIndex]))
			}
			collCounters[key] /= float64(model.Parameters.NElectrons)
		}
	}
	return []string{NormalizedCollisionRateKey, NormalizedCollisionRateCIKey, TotalCollisionsKey}, []any{collisions, collisionsCI, collCounters}, nil
	// datahub.Insert(NormalizedCollisionRateKey, collisions)
	// datahub.Insert(NormalizedCollisionRateCIKey, collisionsCI)
	// datahub.Insert(TotalCollisionsKey, collCounters)
}
