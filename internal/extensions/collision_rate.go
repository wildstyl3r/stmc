package extensions

import (
	"github.com/wildstyl3r/lxgata"
	"github.com/wildstyl3r/stmc/internal/model"
	"github.com/wildstyl3r/stmc/internal/utils"
)

const (
	SingleElectronCollisionRateKey         = "SingleElectronCollisionRate"
	SingleElectronCollisionRateMarginKey   = "SingleElectronCollisionRateMargin"
	SingleElectronCollisionsOutsideKey     = "SingleElectronCollisionsOutside"
	SingleElectronDetailedCollisionRateKey = "SingleElectronDetailedCollisionRate"
	// SingleElectronDetailedCollisionRateLowerMarginKey = "SingleElectronDetailedCollisionRateLowerMargin"
	// SingleElectronDetailedCollisionRateUpperMarginKey = "SingleElectronDetailedCollisionRateUpperMargin"
)

func SingleElectronCollisionRate(model *model.Model) ([]string, []any, error) {
	collisions := map[lxgata.CollisionType]utils.GriddedInterval{}
	collisionsMargin := map[lxgata.CollisionType]utils.GriddedInterval{}
	collisionsOutside := map[lxgata.CollisionType]utils.GriddedInterval{}
	detailedCollisions := map[string]utils.GriddedInterval{}
	// detailedCollisionsLowerMargin := map[string]utils.GriddedInterval{}
	// detailedCollisionsUpperMargin := map[string]utils.GriddedInterval{}
	for collType, distribution := range model.CollisionAtCell {
		if distribution != nil {
			collisionsOfType := utils.GriddedInterval{Values: make([]float64, model.NumCells), Step: model.XStep}
			marginsOfType := utils.GriddedInterval{Values: make([]float64, model.NumCells), Step: model.XStep}
			outside := utils.GriddedInterval{Values: make([]float64, model.NumCells), Step: model.XStep}
			for xIndex := range model.NumCells {
				collisions, margin := distribution[xIndex].MeanWithErrorMargin(0.95)
				collisionsOfType.Values[xIndex], marginsOfType.Values[xIndex] = collisions/model.XStep, margin/model.XStep
			}
			collisionsOfType.Values = append(collisionsOfType.Values, collisionsOfType.Values[len(collisionsOfType.Values)-1])
			collisionsOfType.Values = append(collisionsOfType.Values, collisionsOfType.Values[len(collisionsOfType.Values)-1])
			marginsOfType.Values = append(marginsOfType.Values, marginsOfType.Values[len(marginsOfType.Values)-1])
			marginsOfType.Values = append(marginsOfType.Values, marginsOfType.Values[len(marginsOfType.Values)-1])
			collisions[collType] = collisionsOfType
			collisionsOutside[collType] = outside
			collisionsMargin[collType] = marginsOfType
		}
	}
	for key, val := range model.DetailedCollisionAtCell {
		if val != nil {
			currentCollisions := utils.GriddedInterval{Values: make([]float64, model.NumCells), Step: model.XStep}
			// lowerMargins := utils.GriddedInterval{Values: make([]float64, model.NumCells), Step: model.XStep}
			// upperMargins := utils.GriddedInterval{Values: make([]float64, model.NumCells), Step: model.XStep}
			for xIndex := range model.NumCells {
				xCollisionsMean := val[xIndex].Mean()
				currentCollisions.Values[xIndex] = xCollisionsMean / model.XStep
				// lm, um := utils.PoissonMargin(0.99, val[xIndex])
				// lowerMargins.Values[xIndex] = lm
				// upperMargins.Values[xIndex] = um
				// 1.96 is double-sided quantile for 95% confidence
			}
			currentCollisions.Values = append(currentCollisions.Values, currentCollisions.Values[len(currentCollisions.Values)-1])
			// lowerMargins.Values = append(lowerMargins.Values, lowerMargins.Values[len(lowerMargins.Values)-1])
			detailedCollisions[key] = currentCollisions
			// detailedCollisionsLowerMargin[key] = lowerMargins
			// detailedCollisionsUpperMargin[key] = upperMargins
		}
	}

	return []string{
			SingleElectronCollisionRateKey,
			SingleElectronCollisionRateMarginKey,
			SingleElectronDetailedCollisionRateKey,
			// SingleElectronDetailedCollisionRateLowerMarginKey,
			// SingleElectronDetailedCollisionRateUpperMarginKey,
			SingleElectronCollisionsOutsideKey,
		}, []any{
			collisions,
			collisionsMargin,
			detailedCollisions,
			// detailedCollisionsLowerMargin,
			// detailedCollisionsUpperMargin,
			collisionsOutside,
		}, nil
}
