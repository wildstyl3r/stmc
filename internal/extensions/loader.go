package extensions

import (
	"github.com/wildstyl3r/stmc/internal/model"
)

func LoadExtensions(dh *model.DataHubType) {
	extensions := map[string]model.Extension{
		NormalizedCollisionRateKey:   NormalizedCollisionRate,
		NormalizedCollisionRateCIKey: NormalizedCollisionRate,
		TotalCollisionsKey:           NormalizedCollisionRate,

		GlowDischargeDensityKey:           GlowDischargeDensity,
		SimplifiedGlowDischargeDensityKey: SimplifiedGlowDischargeDensity,
	}
	for name, ext := range extensions {
		dh.Register(name, ext)
	}
}
