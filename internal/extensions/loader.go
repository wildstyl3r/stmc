package extensions

import (
	"github.com/wildstyl3r/stmc/internal/datahub"
)

func Load() {
	extensions := map[string]datahub.Extension{
		NormalizedCollisionRateKey:   NormalizedCollisionRate,
		NormalizedCollisionRateCIKey: NormalizedCollisionRate,
		TotalCollisionsKey:           NormalizedCollisionRate,

		GlowDischargeDensityKey:           GlowDischargeDensity,
		SimplifiedGlowDischargeDensityKey: SimplifiedGlowDischargeDensity,
	}
	for name, ext := range extensions {
		datahub.Register(name, ext)
	}
}
