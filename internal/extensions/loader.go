package extensions

import (
	"github.com/wildstyl3r/stmc/internal/model"
)

func LoadExtensions(dh *model.DataHubType) {
	extensions := map[string]model.Extension{
		NormalizedCollisionRateKey:       NormalizedCollisionRate,
		NormalizedCollisionRateMarginKey: NormalizedCollisionRate,
		TotalCollisionsKey:               NormalizedCollisionRate,

		NormalizedWallLossRateKey:       NormalizedWallLossRate,
		NormalizedWallLossRateMarginKey: NormalizedWallLossRate,

		GlowDischargeDensityKey:           GlowDischargeDensity,
		SimplifiedGlowDischargeDensityKey: SimplifiedGlowDischargeDensity,

		GlowDischargeDensityPPHCKey: GlowDischargeDensityPPHC,
	}
	for name, ext := range extensions {
		dh.Register(name, ext)
	}
}
