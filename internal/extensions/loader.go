package extensions

import (
	"github.com/wildstyl3r/stmc/internal/model"
)

func LoadExtensions(dh *model.DataHubType) {
	extensions := map[string]model.Extension{
		SingleElectronCollisionRateKey:       NormalizedCollisionRate,
		SingleElectronCollisionRateMarginKey: NormalizedCollisionRate,
		TotalCollisionsKey:                   NormalizedCollisionRate,

		NormalizedWallLossRateKey:       NormalizedWallLossRate,
		NormalizedWallLossRateMarginKey: NormalizedWallLossRate,

		AmbipolarDiffusionCoefficientKey: DiffusionConstants,
		CharacteristicDiffusionScaleKey:  DiffusionConstants,

		DensityPerCathodeElectronFluxKey:                DensityPerCathodeElectronFlux,
		GlowDischargeDensityKey:                         GlowDischargeDensity,
		CathodeElectronFluxKey:                          GlowDischargeDensity,
		CathodeElectronFluxMarginKey:                    GlowDischargeDensity,
		CathodeIonFluxKey:                               GlowDischargeDensity,
		CathodeIonFluxMarginKey:                         GlowDischargeDensity,
		FirstDensityPeakIndexKey:                        GlowDischargeDensity,
		SourceIntegralPerCathodeElectronFluxKey:         GlowDischargeDensity,
		SourceIntegralPerCathodeElectronFluxVarianceKey: GlowDischargeDensity,
	}
	for name, ext := range extensions {
		dh.Register(name, ext)
	}
}
