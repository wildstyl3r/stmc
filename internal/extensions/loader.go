package extensions

import (
	"github.com/wildstyl3r/stmc/internal/model"
)

func LoadExtensions(dh *model.DataHubType) {
	extensions := map[string]model.Extension{
		SingleElectronCollisionRateKey:         SingleElectronCollisionRate,
		SingleElectronCollisionRateMarginKey:   SingleElectronCollisionRate,
		SingleElectronCollisionsOutsideKey:     SingleElectronCollisionRate,
		SingleElectronDetailedCollisionRateKey: SingleElectronCollisionRate,
		// SingleElectronDetailedCollisionRateLowerMarginKey: SingleElectronCollisionRate,
		// SingleElectronDetailedCollisionRateUpperMarginKey: SingleElectronCollisionRate,

		NormalizedWallLossRateKey:       NormalizedWallLossRate,
		NormalizedWallLossRateMarginKey: NormalizedWallLossRate,

		AmbipolarDiffusionCoefficientKey: DiffusionConstants,
		CharacteristicDiffusionScaleKey:  DiffusionConstants,

		DensityPerCathodeElectronFluxKey:                DensityPerCathodeElectronFlux,
		CathodeElectronCurrentFractionKey:               GlowDischargeDensity,
		CathodeElectronCurrentFractionMarginKey:         GlowDischargeDensity,
		CathodeIonCurrentFractionKey:                    GlowDischargeDensity,
		CathodeIonCurrentFractionMarginKey:              GlowDischargeDensity,
		FirstDensityPeakIndexKey:                        GlowDischargeDensity,
		SourceIntegralPerCathodeElectronFluxKey:         GlowDischargeDensity,
		SourceIntegralPerCathodeElectronFluxVarianceKey: GlowDischargeDensity,
	}
	for name, ext := range extensions {
		dh.Register(name, ext)
	}
}
