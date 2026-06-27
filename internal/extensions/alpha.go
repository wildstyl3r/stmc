package extensions

import (
	"fmt"
	"math"

	"github.com/schollz/progressbar/v3"
	"github.com/wildstyl3r/stmc/internal/config"
	"github.com/wildstyl3r/stmc/internal/constants"
	"github.com/wildstyl3r/stmc/internal/messages"
	"github.com/wildstyl3r/stmc/internal/model"
	"github.com/wildstyl3r/stmc/internal/utils"
)

type AlphaResult struct {
	utils.Result
	ElectricField                  float64 `csv:"ElectricField" units:"Length:-1,Voltage:1"`
	ElectricFieldPerPressure       float64 `csv:"ElectricFieldPerPressure" units:"Length:-1,Voltage:1,Pressure:-1"`
	ElectricFieldPerDensity        float64 `csv:"ElectricFieldPerPressure" units:"Length:-3,Voltage:1"`
	TownsendAlpha                  float64 `csv:"FirstIonizationCoefficient" units:"Length:-1,Pressure:-1"`
	TownsendAlphaMargin            float64 `csv:"FirstIonizationCoefficientMargin" units:"Length:-1,Pressure:-1"`
	TownsendEta                    float64 `csv:"AttachmentCoefficient" units:"Length:-1,Pressure:-1"`
	TownsendEtaMargin              float64 `csv:"AttachmentCoefficientMargin" units:"Length:-1,Pressure:-1"`
	EffectiveAlpha                 float64 `csv:"EffectiveIonizationCoefficient" units:"Length:-1,Pressure:-1"`
	EffectiveAlphaMargin           float64 `csv:"EffectiveIonizationCoefficientMargin" units:"Length:-1,Pressure:-1"`
	TownsendAlphaPerDensity        float64 `csv:"FirstIonizationCoefficientPerDensity" units:"Length:-4,Pressure:-1"`
	TownsendAlphaPerDensityMargin  float64 `csv:"FirstIonizationCoefficientPerDensityMargin" units:"Length:-4,Pressure:-1"`
	TownsendEtaPerDensity          float64 `csv:"AttachmentCoefficientPerDensity" units:"Length:-4,Pressure:-1"`
	TownsendEtaPerDensityMargin    float64 `csv:"AttachmentCoefficientPerDensityMargin" units:"Length:-4,Pressure:-1"`
	EffectiveAlphaPerDensity       float64 `csv:"EffectiveIonizationCoefficientPerDensity" units:"Length:-4,Pressure:-1"`
	EffectiveAlphaPerDensityMargin float64 `csv:"EffectiveIonizationCoefficientPerDensityMargin" units:"Length:-4,Pressure:-1"`
	DriftVelocity                  float64 `csv:"DriftVelocity" units:"Length:1,Time:-1"`
}

func AlphaCalculation(parameters *config.ModelParameters, outputDir, modelName string, logger messages.Logger) (dataRow, _ utils.ResultInterface, finalModel, _ *model.Model) {
	voltages := []float64{50, 75, 100, 150, 200, 250, 300}
	n := len(voltages)
	progress := progressbar.Default(int64(n))
	l1ps := make([]float64, n)
	gaps := make([]float64, n)
	ionizations := make([]int, n)
	attachments := make([]int, n)
	parameters.SupressSpinner = true
	k := 1.
	parameters.Pressure *= k
	parameters.ConstEField *= k
	for i := range voltages {
		parametersCopy := *parameters
		parametersCopy.GapLength = voltages[i] / -parameters.ConstEField
		gaps[i] = parametersCopy.GapLength
		finalModel = model.NewModel(&parametersCopy)
		LoadExtensions(finalModel.DataHub)
		finalModel.Run(func(m *model.Model) int {
			if m.TotalParticlesEmitted == 0 {
				return m.Parameters.NParticles
			} else {
				return 0
			}
		}, true, logger)
		l1ps[i] = math.Log1p(float64(finalModel.AnodeElectronCounter-finalModel.TotalParticlesEmitted) / (float64(finalModel.TotalParticlesEmitted)))
		ionizations[i] = utils.SumIntSlice(finalModel.IonizationCounters)
		attachments[i] = utils.SumIntSlice(finalModel.AttachmentCounters)
		progress.Describe(fmt.Sprintf("[ emitted: %d; @ anode: %d ; @ cathode: %d ; iz: %d; att: %d; eff alpha estim: %f]",
			finalModel.TotalParticlesEmitted, finalModel.AnodeElectronCounter, finalModel.CathodeElectronCounter, ionizations[i], attachments[i],
			config.FieldToSIeV(l1ps[i]/gaps[i]/k, []config.UnitElement{{Class: config.Length, Power: -1}}, finalModel.Parameters.OutputUnits(), false)))
		progress.Add(1)
	}

	nonNanIxs := utils.NonInfIndicies(l1ps)
	l1ps = utils.Select(l1ps, nonNanIxs)
	gaps = utils.Select(gaps, nonNanIxs)
	ionizations = utils.Select(ionizations, nonNanIxs)
	attachments = utils.Select(attachments, nonNanIxs)
	n = len(nonNanIxs)

	effectiveAlpha, _, effectiveAlphaVar, _ := utils.LinearRegressionMSEInferVariance(gaps, l1ps)

	totalIonizations := utils.SumIntSlice(ionizations)
	totalAttachments := utils.SumIntSlice(attachments)
	var globAlpha, globEta float64
	var varAlpha, varEta float64
	if totalAttachments == 0 {
		globAlpha = effectiveAlpha
		varAlpha = effectiveAlphaVar
	} else if totalIonizations == 0 {
		globEta = -effectiveAlpha
		varEta = effectiveAlphaVar
	} else {
		totalRatio := float64(totalAttachments) / float64(totalIonizations)
		globAlpha = effectiveAlpha / (1 - totalRatio)
		globEta = globAlpha * totalRatio

		jackAlphas := make([]float64, n)
		jackEtas := make([]float64, n)
		jackRatios := make([]float64, n)
		for i := range jackRatios {
			jackRatios[i] = float64(totalAttachments-attachments[i]) / float64(totalIonizations-ionizations[i])
			jackEffectiveAlpha, _ := utils.LinearRegressionMSE(utils.ExcludeView(i, gaps), utils.ExcludeView(i, l1ps))
			jackAlphas[i] = jackEffectiveAlpha / (1 - jackRatios[i])
			jackEtas[i] = jackAlphas[i] * jackRatios[i]
		}
		jackmeanAlpha := utils.Mean(jackAlphas)
		jackmeanEta := utils.Mean(jackEtas)
		for i := range jackAlphas {
			jdAlpha := jackAlphas[i] - jackmeanAlpha
			varAlpha += jdAlpha * jdAlpha
			jdEta := jackEtas[i] - jackmeanEta
			varEta += jdEta * jdEta
		}
		varAlpha *= float64(len(jackAlphas)-1) / float64(len(jackAlphas))
		varEta *= float64(len(jackEtas)-1) / float64(len(jackEtas))
	}

	cr := finalModel.CoreResult()
	cr.ReducedFieldAtCathode = finalModel.Parameters.ConstEField / finalModel.Parameters.GasDensity / constants.Townsend
	return &AlphaResult{
		Result:                         cr,
		ElectricField:                  -parameters.ConstEField,
		ElectricFieldPerPressure:       -parameters.ConstEField / cr.Pressure,
		TownsendAlpha:                  globAlpha / parameters.Pressure,
		TownsendAlphaMargin:            utils.EstimateMargin(0.95, varAlpha, float64(n)) / parameters.Pressure,
		TownsendEta:                    globEta / parameters.Pressure,
		TownsendEtaMargin:              utils.EstimateMargin(0.95, varEta, float64(n)) / parameters.Pressure,
		EffectiveAlpha:                 effectiveAlpha / parameters.Pressure,
		EffectiveAlphaMargin:           utils.EstimateMargin(0.95, effectiveAlphaVar, float64(n)) / parameters.Pressure,
		ElectricFieldPerDensity:        -parameters.ConstEField / parameters.GasDensity,
		TownsendAlphaPerDensity:        globAlpha / parameters.GasDensity,
		TownsendAlphaPerDensityMargin:  utils.EstimateMargin(0.95, varAlpha, float64(n)) / parameters.GasDensity,
		TownsendEtaPerDensity:          globEta / parameters.GasDensity,
		TownsendEtaPerDensityMargin:    utils.EstimateMargin(0.95, varEta, float64(n)) / parameters.GasDensity,
		EffectiveAlphaPerDensity:       effectiveAlpha / parameters.GasDensity,
		EffectiveAlphaPerDensityMargin: utils.EstimateMargin(0.95, effectiveAlphaVar, float64(n)) / parameters.GasDensity,
		DriftVelocity:                  0,
	}, nil, finalModel, nil
}
