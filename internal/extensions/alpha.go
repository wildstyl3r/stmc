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
	ElectricField            float64 `csv:"ElectricField" units:"Length:-1,Voltage:1"`
	ElectricFieldPerPressure float64 `csv:"ElectricFieldPerPressure" units:"Length:-1,Voltage:1,Pressure:-1"`
	TownsendAlpha            float64 `csv:"FirstIonizationCoefficient" units:"Length:-1,Pressure:-1"`
	TownsendAlphaMargin      float64 `csv:"FirstIonizationCoefficientMargin" units:"Length:-1,Pressure:-1"`
	TownsendEta              float64 `csv:"AttachmentCoefficient" units:"Length:-1,Pressure:-1"`
	TownsendEtaMargin        float64 `csv:"AttachmentCoefficientMargin" units:"Length:-1,Pressure:-1"`
	EffectiveAlpha           float64 `csv:"EffectiveIonizationCoefficient" units:"Length:-1,Pressure:-1"`
	EffectiveAlphaMargin     float64 `csv:"EffectiveIonizationCoefficientMargin" units:"Length:-1,Pressure:-1"`
	DriftVelocity            float64 `csv:"DriftVelocity" units:"Length:1,Time:-1"`
}

func AlphaCalculationR(parameters *config.ModelParameters, outputDir, modelName string, logger messages.Logger) (dataRow, _ *AlphaResult, finalModel, _ *model.Model) {
	n := 7
	progress := progressbar.Default(int64(n))
	voltages := []float64{50, 100, 150, 200, 250, 300, 350}
	l1ps := make([]float64, n)
	gaps := make([]float64, n)
	ionizations := make([]int, n)
	attachments := make([]int, n)
	jackRatios := make([]float64, n)
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
			if m.TotalElectronsEmittedOnCathode == 0 {
				return m.Parameters.NElectrons
			} else {
				return 0
			}
		}, logger)
		l1ps[i] = math.Log1p(float64(finalModel.AnodeElectronCounter-finalModel.TotalElectronsEmittedOnCathode) / (float64(finalModel.TotalElectronsEmittedOnCathode)))
		ionizations[i] = utils.SumIntSlice(finalModel.IonizationCounters)
		attachments[i] = utils.SumIntSlice(finalModel.AttachmentCounters)
		progress.Describe(fmt.Sprintf("[ emitted: %d ; arrived at anode: %d ; arrived at cathode: %d ; ionizations: %d; attachments: %d; rough effective alpha estimation: %f]",
			finalModel.TotalElectronsEmittedOnCathode, finalModel.AnodeElectronCounter, finalModel.CathodeElectronCounter, ionizations[i], attachments[i],
			config.FieldToSIeV(l1ps[i]/gaps[i]/k, []config.UnitElement{{Class: config.Length, Power: -1}}, finalModel.Parameters.OutputUnits(), false)))
		progress.Add(1)
	}

	effectiveAlpha, _, effectiveAlphaVar, _ := utils.LinearRegressionMSEInferVariance(gaps, l1ps)

	totalIonizations := utils.SumIntSlice(ionizations)
	totalAttachments := utils.SumIntSlice(attachments)
	totalRatio := float64(totalIonizations) / float64(totalAttachments)
	globeta := effectiveAlpha / (totalRatio - 1)
	globalpha := globeta * totalRatio

	jackAlphas := make([]float64, n)
	jackEtas := make([]float64, n)
	for i := range jackRatios {
		jackRatios[i] = float64(totalIonizations-ionizations[i]) / float64(totalAttachments-attachments[i])
		jackEffectiveAlpha, _ := utils.LinearRegressionMSE(utils.ExcludeView(i, gaps), utils.ExcludeView(i, l1ps))
		jackEtas[i] = jackEffectiveAlpha / (jackRatios[i] - 1)
		jackAlphas[i] = jackEtas[i] * jackRatios[i]
	}
	jackmeanAlpha := utils.Mean(jackAlphas)
	jackmeanEta := utils.Mean(jackEtas)
	var jackvarAlpha, jackvarEta float64
	for i := range jackAlphas {
		jdAlpha := jackAlphas[i] - jackmeanAlpha
		jackvarAlpha += jdAlpha * jdAlpha
		jdEta := jackEtas[i] - jackmeanEta
		jackvarEta += jdEta * jdEta
	}
	jackvarAlpha *= float64(len(jackAlphas)-1) / float64(len(jackAlphas))
	jackvarEta *= float64(len(jackEtas)-1) / float64(len(jackEtas))

	cr := finalModel.CoreResult()
	cr.ReducedFieldAtCathode = finalModel.Parameters.ConstEField / finalModel.Parameters.GasDensity / constants.Townsend
	return &AlphaResult{
		Result:                   cr,
		ElectricField:            -parameters.ConstEField,
		ElectricFieldPerPressure: -parameters.ConstEField / cr.Pressure,
		TownsendAlpha:            globalpha / parameters.Pressure,
		TownsendAlphaMargin:      utils.EstimateMargin(0.95, jackvarAlpha, float64(n)) / parameters.Pressure,
		TownsendEta:              globeta / parameters.Pressure,
		TownsendEtaMargin:        utils.EstimateMargin(0.95, jackvarEta, float64(n)) / parameters.Pressure,
		EffectiveAlpha:           effectiveAlpha / parameters.Pressure,
		EffectiveAlphaMargin:     utils.EstimateMargin(0.95, effectiveAlphaVar, float64(n)) / parameters.Pressure,
		DriftVelocity:            0,
	}, nil, finalModel, nil
}
