package gui

import (
	"maps"
	"slices"

	"fyne.io/fyne/v2"
	"fyne.io/fyne/v2/container"
	"fyne.io/fyne/v2/layout"
	"fyne.io/fyne/v2/widget"
	"github.com/wildstyl3r/stmc/internal/output"
)

func CheckboxizeDataFlags(df *output.DataFlags) fyne.CanvasObject {
	names := slices.Sorted(maps.Keys(df.Sequentials))
	checkBoxes := []fyne.CanvasObject{}
	for _, name := range names {
		check := widget.NewCheck(name, func(b bool) {
			*df.Sequentials[name].SaveFlag = b
		})
		check.Checked = *df.Sequentials[name].SaveFlag
		checkBoxes = append(checkBoxes, check)
	}
	return container.New(layout.NewRowWrapLayout(), checkBoxes...)
}
