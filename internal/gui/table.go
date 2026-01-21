package gui

import (
	"reflect"
	"strconv"

	"fyne.io/fyne/v2"
	"fyne.io/fyne/v2/container"
	"fyne.io/fyne/v2/dialog"
	"fyne.io/fyne/v2/widget"
)

func buildEditableTable(fvTable reflect.Value, valueType string, w fyne.Window) fyne.CanvasObject {
	rowsContainer := container.NewVBox()
	updateTableRows(rowsContainer, fvTable)

	addRowBtn := widget.NewButton("Add "+valueType, func() {
		keyEntry := widget.NewEntry()
		// keyEntry.SetPlaceHolder("Species as in cross-section file")
		valueEntry := widget.NewEntry()
		// valueEntry.SetPlaceHolder("Proportion factor")
		dialog.ShowForm("Add Row", "OK", "Cancel",
			[]*widget.FormItem{
				{Text: valueType, Widget: keyEntry},
				{Text: "Proportion", Widget: valueEntry},
			},
			func(ok bool) {
				if ok {
					key := keyEntry.Text
					val, _ := strconv.ParseFloat(valueEntry.Text, 64)

					fvTable.SetMapIndex(reflect.ValueOf(key), reflect.ValueOf(val))
					updateTableRows(rowsContainer, fvTable)
				}
			}, w)
	})

	return container.NewBorder(
		nil, nil, nil, nil,
		container.NewVBox(
			rowsContainer,
			addRowBtn,
		),
	)
}

func updateTableRows(c *fyne.Container, fvTable reflect.Value) {
	c.Objects = nil // Clear rows

	keys := fvTable.MapKeys()
	for _, key := range keys {
		keyStr := key.String()
		val := fvTable.MapIndex(key).Float()

		// Row: [Key Entry] [Value Entry] [🗑️ Delete]
		keyEntry := widget.NewEntry()
		keyEntry.SetText(keyStr)
		keyEntry.OnChanged = func(s string) {
			oldVal := fvTable.MapIndex(key).Float()
			fvTable.SetMapIndex(reflect.ValueOf(s), reflect.ValueOf(oldVal))
			fvTable.SetMapIndex(key, reflect.Value{}) // Delete old
			key = reflect.ValueOf(s)                  // Update key reference
		}

		valEntry := widget.NewEntry()
		valEntry.SetText(strconv.FormatFloat(val, 'f', 4, 64))
		valEntry.OnChanged = func(s string) {
			if n, err := strconv.Atoi(s); err == nil {
				fvTable.SetMapIndex(key, reflect.ValueOf(n))
			}
		}

		deleteBtn := widget.NewButton("Remove", func() {
			fvTable.SetMapIndex(key, reflect.Value{})
			updateTableRows(c, fvTable)
		})

		row := container.NewGridWithColumns(3,
			keyEntry,
			valEntry,
			deleteBtn,
		)
		c.Add(row)
	}

	c.Refresh()
}
