package gui

import (
	"reflect"

	"fyne.io/fyne/v2"
	"fyne.io/fyne/v2/container"
	"fyne.io/fyne/v2/dialog"
	"fyne.io/fyne/v2/widget"
)

func buildMapAppTabs(fvMap reflect.Value, tabType string, w fyne.Window, meta FieldMeta) fyne.CanvasObject {
	if fvMap.Kind() != reflect.Map || fvMap.Type().Key().Kind() != reflect.String {
		return widget.NewLabel("Unsupported map type")
	}

	tabs := container.NewAppTabs()
	updateTabs(tabs, fvMap, w, meta)

	addBtn := widget.NewButton("Add "+tabType, func() {
		entry := widget.NewEntry()
		entry.SetPlaceHolder("Enter tab key...")

		dialog.ShowForm("New "+tabType, "OK", "Cancel",
			[]*widget.FormItem{
				widget.NewFormItem("Key", entry),
			},
			func(confirmed bool) {
				if confirmed && entry.Text != "" {
					newKey := entry.Text
					newVal := reflect.New(fvMap.Type().Elem().Elem())
					fvMap.SetMapIndex(reflect.ValueOf(newKey), newVal)
					updateTabs(tabs, fvMap, w, meta)
				}
			}, w)
	})

	deleteBtn := widget.NewButton("Delete Selected "+tabType, func() {
		if tabs.Selected() == nil {
			return
		}
		tabText := tabs.Selected().Text

		keyVal := reflect.ValueOf(tabText)
		fvMap.SetMapIndex(keyVal, reflect.Value{})
		updateTabs(tabs, fvMap, w, meta)
	})

	controls := container.NewHBox(addBtn, deleteBtn)

	return container.NewBorder(nil, controls, nil, nil, tabs)
}

func updateTabs(tabs *container.AppTabs, fvMap reflect.Value, w fyne.Window, meta FieldMeta) {
	tabs.Items = nil

	keys := fvMap.MapKeys()
	for _, key := range keys {
		keyStr := key.String()
		if keyStr == "" {
			continue
		}

		subForm := widget.NewForm()
		subVal := fvMap.MapIndex(key)
		buildNestedForm(subForm, subVal, w, meta)

		tabs.Append(container.NewTabItem(keyStr, subForm))
	}

	tabs.Refresh()
}
