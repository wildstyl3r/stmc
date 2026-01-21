package gui

import (
	"maps"
	"reflect"
	"slices"

	"fyne.io/fyne/v2"
	"fyne.io/fyne/v2/container"
	"fyne.io/fyne/v2/dialog"
	"fyne.io/fyne/v2/widget"
)

func buildNestedForm(parent *widget.Form, val reflect.Value, w fyne.Window, meta FieldMeta) {
	val = reflect.Indirect(val)
	typ := val.Type()
	if meta.Sparse {
		zeroFields := getZeroFields(typ, val)
		fancyNames := slices.Sorted(maps.Keys(zeroFields))
		newFieldSelector := widget.NewSelect(fancyNames, nil)
		parent.Append("", container.NewBorder(nil, nil, nil, widget.NewButton("Add field", func() {
			if newFieldSelector.Selected == "" {
				return
			}
			fi := zeroFields[newFieldSelector.Selected]
			tryMakeField(parent, fi, val, w, FieldMeta{Sparse: false}, false)
			delete(zeroFields, newFieldSelector.Selected)
			newFieldSelector.SetOptions(slices.Sorted(maps.Keys(zeroFields)))
			newFieldSelector.SetSelectedIndex(0)
			parent.Refresh()
		}), newFieldSelector))
	}
	for i := 0; i < typ.NumField(); i++ {
		tryMakeField(parent, i, val, w, meta, false)
	}
}

func tryMakeField(parent *widget.Form, fieldIndex int, val reflect.Value, w fyne.Window, meta FieldMeta, prepend bool) {
	typ := val.Type()
	sf := typ.Field(fieldIndex)
	if sf.IsExported() {
		fv := val.Field(fieldIndex)
		tag := sf.Tag.Get("form")
		fieldMeta := parseFormTag(tag)
		if fieldMeta.Hidden {
			return
		}
		label := fieldMeta.Label
		if label == "" {
			label = sf.Name
		}

		var entry fyne.CanvasObject
		switch fv.Kind() {
		case reflect.Struct:
			entry = widget.NewForm()
			buildNestedForm(entry.(*widget.Form), fv, w, fieldMeta)
		case reflect.Map:
			if fv.IsNil() {
				fv.Set(reflect.MakeMap(fv.Type()))
			}
			switch fieldMeta.Widget {
			case "floatmap":
				entry = buildEditableTable(fv, fieldMeta.Element, w)
			default:
				entry = buildMapAppTabs(fv, fieldMeta.Element, w, FieldMeta{Sparse: true})

			}
		default:
			if meta.Sparse && meta.Widget == "" && fv.IsZero() {
				return
			}
			switch fieldMeta.Widget {
			case "checkbox":
				entry = widget.NewCheck("", nil)
				entry.(*widget.Check).SetChecked(fv.Bool())
				entry.(*widget.Check).OnChanged = func(b bool) {
					fv.SetBool(b)
				}

			case "file":
				pathLabel := widget.NewLabel(fv.String())
				pathLabel.Truncation = fyne.TextTruncateOff
				pathLabel.Alignment = fyne.TextAlignLeading

				cancelBtn := widget.NewButton("Clear", func() {})
				cancelBtn.OnTapped = func() {
					pathLabel.SetText("")
					fv.SetString("")
					cancelBtn.Hide()
				}
				cancelBtn.Hide()
				pickBtn := widget.NewButton("Browse filesystem", func() {
					dialog.ShowFileOpen(func(reader fyne.URIReadCloser, err error) {
						if err != nil || reader == nil {
							return
						}
						cancelBtn.Show()
						path := reader.URI().Path()
						fv.SetString(path)
						pathLabel.SetText(path)
						reader.Close()
					}, w)
				})
				entry = container.NewBorder(nil, nil, nil, container.NewHBox(cancelBtn, pickBtn), container.NewHScroll(pathLabel))
			case "folder":
				pathLabel := widget.NewLabel(fv.String())
				pathLabel.Truncation = fyne.TextTruncateOff
				pathLabel.Alignment = fyne.TextAlignLeading

				cancelBtn := widget.NewButton("Clear", func() {})
				cancelBtn.OnTapped = func() {
					pathLabel.SetText("")
					fv.SetString("")
					cancelBtn.Hide()
				}
				cancelBtn.Hide()
				pickBtn := widget.NewButton("Browse filesystem", func() {
					dialog.ShowFolderOpen(func(folderPath fyne.ListableURI, err error) {
						if err != nil {
							dialog.ShowError(err, w)
						}
						if folderPath == nil {
							return
						}
						cancelBtn.Show()
						path := folderPath.Path()
						fv.SetString(path)
						pathLabel.SetText(path)
					}, w)
				})
				entry = container.NewBorder(nil, nil, nil, container.NewHBox(cancelBtn, pickBtn), container.NewHScroll(pathLabel))

			case "text", "":
				if len(fieldMeta.Options) > 0 || fv.Kind() == reflect.Bool {
					if fv.Kind() == reflect.Bool {
						fieldMeta.Options = []string{"false", "true"}
					}
					entry = widget.NewSelect(fieldMeta.Options, nil)
					current := fv.String()
					if current != "" && slices.Contains(fieldMeta.Options, current) {
						entry.(*widget.Select).SetSelected(current)
					} else {
						entry.(*widget.Select).SetSelected(fieldMeta.Options[0])
						makeStr2TypeCallback(fv)(fieldMeta.Options[0])
					}
					entry.(*widget.Select).OnChanged = makeStr2TypeCallback(fv)
				} else {
					entry = widget.NewEntry()
					entry.(*widget.Entry).SetText(valueToStr(fv))
					entry.(*widget.Entry).OnChanged = makeStr2TypeCallback(fv)
				}
			}
		}
		if prepend {
			parent.Items = append([]*widget.FormItem{{Text: label, Widget: entry}}, parent.Items...)
		} else {
			parent.Append(label, entry)
		}

	}
}
