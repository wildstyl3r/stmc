package gui

import (
	"maps"
	"reflect"
	"slices"

	"fyne.io/fyne/v2"
	"fyne.io/fyne/v2/container"
	"fyne.io/fyne/v2/dialog"
	"fyne.io/fyne/v2/layout"
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
			label, entry := tryMakeField(fi, val, w, FieldMeta{Sparse: false}, true)
			parent.Append(label, entry)
			delete(zeroFields, newFieldSelector.Selected)
			newFieldSelector.SetOptions(slices.Sorted(maps.Keys(zeroFields)))
			newFieldSelector.SetSelectedIndex(0)
			parent.Refresh()
		}), newFieldSelector))
	}
	for i := 0; i < typ.NumField(); i++ {
		label, entry := tryMakeField(i, val, w, meta, false)
		if entry != nil {
			parent.Append(label, entry)
		}
	}
}

func tryMakeField(fieldIndex int, val reflect.Value, w fyne.Window, meta FieldMeta, force bool) (string, fyne.CanvasObject) {
	typ := val.Type()
	sf := typ.Field(fieldIndex)
	if sf.IsExported() {
		fv := val.Field(fieldIndex)
		tag := sf.Tag.Get("form")
		fieldMeta := parseFormTag(tag)
		if fieldMeta.Hidden {
			return "", nil
		}
		label := fieldMeta.Label
		if label == "" {
			label = sf.Name
		}

		var entry fyne.CanvasObject
		switch fv.Kind() {
		case reflect.Struct:
			if fieldMeta.Widget == "row" {
				var rowElems []fyne.CanvasObject
				v := reflect.Indirect(fv)
				typ := v.Type()
				// fmt.Printf("typ is %v, nf is %v\n", typ, typ.NumField())
				for i := 0; i < typ.NumField(); i++ {
					l, el := tryMakeField(i, v, w, meta, false)
					if el != nil {
						rowElems = append(rowElems, container.NewHBox(widget.NewLabel(l), el))
					}

				}
				entry = container.New(layout.NewRowWrapLayout(), rowElems...)
			} else {
				entry = widget.NewForm()
				buildNestedForm(entry.(*widget.Form), fv, w, fieldMeta)
				entry = container.NewVBox(entry, widget.NewSeparator())
			}
		case reflect.Map:
			if fv.IsNil() {
				fv.Set(reflect.MakeMap(fv.Type()))
			}
			if meta.Sparse && fv.Len() == 0 && !force {
				return "", nil
			}
			switch fieldMeta.Widget {
			case "floatmap":
				entry = buildEditableTable(fv, fieldMeta.Element, w)
			default:
				entry = container.NewVBox(buildMapAppTabs(fv, fieldMeta.Element, w, FieldMeta{Sparse: true}), widget.NewSeparator())

			}
		default:
			if meta.Sparse && fv.IsZero() && !force {
				return "", nil
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
		return label, entry

	}
	return "", nil
}
