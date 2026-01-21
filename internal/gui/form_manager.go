package gui

import (
	"fmt"
	"io"
	"os"
	"reflect"

	"fyne.io/fyne/v2"
	"fyne.io/fyne/v2/container"
	"fyne.io/fyne/v2/dialog"
	"fyne.io/fyne/v2/widget"
	"github.com/BurntSushi/toml"
	"github.com/wildstyl3r/stmc/internal/config"
	"github.com/wildstyl3r/stmc/internal/output"
	"github.com/wildstyl3r/stmc/internal/router"
)

type FormManager struct {
	data           any
	configFilename *widget.Label
	window         fyne.Window
	form           *widget.Form
	confFlags      *config.Flags
	dataFlags      *output.DataFlags
}

func (fm *FormManager) SetDataContext(data any) {
	fm.data = data
}

func (fm *FormManager) Refresh() {
	fm.form.Items = nil
	rv := reflect.ValueOf(fm.data).Elem()
	buildNestedForm(fm.form, rv, fm.window, FieldMeta{})
	fm.form.Refresh()
}
func (fm *FormManager) GetContent() fyne.CanvasObject {
	configFile := container.NewBorder(
		nil, nil, container.NewHBox(fm.loadButton(), fm.saveButton()),
		fm.runButton(), container.NewScroll(fm.configFilename),
	)
	return container.NewVBox(configFile, fm.form)
}

func (fm *FormManager) loadButton() fyne.CanvasObject {
	return widget.NewButton("Load", func() {
		fmt.Println("load")
		dialog.ShowFileOpen(func(r fyne.URIReadCloser, err error) {
			if err != nil || r == nil {
				return
			}
			defer r.Close()

			data, _ := io.ReadAll(r)
			toml.Unmarshal(data, fm.data)
			fm.configFilename.SetText(r.URI().Path())

			fm.Refresh()
		}, fm.window)
	})
}

func (fm *FormManager) saveButton() fyne.CanvasObject {
	return widget.NewButton("Save as", func() {
		fmt.Println("save as")
		dialog.ShowFileSave(func(writer fyne.URIWriteCloser, err error) {
			if err != nil || writer == nil {
				return
			}
			defer writer.Close()

			data, _ := toml.Marshal(fm.data)
			writer.Write(data)
			fm.configFilename.SetText(writer.URI().Path())
		}, fm.window)
	})
}

func (fm *FormManager) runButton() fyne.CanvasObject {
	return widget.NewButton("Save and run", func() {
		fmt.Println("save&r0")
		if fm.configFilename.Text == "" {
			fmt.Println("save&r1")
			dialog.ShowFileSave(func(writer fyne.URIWriteCloser, err error) {
				if err != nil || writer == nil {
					return
				}
				defer writer.Close()

				data, _ := toml.Marshal(fm.data)
				writer.Write(data)
				fm.configFilename.SetText(writer.URI().Path())
			}, fm.window)
		} else {
			fmt.Println("save&r2")
			f, err := os.Create(fm.configFilename.Text)
			if err != nil {
				fmt.Printf("failed to create file: %v", err)
				return
			}

			if err := toml.NewEncoder(f).Encode(fm.data); err != nil {
				f.Close() // Ensure the file is closed even if encoding fails
				fmt.Printf("failed to encode to TOML: %v", err)
				return
			}

			if err := f.Close(); err != nil {
				fmt.Printf("failed to close file: %v", err)
				return
			}
		}
		fmt.Println("passing config to the runner")
		fm.confFlags.ConfigFileNamePointer = &fm.configFilename.Text
		router.RunConsole(fm.confFlags, fm.dataFlags)
	})
}

func NewFormManager(data any, filename string, w fyne.Window, cf *config.Flags, df *output.DataFlags) *FormManager {
	fm := &FormManager{
		data:           data,
		configFilename: widget.NewLabel(filename),
		window:         w,
		form:           widget.NewForm(),
		confFlags:      cf,
		dataFlags:      df,
	}
	fm.Refresh()
	return fm
}
