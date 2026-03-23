package gui

import (
	"fmt"
	"io"
	"os"
	"reflect"
	"runtime/debug"

	"fyne.io/fyne/v2"
	"fyne.io/fyne/v2/container"
	"fyne.io/fyne/v2/dialog"
	"fyne.io/fyne/v2/storage"
	"fyne.io/fyne/v2/widget"
	"github.com/BurntSushi/toml"
	"github.com/wildstyl3r/stmc/internal/config"
	"github.com/wildstyl3r/stmc/internal/messages"
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
	logger         messages.Logger
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
		saveDialog := dialog.NewFileSave(func(writer fyne.URIWriteCloser, err error) {
			if err != nil || writer == nil {
				return
			}
			uri := writer.URI()
			targetExt := ".toml"

			if uri.Extension() != targetExt {

				writer.Close()
				storage.Delete(uri)

				newName := uri.Name() + targetExt
				parent, _ := storage.Parent(uri)
				newURI, _ := storage.Child(parent, newName)

				exists, _ := storage.Exists(newURI)
				if exists {
					// Ask for confirmation before overwriting
					dialog.ShowConfirm("Overwrite?",
						"A file named "+newName+" already exists. Overwrite?",
						func(confirm bool) {
							if confirm {
								// User said yes, proceed with writing
								writer, err = storage.Writer(newURI)

								if err != nil {
									dialog.ShowError(err, fm.window)
									return
								}

								defer writer.Close()

								data, _ := toml.Marshal(fm.data)
								writer.Write(data)
								fm.configFilename.SetText(writer.URI().Path())
								fm.confFlags.ConfigFileNamePointer = &fm.configFilename.Text
							}
						}, fm.window)
					return
				} else {
					writer, err = storage.Writer(newURI)
					if err != nil {
						dialog.ShowError(err, fm.window)
						return
					}
				}
			}

			defer writer.Close()

			data, _ := toml.Marshal(fm.data)
			writer.Write(data)
			fm.configFilename.SetText(writer.URI().Path())
			fm.confFlags.ConfigFileNamePointer = &fm.configFilename.Text
		}, fm.window)
		saveDialog.SetFilter(storage.NewExtensionFileFilter([]string{".toml"}))
		saveDialog.Show()
	})
}

func (fm *FormManager) runButton() fyne.CanvasObject {
	return widget.NewButton("Save and run", func() {
		if fm.configFilename.Text == "" {
			saveDialog := dialog.NewFileSave(func(writer fyne.URIWriteCloser, err error) {
				if err != nil || writer == nil {
					return
				}
				uri := writer.URI()
				targetExt := ".toml"

				if uri.Extension() != targetExt {
					writer.Close()
					storage.Delete(uri)

					newName := uri.Name() + targetExt
					parent, _ := storage.Parent(uri)
					newURI, _ := storage.Child(parent, newName)

					exists, _ := storage.Exists(newURI)
					if exists {
						// Ask for confirmation before overwriting
						dialog.ShowConfirm("Overwrite?",
							"A file named "+newName+" already exists. Overwrite?",
							func(confirm bool) {
								if confirm {
									// User said yes, proceed with writing
									writer, err = storage.Writer(newURI)

									if err != nil {
										dialog.ShowError(err, fm.window)
										return
									}

									defer writer.Close()

									data, _ := toml.Marshal(fm.data)
									writer.Write(data)
									fm.configFilename.SetText(writer.URI().Path())
									fm.confFlags.ConfigFileNamePointer = &fm.configFilename.Text

									go func() {
										defer func() {
											if r := recover(); r != nil {
												stackTrace := string(debug.Stack())

												fyne.Do(func() {
													dialog.ShowInformation("CRITICAL ERROR", fmt.Sprintf("%v", r), fm.window)
												})

												fmt.Fprintf(os.Stderr, "CRITICAL ERROR: %v\n%s", r, stackTrace)

												// 3. Keep the UI alive to show the error
												// showPanicDialog(r, stackTrace)
											}
										}()
										router.RunConsole(fm.confFlags, fm.dataFlags, fm.logger)
									}()
								}
							}, fm.window)
						return
					} else {
						writer, err = storage.Writer(newURI)
						if err != nil {
							dialog.ShowError(err, fm.window)
							return
						}
					}
				}

				defer writer.Close()

				data, _ := toml.Marshal(fm.data)
				writer.Write(data)
				fm.configFilename.SetText(writer.URI().Path())
				fm.confFlags.ConfigFileNamePointer = &fm.configFilename.Text

				go func() {
					defer func() {
						if r := recover(); r != nil {
							stackTrace := string(debug.Stack())

							fyne.Do(func() {
								dialog.ShowInformation("CRITICAL ERROR", fmt.Sprintf("%v", r), fm.window)
							})

							fmt.Fprintf(os.Stderr, "CRITICAL ERROR: %v\n%s", r, stackTrace)
						}
					}()
					router.RunConsole(fm.confFlags, fm.dataFlags, fm.logger)
				}()
			}, fm.window)

			saveDialog.SetFilter(storage.NewExtensionFileFilter([]string{".toml"}))
			saveDialog.Show()
		} else {
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
			fm.confFlags.ConfigFileNamePointer = &fm.configFilename.Text

			defer func() {
				if r := recover(); r != nil {
					stackTrace := string(debug.Stack())

					dialog.ShowInformation("CRITICAL ERROR", fmt.Sprintf("%v", r), fm.window)
					fmt.Fprintf(os.Stderr, "CRITICAL ERROR: %v\n%s", r, stackTrace)

					// 3. Keep the UI alive to show the error
					// showPanicDialog(r, stackTrace)
				}
			}()
			router.RunConsole(fm.confFlags, fm.dataFlags, fm.logger)
			fmt.Fprintf(os.Stderr, "AFTER ERROR")
		}

	})
}

func NewFormManager(data any, filename string, w fyne.Window, cf *config.Flags, df *output.DataFlags, a fyne.App) *FormManager {
	fm := &FormManager{
		data:           data,
		configFilename: widget.NewLabel(filename),
		window:         w,
		form:           widget.NewForm(),
		confFlags:      cf,
		dataFlags:      df,
		logger:         messages.NewGUILogger(w, a),
	}
	fm.Refresh()
	return fm
}
