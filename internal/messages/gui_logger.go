package messages

import (
	"fmt"

	"fyne.io/fyne/v2"
	"fyne.io/fyne/v2/container"
	"fyne.io/fyne/v2/dialog"
	"fyne.io/fyne/v2/widget"
)

type GUILogger struct {
	window fyne.Window
	app    fyne.App
}

func NewGUILogger(w fyne.Window, a fyne.App) Logger {
	return &GUILogger{
		window: w,
		app:    a,
	}
}

func (gl *GUILogger) Info(format string, args ...any) {
	info := fmt.Sprintf(format, args...)
	dialog.ShowInformation("Status", info, gl.window)
	// dialog.ShowCustomWithoutButtons("Status", container.NewBorder(nil, widget.NewButton("Copy", func() {
	// 	gl.app.Clipboard().SetContent(info)
	// }), nil, nil, container.NewHScroll(widget.NewLabel(info))), gl.window)
	fmt.Printf(format, args...)
}

func (gl *GUILogger) Failure(format string, args ...any) {
	info := fmt.Sprintf(format, args...)
	dialog.ShowInformation("Error", info, gl.window)
	// dialog.ShowCustom("Error", container.NewBorder(nil, widget.NewButton("Copy", func() {
	// 	gl.app.Clipboard().SetContent(info)
	// }), nil, nil, container.NewHScroll(widget.NewLabel(info))), gl.window)
	// panic(fmt.Sprintf(format, args...))
}

func (gl *GUILogger) Progress(descripiton string, max float64) (func(update string), func(increment float64), func()) {
	label := widget.NewLabel("")
	pb := widget.NewProgressBar()
	d := dialog.NewCustomWithoutButtons(descripiton, container.NewVBox(label, pb), gl.window)
	d.Show()
	return func(update string) {
			label.SetText(update)
		}, func(increment float64) {
			pb.SetValue(pb.Value + increment/max)
		}, func() {
			d.Dismiss()
		}
}

func (gl *GUILogger) Busy(descripiton string) func() {
	busy := widget.NewProgressBarInfinite()
	d := dialog.NewCustomWithoutButtons(descripiton, busy, gl.window)
	d.Show()
	return func() {
		busy.Stop()
		d.Dismiss()
	}
}
