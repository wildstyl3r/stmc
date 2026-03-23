package gui

import (
	"bufio"
	"fmt"
	"os"

	"fyne.io/fyne/v2"
	"fyne.io/fyne/v2/container"
	"fyne.io/fyne/v2/widget"
)

type widgetWriter struct {
	widget     *widget.TextGrid
	origStdout *os.File
}

func (w *widgetWriter) Write(p string) (n int, err error) {
	// Append the new data to the existing text
	// The widget methods handle UI updates safely
	// w.widget.SetText(w.widget.Text() + p)
	w.widget.Append(p)
	fmt.Fprintln(w.origStdout, p)
	return len(p), nil
}

func ShowConsoleWindow(a fyne.App) {
	consoleWindow := a.NewWindow("Console Output")
	consoleWindow.Resize(fyne.NewSize(600, 400))

	// Use a MultiLineEntry widget to display text. Disable it so users can't type in it.
	consoleOutput := widget.NewTextGrid()
	// consoleOutput.Disable()

	// Wrap it in a scroll container in case output exceeds window size
	scrollContainer := container.NewScroll(consoleOutput)

	consoleWindow.SetContent(container.NewBorder(nil, widget.NewButton("Copy", func() {
		a.Clipboard().SetContent(consoleOutput.Text())
	}), nil, nil, scrollContainer))

	writer := &widgetWriter{widget: consoleOutput, origStdout: os.Stdout}
	// You could also use io.MultiWriter to still print to the actual console
	// if running from the command line while developing.
	// mw := io.MultiWriter(os.Stdout, writer)

	// origStdout := os.Stdout
	r, w, _ := os.Pipe()

	os.Stdout = w
	os.Stderr = w

	consoleWindow.Show()
	go func() {
		sc := bufio.NewScanner(r)
		for sc.Scan() {
			text := sc.Text()
			fyne.DoAndWait(func() {
				writer.Write(text)
			})

		}
	}()

	// // Example of outputting something to the "console"
	// go func() {
	// 	for {
	// 		log.Println("Logging an update:", time.Now().Format(time.Stamp))
	// 		time.Sleep(2 * time.Second)
	// 	}
	// }()

	consoleWindow.Show()
}
