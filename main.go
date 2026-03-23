package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"runtime"
	"runtime/pprof"
	"time"

	"fyne.io/fyne/v2"
	"fyne.io/fyne/v2/app"
	"fyne.io/fyne/v2/container"
	"fyne.io/fyne/v2/widget"
	"github.com/wildstyl3r/stmc/internal/config"
	"github.com/wildstyl3r/stmc/internal/gui"
	"github.com/wildstyl3r/stmc/internal/messages"
	"github.com/wildstyl3r/stmc/internal/output"
	"github.com/wildstyl3r/stmc/internal/router"
)

func main() {
	console := flag.Bool("no-gui", false, "run as console app (no GUI)")
	startTime := time.Now()
	configFlags := config.NewConfigFlags()
	dataExtractorFlags := output.NewDataFlags()
	flag.Parse()

	if *configFlags.CPUprofile != "" {
		f, err := os.Create(*configFlags.CPUprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	if *console {
		router.RunConsole(&configFlags, &dataExtractorFlags, messages.NewCLILogger())
		return
	}
	runGUI(&configFlags, &dataExtractorFlags)

	fmt.Printf("Total elapsed time: %v\n\n", time.Since(startTime))
	if *configFlags.MEMprofile != "" {
		f, err := os.Create(*configFlags.MEMprofile)
		if err != nil {
			log.Fatal("could not create memory profile: ", err)
		}
		defer f.Close() // error handling omitted for example
		runtime.GC()    // get up-to-date statistics
		// Lookup("allocs") creates a profile similar to go test -memprofile.
		// Alternatively, use Lookup("heap") for a profile
		// that has inuse_space as the default index.
		if err := pprof.Lookup("allocs").WriteTo(f, 0); err != nil {
			log.Fatal("could not write memory profile: ", err)
		}
	}
}

func runGUI(cf *config.Flags, df *output.DataFlags) {
	myApp := app.NewWithID("stmc gui")
	w := myApp.NewWindow("Monte carlo glow discharge model")
	w.Resize(fyne.NewSize(800, 600))
	ac := &config.Config{}
	fm := gui.NewFormManager(ac, "", w, cf, df, myApp)
	gui.ShowConsoleWindow(myApp)

	w.SetContent(container.NewBorder(widget.NewAccordion(widget.NewAccordionItem("Data to save", gui.CheckboxizeDataFlags(df))), nil, nil, nil, container.NewScroll(fm.GetContent())))
	w.ShowAndRun()
}
