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
	"github.com/wildstyl3r/stmc/internal/config"
	"github.com/wildstyl3r/stmc/internal/gui"
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
		router.RunConsole(&configFlags, &dataExtractorFlags)
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

type NestedStruct struct {
	Level  string  `toml:"level" form:"label=Level"`
	Value  float64 `toml:"value" form:"label=Value"`
	Active bool    `toml:"active" form:"label=active,widget=checkbox"`
}

// Rigid TOML structure
type AppConfig struct {
	Name        string                   `toml:"name" form:"label=Name,widget=text"`
	Description string                   `toml:"description" form:"label=Desc,widget=text"`
	Count       int                      `toml:"count" form:"label=Count"`
	Enabled     bool                     `toml:"enabled" form:"label=Enabled,widget=checkbox"`
	LogLevel    string                   `toml:"log_level" form:"label=Log level,options=debug|info|warn|error"` // Dropdown: info, warn, error
	Nested      NestedStruct             `toml:"nested" form:"label=Nested"`
	LogFile     string                   `form:"label=Log File,widget=path" toml:"log_file"`
	Mapped      map[string]*NestedStruct `form:"label=Servers"`
}

func runGUI(cf *config.Flags, df *output.DataFlags) {
	myApp := app.NewWithID("stmc gui")
	w := myApp.NewWindow("Monte carlo glow discharge model")
	w.Resize(fyne.NewSize(800, 600))
	ac := &config.Config{}
	fm := gui.NewFormManager(ac, "", w, cf, df)
	gui.ShowConsoleWindow(myApp)

	w.SetContent(container.NewScroll(fm.GetContent()))
	w.ShowAndRun()
}
