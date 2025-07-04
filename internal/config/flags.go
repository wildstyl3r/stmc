package config

import (
	"flag"
	"runtime"
)

type Flags struct {
	CPUprofile, MEMprofile, ConfigFileNamePointer *string
	Threads                                       *int
	Verbose                                       *bool
	Debug                                         *bool
	RootFindingAlgorithm                          *string
	// dataFlags
}

func NewConfigFlags() Flags {
	return Flags{
		CPUprofile:            flag.String("cpuprofile", "", "write cpu profile to file"),
		MEMprofile:            flag.String("memprofile", "", "write memory profile to `file`"),
		ConfigFileNamePointer: flag.String("i", "inputs/val/BM_He.toml", "model configuration in toml format"),
		Threads:               flag.Int("j", runtime.NumCPU(), "threads to run"),
		Verbose:               flag.Bool("v", true, "verbose"),
		Debug:                 flag.Bool("d", false, "debug"),
		RootFindingAlgorithm:  flag.String("alg", "s", "root finder algorithm for gamma calculation ([s]tochastic approximation, [b]inary search, [t]ernary search)"),
	}
}
