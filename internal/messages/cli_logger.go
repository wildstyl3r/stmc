package messages

import (
	"fmt"

	"github.com/schollz/progressbar/v3"
)

type CLILogger struct {
}

func NewCLILogger() Logger {
	return &CLILogger{}
}

func (cl *CLILogger) Info(format string, args ...any) {
	fmt.Printf(format, args...)
}

func (cl *CLILogger) Failure(format string, args ...any) {
	panic(fmt.Sprintf(format, args...))
}

func (cl *CLILogger) Progress(descripiton string, max float64) (func(update string), func(increment float64), func()) {
	pb := progressbar.Default(int64(max))
	return func(update string) {
			pb.Describe(update)
		}, func(increment float64) {
			pb.Add(int(increment))
		}, func() {
			pb.Finish()
		}
}

func (cl *CLILogger) Busy(description string) func() {
	pb := progressbar.Default(-1)
	return func() {
		pb.Finish()
	}
}
