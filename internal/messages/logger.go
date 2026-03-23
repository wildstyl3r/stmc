package messages

type Logger interface {
	Info(format string, args ...any)
	Failure(format string, args ...any)
	Progress(descripiton string, max float64) (func(update string), func(increment float64), func())
	Busy(descripiton string) func()
}
