package gui

import (
	"reflect"
	"strconv"
)

func getZeroFields(rt reflect.Type, fvStruct reflect.Value) (r map[string]int) {
	r = make(map[string]int)
	for i := 0; i < rt.NumField(); i++ {
		field := rt.Field(i)
		meta := parseFormTag(field.Tag.Get("form"))
		if field.IsExported() && !meta.Hidden {
			fv := fvStruct.Field(i)
			if fv.IsZero() {
				name := field.Name
				if meta.Label != "" {
					name = meta.Label
				}
				r[name] = i
			}
		}
	}
	return r
}

func setNonZeroDefault(fv reflect.Value) {
	switch fv.Kind() {
	case reflect.String:
		fv.SetString("default")
	case reflect.Int:
		fv.SetInt(1)
	case reflect.Bool:
		fv.SetBool(true)
		// Add more defaults as needed
	}
}

func makeStr2TypeCallback(val reflect.Value) func(string) {
	switch val.Kind() {
	case reflect.Int, reflect.Int64:
		return func(s string) {
			if n, err := strconv.Atoi(s); err == nil {
				val.SetInt(int64(n))
			}
		}
	case reflect.Float64:
		return func(s string) {
			if f, err := strconv.ParseFloat(s, 64); err == nil {
				val.SetFloat(f)
			}
		}
	case reflect.Bool:
		return func(s string) {
			if f, err := strconv.ParseBool(s); err == nil {
				val.SetBool(f)
			}
		}
	case reflect.String:
		return func(s string) {
			val.SetString(s)
		}
	default:
		return func(s string) {}
	}
}

func valueToStr(val reflect.Value) string {
	switch val.Kind() {
	case reflect.Int, reflect.Int64:
		return strconv.Itoa(int(val.Int()))
	case reflect.Float64:
		return strconv.FormatFloat(val.Float(), 'f', 2, 64)
	case reflect.Bool:
		return strconv.FormatBool(val.Bool())
	// case reflect.String:
	default:
		return val.String()
	}
}
