package utils

import "reflect"

func GetStructFieldNamesAndAliases(s any) (fieldNames, csvAliases []string) {
	targetReflect := reflect.ValueOf(s)
	targetType := targetReflect.Type()
	for i := 0; i < targetReflect.NumField(); i++ {
		field := targetType.Field(i)
		value := targetReflect.Field(i)
		if field.Anonymous && field.Type.Kind() == reflect.Struct {
			embeddedValue := value
			embeddedType := embeddedValue.Type()
			for j := 0; j < embeddedValue.NumField(); j++ {
				embeddedField := embeddedType.Field(j)
				fieldNames = append(fieldNames, embeddedField.Name)
				csvAliases = append(csvAliases, embeddedField.Tag.Get("csv"))
			}
		} else {
			fieldNames = append(fieldNames, field.Name)
			csvAliases = append(csvAliases, field.Tag.Get("csv"))
		}
	}
	return fieldNames, csvAliases
}
