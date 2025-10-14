package model

type DataHubKeyType string
type DataMapType map[DataHubKeyType]any
type Extension func(*Model) ([]string, []any, error)
type DataHubType struct {
	data       DataMapType
	Extensions map[DataHubKeyType]Extension
}

// var dataHub DataHubType

func (model *Model) calculate(name DataHubKeyType) error {
	if _, exists := model.DataHub.data[name]; !exists {
		if calc, possible := model.DataHub.Extensions[name]; !possible {
			panic("unimplemented")
		} else {
			names, values, err := calc(model)
			if err != nil {
				return err
			}
			for i := range names {
				model.DataHub.data[DataHubKeyType(names[i])] = values[i]
			}
		}
	}
	return nil
}

func (dh *DataHubType) Insert(name DataHubKeyType, value any) {
	dh.data[name] = value
}

func (model *Model) GetMetrics(name DataHubKeyType) any {
	if value, exists := model.DataHub.data[name]; exists {
		return value
	} else {
		model.calculate(name)
		return model.DataHub.data[name]
	}
}

func (model *Model) ApplyToMetrics(name DataHubKeyType, f func(any) any) {
	if value, exists := model.DataHub.data[name]; exists {
		model.DataHub.data[name] = f(value)
	} else {
		model.calculate(name)
		model.DataHub.data[name] = f(value)
	}
}

func NewDataHub() *DataHubType {
	return &DataHubType{
		data:       make(DataMapType),
		Extensions: make(map[DataHubKeyType]Extension),
	}
}

func (dh DataHubType) Register(name string, f Extension) {
	dh.Extensions[DataHubKeyType(name)] = f
}
