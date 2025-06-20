package datahub

import "github.com/wildstyl3r/stmc/internal/model"

type KeyType string
type DataHubType map[KeyType]any
type Extension func(*model.Model) ([]string, []any)

var dataHubExtensions map[KeyType]Extension
var modelPtr *model.Model
var dataHub DataHubType

func (dh DataHubType) calculate(name KeyType, model *model.Model) {
	if _, exists := dataHub[name]; !exists {
		if calc, possible := dataHubExtensions[name]; !possible {
			panic("unimplemented")
		} else {
			names, values := calc(model)
			for i := range names {
				dataHub[KeyType(names[i])] = values[i]
			}
		}
	}
}

func Insert(name KeyType, value any) {
	dataHub[name] = value
}

func Get(name KeyType, model *model.Model) any {
	if modelPtr != model {
		panic("wrong model object for data request")
	}
	if value, exists := dataHub[name]; exists {
		return value
	} else {
		dataHub.calculate(name, model)
		return dataHub[name]
	}
}

func Init() {
	dataHubExtensions = make(map[KeyType]Extension)
}

func Reset(mPtr *model.Model) {
	modelPtr = mPtr
	dataHub = make(DataHubType)
}

func Register(name string, f Extension) {
	dataHubExtensions[KeyType(name)] = f
}
