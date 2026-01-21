package gui

import "strings"

type FieldMeta struct {
	Label      string
	Element    string
	Horizontal bool
	Widget     string
	Options    []string
	Hidden     bool
	Sparse     bool
}

func parseFormTag(tag string) FieldMeta {
	meta := FieldMeta{}
	if tag == "" {
		return meta
	}
	parts := strings.SplitSeq(tag, ",")
	for p := range parts {
		if p == "hide" {
			meta.Hidden = true
			return meta
		}
		if p == "sparse" {
			meta.Sparse = true
			return meta
		}
		kv := strings.SplitN(p, "=", 2)
		if len(kv) != 2 {
			continue
		}
		switch kv[0] {
		case "label":
			meta.Label = kv[1]
		case "widget":
			meta.Widget = kv[1]
		case "options":
			if kv[1] != "" {
				meta.Options = strings.Split(kv[1], "|")
			}
		case "element":
			meta.Element = kv[1]
		case "orientation":
			if kv[1] == "horiz" {
				meta.Horizontal = true
			}
		}
	}
	return meta
}
