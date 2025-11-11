package config

import (
	"slices"
	"testing"
)

func TestTokenizer(t *testing.T) {
	expressions := map[string]struct {
		input  string
		result []string
	}{
		"Trivial":        {"Quantity", []string{"Quantity"}},
		"Simple":         {"A*BCD", []string{"A", "*", "BCD"}},
		"Spaced":         {" Velocity *  mass ", []string{"Velocity", "*", "mass"}},
		"Power":          {"A^3", []string{"A", "^", "3"}},
		"Negative power": {"A^-2", []string{"A", "^", "-2"}},
	}
	for name, test := range expressions {
		t.Run(name, func(t *testing.T) {
			parts := tokenizer(test.input)
			if !slices.Equal(parts, test.result) {
				t.Errorf(`tokenizer(%s) = %#v, want %#v`, test.input, parts, test.result)
			}
		})
	}
}

func TestSanityCheck(t *testing.T) {
	tokenLists := map[string]struct {
		input  []string
		result error
	}{
		"Trivial":                                    {[]string{"A"}, nil},
		"Trivial failing at start":                   {[]string{"*", "A"}, SyntaxError(0)},
		"Trivial failing at end":                     {[]string{"B", "*"}, SyntaxError(1)},
		"Simple expression":                          {[]string{"A", "*", "B"}, nil},
		"Erroneous expression":                       {[]string{"A", "*", "*", "B"}, SyntaxError(2)},
		"Unfinished expression":                      {[]string{"A", "*"}, SyntaxError(1)},
		"Expression with power":                      {[]string{"A", "^", "2"}, nil},
		"Expression with too many powers":            {[]string{"A", "^", "2", "^", "2"}, SyntaxError(3)},
		"Expression with unsupported power":          {[]string{"A", "^", "B"}, SyntaxError(2)},
		"Another expression with unsupported power":  {[]string{"A", "^", "2.718"}, SyntaxError(2)},
		"Expression with negative power":             {[]string{"A", "^", "-2"}, nil},
		"Expression with spaced negative power":      {[]string{"A", "^", "-", "2"}, SyntaxError(2)},
		"Expression with unsupported negative power": {[]string{"A", "^", "-B"}, SyntaxError(2)},
	}
	for name, test := range tokenLists {
		t.Run(name, func(t *testing.T) {
			err := sanityCheck(test.input)
			if err != test.result {
				t.Errorf(`sanityCheck(%#v) = %#v, want %#v`, test.input, err, test.result)
			}
		})
	}
}

func TestDescribe(t *testing.T) {
	expressions := map[string]struct {
		input       string
		description product
		err         error
	}{
		"Trivial": {"a", product{exponent{"a", 1}}, nil},
		"Simple":  {"a*b^2", product{exponent{"a", 1}, exponent{"b", 2}}, nil},
		"Inverse": {"a^-3", product{exponent{"a", -3}}, nil},
	}
	for name, test := range expressions {
		t.Run(name, func(t *testing.T) {
			description, err := expressAsProduct(test.input)
			if err != test.err {
				t.Errorf(`sanityCheck(%#v) = %#v, want %#v`, test.input, err, test.err)
				return
			}
			if !slices.Equal(description, test.description) {
				t.Errorf(`expressAsProduct(%#v) = %#v, want %#v`, test.input, description, test.description)
			}
		})
	}
}

func TestInverseAndConstructName(t *testing.T) {
	expressions := map[string]struct {
		input  string
		output string
		err    error
	}{
		"Trivial":         {"a", "a^-1", nil},
		"Trivial inverse": {"a^-1", "a", nil},
		"Quadratic 1":     {"a^-2", "a^2", nil},
		"Quadratic 2":     {"a^2", "a^-2", nil},
		"Lengthy":         {"a* c^2 *b*d ^ -3*f ^-1", "a^-1*b^-1*c^-2*d^3*f", nil},
	}
	for name, test := range expressions {
		t.Run(name, func(t *testing.T) {
			product, err := expressAsProduct(test.input)
			if err != test.err {
				t.Errorf(`sanityCheck(%#v) = %#v, want %#v`, test.input, err, test.err)
				return
			}
			inverse := product.Inverse().ConstructName()
			if inverse != test.output {
				t.Errorf(`(%#v).Inverse() is %#v, want %#v`, test.input, inverse, test.output)
			}
		})
	}
}
