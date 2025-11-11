package config

import (
	"fmt"
	"math"
	"slices"
	"strconv"
	"strings"

	"github.com/wildstyl3r/stmc/internal/utils"
)

func tokenizer(s string) (t []string) {
	k := -1
	for i := range s {
		if k == -1 && s[i] != ' ' {
			k = i
		}
		if k != -1 {
			if s[i] == ' ' || s[i] == '*' || s[i] == '^' {
				if k != i {
					t = append(t, s[k:i])
				}
				if s[i] == '*' || s[i] == '^' {
					t = append(t, s[i:i+1])
				}
				k = i + 1
			}
		}
	}
	if k != -1 && k != len(s) {
		t = append(t, s[k:])
	}
	return
}

type operation int

const (
	Nothing operation = iota
	Multiplication
	Exponentiation
)

type SyntaxError int

func (se SyntaxError) Error() string {
	return fmt.Sprintf("error during expression check at index %d", int(se))
}

func sanityCheck(t []string) (err error) {
	operatorTokens := map[string]operation{"*": Multiplication, "^": Exponentiation}
	_, firstIsOperator := operatorTokens[t[0]]
	_, lastIsOperator := operatorTokens[t[len(t)-1]]
	if firstIsOperator {
		return SyntaxError(0)
	}
	if lastIsOperator {
		return SyntaxError(len(t) - 1)
	}
	lastOperation := Multiplication
	for i := range len(t) - 1 {
		currentOperation, currentIsOperator := operatorTokens[t[i]]
		_, nextIsOperator := operatorTokens[t[i+1]]
		if (currentIsOperator && nextIsOperator) || (!currentIsOperator && !nextIsOperator) {
			return SyntaxError(i + 1)
		}
		if currentOperation == Exponentiation {
			if lastOperation == Exponentiation {
				return SyntaxError(i)
			}
			if _, err := strconv.Atoi(t[i+1]); err != nil {
				return SyntaxError(i + 1)
			}
		}
		if currentIsOperator {
			lastOperation = currentOperation
		}
	}
	return nil
}

type exponent struct {
	base  string
	power int
}
type product []exponent

func describe(t []string) (p product, err error) {
	for i := 0; i < len(t); i += 2 {
		if len(t) > i+1 && t[i+1] == "^" {
			power, err := strconv.Atoi(t[i+2])
			if err != nil {
				return nil, SyntaxError(i + 2)
			}
			p = append(p, exponent{t[i], power})
			i += 2
		} else {
			p = append(p, exponent{t[i], 1})
		}
	}
	return p, nil
}

func expressAsProduct(s string) (p product, err error) {
	tokens := tokenizer(s)
	err = sanityCheck(tokens)
	if err != nil {
		return nil, err
	}
	return describe(tokens)
}

type value struct {
	description      product
	value            float64
	_ready, constant bool
	compositions     []*_expression
	inverse          *value
}

func (v value) reachableExpressions() (reachable []*_expression) {
	for _, e := range v.compositions {
		if e.CountUnknowns() == 1 {
			reachable = append(reachable, e)
		}
	}
	return reachable
}

type _expression struct {
	operationType                     operation
	leftOperand, rightOperand, result *value
}

func (e *_expression) CountUnknowns() (n int) {
	if !e.leftOperand._ready {
		n++
	}
	if !e.rightOperand._ready {
		n++
	}
	if !e.result._ready {
		n++
	}
	return n
}

func (e *_expression) calculate() (reachable []*_expression) {
	switch e.operationType {
	case Multiplication:
		if e.leftOperand._ready && e.rightOperand._ready {
			e.result.value = e.leftOperand.value * e.rightOperand.value
			e.result._ready = true
			return e.result.reachableExpressions()
		} else if e.leftOperand._ready && e.result._ready {
			e.rightOperand.value = e.result.value / e.leftOperand.value
			e.rightOperand._ready = true
			return e.rightOperand.reachableExpressions()
		} else if e.rightOperand._ready && e.result._ready {
			e.leftOperand.value = e.result.value / e.rightOperand.value
			e.leftOperand._ready = true
			return e.leftOperand.reachableExpressions()
		}
	case Exponentiation:
		if e.leftOperand._ready {
			e.result.value = math.Pow(e.leftOperand.value, e.rightOperand.value)
			e.result._ready = true
			return e.result.reachableExpressions()
		} else if e.result._ready {
			e.leftOperand.value = math.Pow(e.result.value, 1./e.rightOperand.value)
			e.leftOperand._ready = true
			return e.leftOperand.reachableExpressions()
		}
	}
	return reachable
}

func (p product) ConstructName() (r string) {
	productFactors := []string{}
	for factor := range p {
		if p[factor].power == 1 {
			productFactors = append(productFactors, p[factor].base)
		} else {
			productFactors = append(productFactors, p[factor].base+"^"+fmt.Sprintf("%d", p[factor].power))
		}
	}
	slices.Sort(productFactors)
	return strings.Join(productFactors, "*")
}

func (e exponent) Inverse() (r exponent) {
	return exponent{
		e.base,
		-e.power,
	}
}

func (p product) Inverse() (q product) {
	for factor := range p {
		q = append(q, p[factor].Inverse())
	}
	return q
}

type SymbolStorage struct {
	values      map[string]*value
	expressions []*_expression
}

func (st *SymbolStorage) insertConstant(name string, v float64) (ptr *value, err error) {
	if ptr, existed := st.values[name]; existed {
		if ptr.value != v {
			return nil,
				ConstantInsertionError{
					name:             name,
					oldValue:         ptr.value,
					conflictingValue: v,
				}
		}
		return ptr, nil
	}
	ptr = &value{
		description: product{exponent{name, 1}},
		_ready:      true,
		constant:    true,
		value:       v,
	}
	st.values[name] = ptr
	return ptr, nil
}

func NewSymbolStorage() SymbolStorage {
	st := SymbolStorage{
		values:      make(map[string]*value),
		expressions: make([]*_expression, 0),
	}
	st.insertConstant("1", 1)
	st.values["1"].inverse = st.values["1"]
	return st
}

type SymbolStorageInsertionError int

func (err SymbolStorageInsertionError) Error() string {
	return fmt.Sprintf("Too long product (%d)", int(err))
}

func (st *SymbolStorage) getOrInsert(p product) (ptr *value, inserted bool) {
	if ptr, existed := st.values[p.ConstructName()]; existed {
		return ptr, false
	}
	ptr = &value{
		description: p,
	}
	st.values[ptr.description.ConstructName()] = ptr
	return ptr, true
}

type ConstantInsertionError struct {
	name                       string
	oldValue, conflictingValue float64
}

func (ie ConstantInsertionError) Error() string {
	return fmt.Sprintf("Insertion error: redefinition of %#v as %f, previous definition was %f", ie.name, ie.conflictingValue, ie.oldValue)
}

func (st *SymbolStorage) createExpression(op operation, leftOperand, rightOperand, result *value) {
	exprPtr := &_expression{
		operationType: op,
		leftOperand:   leftOperand,
		rightOperand:  rightOperand,
		result:        result,
	}
	leftOperand.compositions = append(leftOperand.compositions, exprPtr)
	rightOperand.compositions = append(rightOperand.compositions, exprPtr)
	result.compositions = append(result.compositions, exprPtr)
}

func (st *SymbolStorage) insertWithInverse(p product) (ptr, inversePtr *value) {
	ptr, inserted := st.getOrInsert(p)
	inversePtr, inverseInserted := st.getOrInsert(p.Inverse())
	if inserted || inverseInserted {
		st.createExpression(Multiplication, ptr, inversePtr, st.values["1"])
		ptr.inverse = inversePtr
		inversePtr.inverse = ptr
	}
	return ptr, inversePtr
}

func (st *SymbolStorage) insertExponentWithBase(e exponent) error {
	power := max(e.power, -e.power)
	ptr, inversePtr := st.insertWithInverse(product{e})
	if power > 1 {
		basePtr, baseInversePtr := st.insertWithInverse(product{exponent{base: e.base, power: 1}})
		powerPtr, err := st.insertConstant(fmt.Sprintf("%d", power), float64(power))
		if err != nil {
			return err
		}
		st.createExpression(Exponentiation, basePtr, powerPtr, ptr)
		st.createExpression(Exponentiation, baseInversePtr, powerPtr, inversePtr)
	}
	return nil
}

func (st *SymbolStorage) Insert(description product) (err error) {
	if len(description) > 8 {
		return fmt.Errorf("too many factors in a product: %#v", description)
	}
	for i := range description {
		if err = st.insertExponentWithBase(description[i]); err != nil {
			return err
		}
	}

	st.insertWithInverse(description)
	upperLimit := min(256, 1<<len(description))
	for i := uint16(1); i < uint16(upperLimit); i++ {
		left, right, err := utils.PartitionByMask(description, i)
		if err != nil {
			return err
		}
		st.insertWithInverse(left)
		st.insertWithInverse(right)
	}
	return nil
}

func fuse(p, q product) (r product) {
	factorPowers := map[string]int{}
	for i := range p {
		factorPowers[p[i].base] += p[i].power
	}
	for i := range q {
		factorPowers[q[i].base] += q[i].power
	}
	for base, power := range factorPowers {
		if power != 0 {
			r = append(r, exponent{base, power})
		}
	}
	if len(r) == 0 {
		r = append(r, exponent{"1", 1})
	}
	return r
}

func (st *SymbolStorage) CrossConnectEverything() {
	pairs := make(map[struct{ first, second *value }]struct{}, len(st.values)*(len(st.values)-1)/2)
	for _, first := range st.values {
		for _, second := range st.values {
			if first != second && first.inverse != second {
				pairs[struct{ first, second *value }{first, second}] = struct{}{}
			}
		}
	}
	for pair := range pairs {
		fused := fuse(pair.first.description, pair.second.description)
		fusedName := fused.ConstructName()
		if _, exists := st.values[fusedName]; exists {
			st.createExpression(Multiplication, pair.first, pair.second, st.values[fusedName])
			st.createExpression(Multiplication, pair.first.inverse, pair.second.inverse, st.values[fusedName].inverse)
		}
	}
}

func (st *SymbolStorage) ClearValues() {
	for _, ptr := range st.values {
		if !ptr.constant {
			ptr._ready = false
			ptr.value = 0
		}
	}
}

type namedCell struct {
	name  string
	value float64
}

func (st *SymbolStorage) FillIfNotReady(row []namedCell) (used []namedCell) {
	for i := range row {
		if !st.values[row[i].name]._ready {
			st.values[row[i].name].value = row[i].value
			st.values[row[i].name]._ready = true
			used = append(used, row[i])
		}
	}
	return
}

func (st *SymbolStorage) Calculate() {
	backlog := make([]*_expression, 0)
	for _, expressionPtr := range st.expressions {
		if expressionPtr.CountUnknowns() == 1 {
			backlog = append(backlog, expressionPtr)
		}
	}
	for len(backlog) > 0 {
		lastIndex := len(backlog) - 1
		backlog = append(backlog[:lastIndex], backlog[lastIndex].calculate()...)
	}
}

type DataLoader struct {
}
