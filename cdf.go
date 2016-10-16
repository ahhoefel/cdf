package cdf

import (
	"fmt"
	"image"
	"image/color"
	"math"
	"math/rand"
	"sort"
)

const (
	basePath = "/Users/hoefel/Development/go/src/github.com/ahhoefel/cdf"
)

type axis bool
type side bool

const (
	xAxis   axis = true
	yAxis   axis = false
	minSide side = true
	maxSide side = false
)

const (
	depth      = 7
	numPoints  = 10000
	showPoints = true
	showBorder = true
	xPixels    = 800
	yPixels    = 800
)

var (
	axisName = map[axis]string{xAxis: "xAxis", yAxis: "yAxis"}
)

const (
	debugPlot = true
)

type Pt struct {
	X, Y float64
}

func (p Pt) get(ax axis) float64 {
	if ax == xAxis {
		return p.X
	}
	return p.Y
}

func randomPoints(n int, scale float64) []Pt {
	var pts []Pt
	for i := 0; i < n; i++ {
		pts = append(pts, Pt{scale * rand.ExpFloat64(), scale * rand.Float64()})
	}
	return pts
}

func randomPointsCorr(n int, scale float64) []Pt {
	var pts []Pt
	for i := 0; i < n; i++ {
		a := rand.ExpFloat64()
		b := rand.Float64()
		pts = append(pts, Pt{scale * (a + b), scale * (a - b)})
	}
	return pts
}

func RandomPointsNorm(n int, scale float64, shift float64) []Pt {
	var pts []Pt
	for i := 0; i < n; i++ {
		a := rand.NormFloat64()*scale + shift
		b := rand.NormFloat64()*scale + shift
		pts = append(pts, Pt{a, b})
	}
	return pts
}

func randomPointsTwoNorm(n int, scale float64, shiftX float64, shiftY float64) []Pt {
	var pts []Pt
	for i := 0; i < n; i++ {
		a := rand.NormFloat64()*scale + shiftX
		b := rand.NormFloat64()*scale + shiftY
		pts = append(pts, Pt{a, b})
	}
	for i := 0; i < n; i++ {
		a := rand.NormFloat64()*scale + shiftX
		b := rand.NormFloat64()*scale + shiftY + shiftY
		pts = append(pts, Pt{a, b})
	}
	return pts
}

type pivot struct {
	Min float64
	Max float64
	med float64
	ax  axis
}

type quad []pivot

type division struct {
	p     pivot
	lower []Pt
	upper []Pt
}

type pend struct {
	pts   []Pt
	ax    axis
	depth int
}

func newQuad(depth int, pts []Pt, ax axis) *quad {
	queue := []pend{{pts, ax, depth}}
	var q quad
	for i := 0; i < len(queue); i++ {
		d := divide(queue[i].pts, queue[i].ax)
		q = append(q, d.p)
		if queue[i].depth > 0 {
			queue = append(queue, pend{d.lower, !queue[i].ax, queue[i].depth - 1})
			queue = append(queue, pend{d.upper, !queue[i].ax, queue[i].depth - 1})
		}
	}
	return &q
}

func divide(pts []Pt, ax axis) division {
	mid := midpoint(pts)
	low := 0
	high := len(pts)
	i := partition(pts, low, high, ax)
	for i != mid {
		if i < mid {
			low = i + 1
			i = partition(pts, low, high, ax)
		} else {
			high = i
			i = partition(pts, low, high, ax)
		}
	}
	var d division
	d.p.Min = ptsMin(pts, ax)
	d.p.Max = ptsMax(pts, ax)
	d.p.med = pts[mid].get(ax)
	d.p.ax = ax
	d.lower = pts[:mid]
	d.upper = pts[mid:]
	return d
}

var midpointToggle int = 0

func midpoint(pts []Pt) int {
	if len(pts)%2 == 0 {
		if midpointToggle == 0 {
			midpointToggle = 1
		} else {
			midpointToggle = 0
		}
		return len(pts)/2 - midpointToggle
	}
	return len(pts) / 2
}

func ptsMin(pts []Pt, ax axis) float64 {
	if len(pts) == 0 {
		return 0
	}
	min := pts[0].get(ax)
	for _, p := range pts {
		if z := p.get(ax); z < min {
			min = z
		}
	}
	return min
}

func ptsMax(pts []Pt, ax axis) float64 {
	if len(pts) == 0 {
		return 0
	}
	min := pts[0].get(ax)
	for _, p := range pts {
		if z := p.get(ax); z > min {
			min = z
		}
	}
	return min
}

// Uses the first element as a pivot for a quicksort partition.
// The result is everything < the pivot on the left of the pivot
// and everything >= pivot are on the right. The index of the pivot
// is returned.
func partition(pts []Pt, low, high int, ax axis) int {
	z := pts[low].get(ax)
	j := low + 1
	k := low + 1
	for k < high {
		if pts[k].get(ax) < z {
			tmp := pts[k]
			pts[k] = pts[j]
			pts[j] = tmp
			j++
		}
		k++
	}
	tmp := pts[j-1]
	pts[j-1] = pts[low]
	pts[low] = tmp
	return j - 1
}

func plotQuad(img *image.RGBA, q *quad, depth int) {
	min := img.Bounds().Min.Y
	max := img.Bounds().Max.Y
	recPlotQuad(img, q, 0, min, max, xAxis, depth)
}

func recPlotQuad(img *image.RGBA, q *quad, index, min, max int, ax axis, depth int) {
	if index >= len(*q) {
		return
	}
	if depth == 0 {
		area := (max - min) * int((*q)[index].Max-(*q)[index].Min)
		brightness := uint8(math.Min(float64(100000/area), 0xff))
		fmt.Printf("area: %v brightness: %v\n", area, brightness)
		plotRec(img, min, max, int((*q)[index].Min), int((*q)[index].Max), !ax, color.RGBA{brightness, brightness, brightness, 0xff})
		//plotLineBounded(img, min, max, int((*q)[index].min), color.RGBA{0xff, 0xff, 0xff, 0xff}, ax)
		//plotLineBounded(img, min, max, int((*q)[index].med), color.RGBA{0x00, 0x00, 0xff, 0xff}, ax)
		//plotLineBounded(img, min, max, int((*q)[index].max), color.RGBA{0xff, 0xff, 0xff, 0xff}, ax)
		//plotLineBounded(img, int((*q)[index].min), int((*q)[index].max), min, color.RGBA{0xff, 0xff, 0xff, 0xff}, !ax)
		//plotLineBounded(img, int((*q)[index].min), int((*q)[index].max), max, color.RGBA{0xff, 0xff, 0xff, 0xff}, !ax)
	}
	recPlotQuad(img, q, index*2+1, int((*q)[index].Min), int((*q)[index].med), !ax, depth-1)
	recPlotQuad(img, q, index*2+2, int((*q)[index].med), int((*q)[index].Max), !ax, depth-1)
}

func plotRec(img *image.RGBA, minA, maxA, minB, maxB int, ax axis, c color.RGBA) {
	if ax == xAxis {
		for x := minA; x <= maxA; x++ {
			for y := minB; y <= maxB; y++ {
				img.Set(x, y, c)
			}
		}
	} else {
		for x := minB; x <= maxB; x++ {
			for y := minA; y <= maxA; y++ {
				img.Set(x, y, c)
			}
		}
	}
}

type bound struct {
	Min float64
	Max float64
}

func (b bound) getSide(s side) float64 {
	if s == minSide {
		return b.Min
	}
	return b.Max
}

type Box struct {
	X        bound
	Y        bound
	pts      []Pt
	subBoxes []Box
	ax       axis
	depth    int
	Weight   float64
}

func (b Box) getBound(ax axis) bound {
	if ax == xAxis {
		return b.X
	}
	return b.Y
}

func (b *Box) setMin(ax axis, val float64) {
	if ax == xAxis {
		b.X.Min = val
	} else {
		b.Y.Min = val
	}
}

func (b *Box) setMax(ax axis, val float64) {
	if ax == xAxis {
		b.X.Max = val
	} else {
		b.Y.Max = val
	}
}

func (b *Box) String() string {
	return fmt.Sprintf("x: %f, %f y: %f, %f", b.X.Min, b.X.Max, b.Y.Min, b.Y.Max)
}

type BoxQuad struct {
	Boxes []Box
	Depth int
}

func (b *BoxQuad) Leaves() []Box {
	return b.Boxes[1<<uint(b.Depth)-1:]
}

func NewBoxQuad(depth int, pts []Pt) *BoxQuad {
	xBound := bound{ptsMin(pts, xAxis), ptsMax(pts, xAxis)}
	yBound := bound{ptsMin(pts, yAxis), ptsMax(pts, yAxis)}
	topBox := Box{xBound, yBound, pts, nil, xAxis, depth, 0}
	boxes := []Box{topBox}
	for i := 0; i < len(boxes); i++ {
		if boxes[i].depth > 0 {
			lower, upper := divideBox(boxes[i])
			boxes = append(boxes, lower)
			boxes = append(boxes, upper)
		}
	}
	q := &BoxQuad{boxes, depth}
	q.fixBounds()
	return q
}

func divideBox(b Box) (Box, Box) {
	mid := midpoint(b.pts)
	low := 0
	high := len(b.pts)
	i := partition(b.pts, low, high, b.ax)
	for i != mid {
		if i < mid {
			low = i + 1
			i = partition(b.pts, low, high, b.ax)
		} else {
			high = i
			i = partition(b.pts, low, high, b.ax)
		}
	}
	lower := Box{b.X, b.Y, b.pts[:mid], nil, !b.ax, b.depth - 1, float64(mid)}
	upper := Box{b.X, b.Y, b.pts[mid:], nil, !b.ax, b.depth - 1, float64(len(b.pts) - mid)}
	lower.setMax(b.ax, b.pts[mid].get(b.ax))
	upper.setMin(b.ax, b.pts[mid].get(b.ax))
	return lower, upper
}

func (q *BoxQuad) fixBounds() {
	q.fixBoundOnSide(xAxis, minSide)
	q.fixBoundOnSide(xAxis, maxSide)
	q.fixBoundOnSide(yAxis, minSide)
	q.fixBoundOnSide(yAxis, maxSide)
}

func (q *BoxQuad) fixBoundOnSide(ax axis, s side) {
	m := make(map[float64][]int)
	leaves := q.Leaves()
	for i, b := range leaves {
		z := b.getBound(ax).getSide(!s)
		m[z] = append(m[z], i)
	}

	for i, b := range leaves {
		z := b.getBound(ax).getSide(s)
		match := false
		indices := m[z]
		for _, j := range indices {
			c := leaves[j]
			if i != j && c.getBound(!ax).Min <= b.getBound(!ax).Max && b.getBound(!ax).Min <= c.getBound(!ax).Max {
				match = true
				break
			}
		}
		if !match {
			leaves[i] = leaves[i].fixBound(ax, s)
		}
	}
}

func (b Box) fixBound(ax axis, s side) Box {
	if ax == xAxis {
		if s == minSide {
			b.X.Min = ptsMin(b.pts, ax)
		} else {
			b.X.Max = ptsMax(b.pts, ax)
		}
	} else {
		if s == minSide {
			b.Y.Min = ptsMin(b.pts, ax)
		} else {
			b.Y.Max = ptsMax(b.pts, ax)
		}
	}
	return b
}

//func (q *BoxQuad) fixWeights() {
//	area := sumArea(q.Leaves())
//	for i, b := range q.Boxes {
//		q.Boxes[i].Weight = b.Area() / area
//	}
//}

func (b Box) Area() float64 {
	return (b.X.Max - b.X.Min) * (b.Y.Max - b.Y.Min)
}

func (b Box) Density() float64 {
	return b.Weight / b.Area()
}

func (b Box) percentBefore(ax axis, z float64) float64 {
	if ax == xAxis {
		if b.X.Max <= z {
			return 1
		}
		if z <= b.X.Min {
			return 0
		}
		return (z - b.X.Min) / (b.X.Max - b.X.Min)
	}
	if b.Y.Max <= z {
		return b.Area()
	}
	if z <= b.Y.Min {
		return 0
	}
	return (z - b.Y.Min) / (b.Y.Max - b.Y.Min)
}

func (b Box) weightBefore(ax axis, z float64) float64 {
	return b.Weight * b.percentBefore(ax, z)
}

func (q *BoxQuad) minArea() float64 {
	min := math.MaxFloat64
	for _, b := range q.Leaves() {
		if a := b.Area(); a < min {
			min = a
		}
	}
	return min
}

func (q *BoxQuad) histogram() map[int]int {
	m := make(map[int]int)
	for _, b := range q.Leaves() {
		m[len(b.pts)]++
	}
	return m
}

func (q *BoxQuad) boxHistogram() map[int]int {
	m := make(map[int]int)
	for _, b := range q.Leaves() {
		m[len(b.subBoxes)]++
	}
	return m
}

func printHistogram(h map[int]int) {
	var keys []int
	for k := range h {
		keys = append(keys, k)
	}
	sort.Ints(keys)
	fmt.Println("len(pts/boxes): #boxes")
	for _, k := range keys {
		fmt.Printf("%d: %d\n", k, h[k])
	}
}

func NewBoxQuadFromBoxes(depth int, b []Box) *BoxQuad {
	xBound := bound{boxMin(b, xAxis), boxMax(b, xAxis)}
	yBound := bound{boxMin(b, yAxis), boxMax(b, yAxis)}
	topBox := Box{xBound, yBound, nil, b, xAxis, depth, sumWeights(b)}
	boxes := []Box{topBox}
	for i := 0; i < len(boxes); i++ {
		if boxes[i].depth > 0 {
			lower, upper := divideBoxOfBoxes(boxes[i])
			boxes = append(boxes, lower)
			boxes = append(boxes, upper)
		}
	}
	q := &BoxQuad{boxes, depth}
	return q
}

func boxMin(boxes []Box, ax axis) float64 {
	if len(boxes) == 0 {
		return math.MaxFloat64
	}
	min := boxes[0].getBound(ax).Min
	for _, b := range boxes {
		if z := b.getBound(ax).Min; z < min {
			min = z
		}
	}
	return min
}

func boxMax(boxes []Box, ax axis) float64 {
	if len(boxes) == 0 {
		return -math.MaxFloat64
	}
	max := boxes[0].getBound(ax).Max
	for _, b := range boxes {
		if z := b.getBound(ax).Max; z > max {
			max = z
		}
	}
	return max
}

func sumArea(boxes []Box) float64 {
	var w float64
	for _, b := range boxes {
		w += b.Area()
	}
	return w
}

func sumWeights(boxes []Box) float64 {
	var w float64
	for _, b := range boxes {
		w += b.Weight
	}
	return w
}

type boxDivision struct {
	boxes   []Box
	indices []int
}

func newBoxDivision(boxes []Box, n, slice int) boxDivision {
	indices := []int{0}
	for i := 0; i < n; i++ {
		if i < slice {
			indices = append(indices, 0)
		} else {
			indices = append(indices, len(boxes))
		}
	}
	return boxDivision{boxes, indices}
}

func (d boxDivision) lengths() []int {
	var lengths []int
	for i := 0; i < len(d.indices)-1; i++ {
		lengths = append(lengths, d.indices[i+1]-d.indices[i])
	}
	return lengths
}

func (d boxDivision) weights() []float64 {
	var weights []float64
	for i := 0; i < len(d.indices)-1; i++ {
		weights = append(weights, sumWeights(d.slice(i)))
	}
	return weights
}

func (d boxDivision) String() string {
	out := "Division\n"
	for i := 0; i < len(d.indices)-1; i++ {
		zeds := []float64{}
		for _, b := range d.slice(i) {
			zeds = append(zeds, b.getBound(xAxis).Min)
		}
		out += fmt.Sprintf("%d: %v\n", i, zeds)
	}
	return out
}

func (d boxDivision) slice(i int) []Box {
	return d.boxes[d.indices[i]:d.indices[i+1]]
}

func (d boxDivision) first(slice int) Box {
	return d.boxes[d.indices[slice]]
}

func (d boxDivision) setFirst(slice int, b Box) {
	d.boxes[d.indices[slice]] = b
}

func (d boxDivision) last(slice int) Box {
	return d.boxes[d.indices[slice+1]-1]
}

func (d boxDivision) setLast(slice int, b Box) {
	d.boxes[d.indices[slice+1]-1] = b
}

func (d boxDivision) length(i int) int {
	return d.indices[i+1] - d.indices[i]
}

type divisionIterator struct {
	d            boxDivision
	slice, index int
}

func (d boxDivision) newIterator(slice int) *divisionIterator {
	return &divisionIterator{d, slice, d.indices[slice]}
}

func (iter *divisionIterator) next() (b Box, ok bool) {
	if iter.index == iter.d.indices[iter.slice+1] {
		return b, false
	}
	b = iter.d.boxes[iter.index]
	iter.index++
	return b, true
}

func (iter *divisionIterator) move(slice int) {
	d := iter.d
	index := iter.index - 1
	b := d.boxes[index]
	if slice < iter.slice {
		d.boxes[index] = d.first(iter.slice)
		for s := iter.slice; s > slice+1; s-- {
			d.setFirst(s, d.first(s-1))
			d.indices[s]++
		}
		d.setFirst(slice+1, b)
		d.indices[slice+1]++
	} else if slice > iter.slice {
		d.boxes[index] = d.last(iter.slice)
		for s := iter.slice; s < slice-1; s++ {
			d.setLast(s, d.last(s+1))
			d.indices[s]--
		}
		d.setLast(slice-1, b)
		d.indices[slice]--
		iter.index--
	}
}

const (
	LOW     = 0
	UNKNOWN = 1
	OVERLAP = 2
	HIGH    = 3
)

func divideBoxOfBoxes(b Box) (Box, Box) {
	midWeight := b.Weight / 2
	fmt.Printf("depth: %d, boxes: %d midweight: %f dims: %s\n", b.depth, len(b.subBoxes), midWeight, b.String())
	div := newBoxDivision(b.subBoxes, 4, UNKNOWN)
	fmt.Printf("Division: %v, weights: %v\n", div.lengths(), div.weights())
	for div.length(UNKNOWN) != 0 {
		zMin := div.first(UNKNOWN).getBound(b.ax).Min
		zMax := div.first(UNKNOWN).getBound(b.ax).Max
		// This could be made more efficient by precomputing the weights.
		wMin, wMax := weightBefore(b.subBoxes, b.ax, zMin, zMax)
		fmt.Printf("zMin, zMax = %f, %f  wMin, wMax = %f, %f", zMin, zMax, wMin, wMax)
		if wMax < midWeight {
			// Everything with max <= zMax goes in low
			fmt.Printf(" sending to low.\n")
			intervalPartitionLow(div, b.ax, zMax)
		} else if wMin > midWeight {
			// Everything with min >= zMin goes in high
			fmt.Printf(" sending to high.\n")
			intervalPartitionHigh(div, b.ax, zMin)
		} else {
			// Everything with min <= zMin and max >= zMax goes in overlap
			// Everything with maxZ <= zMin goes in low
			// Everything with minZ <= zMax goes in high
			fmt.Printf(" splitting.\n")
			intervalPartitionMid(div, b.ax, zMin, zMax)
		}
		fmt.Printf("Division: %v, weights: %v\n", div.lengths(), div.weights())
		//printRegions(div, b.ax)
	}

	// TODO: Partition the overlap. Currently they are in lower.
	// Idea!: ORDER LOW, UNKNOWN, OVERLAP, HIGH. Then when splitting the overlap, leave
	// the lower half in place and append upper half to the high.
	zSplit := findSplit(div, b.ax, midWeight)
	fmt.Printf("Final Z Split: %f.\n", zSplit)
	splitAndAppend(div, zSplit)
	lowerBoxes := div.boxes[:div.indices[HIGH]]
	upperBoxes := div.boxes[div.indices[HIGH]:]
	lower := Box{b.X, b.Y, nil, lowerBoxes, !b.ax, b.depth - 1, sumWeights(lowerBoxes)}
	upper := Box{b.X, b.Y, nil, upperBoxes, !b.ax, b.depth - 1, sumWeights(upperBoxes)}
	lower.setMax(b.ax, math.Max(boxMax(lowerBoxes, b.ax), lower.getBound(b.ax).Min))
	upper.setMin(b.ax, math.Min(boxMin(upperBoxes, b.ax), upper.getBound(b.ax).Max))
	return lower, upper
}

func weightBefore(boxes []Box, ax axis, zMin, zMax float64) (float64, float64) {
	var wMin, wMax float64
	for _, b := range boxes {
		wMin += b.weightBefore(ax, zMin)
		wMax += b.weightBefore(ax, zMax)
	}
	return wMin, wMax
}

// Everything with max <= zMax goes in low
func intervalPartitionLow(d boxDivision, ax axis, zMax float64) {
	iter := d.newIterator(UNKNOWN)
	b, ok := iter.next()
	for ok {
		if b.getBound(ax).Max <= zMax {
			iter.move(LOW)
		}
		b, ok = iter.next()
	}
}

// Everything with min >= zMin goes in high
func intervalPartitionHigh(d boxDivision, ax axis, zMin float64) {
	iter := d.newIterator(UNKNOWN)
	b, ok := iter.next()
	for ok {
		if b.getBound(ax).Min >= zMin {
			iter.move(HIGH)
		}
		b, ok = iter.next()
	}
}

// Everything with min <= zMin and max >= zMax goes in overlap
// Everything with max <= zMin goes in low
// Everything with min <= zMax goes in high
func intervalPartitionMid(d boxDivision, ax axis, zMin, zMax float64) {
	iter := d.newIterator(UNKNOWN)
	b, ok := iter.next()
	for ok {
		if b.getBound(ax).Max <= zMin {
			iter.move(LOW)
		} else if b.getBound(ax).Min >= zMax {
			iter.move(HIGH)
		} else if b.getBound(ax).Min <= zMin && b.getBound(ax).Max >= zMax {
			iter.move(OVERLAP)
		}
		b, ok = iter.next()
	}
}

func findZSplit(d boxDivision, ax axis, targetWeight float64) float64 {
	lowerWeight := sumWeights(d.split(LOW))
	// Define the weight contribution from a box with weight w at z as
	// w(z) = (w/(zMax -zMin))( z - zMin)
	// Thus w(zMax) = w and w(zMin) = 0.
	// Solving for z gives
	// z = (target weight - lowerWeight + sum w zMin /(zMax - zMin) )
	//     / (sum w/(zMax - zMin)
	numer := targetWeight - lowerWeight
	var denom float64
	for _, b := range d.split(OVERLAP) {
		numer += b.Weight * b.Get(ax).Min / (b.Get(ax).Max - b.Get(ax).Min)
		denom += b.Weight / (b.Get(ax).Max - b.Get(ax).Min)
	}
	return numer / denom
}

func splitAndAppend(d boxDivision, ax axis, zSplit float64) []Box {
	boxes = make([]Box, len(d.boxes)+len(d.slice(OVERLAP)))
	for i, b := range d.boxes {
		boxes[i] = b
	}
	for i := 0; 
}

func printBoxes(prefix string, boxes []Box) {
	for i, b := range boxes {
		fmt.Printf("%s\t%d. %s\n", prefix, i, b.String())
	}
}

type basicBox struct {
	x bound
	y bound
}

func distinctBoxes(boxes []Box) int {
	set := make(map[basicBox]bool)
	for _, b := range boxes {
		bb := basicBox{b.X, b.Y}
		set[bb] = true
	}
	return len(set)
}

func printRegions(div boxDivision, ax axis) {
	fmt.Printf("zAxis = %s, distinct boxes %d\n", axisName[ax], distinctBoxes(div.boxes))
	fmt.Printf("overlapMinZ = %f\n", boxMin(div.slice(OVERLAP), ax))
	fmt.Printf("lowerMaxZ = %f\n", boxMax(div.slice(LOW), ax))
	fmt.Printf("upperMinZ = %f\n", boxMin(div.slice(HIGH), ax))
	fmt.Printf("overlapMaxZ = %f\n", boxMax(div.slice(OVERLAP), ax))
	printBoxes("lower", div.slice(LOW))
	printBoxes("upper", div.slice(HIGH))
	printBoxes("overlap", div.slice(OVERLAP))
	fmt.Println()
}
