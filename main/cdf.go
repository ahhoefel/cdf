package main

import (
	"bufio"
	"fmt"
	"image"
	"image/color"
	"image/png"
	"log"
	"math"
	"math/rand"
	"os"
	"path"
	"sort"
	"time"
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

const (
	debugPlot = true
)

func main() {
	rand.Seed(time.Now().Unix())
	main8()
}

func main1() {
	pts := randomPointsCorr(numPoints, 200)
	xSketch := sketchPts(10, pts, xAxis)
	ySketch := sketchPts(10, pts, yAxis)
	fmt.Println("X Sketch")
	fmt.Println(xSketch.String())
	fmt.Println("Y Sketch")
	fmt.Println(ySketch.String())
	img := makePlot(800, 800, pts, xSketch, ySketch)
	writeImage("image.png", img)
}

func makePlot(w, h int, pts []Pt, xSketch *sketch, ySketch *sketch) *image.RGBA {
	img := image.NewRGBA(image.Rect(0, 0, w, h))
	setBackground(img)
	plotSketch(img, xSketch, xAxis)
	plotSketch(img, ySketch, yAxis)
	plotPoints(img, pts)
	return img
}

func main2() {
	pts := randomPointsCorr(numPoints, 200)
	sl := newSlider(10, pts, xAxis)
	img := makePlot2(800, 800, pts, sl)
	writeImage("image.png", img)
}

func makePlot2(w, h int, pts []Pt, s *slider) *image.RGBA {
	img := image.NewRGBA(image.Rect(0, 0, w, h))
	setBackground(img)
	plotSlider(img, s)
	plotPoints(img, pts)
	return img
}

func main3() {
	pts := randomPointsCorr(numPoints, 200)
	q := newQuad(depth, pts, xAxis)
	img := makePlot3(800, 800, pts, q, depth)
	writeImage("image.png", img)
}

func makePlot3(w, h int, pts []Pt, q *quad, depth int) *image.RGBA {
	img := image.NewRGBA(image.Rect(0, 0, w, h))
	setBackground(img)
	plotQuad(img, q, depth)
	plotPoints(img, pts)
	return img
}

func main4() {
	pts := randomPointsTwoNorm(numPoints, 60, 400, 200)
	q := newBoxQuad(depth, pts)
	img := makePlot4(800, 800, pts, q)
	writeImage("image.png", img)
	printHistogram(q.histogram())
}

func makePlot4(w, h int, pts []Pt, q *boxQuad) *image.RGBA {
	img := image.NewRGBA(image.Rect(0, 0, w, h))
	setBackground(img)
	plotBoxQuad(img, q)
	if showPoints {
		plotPoints(img, pts)
	}
	return img
}

func main5() {
	ptsA := randomPointsNorm(numPoints, 60, 300)
	ptsB := randomPointsNorm(numPoints, 60, 400)
	qA := newBoxQuad(depth, ptsA)
	qB := newBoxQuad(depth, ptsB)
	q := newBoxQuadFromBoxes(depth, append(qA.leaves(), qB.leaves()...))
	img := makePlot5(800, 800, ptsA, ptsB, q, qA, qB)
	writeImage("image.png", img)
	printHistogram(q.boxHistogram())
}

func makePlot5(w, h int, ptsA, ptsB []Pt, q, qA, qB *boxQuad) *image.RGBA {
	img := image.NewRGBA(image.Rect(0, 0, w, h))
	setBackground(img)
	//plotBoxQuadByWeight(img, qA)
	//plotBoxQuadByWeight(img, qB)
	plotBoxQuadByWeight(img, q)
	if showPoints {
		plotPoints(img, ptsA)
		plotPoints(img, ptsB)
	}
	return img
}

func testBox(z float64) box {
	return box{bound{z, z}, bound{z, z}, nil, nil, xAxis, 0, 0}
}

func main6() {
	boxes := []box{
		testBox(1),
		testBox(2),
		testBox(3),
		testBox(4),
		testBox(5)}
	div := newBoxDivision(boxes, 4, UNKNOWN)
	iter := div.newIterator(UNKNOWN)
	_, _ = iter.next()
	iter.move(HIGH)
	fmt.Printf("Division: %v\n", div.String())
	_, _ = iter.next()
	iter.move(LOW)
	fmt.Printf("Division: %v\n", div.String())
	_, _ = iter.next()
	iter.move(LOW)
	fmt.Printf("Division: %v\n", div.String())
	_, _ = iter.next()
	iter.move(OVERLAP)
	fmt.Printf("Division: %v\n", div.String())
	_, _ = iter.next()
	iter.move(HIGH)
	fmt.Printf("Division: %v\n", div.String())
}

func main7() {
	makeBox := func(x, y bound) box {
		return box{x, y, nil, nil, xAxis, 0, 1}
	}
	boxes := []box{
		makeBox(bound{0, 10}, bound{0, 20}),
		makeBox(bound{10, 20}, bound{0, 20}),
		makeBox(bound{20, 30}, bound{0, 20}),
	}
	q := newBoxQuadFromBoxes(1, boxes)
	for _, b := range q.leaves() {
		fmt.Printf("x: %f, %f y: %f, %f\n", b.x.min, b.x.max, b.y.min, b.y.max)
	}
}

func main8() {
	fmt.Println("Main 8")
	b := box{bound{50, 100}, bound{50, 300}, nil, nil, xAxis, 0, 0.5}
	img := image.NewRGBA(image.Rect(0, 0, xPixels, yPixels))

	setBackground(img)
	plotBox(img, b, color.RGBA{0xff, 0xff, 0xff, 0x00})
	writeImage("image.png", img)
}

type Pt struct {
	x, y float64
}

func (p Pt) get(ax axis) float64 {
	if ax == xAxis {
		return p.x
	}
	return p.y
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

func randomPointsNorm(n int, scale float64, shift float64) []Pt {
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

func writeImage(fileName string, img image.Image) {
	f, err := os.Create(path.Join(basePath, fileName))
	check(err)
	w := bufio.NewWriter(f)
	check(png.Encode(w, img))
	check(w.Flush())
	check(f.Close())
}

func setBackground(img *image.RGBA) {
	for x := img.Bounds().Min.X; x < img.Bounds().Max.X; x++ {
		for y := img.Bounds().Min.Y; y < img.Bounds().Max.Y; y++ {
			img.Set(x, y, color.RGBA{0x00, 0x00, 0x00, 0xff})
		}
	}
}

func plotPoints(img *image.RGBA, pts []Pt) {
	for _, p := range pts {
		img.Set(int(p.x), int(p.y), color.RGBA{0xff, 0xff, 0xff, 0xff})
	}
}

func plotSketch(img *image.RGBA, s *sketch, a axis) {
	plotLine(img, int(s.min), color.RGBA{0xff, 0x00, 0x00, 0xff}, a)
	plotLine(img, int(s.max), color.RGBA{0x00, 0x00, 0xff, 0xff}, a)
	for _, x := range s.quantiles {
		plotLine(img, int(x), color.RGBA{0xff, 0xff, 0xff, 0xff}, a)
	}
}

func plotLine(img *image.RGBA, z int, c color.RGBA, a axis) {
	if a == xAxis {
		for y := img.Bounds().Min.Y; y <= img.Bounds().Max.Y; y++ {
			img.Set(z, y, c)
		}
	} else {
		for x := img.Bounds().Min.X; x <= img.Bounds().Max.X; x++ {
			img.Set(x, z, c)
		}
	}
}

func check(err error) {
	if err != nil {
		log.Fatal(err)
	}
}

type sketch struct {
	q         int
	n         int
	min       float64
	max       float64
	quantiles []float64
}

func newSketch(q int, a []float64) *sketch {
	sort.Float64s(a)
	var quantiles []float64
	for i := 1; i < q; i++ {
		quantiles = append(quantiles, a[i*len(a)/q])
	}
	return &sketch{q, len(a), a[0], a[len(a)-1], quantiles}
}

func sketchPts(q int, pts []Pt, ax axis) *sketch {
	var a []float64
	if ax == xAxis {
		for _, p := range pts {
			a = append(a, p.x)
		}
	} else {
		for _, p := range pts {
			a = append(a, p.y)
		}
	}
	return newSketch(q, a)
}

func (s sketch) String() string {
	return fmt.Sprintf("Sketch q:%d n:%d min:%f max:%f quantiles:%v", s.q, s.n, s.min, s.max, s.quantiles)
}

type sortablePts struct {
	pts []Pt
	ax  axis
}

func (a sortablePts) Len() int      { return len(a.pts) }
func (a sortablePts) Swap(i, j int) { a.pts[i], a.pts[j] = a.pts[j], a.pts[i] }
func (a sortablePts) Less(i, j int) bool {
	if a.ax == xAxis {
		return a.pts[i].x < a.pts[j].x
	}
	return a.pts[i].y < a.pts[j].y
}

type slider struct {
	ax          axis
	s           *sketch
	subSketches []*sketch
}

func newSlider(q int, pts []Pt, ax axis) *slider {
	var subSketches []*sketch
	sPts := sortablePts{pts, ax}
	sort.Sort(sPts)
	pts = sPts.pts
	s := sketchPts(q, pts, ax)
	var j int
	for _, z := range s.quantiles {
		if ax == xAxis {
			for pts[j].x < z {
				j++
			}
		} else {
			for pts[j].y < z {
				j++
			}
		}
		subSketches = append(subSketches, sketchPts(q, pts[:j], !ax))
		pts = pts[j:]
		j = 0
	}
	subSketches = append(subSketches, sketchPts(q, pts, !ax))
	return &slider{ax, s, subSketches}
}

func plotSlider(img *image.RGBA, sl *slider) {
	plotSketch(img, sl.s, sl.ax)
	prevZ := sl.s.min
	for i, z := range sl.s.quantiles {
		plotSketchBounded(img, int(prevZ), int(z), sl.subSketches[i], !sl.ax)
		prevZ = z
	}
	plotSketchBounded(img, int(prevZ), int(sl.s.max), sl.subSketches[len(sl.subSketches)-1], !sl.ax)
}

func plotSketchBounded(img *image.RGBA, minZ, maxZ int, s *sketch, a axis) {
	plotLineBounded(img, minZ, maxZ, int(s.min), color.RGBA{0xff, 0x00, 0x00, 0xff}, a)
	plotLineBounded(img, minZ, maxZ, int(s.max), color.RGBA{0x00, 0x00, 0xff, 0xff}, a)
	for _, x := range s.quantiles {
		plotLineBounded(img, minZ, maxZ, int(x), color.RGBA{0xff, 0xff, 0xff, 0xff}, a)
	}
}

func plotLineBounded(img *image.RGBA, minA, maxA, z int, c color.RGBA, ax axis) {
	if ax == xAxis {
		for y := minA; y <= maxA; y++ {
			img.Set(z, y, c)
		}
	} else {
		for x := minA; x <= maxA; x++ {
			img.Set(x, z, c)
		}
	}
}

type pivot struct {
	min float64
	max float64
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
	d.p.min = ptsMin(pts, ax)
	d.p.max = ptsMax(pts, ax)
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
		area := (max - min) * int((*q)[index].max-(*q)[index].min)
		brightness := uint8(math.Min(float64(100000/area), 0xff))
		fmt.Printf("area: %v brightness: %v\n", area, brightness)
		plotRec(img, min, max, int((*q)[index].min), int((*q)[index].max), !ax, color.RGBA{brightness, brightness, brightness, 0xff})
		//plotLineBounded(img, min, max, int((*q)[index].min), color.RGBA{0xff, 0xff, 0xff, 0xff}, ax)
		//plotLineBounded(img, min, max, int((*q)[index].med), color.RGBA{0x00, 0x00, 0xff, 0xff}, ax)
		//plotLineBounded(img, min, max, int((*q)[index].max), color.RGBA{0xff, 0xff, 0xff, 0xff}, ax)
		//plotLineBounded(img, int((*q)[index].min), int((*q)[index].max), min, color.RGBA{0xff, 0xff, 0xff, 0xff}, !ax)
		//plotLineBounded(img, int((*q)[index].min), int((*q)[index].max), max, color.RGBA{0xff, 0xff, 0xff, 0xff}, !ax)
	}
	recPlotQuad(img, q, index*2+1, int((*q)[index].min), int((*q)[index].med), !ax, depth-1)
	recPlotQuad(img, q, index*2+2, int((*q)[index].med), int((*q)[index].max), !ax, depth-1)
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
	min float64
	max float64
}

func (b bound) getSide(s side) float64 {
	if s == minSide {
		return b.min
	}
	return b.max
}

type box struct {
	x        bound
	y        bound
	pts      []Pt
	subBoxes []box
	ax       axis
	depth    int
	weight   float64
}

func (b box) getBound(ax axis) bound {
	if ax == xAxis {
		return b.x
	}
	return b.y
}

func (b *box) setMin(ax axis, val float64) {
	if ax == xAxis {
		b.x.min = val
	} else {
		b.y.min = val
	}
}

func (b *box) setMax(ax axis, val float64) {
	if ax == xAxis {
		b.x.max = val
	} else {
		b.y.max = val
	}
}

func (b *box) String() string {
	return fmt.Sprintf("x: %f, %f y: %f, %f", b.x.min, b.x.max, b.y.min, b.y.max)
}

type boxQuad struct {
	boxes []box
	depth int
}

func (b *boxQuad) leaves() []box {
	return b.boxes[1<<uint(b.depth)-1:]
}

func newBoxQuad(depth int, pts []Pt) *boxQuad {
	xBound := bound{ptsMin(pts, xAxis), ptsMax(pts, xAxis)}
	yBound := bound{ptsMin(pts, yAxis), ptsMax(pts, yAxis)}
	topBox := box{xBound, yBound, pts, nil, xAxis, depth, 0}
	boxes := []box{topBox}
	for i := 0; i < len(boxes); i++ {
		if boxes[i].depth > 0 {
			lower, upper := divideBox(boxes[i])
			boxes = append(boxes, lower)
			boxes = append(boxes, upper)
		}
	}
	q := &boxQuad{boxes, depth}
	q.fixBounds()
	q.fixWeights()
	return q
}

func divideBox(b box) (box, box) {
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
	lower := box{b.x, b.y, b.pts[:mid], nil, !b.ax, b.depth - 1, 0}
	upper := box{b.x, b.y, b.pts[mid:], nil, !b.ax, b.depth - 1, 0}
	lower.setMax(b.ax, b.pts[mid].get(b.ax))
	upper.setMin(b.ax, b.pts[mid].get(b.ax))
	return lower, upper
}

func (q *boxQuad) fixBounds() {
	q.fixBoundOnSide(xAxis, minSide)
	q.fixBoundOnSide(xAxis, maxSide)
	q.fixBoundOnSide(yAxis, minSide)
	q.fixBoundOnSide(yAxis, maxSide)
}

func (q *boxQuad) fixBoundOnSide(ax axis, s side) {
	m := make(map[float64][]int)
	leaves := q.leaves()
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
			if i != j && c.getBound(!ax).min <= b.getBound(!ax).max && b.getBound(!ax).min <= c.getBound(!ax).max {
				match = true
				break
			}
		}
		if !match {
			leaves[i] = leaves[i].fixBound(ax, s)
		}
	}
}

func (b box) fixBound(ax axis, s side) box {
	if ax == xAxis {
		if s == minSide {
			b.x.min = ptsMin(b.pts, ax)
		} else {
			b.x.max = ptsMax(b.pts, ax)
		}
	} else {
		if s == minSide {
			b.y.min = ptsMin(b.pts, ax)
		} else {
			b.y.max = ptsMax(b.pts, ax)
		}
	}
	return b
}

func (q *boxQuad) fixWeights() {
	area := sumArea(q.leaves())
	for i, b := range q.boxes {
		q.boxes[i].weight = b.area() / area
	}
}

func plotBoxQuad(img *image.RGBA, q *boxQuad) {
	minArea := q.minArea()
	for _, b := range q.leaves() {
		brightness := uint8(math.Min(float64(2*0xff*minArea/b.area()), 0xff))
		c := color.RGBA{brightness, brightness, brightness, 0xff}
		plotBox(img, b, c)
	}
}

func plotBoxQuadByWeight(img *image.RGBA, q *boxQuad) {
	minArea := q.minArea()
	var brightnessHistogram [8]int
	for _, b := range q.leaves() {
		brightness := uint8(math.Min(float64(0xff*b.weight*minArea/b.area()), 0xff))
		brightnessHistogram[int(brightness>>5)]++
		c := color.RGBA{0xff, 0xff, 0xff, brightness}
		plotBox(img, b, c)
	}
	if debugPlot {
		fmt.Print("Brightness histogram (0->8): ")
		for _, v := range brightnessHistogram {
			fmt.Printf("%d ", v)
		}
		fmt.Print("\n")
	}
}

func scaleByAlpha(c, a uint8) uint8 {
	return uint8(c * a / 0xff)
}

func MinInt(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func uint8Cap(a, b uint8) uint8 {
	return uint8(MinInt(int(a)+int(b), 0xff))
}

func mergeColors(a, b color.RGBA) color.RGBA {
	return color.RGBA{
		uint8Cap(scaleByAlpha(a.R, a.A), scaleByAlpha(b.R, b.A)),
		uint8Cap(scaleByAlpha(a.G, a.A), scaleByAlpha(b.G, b.A)),
		uint8Cap(scaleByAlpha(a.B, a.A), scaleByAlpha(b.B, b.A)),
		uint8Cap(a.A/2, b.A/2)}
}

func plotBox(img *image.RGBA, b box, c color.RGBA) {
	for x := int(b.x.min); x <= int(b.x.max); x++ {
		for y := int(b.y.min); y <= int(b.y.max); y++ {
			r, g, b, a := img.At(x, y).RGBA()
			d := color.RGBA{uint8(r >> 8), uint8(g >> 8), uint8(b >> 8), uint8(a >> 8)}
			//fmt.Printf("d: %v, c: %v merge: %v", d, c, mergeColors(c, d))
			img.Set(x, y, mergeColors(c, d))
			//img.Set(x, y, color.RGBA{0xff, 0xff, 0xff, 0x0})
			//fmt.Printf(" after: %v\n", img.At(x, y))
		}
	}

	if showBorder {
		outlineC := color.RGBA{0xff, 0xff, 0x00, 0x90}
		for y := int(b.y.min); y <= int(b.y.max); y++ {
			img.Set(int(b.x.min), y, outlineC)
			img.Set(int(b.x.max), y, outlineC)
		}
		for x := int(b.x.min); x <= int(b.x.max); x++ {
			img.Set(x, int(b.y.min), outlineC)
			img.Set(x, int(b.y.max), outlineC)
		}
	}
}

func (b box) area() float64 {
	return (b.x.max - b.x.min) * (b.y.max - b.y.min)
}

func (b box) areaBefore(ax axis, z float64) float64 {
	if ax == xAxis {
		if b.x.max <= z {
			return b.area()
		}
		if z <= b.x.min {
			return 0
		}
		return (b.x.max - z) * (b.y.max - b.y.min)
	}
	if b.y.max <= z {
		return b.area()
	}
	if z <= b.y.min {
		return 0
	}
	return (b.x.max - b.x.min) * (b.y.max - z)
}

func (b box) weightBefore(ax axis, z float64) float64 {
	return b.weight * b.areaBefore(ax, z) / b.area()
}

func (q *boxQuad) minArea() float64 {
	min := math.MaxFloat64
	for _, b := range q.leaves() {
		if a := b.area(); a < min {
			min = a
		}
	}
	return min
}

func (q *boxQuad) histogram() map[int]int {
	m := make(map[int]int)
	for _, b := range q.leaves() {
		m[len(b.pts)]++
	}
	return m
}

func (q *boxQuad) boxHistogram() map[int]int {
	m := make(map[int]int)
	for _, b := range q.leaves() {
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

func newBoxQuadFromBoxes(depth int, b []box) *boxQuad {
	xBound := bound{boxMin(b, xAxis), boxMax(b, xAxis)}
	yBound := bound{boxMin(b, yAxis), boxMax(b, yAxis)}
	topBox := box{xBound, yBound, nil, b, xAxis, depth, sumWeights(b)}
	boxes := []box{topBox}
	for i := 0; i < len(boxes); i++ {
		if boxes[i].depth > 0 {
			lower, upper := divideBoxOfBoxes(boxes[i])
			boxes = append(boxes, lower)
			boxes = append(boxes, upper)
		}
	}
	q := &boxQuad{boxes, depth}
	return q
}

func boxMin(boxes []box, ax axis) float64 {
	if len(boxes) == 0 {
		return math.MaxFloat64
	}
	min := boxes[0].getBound(ax).min
	for _, b := range boxes {
		if z := b.getBound(ax).min; z < min {
			min = z
		}
	}
	return min
}

func boxMax(boxes []box, ax axis) float64 {
	if len(boxes) == 0 {
		return -math.MaxFloat64
	}
	max := boxes[0].getBound(ax).max
	for _, b := range boxes {
		if z := b.getBound(ax).max; z > max {
			max = z
		}
	}
	return max
}

func sumArea(boxes []box) float64 {
	var w float64
	for _, b := range boxes {
		w += b.area()
	}
	return w
}

func sumWeights(boxes []box) float64 {
	var w float64
	for _, b := range boxes {
		w += b.weight
	}
	return w
}

type boxDivision struct {
	boxes   []box
	indices []int
}

func newBoxDivision(boxes []box, n, slice int) boxDivision {
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
			zeds = append(zeds, b.getBound(xAxis).min)
		}
		out += fmt.Sprintf("%d: %v\n", i, zeds)
	}
	return out
}

func (d boxDivision) slice(i int) []box {
	return d.boxes[d.indices[i]:d.indices[i+1]]
}

func (d boxDivision) first(slice int) box {
	return d.boxes[d.indices[slice]]
}

func (d boxDivision) setFirst(slice int, b box) {
	d.boxes[d.indices[slice]] = b
}

func (d boxDivision) last(slice int) box {
	return d.boxes[d.indices[slice+1]-1]
}

func (d boxDivision) setLast(slice int, b box) {
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

func (iter *divisionIterator) next() (b box, ok bool) {
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
	OVERLAP = 0
	LOW     = 1
	UNKNOWN = 2
	HIGH    = 3
)

func divideBoxOfBoxes(b box) (box, box) {
	midWeight := b.weight / 2
	fmt.Printf("depth: %d, boxes: %d midweight: %f dims: %s\n", b.depth, len(b.subBoxes), midWeight, b.String())
	div := newBoxDivision(b.subBoxes, 4, UNKNOWN)
	fmt.Printf("Division: %v, weights: %v\n", div.lengths(), div.weights())
	for div.length(UNKNOWN) != 0 {
		zMin := div.first(UNKNOWN).getBound(b.ax).min
		zMax := div.first(UNKNOWN).getBound(b.ax).max
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
	}
	// TODO: Partition the overlap. Currently they are in lower.
	// Idea!: ORDER LOW, UNKNOWN, OVERLAP, HIGH. Then when splitting the overlap, leave
	// the lower half in place and append upper half to the high.
	lowerBoxes := div.boxes[:div.indices[HIGH]]
	upperBoxes := div.boxes[div.indices[HIGH]:]
	lower := box{b.x, b.y, nil, lowerBoxes, !b.ax, b.depth - 1, sumWeights(lowerBoxes)}
	upper := box{b.x, b.y, nil, upperBoxes, !b.ax, b.depth - 1, sumWeights(upperBoxes)}
	fmt.Printf("lower depth: %d, boxes: %d zMin, zMax: %f, %f\n", lower.depth, len(lower.subBoxes), lower.getBound(b.ax).min, lower.getBound(b.ax).max)
	fmt.Printf("upper depth: %d, boxes: %d zMin, zMax: %f, %f\n", upper.depth, len(upper.subBoxes), upper.getBound(b.ax).min, upper.getBound(b.ax).max)
	lower.setMax(b.ax, math.Max(boxMax(lowerBoxes, b.ax), lower.getBound(b.ax).min))
	upper.setMin(b.ax, math.Min(boxMin(upperBoxes, b.ax), upper.getBound(b.ax).max))
	fmt.Printf("lower depth: %d, boxes: %d zMin, zMax: %f, %f\n", lower.depth, len(lower.subBoxes), lower.getBound(b.ax).min, lower.getBound(b.ax).max)
	printBoxes(lowerBoxes)
	fmt.Printf("upper depth: %d, boxes: %d zMin, zMax: %f, %f\n", upper.depth, len(upper.subBoxes), upper.getBound(b.ax).min, upper.getBound(b.ax).max)
	printBoxes(upperBoxes)
	return lower, upper
}

func weightBefore(boxes []box, ax axis, zMin, zMax float64) (float64, float64) {
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
		if b.getBound(ax).max <= zMax {
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
		if b.getBound(ax).min >= zMin {
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
		if b.getBound(ax).max <= zMin {
			iter.move(LOW)
		} else if b.getBound(ax).min >= zMax {
			iter.move(HIGH)
		} else if b.getBound(ax).min <= zMin && b.getBound(ax).max >= zMax {
			iter.move(OVERLAP)
		}
		b, ok = iter.next()
	}
}

func printBoxes(boxes []box) {
	for i, b := range boxes {
		fmt.Printf("%d. %s\n", i, b.String())
	}
}
