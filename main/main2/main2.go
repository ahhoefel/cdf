package main

import "github.com/ahhoefel/cdf/scene"
import "github.com/ahhoefel/cdf"

const (
	width     = 800
	height    = 800
	depth     = 3
	numPoints = 300
)

func main() {
	ptsA := cdf.RandomPointsNorm(numPoints, 60, 300)
	ptsB := cdf.RandomPointsNorm(numPoints, 60, 400)
	q := cdf.NewBoxQuad(depth, append(ptsA, ptsB...))
	s := scene.New()
	s.AddPoints(ptsA...)
	s.AddPoints(ptsB...)
	s.AddBoxes(q.Leaves()...)
	s.Paint("image.png", width, height)
}
