package main

import "github.com/ahhoefel/cdf/scene"
import "github.com/ahhoefel/cdf"

const (
	width     = 800
	height    = 800
	depth     = 5
	numPoints = 1000
)

func main() {
	ptsA := cdf.RandomPointsNorm(numPoints, 60, 300)
	ptsB := cdf.RandomPointsNorm(numPoints, 60, 400)
	qA := cdf.NewBoxQuad(depth, ptsA)
	qB := cdf.NewBoxQuad(depth, ptsB)
	s := scene.New()
	s.AddPoints(ptsA...)
	s.AddPoints(ptsB...)
	s.AddBoxes(append(qA.Leaves(), qB.Leaves()...)...)
	s.Paint("image.png", width, height)
}
