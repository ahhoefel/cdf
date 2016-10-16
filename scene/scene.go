package scene

import (
	"bufio"
	"fmt"
	"image"
	"image/color"
	"image/png"
	"math"
	"os"
	"path"

	"github.com/ahhoefel/cdf"
)

const (
	basePath = "/Users/hoefel/Development/go/src/github.com/ahhoefel/cdf"
)

type Scene struct {
	boxes []cdf.Box
	pts   []cdf.Pt
}

func New() *Scene {
	return &Scene{}
}

func (s *Scene) AddBoxes(b ...cdf.Box) {
	s.boxes = append(s.boxes, b...)
}

func (s *Scene) AddPoints(p ...cdf.Pt) {
	s.pts = append(s.pts, p...)
}

func (s *Scene) Paint(filename string, w, h int) error {
	img := image.NewRGBA(image.Rect(0, 0, w, h))
	paintBackground(img)
	paintBoxes(img, s.boxes)
	paintPoints(img, s.pts)
	printDensities(s.boxes)
	err := writeImage(filename, img)
	return err
}

func writeImage(fileName string, img image.Image) error {
	f, err := os.Create(path.Join(basePath, fileName))
	if err != nil {
		return err
	}
	w := bufio.NewWriter(f)
	err = png.Encode(w, img)
	if err != nil {
		return err
	}
	err = w.Flush()
	if err != nil {
		return err
	}
	err = f.Close()
	if err != nil {
		return err
	}
	return nil
}

func paintBackground(img *image.RGBA) {
	for x := img.Bounds().Min.X; x < img.Bounds().Max.X; x++ {
		for y := img.Bounds().Min.Y; y < img.Bounds().Max.Y; y++ {
			img.Set(x, y, color.RGBA{0x00, 0x00, 0x00, 0xff})
		}
	}
}

func paintBoxes(img *image.RGBA, boxes []cdf.Box) {
	w, h := img.Bounds().Dx(), img.Bounds().Dy()
	d := density(w, h, boxes)
	max := maxDensity(d)
	for _, b := range boxes {
		for x := int(b.X.Min); x < int(b.X.Max); x++ {
			for y := int(b.Y.Min); y < int(b.Y.Max); y++ {
				c := uint8(0xff * d[y*w+x] * d[y*w+x] / (max * max))
				img.Set(x, y, color.RGBA{c, c, c, 0xff})
			}
		}
	}
	for _, b := range boxes {
		paintBoxOutline(img, b)
	}
}

func density(w, h int, boxes []cdf.Box) []float64 {
	d := make([]float64, w*h)
	for _, b := range boxes {
		for x := int(b.X.Min); x < int(b.X.Max); x++ {
			for y := int(b.Y.Min); y < int(b.Y.Max); y++ {
				d[y*w+x] += b.Weight / b.Area()
			}
		}
	}
	return d
}

func maxDensity(d []float64) float64 {
	max := -math.MaxFloat64
	for _, v := range d {
		if v > max {
			max = v
		}
	}
	return max
}

func paintBoxOutline(img *image.RGBA, b cdf.Box) {
	for x := int(b.X.Min); x < int(b.X.Max); x++ {
		img.Set(x, int(b.Y.Min), color.RGBA{0x00, 0x00, 0xff, 0xff})
		img.Set(x, int(b.Y.Max), color.RGBA{0x00, 0x00, 0xff, 0xff})
	}
	for y := int(b.Y.Min); y < int(b.Y.Max); y++ {
		img.Set(int(b.X.Min), y, color.RGBA{0x00, 0x00, 0xff, 0xff})
		img.Set(int(b.X.Max), y, color.RGBA{0x00, 0x00, 0xff, 0xff})
	}
}

func paintPoints(img *image.RGBA, pts []cdf.Pt) {
	for _, p := range pts {
		x := int(p.X)
		y := int(p.Y)
		img.Set(x, y, color.RGBA{0xff, 0xff, 0xff, 0xff})
	}
}

func printDensities(boxes []cdf.Box) {
	for _, b := range boxes {
		fmt.Printf("Weight %f, Area %f, Density %f\n", b.Weight, b.Area(), b.Weight/b.Area())
	}
}
