package main

import (
	"fmt"
	"spline"
)

func main() {
	X := []float64{0.1, 0.4, 1.2, 1.8, 2.0}
	Y := []float64{0.1, 0.7, 0.6, 1.1, 0.9}

	s := spline.Spline{}

	s.Set_points(X, Y, true)

	x := 1.5

	fmt.Printf("spline at %f is %f\n", x, s.Operate(x))

}
