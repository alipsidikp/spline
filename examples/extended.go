package main

import (
	"fmt"
	"github.com/pkelchte/spline"
)

func main() {
	X := []float64{0.1, 0.4, 1.2, 1.8, 2.0}
	Y := []float64{0.1, 0.7, 0.6, 1.1, 0.9}

	s := spline.Spline{}

	s.Set_points(X, Y, true)

	for i := 0; i < len(X); i++ {
		fmt.Printf("%f %f\n", X[i], Y[i])
	}
	fmt.Printf("\n")
	for i := -50; i < 250; i++ {
		x := 0.01 * float64(i)
		fmt.Printf("%f %f\n", x, s.Operate(x))
	}
	fmt.Printf("\n")

}
