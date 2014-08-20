/*
 * Adapted from http://kluge.in-chemnitz.de/opensource/spline/spline.h
 *
 * simple cubic spline interpolation library without external
 * dependencies
 *
 * ---------------------------------------------------------------------
 * Copyright (C) 2014 Pieter Kelchtermans (pieter at kunselhof.be)
 *               2011, 2014 Tino Kluge (ttk448 at gmail.com)
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ---------------------------------------------------------------------
 *
 */

package spline

import (
	"log"
	"math"
	"sort"
)

type band_matrix struct {
	m_upper [][]float64
	m_lower [][]float64
}

func make_band_matrix(dim int, n_u int, n_l int) band_matrix {
	bm := band_matrix{}
	bm.resize(dim, n_u, n_l)
	return bm
}

func (bm *band_matrix) resize(dim int, n_u int, n_l int) {
	if dim <= 0 || n_u < 0 || n_l < 0 {
		log.Fatal()
	}
	bm.m_upper = append(bm.m_upper, make([][]float64, n_u+1)...)
	bm.m_lower = append(bm.m_lower, make([][]float64, n_u+1)...)

	for i := 0; i < len(bm.m_upper); i++ {
		bm.m_upper[i] = append(bm.m_upper[i], make([]float64, dim)...)
	}
	for i := 0; i < len(bm.m_lower); i++ {
		bm.m_lower[i] = append(bm.m_lower[i], make([]float64, dim)...)
	}
}

func (bm *band_matrix) dim() int {
	if len(bm.m_upper) > 0 {
		return len(bm.m_upper[0])
	} else {
		return 0
	}
}

func (bm *band_matrix) num_upper() int {
	return len(bm.m_upper) - 1
}

func (bm *band_matrix) num_lower() int {
	return len(bm.m_lower) - 1
}

// defines the new operator (), so that we can access the elements
// by A(i,j), index going from i=0,...,dim()-1
func (bm *band_matrix) operator(i int, j int) float64 {
	k := j - i // what band is the entry
	if !((i >= 0) && (i < bm.dim()) && (j >= 0) && (j < bm.dim())) {
		log.Fatal()
	}
	if !((-bm.num_lower() <= k) && (k <= bm.num_upper())) {
		log.Fatal()
	}
	// k=0 -> diogonal, k<0 lower left part, k>0 upper right part
	if k >= 0 {
		return bm.m_upper[k][i]
	} else {
		return bm.m_lower[-k][i]
	}
}

func (bm *band_matrix) setOperator(i int, j int, value float64) {
	k := j - i // what band is the entry
	if !((i >= 0) && (i < bm.dim()) && (j >= 0) && (j < bm.dim())) {
		log.Fatal()
	}
	if !((-bm.num_lower() <= k) && (k <= bm.num_upper())) {
		log.Fatal()
	}
	// k=0 -> diogonal, k<0 lower left part, k>0 upper right part
	if k >= 0 {
		bm.m_upper[k][i] = value
	} else {
		bm.m_lower[-k][i] = value
	}
}

// second diag (used in LU decomposition), saved in m_lower
func (bm *band_matrix) saved_diag(i int) float64 {
	if !((i >= 0) && (i < bm.dim())) {
		log.Fatal()
	}
	return bm.m_lower[0][i]
}

func (bm *band_matrix) setSavedDiag(i int, value float64) {
	if !((i >= 0) && (i < bm.dim())) {
		log.Fatal()
	}
	bm.m_lower[0][i] = value
}

// LR-Decomposition of a band matrix
func (bm *band_matrix) lu_decompose() {
	var i_max, j_max int
	var j_min int
	var x float64

	// preconditioning
	// normalize column i so that a_ii=1
	for i := 0; i < bm.dim(); i++ {
		if !(bm.operator(i, i) != 0.0) {
			log.Fatal()
		}
		bm.setSavedDiag(i, 1.0/bm.operator(i, i))
		j_min = int(math.Max(0, float64(i-bm.num_lower())))
		j_max = int(math.Min(float64(bm.dim()-1), float64(i+bm.num_upper())))
		for j := j_min; j <= j_max; j++ {
			bm.setOperator(i, j, bm.operator(i, j)*bm.saved_diag(i))
		}
		bm.setOperator(i, i, 1.0) // prevents rounding errors
	}

	// Gauss LR-Decomposition
	for k := 0; k < bm.dim(); k++ {
		i_max = int(math.Min(float64(bm.dim()-1), float64(k+bm.num_lower()))) // num_lower not a mistake!
		for i := k + 1; i <= i_max; i++ {
			if !(bm.operator(k, k) != 0.0) {
				log.Fatal()
			}
			x = -bm.operator(i, k) / bm.operator(k, k)
			bm.setOperator(i, k, -x) // assembly part of L
			j_max = int(math.Min(float64(bm.dim()-1), float64(k+bm.num_upper())))
			for j := k + 1; j <= j_max; j++ {
				// assembly part of R
				bm.setOperator(i, j, bm.operator(i, j)+x*bm.operator(k, j))
			}
		}
	}
}

// solves Rx=y
func (bm *band_matrix) r_solve(b []float64) []float64 {
	if !(bm.dim() == len(b)) {
		log.Fatal()
	}
	x := make([]float64, bm.dim())
	var j_stop int
	var sum float64
	for i := bm.dim() - 1; i >= 0; i-- {
		sum = 0
		j_stop = int(math.Min(float64(bm.dim()-1), float64(i+bm.num_upper())))
		for j := i + 1; j <= j_stop; j++ {
			sum += bm.operator(i, j) * x[j]
		}
		x[i] = (b[i] - sum) / bm.operator(i, i)
	}
	return x
}

// solves Ly=b
func (bm *band_matrix) l_solve(b []float64) []float64 {
	if !(bm.dim() == len(b)) {
		log.Fatal()
	}
	x := make([]float64, bm.dim())
	var j_start int
	var sum float64
	for i := 0; i < bm.dim(); i++ {
		sum = 0
		j_start = int(math.Max(0, float64(i-bm.num_lower())))
		for j := j_start; j < i; j++ {
			sum += bm.operator(i, j) * x[j]
		}
		x[i] = (b[i] * bm.saved_diag(i)) - sum
	}
	return x
}

func (bm *band_matrix) lu_solve(b []float64, is_lu_decomposed bool) []float64 {
	if !(bm.dim() == len(b)) {
		log.Fatal()
	}
	var x, y []float64
	if is_lu_decomposed == false {
		bm.lu_decompose()
	}
	y = bm.l_solve(b)
	x = bm.r_solve(y)
	return x
}

type Spline struct {
	m_x, m_y           []float64
	m_a, m_b, m_c, m_d []float64
}

func (s *Spline) Set_points(x []float64, y []float64, cubic_spline bool) {
	if !(len(x) == len(y)) {
		log.Fatal()
	}
	s.m_x = x
	s.m_y = y
	n := len(x)

	// TODO sort x and y, rather than returning an error
	for i := 0; i < n-1; i++ {
		if !(s.m_x[i] < s.m_x[i+1]) {
			log.Fatal("input not sorted")
		}
	}

	if cubic_spline == true { // cubic spline interpolation
		// setting up the matrix and right hand side of the equation system
		// for the parameters b[]
		A := make_band_matrix(n, 1, 1)
		rhs := make([]float64, n)
		for i := 1; i < n-1; i++ {
			A.setOperator(i, i-1, 1.0/3.0*(x[i]-x[i-1]))
			A.setOperator(i, i, 2.0/3.0*(x[i+1]-x[i-1]))
			A.setOperator(i, i+1, 1.0/3.0*(x[i+1]-x[i]))
			rhs[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1])
		}
		// boundary conditions, zero curvature b[0]=b[n-1]=0
		A.setOperator(0, 0, 2.0)
		A.setOperator(0, 1, 0.0)
		rhs[0] = 0.0
		A.setOperator(n-1, n-1, 2.0)
		A.setOperator(n-1, n-2, 0.0)
		rhs[n-1] = 0.0

		// solve the equation system to obtain the parameters b[]
		s.m_b = A.lu_solve(rhs, false)

		// calculate parameters a[] and c[] based on b[]
		s.m_a = append(s.m_a, make([]float64, n)...)
		s.m_c = append(s.m_c, make([]float64, n)...)
		for i := 0; i < n-1; i++ {
			s.m_a[i] = 1.0 / 3.0 * (s.m_b[i+1] - s.m_b[i]) / (x[i+1] - x[i])
			s.m_c[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) -
				1.0/3.0*(2.0*s.m_b[i]+s.m_b[i+1])*(x[i+1]-x[i])
		}
	} else { // linear interpolation
		s.m_a = append(s.m_a, make([]float64, n)...)
		s.m_b = append(s.m_b, make([]float64, n)...)
		s.m_c = append(s.m_c, make([]float64, n)...)
		for i := 0; i < n-1; i++ {
			s.m_a[i] = 0.0
			s.m_b[i] = 0.0
			s.m_c[i] = (s.m_y[i+1] - s.m_y[i]) / (s.m_x[i+1] - s.m_x[i])
		}
	}

	// for the right boundary we define
	// f_{n-1}(x) = b*(x-x_{n-1})^2 + c*(x-x_{n-1}) + y_{n-1}
	h := x[n-1] - x[n-2]
	// m_b[n-1] is determined by the boundary condition
	s.m_a[n-1] = 0.0
	s.m_c[n-1] = 3.0*s.m_a[n-2]*h*h + 2.0*s.m_b[n-2]*h + s.m_c[n-2] // = f'_{n-2}(x_{n-1})
}

func (s *Spline) Operate(x float64) float64 {
	n := len(s.m_x)
	// find the closest point m_x[idx] < x, idx=0 even if x<m_x[0]
	if !sort.Float64sAreSorted(s.m_x) {
		log.Fatal("not sorted...")
	}
	it := sort.SearchFloat64s(s.m_x, x)
	idx := int(math.Max(float64(it-1), 0))

	h := x - s.m_x[idx]
	var interpol float64
	if x < s.m_x[0] {
		// extrapolation to the left
		interpol = ((s.m_b[0])*h+s.m_c[0])*h + s.m_y[0]
	} else if x > s.m_x[n-1] {
		// extrapolation to the right
		interpol = ((s.m_b[n-1])*h+s.m_c[n-1])*h + s.m_y[n-1]
	} else {
		// interpolation
		interpol = ((s.m_a[idx]*h+s.m_b[idx])*h+s.m_c[idx])*h + s.m_y[idx]
	}
	return interpol
}
