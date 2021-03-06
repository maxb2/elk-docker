(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

# Minimum and maximum functions that I can differentiate

m_min := (x1, x2) -> my_piecewise3(x1 > x2, x2, x1):
m_max := (x1, x2) -> my_piecewise3(x1 > x2, x1, x2):
m_abs := (x)      -> my_piecewise3(x > 0, x, -x):

# Teach maple to differentiate piecewise functions

`diff/my_piecewise3` :=
    proc(c, x1, x2, x) my_piecewise3(c, diff(x1, x), diff(x2, x)) end proc:

# This is the derivative of xc_E1_scaled = -exp(x)*Ei(-x) = exp(x)E1(x)

`diff/xc_E1_scaled` :=
  proc(y, x) (xc_E1_scaled(y) - 1/y) * diff(y, x) end proc:

# The derivative of xc_erfcx = exp(x^2)*erfc(x) = 2*x*xc_erfcx(x) - 2/sqrt(Pi)
`diff/xc_erfcx` :=
  proc(y, x) (2*y*xc_erfcx(y) - 2/sqrt(Pi)) * diff(y, x) end proc:


# a series of useful definitions

M_C         := 137.0359996287515: (* speed of light *)

X2S         := 1/(2*(6*Pi^2)^(1/3)):
X2S_2D      := 1/(2*(4*Pi)^(1/2)):

X_FACTOR_C    := 3/8*(3/Pi)^(1/3)*4^(2/3):
X_FACTOR_2D_C := 8/(3*sqrt(Pi)):
K_FACTOR_C    := 3/10*(6*Pi^2)^(2/3):

MU_GE       := 10/81:
MU_PBE      := 0.06672455060314922*(Pi^2)/3:
KAPPA_PBE   := 0.8040:

# generic conversion functions

$ifdef xc_dimensions_1d
DIMENSIONS   := 1:
RS_FACTOR    := 1/2:
$elif xc_dimensions_2d
DIMENSIONS   := 2:
RS_FACTOR    := 1/sqrt(Pi):
LDA_X_FACTOR := -X_FACTOR_2D_C:
$else
DIMENSIONS   := 3:
RS_FACTOR    := (3/(4*Pi))^(1/3):
LDA_X_FACTOR := -X_FACTOR_C:
$endif

r_ws       := n  -> RS_FACTOR/n^(1/DIMENSIONS):
n_total    := rs -> (RS_FACTOR/rs)^DIMENSIONS:

n_spin     := (rs, z) -> simplify((1 + z)*n_total(rs)/2):
sigma_spin := (rs, z, xs) -> simplify(xs^2*n_spin(rs, z)^(8/3)):
t_total    := (z, ts0, ts1) ->
  (ts0*((1 + z)/2)^(5/3) + ts1*((1 - z)/2)^(5/3)):
u_total    := (z, us0, us1) -> t_total(z, us0, us1):

# useful formulas that enter several functionals follow

# von Weizsäcker term
t_vw := (z, xt, us0, us1) -> (xt^2 - u_total(z, us0, us1))/8:

# See Eq. (9) of Perdew1992_13244
f_zeta    := z -> ((1 + z)^(4/3) + (1 - z)^(4/3) - 2)/(2^(4/3) - 2):
f_zeta_2d := z -> 1/2*((1 + z)^(3/2) + (1 - z)^(3/2)):

# used in several correlation functionals
mphi := z -> ((1 + z)^(2/3) + (1 - z)^(2/3))/2:
tt   := (rs, z, xt) -> xt/(4*2^(1/3)*mphi(z)*sqrt(rs)):

# in the paper it is beta_a = 0.066725
beta_Hu_Langreth := rs -> 0.066724550603149220*(1 + 0.1*rs)/(1 + 0.1778*rs):

# Generate exchange and kinetic functionals from the expression for the
# enhancement factor
lda_x_spin := (rs, z) -> LDA_X_FACTOR*((1 + z)/2)^(1 + 1/DIMENSIONS)*(RS_FACTOR/rs):
lda_k_spin := (rs, z) -> K_FACTOR_C*((1 + z)/2)^(5/3)*(RS_FACTOR/rs)^2:
if evalb(Polarization = "ferr") then
    gga_exchange_nsp := (func, rs, z, xs0, xs1) ->
             lda_x_spin(rs, 1)*func(rs, 1, xs0):
    gga_exchange := (func, rs, z, xs0, xs1) ->
             lda_x_spin(rs, 1)*func(xs0):
    gga_kinetic := (func, rs, z, xs0, xs1) ->
             lda_k_spin(rs, 1)*func(xs0):

    mgga_exchange_nsp := (func, rs, z, xs0, xs1, u0, u1, t0, t1) ->
             lda_x_spin(rs, 1)*func(rs, 1, xs0, u0, t0):
    mgga_exchange := (func, rs, z, xs0, xs1, u0, u1, t0, t1) ->
             lda_x_spin(rs, 1)*func(xs0, u0, t0):
    mgga_kinetic := (func, rs, z, xs0, xs1, u0, u1) ->
             lda_k_spin(rs, 1)*func(xs0, u0):
else
    gga_exchange_nsp := (func, rs, z, xs0, xs1) ->
             lda_x_spin(rs, z)*func(rs, z, xs0) + lda_x_spin(rs, -z)*func(rs, -z, xs1):
    gga_exchange := (func, rs, z, xs0, xs1) ->
             lda_x_spin(rs, z)*func(xs0) + lda_x_spin(rs, -z)*func(xs1):
    gga_kinetic := (func, rs, z, xs0, xs1) ->
             lda_k_spin(rs, z)*func(xs0) + lda_k_spin(rs, -z)*func(xs1):

    mgga_exchange_nsp := (func, rs, z, xs0, xs1, u0, u1, t0, t1) ->
             lda_x_spin(rs, z)*func(rs, z, xs0, u0, t0) + lda_x_spin(rs, -z)*func(rs, -z, xs1, u1, t1):
    mgga_exchange := (func, rs, z, xs0, xs1, u0, u1, t0, t1) ->
             lda_x_spin(rs, z)*func(xs0, u0, t0) + lda_x_spin(rs, -z)*func(xs1, u1, t1):
    mgga_kinetic := (func, rs, z, xs0, xs1, u0, u1) ->
             lda_k_spin(rs, z)*func(xs0, u0) + lda_k_spin(rs, -z)*func(xs1, u1):
end if:

# This is the Stoll decomposition in our language
lda_stoll_par  := (lda_func, rs, z) ->
  (1 + z)/2 * lda_func(rs*(2/(1 + z))^(1/3), 1):

lda_stoll_perp := (lda_func, rs, z) ->
  + lda_func(rs, z)
  - lda_stoll_par(lda_func, rs,  z)
  - lda_stoll_par(lda_func, rs, -z):

gga_stoll_par  := (gga_func, rs, z, xs, spin) ->
 + gga_func(rs*(2/(1 + z))^(1/3), spin, xs, xs*(1 + spin)/2, xs*(1 - spin)/2)
 * (1 + z)/2:

# Curvature of the Fermi hole (without the current term)
Fermi_D := (xs, ts) -> 1 - xs^2/(8*ts):

# correction to Fermi_D similar to the one found in
#  JCP 127, 214103 (2007); doi: http://dx.doi.org/10.1063/1.2800011
Fermi_D_corrected := (xs, ts) -> (1 - xs^2/(8*ts)) * (1 - exp(-4*ts^2/params_a_Fermi_D_cnst^2)):

# Becke function used in several correlation functionals
b88_R_F := (f_x, rs, z, xs) ->
  1/(2*X_FACTOR_C*n_spin(rs, z)^(1/3)*f_x(xs)):

b88_zss := (css, f_x, rs, z, xs) -> 2*css*b88_R_F(f_x, rs, z, xs):

if evalb(Polarization = "ferr") then
  b88_zab := (cab, f_x, rs, z, xs0, xs1) -> cab*b88_R_F(f_x, rs, 1, xs0):
else
  b88_zab := (cab, f_x, rs, z, xs0, xs1) ->
    cab*(b88_R_F(f_x, rs, z, xs0) + b88_R_F(f_x, rs, -z, xs1)):
end if:

# Power series often used in mggas
mgga_w := t -> (K_FACTOR_C - t)/(K_FACTOR_C + t):
mgga_series_w := (a, n, t) -> add(a[i]*mgga_w(t)^(i-1), i=1..n):

# Used in screened functionals
kF := (rs, z) -> (3*Pi^2*(1 + z))^(1/3) * RS_FACTOR/rs:
nu := (rs, z) -> p_a_cam_omega/kF(rs, z):
