/* 
  This file was generated automatically with ./scripts/maple2c.pl.
  Do not edit this file directly as it can be overwritten!!

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.

  Maple version     : Maple 2016 (X86 64 LINUX)
  Maple source      : ./maple/lda_exc/lda_c_pz.mpl
  Type of functional: lda_exc
*/

#define maple2c_order 4
#define MAPLE2C_FLAGS (XC_FLAGS_I_HAVE_EXC | XC_FLAGS_I_HAVE_VXC | XC_FLAGS_I_HAVE_FXC | XC_FLAGS_I_HAVE_KXC | XC_FLAGS_I_HAVE_LXC)


static inline void
func_unpol(const xc_func_type *p, int order, const double *rho, double *zk, LDA_OUT_PARAMS_NO_EXC(double *))
{

#ifndef XC_DONT_COMPILE_EXC
  double t1, t2, t3, t5, t6, t7, t8, t9;
  double t10, t11, t12, t13, t14, t15, t19, t20;
  double t21, t24, t27, t28, t32, t33, t38;

#ifndef XC_DONT_COMPILE_VXC
  double t42, t44, t47, t49, t50, t54, t56, t68;

#ifndef XC_DONT_COMPILE_FXC
  double t73, t74, t80, t81, t82, t83, t84, t85;
  double t92, t93, t98, t113;

#ifndef XC_DONT_COMPILE_KXC
  double t116, t118, t132, t133, t134, t138, t145, t146;
  double t151, t166;

#ifndef XC_DONT_COMPILE_LXC
  double t172, t178, t190, t207, t208, t227;
#endif

#endif

#endif

#endif

#endif


  lda_c_pz_params *params;

  assert(p->params != NULL);
  params = (lda_c_pz_params * )(p->params);

  t1 = M_CBRT3;
  t2 = 0.1e1 / M_PI;
  t3 = POW_1_3(t2);
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = POW_1_3(rho[0]);
  t8 = 0.1e1 / t7;
  t9 = t6 * t8;
  t10 = t1 * t3 * t9;
  t11 = t10 / 0.4e1;
  t12 = 0.1e1 <= t11;
  t13 = params->gamma[0];
  t14 = params->beta1[0];
  t15 = sqrt(t10);
  t19 = params->beta2[0] * t1;
  t20 = t3 * t6;
  t21 = t20 * t8;
  t24 = 0.1e1 + t14 * t15 / 0.2e1 + t19 * t21 / 0.4e1;
  t27 = params->a[0];
  t28 = log(t11);
  t32 = params->c[0] * t1;
  t33 = t32 * t3;
  t38 = params->d[0] * t1;
  if(zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    zk[0] = my_piecewise3(t12, t13 / t24, t27 * t28 + params->b[0] + t33 * t9 * t28 / 0.4e1 + t38 * t21 / 0.4e1);

#ifndef XC_DONT_COMPILE_VXC

  if(order < 1) return;


  t42 = t24 * t24;
  t44 = t13 / t42;
  t47 = t14 / t15 * t1;
  t49 = 0.1e1 / t7 / rho[0];
  t50 = t20 * t49;
  t54 = -t19 * t50 / 0.12e2 - t47 * t50 / 0.12e2;
  t56 = 0.1e1 / rho[0];
  t68 = my_piecewise3(t12, -t44 * t54, -t27 * t56 / 0.3e1 - t33 * t6 * t49 * t28 / 0.12e2 - t32 * t50 / 0.12e2 - t38 * t50 / 0.12e2);
  if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    vrho[0] = rho[0] * t68 + (my_piecewise3(t12, t13 / t24, t27 * t28 + params->b[0] + t33 * t9 * t28 / 0.4e1 + t38 * t21 / 0.4e1));

#ifndef XC_DONT_COMPILE_FXC

  if(order < 2) return;


  t73 = t13 / t42 / t24;
  t74 = t54 * t54;
  t80 = t1 * t1;
  t81 = t14 / t15 / t10 * t80;
  t82 = t3 * t3;
  t83 = t82 * t5;
  t84 = rho[0] * rho[0];
  t85 = t7 * t7;
  t92 = 0.1e1 / t7 / t84;
  t93 = t20 * t92;
  t98 = -t81 * t83 / t85 / t84 / 0.18e2 + t47 * t93 / 0.9e1 + t19 * t93 / 0.9e1;
  t113 = my_piecewise3(t12, -t44 * t98 + 0.2e1 * t73 * t74, t27 / t84 / 0.3e1 + t33 * t6 * t92 * t28 / 0.9e1 + 0.5e1 / 0.36e2 * t32 * t93 + t38 * t93 / 0.9e1);
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[0] = rho[0] * t113 + 0.2e1 * t68;

#ifndef XC_DONT_COMPILE_KXC

  if(order < 3) return;


  t116 = t42 * t42;
  t118 = t13 / t116;
  t132 = t14 / t15 / t80 / t82 / t5 * t85 / 0.4e1;
  t133 = t84 * t84;
  t134 = 0.1e1 / t133;
  t138 = t84 * rho[0];
  t145 = 0.1e1 / t7 / t138;
  t146 = t20 * t145;
  t151 = -t132 * t2 * t134 / 0.3e1 + 0.2e1 / 0.9e1 * t81 * t83 / t85 / t138 - 0.7e1 / 0.27e2 * t47 * t146 - 0.7e1 / 0.27e2 * t19 * t146;
  t166 = my_piecewise3(t12, -0.6e1 * t118 * t74 * t54 + 0.6e1 * t73 * t54 * t98 - t44 * t151, -0.2e1 / 0.3e1 * t27 / t138 - 0.7e1 / 0.27e2 * t33 * t6 * t145 * t28 - 0.13e2 / 0.36e2 * t32 * t146 - 0.7e1 / 0.27e2 * t38 * t146);
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[0] = rho[0] * t166 + 0.3e1 * t113;

#ifndef XC_DONT_COMPILE_LXC

  if(order < 4) return;


  t172 = t74 * t74;
  t178 = t98 * t98;
  t190 = t133 * rho[0];
  t207 = 0.1e1 / t7 / t133;
  t208 = t20 * t207;
  t227 = my_piecewise3(t12, 0.24e2 * t13 / t116 / t24 * t172 - 0.36e2 * t118 * t74 * t98 + 0.6e1 * t73 * t178 + 0.8e1 * t73 * t54 * t151 - t44 * (-0.5e1 / 0.864e3 * t14 / t15 / t56 / t7 / t190 * t1 * t20 + 0.8e1 / 0.3e1 * t132 * t2 / t190 - 0.80e2 / 0.81e2 * t81 * t83 / t85 / t133 + 0.70e2 / 0.81e2 * t47 * t208 + 0.70e2 / 0.81e2 * t19 * t208), 0.2e1 * t27 * t134 + 0.70e2 / 0.81e2 * t33 * t6 * t207 * t28 + 0.209e3 / 0.162e3 * t32 * t208 + 0.70e2 / 0.81e2 * t38 * t208);
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[0] = rho[0] * t227 + 0.4e1 * t166;

#ifndef XC_DONT_COMPILE_MXC

  if(order < 5) return;


#endif

#endif

#endif

#endif

#endif


}


static inline void
func_ferr(const xc_func_type *p, int order, const double *rho, double *zk, LDA_OUT_PARAMS_NO_EXC(double *))
{

#ifndef XC_DONT_COMPILE_EXC
  double t1, t2, t3, t5, t6, t7, t8, t9;
  double t10, t11, t12, t13, t14, t15, t19, t20;
  double t21, t24, t27, t28, t32, t33, t38;

#ifndef XC_DONT_COMPILE_VXC
  double t42, t44, t47, t49, t50, t54, t56, t68;

#ifndef XC_DONT_COMPILE_FXC
  double t73, t74, t80, t81, t82, t83, t84, t85;
  double t92, t93, t98, t113;

#ifndef XC_DONT_COMPILE_KXC
  double t116, t118, t132, t133, t134, t138, t145, t146;
  double t151, t166;

#ifndef XC_DONT_COMPILE_LXC
  double t172, t178, t190, t207, t208, t227;
#endif

#endif

#endif

#endif

#endif


  lda_c_pz_params *params;

  assert(p->params != NULL);
  params = (lda_c_pz_params * )(p->params);

  t1 = M_CBRT3;
  t2 = 0.1e1 / M_PI;
  t3 = POW_1_3(t2);
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = POW_1_3(rho[0]);
  t8 = 0.1e1 / t7;
  t9 = t6 * t8;
  t10 = t1 * t3 * t9;
  t11 = t10 / 0.4e1;
  t12 = 0.1e1 <= t11;
  t13 = params->gamma[1];
  t14 = params->beta1[1];
  t15 = sqrt(t10);
  t19 = params->beta2[1] * t1;
  t20 = t3 * t6;
  t21 = t20 * t8;
  t24 = 0.1e1 + t14 * t15 / 0.2e1 + t19 * t21 / 0.4e1;
  t27 = params->a[1];
  t28 = log(t11);
  t32 = params->c[1] * t1;
  t33 = t32 * t3;
  t38 = params->d[1] * t1;
  if(zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    zk[0] = my_piecewise3(t12, t13 / t24, t27 * t28 + params->b[1] + t33 * t9 * t28 / 0.4e1 + t38 * t21 / 0.4e1);

#ifndef XC_DONT_COMPILE_VXC

  if(order < 1) return;


  t42 = t24 * t24;
  t44 = t13 / t42;
  t47 = t14 / t15 * t1;
  t49 = 0.1e1 / t7 / rho[0];
  t50 = t20 * t49;
  t54 = -t19 * t50 / 0.12e2 - t47 * t50 / 0.12e2;
  t56 = 0.1e1 / rho[0];
  t68 = my_piecewise3(t12, -t44 * t54, -t27 * t56 / 0.3e1 - t33 * t6 * t49 * t28 / 0.12e2 - t32 * t50 / 0.12e2 - t38 * t50 / 0.12e2);
  if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    vrho[0] = rho[0] * t68 + (my_piecewise3(t12, t13 / t24, t27 * t28 + params->b[1] + t33 * t9 * t28 / 0.4e1 + t38 * t21 / 0.4e1));

#ifndef XC_DONT_COMPILE_FXC

  if(order < 2) return;


  t73 = t13 / t42 / t24;
  t74 = t54 * t54;
  t80 = t1 * t1;
  t81 = t14 / t15 / t10 * t80;
  t82 = t3 * t3;
  t83 = t82 * t5;
  t84 = rho[0] * rho[0];
  t85 = t7 * t7;
  t92 = 0.1e1 / t7 / t84;
  t93 = t20 * t92;
  t98 = -t81 * t83 / t85 / t84 / 0.18e2 + t47 * t93 / 0.9e1 + t19 * t93 / 0.9e1;
  t113 = my_piecewise3(t12, -t44 * t98 + 0.2e1 * t73 * t74, t27 / t84 / 0.3e1 + t33 * t6 * t92 * t28 / 0.9e1 + 0.5e1 / 0.36e2 * t32 * t93 + t38 * t93 / 0.9e1);
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[0] = rho[0] * t113 + 0.2e1 * t68;

#ifndef XC_DONT_COMPILE_KXC

  if(order < 3) return;


  t116 = t42 * t42;
  t118 = t13 / t116;
  t132 = t14 / t15 / t80 / t82 / t5 * t85 / 0.4e1;
  t133 = t84 * t84;
  t134 = 0.1e1 / t133;
  t138 = t84 * rho[0];
  t145 = 0.1e1 / t7 / t138;
  t146 = t20 * t145;
  t151 = -t132 * t2 * t134 / 0.3e1 + 0.2e1 / 0.9e1 * t81 * t83 / t85 / t138 - 0.7e1 / 0.27e2 * t47 * t146 - 0.7e1 / 0.27e2 * t19 * t146;
  t166 = my_piecewise3(t12, -0.6e1 * t118 * t74 * t54 + 0.6e1 * t73 * t54 * t98 - t44 * t151, -0.2e1 / 0.3e1 * t27 / t138 - 0.7e1 / 0.27e2 * t33 * t6 * t145 * t28 - 0.13e2 / 0.36e2 * t32 * t146 - 0.7e1 / 0.27e2 * t38 * t146);
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[0] = rho[0] * t166 + 0.3e1 * t113;

#ifndef XC_DONT_COMPILE_LXC

  if(order < 4) return;


  t172 = t74 * t74;
  t178 = t98 * t98;
  t190 = t133 * rho[0];
  t207 = 0.1e1 / t7 / t133;
  t208 = t20 * t207;
  t227 = my_piecewise3(t12, 0.24e2 * t13 / t116 / t24 * t172 - 0.36e2 * t118 * t74 * t98 + 0.6e1 * t73 * t178 + 0.8e1 * t73 * t54 * t151 - t44 * (-0.5e1 / 0.864e3 * t14 / t15 / t56 / t7 / t190 * t1 * t20 + 0.8e1 / 0.3e1 * t132 * t2 / t190 - 0.80e2 / 0.81e2 * t81 * t83 / t85 / t133 + 0.70e2 / 0.81e2 * t47 * t208 + 0.70e2 / 0.81e2 * t19 * t208), 0.2e1 * t27 * t134 + 0.70e2 / 0.81e2 * t33 * t6 * t207 * t28 + 0.209e3 / 0.162e3 * t32 * t208 + 0.70e2 / 0.81e2 * t38 * t208);
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[0] = rho[0] * t227 + 0.4e1 * t166;

#ifndef XC_DONT_COMPILE_MXC

  if(order < 5) return;


#endif

#endif

#endif

#endif

#endif


}


static inline void
func_pol(const xc_func_type *p, int order, const double *rho, double *zk, LDA_OUT_PARAMS_NO_EXC(double *))
{

#ifndef XC_DONT_COMPILE_EXC
  double t1, t2, t3, t5, t6, t7, t8, t9;
  double t10, t11, t12, t13, t14, t15, t16, t20;
  double t21, t22, t25, t28, t29, t33, t34, t35;
  double t39, t43, t44, t45, t49, t52, t55, t59;
  double t60, t64, t68, t69, t70, t71, t72, t73;
  double t74, t76, t77, t79, t81, t84, t85;

#ifndef XC_DONT_COMPILE_VXC
  double t86, t88, t89, t91, t93, t94, t98, t103;
  double t111, t112, t114, t116, t120, t131, t132, t134;
  double t135, t136, t137, t138, t140, t143, t145, t148;
  double t150, t153, t155;

#ifndef XC_DONT_COMPILE_FXC
  double t158, t159, t163, t164, t168, t170, t171, t172;
  double t173, t174, t177, t181, t182, t187, t193, t201;
  double t204, t205, t209, t216, t228, t229, t231, t233;
  double t234, t235, t236, t237, t240, t241, t242, t244;
  double t247, t248, t249, t252, t255, t257, t261, t262;
  double t265, t268, t271, t274, t276, t280, t281, t285;
  double t288, t291, t294, t296;

#ifndef XC_DONT_COMPILE_KXC
  double t299, t300, t303, t305, t318, t319, t320, t321;
  double t322, t327, t331, t332, t337, t343, t351, t352;
  double t354, t361, t370, t382, t383, t385, t387, t388;
  double t390, t393, t394, t397, t400, t402, t406, t407;
  double t410, t413, t416, t418, t422, t425, t427, t428;
  double t429, t432, t442, t445, t455, t457, t463, t464;
  double t469, t474, t477, t482, t485, t488, t490, t495;
  double t497, t503, t506, t511, t514, t516;

#ifndef XC_DONT_COMPILE_LXC
  double t519, t520, t527, t533, t542, t545, t549, t552;
  double t553, t558, t562, t563, t574, t582, t586, t592;
  double t622, t625, t627, t630, t633, t635, t637, t638;
  double t644, t649, t651, t654, t656, t657, t663, t678;
  double t683, t685, t686, t688, t713, t731, t740, t757;
  double t759, t770, t773, t812, t824, t843, t865, t876;
  double t881, t887, t890, t895;
#endif

#endif

#endif

#endif

#endif


  lda_c_pz_params *params;

  assert(p->params != NULL);
  params = (lda_c_pz_params * )(p->params);

  t1 = M_CBRT3;
  t2 = 0.1e1 / M_PI;
  t3 = POW_1_3(t2);
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = rho[0] + rho[1];
  t8 = POW_1_3(t7);
  t9 = 0.1e1 / t8;
  t10 = t6 * t9;
  t11 = t1 * t3 * t10;
  t12 = t11 / 0.4e1;
  t13 = 0.1e1 <= t12;
  t14 = params->gamma[0];
  t15 = params->beta1[0];
  t16 = sqrt(t11);
  t20 = params->beta2[0] * t1;
  t21 = t3 * t6;
  t22 = t21 * t9;
  t25 = 0.1e1 + t15 * t16 / 0.2e1 + t20 * t22 / 0.4e1;
  t28 = params->a[0];
  t29 = log(t12);
  t33 = params->c[0] * t1;
  t34 = t33 * t3;
  t35 = t10 * t29;
  t39 = params->d[0] * t1;
  t43 = my_piecewise3(t13, t14 / t25, t28 * t29 + params->b[0] + t34 * t35 / 0.4e1 + t39 * t22 / 0.4e1);
  t44 = params->gamma[1];
  t45 = params->beta1[1];
  t49 = params->beta2[1] * t1;
  t52 = 0.1e1 + t45 * t16 / 0.2e1 + t49 * t22 / 0.4e1;
  t55 = params->a[1];
  t59 = params->c[1] * t1;
  t60 = t59 * t3;
  t64 = params->d[1] * t1;
  t68 = my_piecewise3(t13, t44 / t52, t55 * t29 + params->b[1] + t60 * t35 / 0.4e1 + t64 * t22 / 0.4e1);
  t69 = t68 - t43;
  t70 = rho[0] - rho[1];
  t71 = 0.1e1 / t7;
  t72 = t70 * t71;
  t73 = 0.1e1 + t72;
  t74 = POW_1_3(t73);
  t76 = 0.1e1 - t72;
  t77 = POW_1_3(t76);
  t79 = t74 * t73 + t77 * t76 - 0.2e1;
  t81 = M_CBRT2;
  t84 = 0.1e1 / (0.2e1 * t81 - 0.2e1);
  t85 = t69 * t79 * t84;
  if(zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    zk[0] = t43 + t85;

#ifndef XC_DONT_COMPILE_VXC

  if(order < 1) return;


  t86 = t25 * t25;
  t88 = t14 / t86;
  t89 = 0.1e1 / t16;
  t91 = t15 * t89 * t1;
  t93 = 0.1e1 / t8 / t7;
  t94 = t21 * t93;
  t98 = -t20 * t94 / 0.12e2 - t91 * t94 / 0.12e2;
  t103 = t6 * t93 * t29;
  t111 = my_piecewise3(t13, -t88 * t98, -t28 * t71 / 0.3e1 - t34 * t103 / 0.12e2 - t33 * t94 / 0.12e2 - t39 * t94 / 0.12e2);
  t112 = t52 * t52;
  t114 = t44 / t112;
  t116 = t45 * t89 * t1;
  t120 = -t116 * t94 / 0.12e2 - t49 * t94 / 0.12e2;
  t131 = my_piecewise3(t13, -t114 * t120, -t55 * t71 / 0.3e1 - t60 * t103 / 0.12e2 - t59 * t94 / 0.12e2 - t64 * t94 / 0.12e2);
  t132 = t131 - t111;
  t134 = t132 * t79 * t84;
  t135 = t7 * t7;
  t136 = 0.1e1 / t135;
  t137 = t70 * t136;
  t138 = t71 - t137;
  t140 = -t138;
  t143 = 0.4e1 / 0.3e1 * t74 * t138 + 0.4e1 / 0.3e1 * t77 * t140;
  t145 = t69 * t143 * t84;
  if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    vrho[0] = t43 + t85 + t7 * (t111 + t134 + t145);

  t148 = -t71 - t137;
  t150 = -t148;
  t153 = 0.4e1 / 0.3e1 * t74 * t148 + 0.4e1 / 0.3e1 * t77 * t150;
  t155 = t69 * t153 * t84;
  if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    vrho[1] = t43 + t85 + t7 * (t111 + t134 + t155);

#ifndef XC_DONT_COMPILE_FXC

  if(order < 2) return;


  t158 = 0.2e1 * t111;
  t159 = 0.2e1 * t134;
  t163 = t14 / t86 / t25;
  t164 = t98 * t98;
  t168 = 0.1e1 / t16 / t11;
  t170 = t1 * t1;
  t171 = t15 * t168 * t170;
  t172 = t3 * t3;
  t173 = t172 * t5;
  t174 = t8 * t8;
  t177 = t173 / t174 / t135;
  t181 = 0.1e1 / t8 / t135;
  t182 = t21 * t181;
  t187 = -t171 * t177 / 0.18e2 + t91 * t182 / 0.9e1 + t20 * t182 / 0.9e1;
  t193 = t6 * t181 * t29;
  t201 = my_piecewise3(t13, 0.2e1 * t163 * t164 - t88 * t187, t28 * t136 / 0.3e1 + t34 * t193 / 0.9e1 + 0.5e1 / 0.36e2 * t33 * t182 + t39 * t182 / 0.9e1);
  t204 = t44 / t112 / t52;
  t205 = t120 * t120;
  t209 = t45 * t168 * t170;
  t216 = -t209 * t177 / 0.18e2 + t116 * t182 / 0.9e1 + t49 * t182 / 0.9e1;
  t228 = my_piecewise3(t13, -t114 * t216 + 0.2e1 * t204 * t205, t55 * t136 / 0.3e1 + t60 * t193 / 0.9e1 + 0.5e1 / 0.36e2 * t59 * t182 + t64 * t182 / 0.9e1);
  t229 = t228 - t201;
  t231 = t229 * t79 * t84;
  t233 = t132 * t143 * t84;
  t234 = 0.2e1 * t233;
  t235 = t74 * t74;
  t236 = 0.1e1 / t235;
  t237 = t138 * t138;
  t240 = t135 * t7;
  t241 = 0.1e1 / t240;
  t242 = t70 * t241;
  t244 = -0.2e1 * t136 + 0.2e1 * t242;
  t247 = t77 * t77;
  t248 = 0.1e1 / t247;
  t249 = t140 * t140;
  t252 = -t244;
  t255 = 0.4e1 / 0.9e1 * t236 * t237 + 0.4e1 / 0.3e1 * t74 * t244 + 0.4e1 / 0.9e1 * t248 * t249 + 0.4e1 / 0.3e1 * t77 * t252;
  t257 = t69 * t255 * t84;
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[0] = t158 + t159 + 0.2e1 * t145 + t7 * (t201 + t231 + t234 + t257);

  t261 = t132 * t153 * t84;
  t262 = t236 * t148;
  t265 = t74 * t70;
  t268 = t248 * t150;
  t271 = t77 * t70;
  t274 = 0.4e1 / 0.9e1 * t262 * t138 + 0.8e1 / 0.3e1 * t265 * t241 + 0.4e1 / 0.9e1 * t268 * t140 - 0.8e1 / 0.3e1 * t271 * t241;
  t276 = t69 * t274 * t84;
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[1] = t158 + t159 + t145 + t155 + t7 * (t201 + t231 + t233 + t261 + t276);

  t280 = 0.2e1 * t261;
  t281 = t148 * t148;
  t285 = 0.2e1 * t136 + 0.2e1 * t242;
  t288 = t150 * t150;
  t291 = -t285;
  t294 = 0.4e1 / 0.9e1 * t236 * t281 + 0.4e1 / 0.3e1 * t74 * t285 + 0.4e1 / 0.9e1 * t248 * t288 + 0.4e1 / 0.3e1 * t77 * t291;
  t296 = t69 * t294 * t84;
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[2] = t158 + t159 + 0.2e1 * t155 + t7 * (t201 + t231 + t280 + t296);

#ifndef XC_DONT_COMPILE_KXC

  if(order < 3) return;


  t299 = 0.3e1 * t201;
  t300 = 0.3e1 * t231;
  t303 = t86 * t86;
  t305 = t14 / t303;
  t318 = 0.1e1 / t16 / t170 / t172 / t5 * t174 / 0.4e1;
  t319 = t15 * t318;
  t320 = t135 * t135;
  t321 = 0.1e1 / t320;
  t322 = t2 * t321;
  t327 = t173 / t174 / t240;
  t331 = 0.1e1 / t8 / t240;
  t332 = t21 * t331;
  t337 = -t319 * t322 / 0.3e1 + 0.2e1 / 0.9e1 * t171 * t327 - 0.7e1 / 0.27e2 * t91 * t332 - 0.7e1 / 0.27e2 * t20 * t332;
  t343 = t6 * t331 * t29;
  t351 = my_piecewise3(t13, 0.6e1 * t163 * t98 * t187 - 0.6e1 * t305 * t164 * t98 - t88 * t337, -0.2e1 / 0.3e1 * t28 * t241 - 0.7e1 / 0.27e2 * t34 * t343 - 0.13e2 / 0.36e2 * t33 * t332 - 0.7e1 / 0.27e2 * t39 * t332);
  t352 = t112 * t112;
  t354 = t44 / t352;
  t361 = t45 * t318;
  t370 = -t361 * t322 / 0.3e1 + 0.2e1 / 0.9e1 * t209 * t327 - 0.7e1 / 0.27e2 * t116 * t332 - 0.7e1 / 0.27e2 * t49 * t332;
  t382 = my_piecewise3(t13, 0.6e1 * t204 * t120 * t216 - 0.6e1 * t354 * t205 * t120 - t114 * t370, -0.2e1 / 0.3e1 * t55 * t241 - 0.7e1 / 0.27e2 * t60 * t343 - 0.13e2 / 0.36e2 * t59 * t332 - 0.7e1 / 0.27e2 * t64 * t332);
  t383 = t382 - t351;
  t385 = t383 * t79 * t84;
  t387 = t229 * t143 * t84;
  t388 = 0.3e1 * t387;
  t390 = t132 * t255 * t84;
  t393 = 0.1e1 / t235 / t73;
  t394 = t237 * t138;
  t397 = t236 * t138;
  t400 = t70 * t321;
  t402 = 0.6e1 * t241 - 0.6e1 * t400;
  t406 = 0.1e1 / t247 / t76;
  t407 = t249 * t140;
  t410 = t248 * t140;
  t413 = -t402;
  t416 = -0.8e1 / 0.27e2 * t393 * t394 + 0.4e1 / 0.3e1 * t397 * t244 + 0.4e1 / 0.3e1 * t74 * t402 - 0.8e1 / 0.27e2 * t406 * t407 + 0.4e1 / 0.3e1 * t410 * t252 + 0.4e1 / 0.3e1 * t77 * t413;
  t418 = t69 * t416 * t84;
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[0] = t299 + t300 + 0.6e1 * t233 + 0.3e1 * t257 + t7 * (t351 + t385 + t388 + 0.3e1 * t390 + t418);

  t422 = 0.2e1 * t276;
  t425 = t229 * t153 * t84;
  t427 = t132 * t274 * t84;
  t428 = 0.2e1 * t427;
  t429 = t393 * t148;
  t432 = t236 * t70;
  t442 = t406 * t150;
  t445 = t248 * t70;
  t455 = -0.8e1 / 0.27e2 * t429 * t237 + 0.16e2 / 0.9e1 * t432 * t241 * t138 + 0.4e1 / 0.9e1 * t262 * t244 + 0.8e1 / 0.3e1 * t74 * t241 - 0.8e1 * t265 * t321 - 0.8e1 / 0.27e2 * t442 * t249 - 0.16e2 / 0.9e1 * t445 * t241 * t140 + 0.4e1 / 0.9e1 * t268 * t252 - 0.8e1 / 0.3e1 * t77 * t241 + 0.8e1 * t271 * t321;
  t457 = t69 * t455 * t84;
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[1] = t299 + t300 + 0.4e1 * t233 + t257 + t280 + t422 + t7 * (t351 + t385 + 0.2e1 * t387 + t390 + t425 + t428 + t457);

  t463 = t132 * t294 * t84;
  t464 = t393 * t281;
  t469 = t236 * t285;
  t474 = -0.2e1 * t241 - 0.6e1 * t400;
  t477 = t406 * t288;
  t482 = t248 * t291;
  t485 = -t474;
  t488 = -0.8e1 / 0.27e2 * t464 * t138 + 0.16e2 / 0.9e1 * t262 * t242 + 0.4e1 / 0.9e1 * t469 * t138 + 0.4e1 / 0.3e1 * t74 * t474 - 0.8e1 / 0.27e2 * t477 * t140 - 0.16e2 / 0.9e1 * t268 * t242 + 0.4e1 / 0.9e1 * t482 * t140 + 0.4e1 / 0.3e1 * t77 * t485;
  t490 = t69 * t488 * t84;
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[2] = t299 + t300 + t234 + 0.4e1 * t261 + t422 + t296 + t7 * (t351 + t385 + t387 + 0.2e1 * t425 + t428 + t463 + t490);

  t495 = 0.3e1 * t425;
  t497 = t281 * t148;
  t503 = -0.6e1 * t241 - 0.6e1 * t400;
  t506 = t288 * t150;
  t511 = -t503;
  t514 = -0.8e1 / 0.27e2 * t393 * t497 + 0.4e1 / 0.3e1 * t262 * t285 + 0.4e1 / 0.3e1 * t74 * t503 - 0.8e1 / 0.27e2 * t406 * t506 + 0.4e1 / 0.3e1 * t268 * t291 + 0.4e1 / 0.3e1 * t77 * t511;
  t516 = t69 * t514 * t84;
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[3] = t299 + t300 + 0.6e1 * t261 + 0.3e1 * t296 + t7 * (t351 + t385 + t495 + 0.3e1 * t463 + t516);

#ifndef XC_DONT_COMPILE_LXC

  if(order < 4) return;


  t519 = 0.4e1 * t351;
  t520 = 0.4e1 * t385;
  t527 = t164 * t164;
  t533 = t187 * t187;
  t542 = 0.1e1 / t16 / t2 / t71 / 0.48e2;
  t545 = t320 * t7;
  t549 = 0.1e1 / t8 / t545 * t1 * t21;
  t552 = 0.1e1 / t545;
  t553 = t2 * t552;
  t558 = t173 / t174 / t320;
  t562 = 0.1e1 / t8 / t320;
  t563 = t21 * t562;
  t574 = t6 * t562 * t29;
  t582 = my_piecewise3(t13, 0.24e2 * t14 / t303 / t25 * t527 - 0.36e2 * t305 * t164 * t187 + 0.6e1 * t163 * t533 + 0.8e1 * t163 * t98 * t337 - t88 * (-0.5e1 / 0.18e2 * t15 * t542 * t2 * t549 + 0.8e1 / 0.3e1 * t319 * t553 - 0.80e2 / 0.81e2 * t171 * t558 + 0.70e2 / 0.81e2 * t91 * t563 + 0.70e2 / 0.81e2 * t20 * t563), 0.2e1 * t28 * t321 + 0.70e2 / 0.81e2 * t34 * t574 + 0.209e3 / 0.162e3 * t33 * t563 + 0.70e2 / 0.81e2 * t39 * t563);
  t586 = t205 * t205;
  t592 = t216 * t216;
  t622 = my_piecewise3(t13, 0.24e2 * t44 / t352 / t52 * t586 - 0.36e2 * t354 * t205 * t216 + 0.6e1 * t204 * t592 + 0.8e1 * t204 * t120 * t370 - t114 * (-0.5e1 / 0.18e2 * t45 * t542 * t2 * t549 + 0.8e1 / 0.3e1 * t361 * t553 - 0.80e2 / 0.81e2 * t209 * t558 + 0.70e2 / 0.81e2 * t116 * t563 + 0.70e2 / 0.81e2 * t49 * t563), 0.2e1 * t55 * t321 + 0.70e2 / 0.81e2 * t60 * t574 + 0.209e3 / 0.162e3 * t59 * t563 + 0.70e2 / 0.81e2 * t64 * t563);
  t625 = (t622 - t582) * t79 * t84;
  t627 = t383 * t143 * t84;
  t630 = t229 * t255 * t84;
  t633 = t132 * t416 * t84;
  t635 = t73 * t73;
  t637 = 0.1e1 / t235 / t635;
  t638 = t237 * t237;
  t644 = t244 * t244;
  t649 = t70 * t552;
  t651 = -0.24e2 * t321 + 0.24e2 * t649;
  t654 = t76 * t76;
  t656 = 0.1e1 / t247 / t654;
  t657 = t249 * t249;
  t663 = t252 * t252;
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[0] = t519 + t520 + 0.12e2 * t387 + 0.12e2 * t390 + 0.4e1 * t418 + t7 * (t582 + t625 + 0.4e1 * t627 + 0.6e1 * t630 + 0.4e1 * t633 + t69 * (0.40e2 / 0.81e2 * t637 * t638 - 0.16e2 / 0.9e1 * t393 * t237 * t244 + 0.4e1 / 0.3e1 * t236 * t644 + 0.16e2 / 0.9e1 * t397 * t402 + 0.4e1 / 0.3e1 * t74 * t651 + 0.40e2 / 0.81e2 * t656 * t657 - 0.16e2 / 0.9e1 * t406 * t249 * t252 + 0.4e1 / 0.3e1 * t248 * t663 + 0.16e2 / 0.9e1 * t410 * t413 - 0.4e1 / 0.3e1 * t77 * t651) * t84);

  t678 = 0.6e1 * t427;
  t683 = t383 * t153 * t84;
  t685 = t229 * t274 * t84;
  t686 = 0.3e1 * t685;
  t688 = t132 * t455 * t84;
  t713 = 0.32e2 * t265 * t552;
  t731 = 0.32e2 * t271 * t552;
  t740 = -0.8e1 * t432 * t321 * t138 + 0.8e1 * t445 * t321 * t140 + 0.40e2 / 0.81e2 * t637 * t148 * t394 - 0.8e1 / 0.9e1 * t429 * t138 * t244 - 0.16e2 / 0.9e1 * t393 * t70 * t241 * t237 + 0.8e1 / 0.3e1 * t236 * t241 * t138 + 0.8e1 / 0.3e1 * t432 * t241 * t244 + t713 + 0.40e2 / 0.81e2 * t656 * t150 * t407 - 0.8e1 / 0.9e1 * t442 * t140 * t252 + 0.16e2 / 0.9e1 * t406 * t70 * t241 * t249 - 0.8e1 / 0.3e1 * t248 * t241 * t140 - 0.8e1 / 0.3e1 * t445 * t241 * t252 - t731 + 0.4e1 / 0.9e1 * t262 * t402 - 0.16e2 * t74 * t321 + 0.4e1 / 0.9e1 * t268 * t413 + 0.16e2 * t77 * t321;
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[1] = t519 + t520 + 0.9e1 * t387 + 0.6e1 * t390 + t418 + t495 + t678 + 0.3e1 * t457 + t7 * (t69 * t740 * t84 + t582 + t625 + 0.3e1 * t627 + 0.3e1 * t630 + t633 + t683 + t686 + 0.3e1 * t688);

  t757 = t229 * t294 * t84;
  t759 = t132 * t488 * t84;
  t770 = t70 * t70;
  t773 = 0.1e1 / t320 / t135;
  t812 = 0.40e2 / 0.81e2 * t637 * t281 * t237 - 0.64e2 / 0.27e2 * t429 * t138 * t70 * t241 - 0.8e1 / 0.27e2 * t464 * t244 + 0.32e2 / 0.9e1 * t236 * t770 * t773 + 0.16e2 / 0.9e1 * t262 * t241 - 0.16e2 / 0.3e1 * t262 * t400 - 0.8e1 / 0.27e2 * t393 * t285 * t237 + 0.8e1 / 0.9e1 * t236 * t474 * t138 + 0.4e1 / 0.9e1 * t469 * t244 + t713 + 0.40e2 / 0.81e2 * t656 * t288 * t249 + 0.64e2 / 0.27e2 * t442 * t140 * t70 * t241 - 0.8e1 / 0.27e2 * t477 * t252 + 0.32e2 / 0.9e1 * t248 * t770 * t773 - 0.16e2 / 0.9e1 * t268 * t241 + 0.16e2 / 0.3e1 * t268 * t400 - 0.8e1 / 0.27e2 * t406 * t291 * t249 + 0.8e1 / 0.9e1 * t248 * t485 * t140 + 0.4e1 / 0.9e1 * t482 * t252 - t731;
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[2] = t519 + t520 + 0.6e1 * t387 + 0.2e1 * t390 + 0.6e1 * t425 + 0.8e1 * t427 + 0.2e1 * t457 + 0.2e1 * t463 + 0.2e1 * t490 + t7 * (t69 * t812 * t84 + t582 + t625 + 0.2e1 * t627 + t630 + 0.2e1 * t683 + 0.4e1 * t685 + 0.2e1 * t688 + t757 + 0.2e1 * t759);

  t824 = t132 * t514 * t84;
  t843 = 0.12e2 * t321 + 0.24e2 * t649;
  t865 = 0.40e2 / 0.81e2 * t637 * t497 * t138 - 0.16e2 / 0.9e1 * t464 * t242 - 0.8e1 / 0.9e1 * t429 * t285 * t138 + 0.8e1 / 0.3e1 * t432 * t241 * t285 + 0.4e1 / 0.3e1 * t262 * t474 + 0.4e1 / 0.9e1 * t236 * t503 * t138 + 0.4e1 / 0.3e1 * t74 * t843 + 0.40e2 / 0.81e2 * t656 * t506 * t140 + 0.16e2 / 0.9e1 * t477 * t242 - 0.8e1 / 0.9e1 * t442 * t291 * t140 - 0.8e1 / 0.3e1 * t445 * t241 * t291 + 0.4e1 / 0.3e1 * t268 * t485 + 0.4e1 / 0.9e1 * t248 * t511 * t140 - 0.4e1 / 0.3e1 * t77 * t843;
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[3] = t519 + t520 + t388 + 0.9e1 * t425 + t678 + 0.6e1 * t463 + 0.3e1 * t490 + t516 + t7 * (t69 * t865 * t84 + t582 + t625 + t627 + 0.3e1 * t683 + t686 + 0.3e1 * t757 + 0.3e1 * t759 + t824);

  t876 = t281 * t281;
  t881 = t285 * t285;
  t887 = 0.24e2 * t321 + 0.24e2 * t649;
  t890 = t288 * t288;
  t895 = t291 * t291;
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[4] = t519 + t520 + 0.12e2 * t425 + 0.12e2 * t463 + 0.4e1 * t516 + t7 * (t582 + t625 + 0.4e1 * t683 + 0.6e1 * t757 + 0.4e1 * t824 + t69 * (0.40e2 / 0.81e2 * t637 * t876 - 0.16e2 / 0.9e1 * t464 * t285 + 0.4e1 / 0.3e1 * t236 * t881 + 0.16e2 / 0.9e1 * t262 * t503 + 0.4e1 / 0.3e1 * t74 * t887 + 0.40e2 / 0.81e2 * t656 * t890 - 0.16e2 / 0.9e1 * t477 * t291 + 0.4e1 / 0.3e1 * t248 * t895 + 0.16e2 / 0.9e1 * t268 * t511 - 0.4e1 / 0.3e1 * t77 * t887) * t84);

#ifndef XC_DONT_COMPILE_MXC

  if(order < 5) return;


#endif

#endif

#endif

#endif

#endif


}

