/* 
  This file was generated automatically with ./scripts/maple2c.pl.
  Do not edit this file directly as it can be overwritten!!

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.

  Maple version     : Maple 2016 (X86 64 LINUX)
  Maple source      : ./maple/lda_exc/lda_x_rel.mpl
  Type of functional: lda_exc
*/

#define maple2c_order 4
#define MAPLE2C_FLAGS (XC_FLAGS_I_HAVE_EXC | XC_FLAGS_I_HAVE_VXC | XC_FLAGS_I_HAVE_FXC | XC_FLAGS_I_HAVE_KXC | XC_FLAGS_I_HAVE_LXC)


static inline void
func_unpol(const xc_func_type *p, int order, const double *rho, double *zk, LDA_OUT_PARAMS_NO_EXC(double *))
{

#ifndef XC_DONT_COMPILE_EXC
  double t1, t2, t3, t4, t5, t6, t7, t8;
  double t9, t10, t11, t12, t13, t14, t15, t16;
  double t17, t21, t22, t23, t28, t29, t30, t34;
  double t35, t36, t37, t41, t42, t44, t46;

#ifndef XC_DONT_COMPILE_VXC
  double t49, t51, t52, t53, t54, t55, t59, t60;
  double t63, t67, t71, t72;

#ifndef XC_DONT_COMPILE_FXC
  double t76, t83, t88, t92, t95, t97, t98, t101;
  double t107, t111, t112;

#ifndef XC_DONT_COMPILE_KXC
  double t116, t129, t133, t135, t137, t140, t143, t146;
  double t149, t159, t163, t164;

#ifndef XC_DONT_COMPILE_LXC
  double t187, t197, t207, t210, t213, t216;
#endif

#endif

#endif

#endif

#endif



  t1 = M_CBRT3;
  t2 = 0.1e1 / M_PI;
  t3 = POW_1_3(t2);
  t4 = t1 * t3;
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = t4 * t6;
  t8 = M_CBRT2;
  t9 = t8 * t8;
  t10 = POW_1_3(rho[0]);
  t11 = t9 * t10;
  t12 = POW_1_3(0.9e1);
  t13 = t12 * t12;
  t14 = t13 * t1;
  t15 = t3 * t3;
  t16 = 0.1e1 / t15;
  t17 = t10 * t10;
  t21 = 0.1e1 + 0.38075239991386496937e-4 * t14 * t16 * t17;
  t22 = sqrt(t21);
  t23 = t22 * t13;
  t28 = t1 * t1;
  t29 = t12 * t28;
  t30 = 0.1e1 / t3;
  t34 = log(0.35625477770544353752e-2 * t29 * t30 * t10 + sqrt(POW_2(0.35625477770544353752e-2 * t29 * t30 * t10) + 0.1e1));
  t35 = t34 * t12;
  t36 = t28 * t15;
  t37 = 0.1e1 / t17;
  t41 = 0.10396221848752237744e2 * t23 * t4 / t10 - 0.97273285855626056446e3 * t35 * t36 * t37;
  t42 = t41 * t41;
  t44 = 0.1e1 - 0.15e1 * t42;
  t46 = t7 * t11 * t44;
  if(zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    zk[0] = -0.3e1 / 0.16e2 * t46;

#ifndef XC_DONT_COMPILE_VXC

  if(order < 1) return;


  t49 = t10 * rho[0];
  t51 = t49 * t1 * t3;
  t52 = t6 * t9;
  t53 = 0.1e1 / t22;
  t54 = t53 * t12;
  t55 = t28 * t30;
  t59 = 0.1e1 / t49;
  t60 = t4 * t59;
  t63 = t53 * t13;
  t67 = 0.1e1 / t17 / rho[0];
  t71 = 0.11875159256848117917e-2 * t54 * t55 * t37 - 0.34654072829174125813e1 * t23 * t60 - 0.34654072829174125814e1 * t63 * t60 + 0.64848857237084037631e3 * t35 * t36 * t67;
  t72 = t41 * t71;
  if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    vrho[0] = -t46 / 0.4e1 + 0.56250000000000000000e0 * t51 * t52 * t72;

#ifndef XC_DONT_COMPILE_FXC

  if(order < 2) return;


  t76 = t9 * t37;
  t83 = t71 * t71;
  t88 = 0.1e1 / t22 / t21;
  t92 = t55 * t67;
  t95 = rho[0] * rho[0];
  t97 = 0.1e1 / t10 / t95;
  t98 = t4 * t97;
  t101 = t88 * t12;
  t107 = 0.1e1 / t17 / t95;
  t111 = -0.12784227020251018670e-5 * t88 / rho[0] - 0.11875159256848117917e-2 * t54 * t92 + 0.46205430438898834417e1 * t23 * t98 + 0.39583864189493726391e-3 * t101 * t92 + 0.69308145658348251628e1 * t63 * t98 - 0.10808142872847339605e4 * t35 * t36 * t107;
  t112 = t41 * t111;
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[0] = -t7 * t76 * t44 / 0.12e2 + 0.15000000000000000000e1 * t7 * t11 * t72 + 0.56250000000000000000e0 * t51 * t52 * t83 + 0.56250000000000000000e0 * t51 * t52 * t112;

#ifndef XC_DONT_COMPILE_KXC

  if(order < 3) return;


  t116 = t9 * t67;
  t129 = t71 * t111;
  t133 = t21 * t21;
  t135 = 0.1e1 / t22 / t133;
  t137 = t14 * t16;
  t140 = 0.1e1 / t95;
  t143 = t55 * t107;
  t146 = t95 * rho[0];
  t149 = t4 / t10 / t146;
  t159 = 0.1e1 / t17 / t146;
  t163 = 0.48676251190042541751e-10 * t135 * t59 * t137 + 0.25568454040502037340e-5 * t88 * t140 + 0.25069780653346026714e-2 * t54 * t143 - 0.10781267102409728031e2 * t23 * t149 - 0.12784227020251018671e-5 * t135 * t140 - 0.14514083536147699677e-2 * t101 * t143 - 0.20022353190189494915e2 * t63 * t149 + 0.28821714327592905613e4 * t35 * t36 * t159;
  t164 = t41 * t163;
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[0] = t7 * t116 * t44 / 0.18e2 + 0.75000000000000000000e0 * t7 * t76 * t72 + 0.22500000000000000000e1 * t7 * t11 * t83 + 0.22500000000000000000e1 * t7 * t11 * t112 + 0.16875000000000000000e1 * t51 * t52 * t129 + 0.56250000000000000000e0 * t51 * t52 * t164;

#ifndef XC_DONT_COMPILE_LXC

  if(order < 4) return;


  t187 = t111 * t111;
  t197 = 0.1e1 / t22 / t133 / t21;
  t207 = 0.1e1 / t146;
  t210 = t55 * t159;
  t213 = t95 * t95;
  t216 = t4 / t10 / t213;
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[0] = -0.5e1 / 0.54e2 * t7 * t9 * t107 * t44 - 0.66666666666666666667e0 * t7 * t116 * t72 + 0.15000000000000000000e1 * t7 * t76 * t83 + 0.15000000000000000000e1 * t7 * t76 * t112 + 0.90000000000000000000e1 * t7 * t11 * t129 + 0.30000000000000000000e1 * t7 * t11 * t164 + 0.16875000000000000000e1 * t51 * t52 * t187 + 0.22500000000000000000e1 * t51 * t52 * t71 * t163 + 0.56250000000000000000e0 * t51 * t52 * t41 * (-0.27800399189128235230e-13 * t197 * t67 * t29 / t3 / t2 - 0.16225417063347513917e-9 * t135 * t97 * t137 - 0.78125831790422891872e-5 * t88 * t207 - 0.79167728378987452781e-2 * t54 * t210 + 0.35937557008032426770e2 * t23 * t216 + 0.81127085316737569591e-10 * t197 * t97 * t137 + 0.72443953114755772469e-5 * t135 * t207 + 0.61574899850323574387e-2 * t101 * t210 + 0.77009050731498057365e2 * t63 * t216 - 0.10567961920117398725e5 * t35 * t36 / t17 / t213);

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
  double t1, t2, t3, t4, t5, t6, t7, t8;
  double t9, t10, t11, t12, t13, t14, t18, t19;
  double t20, t25, t26, t27, t31, t32, t33, t34;
  double t38, t39, t41, t43;

#ifndef XC_DONT_COMPILE_VXC
  double t46, t47, t48, t49, t50, t51, t52, t56;
  double t57, t60, t64, t68;

#ifndef XC_DONT_COMPILE_FXC
  double t72, t76, t77, t81, t82, t87, t91, t94;
  double t96, t97, t100, t106, t110;

#ifndef XC_DONT_COMPILE_KXC
  double t118, t128, t132, t134, t136, t139, t142, t145;
  double t148, t158, t162;

#ifndef XC_DONT_COMPILE_LXC
  double t187, t196, t206, t209, t212, t215;
#endif

#endif

#endif

#endif

#endif



  t1 = M_CBRT3;
  t2 = 0.1e1 / M_PI;
  t3 = POW_1_3(t2);
  t4 = t1 * t3;
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = POW_1_3(rho[0]);
  t8 = t6 * t7;
  t9 = POW_1_3(0.9e1);
  t10 = t9 * t9;
  t11 = t10 * t1;
  t12 = t3 * t3;
  t13 = 0.1e1 / t12;
  t14 = t7 * t7;
  t18 = 0.1e1 + 0.38075239991386496937e-4 * t11 * t13 * t14;
  t19 = sqrt(t18);
  t20 = t19 * t10;
  t25 = t1 * t1;
  t26 = t9 * t25;
  t27 = 0.1e1 / t3;
  t31 = log(0.35625477770544353752e-2 * t26 * t27 * t7 + sqrt(POW_2(0.35625477770544353752e-2 * t26 * t27 * t7) + 0.1e1));
  t32 = t31 * t9;
  t33 = t25 * t12;
  t34 = 0.1e1 / t14;
  t38 = 0.10396221848752237744e2 * t20 * t4 / t7 - 0.97273285855626056446e3 * t32 * t33 * t34;
  t39 = t38 * t38;
  t41 = 0.1e1 - 0.15e1 * t39;
  t43 = t4 * t8 * t41;
  if(zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    zk[0] = -0.3e1 / 0.8e1 * t43;

#ifndef XC_DONT_COMPILE_VXC

  if(order < 1) return;


  t46 = t7 * rho[0];
  t47 = t46 * t1;
  t48 = t47 * t3;
  t49 = t6 * t38;
  t50 = 0.1e1 / t19;
  t51 = t50 * t9;
  t52 = t25 * t27;
  t56 = 0.1e1 / t46;
  t57 = t4 * t56;
  t60 = t50 * t10;
  t64 = 0.1e1 / t14 / rho[0];
  t68 = 0.11875159256848117917e-2 * t51 * t52 * t34 - 0.34654072829174125813e1 * t20 * t57 - 0.34654072829174125814e1 * t60 * t57 + 0.64848857237084037631e3 * t32 * t33 * t64;
  if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    vrho[0] = -t43 / 0.2e1 + 0.11250000000000000000e1 * t48 * t49 * t68;

#ifndef XC_DONT_COMPILE_FXC

  if(order < 2) return;


  t72 = t6 * t34;
  t76 = t4 * t6;
  t77 = t7 * t38;
  t81 = t3 * t6;
  t82 = t68 * t68;
  t87 = 0.1e1 / t19 / t18;
  t91 = t52 * t64;
  t94 = rho[0] * rho[0];
  t96 = 0.1e1 / t7 / t94;
  t97 = t4 * t96;
  t100 = t87 * t9;
  t106 = 0.1e1 / t14 / t94;
  t110 = -0.12784227020251018670e-5 * t87 / rho[0] - 0.11875159256848117917e-2 * t51 * t91 + 0.46205430438898834417e1 * t20 * t97 + 0.39583864189493726391e-3 * t100 * t91 + 0.69308145658348251628e1 * t60 * t97 - 0.10808142872847339605e4 * t32 * t33 * t106;
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[0] = -t4 * t72 * t41 / 0.6e1 + 0.30000000000000000000e1 * t76 * t77 * t68 + 0.11250000000000000000e1 * t47 * t81 * t82 + 0.11250000000000000000e1 * t48 * t49 * t110;

#ifndef XC_DONT_COMPILE_KXC

  if(order < 3) return;


  t118 = t34 * t38;
  t128 = t6 * t68;
  t132 = t18 * t18;
  t134 = 0.1e1 / t19 / t132;
  t136 = t11 * t13;
  t139 = 0.1e1 / t94;
  t142 = t52 * t106;
  t145 = t94 * rho[0];
  t148 = t4 / t7 / t145;
  t158 = 0.1e1 / t14 / t145;
  t162 = 0.48676251190042541751e-10 * t134 * t56 * t136 + 0.25568454040502037340e-5 * t87 * t139 + 0.25069780653346026714e-2 * t51 * t142 - 0.10781267102409728031e2 * t20 * t148 - 0.12784227020251018671e-5 * t134 * t139 - 0.14514083536147699677e-2 * t100 * t142 - 0.20022353190189494915e2 * t60 * t148 + 0.28821714327592905613e4 * t32 * t33 * t158;
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[0] = t4 * t6 * t64 * t41 / 0.9e1 + 0.15000000000000000000e1 * t76 * t118 * t68 + 0.45000000000000000000e1 * t4 * t8 * t82 + 0.45000000000000000000e1 * t76 * t77 * t110 + 0.33750000000000000000e1 * t48 * t128 * t110 + 0.11250000000000000000e1 * t48 * t49 * t162;

#ifndef XC_DONT_COMPILE_LXC

  if(order < 4) return;


  t187 = t110 * t110;
  t196 = 0.1e1 / t19 / t132 / t18;
  t206 = 0.1e1 / t145;
  t209 = t52 * t158;
  t212 = t94 * t94;
  t215 = t4 / t7 / t212;
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[0] = -0.5e1 / 0.27e2 * t4 * t6 * t106 * t41 - 0.13333333333333333333e1 * t76 * t64 * t38 * t68 + 0.30000000000000000000e1 * t4 * t72 * t82 + 0.30000000000000000000e1 * t76 * t118 * t110 + 0.18000000000000000000e2 * t76 * t7 * t68 * t110 + 0.60000000000000000000e1 * t76 * t77 * t162 + 0.33750000000000000000e1 * t47 * t81 * t187 + 0.45000000000000000000e1 * t48 * t128 * t162 + 0.11250000000000000000e1 * t48 * t49 * (-0.27800399189128235230e-13 * t196 * t64 * t26 / t3 / t2 - 0.16225417063347513917e-9 * t134 * t96 * t136 - 0.78125831790422891872e-5 * t87 * t206 - 0.79167728378987452781e-2 * t51 * t209 + 0.35937557008032426770e2 * t20 * t215 + 0.81127085316737569591e-10 * t196 * t96 * t136 + 0.72443953114755772469e-5 * t134 * t206 + 0.61574899850323574387e-2 * t100 * t209 + 0.77009050731498057365e2 * t60 * t215 - 0.10567961920117398725e5 * t32 * t33 / t14 / t212);

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
  double t1, t2, t3, t4, t5, t6, t7, t8;
  double t9, t10, t11, t12, t13, t14, t15, t17;
  double t18, t20, t21, t22, t23, t24, t25, t26;
  double t27, t28, t32, t33, t34, t39, t40, t41;
  double t45, t46, t47, t48, t52, t53, t55, t56;
  double t58;

#ifndef XC_DONT_COMPILE_VXC
  double t60, t61, t62, t63, t64, t65, t66, t67;
  double t68, t70, t73, t78, t79, t80, t81, t82;
  double t86, t87, t90, t94, t98, t99, t102, t103;
  double t105, t108, t110;

#ifndef XC_DONT_COMPILE_FXC
  double t113, t115, t117, t120, t121, t122, t125, t126;
  double t127, t128, t131, t132, t133, t135, t138, t139;
  double t140, t143, t146, t152, t154, t158, t160, t163;
  double t167, t168, t171, t177, t181, t182, t185, t187;
  double t188, t189, t191, t194, t197, t200, t203, t205;
  double t208, t209, t210, t214, t218, t221, t224, t227;
  double t229;

#ifndef XC_DONT_COMPILE_KXC
  double t233, t235, t238, t240, t242, t244, t247, t248;
  double t251, t252, t255, t258, t260, t261, t264, t267;
  double t268, t269, t271, t275, t276, t279, t282, t285;
  double t291, t295, t298, t300, t303, t304, t306, t308;
  double t313, t318, t328, t332, t333, t336, t339, t340;
  double t343, t353, t356, t366, t368, t376, t377, t378;
  double t380, t381, t383, t384, t386, t388, t389, t391;
  double t392, t398, t400, t405, t410, t413, t418, t421;
  double t424, t426, t429, t430, t431, t440, t446, t449;
  double t454, t457, t459;

#ifndef XC_DONT_COMPILE_LXC
  double t466, t470, t473, t477, t480, t483, t488, t492;
  double t496, t499, t502, t503, t506, t507, t509, t513;
  double t516, t519, t522, t524, t528, t532, t535, t547;
  double t552, t573, t575, t577, t579, t580, t586, t592;
  double t593, t595, t598, t600, t601, t607, t620, t625;
  double t626, t629, t632, t634, t637, t638, t640, t643;
  double t646, t654, t655, t658, t660, t664, t665, t666;
  double t685, t703, t718, t723, t724, t728, t731, t744;
  double t756, t760, t762, t773, t776, t815, t820, t823;
  double t834, t862, t884, t889, t891, t894, t899, t905;
  double t908, t913, t934, t937;
#endif

#endif

#endif

#endif

#endif



  t1 = M_CBRT3;
  t2 = 0.1e1 / M_PI;
  t3 = POW_1_3(t2);
  t4 = t1 * t3;
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = t4 * t6;
  t8 = M_CBRT2;
  t9 = t8 * t8;
  t10 = rho[0] - rho[1];
  t11 = rho[0] + rho[1];
  t12 = 0.1e1 / t11;
  t13 = t10 * t12;
  t14 = 0.1e1 + t13;
  t15 = POW_1_3(t14);
  t17 = 0.1e1 - t13;
  t18 = POW_1_3(t17);
  t20 = t15 * t14 + t18 * t17;
  t21 = t9 * t20;
  t22 = POW_1_3(t11);
  t23 = POW_1_3(0.9e1);
  t24 = t23 * t23;
  t25 = t24 * t1;
  t26 = t3 * t3;
  t27 = 0.1e1 / t26;
  t28 = t22 * t22;
  t32 = 0.1e1 + 0.38075239991386496937e-4 * t25 * t27 * t28;
  t33 = sqrt(t32);
  t34 = t33 * t24;
  t39 = t1 * t1;
  t40 = t23 * t39;
  t41 = 0.1e1 / t3;
  t45 = log(0.35625477770544353752e-2 * t40 * t41 * t22 + sqrt(POW_2(0.35625477770544353752e-2 * t40 * t41 * t22) + 0.1e1));
  t46 = t45 * t23;
  t47 = t39 * t26;
  t48 = 0.1e1 / t28;
  t52 = 0.10396221848752237744e2 * t34 * t4 / t22 - 0.97273285855626056446e3 * t46 * t47 * t48;
  t53 = t52 * t52;
  t55 = 0.1e1 - 0.15e1 * t53;
  t56 = t22 * t55;
  t58 = t7 * t21 * t56;
  if(zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    zk[0] = -0.3e1 / 0.32e2 * t58;

#ifndef XC_DONT_COMPILE_VXC

  if(order < 1) return;


  t60 = t58 / 0.8e1;
  t61 = t22 * t11;
  t62 = t61 * t1;
  t63 = t62 * t3;
  t64 = t6 * t9;
  t65 = t11 * t11;
  t66 = 0.1e1 / t65;
  t67 = t10 * t66;
  t68 = t12 - t67;
  t70 = -t68;
  t73 = 0.4e1 / 0.3e1 * t15 * t68 + 0.4e1 / 0.3e1 * t18 * t70;
  t78 = t3 * t6;
  t79 = t62 * t78;
  t80 = 0.1e1 / t33;
  t81 = t80 * t23;
  t82 = t39 * t41;
  t86 = 0.1e1 / t61;
  t87 = t4 * t86;
  t90 = t80 * t24;
  t94 = 0.1e1 / t28 / t11;
  t98 = 0.11875159256848117917e-2 * t81 * t82 * t48 - 0.34654072829174125813e1 * t34 * t87 - 0.34654072829174125814e1 * t90 * t87 + 0.64848857237084037631e3 * t46 * t47 * t94;
  t99 = t52 * t98;
  t102 = 0.28125000000000000000e0 * t79 * t21 * t99;
  if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    vrho[0] = -t60 - 0.3e1 / 0.32e2 * t63 * t64 * t73 * t55 + t102;

  t103 = -t12 - t67;
  t105 = -t103;
  t108 = 0.4e1 / 0.3e1 * t15 * t103 + 0.4e1 / 0.3e1 * t18 * t105;
  t110 = t64 * t108 * t55;
  if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    vrho[1] = -t60 - 0.3e1 / 0.32e2 * t63 * t110 + t102;

#ifndef XC_DONT_COMPILE_FXC

  if(order < 2) return;


  t113 = t9 * t73;
  t115 = t7 * t113 * t56;
  t117 = t48 * t55;
  t120 = t7 * t21 * t117 / 0.24e2;
  t121 = t4 * t64;
  t122 = t20 * t22;
  t125 = 0.75000000000000000000e0 * t121 * t122 * t99;
  t126 = t15 * t15;
  t127 = 0.1e1 / t126;
  t128 = t68 * t68;
  t131 = t65 * t11;
  t132 = 0.1e1 / t131;
  t133 = t10 * t132;
  t135 = -0.2e1 * t66 + 0.2e1 * t133;
  t138 = t18 * t18;
  t139 = 0.1e1 / t138;
  t140 = t70 * t70;
  t143 = -t135;
  t146 = 0.4e1 / 0.9e1 * t127 * t128 + 0.4e1 / 0.3e1 * t15 * t135 + 0.4e1 / 0.9e1 * t139 * t140 + 0.4e1 / 0.3e1 * t18 * t143;
  t152 = t79 * t113 * t99;
  t154 = t98 * t98;
  t158 = 0.28125000000000000000e0 * t63 * t64 * t20 * t154;
  t160 = 0.1e1 / t33 / t32;
  t163 = t82 * t94;
  t167 = 0.1e1 / t22 / t65;
  t168 = t4 * t167;
  t171 = t160 * t23;
  t177 = 0.1e1 / t28 / t65;
  t181 = -0.12784227020251018670e-5 * t160 * t12 - 0.11875159256848117917e-2 * t81 * t163 + 0.46205430438898834417e1 * t34 * t168 + 0.39583864189493726391e-3 * t171 * t163 + 0.69308145658348251628e1 * t90 * t168 - 0.10808142872847339605e4 * t46 * t47 * t177;
  t182 = t52 * t181;
  t185 = 0.28125000000000000000e0 * t79 * t21 * t182;
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[0] = -t115 / 0.4e1 - t120 + t125 - 0.3e1 / 0.32e2 * t63 * t64 * t146 * t55 + 0.56250000000000000000e0 * t152 + t158 + t185;

  t187 = t22 * t1;
  t188 = t187 * t3;
  t189 = t188 * t110;
  t191 = t127 * t103;
  t194 = t15 * t10;
  t197 = t139 * t105;
  t200 = t18 * t10;
  t203 = 0.4e1 / 0.9e1 * t191 * t68 + 0.8e1 / 0.3e1 * t194 * t132 + 0.4e1 / 0.9e1 * t197 * t70 - 0.8e1 / 0.3e1 * t200 * t132;
  t205 = t64 * t203 * t55;
  t208 = t9 * t108;
  t209 = t208 * t99;
  t210 = t79 * t209;
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[1] = -t115 / 0.8e1 - t120 + t125 - t189 / 0.8e1 - 0.3e1 / 0.32e2 * t63 * t205 + 0.28125000000000000000e0 * t210 + 0.28125000000000000000e0 * t152 + t158 + t185;

  t214 = t103 * t103;
  t218 = 0.2e1 * t66 + 0.2e1 * t133;
  t221 = t105 * t105;
  t224 = -t218;
  t227 = 0.4e1 / 0.9e1 * t127 * t214 + 0.4e1 / 0.3e1 * t15 * t218 + 0.4e1 / 0.9e1 * t139 * t221 + 0.4e1 / 0.3e1 * t18 * t224;
  t229 = t64 * t227 * t55;
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[2] = -t189 / 0.4e1 - t120 + t125 - 0.3e1 / 0.32e2 * t63 * t229 + 0.56250000000000000000e0 * t210 + t158 + t185;

#ifndef XC_DONT_COMPILE_KXC

  if(order < 3) return;


  t233 = t9 * t146;
  t235 = t7 * t233 * t56;
  t238 = t7 * t113 * t117;
  t240 = t73 * t22;
  t242 = t121 * t240 * t99;
  t244 = t94 * t55;
  t247 = t7 * t21 * t244 / 0.36e2;
  t248 = t20 * t48;
  t251 = 0.37500000000000000000e0 * t121 * t248 * t99;
  t252 = t22 * t154;
  t255 = 0.11250000000000000000e1 * t7 * t21 * t252;
  t258 = 0.11250000000000000000e1 * t121 * t122 * t182;
  t260 = 0.1e1 / t126 / t14;
  t261 = t128 * t68;
  t264 = t127 * t68;
  t267 = t65 * t65;
  t268 = 0.1e1 / t267;
  t269 = t10 * t268;
  t271 = 0.6e1 * t132 - 0.6e1 * t269;
  t275 = 0.1e1 / t138 / t17;
  t276 = t140 * t70;
  t279 = t139 * t70;
  t282 = -t271;
  t285 = -0.8e1 / 0.27e2 * t260 * t261 + 0.4e1 / 0.3e1 * t264 * t135 + 0.4e1 / 0.3e1 * t15 * t271 - 0.8e1 / 0.27e2 * t275 * t276 + 0.4e1 / 0.3e1 * t279 * t143 + 0.4e1 / 0.3e1 * t18 * t282;
  t291 = t79 * t233 * t99;
  t295 = t63 * t64 * t73 * t154;
  t298 = t79 * t113 * t182;
  t300 = t98 * t181;
  t303 = 0.84375000000000000000e0 * t79 * t21 * t300;
  t304 = t32 * t32;
  t306 = 0.1e1 / t33 / t304;
  t308 = t25 * t27;
  t313 = t82 * t177;
  t318 = t4 / t22 / t131;
  t328 = 0.1e1 / t28 / t131;
  t332 = 0.48676251190042541751e-10 * t306 * t86 * t308 + 0.25568454040502037340e-5 * t160 * t66 + 0.25069780653346026714e-2 * t81 * t313 - 0.10781267102409728031e2 * t34 * t318 - 0.12784227020251018671e-5 * t306 * t66 - 0.14514083536147699677e-2 * t171 * t313 - 0.20022353190189494915e2 * t90 * t318 + 0.28821714327592905613e4 * t46 * t47 * t328;
  t333 = t52 * t332;
  t336 = 0.28125000000000000000e0 * t79 * t21 * t333;
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[0] = -0.3e1 / 0.8e1 * t235 - t238 / 0.8e1 + 0.22500000000000000000e1 * t242 + t247 + t251 + t255 + t258 - 0.3e1 / 0.32e2 * t63 * t64 * t285 * t55 + 0.84375000000000000000e0 * t291 + 0.84375000000000000000e0 * t295 + 0.84375000000000000000e0 * t298 + t303 + t336;

  t339 = t188 * t205 / 0.4e1;
  t340 = t260 * t103;
  t343 = t127 * t10;
  t353 = t275 * t105;
  t356 = t139 * t10;
  t366 = -0.8e1 / 0.27e2 * t340 * t128 + 0.16e2 / 0.9e1 * t343 * t132 * t68 + 0.4e1 / 0.9e1 * t191 * t135 + 0.8e1 / 0.3e1 * t15 * t132 - 0.8e1 * t194 * t268 - 0.8e1 / 0.27e2 * t353 * t140 - 0.16e2 / 0.9e1 * t356 * t132 * t70 + 0.4e1 / 0.9e1 * t197 * t143 - 0.8e1 / 0.3e1 * t18 * t132 + 0.8e1 * t200 * t268;
  t368 = t64 * t366 * t55;
  t376 = t48 * t1;
  t377 = t376 * t3;
  t378 = t377 * t110;
  t380 = t187 * t78;
  t381 = t380 * t209;
  t383 = t9 * t203;
  t384 = t383 * t99;
  t386 = 0.56250000000000000000e0 * t79 * t384;
  t388 = t64 * t108 * t154;
  t389 = t63 * t388;
  t391 = t208 * t182;
  t392 = t79 * t391;
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[1] = -t235 / 0.8e1 - t339 - 0.3e1 / 0.32e2 * t63 * t368 - t238 / 0.12e2 + 0.15000000000000000000e1 * t242 + t247 + t251 + t255 + t258 + 0.28125000000000000000e0 * t291 + 0.56250000000000000000e0 * t295 + 0.56250000000000000000e0 * t298 + t303 + t336 - t378 / 0.24e2 + 0.75000000000000000000e0 * t381 + t386 + 0.28125000000000000000e0 * t389 + 0.28125000000000000000e0 * t392;

  t398 = t188 * t229;
  t400 = t260 * t214;
  t405 = t127 * t218;
  t410 = -0.2e1 * t132 - 0.6e1 * t269;
  t413 = t275 * t221;
  t418 = t139 * t224;
  t421 = -t410;
  t424 = -0.8e1 / 0.27e2 * t400 * t68 + 0.16e2 / 0.9e1 * t191 * t133 + 0.4e1 / 0.9e1 * t405 * t68 + 0.4e1 / 0.3e1 * t15 * t410 - 0.8e1 / 0.27e2 * t413 * t70 - 0.16e2 / 0.9e1 * t197 * t133 + 0.4e1 / 0.9e1 * t418 * t70 + 0.4e1 / 0.3e1 * t18 * t421;
  t426 = t64 * t424 * t55;
  t429 = t9 * t227;
  t430 = t429 * t99;
  t431 = t79 * t430;
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[2] = -t378 / 0.12e2 - t339 + 0.15000000000000000000e1 * t381 - t238 / 0.24e2 + t247 + t251 + 0.75000000000000000000e0 * t242 + t255 + t258 - t398 / 0.8e1 - 0.3e1 / 0.32e2 * t63 * t426 + 0.28125000000000000000e0 * t431 + t386 + 0.56250000000000000000e0 * t389 + 0.56250000000000000000e0 * t392 + 0.28125000000000000000e0 * t295 + t303 + 0.28125000000000000000e0 * t298 + t336;

  t440 = t214 * t103;
  t446 = -0.6e1 * t132 - 0.6e1 * t269;
  t449 = t221 * t105;
  t454 = -t446;
  t457 = -0.8e1 / 0.27e2 * t260 * t440 + 0.4e1 / 0.3e1 * t191 * t218 + 0.4e1 / 0.3e1 * t15 * t446 - 0.8e1 / 0.27e2 * t275 * t449 + 0.4e1 / 0.3e1 * t197 * t224 + 0.4e1 / 0.3e1 * t18 * t454;
  t459 = t64 * t457 * t55;
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[3] = -t378 / 0.8e1 - 0.3e1 / 0.8e1 * t398 + 0.22500000000000000000e1 * t381 + t247 + t251 + t255 + t258 - 0.3e1 / 0.32e2 * t63 * t459 + 0.84375000000000000000e0 * t431 + 0.84375000000000000000e0 * t389 + 0.84375000000000000000e0 * t392 + t303 + t336;

#ifndef XC_DONT_COMPILE_LXC

  if(order < 4) return;


  t466 = t7 * t233 * t117;
  t470 = t121 * t146 * t22 * t99;
  t473 = t7 * t113 * t244;
  t477 = t121 * t73 * t48 * t99;
  t480 = t7 * t113 * t252;
  t483 = t121 * t240 * t182;
  t488 = 0.5e1 / 0.108e3 * t7 * t21 * t177 * t55;
  t492 = 0.33333333333333333333e0 * t121 * t20 * t94 * t99;
  t496 = 0.75000000000000000000e0 * t7 * t21 * t48 * t154;
  t499 = 0.75000000000000000000e0 * t121 * t248 * t182;
  t502 = 0.45000000000000000000e1 * t121 * t122 * t300;
  t503 = -t466 / 0.4e1 + 0.45000000000000000000e1 * t470 + t473 / 0.9e1 + 0.15000000000000000000e1 * t477 + 0.45000000000000000000e1 * t480 + 0.45000000000000000000e1 * t483 - t488 - t492 + t496 + t499 + t502;
  t506 = 0.15000000000000000000e1 * t121 * t122 * t333;
  t507 = t9 * t285;
  t509 = t79 * t507 * t99;
  t513 = t63 * t64 * t146 * t154;
  t516 = t79 * t233 * t182;
  t519 = t79 * t113 * t300;
  t522 = t79 * t113 * t333;
  t524 = t181 * t181;
  t528 = 0.84375000000000000000e0 * t63 * t64 * t20 * t524;
  t532 = 0.11250000000000000000e1 * t79 * t21 * t98 * t332;
  t535 = 0.1e1 / t33 / t304 / t32;
  t547 = t82 * t328;
  t552 = t4 / t22 / t267;
  t573 = 0.28125000000000000000e0 * t79 * t21 * t52 * (-0.27800399189128235230e-13 * t535 * t94 * t40 / t3 / t2 - 0.16225417063347513917e-9 * t306 * t167 * t308 - 0.78125831790422891872e-5 * t160 * t132 - 0.79167728378987452781e-2 * t81 * t547 + 0.35937557008032426770e2 * t34 * t552 + 0.81127085316737569591e-10 * t535 * t167 * t308 + 0.72443953114755772469e-5 * t306 * t132 + 0.61574899850323574387e-2 * t171 * t547 + 0.77009050731498057365e2 * t90 * t552 - 0.10567961920117398725e5 * t46 * t47 / t28 / t267);
  t575 = t7 * t507 * t56;
  t577 = t14 * t14;
  t579 = 0.1e1 / t126 / t577;
  t580 = t128 * t128;
  t586 = t135 * t135;
  t592 = 0.1e1 / t267 / t11;
  t593 = t10 * t592;
  t595 = -0.24e2 * t268 + 0.24e2 * t593;
  t598 = t17 * t17;
  t600 = 0.1e1 / t138 / t598;
  t601 = t140 * t140;
  t607 = t143 * t143;
  t620 = t506 + 0.11250000000000000000e1 * t509 + 0.16875000000000000000e1 * t513 + 0.16875000000000000000e1 * t516 + 0.33750000000000000000e1 * t519 + 0.11250000000000000000e1 * t522 + t528 + t532 + t573 - t575 / 0.2e1 - 0.3e1 / 0.32e2 * t63 * t64 * (0.40e2 / 0.81e2 * t579 * t580 - 0.16e2 / 0.9e1 * t260 * t128 * t135 + 0.4e1 / 0.3e1 * t127 * t586 + 0.16e2 / 0.9e1 * t264 * t271 + 0.4e1 / 0.3e1 * t15 * t595 + 0.40e2 / 0.81e2 * t600 * t601 - 0.16e2 / 0.9e1 * t275 * t140 * t143 + 0.4e1 / 0.3e1 * t139 * t607 + 0.16e2 / 0.9e1 * t279 * t282 - 0.4e1 / 0.3e1 * t18 * t595) * t55;
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[0] = t503 + t620;

  t625 = t380 * t384;
  t626 = 0.22500000000000000000e1 * t625;
  t629 = t79 * t9 * t366 * t99;
  t632 = t376 * t78 * t209;
  t634 = t380 * t391;
  t637 = t79 * t383 * t182;
  t638 = 0.84375000000000000000e0 * t637;
  t640 = t79 * t208 * t300;
  t643 = t79 * t208 * t333;
  t646 = t499 + t502 + t506 + 0.28125000000000000000e0 * t509 + 0.84375000000000000000e0 * t516 + 0.25312500000000000000e1 * t519 + 0.84375000000000000000e0 * t522 + t532 + t573 + t626 + 0.84375000000000000000e0 * t629 + 0.37500000000000000000e0 * t632 + 0.11250000000000000000e1 * t634 + t638 + 0.84375000000000000000e0 * t640 + 0.28125000000000000000e0 * t643 + 0.22500000000000000000e1 * t470;
  t654 = t377 * t205;
  t655 = t654 / 0.8e1;
  t658 = t94 * t1 * t3 * t110;
  t660 = t188 * t388;
  t664 = t63 * t64 * t203 * t154;
  t665 = 0.84375000000000000000e0 * t664;
  t666 = t188 * t368;
  t685 = 0.32e2 * t194 * t592;
  t703 = 0.32e2 * t200 * t592;
  t718 = 0.40e2 / 0.81e2 * t579 * t103 * t261 - 0.8e1 / 0.9e1 * t340 * t68 * t135 - 0.16e2 / 0.9e1 * t260 * t10 * t132 * t128 + 0.8e1 / 0.3e1 * t127 * t132 * t68 + 0.8e1 / 0.3e1 * t343 * t132 * t135 + t685 + 0.40e2 / 0.81e2 * t600 * t105 * t276 - 0.8e1 / 0.9e1 * t353 * t70 * t143 + 0.16e2 / 0.9e1 * t275 * t10 * t132 * t140 - 0.8e1 / 0.3e1 * t139 * t132 * t70 - 0.8e1 / 0.3e1 * t356 * t132 * t143 - t703 - 0.8e1 * t343 * t268 * t68 + 0.8e1 * t356 * t268 * t70 + 0.4e1 / 0.9e1 * t191 * t271 - 0.16e2 * t15 * t268 + 0.4e1 / 0.9e1 * t197 * t282 + 0.16e2 * t18 * t268;
  t723 = 0.11250000000000000000e1 * t477 + 0.33750000000000000000e1 * t483 - t492 - t466 / 0.8e1 + t473 / 0.12e2 + 0.33750000000000000000e1 * t480 - t488 + t496 + 0.84375000000000000000e0 * t513 + t528 - t575 / 0.8e1 - t655 + t658 / 0.36e2 + 0.11250000000000000000e1 * t660 + t665 - 0.3e1 / 0.8e1 * t666 - 0.3e1 / 0.32e2 * t63 * t64 * t718 * t55;
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[1] = t646 + t723;

  t724 = t380 * t430;
  t728 = t79 * t9 * t424 * t99;
  t731 = t79 * t429 * t182;
  t744 = 0.75000000000000000000e0 * t724 + 0.56250000000000000000e0 * t728 + 0.28125000000000000000e0 * t731 + t499 + t502 + t506 + 0.28125000000000000000e0 * t516 + 0.16875000000000000000e1 * t519 + 0.56250000000000000000e0 * t522 + t532 + t573 + 0.30000000000000000000e1 * t625 + 0.56250000000000000000e0 * t629 + 0.75000000000000000000e0 * t632 + 0.22500000000000000000e1 * t634 + 0.11250000000000000000e1 * t637 + 0.16875000000000000000e1 * t640 + 0.56250000000000000000e0 * t643 + 0.75000000000000000000e0 * t470;
  t756 = t377 * t229;
  t760 = t63 * t64 * t227 * t154;
  t762 = t188 * t426;
  t773 = t10 * t10;
  t776 = 0.1e1 / t267 / t65;
  t815 = 0.40e2 / 0.81e2 * t579 * t214 * t128 - 0.64e2 / 0.27e2 * t340 * t68 * t10 * t132 - 0.8e1 / 0.27e2 * t400 * t135 + 0.32e2 / 0.9e1 * t127 * t773 * t776 + 0.16e2 / 0.9e1 * t191 * t132 - 0.16e2 / 0.3e1 * t191 * t269 - 0.8e1 / 0.27e2 * t260 * t218 * t128 + 0.8e1 / 0.9e1 * t127 * t410 * t68 + 0.4e1 / 0.9e1 * t405 * t135 + t685 + 0.40e2 / 0.81e2 * t600 * t221 * t140 + 0.64e2 / 0.27e2 * t353 * t70 * t10 * t132 - 0.8e1 / 0.27e2 * t413 * t143 + 0.32e2 / 0.9e1 * t139 * t773 * t776 - 0.16e2 / 0.9e1 * t197 * t132 + 0.16e2 / 0.3e1 * t197 * t269 - 0.8e1 / 0.27e2 * t275 * t224 * t140 + 0.8e1 / 0.9e1 * t139 * t421 * t70 + 0.4e1 / 0.9e1 * t418 * t143 - t703;
  t820 = 0.75000000000000000000e0 * t477 + 0.22500000000000000000e1 * t483 - t492 - t466 / 0.24e2 + t473 / 0.18e2 + 0.22500000000000000000e1 * t480 - t488 + t496 + 0.28125000000000000000e0 * t513 + t528 - t654 / 0.6e1 + t658 / 0.18e2 + 0.22500000000000000000e1 * t660 + 0.11250000000000000000e1 * t664 - t666 / 0.4e1 - t756 / 0.24e2 + 0.28125000000000000000e0 * t760 - t762 / 0.4e1 - 0.3e1 / 0.32e2 * t63 * t64 * t815 * t55;
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[2] = t744 + t820;

  t823 = t79 * t9 * t457 * t99;
  t834 = 0.28125000000000000000e0 * t823 + 0.22500000000000000000e1 * t724 + 0.84375000000000000000e0 * t728 + 0.84375000000000000000e0 * t731 + t499 + t502 + t506 + 0.84375000000000000000e0 * t519 + 0.28125000000000000000e0 * t522 + t532 + t573 + t626 + 0.11250000000000000000e1 * t632 + 0.33750000000000000000e1 * t634 + t638 + 0.25312500000000000000e1 * t640 + 0.84375000000000000000e0 * t643;
  t862 = 0.12e2 * t268 + 0.24e2 * t593;
  t884 = 0.40e2 / 0.81e2 * t579 * t440 * t68 - 0.16e2 / 0.9e1 * t400 * t133 - 0.8e1 / 0.9e1 * t340 * t218 * t68 + 0.8e1 / 0.3e1 * t343 * t132 * t218 + 0.4e1 / 0.3e1 * t191 * t410 + 0.4e1 / 0.9e1 * t127 * t446 * t68 + 0.4e1 / 0.3e1 * t15 * t862 + 0.40e2 / 0.81e2 * t600 * t449 * t70 + 0.16e2 / 0.9e1 * t413 * t133 - 0.8e1 / 0.9e1 * t353 * t224 * t70 - 0.8e1 / 0.3e1 * t356 * t132 * t224 + 0.4e1 / 0.3e1 * t197 * t421 + 0.4e1 / 0.9e1 * t139 * t454 * t70 - 0.4e1 / 0.3e1 * t18 * t862;
  t889 = t188 * t459;
  t891 = 0.37500000000000000000e0 * t477 + 0.11250000000000000000e1 * t483 - t492 + t473 / 0.36e2 + 0.11250000000000000000e1 * t480 - t488 + t496 + t528 - t655 + t658 / 0.12e2 + 0.33750000000000000000e1 * t660 + t665 - t756 / 0.8e1 + 0.84375000000000000000e0 * t760 - 0.3e1 / 0.8e1 * t762 - 0.3e1 / 0.32e2 * t63 * t64 * t884 * t55 - t889 / 0.8e1;
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[3] = t834 + t891;

  t894 = t214 * t214;
  t899 = t218 * t218;
  t905 = 0.24e2 * t268 + 0.24e2 * t593;
  t908 = t221 * t221;
  t913 = t224 * t224;
  t934 = t658 / 0.9e1 + 0.15000000000000000000e1 * t632 - 0.3e1 / 0.32e2 * t63 * t64 * (0.40e2 / 0.81e2 * t579 * t894 - 0.16e2 / 0.9e1 * t400 * t218 + 0.4e1 / 0.3e1 * t127 * t899 + 0.16e2 / 0.9e1 * t191 * t446 + 0.4e1 / 0.3e1 * t15 * t905 + 0.40e2 / 0.81e2 * t600 * t908 - 0.16e2 / 0.9e1 * t413 * t224 + 0.4e1 / 0.3e1 * t139 * t913 + 0.16e2 / 0.9e1 * t197 * t454 - 0.4e1 / 0.3e1 * t18 * t905) * t55 + 0.45000000000000000000e1 * t660 + 0.45000000000000000000e1 * t634 + 0.33750000000000000000e1 * t640 + 0.11250000000000000000e1 * t643 - t756 / 0.4e1 + 0.45000000000000000000e1 * t724 + 0.16875000000000000000e1 * t760 + 0.16875000000000000000e1 * t731;
  t937 = -t488 - t492 + t496 + t499 + t502 + t506 + t528 + t532 + t573 - t889 / 0.2e1 + 0.11250000000000000000e1 * t823;
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[4] = t934 + t937;

#ifndef XC_DONT_COMPILE_MXC

  if(order < 5) return;


#endif

#endif

#endif

#endif

#endif


}

