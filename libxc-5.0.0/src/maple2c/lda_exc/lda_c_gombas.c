/* 
  This file was generated automatically with ./scripts/maple2c.pl.
  Do not edit this file directly as it can be overwritten!!

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.

  Maple version     : Maple 2016 (X86 64 LINUX)
  Maple source      : ./maple/lda_exc/lda_c_gombas.mpl
  Type of functional: lda_exc
*/

#define maple2c_order 4
#define MAPLE2C_FLAGS (XC_FLAGS_I_HAVE_EXC | XC_FLAGS_I_HAVE_VXC | XC_FLAGS_I_HAVE_FXC | XC_FLAGS_I_HAVE_KXC | XC_FLAGS_I_HAVE_LXC)


static inline void
func_unpol(const xc_func_type *p, int order, const double *rho, double *zk, LDA_OUT_PARAMS_NO_EXC(double *))
{

#ifndef XC_DONT_COMPILE_EXC
  double t1, t2, t4, t6, t7, t9, t10;

#ifndef XC_DONT_COMPILE_VXC
  double t11, t12, t14, t15, t18, t22, t23, t24;
  double t25;

#ifndef XC_DONT_COMPILE_FXC
  double t32, t33, t35, t36, t39, t40, t44, t47;
  double t48, t49, t51, t52, t53, t54, t56;

#ifndef XC_DONT_COMPILE_KXC
  double t65, t66, t67, t68, t69, t71, t73, t74;
  double t77, t78, t80, t83, t84, t85, t87, t88;
  double t90, t93, t94, t95, t97, t99;

#ifndef XC_DONT_COMPILE_LXC
  double t114, t127, t148, t159;
#endif

#endif

#endif

#endif

#endif



  t1 = POW_1_3(rho[0]);
  t2 = 0.1e1 / t1;
  t4 = 0.1e1 + 0.56200000000000000000e-1 * t2;
  t6 = 0.357e-1 / t4;
  t7 = t2 + 0.239e1;
  t9 = log(t7 * t1);
  t10 = 0.311e-1 * t9;
  if(zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    zk[0] = -t6 - t10;

#ifndef XC_DONT_COMPILE_VXC

  if(order < 1) return;


  t11 = t4 * t4;
  t12 = 0.1e1 / t11;
  t14 = 0.1e1 / t1 / rho[0];
  t15 = t12 * t14;
  t18 = t1 * t1;
  t22 = -0.1e1 / rho[0] / 0.3e1 + t7 / t18 / 0.3e1;
  t23 = 0.1e1 / t7;
  t24 = t22 * t23;
  t25 = t24 * t2;
  if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    vrho[0] = -t6 - t10 + rho[0] * (-0.66877999999999999999e-3 * t15 - 0.311e-1 * t25);

#ifndef XC_DONT_COMPILE_FXC

  if(order < 2) return;


  t32 = 0.1e1 / t11 / t4;
  t33 = rho[0] * rho[0];
  t35 = 0.1e1 / t18 / t33;
  t36 = t32 * t35;
  t39 = 0.1e1 / t1 / t33;
  t40 = t12 * t39;
  t44 = 0.1e1 / t18 / rho[0];
  t47 = 0.2e1 / 0.9e1 / t33 - 0.2e1 / 0.9e1 * t7 * t44;
  t48 = t47 * t23;
  t49 = t48 * t2;
  t51 = t7 * t7;
  t52 = 0.1e1 / t51;
  t53 = t22 * t52;
  t54 = t53 * t44;
  t56 = t24 * t14;
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[0] = -0.13375600000000000000e-2 * t15 - 0.622e-1 * t25 + rho[0] * (-0.25056957333333333333e-4 * t36 + 0.89170666666666666665e-3 * t40 - 0.311e-1 * t49 - 0.10366666666666666667e-1 * t54 + 0.10366666666666666667e-1 * t56);

#ifndef XC_DONT_COMPILE_KXC

  if(order < 3) return;


  t65 = t11 * t11;
  t66 = 0.1e1 / t65;
  t67 = t33 * t33;
  t68 = 0.1e1 / t67;
  t69 = t66 * t68;
  t71 = t33 * rho[0];
  t73 = 0.1e1 / t18 / t71;
  t74 = t32 * t73;
  t77 = 0.1e1 / t1 / t71;
  t78 = t12 * t77;
  t80 = 0.1e1 / t71;
  t83 = 0.10e2 / 0.27e2 * t7 * t35 - 0.10e2 / 0.27e2 * t80;
  t84 = t83 * t23;
  t85 = t84 * t2;
  t87 = t47 * t52;
  t88 = t87 * t44;
  t90 = t48 * t14;
  t93 = 0.1e1 / t51 / t7;
  t94 = t22 * t93;
  t95 = t94 * t80;
  t97 = t53 * t35;
  t99 = t24 * t39;
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[0] = -0.75170871999999999999e-4 * t36 + 0.26751200000000000000e-2 * t40 - 0.933e-1 * t49 - 0.31100000000000000000e-1 * t54 + 0.31100000000000000000e-1 * t56 + rho[0] * (-0.14082010021333333333e-5 * t69 + 0.10022782933333333333e-3 * t74 - 0.20806488888888888888e-2 * t78 - 0.311e-1 * t85 - 0.20733333333333333334e-1 * t88 + 0.20733333333333333334e-1 * t90 - 0.69111111111111111113e-2 * t95 + 0.20733333333333333334e-1 * t97 - 0.13822222222222222223e-1 * t99);

#ifndef XC_DONT_COMPILE_LXC

  if(order < 4) return;


  t114 = t67 * rho[0];
  t127 = 0.1e1 / t1 / t67;
  t148 = t51 * t51;
  t159 = -0.10552119509319111111e-6 / t65 / t4 / t1 / t114 + 0.11265608017066666666e-4 * t66 / t114 - 0.44545701925925925925e-3 * t32 / t18 / t67 + 0.69354962962962962960e-2 * t12 * t127 - 0.311e-1 * (-0.80e2 / 0.81e2 * t7 * t73 + 0.80e2 / 0.81e2 * t68) * t23 * t2 - 0.31100000000000000001e-1 * t83 * t52 * t44 + 0.31100000000000000001e-1 * t84 * t14 - 0.20733333333333333334e-1 * t47 * t93 * t80 + 0.62200000000000000002e-1 * t87 * t35 - 0.41466666666666666668e-1 * t48 * t39 - 0.69111111111111111113e-2 * t22 / t148 * t127 + 0.34555555555555555557e-1 * t94 * t68 - 0.59896296296296296299e-1 * t53 * t73 + 0.32251851851851851854e-1 * t24 * t77;
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[0] = -0.56328040085333333332e-5 * t69 + 0.40091131733333333332e-3 * t74 - 0.83225955555555555555e-2 * t78 - 0.1244e0 * t85 - 0.82933333333333333334e-1 * t88 + 0.82933333333333333334e-1 * t90 - 0.27644444444444444444e-1 * t95 + 0.82933333333333333334e-1 * t97 - 0.55288888888888888890e-1 * t99 + rho[0] * t159;

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
  double t1, t2, t4, t6, t7, t9, t10;

#ifndef XC_DONT_COMPILE_VXC
  double t11, t12, t14, t15, t18, t22, t23, t24;
  double t25;

#ifndef XC_DONT_COMPILE_FXC
  double t32, t33, t35, t36, t39, t40, t44, t47;
  double t48, t49, t51, t52, t53, t54, t56;

#ifndef XC_DONT_COMPILE_KXC
  double t65, t66, t67, t68, t69, t71, t73, t74;
  double t77, t78, t80, t83, t84, t85, t87, t88;
  double t90, t93, t94, t95, t97, t99;

#ifndef XC_DONT_COMPILE_LXC
  double t114, t127, t148, t159;
#endif

#endif

#endif

#endif

#endif



  t1 = POW_1_3(rho[0]);
  t2 = 0.1e1 / t1;
  t4 = 0.1e1 + 0.56200000000000000000e-1 * t2;
  t6 = 0.357e-1 / t4;
  t7 = t2 + 0.239e1;
  t9 = log(t7 * t1);
  t10 = 0.311e-1 * t9;
  if(zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    zk[0] = -t6 - t10;

#ifndef XC_DONT_COMPILE_VXC

  if(order < 1) return;


  t11 = t4 * t4;
  t12 = 0.1e1 / t11;
  t14 = 0.1e1 / t1 / rho[0];
  t15 = t12 * t14;
  t18 = t1 * t1;
  t22 = -0.1e1 / rho[0] / 0.3e1 + t7 / t18 / 0.3e1;
  t23 = 0.1e1 / t7;
  t24 = t22 * t23;
  t25 = t24 * t2;
  if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    vrho[0] = -t6 - t10 + rho[0] * (-0.66877999999999999999e-3 * t15 - 0.311e-1 * t25);

#ifndef XC_DONT_COMPILE_FXC

  if(order < 2) return;


  t32 = 0.1e1 / t11 / t4;
  t33 = rho[0] * rho[0];
  t35 = 0.1e1 / t18 / t33;
  t36 = t32 * t35;
  t39 = 0.1e1 / t1 / t33;
  t40 = t12 * t39;
  t44 = 0.1e1 / t18 / rho[0];
  t47 = 0.2e1 / 0.9e1 / t33 - 0.2e1 / 0.9e1 * t7 * t44;
  t48 = t47 * t23;
  t49 = t48 * t2;
  t51 = t7 * t7;
  t52 = 0.1e1 / t51;
  t53 = t22 * t52;
  t54 = t53 * t44;
  t56 = t24 * t14;
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[0] = -0.13375600000000000000e-2 * t15 - 0.622e-1 * t25 + rho[0] * (-0.25056957333333333333e-4 * t36 + 0.89170666666666666665e-3 * t40 - 0.311e-1 * t49 - 0.10366666666666666667e-1 * t54 + 0.10366666666666666667e-1 * t56);

#ifndef XC_DONT_COMPILE_KXC

  if(order < 3) return;


  t65 = t11 * t11;
  t66 = 0.1e1 / t65;
  t67 = t33 * t33;
  t68 = 0.1e1 / t67;
  t69 = t66 * t68;
  t71 = t33 * rho[0];
  t73 = 0.1e1 / t18 / t71;
  t74 = t32 * t73;
  t77 = 0.1e1 / t1 / t71;
  t78 = t12 * t77;
  t80 = 0.1e1 / t71;
  t83 = 0.10e2 / 0.27e2 * t7 * t35 - 0.10e2 / 0.27e2 * t80;
  t84 = t83 * t23;
  t85 = t84 * t2;
  t87 = t47 * t52;
  t88 = t87 * t44;
  t90 = t48 * t14;
  t93 = 0.1e1 / t51 / t7;
  t94 = t22 * t93;
  t95 = t94 * t80;
  t97 = t53 * t35;
  t99 = t24 * t39;
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[0] = -0.75170871999999999999e-4 * t36 + 0.26751200000000000000e-2 * t40 - 0.933e-1 * t49 - 0.31100000000000000000e-1 * t54 + 0.31100000000000000000e-1 * t56 + rho[0] * (-0.14082010021333333333e-5 * t69 + 0.10022782933333333333e-3 * t74 - 0.20806488888888888888e-2 * t78 - 0.311e-1 * t85 - 0.20733333333333333334e-1 * t88 + 0.20733333333333333334e-1 * t90 - 0.69111111111111111113e-2 * t95 + 0.20733333333333333334e-1 * t97 - 0.13822222222222222223e-1 * t99);

#ifndef XC_DONT_COMPILE_LXC

  if(order < 4) return;


  t114 = t67 * rho[0];
  t127 = 0.1e1 / t1 / t67;
  t148 = t51 * t51;
  t159 = -0.10552119509319111111e-6 / t65 / t4 / t1 / t114 + 0.11265608017066666666e-4 * t66 / t114 - 0.44545701925925925925e-3 * t32 / t18 / t67 + 0.69354962962962962960e-2 * t12 * t127 - 0.311e-1 * (-0.80e2 / 0.81e2 * t7 * t73 + 0.80e2 / 0.81e2 * t68) * t23 * t2 - 0.31100000000000000001e-1 * t83 * t52 * t44 + 0.31100000000000000001e-1 * t84 * t14 - 0.20733333333333333334e-1 * t47 * t93 * t80 + 0.62200000000000000002e-1 * t87 * t35 - 0.41466666666666666668e-1 * t48 * t39 - 0.69111111111111111113e-2 * t22 / t148 * t127 + 0.34555555555555555557e-1 * t94 * t68 - 0.59896296296296296299e-1 * t53 * t73 + 0.32251851851851851854e-1 * t24 * t77;
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[0] = -0.56328040085333333332e-5 * t69 + 0.40091131733333333332e-3 * t74 - 0.83225955555555555555e-2 * t78 - 0.1244e0 * t85 - 0.82933333333333333334e-1 * t88 + 0.82933333333333333334e-1 * t90 - 0.27644444444444444444e-1 * t95 + 0.82933333333333333334e-1 * t97 - 0.55288888888888888890e-1 * t99 + rho[0] * t159;

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
  double t1, t2, t3, t5, t7, t8, t10, t11;

#ifndef XC_DONT_COMPILE_VXC
  double t12, t13, t15, t16, t19, t23, t24, t25;
  double t26;

#ifndef XC_DONT_COMPILE_FXC
  double t33, t34, t36, t37, t40, t41, t45, t48;
  double t49, t50, t52, t53, t54, t55, t57;

#ifndef XC_DONT_COMPILE_KXC
  double t66, t67, t68, t69, t70, t72, t74, t75;
  double t78, t79, t81, t84, t85, t86, t88, t89;
  double t91, t94, t95, t96, t98, t100;

#ifndef XC_DONT_COMPILE_LXC
  double t115, t128, t149, t160;
#endif

#endif

#endif

#endif

#endif



  t1 = rho[0] + rho[1];
  t2 = POW_1_3(t1);
  t3 = 0.1e1 / t2;
  t5 = 0.1e1 + 0.56200000000000000000e-1 * t3;
  t7 = 0.357e-1 / t5;
  t8 = t3 + 0.239e1;
  t10 = log(t8 * t2);
  t11 = 0.311e-1 * t10;
  if(zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    zk[0] = -t7 - t11;

#ifndef XC_DONT_COMPILE_VXC

  if(order < 1) return;


  t12 = t5 * t5;
  t13 = 0.1e1 / t12;
  t15 = 0.1e1 / t2 / t1;
  t16 = t13 * t15;
  t19 = t2 * t2;
  t23 = -0.1e1 / t1 / 0.3e1 + t8 / t19 / 0.3e1;
  t24 = 0.1e1 / t8;
  t25 = t23 * t24;
  t26 = t25 * t3;
  if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    vrho[0] = -t7 - t11 + t1 * (-0.66877999999999999999e-3 * t16 - 0.311e-1 * t26);

  if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    vrho[1] = vrho[0];

#ifndef XC_DONT_COMPILE_FXC

  if(order < 2) return;


  t33 = 0.1e1 / t12 / t5;
  t34 = t1 * t1;
  t36 = 0.1e1 / t19 / t34;
  t37 = t33 * t36;
  t40 = 0.1e1 / t2 / t34;
  t41 = t13 * t40;
  t45 = 0.1e1 / t19 / t1;
  t48 = 0.2e1 / 0.9e1 / t34 - 0.2e1 / 0.9e1 * t8 * t45;
  t49 = t48 * t24;
  t50 = t49 * t3;
  t52 = t8 * t8;
  t53 = 0.1e1 / t52;
  t54 = t23 * t53;
  t55 = t54 * t45;
  t57 = t25 * t15;
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[0] = -0.13375600000000000000e-2 * t16 - 0.622e-1 * t26 + t1 * (-0.25056957333333333333e-4 * t37 + 0.89170666666666666665e-3 * t41 - 0.311e-1 * t50 - 0.10366666666666666667e-1 * t55 + 0.10366666666666666667e-1 * t57);

  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[1] = v2rho2[0];

  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[2] = v2rho2[1];

#ifndef XC_DONT_COMPILE_KXC

  if(order < 3) return;


  t66 = t12 * t12;
  t67 = 0.1e1 / t66;
  t68 = t34 * t34;
  t69 = 0.1e1 / t68;
  t70 = t67 * t69;
  t72 = t34 * t1;
  t74 = 0.1e1 / t19 / t72;
  t75 = t33 * t74;
  t78 = 0.1e1 / t2 / t72;
  t79 = t13 * t78;
  t81 = 0.1e1 / t72;
  t84 = 0.10e2 / 0.27e2 * t8 * t36 - 0.10e2 / 0.27e2 * t81;
  t85 = t84 * t24;
  t86 = t85 * t3;
  t88 = t48 * t53;
  t89 = t88 * t45;
  t91 = t49 * t15;
  t94 = 0.1e1 / t52 / t8;
  t95 = t23 * t94;
  t96 = t95 * t81;
  t98 = t54 * t36;
  t100 = t25 * t40;
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[0] = -0.75170871999999999999e-4 * t37 + 0.26751200000000000000e-2 * t41 - 0.933e-1 * t50 - 0.31100000000000000000e-1 * t55 + 0.31100000000000000000e-1 * t57 + t1 * (-0.14082010021333333333e-5 * t70 + 0.10022782933333333333e-3 * t75 - 0.20806488888888888888e-2 * t79 - 0.311e-1 * t86 - 0.20733333333333333334e-1 * t89 + 0.20733333333333333334e-1 * t91 - 0.69111111111111111113e-2 * t96 + 0.20733333333333333334e-1 * t98 - 0.13822222222222222223e-1 * t100);

  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[1] = v3rho3[0];

  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[2] = v3rho3[1];

  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[3] = v3rho3[2];

#ifndef XC_DONT_COMPILE_LXC

  if(order < 4) return;


  t115 = t68 * t1;
  t128 = 0.1e1 / t2 / t68;
  t149 = t52 * t52;
  t160 = -0.10552119509319111111e-6 / t66 / t5 / t2 / t115 + 0.11265608017066666666e-4 * t67 / t115 - 0.44545701925925925925e-3 * t33 / t19 / t68 + 0.69354962962962962960e-2 * t13 * t128 - 0.311e-1 * (-0.80e2 / 0.81e2 * t8 * t74 + 0.80e2 / 0.81e2 * t69) * t24 * t3 - 0.31100000000000000001e-1 * t84 * t53 * t45 + 0.31100000000000000001e-1 * t85 * t15 - 0.20733333333333333334e-1 * t48 * t94 * t81 + 0.62200000000000000002e-1 * t88 * t36 - 0.41466666666666666668e-1 * t49 * t40 - 0.69111111111111111113e-2 * t23 / t149 * t128 + 0.34555555555555555557e-1 * t95 * t69 - 0.59896296296296296299e-1 * t54 * t74 + 0.32251851851851851854e-1 * t25 * t78;
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[0] = -0.56328040085333333332e-5 * t70 + 0.40091131733333333332e-3 * t75 - 0.83225955555555555555e-2 * t79 - 0.1244e0 * t86 - 0.82933333333333333334e-1 * t89 + 0.82933333333333333334e-1 * t91 - 0.27644444444444444444e-1 * t96 + 0.82933333333333333334e-1 * t98 - 0.55288888888888888890e-1 * t100 + t1 * t160;

  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[1] = v4rho4[0];

  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[2] = v4rho4[1];

  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[3] = v4rho4[2];

  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[4] = v4rho4[3];

#ifndef XC_DONT_COMPILE_MXC

  if(order < 5) return;


#endif

#endif

#endif

#endif

#endif


}
