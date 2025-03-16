#include <Rcpp.h>
#include <unordered_set>
using namespace Rcpp;

// unique()
// [[Rcpp::export]]
NumericVector unique_cpp(NumericVector x) {
  std::unordered_set<double> seen;
  NumericVector result;
  for (auto& val : x) {
    if (seen.insert(val).second) {
      result.push_back(val);
    }
  }
  return result;
}










// XWX_XWSWX_equicorr_cpp()
// [[Rcpp::export]]
List XWX_XWSWX_equicorr_cpp(double rho, int num_leaves, int I, NumericVector n_i, NumericVector nodesis, NumericVector epsilon) {
  double theta = rho/(1-rho);
  NumericVector nodesis_shift = nodesis-1;
  NumericMatrix XWSWX(num_leaves, num_leaves);
  NumericMatrix XWX(num_leaves, num_leaves);

  int index_tracker1 = 0;
  int index_tracker2 = 0;
  int node_gg1;
  int node_gg2;

  for(int i = 0; i < I; i++) {
    int n = n_i(i);
    double phi = theta/(1+n*theta);
    NumericVector nodes_i(n);
    NumericVector tilde_blank_i(n);
    NumericVector epsilon_i(n);
    NumericVector epsilon_pert_i(n);

    index_tracker2 = index_tracker1+n-1;
    // Extract node numbers and epsilon vector of a single group (index i)
    // NTS: Could do the below in one single loop if updating epsilon_i at
    //      gg+1 extra along each line (this format an easier baseline for
    //      adapting to other weights though).
    for (int gg = 0; gg < n; gg++) {
      nodes_i(gg) = nodesis_shift(index_tracker1 + gg);
      epsilon_i(gg) = epsilon(index_tracker1 + gg);
    }
    double phi_sum_eps_i = phi * sum(epsilon_i);
    for (int gg = 0; gg < n; gg++) {
      epsilon_pert_i(gg) = epsilon_i(gg) - phi_sum_eps_i;
    }
    index_tracker1 = index_tracker2+1;


    for (int gg1 = 0; gg1 < n; gg1++) {
      node_gg1 = nodes_i(gg1);
      XWX(node_gg1, node_gg1) += 1;
      for (int gg2 = 0; gg2 < n; gg2++) {
        node_gg2 = nodes_i(gg2);
        XWX(node_gg1, node_gg2) -= phi;
        XWSWX(node_gg1, node_gg2) += epsilon_pert_i(gg1) * epsilon_pert_i(gg2);
      }
    }

  }

  List list_XWX_XWSWX = List::create(XWX, XWSWX);

  return list_XWX_XWSWX;
}

// XWX_XWSWX_XX_equicorr_cpp()
// [[Rcpp::export]]
List XWX_XWSWX_XX_equicorr_cpp(double rho, int num_leaves, int I, NumericVector n_i, NumericVector nodesis, NumericVector epsilon) {
  double theta = rho/(1-rho);
  NumericVector nodesis_shift = nodesis-1;
  NumericMatrix XWSWX(num_leaves, num_leaves);
  NumericMatrix XWX(num_leaves, num_leaves);
  NumericMatrix XX(num_leaves, num_leaves);

  int index_tracker1 = 0;
  int index_tracker2 = 0;
  int node_gg1;
  int node_gg2;

  for(int i = 0; i < I; i++) {
    int n = n_i(i);
    double phi = theta/(1+n*theta);
    NumericVector nodes_i(n);
    NumericVector tilde_blank_i(n);
    NumericVector epsilon_i(n);
    NumericVector epsilon_pert_i(n);

    index_tracker2 = index_tracker1+n-1;
    // Extract node numbers and epsilon vector of a single group (index i)
    // NTS: Could do the below in one single loop if updating epsilon_i at
    //      gg+1 extra along each line (this format an easier baseline for
    //      adapting to other weights though).
    for (int gg = 0; gg < n; gg++) {
      nodes_i(gg) = nodesis_shift(index_tracker1 + gg);
      epsilon_i(gg) = epsilon(index_tracker1 + gg);
    }
    double phi_sum_eps_i = phi * sum(epsilon_i);
    for (int gg = 0; gg < n; gg++) {
      epsilon_pert_i(gg) = epsilon_i(gg) - phi_sum_eps_i;
    }
    index_tracker1 = index_tracker2+1;


    for (int gg1 = 0; gg1 < n; gg1++) {
      node_gg1 = nodes_i(gg1);
      XWX(node_gg1, node_gg1) += 1;
      XX(node_gg1, node_gg1) += 1;
      for (int gg2 = 0; gg2 < n; gg2++) {
        node_gg2 = nodes_i(gg2);
        XWX(node_gg1, node_gg2) -= phi;
        XWSWX(node_gg1, node_gg2) += epsilon_pert_i(gg1) * epsilon_pert_i(gg2);
      }
    }

  }

  List list_XWX_XWSWX_XX = List::create(XWX, XWSWX, XX);

  return list_XWX_XWSWX_XX;
}

// XWX_XWY_equicorr_cpp()
// [[Rcpp::export]]
List XWX_XWY_equicorr_cpp(double rho, int num_leaves, int I, NumericVector n_i, NumericVector nodesis, NumericVector Y) {
  double theta = rho/(1-rho);
  NumericVector nodesis_shift = nodesis-1;
  NumericMatrix XWX(num_leaves, num_leaves);
  NumericVector XWY(num_leaves);

  int index_tracker1 = 0;
  int index_tracker2 = 0;
  int node_gg1;
  int node_gg2;

  for(int i = 0; i < I; i++) {
    int n = n_i(i);
    double phi = theta/(1+n*theta);
    NumericVector nodes_i(n);
    NumericVector y_i(n);

    index_tracker2 = index_tracker1+n-1;
    for (int gg = 0; gg < n; gg++) {
      nodes_i(gg) = nodesis_shift(index_tracker1 + gg);
      y_i(gg) = Y(index_tracker1 + gg);
    }
    index_tracker1 = index_tracker2+1;
    double phi_sum_y_i = phi * sum(y_i);

    // Comp. improvements over commented code (for inv equicorr structure)
    for (int gg1 = 0; gg1 < n; gg1++) {
      node_gg1 = nodes_i(gg1);
      XWX(node_gg1, node_gg1) += 1;
      XWY(node_gg1) += y_i(gg1) - phi_sum_y_i;
      for (int gg2 = 0; gg2 < n; gg2++) {
        node_gg2 = nodes_i(gg2);
        XWX(node_gg1, node_gg2) -= phi;
      }
    }

  }

  List list_XWX_XWY = List::create(XWX, XWY);

  return list_XWX_XWY;
}

// XWX_XWSWX_XX_equicorr_cpp()
// [[Rcpp::export]]
List XWX_XWSWX_XtestX_equicorr_cpp(double rho, int num_leaves, int I, NumericVector n_i, NumericVector nodesis, NumericVector epsilon, int N_test, NumericVector nodesis_test) {
  double theta = rho/(1-rho);
  NumericVector nodesis_shift = nodesis-1;
  NumericVector nodesis_shift_test = nodesis_test-1;
  // Rcpp::Rcout << nodesis << std::endl;
  NumericMatrix XWSWX(num_leaves, num_leaves);
  NumericMatrix XWX(num_leaves, num_leaves);
  NumericMatrix XX(num_leaves, num_leaves);

  int index_tracker1 = 0;
  int index_tracker2 = 0;

  for(int i = 0; i < I; i++) {
    int n = n_i(i);
    double phi = theta/(1+theta*n);
    index_tracker2 = index_tracker1+n-1;
    NumericVector nodes_i(n);
    NumericVector epsilon_i(n);
    NumericVector vec_sub(num_leaves);
    NumericVector eps_sub(num_leaves);
    for (int gg = 0; gg < n; gg++) {
      nodes_i(gg) = nodesis_shift(index_tracker1 + gg);
      epsilon_i(gg) = epsilon(index_tracker1 + gg);
      vec_sub(nodes_i(gg)) += 1;
      eps_sub(nodes_i(gg)) += epsilon_i(gg);
    }
    double sum_epsilon_i = sum(epsilon_i);
    index_tracker1 = index_tracker2+1;


    NumericVector unique_nodes_i;
    unique_nodes_i = unique_cpp(nodes_i);

    for(int k : unique_nodes_i) {
      for(int kk : unique_nodes_i) {
        if (k == kk) {
          XWX(k,kk) += vec_sub(k);
        }
        XWX(k,kk) -= phi*vec_sub(k)*vec_sub(kk);
        XWSWX(k,kk) += eps_sub(k)*eps_sub(kk) - phi*sum_epsilon_i*(eps_sub(k)*vec_sub(kk)+vec_sub(k)*eps_sub(kk)) + pow(phi,2)*pow(sum_epsilon_i,2)*vec_sub(k)*vec_sub(kk);
      }
    }

  }

  for(int i = 0; i < N_test; i++) {
    int loc = nodesis_shift_test(i);
    XX(loc,loc) += 1;
  }

  List list_XWX_XWSWX_XX = List::create(XWX, XWSWX, XX);

  return list_XWX_XWSWX_XX;
}










// XWX_XWSWX_ar1_cpp()
// [[Rcpp::export]]
List XWX_XWSWX_ar1_cpp(double rho, int num_leaves, int I, NumericVector n_i, NumericVector nodesis, NumericVector epsilon) {
  double rho2 = pow(rho,2);
  NumericVector nodesis_shift = nodesis-1;
  // Rcpp::Rcout << nodesis << std::endl;
  NumericMatrix XWSWX(num_leaves, num_leaves);
  NumericMatrix XWX(num_leaves, num_leaves);

  int index_tracker1 = 0;
  int index_tracker2 = 0;
  int node_gg1;
  int node_gg2;

  for(int i = 0; i < I; i++) {
    int n = n_i(i);
    NumericVector nodes_i(n);
    NumericVector tilde_blank_i(n);
    NumericVector epsilon_i(n);
    NumericVector tilde_epsilon_i(n);
    NumericVector bar_epsilon_i(n);
    NumericVector tmrb_epsilon_i(n); // tilde_eps - rho * bar_eps

    index_tracker2 = index_tracker1+n-1;
    // Extract node numbers and epsilon vector of a single group (index i)
    // NTS: Could do the below in one single loop if updating epsilon_i at
    //      gg+1 extra along each line (this format an easier baseline for
    //      adapting to other ARMA processes etc though).
    for (int gg = 0; gg < n; gg++) {
      nodes_i(gg) = nodesis_shift(index_tracker1 + gg);
      epsilon_i(gg) = epsilon(index_tracker1 + gg);
    }
    for (int gg = 0; gg < n; gg++) {
      if (gg == 0) {
        tilde_epsilon_i(gg) = epsilon_i(gg);
        bar_epsilon_i(gg) = epsilon_i(gg+1);
        tilde_blank_i(gg) = 1;
      } else if (gg == n-1) {
        tilde_epsilon_i(gg) = epsilon_i(gg);
        bar_epsilon_i(gg) = epsilon_i(gg-1);
        tilde_blank_i(gg) = 1;
      } else {
        tilde_epsilon_i(gg) = (1+rho2)*epsilon_i(gg);
        bar_epsilon_i(gg) = epsilon_i(gg+1) + epsilon_i(gg-1);
        tilde_blank_i(gg) = 1+rho2;
      }
      tmrb_epsilon_i(gg) = tilde_epsilon_i(gg) - rho * bar_epsilon_i(gg);
    }
    index_tracker1 = index_tracker2+1;

    // Old run for calculating XWX
    //    for (int gg1 = 0; gg1 < n; gg1++) {
    //      for (int gg2 = 0; gg2 < n; gg2++) {
    //        node_gg1 = nodes_i(gg1);
    //        node_gg2 = nodes_i(gg2);
    //        if (gg1 == gg2) {
    //          if (gg1 == 0 || gg1 == n-1) {
    //            XWX(node_gg1, node_gg2) += 1;
    //          } else {
    //            XWX(node_gg1, node_gg2) += 1+rho2;
    //          }
    //        } else if (abs(gg2-gg1) == 1) {
    //          XWX(node_gg1, node_gg2) -= rho;
    //        }
    //      }
    //    }

    for (int gg1 = 0; gg1 < n; gg1++) {
      node_gg1 = nodes_i(gg1);
      XWX(node_gg1, node_gg1) += tilde_blank_i(gg1);
      for (int gg2 = 0; gg2 < n; gg2++) {
        node_gg2 = nodes_i(gg2);

    //   XWSWX(node_gg1, node_gg2) += tilde_epsilon_i(gg1) * tilde_epsilon_i(gg2) +
    //                                 rho2 * bar_epsilon_i(gg1) * bar_epsilon_i(gg2) -
    //                                 rho * tilde_epsilon_i(gg1) * bar_epsilon_i(gg2) -
    //                                 rho * bar_epsilon_i(gg1) * tilde_epsilon_i(gg2);

        XWSWX(node_gg1, node_gg2) += tmrb_epsilon_i(gg1) * tmrb_epsilon_i(gg2);

        if (abs(gg2-gg1) == 1) {
          XWX(node_gg1, node_gg2) -= rho;
        }

      }
    }

  }

  List list_XWX_XWSWX = List::create(XWX, XWSWX);

  return list_XWX_XWSWX;
}

// XWX_XWSWX_XX_ar1_cpp()
// [[Rcpp::export]]
List XWX_XWSWX_XX_ar1_cpp(double rho, int num_leaves, int I, NumericVector n_i, NumericVector nodesis, NumericVector epsilon) {
  double rho2 = pow(rho,2);
  NumericVector nodesis_shift = nodesis-1;
  // Rcpp::Rcout << nodesis << std::endl;
  NumericMatrix XWSWX(num_leaves, num_leaves);
  NumericMatrix XWX(num_leaves, num_leaves);
  NumericMatrix XX(num_leaves, num_leaves);

  int index_tracker1 = 0;
  int index_tracker2 = 0;
  int node_gg1;
  int node_gg2;

  for(int i = 0; i < I; i++) {
    int n = n_i(i);
    NumericVector nodes_i(n);
    NumericVector tilde_blank_i(n);
    NumericVector epsilon_i(n);
    NumericVector tilde_epsilon_i(n);
    NumericVector bar_epsilon_i(n);
    NumericVector tmrb_epsilon_i(n); // tilde_eps - rho * bar_eps

    index_tracker2 = index_tracker1+n-1;
    // Extract node numbers and epsilon vector of a single group (index i)
    // NTS: Could do the below in one single loop if updating epsilon_i at
    //      gg+1 extra along each line (this format an easier baseline for
    //      adapting to other ARMA processes etc though).
    for (int gg = 0; gg < n; gg++) {
      nodes_i(gg) = nodesis_shift(index_tracker1 + gg);
      epsilon_i(gg) = epsilon(index_tracker1 + gg);
    }
    for (int gg = 0; gg < n; gg++) {
      if (gg == 0) {
        tilde_epsilon_i(gg) = epsilon_i(gg);
        bar_epsilon_i(gg) = epsilon_i(gg+1);
        tilde_blank_i(gg) = 1;
      } else if (gg == n-1) {
        tilde_epsilon_i(gg) = epsilon_i(gg);
        bar_epsilon_i(gg) = epsilon_i(gg-1);
        tilde_blank_i(gg) = 1;
      } else {
        tilde_epsilon_i(gg) = (1+rho2)*epsilon_i(gg);
        bar_epsilon_i(gg) = epsilon_i(gg+1) + epsilon_i(gg-1);
        tilde_blank_i(gg) = 1+rho2;
      }
      tmrb_epsilon_i(gg) = tilde_epsilon_i(gg) - rho * bar_epsilon_i(gg);
    }
    index_tracker1 = index_tracker2+1;

    // Old run for calculating XWX
    //    for (int gg1 = 0; gg1 < n; gg1++) {
    //      for (int gg2 = 0; gg2 < n; gg2++) {
    //        node_gg1 = nodes_i(gg1);
    //        node_gg2 = nodes_i(gg2);
    //        if (gg1 == gg2) {
    //          if (gg1 == 0 || gg1 == n-1) {
    //            XWX(node_gg1, node_gg2) += 1;
    //          } else {
    //            XWX(node_gg1, node_gg2) += 1+rho2;
    //          }
    //        } else if (abs(gg2-gg1) == 1) {
    //          XWX(node_gg1, node_gg2) -= rho;
    //        }
    //      }
    //    }

    for (int gg1 = 0; gg1 < n; gg1++) {
      node_gg1 = nodes_i(gg1);
      XWX(node_gg1, node_gg1) += tilde_blank_i(gg1);
      XX(node_gg1, node_gg1) += 1;
      for (int gg2 = 0; gg2 < n; gg2++) {
        node_gg2 = nodes_i(gg2);

        //   XWSWX(node_gg1, node_gg2) += tilde_epsilon_i(gg1) * tilde_epsilon_i(gg2) +
        //                                 rho2 * bar_epsilon_i(gg1) * bar_epsilon_i(gg2) -
        //                                 rho * tilde_epsilon_i(gg1) * bar_epsilon_i(gg2) -
        //                                 rho * bar_epsilon_i(gg1) * tilde_epsilon_i(gg2);

        XWSWX(node_gg1, node_gg2) +=  (gg1) * tmrb_epsilon_i(gg2);

        if (abs(gg2-gg1) == 1) {
          XWX(node_gg1, node_gg2) -= rho;
        }

      }
    }

  }

  List list_XWX_XWSWX_XX = List::create(XWX, XWSWX, XX);

  return list_XWX_XWSWX_XX;
}

// XWX_XWY_ar1_cpp()
// [[Rcpp::export]]
List XWX_XWY_ar1_cpp(double rho, int num_leaves, int I, NumericVector n_i, NumericVector nodesis, NumericVector Y) {
  double rho2 = pow(rho,2);
  NumericVector nodesis_shift = nodesis-1;
  // Rcpp::Rcout << nodesis << std::endl;
  NumericMatrix XWX(num_leaves, num_leaves);
  NumericVector XWY(num_leaves);

  int index_tracker1 = 0;
  int index_tracker2 = 0;
  int node_gg1;
  // int node_gg2;

  for(int i = 0; i < I; i++) {
    int n = n_i(i);
    NumericVector nodes_i(n);
    NumericVector y_i(n);
    NumericVector tilde_blank_i(n);

    index_tracker2 = index_tracker1+n-1;
    for (int gg = 0; gg < n; gg++) {
      nodes_i(gg) = nodesis_shift(index_tracker1 + gg);
      y_i(gg) = Y(index_tracker1 + gg);
      if (gg == 0 || gg == n-1) {
        tilde_blank_i(gg) = 1;
      } else {
        tilde_blank_i(gg) = 1+rho2;
      }
    }
    index_tracker1 = index_tracker2+1;

    // Comp. improvements over commented code (for inv AR(1) structure)
    // NTS: Commented code more generalisable for other working weights
    for (int gg1 = 0; gg1 < n; gg1++) {
      node_gg1 = nodes_i(gg1);
      XWX(node_gg1, node_gg1) += tilde_blank_i(gg1);
      if ((gg1 != 0) & (gg1 != n-1)) {
        XWX(node_gg1, nodes_i(gg1-1)) -= rho;
        XWX(node_gg1, nodes_i(gg1+1)) -= rho;
        XWY(node_gg1) += (1+rho2)*y_i(gg1) - rho*(y_i(gg1-1)+y_i(gg1+1));
      } else if (gg1 == 0) {
        XWX(node_gg1, nodes_i(gg1+1)) -= rho;
        XWY(node_gg1) += y_i(gg1) - rho*y_i(gg1+1);
      } else if (gg1 == n-1) {
        XWX(node_gg1, nodes_i(gg1-1)) -= rho;
        XWY(node_gg1) += y_i(gg1) - rho*y_i(gg1-1);
      }
    }

    //    for (int gg1 = 0; gg1 < n; gg1++) {
    //      for (int gg2 = 0; gg2 < n; gg2++) {
    //        node_gg1 = nodes_i(gg1);
    //        node_gg2 = nodes_i(gg2);
    //        if (gg1 == gg2) {
    //          if (gg1 == 0 || gg1 == n-1) {
    //            XWX(node_gg1, node_gg2) += 1;
    //          } else {
    //            XWX(node_gg1, node_gg2) += 1+rho2;
    //          }
    //        } else if (abs(gg2-gg1) == 1) {
    //          XWX(node_gg1, node_gg2) -= rho;
    //        }
    //      }
    //    }







//    for (int gg1 = 0; gg1 < n; gg1++) {
//      node_gg1 = nodes_i(gg1);
//      for (int gg2 = 0; gg2 < n; gg2++) {
//        node_gg2 = nodes_i(gg2);
//        if (gg1 == gg2) {
//          if (gg1 == 0 || gg1 == n-1) {
//            XWX(node_gg1, node_gg2) += 1;
//          } else {
//            XWX(node_gg1, node_gg2) += 1+rho2;
//          }
//        } else if (abs(gg2-gg1) == 1) {
//          XWX(node_gg1, node_gg2) -= rho;
//        }
//      }
//
//      if (gg1 == 0) {
//        XWY(node_gg1) += y_i(gg1) - rho*y_i(gg1+1);
//      } else if (gg1 == n-1) {
//        XWY(node_gg1) += y_i(gg1) - rho*y_i(gg1-1);
//      } else {
//        XWY(node_gg1) += (1+rho2)*y_i(gg1) - rho*(y_i(gg1-1)+y_i(gg1+1));
//      }
//    }








//    for (int gg1 = 0; gg1 < n; gg1++) {
//      node_gg1 = nodes_i(gg1);
//      for (int gg2 = 0; gg2 < n; gg2++) {
//        node_gg2 = nodes_i(gg2);
//        if (gg1 == gg2) {
//          if (gg1 == 0 || gg1 == n-1) {
//            XWX(node_gg1, node_gg2) += 1;
//          } else {
//            XWX(node_gg1, node_gg2) += 1+rho2;
//          }
//        } else if (abs(gg2-gg1) == 1) {
//          XWX(node_gg1, node_gg2) -= rho;
//        }
//      }
//
//      if (gg1 == 0) {
//        XWY(node_gg1) += y_i(gg1) - rho*y_i(gg1+1);
//      } else if (gg1 == n-1) {
//        XWY(node_gg1) += y_i(gg1) - rho*y_i(gg1-1);
//      } else {
//        XWY(node_gg1) += (1+rho2)*y_i(gg1) - rho*(y_i(gg1-1)+y_i(gg1+1));
//      }
//
//    }

  }

  List list_XWX_XWY = List::create(XWX, XWY);

  return list_XWX_XWY;
}


// XWX_XWSWX_XX_ar1_cpp()
// [[Rcpp::export]]
List XWX_XWSWX_XX_OLD_ar1_cpp(double rho, int num_leaves, int I, NumericVector n_i, NumericVector nodesis, NumericVector epsilon) {
  double rho2 = pow(rho,2);
  NumericVector nodesis_shift = nodesis-1;
  // Rcpp::Rcout << nodesis << std::endl;
  NumericMatrix XWSWX(num_leaves, num_leaves);
  NumericMatrix XWX(num_leaves, num_leaves);
  NumericMatrix XX(num_leaves, num_leaves);

  int index_tracker1 = 0;
  int index_tracker2 = 0;
  int node_gg1;
  int node_gg2;

  for(int i = 0; i < I; i++) {
    int n = n_i(i);
    NumericVector nodes_i(n);
    NumericVector tilde_blank_i(n);
    NumericVector epsilon_i(n);
    NumericVector tilde_epsilon_i(n);
    NumericVector bar_epsilon_i(n);
    NumericVector tmrb_epsilon_i(n); // tilde_eps - rho * bar_eps

    index_tracker2 = index_tracker1+n-1;
    // Extract node numbers and epsilon vector of a single group (index i)
    // NTS: Could do the below in one single loop if updating epsilon_i at
    //      gg+1 extra along each line (this format an easier baseline for
    //      adapting to other ARMA processes etc though).
    for (int gg = 0; gg < n; gg++) {
      nodes_i(gg) = nodesis_shift(index_tracker1 + gg);
      epsilon_i(gg) = epsilon(index_tracker1 + gg);
    }
    for (int gg = 0; gg < n; gg++) {
      if (gg == 0) {
        tilde_epsilon_i(gg) = epsilon_i(gg);
        bar_epsilon_i(gg) = epsilon_i(gg+1);
        tilde_blank_i(gg) = 1;
      } else if (gg == n-1) {
        tilde_epsilon_i(gg) = epsilon_i(gg);
        bar_epsilon_i(gg) = epsilon_i(gg-1);
        tilde_blank_i(gg) = 1;
      } else {
        tilde_epsilon_i(gg) = (1+rho2)*epsilon_i(gg);
        bar_epsilon_i(gg) = epsilon_i(gg+1) + epsilon_i(gg-1);
        tilde_blank_i(gg) = 1+rho2;
      }
      tmrb_epsilon_i(gg) = tilde_epsilon_i(gg) - rho * bar_epsilon_i(gg);
    }
    index_tracker1 = index_tracker2+1;

    // Old run for calculating XWX
    //    for (int gg1 = 0; gg1 < n; gg1++) {
    //      for (int gg2 = 0; gg2 < n; gg2++) {
    //        node_gg1 = nodes_i(gg1);
    //        node_gg2 = nodes_i(gg2);
    //        if (gg1 == gg2) {
    //          if (gg1 == 0 || gg1 == n-1) {
    //            XWX(node_gg1, node_gg2) += 1;
    //          } else {
    //            XWX(node_gg1, node_gg2) += 1+rho2;
    //          }
    //        } else if (abs(gg2-gg1) == 1) {
    //          XWX(node_gg1, node_gg2) -= rho;
    //        }
    //      }
    //    }

    for (int gg1 = 0; gg1 < n; gg1++) {
      node_gg1 = nodes_i(gg1);
      XWX(node_gg1, node_gg1) += tilde_blank_i(gg1);
      for (int gg2 = 0; gg2 < n; gg2++) {
        node_gg2 = nodes_i(gg2);

        //   XWSWX(node_gg1, node_gg2) += tilde_epsilon_i(gg1) * tilde_epsilon_i(gg2) +
        //                                 rho2 * bar_epsilon_i(gg1) * bar_epsilon_i(gg2) -
        //                                 rho * tilde_epsilon_i(gg1) * bar_epsilon_i(gg2) -
        //                                 rho * bar_epsilon_i(gg1) * tilde_epsilon_i(gg2);

        XWSWX(node_gg1, node_gg2) += tmrb_epsilon_i(gg1) * tmrb_epsilon_i(gg2);

        if (abs(gg2-gg1) == 1) {
          XWX(node_gg1, node_gg2) -= rho;
        }

      }
    }

  }

  List list_XWX_XWSWX = List::create(XWX, XWSWX);

  return list_XWX_XWSWX;
}


// XWX_XWSWX_XX_ar1_cpp()
// [[Rcpp::export]]
List XWX_XWSWX_XtestX_ar1_cpp(double rho, int num_leaves, int I, NumericVector n_i, NumericVector nodesis, NumericVector epsilon, int N_test, NumericVector nodesis_test) {
  double rho2 = pow(rho,2);
  NumericVector nodesis_shift = nodesis-1;
  NumericVector nodesis_shift_test = nodesis_test-1;
  // Rcpp::Rcout << nodesis << std::endl;
  NumericMatrix XWSWX(num_leaves, num_leaves);
  NumericMatrix XWX(num_leaves, num_leaves);
  NumericMatrix XX(num_leaves, num_leaves);

  int index_tracker1 = 0;
  int index_tracker2 = 0;
  int node_gg1;
  int node_gg2;

  for(int i = 0; i < I; i++) {
    int n = n_i(i);
    NumericVector nodes_i(n);
    NumericVector tilde_blank_i(n);
    NumericVector epsilon_i(n);
    NumericVector tilde_epsilon_i(n);
    NumericVector bar_epsilon_i(n);
    NumericVector tmrb_epsilon_i(n); // tilde_eps - rho * bar_eps

    index_tracker2 = index_tracker1+n-1;
    // Extract node numbers and epsilon vector of a single group (index i)
    // NTS: Could do the below in one single loop if updating epsilon_i at
    //      gg+1 extra along each line (this format an easier baseline for
    //      adapting to other ARMA processes etc though).
    for (int gg = 0; gg < n; gg++) {
      nodes_i(gg) = nodesis_shift(index_tracker1 + gg);
      epsilon_i(gg) = epsilon(index_tracker1 + gg);
    }
    for (int gg = 0; gg < n; gg++) {
      if (gg == 0) {
        tilde_epsilon_i(gg) = epsilon_i(gg);
        bar_epsilon_i(gg) = epsilon_i(gg+1);
        tilde_blank_i(gg) = 1;
      } else if (gg == n-1) {
        tilde_epsilon_i(gg) = epsilon_i(gg);
        bar_epsilon_i(gg) = epsilon_i(gg-1);
        tilde_blank_i(gg) = 1;
      } else {
        tilde_epsilon_i(gg) = (1+rho2)*epsilon_i(gg);
        bar_epsilon_i(gg) = epsilon_i(gg+1) + epsilon_i(gg-1);
        tilde_blank_i(gg) = 1+rho2;
      }
      tmrb_epsilon_i(gg) = tilde_epsilon_i(gg) - rho * bar_epsilon_i(gg);
    }
    index_tracker1 = index_tracker2+1;

    // Old run for calculating XWX
    //    for (int gg1 = 0; gg1 < n; gg1++) {
    //      for (int gg2 = 0; gg2 < n; gg2++) {
    //        node_gg1 = nodes_i(gg1);
    //        node_gg2 = nodes_i(gg2);
    //        if (gg1 == gg2) {
    //          if (gg1 == 0 || gg1 == n-1) {
    //            XWX(node_gg1, node_gg2) += 1;
    //          } else {
    //            XWX(node_gg1, node_gg2) += 1+rho2;
    //          }
    //        } else if (abs(gg2-gg1) == 1) {
    //          XWX(node_gg1, node_gg2) -= rho;
    //        }
    //      }
    //    }

    for (int gg1 = 0; gg1 < n; gg1++) {
      node_gg1 = nodes_i(gg1);
      XWX(node_gg1, node_gg1) += tilde_blank_i(gg1);
      XX(node_gg1, node_gg1) += 1;
      for (int gg2 = 0; gg2 < n; gg2++) {
        node_gg2 = nodes_i(gg2);

        //   XWSWX(node_gg1, node_gg2) += tilde_epsilon_i(gg1) * tilde_epsilon_i(gg2) +
        //                                 rho2 * bar_epsilon_i(gg1) * bar_epsilon_i(gg2) -
        //                                 rho * tilde_epsilon_i(gg1) * bar_epsilon_i(gg2) -
        //                                 rho * bar_epsilon_i(gg1) * tilde_epsilon_i(gg2);

        XWSWX(node_gg1, node_gg2) += tmrb_epsilon_i(gg1) * tmrb_epsilon_i(gg2);

        if (abs(gg2-gg1) == 1) {
          XWX(node_gg1, node_gg2) -= rho;
        }

      }
    }

  }

  for(int i = 0; i < N_test; i++) {
    int loc = nodesis_shift_test(i);
    XX(loc,loc) += 1;
  }

  List list_XWX_XWSWX_XX = List::create(XWX, XWSWX, XX);

  return list_XWX_XWSWX_XX;
}






