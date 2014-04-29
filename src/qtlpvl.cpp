// [[Rcpp::depends(RcppEigen)]]

#include "../inst/include/qtlpvl.h"

inline MatrixXd AtA(const MatrixXd& A) {
  int n(A.cols());
  return MatrixXd(n,n).setZero().selfadjointView<Lower>()
    .rankUpdate(A.adjoint());
}

inline ArrayXd Dplus(const ArrayXd& d, const double threshold) {
  ArrayXd di(d.size());
  double comp(d.maxCoeff() * threshold);
  for (int j = 0; j < d.size(); ++j) di[j] = (d[j] < comp) ? 0. : 1./d[j];
  return di;
}

//' Calculate the residual of lm(Y ~ X)
//'
//' An implementation for calculate residual of linear
//' regression. It's ' not as fast as lm_resid_llt (about 50% slower),
//' but is rank revealing, so should be safe for linear correlated X.
//'
//' @param X A model matrix
//' @param Y The response matrix
//' @param threshold Eigen decomposition is used in calculation. An
//' eigen value smaller than threshold * largest eigen value is
//' considered zero.
//' @return The residual matrix
// [[Rcpp::export]]
MatrixXd lm_resid_svd(const MapMatd& X,
		     const MapMatd& Y,
		     const double threshold=1e-7){
  JacobiSVD<MatrixXd>  UDV(X.jacobiSvd(ComputeThinU|ComputeThinV));
  MatrixXd             VDi(UDV.matrixV() * 
			   Dplus(UDV.singularValues().array(), threshold).matrix().asDiagonal());
  return Y - X * VDi * UDV.matrixU().adjoint() * Y;
}

//' Calculate the residual of lm(Y ~ X)
//' @inheritParams lm_resid_svd
//' @return The residual matrix
// [[Rcpp::export]]
MatrixXd lm_resid_llt(const MapMatd& X,
                     const MapMatd& Y){
  const LLT<MatrixXd> llt(AtA(X));
  return Y - X * llt.solve(X.adjoint() * Y);
}


//' Calculate the residual of lm(Y ~ X)
//' @inheritParams lm_resid_svd
//' @return The residual matrix
// [[Rcpp::export]]
MatrixXd lm_resid_qr(const MapMatd& X,
		    const MapMatd& Y){
  Eigen::HouseholderQR<MatrixXd> QR(X);
  return Y - X * QR.solve(Y);
}


//' Calculate the residual of lm(Y ~ X)
//' @inheritParams lm_resid_svd
//' @return The residual matrix
// [[Rcpp::export]]
MatrixXd lm_resid_symmEigen(const MapMatd& X,
			   const MapMatd& Y,
			   const double threshold = 1e-7){
  SelfAdjointEigenSolver<MatrixXd> eig(AtA(X).selfadjointView<Lower>());
  MatrixXd   VDi(eig.eigenvectors() *
		 Dplus(eig.eigenvalues(), threshold).sqrt().matrix().asDiagonal());
  return Y - X * VDi * VDi.adjoint() * X.adjoint() * Y;
}

////////////////////////////////////////////////////////////////////////////
template<typename MatrixType, unsigned int UpLo>
inline double det_selfadjoint(const Eigen::SelfAdjointView< MatrixType, UpLo >& X,
                              bool logarithm = true){
  const VectorXd Dvec(X.ldlt().vectorD());
  double res = Dvec.array().log().sum();
  if (logarithm) {
    return res;
  } else {
    return exp(res);
  }
}

// It's the caller's responsibility to ensure X is self adjoint matrix
inline double det_selfadjoint(const MatrixXd& X, bool logarithm = true){
  return det_selfadjoint(X.selfadjointView<Lower>(),
			 logarithm);
}


//' Sample covariance matrix of residual for linear model
//' @inheritParams lm_resid_svd
//' @return The residual covariance matrix
// [[Rcpp::export]]
MatrixXd lm_resid_cov(const MapMatd& X,
                     const MapMatd& Y,
                     const double threshold = 1e-7){
  const int d(Y.cols());
  const Eigen::SelfAdjointEigenSolver<MatrixXd> eig(AtA(X));
  const ArrayXd Dp(Dplus(eig.eigenvalues(), threshold).sqrt());
  const MatrixXd VDp(eig.eigenvectors() * Dp.matrix().asDiagonal());
  const MatrixXd fitsqrt(VDp.adjoint() * X.adjoint() * Y);   // square root of Yhat'Yhat
  return MatrixXd(d, d).setZero().selfadjointView<Lower>()
    .rankUpdate(Y.adjoint())
    .rankUpdate(fitsqrt.adjoint(), -1.0);
}

//' Determinant of covariance matrix of residual for linear model
//' @inheritParams lm_resid_svd
//' @param logarithm A bool value indicating whether determinant or
//' log determinate should be returned.
//' @return The log determinant of residual covariance matrix
// [[Rcpp::export]]
double lm_resid_cov_det(const MapMatd& X,
		       const MapMatd& Y,
		       const bool logarithm = true,
		       const double threshold = 1e-7){
  MatrixXd ete = lm_resid_cov(X, Y, threshold);
  const VectorXd Dvec(ete.selfadjointView<Lower>().ldlt().vectorD());
  double res = Dvec.array().log().sum();
  if (logarithm) {
    return res;
  } else {
    return exp(res);
  }
}

//' Determinant of self cross product
//'
//' Calculate the determinate of X^T X. By using the symmetric property
//' of X^T X, this implementation is much faster than calling the
//' general function `determinant` on crossproduct(X)
//'
//' @param X matrix
//' @param logarithm A bool value indicating whether determinant or
//' log determinate should be returned.
//' @return The determinant or log determinant
// [[Rcpp::export]]
double det_AtA(const MapMatd& X, bool logarithm = true){
  const MatrixXd ata(AtA(X));
  return det_selfadjoint(ata, logarithm);
}
