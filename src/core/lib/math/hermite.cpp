//==================================================================================
// BSD 2-Clause License
//
// Copyright (c) 2025, CEA-List
//
// All rights reserved.
//
// Author TPOC: jules.dumezy@cea.fr
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//==================================================================================

#include "math/hermite.h"
#include "utils/exception.h"

#include <cmath>
#include <functional>
#include <vector>

namespace lbcrypto {

std::vector<std::complex<double>> hermiteInterpOrder1Coeff(std::function<double(double)> f, int p) {
  if (!p) {
    OPENFHE_THROW("The number of points of interest cannot be zero!");
  }
  std::vector<double> y(p);
  for (int i = 0; i < p; ++i) {
    y[i] = f(i);
  }

  std::vector<std::complex<double>> alpha(p, std::complex<double>(0.0, 0.0));

  double sum_y = 0;
  for (int i = 0; i < p; ++i) {
      sum_y += y[i];
  }
  alpha[0] = std::complex<double>(1.0 / p * sum_y, 0);

  for (int k = 1; k < p; ++k) {
    for (int l = 0; l < p; ++l) {
      std::complex<double> exp_term = std::exp(std::complex<double>(0, -2 * M_PI * k * l / p));
      alpha[k] += y[l] * exp_term;
    }

    alpha[k] *= 2.0 * (p - k) / (p * p);
  }

  return alpha;
}



std::vector<std::complex<double>> hermiteInterpOrder2Coeff(std::function<double(int)> f, int p) {
  if (!p) {
    OPENFHE_THROW("The number of points of interest cannot be zero!");
  }
  std::vector<double> y(p);
  for (int i = 0; i < p; ++i) {
    y[i] = f(i);
  }

  std::vector<std::complex<double>> alpha(p, std::complex<double>(0.0, 0.0));
  std::vector<std::complex<double>> beta(p/2, std::complex<double>(0.0, 0.0));
  std::vector<std::complex<double>> delta(p/2, std::complex<double>(0.0, 0.0));
  std::vector<std::complex<double>> theta(p/2, std::complex<double>(0.0, 0.0));

  double sum_y = 0;
  for (int i = 0; i < p; ++i) {
    sum_y += y[i];
  }
  alpha[0] = std::complex<double>(1.0 / p * sum_y, 0);

  for (int k = 1; k < p; ++k) {
    for (int l = 0; l < p; ++l) {
      std::complex<double> exp_term = std::exp(std::complex<double>(0, -2 * M_PI * k * l / p));
      alpha[k] += y[l] * exp_term;

      if (k > 0 && k <= p / 2) {
        beta[k - 1] += y[l] * exp_term;

        exp_term = std::exp(std::complex<double>(0, -2 * M_PI * (p + k) * l / p));
        delta[k - 1] += y[l] * exp_term;

        exp_term = std::exp(std::complex<double>(0, -2 * M_PI * (p - k) * l / p));
        theta[k - 1] += y[l] * exp_term;
      }
    }

    alpha[k] *= 2. * (p - k) / (p * p);

    if (k > 0 && k <= p / 2) {
      int gamma = (p % 2 == 0 && k == p / 2) ? 1. : 0.;
      double scale = ((2. - gamma) * k * (p - k)) / (p * p * p);
      beta[k - 1] *= scale;
      delta[k - 1] *= scale;
      theta[k - 1] *= scale;
    }
  }

  std::vector<std::complex<double>> coeffs(3*p/2, std::complex<double>(0.0, 0.0));
  coeffs[0] = alpha[0];
  for (int i = 1; i < 3*p/2; ++i) {
    if (i < p/2)
      coeffs[i] = alpha[i] + beta[i];
    else if (i == p/2)
      coeffs[i] = alpha[i] + beta[i] - (1 - p % 2) * 0.5 * theta[p - i - (p % 2)];
    else if (i < p)
      coeffs[i] = alpha[i] - 0.5 * theta[p - i]; 
    else if (i > p)
      coeffs[i] = - 0.5 * delta[i - p];
  }
  return coeffs;
}



std::vector<std::complex<double>> hermiteInterpOrder3Coeff(std::function<double(int)> f, int p) {
  if (!p) {
    OPENFHE_THROW("The number of points of interest cannot be zero!");
  }
  std::vector<double> y(p);
  for (int i = 0; i < p; ++i) {
    y[i] = f(i);
  }

  std::vector<std::complex<double>> alpha(p, std::complex<double>(0.0, 0.0));
  std::vector<std::complex<double>> beta(p, std::complex<double>(0.0, 0.0));
  std::vector<std::complex<double>> delta(p, std::complex<double>(0.0, 0.0));
  std::vector<std::complex<double>> theta(p, std::complex<double>(0.0, 0.0));

  double sum_y = 0;
  for (int i = 0; i < p; ++i) {
    sum_y += y[i];
  }
  alpha[0] = std::complex<double>(1.0 / p * sum_y, 0);

  for (int k = 1; k < p; ++k) {
    for (int l = 0; l < p; ++l) {
      std::complex<double> exp_term = std::exp(std::complex<double>(0, -2 * M_PI * k * l / p));
      alpha[k] += y[l] * exp_term;

      beta[k] += y[l] * exp_term;

      exp_term = std::exp(std::complex<double>(0, -2 * M_PI * (p + k) * l / p));
      delta[k] += y[l] * exp_term;

      exp_term = std::exp(std::complex<double>(0, -2 * M_PI * (p - k) * l / p));
      theta[k] += y[l] * exp_term;
    }

    alpha[k] *= 2. * (p - k) / (p * p);

    double scale = (2. * k * (p - k) * (2. * p - k)) / (3. * std::pow(p, 4));
    beta[k] *= scale;
    delta[k] *= scale;
    theta[k] *= scale;
  }

  std::vector<std::complex<double>> coeffs(2*p, std::complex<double>(0.0, 0.0));
  coeffs[0] = alpha[0];
  for (int i = 1; i < 2*p; ++i) {
    if (i < p)
      coeffs[i] = alpha[i] + beta[i] - 0.5 * theta[p - i];
    if (i > p)
      coeffs[i] = - 0.5 * delta[i - p];
  }
  return coeffs;
}

}
