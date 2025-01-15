//==================================================================================
// BSD 2-Clause License
//
// This file has been modified from the original version.
// Changes made by Jules Dumezy at CEA-List in 2025.
//
// Copyright (c) 2025, CEA-List
//
// Author TPOC: jules.dumezy@cea.fr
//
// Copyright (c) 2014-2022, NJIT, Duality Technologies Inc. and other contributors
//
// All rights reserved.
//
// Author TPOC: contact@openfhe.org
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

/*
  This code provides Chebyshev approximation utilities
 */

#include "math/chebyshev.h"
#include "utils/exception.h"

#include <cmath>
#include <cstdint>
#include <functional>
#include <vector>

namespace lbcrypto {

std::vector<double> EvalChebyshevCoefficients(std::function<double(double)> func, double a, double b, uint32_t degree) {
    if (!degree) {
        OPENFHE_THROW("The degree of approximation can not be zero");
    }
    // the number of coefficients to be generated should be degree+1 as zero is also included
    size_t coeffTotal{degree + 1};
    double bMinusA = 0.5 * (b - a);
    double bPlusA  = 0.5 * (b + a);
    double PiByDeg = M_PI / static_cast<double>(coeffTotal);
    std::vector<double> functionPoints(coeffTotal);
    for (size_t i = 0; i < coeffTotal; ++i)
        functionPoints[i] = func(std::cos(PiByDeg * (i + 0.5)) * bMinusA + bPlusA);

    double multFactor = 2.0 / static_cast<double>(coeffTotal);
    std::vector<double> coefficients(coeffTotal);
    for (size_t i = 0; i < coeffTotal; ++i) {
        for (size_t j = 0; j < coeffTotal; ++j)
            coefficients[i] += functionPoints[j] * std::cos(PiByDeg * i * (j + 0.5));
        coefficients[i] *= multFactor;
    }

    return coefficients;
}

std::vector<std::complex<double>> EvalChebyshevCoefficients(std::function<std::complex<double>(double)> func, double a, double b, uint32_t degree) {
    if (!degree) {
        OPENFHE_THROW("The degree of approximation can not be zero");
    }

    size_t coeffTotal = degree + 1;
    double bMinusA = 0.5 * (b - a);
    double bPlusA = 0.5 * (b + a);
    double PiByDeg = M_PI / static_cast<double>(coeffTotal);
    std::vector<std::complex<double>> functionPoints(coeffTotal);
    for (size_t i = 0; i < coeffTotal; ++i) {
        functionPoints[i] = func(std::cos(PiByDeg * (i + 0.5)) * bMinusA + bPlusA);
    }

    double multFactor = 2.0 / static_cast<double>(coeffTotal);
    std::vector<std::complex<double>> coefficients(coeffTotal, std::complex<double>(0.0, 0.0));
    for (size_t i = 0; i < coeffTotal; ++i) {
        for (size_t j = 0; j < coeffTotal; ++j) {
            coefficients[i] += functionPoints[j] * std::cos(PiByDeg * i * (j + 0.5));
        }
        coefficients[i] *= multFactor;
    }

    return coefficients;
}

}  // namespace lbcrypto
