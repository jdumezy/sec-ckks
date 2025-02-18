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

#define PROFILE

#include "openfhe.h"

using namespace lbcrypto;

void FuncBootstrapExample(std::function<double(double)> f, int bits, int p);

int main(int argc, char* argv[]) {
    int bits = 4;
    // p will determine the plaintext space, which is [0..p-1]
    int p = pow(2, bits);

    // f is the function that is evaluated during bootstrapping
    auto f = [](double x) -> double { return (double) ((int)x % 2); };
    /*auto f = [](double x) -> double { return x * x; };*/

    FuncBootstrapExample(f, bits, p);
}

void FuncBootstrapExample(std::function<double(double)> f, int bits, int p) {
    CCParams<CryptoContextCKKSRNS> parameters;

    SecretKeyDist secretKeyDist = SPARSE_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);

    // For secure computations, use parameters.SetSecurityLevel(HEStd_128_classic);
    // and remove the two lines below.
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetRingDim(1 << 12);

    parameters.SetNumLargeDigits(3);
    parameters.SetKeySwitchTechnique(HYBRID);

#if NATIVEINT == 128
    OPENFHE_THROW("NATIVEINT 128 is not supported, use 64.");
#endif
    // Currently, only FIXEDMANUAL mode is supported for functional bootstrapping.
    ScalingTechnique rescaleTech = FIXEDMANUAL;
    usint dcrtBits               = 48;
    usint firstMod               = 49;


    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);

    std::vector<uint32_t> levelBudget = {4, 2};

    std::vector<uint32_t> bsgsDim = {0, 0};

    /*usint depth = levelBudget[0] + levelBudget[1] + 12 + std::log2(p);*/
    usint depth = levelBudget[0] + levelBudget[1] + 1 + 5 + 4 + bits + 2;
    parameters.SetMultiplicativeDepth(depth);

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(FHE);
    
    cc->Enable(FBTS);

    usint ringDim = cc->GetRingDimension();
    std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl << std::endl;

    // Precomputations for bootstrapping
    int numSlots = ringDim/2;
    cc->EvalFuncBootstrapSetup(levelBudget, bsgsDim, numSlots, bits);

    // Key Generation
    auto keyPair = cc->KeyGen();
    cc->EvalMultKeyGen(keyPair.secretKey);
    cc->EvalBootstrapKeyGen(keyPair.secretKey, numSlots);

    // Input vector [0, 1, ..., p-2, p-1]
    std::vector<std::complex<double>> x;
    for (int i = 0; i < p; ++i) {
        x.push_back((std::complex<double>(i, p - i - 1)));
    }

    Plaintext ptxt =  cc->MakeCKKSPackedPlaintext(x, 1, depth - levelBudget[1], nullptr, numSlots);
    ptxt->SetLength(numSlots);
    std::cout << "Input = " << ptxt << std::endl;

    Ciphertext<DCRTPoly> ctxt = cc->Encrypt(keyPair.publicKey, ptxt);

    // Order of the hermite interpolation, should be between 1 and 3
    int hermite_order = 1;

    // We evaluate f(ctxt) while bootstrapping
    auto ctxtResult = cc->EvalFuncBootstrap(ctxt, f, p, hermite_order);

    Plaintext result;
    cc->Decrypt(keyPair.secretKey, ctxtResult, &result);
    result->SetLength(p);
    // Output vector [f(0), f(1), ..., f(p-2), f(p-1)]
    std::cout << "Output = " << result << std::endl;
}
