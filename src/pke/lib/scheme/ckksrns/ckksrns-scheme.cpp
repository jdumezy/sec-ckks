//==================================================================================
// BSD 2-Clause License
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
CKKS implementation. See https://eprint.iacr.org/2020/1118 for details.
 */

#include <memory>
#define PROFILE

#include "scheme/ckksrns/ckksrns-scheme.h"

namespace lbcrypto {

void SchemeCKKSRNS::Enable(PKESchemeFeature feature) {
    switch (feature) {
        case PKE:
            if (m_PKE == nullptr)
                m_PKE = std::make_shared<PKECKKSRNS>();
            break;
        case KEYSWITCH:
            // m_KeySwitch must be initialized later by calling SetKeySwitchingTechnique() with the value of key switching technique from cryptoparams
            break;
        case PRE:
            if (m_PRE == nullptr)
                m_PRE = std::make_shared<PRECKKSRNS>();
            break;
        case LEVELEDSHE:
            if (m_LeveledSHE == nullptr)
                m_LeveledSHE = std::make_shared<LeveledSHECKKSRNS>();
            break;
        case MULTIPARTY:
            if (m_Multiparty == nullptr)
                m_Multiparty = std::make_shared<MultipartyCKKSRNS>();
            break;
        case ADVANCEDSHE:
            if (m_AdvancedSHE == nullptr)
                m_AdvancedSHE = std::make_shared<AdvancedSHECKKSRNS>();
            break;
        case FHE:
            if (m_FHE == nullptr)
                m_FHE = std::make_shared<FHECKKSRNS>();
            break;
        case SCHEMESWITCH:
            if (m_SchemeSwitch == nullptr)
                m_SchemeSwitch = std::make_shared<SWITCHCKKSRNS>();
            break;
        case FBTS:
            m_Fbts = true;
            break;
        default:
            std::stringstream ss;
            ss << feature;
            OPENFHE_THROW(std::string("This feature [") + ss.str() + "] is not supported for CKKSRNS scheme");
    }
}

}  // namespace lbcrypto
