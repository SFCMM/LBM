// SPDX-License-Identifier: BSD-3-Clause

#ifndef sfcmm_common_H
#define sfcmm_common_H
namespace sfcmm {
static constexpr int MAX_DIM = 4;
}

#include "common/boundingbox.h"
#include "common/compiler_config.h"
#include "common/constants.h"
#include "common/globalmpi.h"
#include "common/log.h"
#include "common/macros.h"
#include "common/randxor.h"
#include "common/random_special.h"
#include "common/sfcmm_types.h"
#include "common/timer.h"

#include "common/geometry/circle.h"

#include "common/util/backtrace.h"
#include "common/util/base64.h"
#include "common/util/binary.h"
#include "common/util/eigen.h"
#include "common/util/string_helper.h"
#include "common/util/sys.h"

#include "common/math/cartesian.h"
#include "common/math/hilbert.h"
#include "common/math/integration.h"
#include "common/math/mathfunctions.h"

#endif