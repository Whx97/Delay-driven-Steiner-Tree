#pragma once

#include <boost/functional/hash.hpp>
#include <unordered_map>

#include "salt/base/eval.h"

namespace salt {

class CriticalBuilder {
public:
    void Run(const Net& net, Tree& tree);
};

}  // namespace salt