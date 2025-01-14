#pragma once

#include "bonn.h"
#include "salt/base/eval.h"
#include "salt/base/tree.h"
#include "salt/utils/enum.h"
// #include "myself/self.h"
#include "delay_tree/DelayTreeBuilder.h"
// Note: need to set salt::BonnBuilder::unitCap

BETTER_ENUM(Method, int, FLUTE, RSA, ESRSA, MST, BRBC, KRYS, KRY, PD, ES, BONN, SALT_R0, SALT_R1, SALT_R2, SALT_R3, DelayTree_R1, DelayTree_R2, Critical);

void GetATree(const salt::Net& net, salt::Tree& tree, Method type, double eps, bool checkTree = true);