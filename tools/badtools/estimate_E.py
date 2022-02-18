#!/usr/bin/env python
# coding: utf-8

import math

#kBT前の比例定数
const = 0.1

#kB(eV/K) : ボルツマン定数
kb = 8.617333 * (10**(-5))

#T(K) : 絶対温度
T = 600

#q(無次元) : イオンの価数
q = 1

#d(Å) : ジャンプの平均距離
d = 1.5

#E(eV/Å) : かける電場の強さ
E = (const * kb * T) / (q * d)
logE = math.log10(E)

print("logE = ", logE)
