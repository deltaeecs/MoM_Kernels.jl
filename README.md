# MoM_Kernels


[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://deltaeecs.github.io/MoM_Kernels.jl/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://deltaeecs.github.io/MoM_Kernels.jl/dev/)
[![Build Status](https://github.com/deltaeecs/MoM_Kernels.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/deltaeecs/MoM_Kernels.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/deltaeecs/MoM_Kernels.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/deltaeecs/MoM_Kernels.jl)

[![CEM_MoMs](https://img.shields.io/badge/CEM_MoMs-github-orange.svg)](https://github.com/deltaeecs/CEM_MoMs.git)
[![CEM_MoMs](https://img.shields.io/badge/CEM_MoMs-gitee-orange.svg)](https://gitee.com/deltaeecs/CEM_MoMs.git)

#### 介绍

提供 CEM_MoMs (主页：[github](https://github.com/deltaeecs/CEM_MoMs.git), [gitee](https://gitee.com/deltaeecs/CEM_MoMs.git)) 包的核心算法，包括各类基函数的矩阵元计算、多层快速多极子（MLFMA）各种类型和函数、求解器、后处理（RCS、远场计算等功能）。但均不包含绘图函数。 CEM_MoMs 本身被拆分为几个独立的包以方便开发时快读编译，同时避免在无图形化界面使用时调入绘图相关包而导致报错。