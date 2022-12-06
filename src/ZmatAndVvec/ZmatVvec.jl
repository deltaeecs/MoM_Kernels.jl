## 导入需要的阻抗矩阵

# 矩阵元的计算涉及的奇异积分处理函数
include("Singularity.jl")

# EFIE 基函数相关矩阵元素\激励向量计算相关函数
include("EFIE/EFIERWGTri.jl")
include("EFIE/EFIEPWCTetra.jl")
include("EFIE/EFIESWGTetra.jl")
include("EFIE/EFIEPWCHexa.jl")
include("EFIE/EFIERBFHexa.jl")
include("EFIE/EFIEExcitedVectors.jl")
include("EFIE/EFIEPWCTetraHexa.jl")
# 体表面积分方程
include("EFIE/EFIEVSIE.jl")

# MFIE 基函数相关矩阵元素\激励向量计算相关函数
include("MFIE/MFIERWGTri.jl")
include("MFIE/MFIEExcitedVectors.jl")

# CFIE 基函数相关矩阵元素\激励向量计算相关函数
include("CFIE/CFIERWGTri.jl")
include("CFIE/CFIEExcitedVectors.jl")


##