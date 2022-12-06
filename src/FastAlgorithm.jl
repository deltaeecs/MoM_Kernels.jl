## 本文件采用MLFMA进行快速计算
# 导入相关文件
include("MLFMA/MLFMA.jl")
include("MLFMA/Znear/Znear.jl")
include("MLFMA/Precondition/Precondition.jl")

"""
根据几何信息与基函数数量，计算阻抗矩阵算子和激励向量
输入：
geosInfo::  几何信息，三角形、四面体、六面体的向量
nbf::       基函数数量
source::    激励源
返回：
Zopt::      阻抗矩阵算子，由近场稀疏矩阵和远场八叉树聚合、转移、解聚组成
V::         激励向量
Octree::    八叉树
ZnearCSC::  阻抗矩阵近场元
"""
function getImpedanceOptAndExciteVOctree(geosInfo, bfsInfo, source)
    # 根据几何信息与基函数数量，计算阻抗矩阵算子
    Zopt    =   getImpedanceOpt(geosInfo, bfsInfo)
    # 计算激励向量
    V       =   getExcitationVector(geosInfo, nbf, source);
    # 返回
    Zopt, V
end

"""
根据几何信息与基函数数量，计算阻抗矩阵算子
输入：
geosInfo::  几何信息，三角形、四面体、六面体的向量
nbf::       基函数数量
source::    激励源
返回：
Zopt::      阻抗矩阵算子，由近场稀疏矩阵和远场八叉树聚合、转移、解聚组成
V::         激励向量
Octree::    八叉树
ZnearCSC::  阻抗矩阵近场元
"""
function getImpedanceOpt(geosInfo, bfsInfo)
    # 计算八叉树
    nLevels, octree     =   getOctreeAndReOrderBFs!(geosInfo, bfsInfo)
    # 叶层
    leafLevel   =   octree.levels[nLevels]
    # 计算近场矩阵CSC
    ZnearCSC    =   calZnearCSC(leafLevel, geosInfo, bfsInfo)
    # 构建矩阵向量乘积算子
    Zopt    =   MLMFAIterator(ZnearCSC, octree, geosInfo, bfsInfo)
    # 返回
    Zopt
end
