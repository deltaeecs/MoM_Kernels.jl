
## 导入相关包
include("EFIEVSIERWGSWG.jl")
include("EFIEVSIERWGRBF.jl")
include("EFIEVSIERWGPWCHexa.jl")
include("EFIEVSIERWGPWCTetra.jl")

"""
计算VSIE的矩阵
"""
function impedancemat4VSIE(geosInfo::Vector{VT}, nbf::Integer) where {VT<:AbstractVector}
    # 计算 混合基函数 下的  RWG + SWG/PWC/RBF 下的阻抗矩阵
    Zmat = begin
        # RWG + SWG
        if (VSBFTypes.sbfType <: RWG) && (VSBFTypes.vbfType <: SWG)
            impedancemat4VSIERWGSWG(geosInfo, nbf)
        # RWG + RBF
        elseif (VSBFTypes.sbfType <: RWG) && (VSBFTypes.vbfType <: RBF)
            impedancemat4VSIERWGRBF(geosInfo, nbf)
        # RWG + PWC
        elseif (VSBFTypes.sbfType <: RWG) && (VSBFTypes.vbfType <: PWC)
            impedancemat4VSIERWGPWC(geosInfo, nbf)
        else
            thorw("检查体面积分的基函数类型，当前组合不支持！")
        end
    end

    Zmat
end