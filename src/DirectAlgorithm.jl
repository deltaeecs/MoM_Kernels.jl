## 计算求解

"""
根据几何信息与基函数数量，计算阻抗矩阵
输入：
geosInfo::  几何信息，三角形、四面体、六面体的向量
nbf::       基函数数量
返回：
Zmat
"""
function getImpedanceMatrix(geosInfo::Vector{ST}, nbf::Integer) where {ST<:SurfaceCellType}
    if SimulationParams.ieT == :EFIE
        # 计算 RWG下 的 EFIE 阻抗矩阵
        Zmat =   impedancemat4EFIE4PEC(geosInfo, nbf, VSBFTypes.sbfType)
    elseif SimulationParams.ieT == :MFIE
        # 计算 RWG下 的 MFIE 阻抗矩阵
        Zmat =   impedancemat4MFIE4PEC(geosInfo, nbf, VSBFTypes.sbfType)
    elseif SimulationParams.ieT == :CFIE
        # 计算 RWG下 的 CFIE 阻抗矩阵
        Zmat =   impedancemat4CFIE4PEC(geosInfo, nbf, VSBFTypes.sbfType)
    end
    return Zmat
end

function getImpedanceMatrix(geosInfo::Vector{VT}, nbf::Integer) where {VT<:VolumeCellType}
    # 计算 SWG/PWC/RBF 下的 EFIE 阻抗矩阵
    Zmat =   impedancemat4VIE(geosInfo, nbf, VSBFTypes.vbfType)
    return Zmat
end

function getImpedanceMatrix(geosInfo::Vector{VT}, nbf::Integer) where {VT<:AbstractVector}
    if eltype(geosInfo[1]) <: SurfaceCellType
        # 计算 混合基函数下带有面网格的SWG/PWC/RBF 下的阻抗矩阵
        Zmat =   impedancemat4VSIE(geosInfo, nbf)
    else
        Zmat =   impedancemat4VIE(geosInfo, nbf, VSBFTypes.vbfType)
    end
    return Zmat
end


"""
根据几何信息与基函数数量，计算激励向量
输入：
geosInfo::  几何信息，三角形、四面体、六面体的向量
nbf::       基函数数量
source::    激励源
返回：
V::         激励向量
"""
function getExcitationVector(geosInfo::Vector{VST}, nbf, source) where {VST <: SurfaceCellType}
    if SimulationParams.ieT == :EFIE
        # 计算 RWG下 的 EFIE 激励向量
        V = excitationVectorEFIE(source, geosInfo, nbf)
    elseif SimulationParams.ieT == :MFIE
        # 计算 RWG下 的 MFIE 激励向量
        V = excitationVectorMFIE(source, geosInfo, nbf)
    elseif SimulationParams.ieT == :CFIE
        # 计算 RWG下 的 CFIE 激励向量
        V = excitationVectorCFIE(source, geosInfo, nbf)
    end
    return V
end

"""
根据几何信息与基函数数量，计算激励向量
输入：
geosInfo::  几何信息，三角形、四面体、六面体的向量
nbf::       基函数数量
source::    激励源
返回：
V::         激励向量
"""
function getExcitationVector(geosInfo::Vector{VST}, nbf, source) where {VST<:VolumeCellType}
    # 计算 体积分的 激励向量
    V = excitationVectorEFIE(source, geosInfo, nbf)
    return V
end

"""
根据几何信息与基函数数量，计算激励向量
输入：
geosInfo::  几何信息，三角形、四面体、六面体的向量
nbf::       基函数数量
source::    激励源
返回：
V::         激励向量
"""
function getExcitationVector(geosInfo::Vector{VT}, nbf, source) where {VT<:AbstractVector}
    if SimulationParams.ieT == :EFIE
        # 计算 RWG下 的 EFIE 激励向量
        V = excitationVectorEFIE(source, geosInfo, nbf)
    elseif SimulationParams.ieT == :MFIE
        # 计算 RWG下 的 MFIE 激励向量
        V = excitationVectorMFIE(source, geosInfo, nbf)
    elseif SimulationParams.ieT == :CFIE
        # 计算 RWG下 的 CFIE 激励向量
        V = excitationVectorCFIE(source, geosInfo, nbf)
    end
    return V
end


"""
根据几何信息与基函数数量，计算阻抗矩阵和激励向量
输入：
geosInfo::  几何信息，三角形、四面体、六面体的向量
nbf::       基函数数量
source::    激励源
返回：
Zmat::      阻抗矩阵
V::         激励向量
"""
function getImpedanceMatAndExciteV(geosInfo, nbf::Integer, source)
    Zmat    =   getImpedanceMatrix(geosInfo, nbf)
    V       =   getExcitationVector(geosInfo, nbf, source)
    Zmat, V
end

"""
根据几何信息与基函数数量，计算阻抗矩阵和激励向量
输入：
geosInfo::  几何信息，三角形、四面体、六面体的向量
bfsInfo::   基函数信息
source::    激励源
返回：
Zmat::      阻抗矩阵
V::         激励向量
"""
function getImpedanceMatAndExciteV(geosInfo, bfsInfo::Vector, source)
    nbf = length(bfsInfo)
    getImpedanceMatAndExciteV(geosInfo, nbf, source)
end
