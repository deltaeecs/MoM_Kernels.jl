include("AggOnBF/AggEFIE.jl")
include("AggOnBF/AggMFIE.jl")
include("AggOnBF/AggCFIE.jl")

"""
根据积分方程类型计算基层聚合项
"""
function getAggSBFOnLevel(level, geosInfo::Vector{ST}, 
        bfsInfo::Vector{BFT}) where {ST<:SurfaceCellType, BFT<:BasisFunctionType}
    if SimulationParams.ieT == :EFIE
        return aggSBFOnLevelEFIE(level, geosInfo, bfsInfo)
    elseif SimulationParams.ieT == :MFIE
        return aggSBFOnLevelMFIE(level, geosInfo, bfsInfo)
    elseif SimulationParams.ieT == :CFIE
        return aggSBFOnLevelCFIE(level, geosInfo, bfsInfo)
    end
end

"""
根据积分方程类型计算基层聚合项
"""
function getAggSBFOnLevel(level, geosInfo::Vector{VT}, 
        bfsInfo::Vector{BFT}) where {VT<:VolumeCellType, BFT<:BasisFunctionType}
    return aggSBFOnLevel(level, geosInfo, bfsInfo)
end


"""
根据积分方程类型计算基层聚合项
"""
function getAggSBFOnLevel(level, geosInfoV::Vector{VT1}, 
        bfsInfo::Vector{VT2}) where {VT1<:AbstractVector, VT2<:AbstractVector}
    FT  =   Precision.FT
    # 层采样点
    polesr̂sθsϕs =   level.poles.r̂sθsϕs
    # poles索引
    polesIndices    =   eachindex(polesr̂sθsϕs)
    # 采样多极子数量
    nPoles  =   polesIndices.stop
    # 预分配内存
    nbf     =   getNUnknown(bfsInfo)
    aggSBF  =   zeros(Complex{FT}, nPoles, 2, nbf)
    disaggSBF   =   zeros(Complex{FT}, nPoles, 2, nbf)

    if eltype(geosInfoV[1]) <: SurfaceCellType    
        # 面元、体元
        geoSs   =   geosInfoV[1]
        geoVs   =   geosInfoV[2]
        # 面基函数，体基函数
        bfSs    =   bfsInfo[1]
        bfVs    =   bfsInfo[2]
        # 面部分
        if SimulationParams.ieT == :EFIE
            aggSBFOnLevelEFIE!(aggSBF, disaggSBF, level, geoSs, eltype(bfSs))
        elseif SimulationParams.ieT == :MFIE
            aggSBFOnLevelMFIE!(aggSBF, disaggSBF, level, geoSs, eltype(bfSs))
        elseif SimulationParams.ieT == :CFIE
            aggSBFOnLevelCFIE!(aggSBF, disaggSBF, level, geoSs, eltype(bfSs))
        end
        # 体部分
        aggSBFOnLevel!(aggSBF, disaggSBF, level, geoVs, eltype(bfVs); setzero = false)
    else
        for i in 1:length(geosInfoV)
            geosInfo = geosInfoV[i]
            aggSBFOnLevel!(aggSBF, disaggSBF, level, geosInfo, eltype(bfsInfo[i]); setzero = false)
        end
    end

    return aggSBF, disaggSBF

end
