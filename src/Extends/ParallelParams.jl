# 本文件用于扩展程序如 MPI 等的并行实现，不直接在本包使用

Base.@kwdef mutable struct ParallelParamsType
    nprocs::Int
end

"""
保存并行参数的实例。
"""
ParallelParams = ParallelParamsType(0)

"""
    set_nprocs!([;nprocs=1, np=nprocs])

设置并行核心数量为 `nprocs` 。
"""
function set_nprocs!(;nprocs=1, np=nprocs)
    ParallelParams.nprocs = np
    nothing
end

"""
    GeosIntervalType{T}

保存网格数据区间的类。
"""
Base.@kwdef mutable struct GeosIntervalType{T}
    tri::T
    tetra::T
    hexa::T
end

"""
保存网格数据区间的实例。
"""
GeosInterval = GeosIntervalType{UnitRange{Int}}(0:0, 0:0, 0:0)

"""
    set_geosInterval!(fn)

通过文件 `fn` 设置网格数据区间。
"""
function set_geosInterval!(fn)
    data = loadGeoInterval(fn)
    GeosInterval.tri = data.tri
    GeosInterval.tetra = data.tetra
    GeosInterval.hexa = data.hexa
    nothing
end

"""
    set_geosInterval!(fn)

载入文件 `fn` 读取网格数据区间。
"""
function loadGeoInterval(fn)
    load(fn, "data")
end
