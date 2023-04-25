
"""
getGeoIDsInCubeChunk(cubes, ckunkIndice)

获取 ckunkIndice 内的所有 cube 的 geo ID ， 返回为 Tuple 形式以适应数组索引相关API

"""
function getGeoIDsInCubeChunk(cubes, chunkIndice::Tuple)

    geoIDs = reduce(vcat, cubes[i].geoIDs for i in chunkIndice[1])

    return (unique!(sort!(geoIDs)), )

end

"""
getGeoIDsInCubeChunk(cubes, ckunkIndice)

获取 ckunkIndice 内的所有 cube 的 geo ID ， 返回为 Tuple 形式以适应数组索引相关API

"""
function getGeoIDsInCubeChunk(cubes, chunkIndice::UnitRange)

    geoIDs = reduce(vcat, cubes[i].geoIDs for i in chunkIndice)

    return (unique!(sort!(geoIDs)), )

end

"""
getNeighborCubeIDs(cubes, chunkIndice)

获取 ckunkIndice 内的所有 cube 的 邻盒子ID， 返回为 Tuple 形式以适应数组索引相关API

TBW
"""
function getNeighborCubeIDs(cubes, chunkIndice::Tuple)

    neighborCubeIDs = reduce(vcat, cubes[i].neighbors for i in chunkIndice[1])

    return (unique!(sort!(neighborCubeIDs)), )

end

function getNeighborCubeIDs(cubes, chunkIndice::AbstractVector)

    neighborCubeIDs = reduce(vcat, cubes[i].neighbors for i in chunkIndice)

    return (unique!(sort!(neighborCubeIDs)), )

end


"""
saveGeosInfoChunks(geos::AbstractVector, cubes, name::AbstractString, nchunk::Int; dir = "")

将几何信息保存在
TBW
"""
function saveGeosInfoChunks(geos::AbstractVector, cubes, name::AbstractString, nchunk::Int; dir = "", cubes_ChunksIndices =   sizeChunks2idxs(length(cubes), nchunk))
    # 拿到各块的包含邻盒子的id
    cubesNeighbors_ChunksIndices    =   ThreadsX.mapi(chunkIndice -> getNeighborCubeIDs(cubes, chunkIndice), cubes_ChunksIndices)
    # 拿到包含邻盒子内的该块的所有几何信息 id
    geoInfo_chunks_indices  =   ThreadsX.mapi(chunkIndice -> getGeoIDsInCubeChunk(cubes, chunkIndice), cubesNeighbors_ChunksIndices)
    # 保存
    saveVec2Chunks(geos, name, geoInfo_chunks_indices; dir = dir, showpmeter = true, message = "GeosInfo")

    nothing

end

function getMeshDataSaveGeosInterval(filename; meshUnit=:mm, dir = "temp/GeosInfo");
    meshData, εᵣs   =  getMeshData(filename; meshUnit=meshUnit);
    saveGeoInterval(meshData; dir = dir)
    return meshData, εᵣs
end

function saveGeoInterval(meshData; dir = "temp/GeosInfo")
    !ispath(dir) && mkpath(dir)
    data = (tri = 1:meshData.trinum, tetra = (meshData.trinum + 1):(meshData.trinum + meshData.tetranum),
            hexa = (meshData.trinum + meshData.tetranum + 1):meshData.geonum,)
    jldsave(joinpath(dir, "geoInterval.jld2"), data = data)
    nothing
end


"""
    saveVec2Chunks(y::AbstractVector, name::AbstractString, nchunk::Int; dir = "")

    把向量分块保存。
TBW
"""
function saveVec2Chunks(y::AbstractVector, name::AbstractString, nchunk::Int; dir = "")

	indices = sizeChunks2idxs(length(y), nchunk)

	saveVec2Chunks(y, name, indices; dir = dir)

	nothing

end

"""
    saveVec2Chunks(y::AbstractVector, name::AbstractString, indices; dir = "")

    把向量分块保存。
TBW
"""
function saveVec2Chunks(y::AbstractVector, name::AbstractString, indices; dir = "", showpmeter = false, message = "")

	!ispath(dir) && mkpath(dir)

    pmeter = Progress(length(indices); desc = "Saving $message...", enabled = showpmeter)
	@floop for (i, indice) in enumerate(indices)
		jldsave(joinpath(dir, "$(name)_part_$i.jld2"), data = y[indice...], size = (length(y), ), indice = indice)
        next!(pmeter)
	end

	nothing

end


"""
    saveOctree(octree; dir="")

    保存八叉树。
TBW
"""
function saveOctree(octree; dir="")

    !ispath(dir) && mkpath(dir)

    # data = Dict{Symbol, Any}()
    ks = Symbol[]
    vals = []

    fieldsKeept = (:nLevels, :leafCubeEdgel, :bigCubeLowerCoor)

    @floop for k in fieldsKeept
        # data[k] = getfield(octree, k)
        push!(ks, k)
        push!(vals, getfield(octree, k))
    end

    nLevels = octree.nLevels
    levels = octree.levels
    kcubeIndices = nothing
    for iLevel in nLevels:-1:1
        @info "Saving level $(nLevels - iLevel + 1) / $nLevels."
        level = levels[iLevel]
        kcubeIndices = saveLevel(level; dir=dir, kcubeIndices = kcubeIndices)
    end

    # data[:levelsname] = joinpath(dir, "Level")
    push!(ks, :levelsname)
    push!(vals, joinpath(dir, "Level"))

    jldsave(joinpath(dir, "Octree.jld2"), data = (; zip(ks, vals)...))#data)

end

"""
    get_partition(nCubes, sizePoles, np)

    根据给定的盒子数、多极子数、进程数返回该层 pattern 的三个维度的划分数量。

"""
function get_partition(nCubes, sizePoles, np)

    # partition = if nCubes > 3np
    #     (1, 1, np)
    # else
    #     aggSize = (sizePoles, 2, nCubes)
    #     # 分区
    #     slicedim2mpi(aggSize, np)
    # end

    temp = slicedim2mpi((sizePoles, nCubes), np)
    partition = (temp[1], 1, temp[2])

    return partition

end


"""
    saveLevel(level, np = ParallelParams.nprocs; dir="", kcubeIndices = nothing)

    保存层信息。
TBW
"""
function saveLevel(level, np = ParallelParams.nprocs; dir="", kcubeIndices = nothing)

    !ispath(dir) && mkpath(dir)

    # cube要单独处理
    cubes = level.cubes
    level.cubes = eltype(cubes)[]

    # 多极子数
    sizePoles   =   length(level.poles.r̂sθsϕs)

    partition =   get_partition(length(cubes), sizePoles, np)

    indices = saveCubes(cubes, partition; name = "Level_$(level.ID)_Cubes", dir=dir, kcubeIndices = kcubeIndices)

    # 保存
    jldsave(joinpath(dir, "Level_$(level.ID).jld2"), data = level)

    level.cubes = cubes

    return indices

end

"""
这六个函数用于寻找盒子的子盒子区间内的比较函数，多重分派以实现
"""
func4Cube1stkInterval(cube::CubeInfo) = first(cube.kidsInterval)
func4Cube1stkInterval(i::T) where T <: Integer = i
func4Cube1stkInterval(interval::T) where T <: UnitRange = first(interval)
func4CubelastkInterval(cube::CubeInfo) = last(cube.kidsInterval)
func4CubelastkInterval(i::T) where T <: Integer = i
func4CubelastkInterval(interval::T) where T <: UnitRange = last(interval)

"""
根据 partition 计算在盒子方向本层所有 rank 到子层所有 rank 的 map。
"""
function get_partition_map(partition, kcubeIndices)
    collect(eachcol(reshape(repeat(1:length(kcubeIndices), inner = prod(partition) ÷ length(kcubeIndices)), :, partition[3])))
end

"""
    saveCubes(cubes, nchunk = ParallelParams.nprocs; name, dir="", kcubeIndices = nothing)

    保存盒子。
TBW
"""
function saveCubes(cubes, partition; name, dir="", kcubeIndices = nothing)
    nchunk = partition[3]
	!ispath(dir) && mkpath(dir)
    # 对盒子按 盒子数 和块数分块
	indices =   sizeChunks2idxs(length(cubes), nchunk)
    # 拿到各块的包含邻盒子的id
    cubesFarNeighbors_ChunksIndices    =   ThreadsX.mapi(chunkIndice -> getNeiFarNeighborCubeIDs(cubes, chunkIndice), indices)

    pmeter  =  Progress(length(indices), "Saving cubes...")

    # 重新组合 kcubeIndices，以避免层间分区不在同一维度时 子盒子区间计算错误
    kindices =  isnothing(kcubeIndices) ? nothing : get_partition_map(partition, kcubeIndices)
    @floop for (i, indice) in enumerate(indices)
        data = OffsetVector(cubes[indice...], indice...)
        idcs = indice[1]
        ghostindices::Vector{Int} = setdiff(cubesFarNeighbors_ChunksIndices[i][1], indice[1])
        # 子盒子与本层盒子区间错位时也会产生 ghost 数据需要保存在本地
        !isnothing(kcubeIndices) && begin
            tCubesInterval  =   last(searchsorted(cubes, first(kcubeIndices[first(kindices[i])]); by = func4Cube1stkInterval)):first(searchsorted(cubes, last(kcubeIndices[last(kindices[i])]); by = func4CubelastkInterval))
            otherGhostIdcs  =   setdiff(tCubesInterval, indice[1])
            # 将此部分 idcs 补充进来
            unique!(sort!(append!(ghostindices, otherGhostIdcs)))
        end

        ghostdata = sparsevec(ghostindices, cubes[ghostindices])
        cubes_i = PartitionedVector{eltype(cubes)}(length(cubes), data, idcs, ghostdata, ghostindices)
		jldsave(joinpath(dir, "$(name)_part_$i.jld2"), data = cubes_i)

        next!(pmeter)
	end

    return indices
    
end


"""
    getNeiFarNeighborCubeIDs(cubes, chunkIndice::Tuple)


    getFarNeighborCubeIDs(cubes, chunkIndice)

    获取 ckunkIndice 内的所有 cube 的 远亲盒子ID， 返回为 Tuple 形式以适应数组索引相关API

TBW
"""
function getNeiFarNeighborCubeIDs(cubes, chunkIndice::Tuple)

    neighborCubeIDs     = reduce(vcat, cubes[i].neighbors for i in chunkIndice[1])
    farneighborCubeIDs  = reduce(vcat, cubes[i].farneighbors for i in chunkIndice[1])

    return (unique!(sort!(vcat(neighborCubeIDs, farneighborCubeIDs))), )

end

