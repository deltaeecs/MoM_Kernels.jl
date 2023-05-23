using Primes

"""
    slicedim2bounds(sz::Int, nc::Int)

将区间 `1:sz` 划分为 `nc` 个区间并返回区间上下界。
*从 [MPIArray4MoMs](https://github.com/deltaeecs/MPIArray4MoMs.jl) 借的！
为的是避免提前引入 MPI 导致在集群上的 bug。因此该函数的修改必须与 
[MPIArray4MoMs](https://github.com/deltaeecs/MPIArray4MoMs.jl) 同步。*
"""
function slicedim2bounds(sz::Int, nc::Int)
    if sz >= nc
        chunk_size = div(sz,nc)
        remainder = rem(sz,nc)
        grid = zeros(Int64, nc+1)
        for i = 1:(nc+1)
            grid[i] += (i-1)*chunk_size + 1
            if i<= remainder
                grid[i] += i-1
            else
                grid[i] += remainder
            end
        end
        return grid
    else
        return [[1:(sz+1);]; zeros(Int, nc-sz)]
    end
end

"""
    slicedim2bounds(dims, nc::Int)

将区间 `dims` 划分为 `nc` 个区间并返回区间上下界。
*从 [MPIArray4MoMs](https://github.com/deltaeecs/MPIArray4MoMs.jl) 借的！
为的是避免提前引入 MPI 导致在集群上的 bug。因此该函数的修改必须与 
[MPIArray4MoMs](https://github.com/deltaeecs/MPIArray4MoMs.jl) 同步。*
"""
function slicedim2partition(dims, nc::Int)
    dims = [dims...]
    chunks = ones(Int, length(dims))
    f = sort!(collect(keys(factor(nc))), rev=true)
    k = 1
    while nc > 1
        # repeatedly allocate largest factor to largest dim
        if nc % f[k] != 0
            k += 1
            if k > length(f)
                break
            end
        end
        fac = f[k]
        (d, dno) = findmax(dims)
        # resolve ties to highest dim
        dno = findlast(isequal(d), dims)
        if dims[dno] >= fac
            dims[dno] = div(dims[dno], fac)
            chunks[dno] *= fac
        end
        nc = div(nc, fac)
    end
    return Tuple(chunks)
end


@doc """
    sizeChunks2cuts(Asize, chunks)
    sizeChunks2cuts(Asize::Int, chunks)
    sizeChunks2cuts(Asize, chunks::Int)
    sizeChunks2cuts(Asize::Int, chunks::Int)

将数组大小 `Asize` 按 `chunks` 进行分块。
*从 [MPIArray4MoMs](https://github.com/deltaeecs/MPIArray4MoMs.jl) 借的！
为的是避免提前引入 MPI 导致在集群上的 bug。因此该函数的修改必须与 
[MPIArray4MoMs](https://github.com/deltaeecs/MPIArray4MoMs.jl) 同步。*
"""
function sizeChunks2cuts(Asize, chunks)
    map(slicedim2bounds, Asize, chunks)
end
function sizeChunks2cuts(Asize::Int, chunks)
    map(slicedim2bounds, (Asize, ), chunks)
end
function sizeChunks2cuts(Asize, chunks::Int)
    map(slicedim2bounds, Asize, (chunks, ))
end
function sizeChunks2cuts(Asize::Int, chunks::Int)
    map(slicedim2bounds, (Asize, ), (chunks, ))
end

@doc """
    sizeChunksCuts2indices(Asize, nchunk, cuts::Tuple)
    sizeChunksCuts2indices(Asize, nchunk, cuts::Vector{I}) where{I<:Integer}

根据数组大小 `Asize` 分块数量 `nchunk` 以及各块索引区间 `cuts` 计算各块的索引。
*从 [MPIArray4MoMs](https://github.com/deltaeecs/MPIArray4MoMs.jl) 借的！
为的是避免提前引入 MPI 导致在集群上的 bug。因此该函数的修改必须与 
[MPIArray4MoMs](https://github.com/deltaeecs/MPIArray4MoMs.jl) 同步。*
"""
function sizeChunksCuts2indices(Asize, nchunk, cuts::Tuple)
    n = length(Asize)
    idxs = Array{NTuple{n,UnitRange{Int}}, n}(undef, nchunk...)
    for cidx in CartesianIndices(tuple(nchunk...))
        if n > 0
            idxs[cidx.I...] = ntuple(i -> (cuts[i][cidx[i]]:cuts[i][cidx[i] + 1] - 1), n)
        else
            throw("0 dim array not supported.")
        end
    end
    return idxs
end
function sizeChunksCuts2indices(Asize, nchunk, cuts::Vector{I}) where{I<:Integer}
    n = length(Asize)
    idxs = Array{NTuple{n,UnitRange{Int}}, n}(undef, nchunk...)
    for cidx in CartesianIndices(tuple(nchunk...))
        idxs[cidx.I...] = (cuts[cidx[1]]:cuts[cidx[1] + 1] - 1, )
    end
    return idxs
end

"""
    sizeChunks2idxs(Asize, nchunk)

Borrowed form DistributedArray.jl, get the slice of matrix
size Asize on each dimension with nchunk.
*从 [MPIArray4MoMs](https://github.com/deltaeecs/MPIArray4MoMs.jl) 借的！
为的是避免提前引入 MPI 导致在集群上的 bug。因此该函数的修改必须与 
[MPIArray4MoMs](https://github.com/deltaeecs/MPIArray4MoMs.jl) 同步。*
"""
function sizeChunks2idxs(Asize, nchunk)
    cuts = sizeChunks2cuts(Asize, nchunk)
    return sizeChunksCuts2indices(Asize, nchunk, cuts)
end
