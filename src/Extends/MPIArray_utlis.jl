using Primes

"""
slicedim2mpi(sz::Int, nc::Int)

从 MPIArray4MoMs借的！为的是避免提前引入 MPI 导致在集群上的 bug。
TBW
"""
function slicedim2mpi(sz::Int, nc::Int)
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
slicedim2mpi(dims, nc::Int)

从 MPIArray4MoMs借的！为的是避免提前引入 MPI 导致在集群上的 bug。
TBW
"""
function slicedim2mpi(dims, nc::Int)
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



function sizeChunks2cuts(Asize, chunks)
    map(slicedim2mpi, Asize, chunks)
end

function sizeChunks2cuts(Asize, chunks::Int)
    map(slicedim2mpi, Asize, (chunks, ))
end

function sizeChunks2cuts(Asize::Int, chunks)
    map(slicedim2mpi, (Asize, ), chunks)
end

function sizeChunks2cuts(Asize::Int, chunks::Int)
    map(slicedim2mpi, (Asize, ), (chunks, ))
end

function sizeChunksCuts2indices(Asize, chunks, cuts::Tuple)
    n = length(Asize)
    idxs = Array{NTuple{n,UnitRange{Int}}, n}(undef, chunks...)
    for cidx in CartesianIndices(tuple(chunks...))
        if n > 0
            idxs[cidx.I...] = ntuple(i -> (cuts[i][cidx[i]]:cuts[i][cidx[i] + 1] - 1), n)
        else
            throw("0 dim array not supported.")
        end
    end

    return idxs
end

function sizeChunksCuts2indices(Asize, chunks, cuts::Vector{I}) where{I<:Integer}
    n = length(Asize)
    idxs = Array{NTuple{n,UnitRange{Int}}, n}(undef, chunks...)
    for cidx in CartesianIndices(tuple(chunks...))
        idxs[cidx.I...] = (cuts[cidx[1]]:cuts[cidx[1] + 1] - 1, )
    end

    return idxs
end


"""
    sizeChunks2idxs(Asize, chunks)

    Borrowed form DistributedArray.jl, get the slice of matrix
    size Asize on each dimension with chunks.

TBW
"""
function sizeChunks2idxs(Asize, chunks)
    cuts = sizeChunks2cuts(Asize, chunks)
    return sizeChunksCuts2indices(Asize, chunks, cuts)
end
