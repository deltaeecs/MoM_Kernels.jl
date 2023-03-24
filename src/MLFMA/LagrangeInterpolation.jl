
"""
计算分别采用高斯求积、中点求积计算 θ,ϕ 方向的采样点的坐标、权重
lb::FT，  积分区域下界
hb::FT,   积分区域上界
Nsample::IT, 采样点数
mod::Symbol， 模式，接受 :uni, 均值积分 :glq, 高斯-勒让德积分 两种模式
"""
function integral1DXW(lb::FT, hb::FT, Nsample::IT, mod::Symbol) where{IT<:Integer, FT<:Real}
    
    X, W = zeros(FT, Nsample), zeros(FT, Nsample)

    if mod == :uni # 均值积分
        # 采样间隔
        dl = (hb-lb)/Nsample
        # 采样点、权重计算
        for j in 1:Nsample
            X[j] = lb + (j-1)*dl + dl/2
            W[j] = dl
        end
    elseif mod == :glq # 高斯-勒让德积分
        Dx      =   0.5 * (hb - lb)
        center  =   0.5 * (hb + lb)
        XGL, WGL    =   gausslegendre(Nsample)
        X   .=   center .+ Dx .* XGL
        W   .=   abs(Dx) .* WGL
    else
        throw("仅接受 :uni, 均值积分 :glq, 高斯-勒让德积分 两种模式")
    end
    return X, W
    
end


"""
计算各层八叉树求积坐标、求积权重
lb::FT，  积分区域下界
hb::FT,   积分区域上界
nlevels::IT, 八叉树叶层ID
mod::Symbol， 模式，接受 :uni, 均值积分(ϕ方向) :glq, 高斯-勒让德积分(θ方向) 两种模式
"""
function octreeXWNCal(lb::FT, hb::FT, L::IT, mod::Symbol) where{IT<:Integer, FT<:Real}

    N = if mod == :uni  # 均值积分
        2*(L + 1)
    elseif mod == :glq  # 高斯-勒让德积分
        L + 1
    else
        throw("仅接受 :uni, 均值积分 :glq, 高斯-勒让德积分 两种模式")
    end

    Xs, Ws    =   integral1DXW(lb, hb, N, mod)

    return Xs, Ws
end

"""
    gq_xsws_on_sphere(L)

    计算单位球面 2(L+1) 阶高斯求积的采样点坐标权重

TBW
"""
function gq_xsws_on_sphere(L; FT = Precision.FT)

    ## 积分点和求积权重数据
    # θ方向
    Xcosθs, Wθs   =   octreeXWNCal(one(FT), -one(FT), L, :glq)
    # 将θ方向高斯-勒让德求积坐标从 [1.,-1.] 转换到 [0,π]
    Xθs      =   acos.(Xcosθs)
    # ϕ方向
    Xϕs, Wϕs = octreeXWNCal(zero(FT), convert(FT, 2π), L, :uni)
    
    # 所有采样点直角坐标
    nodes = reduce(hcat, [sphere2cart(1, θ, ϕ) for ϕ in Xϕs for θ in Xθs ])

    # 所有采样点权重
    ws =  [Wθ * Wϕ for Wϕ in Wϕs for Wθ in Wθs]

    return nodes, ws

end

"""
多极子的极信息，即角谱空间采样信息
Xθs::Vector{FT}，θ方向的采样点坐标（rad单位），高斯-勒让德求积
Xϕs::Vector{FT}，ϕ方向的采样点坐标（rad单位），均值求积
Wθϕs::Vector{FT}，采样点权重，用于积分时使用，在MLFMA中直接乘在转移项
"""
struct GLPolesInfo{FT<:Real} <:PolesInfo{FT}
    Xθs::Vector{FT}
    Xϕs::Vector{FT}
    Wθϕs::Vector{FT}
    r̂sθsϕs::Vector{r̂θϕInfo{FT}}
end

"""
计算八叉树的积分相关信息，包括截断项、各层积分点和求积权重数据
输入:
levelCubeEdgel::FT,  层盒子边长, 一般叶层为0.25λ，其中 λ 为区域局部波长。
返回值
L           ::IT， 层 截断项
levelsPoles ::Vector{GLPolesInfo{FT}}，从叶层到第 “2” 层的角谱空间采样信息
"""
function levelIntegralInfoCal(levelCubeEdgel::FT, ::Union{Val{:Lagrange2Step}, Val{:Lagrange1Step}}) where{FT<:Real}
    ## 计算截断项
    L = truncationLCal(levelCubeEdgel)
    
    ## 积分点和求积权重数据
    # θ方向
    Xcosθs, Wθs   =   octreeXWNCal(one(FT), -one(FT), L, :glq)
    # 将θ方向高斯-勒让德求积坐标从 [1.,-1.] 转换到 [0,π]
    Xθs      =   acos.(Xcosθs)
    # ϕ方向
    Xϕs, Wϕs = octreeXWNCal(zero(FT), convert(FT, 2π), L, :uni)
    
    # 将数据保存在 levelsPoles 中，按照 θ 方向连续的顺序，将所有采样点信息保存为一向量
    # 计算所有极子的信息
    r̂sθsϕs =  [r̂θϕInfo{FT}(θ, ϕ) for ϕ in Xϕs for θ in Xθs]
    # 所有采样点权重
    Wθϕs   =  [Wθ * Wϕ for Wϕ in Wϕs for Wθ in Wθs]

    # 创建Poles实例保存
    Poles  =   GLPolesInfo{FT}(Xθs, Xϕs, Wθϕs, r̂sθsϕs)
    
    return L, Poles
end

"""
保存 θ, ϕ 两个方向的稀疏插值矩阵，
θ方向为 (npXθs, ntXθs) 稀疏矩阵, 用于左乘本层多极子矩阵，在 θ 方向插值
ϕ方向为 (ntXϕs, ntXϕs) 稀疏矩阵, 用于右乘本层多极子矩阵，在 ϕ 方向插值
θCSCT   ::SparseMatrixCSC{FT} 稀疏矩阵, θ 方向插值矩阵的转置，用于左乘本层多极子矩阵，在 θ 方向反插值
ϕCSCT   ::SparseMatrixCSC{FT} 稀疏矩阵, ϕ 方向插值矩阵的转置，用于左乘本层多极子矩阵，在 ϕ 方向反插值
"""
mutable struct LagrangeInterpInfo{IT<:Integer, FT<:Real} <: InterpInfo{IT, FT}
    θCSC    ::SparseMatrixCSC{FT, IT}
    ϕCSC    ::SparseMatrixCSC{FT, IT}
    θCSCT   ::SparseMatrixCSC{FT, IT}
    ϕCSCT   ::SparseMatrixCSC{FT, IT}
    LagrangeInterpInfo{IT, FT}() where {IT<:Integer, FT<:Real}  =   new{IT, FT}()
    LagrangeInterpInfo{IT, FT}(θCSC, ϕCSC, θCSCT, ϕCSCT) where {IT<:Integer, FT<:Real}  =   new{IT, FT}(θCSC, ϕCSC, θCSCT, ϕCSCT)
end

"""
带参数的构造函数
"""
function LagrangeInterpInfo(θCSC::SparseMatrixCSC{FT, IT}, ϕCSC::SparseMatrixCSC{FT, IT}) where {FT<:Real, IT}

    θCSC     =   θCSC
    ϕCSC     =   ϕCSC
    θCSCT    =   sparse(transpose(θCSC))
    ϕCSCT    =   sparse(transpose(ϕCSC))

    return LagrangeInterpInfo{IT, FT}(θCSC, ϕCSC, θCSCT, ϕCSCT)

end


"""
拉格朗日分步插值
"""
function interpolate(weights::LagrangeInterpInfo{IT, FT}, data::AbstractArray) where {IT, FT}
    target  =   zeros(eltype(data), size(weights.θCSC, 1), size(data, 2))
    interpolate!(target, weights, data)
end
"""
拉格朗日分步反插值
"""
function anterpolate(weights::LagrangeInterpInfo{IT, FT}, data::AbstractArray) where {IT, FT}
    target  =   zeros(eltype(data), size(weights.ϕCSCT, 1), size(data, 2))
    anterpolate!(target, weights, data)
end

"""
拉格朗日分步插值
"""
function interpolate!(target::AbstractArray, weights::LagrangeInterpInfo{IT, FT}, data::AbstractArray) where {IT, FT}
    # target .= weights.θCSC * (weights.ϕCSC * data)
    mul!(target, weights.θCSC, weights.ϕCSC * data)
end
"""
拉格朗日分步反插值
"""
function anterpolate!(target::AbstractArray, weights::LagrangeInterpInfo{IT, FT}, data::AbstractArray) where {IT, FT}
    #target .= weights.ϕCSCT * (weights.θCSCT * data)
    mul!(target, weights.ϕCSCT, weights.θCSCT * data)
end


"""
保存总的 稀疏插值矩阵，用于单步插值，根据稀疏度决定保存稀疏阵或是稠密阵
θϕCSC       ::AbstractMatrix{FT} 稀疏矩阵, θ 方向插值矩阵的转置，用于左乘本层多极子矩阵，在 θ 方向反插值
θϕCSCT      ::AbstractMatrix{FT} 稀疏矩阵, ϕ 方向插值矩阵的转置，用于左乘本层多极子矩阵，在 ϕ 方向反插值
"""
mutable struct LagrangeInterp1StepInfo{IT, FT<:Real} <: InterpInfo{IT, FT}
    θϕCSC    ::SparseMatrixCSC{FT, IT}
    θϕCSCT   ::SparseMatrixCSC{FT, IT}
    LagrangeInterp1StepInfo{IT, FT}() where {IT, FT<:Real}  =   new{IT, FT}()
    LagrangeInterp1StepInfo{IT, FT}(θϕCSC, θϕCSCT) where {IT, FT<:Real}  =   new{IT, FT}(θϕCSC, θϕCSCT)
end


"""
带参数的构造函数
"""
function LagrangeInterp1StepInfo(θCSC::SparseMatrixCSC{FT, IT}, ϕCSC::SparseMatrixCSC{FT, IT}) where {FT<:Real, IT}

    θϕCSC    =   θCSC*ϕCSC
    θϕCSCT   =   sparse(transpose(θϕCSC))

    return LagrangeInterp1StepInfo{IT, FT}(θϕCSC, θϕCSCT)

end

function Base.convert(::Type{LagrangeInterp1StepInfo}, interp::LagrangeInterpInfo{IT, FT}) where {IT, FT}
    LagrangeInterp1StepInfo(interp.θCSC, interp.ϕCSC)
end

"""
拉格朗日单步插值
"""
function interpolate(weights::LagrangeInterp1StepInfo{IT, FT}, data::AbstractArray) where {IT, FT}
    target  =   zeros(eltype(data), size(weights.θϕCSC, 1), size(data, 2))
    interpolate!(target, weights, data)
end

"""
拉格朗日单步反插值
"""
function anterpolate(weights::LagrangeInterp1StepInfo{IT, FT}, data::AbstractArray) where {IT, FT}
    target  =   zeros(eltype(data), size(weights.θϕCSCT, 1), size(data, 2))
    anterpolate!(target, weights, data)
end

"""
拉格朗日单步插值
"""
function interpolate!(target::AbstractArray, weights::LagrangeInterp1StepInfo{IT, FT}, data::AbstractArray) where {IT, FT}
    #target .= weights.θϕCSC * data
    mul!(target, weights.θϕCSC, data)
end

"""
拉格朗日单步反插值
"""
@inline function anterpolate!(target::AbstractArray, weights::LagrangeInterp1StepInfo{IT, FT}, data::AbstractArray) where {IT, FT}
    #target .= weights.θϕCSCT * data
    mul!(target, weights.θϕCSCT, data)
end

"""
采用局部插值、两步插值法，计算局部坐标到全局坐标的稀疏插值矩阵，Julia数据存储为列主的，因此使用 压缩稀疏列(Compressed Sparse Column, CSC)
pLevelPoles::GLPolesInfo{FT}， 父层多极子信息
tLevelPoles::GLPolesInfo{FT}， 本层多极子信息
"""
function interpolationCSCMatCal(pLevelPoles::GLPolesInfo{FT}, tLevelPoles::GLPolesInfo{FT}, nlocalInterp::IT) where {IT<:Integer, FT<:Real}
    
    # 父层本层多极子数据
    # 父层本层的 θ 坐标
    pXθs    =   pLevelPoles.Xθs
    tXθs    =   tLevelPoles.Xθs
    # θ方向的父层本层插值点数
    npXθs   =   length(pXθs)
    ntXθs   =   length(tXθs)
    # 父层本层的ϕ坐标
    pXϕs    =   pLevelPoles.Xϕs
    tXϕs    =   tLevelPoles.Xϕs
    # ϕ方向的父层本层插值点数
    npXϕs   =   length(pXϕs)
    ntXϕs   =   length(tXϕs)


    ################################################################
    # 先计算 θ 方向

    # 判断插值点是否过多
    nlocalInterp > npXθs && throw("插值点数设置过多！(一般是频率与网格不匹配，输入频率太低导致。)")
    # θ方向父层坐标在本层坐标的区间位置pInt
    pθsIntθs    =   cooraInCoorb(pXθs, tXθs)

    # 初始化权重矩阵
    interWθs    =   ones(FT, (nlocalInterp, npXθs))
    interIDθs   =   zeros(IT, (nlocalInterp, npXθs))
    # 插值点的相对坐标（[-nlocalInterp/2 : nlocalInterp/2 ]）
    RelativeOffsets =   (1:nlocalInterp) .- nlocalInterp ÷ 2
    # 对目标插值点循环计算其插值点
    @inbounds for ipXθs in 1:npXθs
        pθIntθ  =   pθsIntθs[ipXθs]
        # 插值点坐标
        iInterIDθs   =   pθIntθ .+  RelativeOffsets
        # 当插值点越过边界（超出 θ 的取值范围），要调整对应的 θ 值到负值或超出 π 的值
        ###### 采用越界插值，不在此处调整插值点
        # if iInterIDθs[1] < 1
        #     iInterIDθs   =  iInterIDθs .- (iInterIDθs[1] - 1)
        # elseif iInterIDθs[end] > ntXθs
        #     iInterIDθs   =  iInterIDθs .- (iInterIDθs[end] - ntXθs)
        # end
        # # 写入数据
        interIDθs[:,ipXθs] .=  iInterIDθs
    end

    # 利用向量化直接计算所有插值权重，此处采用高斯插值
    @inbounds for i in 1:nlocalInterp
        # 所有待插值点处的第i个插值点
        θsInterp    =   [pickθ(idx, tXθs) for idx in interIDθs[i, :]]
        sinHalfDiffθpLevel  =   sin.((pXθs .- θsInterp) ./ 2)
        for j in 1:nlocalInterp
            if i != j
                θsInterpLocal       =   [pickθ(idx, tXθs) for idx in interIDθs[j, :]]
                sinHalfDiffθtLevel  =   sin.((θsInterpLocal .- θsInterp) ./ 2)
                interWθs[j, :]     .*=   sinHalfDiffθpLevel ./ sinHalfDiffθtLevel
            end #if
        end #for
    end #for

    # 计算误差导致权重之和不为1，此处进行修正
    @inbounds for i in 1:npXθs
        interWθs[:, i] ./=  sum(interWθs[:, i])
    end #for

    ################################################################
    # 计算ϕ方向

    # 判断插值点是否过多
    nlocalInterp > npXϕs && throw("插值点数设置过多！")
    # ϕ方向父层坐标在本层坐标的区间位置pInt
    pϕsIntϕs    =   cooraInCoorb(pXϕs, tXϕs)

    # 初始化权重矩阵
    interWϕs    =   ones(FT, (nlocalInterp, npXϕs))
    interIDϕs   =   zeros(IT, (nlocalInterp, npXϕs))
    # 插值点的相对坐标（[-nlocalInterp/2 : nlocalInterp/2 ]）
    RelativeOffsets =   (1:nlocalInterp) .- nlocalInterp ÷ 2
    # 对目标插值点循环计算其插值点
    @inbounds for ipXϕs in 1:npXϕs
        pϕIntϕ  =   pϕsIntϕs[ipXϕs]
        # 插值点坐标
        iInterIDϕs   =   collect(pϕIntϕ .+  RelativeOffsets)
        # 当插值点在边界附近时要修正，与θ方向不同， ϕ方向的具有周期性，因此直接在索引角度值时进行修改
        # 写入数据
        interIDϕs[:,ipXϕs] .=  iInterIDϕs
    end #ipXϕs

    # 利用向量化直接计算所有插值权重，此处采用高斯插值
    @inbounds for i in 1:nlocalInterp
        # 所有待插值点处的第i个插值点
        ϕsInterp    =   [pickϕ(idx, tXϕs) for idx in interIDϕs[i,:]]
        sinHalfDiffϕpLevel  =   sin.((pXϕs .- ϕsInterp) ./ 2)
        for j in 1:nlocalInterp
            if i != j
                ϕsInterpLocal       =   [pickϕ(idx, tXϕs) for idx in interIDϕs[j,:]]
                sinHalfDiffϕtLevel  =   sin.((ϕsInterpLocal .- ϕsInterp) ./ 2)
                interWϕs[j, :]     .*=   sinHalfDiffϕpLevel ./ sinHalfDiffϕtLevel
            end #if
        end #for
    end #for

    # 计算误差导致权重之和不为1，此处进行修正
    @inbounds for i in 1:npXϕs
        interWϕs[:, i] ./=  sum(view(interWϕs, :, i))
    end #for


    ####### 放弃以下部分，原用于两侧插值，现用单侧插值，即将所有采样点排列成一个向量 ######
    # # 计算完权重后要将插值点索引的超出边界的值重新计算，如 [0 → end], [end + 1 → 1]等, 用于构建插值矩阵
    # @inbounds for i in eachindex(interIDϕs)
    #     if  interIDϕs[i] < 1
    #         interIDϕs[i] += ntXϕs
    #     elseif interIDϕs[i] > ntXϕs
    #         interIDϕs[i] -= ntXϕs
    #     end
    # end

    # #######################
    # # 将 θ 方向的计算结果保存为大小 (npXθs, ntXθs) 的稀疏化矩阵
    # # 给定的插值点为列坐标，要计算出行坐标
    # rawIDθs     =   repeat(collect(IT, 1:size(interIDθs, 2)); inner = size(interIDθs, 1))
    # """构建 CSC (npXθs, ntXθs) 稀疏矩阵, 用于左乘本层多极子矩阵，在 θ 方向插值"""
    # interpθCSC  =   sparse(rawIDθs, view(interIDθs, :), view(interWθs, :))
    # #######################
    # # 将 ϕ 方向的计算结果保存为大小 (ntXϕs, ntXϕs) 的稀疏化矩阵
    # # 给定的插值点为列坐标，要计算出行坐标
    # rawIDϕs     =   repeat(collect(IT, 1:size(interIDϕs, 2)); inner = size(interIDϕs, 1))
    # """构建 CSC (ntXϕs, ntXϕs) 稀疏矩阵, 用右乘本层多极子矩阵，在 ϕ 方向插值"""
    # interpϕCSC  =   sparse(view(interIDϕs, :), rawIDϕs, view(interWϕs, :))

    # #############新方法，将所有采样点视为 1 个向量，进行插值
    # 采用两步插值:
    # 第单步在 ϕ 方向插值，将采样点数量从本层的 ntXθs*ntXϕs 插值到临时的 ntXθs*npXϕs
    # 第二步在 θ 方向插值，将采样点数量从临时的 ntXθs*npXϕs 插值到父层的 npXθs*npXϕs
    # 父层、本层、插值临时总采样点数
    npSample    =   npXθs*npXϕs
    ntSample    =   ntXθs*ntXϕs
    ntempSample =   ntXθs*npXϕs
    
    ### 第单步 ϕ 方向
    # 将本层所有采样点按 列主排列（列为 θ 方向，行为 ϕ 方向）， 构建本层采样点索引矩阵
    tSampleIndexes  =   reshape(collect(1:ntSample), ntXθs, ntXϕs)
    # 所有临时采样点（待插值点）在本层采样点（插值用到的格点）全局id的向量
    interIDGlobalϕs =   repeat(interIDϕs, inner = (1, ntXθs))
    # 对 θ 方向的所有 ϕ 插值点循环计算其全局插值id
    for itθ in 1:ntXθs
        for ipϕ in 1:npXϕs
            interIDGlobalϕs[:, (ipϕ - 1)*ntXθs + itθ]    .=  [pickCycleVec(interIDϕ, tSampleIndexes[itθ,:])  for interIDϕ in interIDϕs[:,ipϕ]]
        end
    end
    # 上面的全局插值id为列号，下面计算每个元素的行号，从而构建CSC稀疏矩阵
    rawIDϕs =   repeat(collect(IT, 1:ntempSample); inner = nlocalInterp)
    # 权重也要在 tθ 方向重复以匹配所有 θ 处的采样点
    interWGlobalϕs  =   repeat(interWϕs, inner = (1, ntXθs))
    # 构建 稀疏 插值矩阵
    interpϕCSC  =   sparse(rawIDϕs, view(interIDGlobalϕs, :), view(interWGlobalϕs, :))


    #### 第二步 θ 方向插值
    # θ 方向的插值要跨过极点，跨过极点的 真实采样点 与计算权重时的 虚拟采样点 算出来的 θ̂  ϕ̂  是相反的！，因此，为使用真实点替代虚拟点进行插值计算，这些虚拟点对应的权重要取反。
    # 修正越过极点区域的插值权重
    for i in eachindex(interIDθs)
        ((interIDθs[i] < 1) | (interIDθs[i] > ntXθs)) && (interWθs[i] *= -1) 
    end
    # 将 ϕ 插值后的所有临时采样点按 列主排列（列为 θ 方向，行为 ϕ 方向）， 构建本层采样点索引矩阵
    tempSampleIndexes  =   reshape(collect(1:ntempSample), ntXθs, npXϕs)
    # 所有父层采样点（待插值点）在 临时采样点（插值用到的格点）全局id的向量
    interIDGlobalθs =   repeat(interIDθs, outer = (1, npXϕs))
    # 对所有父层采样点（待插值点）循环计算其全局插值id
    halfnpϕ = npXϕs ÷ 2
    for ipϕ in 1:npXϕs
        for ipθ in 1:npXθs
            # 创建包含所有待计算的 nlocalInterp 个插值点id的向量
            inGlobalIDs     =   zeros(IT, nlocalInterp)
            for jInter in 1:nlocalInterp
                # θ方向的插值点id
                interIDθ    =   interIDθs[jInter, ipθ]
                # 待计算的插值点id在所有id数组tempSampleIndexes中的位置
                targetIdxInTempSampleIndexes  =  [interIDθ, ipϕ] 
                # 插值点越过极点的修正
                if (interIDθ < 1) | (interIDθ > ntXθs)
                    # 插值点越过“北”极点
                    if interIDθ < 1
                        targetIdxInTempSampleIndexes[1] = -interIDθ + 1
                    # 插值点越过“南”极点
                    elseif  interIDθ > ntXθs
                        targetIdxInTempSampleIndexes[1] = 2ntXθs + 1 -interIDθ
                    end
                    ## 越过极点后其插值点在另一半球
                    # 待插值点在的 ϕ 值在 [0, π)
                    if ipϕ <= halfnpϕ
                        targetIdxInTempSampleIndexes[2] =  ipϕ + halfnpϕ
                    # 待插值点在的 ϕ 值在 [π, 2π)
                    else
                        targetIdxInTempSampleIndexes[2] =  ipϕ - halfnpϕ
                    end
                end
                inGlobalIDs[jInter] = tempSampleIndexes[targetIdxInTempSampleIndexes...]
            end
            # 将结果写入数组
            interIDGlobalθs[:, (ipϕ - 1)*npXθs + ipθ]    .=  inGlobalIDs
        end
    end
    # 上面的全局插值id为列号，下面计算每个元素的行号，从而构建CSC稀疏矩阵
    rawIDθs =   repeat(collect(IT, 1:npSample); inner = nlocalInterp)
    # 权重也要在 tϕ 方向重复以匹配所有 θ 处的采样点
    interWGlobalθs  =   repeat(interWθs; outer = (1, npXϕs))
    # 构建 稀疏 插值矩阵
    interpθCSC  =   sparse(rawIDθs, view(interIDGlobalθs, :), view(interWGlobalθs, :))

    # 去除插值权重为 0 的部分
    dropzeros!(interpθCSC)
    dropzeros!(interpϕCSC)

    # 将权重保存插值类并返回
    if MLFMAParams.InterpolationMethod == :Lagrange1Step
        return LagrangeInterp1StepInfo(interpθCSC, interpϕCSC)
    else
        return LagrangeInterpInfo(interpθCSC, interpϕCSC)
    end

end

"""
计算一维坐标coora在坐标corrb中的位置
"""
function cooraInCoorb(coora::Vector{T}, coorb::Vector{T}) where {T<:Number}
    
    targetIDs =    fill!(similar(coora, Int), 0)
    @inbounds for i in eachindex(coora)
        a = coora[i]
        for j in eachindex(coorb)
            b = coorb[j]
            if a >= b
                targetIDs[i] = j
                continue
            end
        end
    end

    return targetIDs
end

"""
利用 ϕ 的周期性索引超出上下界的 ϕ 值
"""
function pickϕ(index::Integer, ϕs::Vector{TT}) where {TT<:Real}
    re = zero(TT)
    ln = length(ϕs)
    if  index < 1
        re = (ϕs[index + ln] - 2π)
    elseif index > ln
        re = ϕs[index - ln] + 2π
    else
        re = ϕs[index]
    end
    return re
end


"""
根据循环向量 cycleVec 的周期性索引超出上下界的 index 对应的值
"""
function pickCycleVec(index::Integer, cycleVec::Vector{T}) where {T<:Real}

    ln = length(cycleVec)
    if  index < 1
        re = cycleVec[index + ln]
    elseif index > ln
        re = cycleVec[index - ln]
    else
        re = cycleVec[index]
    end
    return re
end

"""
利用 θ 在极点附近的对称性计算索引超出上下界的 θ 值
"""
function pickθ(index::Integer, θs::Vector{T}) where {T<:Real}
    re = zero(T)
    ln = length(θs)
    if  index < 1
        re = -θs[-index + 1]
    elseif index > ln
        re = 2pi - θs[2ln + 1 - index]
    else
        re = θs[index]
    end
    return re
end

