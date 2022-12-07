"""
计算三角形面片上的加权电流。
电流基函数公式为：Jₙ = Iₙfₙ
同一个三角形面片上存在三个基函数，因此
Jₜ = ∑ₜₙ₌₁³ Iₜₙfₜₙ
输入：
ICoeff          ::Vector{Complex{FT}}  计算得到的电流系数
trianglesInfo   ::Vector{TriangleInfo{IT, FT}}，三角形信息
输出值:
Jtri         ::Marrix{Complex{FT}}, 三角形上加权后的电流
"""
function geoElectricJCal(ICoeff::Vector{CT}, trianglesInfo::AbstractVector{TriangleInfo{IT, FT}}, ::Type{BFT}= VSBFTypes.sbfType) where {IT<:Integer, FT<:Real, CT<:Complex{FT}, BFT<:RWG}
    # 三角形数
    ntri    =   length(trianglesInfo)
    # 给结果预分配内存
    Jtris   =   zeros(CT, 3, ntri)
    # 是否为偏置数组
    isoffset    =   isa(trianglesInfo, OffsetVector)
    geoIdx      =   eachindex(trianglesInfo)
    # 对六面体循环计算
    for it in 1:ntri
        # 第it个六面体
        idx =   isoffset ? (geoIdx.offset + it) : it
        tri =   trianglesInfo[idx]
        # 该三角形所在的三个基函数id
        tms     =   tri.inBfsID
        # 相关的电流
        Jtri =   @view Jtris[:, it]
        # 对高斯求积点循环
        for gi in 1:GQPNTri
            # 初始化电流
            Jtritemp   =   zero(MVec3D{CT})
            # 采样点
            rgi =   getGQPTri(tri, gi)
            for mi in 1:3
                # 判断边是不是基函数
                m = tms[mi] 
                m == 0 && continue
                # 自由点到采样点向量ρ
                ρgm =   rgi - tri.vertices[:, mi]
                # 将本基函数（边）的电流累加到结果
                Jtritemp  +=   ICoeff[m]*tri.edgel[mi]*ρgm
            end
            Jtritemp  *=   TriGQInfo.weight[gi]
            # 结果累加
            Jtri .+=   Jtritemp
        end #for gi
        # 结果修正
        Jtri     ./=   2tri.area
    end #for it

    return Jtris
end # function


"""
计算四面体上的电流。
分片常数基 PWC 基函数
电流基函数公式为：Jₙ = κₙIₙfₙ
同一个四面体面片上存在 x̂, ŷ, ẑ 方向的三个基函数，因此
Jₜ = κₜ ∑ₜₙ₌₁³Iₜₙfₜₙ
输入：
ICoeff          ::Vector{Complex{FT}}  计算得到的电流系数
tetrasInfo      ::Vector{TetrahedraInfo{IT, FT, CT}}, 
输出值:
Jtetra          ::Marrix{Complex{FT}}, 四面体上的电流
"""
function geoElectricJCal(ICoeff::Vector{CT}, tetrasInfo::AbstractVector{TetrahedraInfo{IT, FT, CT}}, ::Type{BFT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}, BFT<:PWC}
    # 四面体数
    ntetra  =   length(tetrasInfo)
    # 给结果预分配内存
    Jtetras =   zeros(CT, 3, ntetra)
    # 是否为偏置数组
    isoffset    =   isa(tetrasInfo, OffsetVector)
    geoIdx      =   eachindex(tetrasInfo)
    # 判断体电流的离散方式
    discreteJ::Bool =   (SimulationParams.discreteVar === "J")
    # 对四面体循环计算
    for it in 1:ntetra
        # 第it个六面体
        idx     =   isoffset ? (geoIdx.offset + it) : it
        tetra   =   tetrasInfo[idx]
        # 该四面体所在的三个基函数id
        tms     =   tetra.inBfsID
        # 将本基函数的电流累加到结果
        if discreteJ
            Jtetras[:, it] .+=   view(ICoeff, tms)
        else
            Jtetras[:, it] .+=   tetra.κ .* view(ICoeff, tms)
        end
    end #for ti

    return Jtetras
end # function

"""
计算四面体上的电流。
分片常数基 SWG 基函数
电流基函数公式为：Jₙ = κₙIₙfₙ
同一个四面体面片上存在四个基函数，因此
Jₜ = ∑ₜₙ₌₁⁴ Iₜₙfₜₙ
输入：
ICoeff          ::Vector{Complex{FT}}  计算得到的电流系数
tetrasInfo      ::Vector{TetrahedraInfo{IT, FT, CT}}, 
输出值:
Jtetra          ::Marrix{Complex{FT}}, 四面体上加权后的电流
"""
function geoElectricJCal(ICoeff::Vector{CT}, tetrasInfo::AbstractVector{TetrahedraInfo{IT, FT, CT}}, ::Type{BFT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}, BFT<:SWG}
    # 四面体数
    ntetra  =   length(tetrasInfo)
    # 给结果预分配内存
    Jtetras =   zeros(CT, 3, ntetra)
    # 是否为偏置数组
    isoffset    =   isa(tetrasInfo, OffsetVector)
    geoIdx      =   eachindex(tetrasInfo)
    # 判断体电流的离散方式
    discreteJ::Bool =   (SimulationParams.discreteVar === "J")
    # 对六面体循环计算
    for it in 1:ntetra
        # 第it个六面体
        idx = isoffset ? (geoIdx.offset + it) : it
        tetra   =   tetrasInfo[idx]
        # 该四面体所在的三个基函数id
        tms     =   tetra.inBfsID
        # 相关的电流
        Jtetra  =   @view Jtetras[:, it]
        # 对高斯求积点循环
        for gi in 1:GQPNTetra
            # 初始化电流
            Jtemp   =   zero(MVec3D{CT})
            # 采样点
            rgi     =   getGQPTetra(tetra, gi)
            for mi in 1:4
                # 判断边是不是基函数
                m   =   tms[mi]
                # 自由点到采样点向量ρ
                ρgm =   rgi - tetra.vertices[:, mi]
                # 将本基函数（边）的电流累加到结果
                Jtemp   +=   ICoeff[m]*tetra.facesArea[mi]*ρgm
            end
            Jtemp    *=   TetraGQInfo.weight[gi]
            # 结果累加
            Jtetra  .+=   Jtemp
        end #for gi
        # 结果修正
        if discreteJ
            Jtetra      ./=   3tetra.volume
        else
            Jtetra      .*=   tetra.κ/(3tetra.volume)
        end
    end #for ti

    return Jtetras
end # function

"""
计算六面体上的电流。
分片常数基 PWC 基函数
电流基函数公式为：Jₙ = κₙIₙfₙ
同一个六面体面片上存在 x̂, ŷ, ẑ 方向的三个基函数，因此
Jₜ = κₜ ∑ₜₙ₌₁³Iₜₙfₜₙ
输入：
ICoeff          ::Vector{Complex{FT}}  计算得到的电流系数
hexasInfo      ::Vector{HexahedraInfo{IT, FT, CT}}, 
输出值:
Jhexa          ::Marrix{Complex{FT}}, 六面体上的电流
"""
function geoElectricJCal(ICoeff::Vector{CT}, hexasInfo::AbstractVector{HexahedraInfo{IT, FT, CT}}, ::Type{BFT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}, BFT<:PWC}
    # 六面体数
    nhexa  =   length(hexasInfo)
    # 给结果预分配内存
    Jhexas =   zeros(CT, 3, nhexa)
    # 是否为偏置数组
    isoffset    =   isa(hexasInfo, OffsetVector)
    geoIdx      =   eachindex(hexasInfo)
    # 判断体电流的离散方式
    discreteJ::Bool =   (SimulationParams.discreteVar === "J")
    # 对六面体循环计算
    for it in 1:nhexa
        # 第it个六面体
        idx = isoffset ? (geoIdx.offset + it) : it
        hexa   =   hexasInfo[idx]
        # 该六面体所在的三个基函数id
        tms     =   hexa.inBfsID
        # 将本基函数的电流累加到结果
        if discreteJ
            Jhexas[:, it] .+=   view(ICoeff, tms)
        else
            Jhexas[:, it] .+=   hexa.κ .* view(ICoeff, tms)
        end
    end #for it

    return Jhexas
end # function


"""
计算六面体上的电流。
分片常数基 RBF 基函数
电流基函数公式为：Jₙ = κₙIₙfₙ
同一个六面体面片上存在 6 或 3 个基函数，因此
Jₜ = ∑ₜₙ₌₁ Iₜₙfₜₙ
输入：
ICoeff          ::Vector{Complex{FT}}  计算得到的电流系数
hexasInfo       ::Vector{HexahedraInfo{IT, FT, CT}}, 
输出值:
Jhexa          ::Marrix{Complex{FT}}, 六面体上加权后的电流
"""
function geoElectricJCal(ICoeff::Vector{CT}, hexasInfo::AbstractVector{HexahedraInfo{IT, FT, CT}}, ::Type{BFT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}, BFT<:RBF}
    # 六面体数
    nhexa   =   length(hexasInfo)
    # 给结果预分配内存
    Jhexas  =   zeros(CT, 3, nhexa)
    # 是否为偏置数组
    isoffset    =   isa(hexasInfo, OffsetVector)
    geoIdx      =   eachindex(hexasInfo)
    # 判断体电流的离散方式
    discreteJ::Bool =   (SimulationParams.discreteVar === "J")
    # 对六面体循环计算
    for it in 1:nhexa
        # 第it个六面体
        idx = isoffset ? (geoIdx.offset + it) : it
        hexa   =   hexasInfo[idx]
        # 该六面体所在的三个基函数id
        tms    =   hexa.inBfsID
        # 相关的电流
        Jhexa  =   @view Jhexas[:, it]
        # 自由端
        freeVms =  Vector{MMatrix{3, GQPNQuad, FT, 3*GQPNQuad}}(undef, 6)
        @inbounds for mi in 1:6
            freeVms[mi] =   getFreeVns(hexa, mi)
        end
        # 对高斯求积点循环
        for gi in 1:GQPNHexa
            # 初始化电流
            Jtemp   =   zero(MVec3D{CT})
            # 采样点
            rgi     =   getGQPHexa(hexa, gi)
            idm3D   =   GQ1DID2GQ3DIDVector[gi]
            for mi in 1:6
                # 判断边是不是基函数
                m   =   tms[mi]
                idm =   getFreeVIDFromGQ3DID(idm3D, mi)
                # 自由点到采样点向量ρ
                ρgm =   rgi - view(freeVms[mi], :, idm)
                # 将本基函数（边）的电流累加到结果
                Jtemp   +=   ICoeff[m]*hexa.facesArea[mi]*ρgm
            end
            Jtemp    *=   HexaGQInfo.weight[gi]
            # 结果累加
            Jhexa  .+=   Jtemp
        end #for gi
        # 结果修正
        if discreteJ
            Jhexa      ./=   hexa.volume
        else
            Jhexa      .*=   hexa.κ/(hexa.volume)
        end
    end #for ti

    return Jhexas
end # function

"""
计算给定三角形面片位置 r 处的电流
电流基函数公式为：Jₙ = Iₙfₙ
同一个三角形面片上存在三个基函数，因此
Jₜ = ∑ₜₙ₌₁³ Iₜₙfₜₙ
输入：
r               ::Vec3D{FT}
ICoeff          ::Vec3D{CT}  三角形上的三个基函数的计算得到的电流系数
triangleInfo    ::TriangleInfo{IT, FT}，三角形信息
输出值:
Jtrir           ::Complex{FT}, 三角形上加权后的电流
"""
function electricJCal(r::Vec3D{FT}, ICoefftri::Vec3D{CT}, tri::TriangleInfo{IT, FT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}}
    # 给结果预分配内存
    Jr   =   zero(MVec3D{CT})
    for mi in 1:3
        # 自由点到采样点向量ρ
        ρgm     =   r - tri.vertices[:, mi]
        # 将本基函数（边）的电流累加到结果
        Jr     +=   ICoefftri[mi]*tri.edgel[mi]*ρgm
    end
    # 结果修正
    Jr ./= 2tri.area

    return Jr
end # function


"""
计算所有三角形上的高斯求积点电流权重乘积 JₙᵢWᵢ
电流基函数公式为：Jₙ = Iₙfₙ
同一个三角形面片上存在三个基函数，因此
JₙᵢWᵢ = ∑ₜₙ₌₁³ Iₜₙlₜₙ/2Sₜₙ
输入：
ICoeff          ::Vector{Complex{FT}}  计算得到的电流系数
trianglesInfo   ::Vector{TriangleInfo{IT, FT}}，三角形信息
输出值:
Jtri         ::Marrix{Complex{FT}}, 三角形上加权后的电流
"""
@fastmath function electricJCal(ICoeff::Vector{CT}, trianglesInfo::AbstractVector{TriangleInfo{IT, FT}}
                        ) where {IT<:Integer, FT<:Real, CT<:Complex{FT}}

    # 三角形数
    ntri    =   length(trianglesInfo)
    # 给结果预分配内存
    Jtris   =   zeros(CT, 3, GQPNTri, ntri)

    # Progress Meter
    pmeter  =   Progress(ntri, "Calculating J on triangles' gaussquad points ($GQPNTri × $ntri)")

    # 对三角形循环计算
    @threads for ti in 1:ntri
        # 第ti个三角形
        @inbounds tri     =   trianglesInfo[ti]
        # 该三角形所在的三个基函数id
        tms     =   tri.inBfsID
        # 相关的电流
        @inbounds Jtri    =   @view Jtris[:, :, ti]
        # 对高斯求积点循环
        @inbounds for gi in 1:GQPNTri
            # 初始化电流
            Jtritemp   =   zero(MVec3D{CT})
            # 采样点
            rgi =   getGQPTri(tri, gi)
            for mi in 1:3
                # 判断边是不是基函数
                m = tms[mi] 
                m == 0 && continue
                # 自由点到采样点向量ρ
                ρgm = rgi - tri.vertices[:, mi]
                # 将本基函数（边）的电流累加到结果
                Jtritemp +=  ICoeff[m]*tri.edgel[mi]*ρgm
            end
            Jtritemp    /=  2tri.area
            # 结果累加
            Jtri[:, gi] .=  Jtritemp
        end #for gi

        # 更新进度条
        next!(pmeter)

    end #for ti

    return Jtris

end
