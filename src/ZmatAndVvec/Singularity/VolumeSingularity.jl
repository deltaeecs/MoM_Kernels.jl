@doc raw"""
    volumeSingularityIgIvecg(rtveclc::AbstractVector{FT}, volumeCell::TetrahedraInfo{IT, FT, CT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}}
    volumeSingularityIgIvecg(rtveclc::AbstractVector{FT}, volumeCell::HexahedraInfo{IT, FT, CT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}}

计算场点`rgt`在体网格`volumeCell`上的奇异性。
计算结果为：
```math
\begin{aligned}
I_{gV}  &= \int{g(R)dV'}\\
        &= -\sum_{S_i}{d_i\sum_{n}^{SglrOrder}{\frac{coeffgreen(n)}{n+2}I_{RS}^{n-1}}}\\
\boldsymbol{I}_{gV}  &= \int{\boldsymbol{R}g(R)dV'}\\
        &= -\sum_{S_i}{\hat{\bm{n}}_i \sum_{n=0}^{SglrOrder}{\frac{coeffgreen(n)}{n+1}I^{n+1}_{RS}}}\\
\end{aligned}
```
"""
function volumeSingularityIgIvecg(rtveclc::AbstractVector{FT}, volumeCell::TetrahedraInfo{IT, FT, CT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}}
    # 体元的面
    faces   =   volumeCell.faces
    # 体元的面数
    nFaces  =   length(faces)
    # 预分配 Iᵣⁿ、 ISᵣ 数组
    ISᵣ     =   OffsetArray(zero(MVector{SglrOrder, FT}), -2)
    Ilᵣ     =   OffsetArray(zero(MVector{SglrOrder, FT}), -2)
    # 结果
    IVg     =   zero(CT)
    IvecVg  =   zero(MVec3D{CT})

    # 预分配临时变量以加速
    r0      =   zero(MVec3D{FT})
    p02ivec =   zero(MVec3D{FT})
    
    # 对高斯求积点循环
    for iface in 1:nFaces
        # 第 iface 个面
        face    =   faces[iface]
        # 置零
        ISᵣ .=  0
        # 该积分点在到源三角形上的距离（带正负，以源三角形法向为正向）
        # dts     =   view(volumeCell.facesn̂, :, iface) ⋅ (rtveclc .- view(face.vertices, :, 1))
        dts  =   zero(FT)
        for ii in 1:3
            dts += volumeCell.facesn̂[ii, iface] * (rtveclc[ii] - face.vertices[ii, 1])
        end
        # 投影点
        # r0     .=   rtveclc .- dts .* view(volumeCell.facesn̂, :, iface)
        for ii in 1:3
            r0[ii]   =   rtveclc[ii] .- dts .* volumeCell.facesn̂[ii, iface]
        end
        
        # 距离的绝对值、平方
        dts²    =   dts*dts
        dtsAbs  =   abs(dts)
        
        # 边的信息
        edgel   =   face.edgel
        edgev̂   =   face.edgev̂
        edgen̂   =   face.edgen̂
        # 边数
        nEdge   =   length(edgel)
        # 对 nEdge 个边循环
        for edgei in 1:nEdge
            # 置零避免累加错误
            Ilᵣ .=  0
            # 构成该边的第一个点
            edgeNodei⁻  =   view(face.vertices, :, MoM_Basics.EDGEVmINTriVsID[edgei])
            # 该边边长
            lj      =   edgel[edgei]
            # 每个边的局部坐标下的 ljⱼ⁺ \lⱼ⁻
            # lⱼ⁻     =   (edgeNodei⁻[1] - r0[1]) * edgev̂[1, edgei] + (edgeNodei⁻[2] - r0[2]) * edgev̂[2, edgei] + (edgeNodei⁻[3] - r0[3]) * edgev̂[3, edgei]
            lⱼ⁻  =   zero(FT)
            for ii in 1:3
                lⱼ⁻ += (edgeNodei⁻[ii] - r0[ii]) * edgev̂[ii, edgei]
            end
            lⱼ⁺ =   lⱼ⁻ + lj
            # 投影点 r0 到各个边的垂足的向量 p02ivec
            # p02ivec    .=   edgeNodei⁻ .- lⱼ⁻ .* view(edgev̂, :, edgei) .- r0
            for ii in 1:3
                p02ivec[ii] = edgeNodei⁻[ii] - lⱼ⁻ * edgev̂[ii, edgei] - r0[ii]
            end
            # 该向量长的平方、该向量长
            @views p02jl   =   p02ivec ⋅ edgen̂[:, edgei]
            p02jl²  =   abs2(p02jl)
            # 将p02ivec化为单位向量
            p02ivec ./=   p02jl
            # 场点到在该边上投影点的距离
            R0²     =   p02jl² + dts²
            # 场点到在该边上正、负端点的距离
            R⁺      =   sqrt(abs2(lⱼ⁺) + R0²)
            R⁻      =   sqrt(abs2(lⱼ⁻) + R0²)
            
            # 此处为避免数值误差，将阈值 ϵl 设定为 1e-2 边长
            let fⱼ = zero(FT), βⱼ =  zero(FT), ϵl = 1e-3lj
                if abs(p02jl) < ϵl
                    # 投影点过于靠近该边即 p02jl = 0 时， βⱼ = 0
                    if dtsAbs < ϵl
                        # 视为积分点与边重合，此时 βⱼ = 0， Kᵣ不变
                        continue
                    else
                        # 视为投影点与边重合，此时 βⱼ = 0， Kᵣ不变
                        fⱼ      =   log((lⱼ⁺ + R⁺)/(lⱼ⁻ + R⁻))
                    end
                else
                    fⱼ      =   log((lⱼ⁺ + R⁺)/(lⱼ⁻ + R⁻))
                    if dtsAbs < ϵl
                        # 视为积分点与源三角形面重合，此时 βⱼ ≠ 0
                        βⱼ      =   atan((p02jl*lⱼ⁺)/R0²) - atan((p02jl*lⱼ⁻)/R0²)
                        # ISᵣ累加
                        ISᵣ[-1] +=   p02jl*fⱼ
                    else
                        # 正常处理
                        βⱼ      =   atan((p02jl*lⱼ⁺)/(R0² + dtsAbs*R⁺)) - atan((p02jl*lⱼ⁻)/(R0² + dtsAbs*R⁻))
                        # ISᵣ累加
                        ISᵣ[-1] +=   p02jl*fⱼ - dtsAbs*βⱼ
                    end
                end #if p02jl
                # 计算 Ilᵣ 所有项
                Ilᵣ[-1]  =   fⱼ
                Ilᵣ[0]   =   lj
                R⁺ⁿ =   one(FT)
                R⁻ⁿ =   one(FT) 
                for n in 1:(SglrOrder-2)
                    R⁺ⁿ    *=   R⁺
                    R⁻ⁿ    *=   R⁻
                    Ilᵣ[n]  =   (lⱼ⁺*R⁺ⁿ - lⱼ⁻*R⁻ⁿ  + n*R0²*Ilᵣ[n-2])/(n+1)
                end # n
                # 计算 ISᵣ 中与边相关的项
                for n in 1:(SglrOrder-2)
                    ISᵣ[n]  +=   p02jl * Ilᵣ[n]
                end # n
            end #let
        end # edgei
        # 其它 ISᵣ 相关项
        ISᵣ[0]   =   abs(volumeCell.facesArea[iface])
        for n in 1:(SglrOrder-2)
            ISᵣ[n]  +=  n*dts²*ISᵣ[n-2]
            ISᵣ[n]  /=  n + 2
        end
        # 累加 IVg 结果
        ISg      =   zero(CT)
        for n in 0:(SglrOrder-1)
            ISg -=  SSCgdivnp2[n]*ISᵣ[n-1]
        end
        IVg     +=  dts*ISg

        # 累加 IvecVg 结果
        IvecVgi  =  zero(CT)
        for n in 0:(SglrOrder-3)
            IvecVgi -=  SSCgdivnp1[n]*ISᵣ[n+1]
        end
        @views IvecVg .+=  volumeCell.facesn̂[:, iface] * IvecVgi

    end # iface

    return IVg, IvecVg

end
function volumeSingularityIgIvecg(rtveclc::AbstractVector{FT}, volumeCell::HexahedraInfo{IT, FT, CT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}}
    # 体元的面
    faces   =   volumeCell.faces
    # 体元的面数
    nFaces  =   length(faces)
    # 预分配 Iᵣⁿ、 ISᵣ 数组
    ISᵣ     =   OffsetArray(zero(MVector{SglrOrder, FT}), -2)
    Ilᵣ     =   OffsetArray(zero(MVector{SglrOrder, FT}), -2)
    # 结果
    IVg     =   zero(CT)
    IvecVg  =   zero(MVec3D{CT})

    # 预分配临时变量以加速
    p02ivec         =   zero(MVec3D{FT})
    r0  =   zero(MVec3D{FT})
    
    # 对高斯求积点循环
    for iface in 1:nFaces
        # 第 iface 个面
        face    =   faces[iface]
        # 置零
        ISᵣ .=  0
        # 该积分点在到源三角形上的距离（带正负，以源三角形法向为正向）
        # dts     =   view(volumeCell.facesn̂, :, iface) ⋅ (rtveclc .- view(face.vertices, :, 1))
        dts  =   zero(FT)
        for ii in 1:3
            dts += volumeCell.facesn̂[ii, iface] * (rtveclc[ii] - face.vertices[ii, 1])
        end
        # 投影点
        # r0     .=   rtveclc .- dts .* view(volumeCell.facesn̂, :, iface)
        # @views r0   =   rtveclc .- dts .* volumeCell.facesn̂[:, iface]
        for ii in 1:3
            r0[ii]   =   rtveclc[ii] .- dts .* volumeCell.facesn̂[ii, iface]
        end
        
        # 距离的绝对值、平方
        dts²    =   dts*dts
        dtsAbs  =   abs(dts)
        
        # 边的信息
        edgel   =   face.edgel
        edgev̂   =   face.edgev̂
        edgen̂   =   face.edgen̂
        # 边数
        nEdge   =   length(edgel)
        # 对 nEdge 个边循环
        for edgei in 1:nEdge
            # 置零避免累加错误
            Ilᵣ .=  0
            # 构成该边的第一个点
            edgeNodei⁻  =   view(face.vertices, :, edgei)
            # 该边边长
            lj      =   edgel[edgei]
            # 每个边的局部坐标下的 ljⱼ⁺ \lⱼ⁻
            lⱼ⁻  =   zero(FT)
            for ii in 1:3
                lⱼ⁻ += (edgeNodei⁻[ii] - r0[ii]) * edgev̂[ii, edgei]
            end
            lⱼ⁺ =   lⱼ⁻ + lj
            # 投影点 r0 到各个边的垂足的向量 p02ivec
            # @views p02ivec .=   edgeNodei⁻ .- lⱼ⁻ * edgev̂[:, edgei] .- r0
            for ii in 1:3
                p02ivec[ii] = edgeNodei⁻[ii] - lⱼ⁻ * edgev̂[ii, edgei] - r0[ii]
            end
            # 该向量长的平方、该向量长
            @views p02jl   =   p02ivec ⋅ edgen̂[:, edgei]
            p02jl²  =   abs2(p02jl)
            # 将p02ivec化为单位向量
            p02ivec ./=   p02jl
            # 场点到在该边上投影点的距离
            R0²     =   p02jl² + dts²
            # 场点到在该边上正、负端点的距离
            R⁺      =   sqrt(abs2(lⱼ⁺) + R0²)
            R⁻      =   sqrt(abs2(lⱼ⁻) + R0²)
            
            # 此处为避免数值误差，将阈值 ϵl 设定为 1e-2 边长
            let fⱼ = zero(FT), βⱼ =  zero(FT), ϵl = 1e-3lj
                if abs(p02jl) < ϵl
                    # 投影点过于靠近该边即 p02jl = 0 时， βⱼ = 0
                    if dtsAbs < ϵl
                        # 视为积分点与边重合，此时 βⱼ = 0， Kᵣ不变
                        continue
                    else
                        # 视为投影点与边重合，此时 βⱼ = 0， Kᵣ不变
                        fⱼ      =   log((lⱼ⁺ + R⁺)/(lⱼ⁻ + R⁻))
                    end
                else
                    fⱼ      =   log((lⱼ⁺ + R⁺)/(lⱼ⁻ + R⁻))
                    if dtsAbs < ϵl
                        # 视为积分点与源三角形面重合，此时 βⱼ ≠ 0
                        βⱼ      =   atan((p02jl*lⱼ⁺)/R0²) - atan((p02jl*lⱼ⁻)/R0²)
                        # ISᵣ累加
                        ISᵣ[-1] +=   p02jl*fⱼ
                    else
                        # 正常处理
                        βⱼ      =   atan((p02jl*lⱼ⁺)/(R0² + dtsAbs*R⁺)) - atan((p02jl*lⱼ⁻)/(R0² + dtsAbs*R⁻))
                        # ISᵣ累加
                        ISᵣ[-1] +=   p02jl*fⱼ - dtsAbs*βⱼ
                    end
                end #if p02jl
                # 计算 Ilᵣ 所有项
                Ilᵣ[-1]  =   fⱼ
                Ilᵣ[0]   =   lj
                R⁺ⁿ =   one(FT)
                R⁻ⁿ =   one(FT) 
                for n in 1:(SglrOrder-2)
                    R⁺ⁿ    *=   R⁺
                    R⁻ⁿ    *=   R⁻
                    Ilᵣ[n]  =   (lⱼ⁺*R⁺ⁿ - lⱼ⁻*R⁻ⁿ  + n*R0²*Ilᵣ[n-2])/(n+1)
                end # n
                # 计算 ISᵣ 中与边相关的项
                for n in 1:(SglrOrder-2)
                    ISᵣ[n]  +=   p02jl * Ilᵣ[n]
                end # n
            end #let
        end # edgei
        # 其它 ISᵣ 相关项
        ISᵣ[0]   =   abs(volumeCell.facesArea[iface])
        for n in 1:(SglrOrder-2)
            ISᵣ[n]  +=  n*dts²*ISᵣ[n-2]
            ISᵣ[n]  /=  n + 2
        end
        # 累加 IVg 结果
        ISg      =   zero(CT)
        for n in 0:(SglrOrder-1)
            ISg -=  SSCgdivnp2[n]*ISᵣ[n-1]
        end
        IVg     +=  dts*ISg

        # 累加 IvecVg 结果
        IvecVgi  =  zero(CT)
        for n in 0:(SglrOrder-3)
            IvecVgi -=  SSCgdivnp1[n]*ISᵣ[n+1]
        end
        @views IvecVg .+=  volumeCell.facesn̂[:, iface] * IvecVgi

    end # iface

    return IVg, IvecVg

end

@doc raw"""
    volumeSingularityIg(rtveclc::AbstractVector{FT}, volumeCell::HexahedraInfo{IT, FT, CT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}}

计算场点`rgt`在体网格`volumeCell`上的奇异性。
计算结果为：
```math
\begin{aligned}
I_{gV}  &= \int{g(R)dV'}\\
        &= -\sum_{S_i}{d_i\sum_{n}^{SglrOrder}{\frac{coeffgreen(n)}{n+2}I_{RS}^{n-1}}}\\
\end{aligned}
```
"""
function volumeSingularityIg(rtveclc::AbstractVector{FT}, volumeCell::HexahedraInfo{IT, FT, CT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}}
    # 体元的面
    faces   =   volumeCell.faces
    # 体元的面数
    nFaces  =   length(faces)
    # 预分配 Iᵣⁿ、 ISᵣ 数组
    ISᵣ     =   OffsetArray(zero(MVector{SglrOrder, FT}), -2)
    Ilᵣ     =   OffsetArray(zero(MVector{SglrOrder, FT}), -2)
    # 结果
    IVg     =   zero(CT)

    # 预分配临时变量以加速
    p02ivec         =   zero(MVec3D{FT})
    r0  =   zero(MVec3D{FT})
    
    # 对高斯求积点循环
    for iface in 1:nFaces
        # 第 iface 个面
        face    =   faces[iface]
        # 置零
        ISᵣ .=  0
        # 该积分点在到源三角形上的距离（带正负，以源三角形法向为正向）
        dts  =   zero(FT)
        for ii in 1:3
            dts += volumeCell.facesn̂[ii, iface] * (rtveclc[ii] - face.vertices[ii, 1])
        end
        # 投影点
        # @views r0   =   rtveclc .- dts .* volumeCell.facesn̂[:, iface]
        for ii in 1:3
            r0[ii]   =   rtveclc[ii] .- dts .* volumeCell.facesn̂[ii, iface]
        end
        
        # 距离的绝对值、平方
        dts²    =   dts*dts
        dtsAbs  =   abs(dts)
        
        # 边的信息
        edgel   =   face.edgel
        edgev̂   =   face.edgev̂
        edgen̂   =   face.edgen̂
        # 边数
        nEdge   =   length(edgel)
        # 对 nEdge 个边循环
        for edgei in 1:nEdge
            # 置零避免累加错误
            Ilᵣ .=  0
            # 构成该边的第一个点
            edgeNodei⁻  =   view(face.vertices, :, edgei)
            # 该边边长
            lj      =   edgel[edgei]
            # 每个边的局部坐标下的 ljⱼ⁺ \lⱼ⁻
            lⱼ⁻  =   zero(FT)
            for ii in 1:3
                lⱼ⁻ += (edgeNodei⁻[ii] - r0[ii]) * edgev̂[ii, edgei]
            end
            lⱼ⁺ =   lⱼ⁻ + lj
            # 投影点 r0 到各个边的垂足的向量 p02ivec
            for ii in 1:3
                p02ivec[ii] = edgeNodei⁻[ii] - lⱼ⁻ * edgev̂[ii, edgei] - r0[ii]
            end
            # 该向量长的平方、该向量长
            @views p02jl   =   p02ivec ⋅ edgen̂[:, edgei]
            p02jl²  =   abs2(p02jl)
            # 将p02ivec化为单位向量
            p02ivec ./=   p02jl
            # 场点到在该边上投影点的距离
            R0²     =   p02jl² + dts²
            # 场点到在该边上正、负端点的距离
            R⁺      =   sqrt(abs2(lⱼ⁺) + R0²)
            R⁻      =   sqrt(abs2(lⱼ⁻) + R0²)
            
            # 此处为避免数值误差，将阈值 ϵl 设定为 1e-2 边长
            let fⱼ = zero(FT), βⱼ =  zero(FT), ϵl = 1e-3lj
                if abs(p02jl) < ϵl
                    # 投影点过于靠近该边即 p02jl = 0 时， βⱼ = 0
                    if dtsAbs < ϵl
                        # 视为积分点与边重合，此时 βⱼ = 0， Kᵣ不变
                        continue
                    else
                        # 视为投影点与边重合，此时 βⱼ = 0， Kᵣ不变
                        fⱼ      =   log((lⱼ⁺ + R⁺)/(lⱼ⁻ + R⁻))
                    end
                else
                    fⱼ      =   log((lⱼ⁺ + R⁺)/(lⱼ⁻ + R⁻))
                    if dtsAbs < ϵl
                        # 视为积分点与源三角形面重合，此时 βⱼ ≠ 0
                        βⱼ      =   atan((p02jl*lⱼ⁺)/R0²) - atan((p02jl*lⱼ⁻)/R0²)
                        # ISᵣ累加
                        ISᵣ[-1] +=   p02jl*fⱼ
                    else
                        # 正常处理
                        βⱼ      =   atan((p02jl*lⱼ⁺)/(R0² + dtsAbs*R⁺)) - atan((p02jl*lⱼ⁻)/(R0² + dtsAbs*R⁻))
                        # ISᵣ累加
                        ISᵣ[-1] +=   p02jl*fⱼ - dtsAbs*βⱼ
                    end
                end #if p02jl
                # 计算 Ilᵣ 所有项
                Ilᵣ[-1]  =   fⱼ
                Ilᵣ[0]   =   lj
                R⁺ⁿ =   one(FT)
                R⁻ⁿ =   one(FT) 
                for n in 1:(SglrOrder-2)
                    R⁺ⁿ    *=   R⁺
                    R⁻ⁿ    *=   R⁻
                    Ilᵣ[n]  =   (lⱼ⁺*R⁺ⁿ - lⱼ⁻*R⁻ⁿ  + n*R0²*Ilᵣ[n-2])/(n+1)
                end # n
                # 计算 ISᵣ 中与边相关的项
                for n in 1:(SglrOrder-2)
                    ISᵣ[n]  +=   p02jl * Ilᵣ[n]
                end # n
            end #let
        end # edgei
        # 其它 ISᵣ 相关项
        ISᵣ[0]   =   abs(volumeCell.facesArea[iface])
        for n in 1:(SglrOrder-2)
            ISᵣ[n]  +=  n*dts²*ISᵣ[n-2]
            ISᵣ[n]  /=  n + 2
        end
        # 累加 IVg 结果
        ISg      =   zero(CT)
        for n in 0:(SglrOrder-1)
            ISg -=  SSCgdivnp2[n]*ISᵣ[n-1]
        end
        IVg     +=  dts*ISg

    end # iface

    return IVg

end

@doc raw"""
    volumeSingularityLOpDyad(rtveclc::AbstractVector{FT}, volumeCell::TetrahedraInfo{IT, FT, CT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}}
    volumeSingularityLOpDyad(rtveclc::AbstractVector{FT}, volumeCell::HexahedraInfo{IT, FT, CT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}}
    
计算场点`rtveclc`在体网格`volumeCell`上的并矢格林函数奇异性。
计算结果为：
```math
\begin{aligned}
\overline{I}_{V}  &= \int{(k^2 I + ∇∇)G(R) dV'}
\end{aligned}
```
"""
function volumeSingularityLOpDyad(rtveclc::AbstractVector{FT}, volumeCell::TetrahedraInfo{IT, FT, CT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}}
    # 体元的面
    faces   =   volumeCell.faces
    # 体元的面数
    nFaces  =   length(faces)
    # 常数项
    k²      =   Params.k²
    # 预分配 Iᵣⁿ、 ISᵣ, K̂ᵣⁿ 结果数组
    ISᵣ     =   OffsetArray(zero(MVector{SglrOrder+1, FT}), -2)
    Ilᵣ     =   OffsetArray(zero(MVector{SglrOrder+1, FT}), -2)
    # 一些要计算的量
    dC₁K::CT    =   zero(CT)
    ûC₂Ipn̂dC₂Kvec   =   zero(MVec3D{CT})
    n̂ûC₂Ipn̂dC₂KDyad =   zero(MMatrix{3, 3, CT})
    re              =   zero(MMatrix{3, 3, CT})
    # 预分配临时变量以加速
    p02ivec         =   zero(MVec3D{FT})
    r0              =   zero(MVec3D{FT})
    
    # 对高斯求积点循环
    for iface in 1:nFaces
        # 第 iface 个面
        face    =   faces[iface]
        # 置零
        ISᵣ .=  0
        # 该积分点在到源三角形上的距离（带正负，以源三角形法向为正向）
        # @views dts  =   volumeCell.facesn̂[:, iface] ⋅ (rtveclc .- face.vertices[:, 1])
        dts  =   zero(FT)
        for ii in 1:3
            dts += volumeCell.facesn̂[ii, iface] * (rtveclc[ii] - face.vertices[ii, 1])
        end
        # 投影点
        # @views r0   =   rtveclc .- dts .* volumeCell.facesn̂[:, iface]
        for ii in 1:3
            r0[ii]   =   rtveclc[ii] .- dts .* volumeCell.facesn̂[ii, iface]
        end
        
        # 距离的绝对值、平方
        dts²    =   dts*dts
        dtsAbs  =   abs(dts)
        
        # 边的信息
        edgel   =   face.edgel
        edgev̂   =   face.edgev̂
        edgen̂   =   face.edgen̂
        # 边数
        nEdge   =   length(edgel)

        # 初始化置零
        ûC₂Ipn̂dC₂Kvec   .=  0
        # Kᵣm3 
        local diKᵣm3  =   zero(FT)
        # 对 nEdge 个边循环
        for edgei in 1:nEdge
            # 置零避免累加错误
            Ilᵣ .=  0
            # 构成该边的第一个点
            @views edgeNodei⁻  =    face.vertices[:, MoM_Basics.EDGEVmINTriVsID[edgei]]
            # 该边边长
            lj          =   edgel[edgei]
            # 每个边的局部坐标下的 ljⱼ⁺ \lⱼ⁻
            lⱼ⁻  =   zero(FT)
            for ii in 1:3
                lⱼ⁻ += (edgeNodei⁻[ii] - r0[ii]) * edgev̂[ii, edgei]
            end
            lⱼ⁺ =   lⱼ⁻ + lj
            # 投影点 r0 到各个边的垂足的向量 p02ivec
            # @views p02ivec .=   edgeNodei⁻ .- lⱼ⁻ * edgev̂[:, edgei] .- r0
            for ii in 1:3
                p02ivec[ii] = edgeNodei⁻[ii] - lⱼ⁻ * edgev̂[ii, edgei] - r0[ii]
            end
            # 该向量长的平方、该向量长
            @views p02jl   =   p02ivec ⋅ edgen̂[:, edgei]
            p02jl²  =   abs2(p02jl)
            # 将p02ivec化为单位向量
            p02ivec ./=   p02jl
            # 场点到在该边上投影点的距离
            R0²     =   p02jl² + dts²
            # 场点到在该边上正、负端点的距离
            R⁺      =   sqrt(abs2(lⱼ⁺) + R0²)
            R⁻      =   sqrt(abs2(lⱼ⁻) + R0²)
            
            # 此处为避免数值误差，将阈值 ϵl 设定为 1e-2 边长
            let fⱼ = zero(FT), βⱼ =  zero(FT), ϵl = 1e-2lj
                if abs(p02jl) < ϵl
                    # 投影点过于靠近该边即 p02jl = 0 时， βⱼ = 0
                    if dtsAbs < ϵl
                        # 视为积分点与边重合，此时 βⱼ = 0， Kᵣ不变
                        continue
                    else
                        # 视为投影点与边重合，此时 βⱼ = 0， Kᵣ不变
                        fⱼ      =   log((lⱼ⁺ + R⁺)/(lⱼ⁻ + R⁻))
                    end
                else
                    fⱼ      =   log((lⱼ⁺ + R⁺)/(lⱼ⁻ + R⁻))
                    if dtsAbs < ϵl
                        # 视为积分点与源三角形面重合，此时 βⱼ ≠ 0
                        βⱼ      =   atan((p02jl*lⱼ⁺)/R0²) - atan((p02jl*lⱼ⁻)/R0²)
                        # Kᵣ累加
                        ISᵣ[-1] +=   p02jl*fⱼ
                        # diKᵣm3项累加
                        diKᵣm3  +=   βⱼ
                    else
                        # 正常处理
                        βⱼ      =   atan((p02jl*lⱼ⁺)/(R0² + dtsAbs*R⁺)) - atan((p02jl*lⱼ⁻)/(R0² + dtsAbs*R⁻))
                        # Kᵣ累加
                        ISᵣ[-1] +=   p02jl*fⱼ - dtsAbs*βⱼ
                        # diKᵣm3项累加
                        diKᵣm3  +=   βⱼ
                    end
                end #if p02jl
                # 计算 Ilᵣ 所有项
                Ilᵣ[-1]  =   fⱼ
                Ilᵣ[0]   =   lj
                R⁺ⁿ =   one(FT)
                R⁻ⁿ =   one(FT) 
                for n in 1:(SglrOrder-2)
                    R⁺ⁿ    *=   R⁺
                    R⁻ⁿ    *=   R⁻
                    Ilᵣ[n]  =   (lⱼ⁺*R⁺ⁿ - lⱼ⁻*R⁻ⁿ  + n*R0²*Ilᵣ[n-2])/(n+1)
                end # n
                # 计算 ISᵣ 中与边相关的项
                for n in 1:(SglrOrder-1)
                    ISᵣ[n]  +=   p02jl * Ilᵣ[n]
                end # n
                # 计算 C2I
                # C2I     =    CT(Ilᵣ[-1])
                # for n in eachindex(VSC₂ⁿ)
                #     C2I -=  VSC₂ⁿ[n] * Ilᵣ[n + 2]/(n + 2)
                # end
                C2I =  Ilᵣ ⋅ VSC₃ⁿ

                # 累加 ûC₂Ipn̂dC₂Kvec 与边有关的项
                @views ûC₂Ipn̂dC₂Kvec   .+=  edgen̂[:, edgei] .* C2I
            end #let
        end # edgei
        # 其它 ISᵣ 相关项
        diKᵣm3  *=   sign(dts)
        ISᵣ[0]   =   abs(volumeCell.facesArea[iface])
        for n in 1:(SglrOrder-1)
            ISᵣ[n]  +=  n*dts²*ISᵣ[n-2]
            ISᵣ[n]  /=  n + 2
        end
        #= 累加面 格林函数在 面上的 奇异性计算结果, 注意此处点乘不能改变顺序
        Julia 程序处理 点乘时，若前面一个向量为复数，则取其共轭，因此，要将 浮点数 向量 ISᵣ 放在前面 =#
        dC₁K    +=  dts*div4π*(ISᵣ ⋅ VSC₁ⁿ)

        # 累加 ûC₂Ipn̂dC₂Kvec 与面有关的项
        @views ûC₂Ipn̂dC₂Kvec    .+=  (diKᵣm3 + dts*(ISᵣ[-1:(SglrOrder-3)] ⋅ VSC₂ⁿ)) .* volumeCell.facesn̂[:, iface]
        ûC₂Ipn̂dC₂Kvec .*= div4π
        
        @views n̂ûC₂Ipn̂dC₂KDyad  .+=  volumeCell.facesn̂[:, iface] * conj(ûC₂Ipn̂dC₂Kvec')
    end # iface

    re  .=   ( I * (-dC₁K*k²)) + n̂ûC₂Ipn̂dC₂KDyad

    return re

end
function volumeSingularityLOpDyad(rtveclc::AbstractVector{FT}, volumeCell::HexahedraInfo{IT, FT, CT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}}
    # 体元的面
    faces   =   volumeCell.faces
    # 体元的面数
    nFaces  =   length(faces)
    # 常数项
    k²      =   Params.k²
    # 预分配 Iᵣⁿ、 ISᵣ, K̂ᵣⁿ 结果数组
    ISᵣ     =   OffsetArray(zero(MVector{SglrOrder+1, FT}), -2)
    Ilᵣ     =   OffsetArray(zero(MVector{SglrOrder+1, FT}), -2)
    # 一些要计算的量
    dC₁K::CT    =   zero(CT)
    ûC₂Ipn̂dC₂Kvec   =   zero(MVec3D{CT})
    n̂ûC₂Ipn̂dC₂KDyad =   zero(MMatrix{3, 3, CT})
    re              =   zero(MMatrix{3, 3, CT})
    p02ivec         =   zero(MVec3D{FT})
    r0              =   zero(MVec3D{FT})

    # 对高斯求积点循环
    for iface in 1:nFaces
        # 第 iface 个面
        face    =   faces[iface]
        # 置零
        ISᵣ .=  0
        # 该积分点在到源三角形上的距离（带正负，以源三角形法向为正向）
        @views dts  =   volumeCell.facesn̂[:, iface] ⋅ (rtveclc .- face.vertices[:, 1])
        # 投影点
        # @views r0   =   rtveclc .- dts .* volumeCell.facesn̂[:, iface]
        for ii in 1:3
            r0[ii]   =   rtveclc[ii] .- dts .* volumeCell.facesn̂[ii, iface]
        end
        # 距离的绝对值、平方
        dts²    =   dts*dts
        dtsAbs  =   abs(dts)
        
        # 边的信息
        edgel   =   face.edgel
        edgev̂   =   face.edgev̂
        edgen̂   =   face.edgen̂
        # 边数
        nEdge   =   length(edgel)

        # 初始化置零
        ûC₂Ipn̂dC₂Kvec   .=  0
        # Kᵣm3 
        local diKᵣm3  =   zero(FT)
        # 对 nEdge 个边循环
        for edgei in 1:nEdge
            # 置零避免累加错误
            Ilᵣ .=  0
            # 构成该边的第一个点
            @views edgeNodei⁻  =    face.vertices[:, edgei]
            # 该边边长
            lj          =   edgel[edgei]
            # 每个边的局部坐标下的 ljⱼ⁺ \lⱼ⁻
            # @views lⱼ⁻  =   (edgeNodei⁻ .- r0) ⋅ edgev̂[:, edgei]
            @views lⱼ⁻  =   zero(FT)
            for ii in 1:3
                lⱼ⁻ += (edgeNodei⁻[ii] - r0[ii]) * edgev̂[ii, edgei]
            end
            lⱼ⁺ =   lⱼ⁻ + lj
            # 投影点 r0 到各个边的垂足的向量 p02ivec
            # @views p02ivec .=   edgeNodei⁻ .- lⱼ⁻ * edgev̂[:, edgei] .- r0
            for ii in 1:3
                p02ivec[ii] = edgeNodei⁻[ii] - lⱼ⁻ * edgev̂[ii, edgei] - r0[ii]
            end
            # 该向量长的平方、该向量长
            @views p02jl  =   p02ivec ⋅ edgen̂[:, edgei]
            p02jl²  =   abs2(p02jl)
            # 将p02ivec化为单位向量
            p02ivec ./=   p02jl
            # 场点到在该边上投影点的距离
            R0²     =   p02jl² + dts²
            # 场点到在该边上正、负端点的距离
            R⁺      =   sqrt(abs2(lⱼ⁺) + R0²)
            R⁻      =   sqrt(abs2(lⱼ⁻) + R0²)
            
            # 此处为避免数值误差，将阈值 ϵl 设定为 1e-2 边长
            let fⱼ = zero(FT), βⱼ =  zero(FT), ϵl = lj/100
                if abs(p02jl) < ϵl
                    # 投影点过于靠近该边即 p02jl = 0 时， βⱼ = 0
                    if dtsAbs < ϵl
                        # 视为积分点与边重合，此时 βⱼ = 0， Kᵣ不变
                        continue
                    else
                        # 视为投影点与边重合，此时 βⱼ = 0， Kᵣ不变
                        fⱼ      =   log((lⱼ⁺ + R⁺)/(lⱼ⁻ + R⁻))
                    end
                else
                    fⱼ      =   log((lⱼ⁺ + R⁺)/(lⱼ⁻ + R⁻))
                    if dtsAbs < ϵl
                        # 视为积分点与源三角形面重合，此时 βⱼ ≠ 0
                        βⱼ      =   atan((p02jl*lⱼ⁺)/R0²) - atan((p02jl*lⱼ⁻)/R0²)
                        # Kᵣ累加
                        ISᵣ[-1] +=   p02jl*fⱼ
                        # diKᵣm3项累加
                        diKᵣm3  +=   βⱼ
                    else
                        # 正常处理
                        βⱼ      =   atan((p02jl*lⱼ⁺)/(R0² + dtsAbs*R⁺)) - atan((p02jl*lⱼ⁻)/(R0² + dtsAbs*R⁻))
                        # Kᵣ累加
                        ISᵣ[-1] +=   p02jl*fⱼ - dtsAbs*βⱼ
                        # diKᵣm3项累加
                        diKᵣm3  +=   βⱼ
                    end
                end #if p02jl
                # 计算 Ilᵣ 所有项
                Ilᵣ[-1]  =   fⱼ
                Ilᵣ[0]   =   lj
                R⁺ⁿ =   one(FT)
                R⁻ⁿ =   one(FT) 
                for n in 1:(SglrOrder-2)
                    R⁺ⁿ    *=   R⁺
                    R⁻ⁿ    *=   R⁻
                    Ilᵣ[n]  =   (lⱼ⁺*R⁺ⁿ - lⱼ⁻*R⁻ⁿ  + n*R0²*Ilᵣ[n-2])/(n+1)
                end # n
                # 计算 ISᵣ 中与边相关的项
                for n in 1:(SglrOrder-1)
                    ISᵣ[n]  +=   p02jl * Ilᵣ[n]
                end # n
                # 计算 C2I
                # C2I     =    CT(Ilᵣ[-1])
                # for n in eachindex(VSC₂ⁿ)
                #     C2I -=  VSC₂ⁿ[n] * Ilᵣ[n + 2]/(n + 2)
                # end
                C2I =  Ilᵣ ⋅ VSC₃ⁿ

                # 累加 ûC₂Ipn̂dC₂Kvec 与边有关的项
                @views ûC₂Ipn̂dC₂Kvec   .+=  edgen̂[:, edgei] .* C2I
            end #let
        end # edgei
        # 其它 ISᵣ 相关项
        diKᵣm3  *=   sign(dts)
        ISᵣ[0]   =   abs(volumeCell.facesArea[iface])
        for n in 1:(SglrOrder-1)
            ISᵣ[n]  +=  n*dts²*ISᵣ[n-2]
            ISᵣ[n]  /=  n + 2
        end
        #= 累加面 格林函数在 面上的 奇异性计算结果, 注意此处点乘不能改变顺序
        Julia 程序处理 点乘时，若前面一个向量为复数，则取其共轭，因此，要将 浮点数 向量 ISᵣ 放在前面 =#
        dC₁K    +=  dts*div4π*(ISᵣ ⋅ VSC₁ⁿ)

        # 累加 ûC₂Ipn̂dC₂Kvec 与面有关的项
        @views ûC₂Ipn̂dC₂Kvec    .+=  (diKᵣm3 + dts*(ISᵣ[-1:(SglrOrder-3)] ⋅ VSC₂ⁿ)) .* volumeCell.facesn̂[:, iface]
        ûC₂Ipn̂dC₂Kvec .*= div4π
        
        @views n̂ûC₂Ipn̂dC₂KDyad  .+=  volumeCell.facesn̂[:, iface] * conj(ûC₂Ipn̂dC₂Kvec')

    end # iface

    re  .=   ( I * (-dC₁K*k²)) + n̂ûC₂Ipn̂dC₂KDyad

    return re

end