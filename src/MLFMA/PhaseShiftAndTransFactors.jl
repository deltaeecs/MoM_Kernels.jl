## 本文件用于预计算子层盒子向父层的相移因子和本层盒子的远亲之间的转移因子

"""
本函数用于给输入的本(level)层的盒子与其子盒子之间计算相移因子，
由盒子排列的规律性和相移因子的对称性，可知：
只需要计算8个相移因子，即可用于所有盒子到其子盒子的相移，
且这八个盒子关于原点对称的两两之间的相移因子为共轭关系
计算完成直接保存在 level 不再返回
"""
function setLevelsShiftFactor!(nLevels::Int, levels::Dict{Int, LV}) where{LV<:AbstractLevel}

    FT = typeof(levels[nLevels].cubeEdgel)
    # 常数
    JK_0 = Params.JK_0
    # 从 第“2”层开始计算到 nLevels-1 层计算子层到本层的相移因子
    for iLevel in 2:(nLevels - 1)
        # 该层
        level = levels[iLevel]

        # z轴以下4个子盒子的偏置，其他四个通过共轭对称计算
        kidsOffsets = [  -1 -1 -1;  1  -1 -1; -1  1 -1; 1  1 -1;]'
        # 本层盒子中心到子层盒子中心之间的 对应的 4个 偏置向量
        ΔCt2Cks  =   SMatrix{3, 4, FT}(level.cubeEdgel/4 .* kidsOffsets)

        # 多极子采样信息
        polesr̂sθsϕs =   level.poles.r̂sθsϕs
        nPoles  =   length(polesr̂sθsϕs)

        # 预分配内存
        phaseShift2Kids     =   zeros(Complex{FT}, nPoles, 8)
        
        # 开始循环计算
        for iKid in 1:4
            # 到子盒子的偏置
            ΔCt2Ck  =   ΔCt2Cks[:,iKid]
            # 开始对子盒子循环
            for iPole in 1:nPoles
                # 多极子采样点信息
                r̂θϕ =  polesr̂sθsϕs[iPole]
                # k̂
                k̂   =  r̂θϕ.r̂
                # 相移因子
                phaseShift2Kids[iPole, iKid] = exp(-JK_0*(k̂ ⋅ ΔCt2Ck))
            end # iPole
        end # iKid
        # 对称性计算另外四个子盒子
        phaseShift2Kids[:, 8:-1:5]   .=   conj(phaseShift2Kids[:, 1:4])

        # 取共轭计算从子层盒子到本盒子的相移因子
        phaseShiftFromKids  =   conj(phaseShift2Kids)

        # 保存在层信息中
        level.phaseShift2Kids      =   phaseShift2Kids
        level.phaseShiftFromKids   =   phaseShiftFromKids

    end # iLevel
    return nothing
end

"""
第一类球汉克尔函数，使用 SpecialFunctions.jl，
适用于非整数阶、复数变量，算的较慢，只在计算有耗介质（复数波矢）时调用
"""
function spherical_h1l(l, x::T) where T
    sphericalbesselj(l,x) + im*sphericalbessely(l,x)
end

"""
第一类球汉克尔函数，使用GSL.jl(GNU Scientific Library)，适用于 l 为整数，x 为浮点数时算的更快
"""
spherical_h1l(l::Integer,x::Real) = sf_bessel_jl(l, x) + im*sf_bessel_yl(l, x)

"""
一次计算 0:lmax 的多阶第一类球汉克尔函数， 保存在数组里
"""
spherical_h1l_array(l, x::T) where T =   [spherical_h1l(l, x::T) for l in 0:lmax]

"""
一次计算 0:lmax 的多阶第一类球汉克尔函数， 保存在数组里
"""
spherical_h1l_array(lmax::Integer,x::T)  where {T<:Real}   =   Complex{T}.(sf_bessel_jl_array(lmax, x) .+ im*sf_bessel_yl_array(lmax, x))


"""
第二类球汉克尔函数，使用 SpecialFunctions.jl，
适用于非整数阶、复数变量，算的较慢，只在计算有耗介质（复数波矢）时调用
"""
function spherical_h2l(l, x::T) where T
    sphericalbesselj(l,x) - im*sphericalbessely(l,x)
end

"""
第二类球汉克尔函数，使用GSL.jl(GNU Scientific Library)，，适用于 l 为整数，x 为浮点数时算的更快
"""
spherical_h2l(l::Integer,x::Real) = sf_bessel_jl(l, x) - im*sf_bessel_yl(l, x)

"""
一次计算 0:lmax 的多阶第二类球汉克尔函数， 保存在数组里
"""
spherical_h2l_array(lmax, x::T) where T =   [spherical_h2l(l, x::T) for l in 0:lmax]


"""
一次计算 0:lmax 的多阶第二类球汉克尔函数， 保存在数组里
"""
spherical_h2l_array(lmax::Integer,x::T) where {T<:Real}     =   Complex{T}.(sf_bessel_jl_array(lmax, x) .- im*sf_bessel_yl_array(lmax, x))


"""
计算 第“2”层 到 叶 层的转移因子，
转移因子只存在于远亲组，每层远亲组最多有 7^3 - 3^3 = 316种结果
"""
function setLevelTransFactor!(nLevels::Int, levels::Dict{Int, LV}) where{LV<:AbstractLevel}

    # 预分配所有的316个远亲的id
    all316FarNeighID    =   zeros(Int32, 3, 316)
    # 所有343个 7*7*7 盒子id到这316个盒子的索引
    all343InFar316      =   OffsetArray(zeros(Int32, 7, 7, 7), -3:3, -3:3, -3:3)

    # 开始找出远亲组
    indexfar316 =   0
    for k in -3:3
        for j in -3:3
            for i in -3:3
                if (abs(i) > 1) | (abs(j) > 1) | (abs(k) > 1)
                    indexfar316 += 1
                    all343InFar316[i, j, k]             =   indexfar316
                    all316FarNeighID[:, indexfar316]   .=   [i, j, k]
                end
            end # for i
        end # for j
    end # for k

    
    # 进度条
    pmeter = Progress(length(2:nLevels); desc = "Calculating translation factors...")

    # 从 第“2”层开始计算到 叶 层计算子层到本层的转移因子
    for iLevel in 2:nLevels
        # 该层
        level = levels[iLevel]
        # 该层截断项
        truncL  =   level.L
    
        calαTransOnLevel!(level, truncL, all316FarNeighID, all343InFar316)
        
        # 更新进度条
        next!(pmeter)

    end # for iLevel

    return nothing

end

"""
计算 level 层的转移因子，
转移因子只存在于远亲组，每层远亲组最多有 7^3 - 3^3 = 316种结果
"""
function calαTransOnLevel!(level, truncL, all316FarNeighID, all343InFar316)

    # 常数
    mjKdiv16π²  =  -Params.JK_0/(4π)^2
    K_0         =   Params.K_0
    # 多极子采样信息
    polesr̂sθsϕs =   level.poles.r̂sθsϕs
    # 权重Weights
    Wθϕs        =   level.poles.Wθϕs
    # 多极子数
    nPoles  =   length(polesr̂sθsϕs)

    # 浮点数类型
    FT = typeof(level.cubeEdgel)

    # 预分配内存
    αTrans  =   zeros(Complex{FT}, nPoles, 316)

    pmeter  =   Progress(316; desc = "Calculating translation factors on level $(level.ID)...")
    
    @floop WorkStealingEx() for iFarNei in 1:316
        # 本层盒子中心到远亲盒子中心之间的 对应的 316 个 偏置向量
        RabVec  =   SVec3D{FT}(level.cubeEdgel .* all316FarNeighID[:, iFarNei])
        # 本层盒子中心到远亲盒子中心之间的距离
        Rab     =   norm(RabVec)
        # R̂ab
        R̂ab     =   RabVec / Rab
        # 预分配汉克尔函数内存
        h2lxs   =   zeros(Complex{FT}, truncL+1)
        h2lxsOffset =   OffsetArray(h2lxs, 0:truncL)

        for iPole in 1:nPoles
            # 多极子采样点信息
            r̂θϕ =  polesr̂sθsϕs[iPole]
            # k̂
            k̂   =  r̂θϕ.r̂
            # k̂ ⋅ R̂ab
            cosϕ = clamp(k̂ ⋅ R̂ab, -1., 1.)

            # 勒让德多项式在每一次的累加循环里计算，此处直接计算出所有的 (0-truncL) 为索引
            legendrePls =   collectPl(cosϕ, lmax = truncL)
            # 计算出所有的汉克尔函数值，并保存在从 0 开始索引的偏置数组
            h2lxs      .=   spherical_h2l_array(truncL, K_0*Rab)
            # 置零用于累加
            αTransTemp  =   zero(Complex{FT})
            # 累加计算
            jˡ = im
            for l in 0:truncL
                jˡ *= -im
                αTransTemp +=  jˡ*(2l + 1)*h2lxsOffset[l]*legendrePls[l]
            end #for l
            # 系数修正
            αTransTemp *=  mjKdiv16π² * Wθϕs[iPole]
            # 记录结果
            αTrans[iPole, iFarNei]  =   αTransTemp
        end # for iPole
        # 更新进度条
        next!(pmeter)
    end # for iFarNei

    # 写入层信息
    level.αTrans        =   αTrans
    level.αTransIndex   =   all343InFar316

    return nothing

end


function calαTrans(truncL, θ, ϕ, RVec)
    # k̂
    k̂   =   r̂func(θ, ϕ)
    Rab =   norm(RVec)
    R̂ab =   RVec ./ Rab
    # k̂ ⋅ R̂ab
    cosϕ    =   k̂ ⋅ R̂ab

    # 勒让德多项式在每一次的累加循环里计算，此处直接计算出所有的 (0-truncL) 为索引
    legendrePls =   collectPl(cosϕ, lmax = truncL)
    # 计算出所有的汉克尔函数值，并保存在从 0 开始索引的偏置数组
    h2lxs       =   OffsetArray(spherical_h2l_array(truncL, K_0*Rab), 0:truncL)
    # 置零用于累加
    αTransTemp  =   zero(ComplexF64)
    # 累加计算
    jˡ = im
    for l in 0:truncL
        jˡ *= -im
        αTransTemp +=  jˡ*(2l + 1)*h2lxs[l]*legendrePls[l]
    end #for l
    # 系数修正
    αTransTemp  *=  -Params.JK_0/(4π)^2

    return αTransTemp
end