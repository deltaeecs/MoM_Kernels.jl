
"""
在球坐标为r̂θϕ处计算辐射积分，采用RWG基函数时，三角形上没有统一的电流值，每一点上都是三边电流的叠加，
此时:
N(θ, ϕ) =   ∑ₙ(∫ₛ Jˢ exp( jkr̂(θ, ϕ)⋅rₙ ) dS)
        =   ∑ₙ(∫ₛ (∑ₜₙ₌₁³ Iₙfₙ)exp(jkr̂(θ, ϕ)⋅rₙ) dS)
        =   ∑ₙ(Sₜ (∑ₜₙ₌₁³ Iₙlₙρₙ/(2Sₙ))exp(jkr̂(θ, ϕ)⋅rₙ) )
        =   ∑ₙ(∑ᵢWᵢ(∑ₜₙ₌₁³ Iₙlₙρₙ/2)exp(jkr̂(θ, ϕ)⋅rₙ) )
"""
function radarCrossSection(θs_obs, ϕs_obs, ICoeff::Vector{CT}, 
    trianglesInfo::Vector{ST}, ::Type{BFT} = VSBFTypes.sbfType; str = "") where {CT<:Complex, ST<:TriangleInfo, BFT<:RWG}
    FT = Precision.FT
    # 高斯求积点电流权重乘积
    Jtris       =   electricJCal(ICoeff, trianglesInfo)
    # 常数
    K_0         =   Params.K_0
    Nθ_obs      =   length(θs_obs)
    Nϕ_obs      =   length(ϕs_obs)

    # 观测角度信息
    nobs        =   Nθ_obs * Nϕ_obs
    θsobsInfo   =   [∠Info{FT}(θ_obs) for θ_obs in θs_obs]
    ϕsobsInfo   =   [∠Info{FT}(ϕ_obs) for ϕ_obs in ϕs_obs]
    r̂θsϕs       =   [r̂θϕInfo(θobsInfo, ϕobsInfo) for θobsInfo in θsobsInfo, ϕobsInfo in ϕsobsInfo]
    
    # 预分配RCS内存
    RCSθsϕs     =   zeros(FT, (2, length(θsobsInfo), length(ϕsobsInfo)))
    RCSθsϕsrsp  =   reshape(RCSθsϕs, (2, nobs))
    
    # 进度条
    pmeter  =   Progress(nobs, "Calculating RCS ($Nθ_obs × $Nϕ_obs)")
    # 计算RCS
    @threads for ii in 1:nobs
        # 辐射积分
        Nθϕ     =   raditionalIntegralNθϕCal(r̂θsϕs[ii], trianglesInfo, Jtris)
        # RCS
        RCSθϕ   =   abs2.(Nθϕ)  
        RCSθsϕsrsp[:, ii]   .=  (K_0*η_0)^2/4π*RCSθϕ
        # 更新进度条
        next!(pmeter)
    end #for ii
    # dB形式
    RCSθsϕsdB   =   10log10.(RCSθsϕs)
    # 总的RCS
    RCS         =   RCSθsϕs[1, :, :] + RCSθsϕs[2, :, :]
    RCSdB       =   10log10.(RCS)
    # 保存数据
    save_RCS(θs_obs, ϕs_obs, RCS; str=str)
    # 绘图
    RCSPlot(θs_obs, ϕs_obs, RCS, RCSdB)
    # 返回
    return RCSθsϕs, RCSθsϕsdB, RCS, RCSdB

end # end function


"""
在球坐标为r̂θϕ处计算辐射积分，采用SWG基函数时，四面体上没有统一的电流值，每一点上都是四个SWG基函数电流的叠加，
此时:
N(θ, ϕ) =   ∑ₙ(∫ₛ Jˢ exp( jkr̂(θ, ϕ)⋅rₙ ) dS)
        =   ∑ₙ(∫ₛ (∑ₜₙ₌₁³ Iₙfₙ)exp(jkr̂(θ, ϕ)⋅rₙ) dS)
        =   ∑ₙ(Sₜ (∑ₜₙ₌₁³ Iₙsₙρₙ/(3Vₙ))exp(jkr̂(θ, ϕ)⋅rₙ) )
        =   ∑ₙ(∑ᵢWᵢ(∑ₜₙ₌₁³ Iₙsₙρₙ/3)exp(jkr̂(θ, ϕ)⋅rₙ) )
"""
function radarCrossSection(θs_obs, ϕs_obs, ICoeff::Vector{CT}, geosInfo::Vector{VT}, 
                            bfT::Type{BFT} = VSBFTypes.vbfType; str = "") where {VT<:VolumeCellType, CT<:Complex, BFT<:BasisFunctionType}
    FT = Precision.FT
    # 高斯求积点电流权重乘积
    Jgeos       =   geoElectricJCal(ICoeff, geosInfo, bfT)
    # 常数
    K_0         =   Params.K_0

    # 观测角度信息
    nobs        =   length(θs_obs) * length(ϕs_obs)
    θsobsInfo   =   [∠Info{FT}(θ_obs) for θ_obs in θs_obs]
    ϕsobsInfo   =   [∠Info{FT}(ϕ_obs) for ϕ_obs in ϕs_obs]
    r̂θsϕs       =   [r̂θϕInfo(θobsInfo, ϕobsInfo) for θobsInfo in θsobsInfo, ϕobsInfo in ϕsobsInfo]
    
    # 预分配RCS内存
    RCSθsϕs     =   zeros(FT, (2, length(θsobsInfo), length(ϕsobsInfo)))
    RCSθsϕsrsp  =   reshape(RCSθsϕs, (2, nobs))
    
    # 进度条
    pmeter  =   Progress(nobs, "Calculating RCS ($(length(θs_obs)) × $(length(ϕs_obs))))")
    # 计算RCS
    @threads for ii in 1:nobs
        # 辐射积分
        Nθϕ     =   raditionalIntegralNθϕCal(r̂θsϕs[ii], geosInfo, Jgeos)
        # RCS
        RCSθϕ   =   abs2.(Nθϕ)
        RCSθsϕsrsp[:, ii]   .=  (K_0*η_0)^2/4π*RCSθϕ
        # 更新进度条
        next!(pmeter)
    end #for ii
    # dB形式
    RCSθsϕsdB   =   10log10.(RCSθsϕs)
    # 总的RCS
    @views RCS  =   RCSθsϕs[1, :, :] .+ RCSθsϕs[2, :, :]
    RCSdB       =   10log10.(RCS)
    # 保存数据
    save_RCS(θs_obs, ϕs_obs, RCS; str=str)
    # 绘图
    RCSPlot(θs_obs, ϕs_obs, RCS, RCSdB)
    # 返回
    return RCSθsϕs, RCSθsϕsdB, RCS, RCSdB

end # end function


"""
在球坐标为r̂θϕ处计算辐射积分，采用SWG基函数时，四面体上没有统一的电流值，每一点上都是四个SWG基函数电流的叠加，
此时:
N(θ, ϕ) =   ∑ₙ(∫ₛ Jˢ exp( jkr̂(θ, ϕ)⋅rₙ ) dS)
        =   ∑ₙ(∫ₛ (∑ₜₙ₌₁³ Iₙfₙ)exp(jkr̂(θ, ϕ)⋅rₙ) dS)
        =   ∑ₙ(Sₜ (∑ₜₙ₌₁³ Iₙsₙρₙ/(3Vₙ))exp(jkr̂(θ, ϕ)⋅rₙ) )
        =   ∑ₙ(∑ᵢWᵢ(∑ₜₙ₌₁³ Iₙsₙρₙ/3)exp(jkr̂(θ, ϕ)⋅rₙ) )
"""
function radarCrossSection(θs_obs, ϕs_obs,
    ICoeff::Vector{CT}, geosInfo::Vector{VT}; str = "") where {CT<:Complex, VT<:AbstractVector}
    FT = Precision.FT
    # 面网格、体网格
    tris    =   geosInfo[1]
    geosV   =   geosInfo[2]
    sbfT    =   getBFTfromCellT(eltype(tris))
    vbfT    =   getBFTfromCellT(eltype(geosV))
    # 高斯求积点电流权重乘积
    Jtris       =   geoElectricJCal(ICoeff, tris,  sbfT)
    JgeoVs      =   geoElectricJCal(ICoeff, geosV, vbfT)
    # 常数
    K_0         =   Params.K_0

    # 观测角度信息
    nobs        =   length(θs_obs) * length(ϕs_obs)
    θsobsInfo   =   [∠Info{FT}(θ_obs) for θ_obs in θs_obs]
    ϕsobsInfo   =   [∠Info{FT}(ϕ_obs) for ϕ_obs in ϕs_obs]
    r̂θsϕs       =   [r̂θϕInfo(θobsInfo, ϕobsInfo) for θobsInfo in θsobsInfo, ϕobsInfo in ϕsobsInfo]
    
    # 预分配RCS内存
    RCSθsϕs     =   zeros(FT, (2, length(θsobsInfo), length(ϕsobsInfo)))
    RCSθsϕsrsp  =   reshape(RCSθsϕs, (2, nobs))
    
    # 进度条
    pmeter  =   Progress(nobs, "Calculating RCS ($(length(θs_obs)) × $(length(ϕs_obs))))")
    # 计算RCS
    @threads for ii in 1:nobs
        # 辐射积分
        Nθϕ     =   raditionalIntegralNθϕCal(r̂θsϕs[ii], tris, Jtris)
        Nθϕ   .+=   raditionalIntegralNθϕCal(r̂θsϕs[ii], geosV, JgeoVs)
        # RCS
        RCSθϕ   =   abs2.(Nθϕ)
        RCSθsϕsrsp[:, ii]   .=  (K_0*η_0)^2/4π*RCSθϕ
        # 更新进度条
        next!(pmeter)
    end #for ii
    # dB形式
    RCSθsϕsdB   =   10log10.(RCSθsϕs)
    # 总的RCS
    @views RCS  =   RCSθsϕs[1, :, :] .+ RCSθsϕs[2, :, :]
    RCSdB       =   10log10.(RCS)
    # 保存数据
    save_RCS(θs_obs, ϕs_obs, RCS; str=str)
    # 绘图
    RCSPlot(θs_obs, ϕs_obs, RCS, RCSdB)

    # 返回
    return RCSθsϕs, RCSθsϕsdB, RCS, RCSdB

end # end function

function save_RCS(θs_obs, ϕs_obs, RCS; str="")
    # 保存数据
    open(SimulationParams.resultDir*"farEm2$str.txt", "w") do io
        θs_obs_deg  =   θs_obs/pi*180
        ϕs_obs_deg  =   ϕs_obs/pi*180
        for θii in 1:length(θs_obs_deg)
            write(io, "$(θs_obs_deg[θii])\t")
            for ϕjj in 1:length(ϕs_obs_deg)
                write(io, "$(RCS[θii, ϕjj])\t")
            end
            write(io, "\n")
        end
    end
end

"""
RCS 绘图
"""
function RCSPlot(θs_obs, ϕs_obs, RCS::Matrix{FT}, RCSdB::Matrix{FT};  str = "") where{FT<:Real}
    
    θs_obs_deg  =   θs_obs/pi*180
    ϕs_obs_deg  =   ϕs_obs/pi*180 #[2:end-1]
    # 标签
    labels  =   reshape(["ϕ = $(ϕ_obs_deg)°" for ϕ_obs_deg in ϕs_obs_deg], (1, length(ϕs_obs_deg)))
    # 绘图
    θs_obs_deg  =   θs_obs/pi*180
    ϕs_obs_deg  =   ϕs_obs/pi*180
    # 标签
    labels  =   reshape(["ϕ = $(ϕ_obs_deg)°" for ϕ_obs_deg in ϕs_obs_deg], (1, length(ϕs_obs_deg)))
    # 绘图
    figRCSm²  = lineplot(θs_obs_deg, RCS; name = labels, xlim = extrema(θs_obs_deg), ylim = extrema(RCS), xlabel="θ", ylabel="m²", title = "RCS(m²)(θϕ)$str")
    figRCSdB  = lineplot(θs_obs_deg, RCSdB; name = labels, xlim = extrema(θs_obs_deg), ylim = extrema(RCSdB), xlabel="θ", ylabel="dB", title = "RCS(dB)(θϕ)$str")

    SimulationParams.SHOWIMAGE && display(figRCSm²)
    SimulationParams.SHOWIMAGE && display(figRCSdB)

    return figRCSdB

end