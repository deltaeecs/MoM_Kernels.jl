function test_mat_solving(geosInfo, nbf)
    
    # 阻抗矩阵
    Zmat = getImpedanceMatrix(geosInfo, nbf)

    # 源
    source  =   PlaneWave(π, 0, 0f0, 1f0)
    # 激励向量
    V    =   getExcitationVector(geosInfo, nbf, source);
    @test true

    # 迭代求解
    ICoeff, _   =   solve(Zmat, V; solverT = :direct);
    @test true

    ICoeff

end

function test_opt_solving(geosInfo, bfsInfo; source = nothing)
    ## 快速算法
    # 计算阻抗矩阵（MLFMA计算），注意此处根据基函数在八叉树的位置信息改变了基函数顺序
    # Zopt, octree, ZnearCSC  =   getImpedanceOpt(geosInfo, bfsInfo);
    nLevels, octree     =   getOctreeAndReOrderBFs!(geosInfo, bfsInfo; nInterp = 4)
    @test true

    # 叶层
    leafLevel   =   octree.levels[nLevels];
    # 计算近场矩阵CSC
    ZnearCSC     =   calZnearCSC(leafLevel, geosInfo, bfsInfo);
    @test true

    # 构建矩阵向量乘积算子
    Zopt  =   MLMFAIterator(ZnearCSC, octree, geosInfo, bfsInfo);
    @test true

    ## 根据近场矩阵和八叉树计算 SAI 左预条件
    Zprel    =   sparseApproximateInversePl(ZnearCSC, leafLevel)
    @test true

    # 源
    isnothing(source) && begin
        source  =   MagneticDipole(;Iml = 1., phase = 0., orient = (π/2, π/2, 0))#PlaneWave(π, 0, 0f0, 1f0)
    end
    # 激励向量
    V    =   getExcitationVector(geosInfo, size(ZnearCSC, 1), source);
    @test true

    # 迭代求解
    ICoeff, ch   =   solve(Zopt, V; solverT = :gmres, Pl = Zprel);
    @test true

    # 原地迭代求解
    fill!(ICoeff, 0)
    ch   =   solve!(Zopt, ICoeff, V; solverT = :gmres!, Pl = Zprel);
    @test true

    ICoeff

end