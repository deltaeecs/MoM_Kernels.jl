function test_matrix_solving(geosInfo, bfsInfo)
    ## 快速算法
    # 计算阻抗矩阵（MLFMA计算），注意此处根据基函数在八叉树的位置信息改变了基函数顺序
    # Zopt, octree, ZnearCSC  =   getImpedanceOpt(geosInfo, bfsInfo);
    nLevels, octree     =   getOctreeAndReOrderBFs!(geosInfo, bfsInfo)
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

    source  =   PlaneWave(π, 0, 0f0, 1f0)

    V    =   getExcitationVector(geosInfo, size(ZnearCSC, 1), source);
    @test true

    ICoeff, ch   =   solve(Zopt, V; solverT = :gmres, Pl = Zprel);
    @test true

    fill!(ICoeff, 0)
    ch   =   solve!(Zopt, ICoeff, V; solverT = :gmres!, Pl = Zprel);
    @test true

    ICoeff

end