using ManifoldMeshes
using Test

@testset "CSRMapping basics" begin
    # Uniform CSR
    csr = ManifoldMeshes.CSRMapping(3, 4)
    @test length(csr) == 3
    @test ManifoldMeshes.n_neighbors(csr, 1) == 4
    @test ManifoldMeshes.n_neighbors(csr, 2) == 4

    # Fill and read back via getindex_fixed
    csr.values[1:4] = [10, 20, 30, 40]
    csr.values[5:8] = [50, 60, 70, 80]
    csr.values[9:12] = [90, 100, 110, 120]

    @test ManifoldMeshes.getindex_fixed(csr, 1, Val(4)) == (10, 20, 30, 40)
    @test ManifoldMeshes.getindex_fixed(csr, 2, Val(4)) == (50, 60, 70, 80)

    # Variable CSR
    counts = [2, 3, 1]
    csr2, ptrs = ManifoldMeshes.CSRMapping(3, counts)
    @test length(csr2) == 3
    @test ManifoldMeshes.n_neighbors(csr2, 1) == 2
    @test ManifoldMeshes.n_neighbors(csr2, 2) == 3
    @test ManifoldMeshes.n_neighbors(csr2, 3) == 1

    # Write via ptrs
    csr2.values[ptrs[1]] = 1;
    ptrs[1] += 1
    csr2.values[ptrs[1]] = 2;
    ptrs[1] += 1
    csr2.values[ptrs[2]] = 3;
    ptrs[2] += 1
    csr2.values[ptrs[2]] = 4;
    ptrs[2] += 1
    csr2.values[ptrs[2]] = 5;
    ptrs[2] += 1
    csr2.values[ptrs[3]] = 6;
    ptrs[3] += 1

    @test collect(csr2[1]) == [1, 2]
    @test collect(csr2[2]) == [3, 4, 5]
    @test collect(csr2[3]) == [6]

    # Empty neighbor
    counts0 = [0, 2]
    csr0, _ = ManifoldMeshes.CSRMapping(2, counts0)
    @test length(csr0[1]) == 0
    @test length(csr0[2]) == 2
end

@testset "CSR getindex returns view" begin
    csr = ManifoldMeshes.CSRMapping(3, 2)
    csr.values[1:2] = [10, 20]
    csr.values[3:4] = [30, 40]
    csr.values[5:6] = [50, 60]

    v = csr[2]
    @test isa(v, SubArray)  # view, not copy
    @test v == [30, 40]

    # Mutating the view should mutate the underlying CSR
    v[1] = 999
    @test csr.values[3] == 999
end

@testset "CSR uniform permutation" begin
    csr = ManifoldMeshes.CSRMapping(3, 2)
    csr.values[1:2] = [10, 20]
    csr.values[3:4] = [30, 40]
    csr.values[5:6] = [50, 60]

    perm = [3, 1, 2]  # new_id -> old_id
    perm_csr = ManifoldMeshes._permute_uniform_csr(csr, perm, Val(2))

    @test perm_csr[1] == [50, 60]  # was old 3
    @test perm_csr[2] == [10, 20]  # was old 1
    @test perm_csr[3] == [30, 40]  # was old 2
end
